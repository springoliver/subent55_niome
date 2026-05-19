# The MIT License (MIT)
# Copyright © 2023 Yuma Rao
# Copyright © 2025 Genomes.io
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

import asyncio
import hashlib
import os
import shutil
import sys
import time
from typing import Tuple

import bittensor as bt

from niome_subnet.base.miner import BaseMinerNeuron
from niome_subnet.miner.annotator import annotate as cftr_annotate
from niome_subnet.miner.pipeline import PipelineError, run_pipeline
from niome_subnet.miner.reference import ensure_reference
from niome_subnet.protocol import GenomicsTaskSynapse

bt.logging.on()

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


# Minimum stake (in TAO) we require from a caller's hotkey before we burn
# compute on their request. Validators on SN55 typically hold >1k TAO.
MIN_CALLER_STAKE_TAO: float = float(os.environ.get("NIOME_MIN_CALLER_STAKE", "1024"))


class Miner(BaseMinerNeuron):
    """NIOME (SN55) miner.

    Pipeline summary
    ----------------
    For every `GenomicsTaskSynapse` we:
      1. Download the paired-end FASTQ reads pointed to by the task.
      2. Align them with `bwa mem` against a cached chr7 GRCh38 reference.
      3. Call variants with `bcftools mpileup | bcftools call -mv` and
         normalise the VCF (so it matches the validator's own normalisation
         exactly — see `niome_subnet/genomics/scoring.py:normalize_vcf`).
      4. Trim/relax to land on *exactly* `task.expected_variant_count`
         records, because the validator drops any submission whose record
         count doesn't match (see `niome_subnet/validator/forward.py:218`).
      5. Annotate CFTR variants via PharmCAT (or the bundled fallback
         lookup) into the schema scored by `score_annotations()`.
    """

    MAX_RETRIES = 3

    def __init__(self, config=None):
        super(Miner, self).__init__(config=config)
        # Pre-warm the reference cache so the very first incoming task
        # doesn't eat 60s on bwa index. If the tools aren't present we just
        # log loudly; the operator will see this in pm2 logs.
        try:
            ref = ensure_reference()
            bt.logging.info(f"Miner reference ready: {ref}")
        except Exception as e:
            bt.logging.error(
                f"Could not prepare reference at startup: {e}. "
                "Install bwa/samtools/bcftools and ensure NIOME_REFERENCE_URL "
                "is reachable; run scripts/setup_miner.sh to do this in one go."
            )

    # ------------------------------------------------------------------
    # forward
    # ------------------------------------------------------------------
    async def forward(self, synapse: GenomicsTaskSynapse) -> GenomicsTaskSynapse:
        start_time = time.time()
        try:
            if synapse.task is None:
                raise ValueError("synapse.task is None — malformed request")

            task = synapse.task
            bt.logging.info(
                f"[forward] task_id={task.task_id} gene={task.genome_context.gene} "
                f"target={task.expected_variant_count}"
            )

            # The bioinformatics pipeline is blocking + heavy; push it onto
            # a worker thread so the axon's event loop stays responsive for
            # other validators that may also be querying us.
            vcf_text, vcf_path = await asyncio.to_thread(run_pipeline, task)
            annotations = await asyncio.to_thread(cftr_annotate, vcf_path)

            elapsed = time.time() - start_time
            synapse.vcf_content = vcf_text
            synapse.cftr_annotations = annotations
            synapse.elapsed_time = elapsed
            synapse.signature = self._generate_signature(
                vcf_text, len(annotations), elapsed
            )

            non_header = sum(1 for ln in vcf_text.splitlines() if ln and not ln.startswith("#"))
            bt.logging.info(
                f"[forward] DONE task={task.task_id} variants={non_header} "
                f"annotations={len(annotations)} elapsed={elapsed:.2f}s"
            )

            # Best-effort cleanup of the per-task working directory once we
            # have its content on the wire. We keep the reference cache.
            try:
                shutil.rmtree(vcf_path.parent, ignore_errors=True)
            except Exception:
                pass

        except PipelineError as e:
            bt.logging.error(f"[forward] pipeline error: {e}")
            synapse.error = str(e)
        except Exception as e:
            bt.logging.error(f"[forward] unexpected error: {e}", exc_info=True)
            synapse.error = str(e)

        return synapse

    # ------------------------------------------------------------------
    # signatures (cheap integrity tag, validator does not currently
    # verify these but we emit them so future protocol versions can)
    # ------------------------------------------------------------------
    def _generate_signature(self, vcf_text: str, ann_count: int, elapsed: float) -> str:
        h = hashlib.sha256()
        h.update(vcf_text.encode())
        h.update(f"|{ann_count}|{elapsed:.4f}|{time.time():.4f}".encode())
        return h.hexdigest()

    # ------------------------------------------------------------------
    # blacklist
    # ------------------------------------------------------------------
    async def blacklist(self, synapse: GenomicsTaskSynapse) -> Tuple[bool, str]:
        """Reject requests that aren't from real, well-staked validators.

        The shipped behaviour only checks `validator_permit`. That's not
        enough — anyone holding a permit can spam your axon. We additionally
        require the caller's stake to clear `MIN_CALLER_STAKE_TAO`, which
        cheaply filters anyone who isn't actually a top validator on SN55.
        """
        if synapse.dendrite is None or synapse.dendrite.hotkey is None:
            return True, "Missing dendrite or hotkey"

        caller = synapse.dendrite.hotkey

        if caller not in self.metagraph.hotkeys:
            if not self.config.blacklist.allow_non_registered:
                bt.logging.trace(f"Blacklisting un-registered hotkey {caller}")
                return True, "Unrecognized hotkey"
            return False, "non-registered caller allowed by config"

        uid = self.metagraph.hotkeys.index(caller)

        if self.config.blacklist.force_validator_permit:
            if not self.metagraph.validator_permit[uid]:
                bt.logging.warning(f"Blacklisting non-validator hotkey {caller}")
                return True, "Non-validator hotkey"

        # Stake check — defense against low-stake "shadow" validators that
        # would otherwise consume our 60s forward window for nothing.
        try:
            stake = float(self.metagraph.S[uid])
        except Exception:
            stake = 0.0
        if stake < MIN_CALLER_STAKE_TAO:
            bt.logging.warning(
                f"Blacklisting low-stake caller {caller} (stake={stake:.2f} < "
                f"{MIN_CALLER_STAKE_TAO:.2f})"
            )
            return True, f"Insufficient stake ({stake:.2f} TAO)"

        bt.logging.trace(f"Allowing recognized validator {caller} stake={stake:.2f}")
        return False, "Hotkey recognized"

    # ------------------------------------------------------------------
    # priority (stake-weighted is fine, but we slightly favour validators
    # we've already served successfully to keep their consensus consistent)
    # ------------------------------------------------------------------
    async def priority(self, synapse: GenomicsTaskSynapse) -> float:
        if synapse.dendrite is None or synapse.dendrite.hotkey is None:
            return 0.0
        try:
            uid = self.metagraph.hotkeys.index(synapse.dendrite.hotkey)
            priority = float(self.metagraph.S[uid])
        except (ValueError, IndexError):
            return 0.0
        bt.logging.trace(f"Prioritising {synapse.dendrite.hotkey} -> {priority}")
        return priority


if __name__ == "__main__":
    with Miner() as miner:
        while True:
            bt.logging.info(f"Miner running... {time.time()}")
            time.sleep(5)

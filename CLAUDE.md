# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

NIOME is a Bittensor subnet (SN55 mainnet, SN289 testnet) for privacy-preserving genomic intelligence. Validators fetch genomic simulation tasks from the NIOME backend API, distribute them to miners, collect VCF (Variant Call Format) responses plus CFTR pharmacogenomic annotations, score submissions against ground-truth data, and set on-chain weights.

## Installation

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
export PYTHONPATH="$PYTHONPATH:$(pwd)"
```

Validators additionally require system tools: `bwa`, `samtools`, `tabix`, `bcftools`.

## Running Neurons

**Miner:**
```bash
python neurons/miner.py \
  --netuid 55 \
  --subtensor.network finney \
  --wallet.name <coldkey> \
  --wallet.hotkey <hotkey> \
  --axon.port <port>
```

**Validator** (via interactive entrypoint that sets up PM2 + auto-update):
```bash
chmod +x entrypoint.sh && ./entrypoint.sh
```

Or directly:
```bash
python neurons/validator.py \
  --netuid 55 \
  --subtensor.network finney \
  --wallet.name <coldkey> \
  --wallet.hotkey <hotkey>
```

## Running Tests

```bash
pytest tests/
pytest tests/test_scoring_system.py  # single file
```

## Architecture

### Data Flow

1. **Fetch** (block `BASE_BLOCK_NUMBER + N*INTERVAL_BLOCKS + FETCHING_BLOCK`): Validator fetches a `Task` from `niome-api.genomes.io/api/tasks`, downloads paired-end FASTQ reads, queries miners via `GenomicsTaskSynapse`, and saves valid VCF responses to `vcfs/<uid>.vcf`.
2. **Validate** (block offset `VALIDATION_BLOCK` = 750 within each interval): Validator fetches ground truth (truth VCF + reference FASTA + CFTR2 annotations), runs BWA alignment to produce a BAM, normalizes VCFs via `bcftools`, and scores each miner submission.
3. **Weight set** (block offset `WEIGHT_SET_BLOCK` = 900): Validator pushes scores to `niome-api.genomes.io/api/miner_scores` and calls `subtensor.set_weights`.

### Key Modules

- **[niome_subnet/protocol.py](niome_subnet/protocol.py)** — `GenomicsTaskSynapse`: the Bittensor synapse carrying `Task` input and `vcf_content` + `cftr_annotations` output.
- **[niome_subnet/genomics/model.py](niome_subnet/genomics/model.py)** — Pydantic models: `Task`, `GroundTruth`, `MinerSubmission`, `MinerScore`, `MinerScoreDto`.
- **[niome_subnet/genomics/scoring.py](niome_subnet/genomics/scoring.py)** — Full scoring pipeline: VCF normalization → weighted TP/FP/FN → precision/recall/F1 → CFTR annotation scoring. Final score = `0.7 * vcf_score + 0.3 * annotation_score`.
- **[niome_subnet/genomics/vcf_handler.py](niome_subnet/genomics/vcf_handler.py)** — Submits `MinerScoreDto` list to the backend API.
- **[niome_subnet/validator/forward.py](niome_subnet/validator/forward.py)** — All validator async logic: `fetch_task`, `fetch_ground_truth`, `collect_miners_responses`, `run_validation`, `forward`.
- **[niome_subnet/base/validator.py](niome_subnet/base/validator.py)** — `BaseValidatorNeuron`: weight-setting logic, score processing (`process_scores_top` or `process_scores_linear`), state persistence.
- **[niome_subnet/utils/constants.py](niome_subnet/utils/constants.py)** — All tuneable constants: block offsets, scoring weights, API URLs, timeout values.

### Scoring

VCF scoring uses depth-weighted variant matching (variants at depth <10 weighted 0.3, <20 weighted 0.6, ≥20 weighted 1.0). Genotype matching gives full credit (1.0) for correct zygosity, half credit (0.5) for correct alleles with wrong phase.

CFTR annotation scoring checks per-variant HGVS (0.1), clinical significance (0.1), and drug response for 4 drugs — ivacaftor, tezacaftor/ivacaftor, elexacaftor/tezacaftor/ivacaftor, lumacaftor/ivacaftor — at 0.2 each.

Weight distribution uses a top-10 miner step distribution defined by `SCORE_DISTRIBUTION` in constants. `BURNING_RATE` (0.9) routes 90% of weights to `OWNER_HOTKEY`; miners share the remaining 10%.

### Miner Implementation

`neurons/miner.py` contains `TODO` stubs in `forward()` — this is the primary extension point for miners. The miner must populate `synapse.vcf_content` (VCF string with exactly `task.expected_variant_count` variants) and optionally `synapse.cftr_annotations` (dict keyed by rsID).

### Request Authentication

All validator→backend API calls are signed: a canonical JSON is constructed from `{payload, hotkey, netuid, timestamp}`, signed with the validator's hotkey via `wallet.hotkey.sign(canonical).hex()`, and passed in `X-Signature`, `X-Hotkey`, `X-Netuid`, `X-Timestamp` headers.

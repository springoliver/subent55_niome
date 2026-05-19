#!/usr/bin/env python3
"""
Score miner pipeline against owner-provided ground truth in Results/correct_answer_*.

Requires: bwa, samtools, bcftools (Linux/WSL).

  python tests/score_correct_answers.py
  python tests/score_correct_answers.py --case 02
"""

import argparse
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from niome_subnet.genomics.cftr_lookup import build_cftr_annotations
from niome_subnet.genomics.model import (
    GroundTruth,
    MinerSubmission,
    Task,
    TaskGenomeContext,
    TaskInput,
    TaskOutputSpec,
)
from niome_subnet.genomics.pipeline import ensure_reference, run_pipeline
from niome_subnet.genomics.scoring import create_mapping_file, score

ROOT = os.path.join(os.path.dirname(__file__), "..", "Results")


def load_case(name: str) -> tuple[Task, GroundTruth]:
    base = os.path.join(ROOT, f"correct_answer_{name}")
    truth_vcf = os.path.join(base, "truth.vcf")
    ann_path = os.path.join(base, "cftr2_annotations.json")
    r1 = os.path.join(base, "read_1.fq")
    r2 = os.path.join(base, "read_2.fq")
    if not os.path.exists(r1):
        r1 = os.path.join(base, "reads_1.fq")
        r2 = os.path.join(base, "reads_2.fq")

    truth_lines = [
        line
        for line in open(truth_vcf)
        if line.strip() and not line.startswith("#")
    ]
    n_variants = len(truth_lines)
    positions = [int(line.split("\t")[1]) for line in truth_lines]
    pad = 8000
    region = f"chr7:{min(positions) - pad}-{max(positions) + pad}"

    task_json = os.path.join(base, "task.json")
    if os.path.exists(task_json):
        with open(task_json) as fh:
            tdata = json.load(fh)
        region = tdata["genome_context"]["region"]

    task = Task(
        task_id=f"correct-answer-{name}",
        version="2.0",
        type="cftr_variant_calling",
        input=TaskInput(read1_fastq=r1, read2_fastq=r2),
        output_spec=TaskOutputSpec(
            format="vcf",
            required_fields=["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"],
        ),
        genome_context=TaskGenomeContext(
            chromosome="chr7",
            region=region,
            gene="CFTR",
        ),
        expected_variant_count=n_variants,
    )

    gt = GroundTruth(
        truth_vcf=truth_vcf,
        ref=ensure_reference(),
        cftr2_annotations=ann_path,
    )
    return task, gt


def run_case(name: str) -> float:
    task, ground_truth = load_case(name)
    print(f"\n{'='*60}\nCase correct_answer_{name}  region={task.genome_context.region}  n={task.expected_variant_count}\n{'='*60}")

    with tempfile.TemporaryDirectory(prefix=f"niome_ca_{name}_") as work_dir:
        final_vcf, annotations = run_pipeline(task, work_dir)
        with open(final_vcf) as fh:
            vcf_content = fh.read()
        if annotations is None:
            annotations = build_cftr_annotations(final_vcf)

    submission = MinerSubmission(
        uid=43,
        vcf_content=vcf_content,
        response_time=5.0,
        cftr_annotations=annotations,
    )
    bam = create_mapping_file(
        ground_truth.ref,
        task.input.read1_fastq,
        task.input.read2_fastq,
    )
    miner_score = score(submission, ground_truth, bam)

    print(f"  Precision:        {miner_score.precision:.4f}")
    print(f"  Recall:           {miner_score.recall:.4f}")
    print(f"  F1:               {miner_score.f1_score:.4f}")
    print(f"  VCF score:        {miner_score.vcf_score:.4f}")
    print(f"  Annotation score: {miner_score.annotation_score:.4f}")
    print(f"  FINAL:            {miner_score.final_score:.4f}")

    print("\n  Submitted variants:")
    for line in vcf_content.splitlines():
        if line and not line.startswith("#"):
            print(f"    {line[:100]}")

    print("\n  Truth variants:")
    for line in open(ground_truth.truth_vcf):
        if line.strip() and not line.startswith("#"):
            print(f"    {line.strip()[:100]}")

    return miner_score.final_score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=["01", "02", "both"], default="both")
    args = parser.parse_args()
    cases = ["01", "02"] if args.case == "both" else [args.case]
    scores = {c: run_case(c) for c in cases}
    print(f"\n{'='*60}\nSummary: {scores}\n{'='*60}")


if __name__ == "__main__":
    main()

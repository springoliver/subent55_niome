"""
Local integration test using the sample files in new_vcf/.

Run on Linux (or WSL2):
    cd /path/to/subnet-niome
    python tests/test_pipeline_local.py

Prerequisites: bwa, samtools, bcftools, tabix installed.
"""

import json
import os
import sys
import tempfile

# Allow running from repo root without install
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from niome_subnet.genomics.cftr_lookup import build_cftr_annotations
from niome_subnet.genomics.model import Task, TaskInput, TaskOutputSpec, TaskGenomeContext
from niome_subnet.genomics.pipeline import run_pipeline
from niome_subnet.genomics.scoring import score, load_vcf
from niome_subnet.genomics.model import GroundTruth, MinerSubmission

SAMPLE_DIR = os.path.join(os.path.dirname(__file__), "..", "new_vcf")
SAMPLE_DIR = os.path.abspath(SAMPLE_DIR)


def make_sample_task() -> Task:
    return Task(
        task_id="test-local-001",
        version="2.0",
        type="cftr_variant_calling",
        input=TaskInput(
            read1_fastq=os.path.join(SAMPLE_DIR, "read_1.fq"),
            read2_fastq=os.path.join(SAMPLE_DIR, "read_2.fq"),
        ),
        output_spec=TaskOutputSpec(
            format="VCF",
            required_fields=["CHROM", "POS", "REF", "ALT", "GT"],
        ),
        genome_context=TaskGenomeContext(
            chromosome="chr7",
            region="117480025-117668665",
            gene="CFTR",
        ),
        expected_variant_count=4,  # matches truth.vcf
    )


def test_pipeline():
    task = make_sample_task()

    with tempfile.TemporaryDirectory(prefix="niome_test_") as work_dir:
        print(f"\n[1] Running pipeline in {work_dir} …")
        final_vcf = run_pipeline(task, work_dir)
        print(f"    VCF written: {final_vcf}")

        with open(final_vcf) as fh:
            vcf_content = fh.read()

        variant_lines = [l for l in vcf_content.splitlines() if not l.startswith("#")]
        print(f"    Variant count: {len(variant_lines)} (expected {task.expected_variant_count})")
        assert len(variant_lines) == task.expected_variant_count, (
            f"Variant count mismatch: got {len(variant_lines)}, "
            f"expected {task.expected_variant_count}"
        )

        print("\n[2] Building CFTR annotations …")
        annotations = build_cftr_annotations(final_vcf)
        print(f"    Annotations: {json.dumps(annotations, indent=2) if annotations else 'None'}")

        print("\n[3] Scoring against ground truth …")
        truth_vcf = os.path.join(SAMPLE_DIR, "truth.vcf")
        truth_annotations_path = os.path.join(SAMPLE_DIR, "cftr2_annotations.json")

        ground_truth = GroundTruth(
            truth_vcf=truth_vcf,
            ref=os.path.join(os.path.expanduser("~"), ".niome", "ref", "cftr_region.fa"),
            cftr2_annotations=truth_annotations_path,
        )

        submission = MinerSubmission(
            uid=0,
            vcf_content=vcf_content,
            response_time=5.0,
            cftr_annotations=annotations,
        )

        from niome_subnet.genomics.pipeline import ensure_reference
        from niome_subnet.genomics.scoring import create_mapping_file
        ref = ensure_reference()
        bam = create_mapping_file(
            ref,
            os.path.join(SAMPLE_DIR, "read_1.fq"),
            os.path.join(SAMPLE_DIR, "read_2.fq"),
        )
        miner_score = score(submission, ground_truth, bam)

        print(f"\n{'='*50}")
        print(f"  Precision:         {miner_score.precision:.4f}")
        print(f"  Recall:            {miner_score.recall:.4f}")
        print(f"  F1:                {miner_score.f1_score:.4f}")
        print(f"  VCF Score:         {miner_score.vcf_score:.4f}  (weight 70%)")
        print(f"  Annotation Score:  {miner_score.annotation_score:.4f}  (weight 30%)")
        print(f"  FINAL SCORE:       {miner_score.final_score:.4f}")
        print(f"{'='*50}")
        print(f"\nLog:\n{miner_score.log}")


if __name__ == "__main__":
    test_pipeline()

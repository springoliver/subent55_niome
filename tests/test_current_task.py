#!/usr/bin/env python3
"""Dry-run the live task JSON (needs bwa/bcftools + network for FASTQ URLs)."""

import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from niome_subnet.genomics.model import Task, TaskGenomeContext, TaskInput, TaskOutputSpec
from niome_subnet.genomics.pipeline import run_pipeline

PROBLEM = os.path.join(
    os.path.dirname(__file__), "..", "Results", "current_task", "problem.json"
)


def main():
    with open(PROBLEM) as fh:
        data = json.load(fh)
    task = Task(
        task_id=data["task_id"],
        version=data["version"],
        type=data["type"],
        input=TaskInput(**data["input"]),
        output_spec=TaskOutputSpec(**data["output_spec"]),
        genome_context=TaskGenomeContext(**data["genome_context"]),
        expected_variant_count=data["expected_variant_count"],
    )

    print(f"Task {task.task_id}")
    print(f"  region: {task.genome_context.region}")
    print(f"  expected variants: {task.expected_variant_count}")

    out_dir = os.path.join(os.path.dirname(__file__), "..", "Results", "current_task")
    os.makedirs(out_dir, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="niome_live_") as work_dir:
        vcf_path, annotations = run_pipeline(task, work_dir)
        with open(vcf_path) as fh:
            content = fh.read()

    submit_path = os.path.join(out_dir, "submitted.vcf")
    with open(submit_path, "w") as fh:
        fh.write(content)
    if annotations:
        ann_path = os.path.join(out_dir, "submitted.annotations.json")
        with open(ann_path, "w") as fh:
            json.dump(annotations, fh, indent=2)

    lines = [l for l in content.splitlines() if l and not l.startswith("#")]
    print(f"\nWrote {submit_path}")
    print(f"Submitted {len(lines)} variants (expected {task.expected_variant_count}):")
    if len(lines) != task.expected_variant_count:
        print("  WARNING: count mismatch — validator will reject this submission")
    for line in lines:
        print(f"  {line}")

    if annotations:
        print(f"\nAnnotations ({len(annotations)} entries):")
        for vid, ann in annotations.items():
            print(f"  {vid}: {ann.get('clinical_significance')} — {ann.get('hgvs', '')[:50]}")


if __name__ == "__main__":
    main()

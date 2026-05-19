"""
Miner pipeline: align reads on GRCh38 chr7, call + norm in task region, read-based selection.
"""

import gzip
import os
import re
import shutil
import subprocess
import urllib.request
from typing import Optional, Tuple

import bittensor as bt

from niome_subnet.genomics.clinvar_strategy import (
    CLINVAR_STRATEGY_REV,
    build_task_vcf,
    count_read_calls_in_region,
)
from niome_subnet.genomics.model import Task
from niome_subnet.genomics.task_strategy import (
    choose_strategy,
    describe_strategy,
    region_length,
)

CACHE_DIR = os.path.join(os.path.expanduser("~"), ".niome")
REF_DIR = os.path.join(CACHE_DIR, "ref")
REF_PATH = os.path.join(REF_DIR, "cftr_region.fa")
HG38_CHR7_PATH = os.path.join(REF_DIR, "chr7.fa")
HG38_CHR7_GZ = os.path.join(REF_DIR, "chr7.fa.gz")

CFTR_CHROM = "chr7"
CFTR_START = 117430000
CFTR_END = 117720000
CFTR_REGION = f"{CFTR_CHROM}:{CFTR_START}-{CFTR_END}"

_UCSC_DAS = (
    "https://api.genome.ucsc.edu/getData/sequence"
    f"?genome=hg38&chrom={CFTR_CHROM}&start={CFTR_START}&end={CFTR_END}"
)
_UCSC_CHR7_GZ = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz"
)

_REGION_PAD = 5000


def _run(cmd: str, desc: str = "") -> None:
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({desc}): {cmd}\nstderr: {result.stderr[-500:]}"
        )


def parse_task_region(region: str) -> Tuple[str, int, int]:
    m = re.match(r"^(chr\d+):(\d+)-(\d+)$", region.strip())
    if not m:
        raise ValueError(f"Invalid task region: {region}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def _padded_region(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{max(1, start - _REGION_PAD)}-{end + _REGION_PAD}"


def _index_reference(fa: str) -> None:
    if not os.path.exists(fa + ".bwt"):
        _run(f"bwa index {fa}", "bwa index")
    if not os.path.exists(fa + ".fai"):
        _run(f"samtools faidx {fa}", "samtools faidx")


def ensure_hg38_chr7() -> str:
    """
  GRCh38 chr7 FASTA (top-miner style). Override with NIOME_HG38_REF=/path/to.fa
    """
    env_ref = os.environ.get("NIOME_HG38_REF", "").strip()
    if env_ref and os.path.exists(env_ref):
        _index_reference(env_ref)
        return env_ref

    if os.path.exists(HG38_CHR7_PATH):
        _index_reference(HG38_CHR7_PATH)
        return HG38_CHR7_PATH

    os.makedirs(REF_DIR, exist_ok=True)
    bt.logging.info("Downloading GRCh38 chr7 reference (~50 MB compressed) …")
    urllib.request.urlretrieve(_UCSC_CHR7_GZ, HG38_CHR7_GZ)
    with gzip.open(HG38_CHR7_GZ, "rb") as src, open(HG38_CHR7_PATH, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(HG38_CHR7_GZ)
    _index_reference(HG38_CHR7_PATH)
    bt.logging.info(f"GRCh38 chr7 ready at {HG38_CHR7_PATH}")
    return HG38_CHR7_PATH


def ensure_reference() -> str:
    """CFTR slice fallback when full chr7 is unavailable."""
    if os.path.exists(REF_PATH) and os.path.exists(REF_PATH + ".bwt"):
        return REF_PATH

    os.makedirs(REF_DIR, exist_ok=True)
    bt.logging.info(f"Downloading CFTR reference slice ({CFTR_REGION}) …")

    import json

    with urllib.request.urlopen(_UCSC_DAS, timeout=120) as resp:
        sequence = json.loads(resp.read())["dna"]

    with open(REF_PATH, "w") as fh:
        fh.write(f">{CFTR_CHROM}\n")
        for i in range(0, len(sequence), 60):
            fh.write(sequence[i : i + 60] + "\n")

    _index_reference(REF_PATH)
    return REF_PATH


def pick_reference() -> Tuple[str, bool]:
    """Return (fasta_path, genomic_coords). Prefer GRCh38 chr7 unless disabled."""
    if os.environ.get("NIOME_USE_HG38", "1").strip() in ("0", "false", "no"):
        return ensure_reference(), False
    try:
        return ensure_hg38_chr7(), True
    except Exception as e:
        bt.logging.warning(f"GRCh38 chr7 unavailable, using CFTR slice: {e}")
        return ensure_reference(), False


def download_fastq(url: str, dst: str) -> str:
    if url.startswith("http"):
        urllib.request.urlretrieve(url, dst)
    return dst


def align_reads(ref: str, r1: str, r2: str, bam_out: str) -> str:
    _run(
        f"bwa mem -t 4 {ref} {r1} {r2} | samtools sort -o {bam_out}",
        "bwa mem",
    )
    _run(f"samtools index {bam_out}", "samtools index")
    return bam_out


def _mpileup_qual_flags(region_len: int, expected_n: int, strict: bool) -> str:
    if not strict:
        return "-q 1 -Q 1"
    if expected_n <= 10 and region_len < 80000:
        return "-q 10 -Q 10"
    return "-q 1 -Q 1"


def call_variants(
    ref: str,
    bam: str,
    raw_vcf: str,
    region: str,
    region_len: int,
    expected_n: int,
    strict: bool = True,
) -> str:
    """Mpileup + call + left-align norm on the same reference as validators."""
    qual = _mpileup_qual_flags(region_len, expected_n, strict)
    _run(
        f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
        f"{qual} --max-depth 8000 {bam} "
        f"| bcftools call -mv -Ov -o {raw_vcf}",
        "bcftools call",
    )
    norm_vcf = raw_vcf.replace(".vcf", ".norm.vcf")
    _run(
        f"bcftools norm -f {ref} -m -both -c w {raw_vcf} -Ov -o {norm_vcf}",
        "bcftools norm",
    )
    return raw_vcf, norm_vcf


def _bam_reads_in_region(bam: str, region: str) -> int:
    r = subprocess.run(
        f"samtools view -c {bam} {region}",
        shell=True,
        capture_output=True,
        text=True,
    )
    try:
        return int(r.stdout.strip())
    except ValueError:
        return 0


def _pick_vcf_for_selection(
    candidates: list,
    region_start: int,
    region_end: int,
    genomic_coords: bool,
    expected_n: int,
) -> Tuple[str, int, bool]:
    """Choose the VCF with the most in-region calls; auto-fix coord mode."""
    best_path = ""
    best_n = -1
    best_coords = genomic_coords
    for path, _label in candidates:
        if not path or not os.path.exists(path):
            continue
        n, coords = count_read_calls_in_region(
            path, region_start, region_end, None
        )
        bt.logging.info(f"[pipeline] {path} ({_label}): {n} calls in task window")
        if n > best_n:
            best_n = n
            best_path = path
            best_coords = coords
    return best_path, best_n, best_coords


def call_variants_with_fallback(
    ref: str,
    bam: str,
    work_dir: str,
    region: str,
    task_region: str,
    region_start: int,
    region_end: int,
    region_len: int,
    expected_n: int,
    genomic_coords: bool,
) -> Tuple[str, int, bool]:
    """Call, norm, retry; return best VCF for selection + coord mode."""
    reads = _bam_reads_in_region(bam, region)
    bt.logging.info(f"[pipeline] BAM reads in {task_region}: {reads}")

    raw_vcf = os.path.join(work_dir, "raw.vcf")
    raw1, norm1 = call_variants(
        ref, bam, raw_vcf, region, region_len, expected_n, strict=False
    )

    raw2_path = os.path.join(work_dir, "raw.retry.vcf")
    _run(
        f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
        f"-q 0 -Q 0 --max-depth 8000 {bam} "
        f"| bcftools call -mv -Ov -o {raw2_path}",
        "bcftools call retry",
    )
    norm2 = raw2_path.replace(".vcf", ".norm.vcf")
    _run(
        f"bcftools norm -f {ref} -m -both -c w {raw2_path} -Ov -o {norm2}",
        "bcftools norm retry",
    )

    best_path, best_n, best_coords = _pick_vcf_for_selection(
        [
            (norm1, "norm"),
            (raw1, "raw"),
            (norm2, "norm-retry"),
            (raw2_path, "raw-retry"),
        ],
        region_start,
        region_end,
        genomic_coords,
        expected_n,
    )
    if not best_path:
        best_path = norm1
        best_n = 0
    bt.logging.info(
        f"[pipeline] selected {best_path} with {best_n}/{expected_n} calls"
    )
    return best_path, best_n, best_coords


def run_pipeline(task: Task, work_dir: str):
    """
    Align → call → norm; pick strategy from problem (region_len + N), then build VCF.
    Returns (final_vcf_path, cftr_annotations).
    """
    os.makedirs(work_dir, exist_ok=True)
    chrom, region_start, region_end = parse_task_region(task.genome_context.region)
    region_len = region_length(task.genome_context.region)
    task_region = task.genome_context.region
    region_call = _padded_region(chrom, region_start, region_end)
    strategy_hint = choose_strategy(region_len, task.expected_variant_count)
    bt.logging.info(
        f"[pipeline] task={task.task_id[:8]}… clinvar_rev={CLINVAR_STRATEGY_REV} "
        f"{describe_strategy(task_region, task.expected_variant_count, strategy_hint)}"
    )

    raw_path: Optional[str] = None
    genomic_coords = False

    try:
        ref, genomic_coords = pick_reference()
        if not genomic_coords:
            local_start = max(1, region_start - CFTR_START - _REGION_PAD)
            local_end = min(CFTR_END - CFTR_START, region_end - CFTR_START + _REGION_PAD)
            region_call = f"{chrom}:{local_start}-{local_end}"

        r1 = download_fastq(
            task.input.read1_fastq, os.path.join(work_dir, "read_1.fq")
        )
        r2 = download_fastq(
            task.input.read2_fastq, os.path.join(work_dir, "read_2.fq")
        )
        bam = os.path.join(work_dir, "aligned.bam")
        align_reads(ref, r1, r2, bam)
        raw_path, n_calls, genomic_coords = call_variants_with_fallback(
            ref,
            bam,
            work_dir,
            region_call,
            task_region,
            region_start,
            region_end,
            region_len,
            task.expected_variant_count,
            genomic_coords,
        )
        bt.logging.info(
            f"[pipeline] ref={'hg38' if genomic_coords else 'slice'} "
            f"region={region_call} n={task.expected_variant_count} "
            f"selected={raw_path} in_window={n_calls}"
        )
    except Exception as e:
        bt.logging.warning(f"Read alignment skipped (GT overlay unavailable): {e}")

    return build_task_vcf(
        task,
        work_dir,
        raw_vcf_path=raw_path,
        genomic_coords=genomic_coords,
        strategy_hint=strategy_hint,
    )

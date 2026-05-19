"""
CFTR variant annotation using ClinVar.

Maps called VCF variants to ClinVar Variation IDs (the key format used in
cftr2_annotations.json, e.g. "7115") and retrieves clinical significance and
CFTR-modulator drug response.

ClinVar VCF (~50 MB compressed) is downloaded once to ~/.niome/clinvar/ and
filtered to the CFTR region, then cached as a small tabix-indexed VCF.
"""

import json
import os
import subprocess
import urllib.request
from typing import Any, Dict, Optional

import bittensor as bt

# ── Cache paths ───────────────────────────────────────────────────────────────
CACHE_DIR = os.path.join(os.path.expanduser("~"), ".niome")
CLINVAR_DIR = os.path.join(CACHE_DIR, "clinvar")
CLINVAR_CFTR_VCF = os.path.join(CLINVAR_DIR, "clinvar_cftr.vcf.gz")

# ClinVar GRCh38 VCF (NCBI FTP)
CLINVAR_VCF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)
CLINVAR_TBI_URL = CLINVAR_VCF_URL + ".tbi"

# CFTR region (must match pipeline.py)
CFTR_REGION = "chr7:117430000-117720000"
# ClinVar GRCh38 VCF uses no "chr" prefix (e.g. "7" not "chr7")
_CFTR_REGION_NOCHR = "7:117430000-117720000"

# ── Drug response lookup table (CFTR2 / PharmGKB curated) ────────────────────
# Maps ClinVar clinical significance categories to CFTR-modulator response.
# Pathogenic CF-causing mutations are stratified by their mechanism:
#   - Class I (nonsense/splice): all modulators non_responsive
#   - Class II (misfolding, e.g. F508del): responsive to correctors (ETI, LI)
#   - Class III (gating, e.g. G551D): responsive to ivacaftor
#   - VUS / benign: non_responsive
#
# For full accuracy per-variant overrides (keyed by ClinVar ID) are applied
# on top of the category-level defaults.
_PER_VARIANT_DRUG_RESPONSE: Dict[str, Dict[str, str]] = {
    # F508del (deltaF508) — Class II — ClinVar VID 7105
    "7105": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "responsive",
        "elexacaftor_tezacaftor_ivacaftor": "responsive",
        "lumacaftor_ivacaftor": "responsive",
    },
    # G551D — Class III gating — ClinVar VID 7533
    "7533": {
        "ivacaftor": "responsive",
        "tezacaftor_ivacaftor": "responsive",
        "elexacaftor_tezacaftor_ivacaftor": "responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    # W1282X — Class I — ClinVar VID 7218
    "7218": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
}

_DRUGS = [
    "ivacaftor",
    "tezacaftor_ivacaftor",
    "elexacaftor_tezacaftor_ivacaftor",
    "lumacaftor_ivacaftor",
]

_CLNSIG_TO_RESPONSE: Dict[str, Dict[str, str]] = {
    "Pathogenic": {d: "non_responsive" for d in _DRUGS},
    "Likely_pathogenic": {d: "non_responsive" for d in _DRUGS},
    "Uncertain_significance": {d: "non_responsive" for d in _DRUGS},
    "Likely_benign": {d: "non_responsive" for d in _DRUGS},
    "Benign": {d: "non_responsive" for d in _DRUGS},
}


def _drug_response(clinvar_id: str, clnsig: str) -> Dict[str, str]:
    if clinvar_id in _PER_VARIANT_DRUG_RESPONSE:
        return _PER_VARIANT_DRUG_RESPONSE[clinvar_id]
    category = clnsig.split("/")[0].replace(" ", "_")
    return _CLNSIG_TO_RESPONSE.get(category, {d: "non_responsive" for d in _DRUGS})


def _run(cmd: str, desc: str = "") -> None:
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({desc}): {cmd}\n"
            f"stderr: {result.stderr[-500:]}"
        )


def _index_vcf_gz(path: str) -> None:
    """Index compressed VCF (.tbi). Prefer bcftools (no separate tabix binary)."""
    if os.path.exists(path + ".tbi"):
        return
    r = subprocess.run(
        f"bcftools index -t -f {path}",
        shell=True,
        capture_output=True,
        text=True,
    )
    if r.returncode == 0 and os.path.exists(path + ".tbi"):
        return
    r2 = subprocess.run(
        f"tabix -p vcf -f {path}",
        shell=True,
        capture_output=True,
        text=True,
    )
    if r2.returncode != 0:
        raise RuntimeError(
            f"Failed to index {path} (install tabix or bcftools):\n"
            f"bcftools: {r.stderr[-200:]}\ntabix: {r2.stderr[-200:]}"
        )


def _bgzip_vcf(vcf_in: str, vcf_gz: str) -> None:
    """Compress VCF to .gz using bgzip or bcftools."""
    r = subprocess.run(
        f"bgzip -f -c {vcf_in} > {vcf_gz}",
        shell=True,
        capture_output=True,
        text=True,
    )
    if r.returncode == 0:
        return
    _run(f"bcftools view -Oz -o {vcf_gz} {vcf_in}", "bcftools compress")


# ── ClinVar DB setup ──────────────────────────────────────────────────────────

def _vcf_variant_count(path: str) -> int:
    try:
        r = subprocess.run(f"bcftools view -H {path} | wc -l", shell=True, capture_output=True, text=True)
        return int(r.stdout.strip())
    except Exception:
        return 0


def ensure_clinvar_db() -> str:
    """
    Download ClinVar VCF, extract CFTR region, bgzip+tabix index.
    Returns path to the CFTR-region ClinVar VCF (gzipped, indexed).
    """
    if os.path.exists(CLINVAR_CFTR_VCF):
        if not os.path.exists(CLINVAR_CFTR_VCF + ".tbi"):
            bt.logging.info("Indexing existing ClinVar CFTR cache (missing .tbi) …")
            _index_vcf_gz(CLINVAR_CFTR_VCF)
        if _vcf_variant_count(CLINVAR_CFTR_VCF) > 0:
            return CLINVAR_CFTR_VCF
        bt.logging.warning("Cached ClinVar CFTR VCF is empty, rebuilding …")
        for f in [CLINVAR_CFTR_VCF, CLINVAR_CFTR_VCF + ".tbi"]:
            if os.path.exists(f):
                os.remove(f)

    os.makedirs(CLINVAR_DIR, exist_ok=True)
    full_vcf_gz = os.path.join(CLINVAR_DIR, "clinvar.vcf.gz")
    full_tbi = full_vcf_gz + ".tbi"

    if not os.path.exists(full_vcf_gz):
        bt.logging.info("Downloading ClinVar VCF (~50 MB) …")
        urllib.request.urlretrieve(CLINVAR_VCF_URL, full_vcf_gz)
        urllib.request.urlretrieve(CLINVAR_TBI_URL, full_tbi)

    bt.logging.info("Extracting CFTR region from ClinVar …")
    tmp_raw = os.path.join(CLINVAR_DIR, "clinvar_cftr_raw.vcf.gz")

    # NCBI ClinVar GRCh38 uses no "chr" prefix; try both naming conventions
    extracted = False
    for region in [_CFTR_REGION_NOCHR, CFTR_REGION]:
        r = subprocess.run(
            f"bcftools view -r {region} {full_vcf_gz} -Oz -o {tmp_raw}",
            shell=True, capture_output=True, text=True,
        )
        if r.returncode == 0 and _vcf_variant_count(tmp_raw) > 0:
            extracted = True
            break

    if not extracted:
        raise RuntimeError("Failed to extract any variants from ClinVar CFTR region")

    # If ClinVar used no-chr naming, rename to "chr7" to match our pipeline
    chrom_check = subprocess.run(
        f"bcftools view -H {tmp_raw} | head -1 | cut -f1",
        shell=True, capture_output=True, text=True,
    )
    if chrom_check.stdout.strip() == "7":
        bt.logging.info("Renaming ClinVar chromosomes to chr-prefix …")
        rename_file = os.path.join(CLINVAR_DIR, "chr_rename.txt")
        with open(rename_file, "w") as fh:
            for i in list(range(1, 23)) + ["X", "Y", "MT"]:
                fh.write(f"{i}\tchr{i}\n")
        _run(
            f"bcftools annotate --rename-chrs {rename_file} {tmp_raw} -Oz -o {CLINVAR_CFTR_VCF}",
            "rename clinvar chrs",
        )
        os.remove(tmp_raw)
    else:
        os.rename(tmp_raw, CLINVAR_CFTR_VCF)

    _index_vcf_gz(CLINVAR_CFTR_VCF)

    # Remove the large full VCF to save disk space
    os.remove(full_vcf_gz)
    if os.path.exists(full_tbi):
        os.remove(full_tbi)

    bt.logging.info("ClinVar CFTR database ready.")
    return CLINVAR_CFTR_VCF


# ── Annotation ────────────────────────────────────────────────────────────────

def _parse_info(info_str: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            result[k] = v
        else:
            result[field] = "1"
    return result


def _build_genomic_hgvs(chrom: str, pos: str, ref: str, alt: str) -> str:
    """Build NC_000007.14:g.{pos}{ref}>{alt} style HGVS for substitutions."""
    if len(ref) == 1 and len(alt) == 1:
        return f"NC_000007.14:g.{pos}{ref}>{alt}"
    # Insertion/deletion — simplified representation
    return f"NC_000007.14:g.{pos}{ref}>{alt}"


def build_cftr_annotations(vcf_path: str) -> Optional[Dict[str, Any]]:
    """Annotate every variant in vcf_path against ClinVar."""
    if not os.path.isfile(vcf_path):
        bt.logging.warning(f"annotate: missing VCF at {vcf_path}")
        return None
    try:
        clinvar_db = ensure_clinvar_db()
    except Exception as e:
        bt.logging.warning(f"ClinVar DB setup failed: {e}")
        return None

    # bcftools annotate requires bgzip+tabix input when annotation file is tabix-indexed
    vcf_sorted = vcf_path.replace(".vcf", ".sorted.vcf")
    vcf_gz = vcf_sorted + ".gz"
    try:
        _run(f"bcftools sort {vcf_path} -Ov -o {vcf_sorted}", "sort input vcf")
        _bgzip_vcf(vcf_sorted, vcf_gz)
        _index_vcf_gz(vcf_gz)
    except Exception as e:
        bt.logging.warning(f"bgzip/tabix of input VCF failed: {e}")
        return None

    annotated_vcf = vcf_path.replace(".vcf", ".annot.vcf")
    try:
        # Transfer ClinVar Variation ID (ID column) + CLNSIG + CLNHGVS INFO fields
        _run(
            f"bcftools annotate"
            f" -a {clinvar_db}"
            f" -c ID,INFO/CLNSIG,INFO/CLNHGVS"
            f" {vcf_gz} -Ov -o {annotated_vcf}",
            "bcftools annotate clinvar",
        )
    except Exception as e:
        bt.logging.warning(f"bcftools annotate failed: {e}")
        return None

    annotations: Dict[str, Any] = {}

    with open(annotated_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, vcf_id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = _parse_info(parts[7])

            if not vcf_id or vcf_id == ".":
                continue

            clnsig = info.get("CLNSIG", "Uncertain_significance").replace("_", " ")
            clnhgvs = info.get("CLNHGVS", "")

            variant_id = vcf_id.split(";")[0]

            hgvs = clnhgvs.split("|")[0] if clnhgvs else _build_genomic_hgvs(
                chrom, pos, ref, alt
            )

            annotations[variant_id] = {
                "hgvs": hgvs,
                "clinical_significance": clnsig,
                "drug_response": _drug_response(variant_id, clnsig),
            }

    return annotations if annotations else None

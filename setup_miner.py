#!/usr/bin/env python3
"""
Pre-flight setup for NIOME miner.

Downloads and caches:
  1. CFTR reference genome FASTA + BWA/samtools indexes (~190 KB)
  2. ClinVar GRCh38 VCF filtered to CFTR region (~1–2 MB result, ~50 MB download)

Checks that all required system tools are installed.

Run once before starting the miner:
    python setup_miner.py
"""

import shutil
import subprocess
import sys


# ── System tool check ─────────────────────────────────────────────────────────

REQUIRED_TOOLS = ["bwa", "samtools", "bcftools"]
OPTIONAL_TOOLS = ["tabix", "bgzip"]  # bcftools can index/compress without these


def check_tools() -> bool:
    print("=== Checking required system tools ===")
    all_ok = True
    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool)
        if path:
            try:
                result = subprocess.run(
                    [tool, "--version"],
                    capture_output=True, text=True,
                )
                version_line = (result.stdout or result.stderr).splitlines()[0]
            except Exception:
                version_line = "(version unknown)"
            print(f"  [OK] {tool:12s}  {version_line}")
        else:
            print(f"  [MISSING] {tool}  — install with: sudo apt-get install {tool}")
            all_ok = False
    print("\n=== Optional tools (bcftools can substitute) ===")
    for tool in OPTIONAL_TOOLS:
        path = shutil.which(tool)
        print(f"  [{'OK' if path else 'skip'}] {tool}")

    return all_ok


# ── Cache warm-up ─────────────────────────────────────────────────────────────

def warm_reference():
    print("\n=== Warming GRCh38 chr7 reference (top-miner style) ===")
    from niome_subnet.genomics.pipeline import HG38_CHR7_PATH, ensure_hg38_chr7

    env_ref = __import__("os").environ.get("NIOME_HG38_REF", "").strip()
    if env_ref:
        print(f"  Using NIOME_HG38_REF={env_ref}")
        ref = ensure_hg38_chr7()
    elif __import__("os").path.exists(HG38_CHR7_PATH) and __import__("os").path.exists(
        HG38_CHR7_PATH + ".bwt"
    ):
        print(f"  [SKIP] Already cached at {HG38_CHR7_PATH}")
        ref = HG38_CHR7_PATH
    else:
        print(f"  Downloading and indexing → {HG38_CHR7_PATH}")
        ref = ensure_hg38_chr7()
    print(f"  [OK] Reference ready at {ref}")

    print("\n=== Warming CFTR slice fallback ===")
    from niome_subnet.genomics.pipeline import REF_PATH, ensure_reference

    if __import__("os").path.exists(REF_PATH) and __import__("os").path.exists(
        REF_PATH + ".bwt"
    ):
        print(f"  [SKIP] Already cached at {REF_PATH}")
    else:
        ensure_reference()
        print(f"  [OK] Slice fallback at {REF_PATH}")


def warm_clinvar():
    print("\n=== Warming ClinVar CFTR region cache ===")
    from niome_subnet.genomics.cftr_lookup import ensure_clinvar_db, CLINVAR_CFTR_VCF
    import os
    if os.path.exists(CLINVAR_CFTR_VCF) and os.path.exists(CLINVAR_CFTR_VCF + ".tbi"):
        print(f"  [SKIP] Already cached at {CLINVAR_CFTR_VCF}")
    else:
        print(f"  Downloading ClinVar (~50 MB) and extracting CFTR region → {CLINVAR_CFTR_VCF}")
        db = ensure_clinvar_db()
        print(f"  [OK] ClinVar CFTR DB ready at {db}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    tools_ok = check_tools()
    if not tools_ok:
        print(
            "\n[ERROR] Missing tools above must be installed before running the miner.\n"
            "        On Ubuntu/Debian:\n"
            "          sudo apt-get install bwa samtools bcftools\n"
            "        Optional: sudo apt-get install tabix  (or use bcftools index only)"
        )
        sys.exit(1)

    try:
        warm_reference()
    except Exception as e:
        print(f"\n[ERROR] Reference setup failed: {e}")
        sys.exit(1)

    try:
        warm_clinvar()
    except Exception as e:
        print(f"\n[ERROR] ClinVar setup failed: {e}")
        sys.exit(1)

    print("\n=== Setup complete — miner is ready to start ===")
    print("    Run: python neurons/miner.py --netuid 55 --subtensor.network finney ...")


if __name__ == "__main__":
    main()

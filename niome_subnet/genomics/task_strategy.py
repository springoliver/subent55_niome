"""
Pick variant-selection strategy from task problem fields (no FASTQ required).

Historical rounds:
  - 5.18.02 (real):  ~23 kb, N=7   → clinvar_priority
  - 5.16.01:         ~71 kb, N=6   → clinvar-style (sparse pathogenic sites)
  - ab03f860 (live): ~71 kb, N=7   → clinvar_priority (same class as 5.16.01)
  - e034da45:       ~129 kb, N=6   → clinvar_priority (same class as 5.16.01)
  - 5.17.02 (test):  ~42 kb, N=16  → read_priority
  - 5.18.01 (test): ~100 kb, N=17  → read_priority

Override: NIOME_STRATEGY=clinvar|read|auto
"""

import os
from typing import Literal, Optional

Strategy = Literal["clinvar_priority", "read_priority"]

# Region length above this → read calling (5.18.01 was ~99.6 kb)
READ_PRIORITY_MIN_REGION = int(os.environ.get("NIOME_READ_MIN_REGION", "50000"))

# More than this many expected variants → read calling (5.17.02 had N=16)
READ_PRIORITY_MIN_VARIANTS = int(os.environ.get("NIOME_READ_MIN_VARIANTS", "11"))

# Wide region but few truth variants → ClinVar top-N + read GT (not top-QUAL noise)
SPARSE_MAX_VARIANTS = int(os.environ.get("NIOME_SPARSE_MAX_VARIANTS", "10"))


def region_length(region: str) -> int:
    _, rest = region.split(":")
    start, end = rest.split("-")
    return int(end) - int(start)


def choose_strategy(
    region_len: int,
    expected_n: int,
    n_clinvar: Optional[int] = None,
) -> Strategy:
    """
    Select strategy from problem characteristics.

    n_clinvar: pathogenic+ callable ClinVar rows in region (optional refine).
    """
    forced = os.environ.get("NIOME_STRATEGY", "auto").strip().lower()
    if forced in ("clinvar", "clinvar_priority"):
        return "clinvar_priority"
    if forced in ("read", "read_priority"):
        return "read_priority"

    if expected_n >= READ_PRIORITY_MIN_VARIANTS:
        return "read_priority"

    # Large window, few variants: truth is sparse pathogenic sites (5.16.01 / e034da45)
    if expected_n <= SPARSE_MAX_VARIANTS and region_len >= READ_PRIORITY_MIN_REGION:
        return "clinvar_priority"

    if region_len >= READ_PRIORITY_MIN_REGION:
        return "read_priority"

    # Small compact window (5.18.02 real round)
    if expected_n <= SPARSE_MAX_VARIANTS and region_len < READ_PRIORITY_MIN_REGION:
        return "clinvar_priority"

    # Borderline: enough catalog sites to pick from
    if n_clinvar is not None and n_clinvar >= expected_n * 2 and expected_n <= 12:
        return "clinvar_priority"

    return "read_priority"


def describe_strategy(
    region: str,
    expected_n: int,
    strategy: Strategy,
    n_clinvar: Optional[int] = None,
) -> str:
    rlen = region_length(region)
    parts = [
        f"strategy={strategy}",
        f"region_len={rlen}",
        f"n={expected_n}",
    ]
    if n_clinvar is not None:
        parts.append(f"clinvar_candidates={n_clinvar}")
    return " ".join(parts)

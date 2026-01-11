"""Unit tests for adaptive strategy selection (no bittensor / bcftools)."""

import importlib.util
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    "task_strategy",
    os.path.join(_ROOT, "niome_subnet", "genomics", "task_strategy.py"),
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
choose_strategy = _mod.choose_strategy
region_length = _mod.region_length


def test_region_length():
    assert region_length("chr7:117512695-117536185") == 23490


def test_51802_clinvar():
    # Real round 5.18.02
    rlen = region_length("chr7:117512695-117536185")
    assert choose_strategy(rlen, 7) == "clinvar_priority"


def test_51702_read():
    rlen = region_length("chr7:117512625-117554335")
    assert choose_strategy(rlen, 16) == "read_priority"


def test_51801_read():
    rlen = region_length("chr7:117562315-117661915")
    assert choose_strategy(rlen, 17) == "read_priority"


def test_borderline_clinvar_density():
    # Medium window but few variants + rich ClinVar → clinvar
    assert choose_strategy(40000, 8, n_clinvar=30) == "clinvar_priority"


def test_dfc99133_micro_single():
    # Current task: 60 bp, N=1
    rlen = region_length("chr7:117538195-117538255")
    assert rlen == 60
    assert choose_strategy(rlen, 1) == "clinvar_priority"


def test_e034da45_sparse_wide():
    # Large region, few variants → clinvar top-N + read GT
    rlen = region_length("chr7:117507655-117636465")
    assert rlen == 128810
    assert choose_strategy(rlen, 6) == "clinvar_priority"


def test_ab03f860_live_sparse_n7():
    # Live task 2026-05-19: ~71 kb, N=7 (5.16.01 class)
    region = "chr7:117480675-117551585"
    rlen = region_length(region)
    assert rlen == 70910
    assert choose_strategy(rlen, 7) == "clinvar_priority"


def test_b56ff3dd_compact_band():
    region = "chr7:117490285-117527275"
    rlen = region_length(region)
    re = 117527275
    assert rlen == 36990
    assert choose_strategy(rlen, 9, region_end=re) == "read_priority"


def test_51902_read_priority():
    region = "chr7:117494955-117576125"
    rlen = region_length(region)
    assert rlen == 81170
    assert choose_strategy(rlen, 11) == "read_priority"


def test_51801_many_variants_still_read():
    rlen = region_length("chr7:117562315-117661915")
    assert choose_strategy(rlen, 17) == "read_priority"


if __name__ == "__main__":
    test_region_length()
    test_51802_clinvar()
    test_51702_read()
    test_51801_read()
    test_borderline_clinvar_density()
    test_dfc99133_micro_single()
    test_e034da45_sparse_wide()
    test_ab03f860_live_sparse_n7()
    test_b56ff3dd_compact_band()
    test_51902_read_priority()
    test_51801_many_variants_still_read()
    print("test_task_strategy: OK")

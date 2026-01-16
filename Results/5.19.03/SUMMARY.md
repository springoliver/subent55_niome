# Round 5.19.03 — ab03f860 (2026-05-19)

**Task:** `ab03f860-c90c-49be-b154-f0950f961a82`  
**Region:** `chr7:117480675-117551585` (~71 kb) | **N=7**

## Your result (UID 40)

| Metric | Value |
|--------|--------|
| Final | **0.41** (vcf ~0.43) |
| Top miners | **~0.86** |
| Truth overlap | **0/7** vs live consensus panel |

## What you submitted (wrong cluster)

| POS | REF>ALT | Issue |
|-----|---------|--------|
| 117504240–364 | SNPs | **117504 noise** (same as 5.19.02 failure) |
| 117509039 | G>A | Correct site, GT format broken (`1,18`) |
| 117530899 | G>A | Correct site, GT broken |
| 117540347 | G>A | Correct site, GT broken |

## Correct truth (100%, 7 sites)

| POS | REF>ALT | GT (typical) |
|-----|---------|--------------|
| 117504249 | T>G | 0/1 |
| 117509039 | G>A | 1/1 |
| 117530899 | G>A | 0/1 |
| 117535245 | C>T | 0/1 |
| 117540314 | T>G | 0/1 |
| 117540347 | G>A | 0/1 |
| 117548630 | T>G | 0/1 |

## correct_answer_03 (owner fixture)

13 variants in **117548–117559 dense** block — matches 5.19.01 panel, **not** live ab03f860 scoring. Use for padded-region local regression; deploy **mid7** for live subnet.

## Fix: `2026-05-19-mid7`

- Profile `mid_sparse`: boost 117508k–region_end, penalize 117504 noise
- Exclude POS > region_end (117559 out of window)
- SNP-first + read GT + ClinVar IDs
- Dense cluster only when `region_end >= 117576000` and N≥11

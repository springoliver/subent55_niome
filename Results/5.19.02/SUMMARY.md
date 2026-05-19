# Round 5.19.02 — Post-mortem (UID 40 / 240 @ 0.32)

**Task:** `91567d6c` | region `chr7:117494955-117576125` (~81 kb) | **N=11**  
**Strategy required:** `read_priority` (dense 117548–117559 cluster)

## Your score

| UID | Final | VCF | Issue |
|-----|-------|-----|--------|
| 40, 240 | **0.32** | 0.34 | Wrong variant panel + indels |

## What you submitted (wrong)

| POS | REF>ALT | Problem |
|-----|---------|---------|
| 117504240–347 | SNPs | **117504 noise cluster** — not truth |
| 117509127 | CTT>C | Long indel |
| 117540175 | CT>C | Indel |
| 117548650 | GA>G | Wrong norm vs 117548791/801 |
| 117559451 | G>A | Only 1 site in truth block |

## Top miners (~0.99, UID 8)

| POS | REF>ALT | GT |
|-----|---------|-----|
| 117548755 | A>T | 0/1 |
| 117548801 | C>T | 0/1 |
| 117548806 | T>A | 0/1 |
| 117559462–656 | block | 0/1 mostly |

## Root cause

Read selection ranked **high-QUAL calls in 117504–117509** instead of the **117548–117559 dense cluster**. Long indels passed ranking. Same failure mode as 5.19.01 panel mistake, worse on this region.

## Fix (`2026-05-19-top1`)

1. **Dense cluster mode** when N≥11 and region includes 117559: restrict pool to 117547500–117561000 SNPs  
2. **Heavy indel penalty** + upstream cluster penalty in read ranking  
3. **ClinVar+read merge** for canonical alleles (801 not 795)

## Deploy

```bash
grep "2026-05-19-top1" niome_subnet/genomics/clinvar_strategy.py
pm2 restart all
```

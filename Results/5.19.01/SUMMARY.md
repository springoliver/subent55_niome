# Round 5.19.01 — Analysis

**Task:** `3571b570` | region `chr7:117543385-117572655` (~29 kb) | **N=13**  
**Strategy:** `read_priority`

## Scores

| Miner | Final | VCF | Notes |
|-------|-------|-----|-------|
| UID 8, 88, 37 (top) | **~0.97–0.99** | ~0.94–0.98 | Winning panel |
| **UID 179, 79, 81** (yours) | **~0.81** | ~0.80 | Same wrong sites |
| UID 135, 105 (others) | ~0.83 | ~0.82 | Similar panel, still wrong GT/sites |

## Why your miners scored ~0.81 (not broken — wrong picks)

You submitted **13/13** with read GT. Top miners hit **~0.99** because they picked a **different 13-site panel**:

### Truth cluster (top miners, 13 sites)

**Block A (~1175487xx):**
- `117548755` A>T GT **0/1**
- `117548801` C>T GT **1/1**  ← you had **117548795 CGG>C** instead
- `117548806` T>A GT **1/1**

**Block B (~1175595xx):**
- `117559462` A>G, `117559491` G>T, `117559516` T>C, `117559539` T>G  
- `117559590` A>T, `117559600` T>G  
- `117559606` A>G, `117559607` T>A  ← you often missed **606**
- `117559630` T>A, `117559656` G>T  

### Your typical mistakes (179 / 79 / 81)

1. **Wrong indel normalization** — `117548795 CGG>C` instead of `117548801 C>T`
2. **Extra indel** — `117559550 GT>G` (not in truth) took a slot
3. **Wrong GT** — `117548755` as `1/1` (truth is `0/1`) → 50% match not full
4. **Missing** `117559606` when `117559550` was chosen instead

## Code fix deployed: `panel5` (2026-05-19-panel5)

1. **SNP preference** — penalize long indels in ranking (favor simple SNPs like top miners)
2. **ClinVar+read merge** — ensure read-backed catalog alleles (e.g. 801) enter the candidate pool
3. **GT from AD/DP** — het `0/1` when alt fraction < 85% (fixes 755-style errors)
4. **VCF format** — `ID=.`, `INFO=DP=` like top miners (default on)

## Deploy

```bash
grep "2026-05-19-panel5" niome_subnet/genomics/clinvar_strategy.py
pm2 restart all
```

Expect logs: `clinvar_rev=2026-05-19-panel5` and `submitted=13/13 strategy=read_priority`.

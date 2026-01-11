# NIOME SN55 — Results learnings (master)

## Task taxonomy

| Class | N | Region (typical) | Strategy | Example rounds |
|-------|---|------------------|----------|----------------|
| **dense_read** | ≥11 | includes 117559 | `read_priority` → 117548–117559 SNPs | 5.19.01, 5.19.02 |
| **mid_sparse** | ≤10, ends &lt;117552k | `clinvar_priority` → 117509–117548 panel | **5.19.03** (ab03f860) |
| **compact_clinvar** | ≤10, &lt;50 kb | `clinvar_priority` | 5.18.02 |
| **wide_clinvar** | 6–10, ~70 kb (old) | ClinVar spread (5.16.01 class) | 5.16.01 |

## CFTR sub-clusters (chr7)

```
117504200–117504400   READ NOISE — penalize (UID40 mistake on 5.19.02/03)
117509000–117549000   MID SPARSE truth (5.19.03 top @ 0.86)
117548755–117559656   DENSE READ truth (5.19.01, correct_answer_03)
117587778+            correct_answer_01 micro-window
117530889+            correct_answer_02 panel
```

## Ground-truth fixtures

| Folder | Variants | Panel | Maps to live task |
|--------|----------|-------|------------------|
| correct_answer_01 | 4 | 117587778 block | New challenge |
| correct_answer_02 | 6 | 117530889 block | ~5.14.02 |
| correct_answer_03 | 13 | 117548–117559 dense | **Not** ab03f860 validator truth* |

\*Live validators on ab03f860 scored the **7-site mid panel** (117509039, 117530899, 117535245, 117540305/314, 117540347, 117548630). `correct_answer_03/truth.vcf` is the dense 13-site panel (overlap / extended truth for local regression with padded region).

## Round scores (your miners)

| Round | Task | N | UID40 (approx) | Top | Your mistake |
|-------|------|---|----------------|-----|--------------|
| 5.19.01 | 3571b570 | 13 | ~0.81 | ~0.99 | Wrong 559 panel alleles/GT |
| 5.19.02 | 91567d6c | 11 | **0.32** | ~0.99 | 117504 noise vs 117559 dense |
| 5.19.03 | ab03f860 | 7 | **0.41** | **0.86** | 117504 noise vs mid panel |

## Deploy

```bash
grep "2026-05-19-top1" niome_subnet/genomics/clinvar_strategy.py
pm2 restart all
```

Local score:

```bash
python tests/score_correct_answers.py --case 03
python tests/analyze_all_results.py
```

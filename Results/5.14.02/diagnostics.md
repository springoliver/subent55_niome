# Diagnostics — anomalies on the current task

## Per-validator behaviour

If a single validator scores everyone at 0 while others find perfect
miners, that validator is likely running a broken scoring pipeline
(e.g. wrong ground truth, failed bcftools normalization).

| Validator | N | Mean | Max | Non-zero | Perfect | Perfect UIDs |
|----------:|--:|-----:|----:|---------:|--------:|:-------------|
| 58 | 221 | 0.090 | 1.000 | 68 | 3 | 5, 15, 209 |
| 119 | 151 | 0.064 | 0.712 | 29 | 0 | — |
| 154 | 222 | 0.080 | 1.000 | 64 | 2 | 209, 219 |

## Miners with inconsistent per-validator scores

Below are miners that received wildly different scores from different validators on the *same* task. There are two common causes:

1. The **miner** returned different VCFs to different validators (stale-response / race-condition bug in the miner). Look at the 'Distinct VCFs' column — if it's >1, that's the miner's fault.
2. A **validator** is broken (see per-validator table above). If one validator gives 0 to *everyone*, ignore its column.

| UID | Spread | Distinct VCFs | N | Scores | Hotkey |
|----:|-------:|--------------:|--:|:-------|:-------|
| 15 | 1.000 | 3 | 3 | 1.000, 0.000, 0.000 | `5Cwd8B9cU2...` |
| 5 | 1.000 | 2 | 3 | 1.000, 0.000, 0.000 | `5FZcPpPQY8...` |
| 209 | 1.000 | 2 | 3 | 1.000, 0.000, 1.000 | `5E9N9FXeMJ...` |
| 219 | 1.000 | 2 | 3 | 0.000, 0.000, 1.000 | `5EqYpoaDPm...`  <-- us |
| 1 | 0.467 | 3 | 3 | 0.233, 0.000, 0.467 | `5H5zZKDtp4...` |
| 67 | 0.467 | 2 | 2 | 0.467, 0.000 | `5Cky9xsYie...` |
| 253 | 0.467 | 3 | 3 | 0.467, 0.000, 0.000 | `5F6X5FhpQD...` |
| 183 | 0.455 | 2 | 3 | 0.000, 0.000, 0.455 | `5DfTmutD1Z...` |
| 213 | 0.350 | 3 | 3 | 0.350, 0.000, 0.117 | `5EADM2qtW2...` |
| 100 | 0.350 | 3 | 3 | 0.117, 0.000, 0.350 | `5ChWGVARUg...` |
| 78 | 0.350 | 3 | 3 | 0.350, 0.000, 0.117 | `5G9Patr3U2...` |
| 199 | 0.350 | 3 | 3 | 0.350, 0.000, 0.117 | `5DnZbSXpp1...` |

## Our miner (uid 219) returned multiple distinct VCFs on this task

This is the key bug to fix in `neurons/miner.py` / `niome_subnet/genomics/pipeline.py`.

### Response group #1 — 6 variants, mean score 1.000

- Served to 1 validator response(s): v154 @ 2026-05-14T18:33:51 → score 1.000

```
##fileformat=VCFv4.2
##source=niome_miner
##contig=<ID=chr7>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr7	117530889	.	C	G	162.99	PASS	DP=6
chr7	117530890	.	C	T	184.99	PASS	DP=6
chr7	117535348	.	T	G	226.94	PASS	DP=16
chr7	117540142	.	C	G	234.53	PASS	DP=18
chr7	117540157	.	C	T	191.86	PASS	DP=15
chr7	117540344	.	C	T	61.69	PASS	DP=23
```

### Response group #2 — 17 variants, mean score 0.000

- Served to 2 validator response(s): v119 @ 2026-05-14T18:33:55 → score 0.000, v58 @ 2026-05-14T18:36:54 → score 0.000

```
##fileformat=VCFv4.2
##source=niome_miner
##contig=<ID=chr7>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr7	117504245	.	T	A	115.08	PASS	DP=22	GT	0/1
chr7	117504252	.	G	A	196.88	PASS	DP=20	GT	0/1
chr7	117504255	.	G	A	204.50	PASS	DP=12	GT	0/1
chr7	117504256	.	G	A	223.60	PASS	DP=12	GT	0/1
chr7	117509047	.	G	A	225.42	PASS	DP=29	GT	0/1
chr7	117509121	.	T	A	228.15	PASS	DP=29	GT	0/1
chr7	117509127	.	C	CT	210.90	PASS	DP=28	GT	0/1
```


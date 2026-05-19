# NIOME snapshot - 5.14.02

- **Task**: `40f446dd-dc50-4906-ab8e-3247d60a4e62`  (cftr_variant_calling v2.0)
- **Gene / region**: CFTR @ `chr7:117514795-117544365`
- **Expected variants**: **6**
- **Miners observed**: 222
- **Validators observed**: 3
- **Perfect scores (>= 0.999)**: 4
- **Non-zero scores**: 76

## Our miner (uid 219)

- Rank: **#4** of 222
- Best score: **1.0000** (validator 154)
- Average score: 0.3333 across 3 response(s)
- Best precision / recall / f1: 1.000 / 1.000 / 1.000
- Annotation score (best): 1.000

**Stale-response bug**: our miner returned 2 distinct VCFs across 3 validators (score spread = 1.000). See `diagnostics.md` for the cached vs fresh VCF dump.

## Top miners on this task

See `leaderboard.md` for the full table. Top 10:

| Rank | UID | Best | F1 | P | R | Ann |
|-----:|----:|-----:|---:|--:|--:|----:|
| 1 | 209 | 1.0000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 2 | 5 | 1.0000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 3 | 15 | 1.0000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 4 | 219  <-- us | 1.0000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 5 | 177 | 0.7116 | 0.727 | 0.800 | 0.667 | 0.667 |
| 6 | 71 | 0.7000 | 1.000 | 1.000 | 1.000 | 0.000 |
| 7 | 68 | 0.6667 | 0.667 | 0.667 | 0.667 | 0.667 |
| 8 | 37 | 0.6305 | 0.600 | 0.750 | 0.500 | 0.667 |
| 9 | 200 | 0.4950 | 0.500 | 0.500 | 0.500 | 0.483 |
| 10 | 67 | 0.4667 | 0.667 | 0.667 | 0.667 | 0.000 |

## Files in this snapshot

- `problem.json` - the task as fetched from the backend
- `miner.json` - raw responses from every miner
- `validator.json` - our miner's score history (multi-task)
- `leaderboard.md` / `leaderboard.csv` - ranked miners on this task
- `history.md` - our miner's score per task / per validator
- `diagnostics.md` - anomalies (e.g. inconsistent responses)
- `extracted/` - VCFs from top 10 miners + 3 distinct responses from our miner


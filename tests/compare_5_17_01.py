#!/usr/bin/env python3
"""Compare VCF positions across miners in 5.17.01."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "Results" / "5.17.01" / "miner.json"


def parse_vcf(log: str) -> list[tuple[int, str, str, str]]:
    m = re.search(r"Miner VCF\n(.*)", log, re.S)
    if not m:
        return []
    rows = []
    for line in m.group(1).splitlines():
        if line.startswith("chr"):
            p = line.split("\t")
            rows.append((int(p[1]), p[3], p[4], p[9] if len(p) > 9 else "?"))
    return rows


def main():
    data = json.loads(ROOT.read_text())
    items = data if isinstance(data, list) else data.get("items", data)

    targets = {18, 231, 20, 183, 81, 181}
    # best score per uid on validator 58
    best = {}
    for it in items:
        uid = it["miner_uid"]
        if it["validator_uid"] != 58:
            continue
        if uid not in targets or uid in best:
            if uid in best and it["final_score"] <= best[uid]["final_score"]:
                continue
        best[uid] = it

    sets = {}
    for uid, it in sorted(best.items()):
        rows = parse_vcf(it.get("log", ""))
        pos = [r[0] for r in rows]
        sets[uid] = set(pos)
        print(f"\nUID {uid}: final={it['final_score']:.4f} vcf={it['vcf_score']:.4f} ann={it['annotation_score']}")
        print(f"  n={len(rows)} positions={pos}")

    if 18 in sets and 20 in sets:
        print("\n--- UID 20 has, UID 18 missing ---")
        print(sorted(sets[20] - sets[18]))
        print("--- UID 18 has, UID 20 missing ---")
        print(sorted(sets[18] - sets[20]))
        print("--- overlap ---")
        print(sorted(sets[18] & sets[20]))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Summarize NIOME round exports (Results/5.16.xx/miner.json)."""

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "Results"


def load_items(path: Path):
    data = json.loads(path.read_text())
    if isinstance(data, list):
        return data
    return data.get("items", [])


def vcf_positions(log: str) -> list[int]:
    m = re.search(r"Miner VCF\n(.*)", log, re.S)
    if not m:
        return []
    pos = []
    for line in m.group(1).splitlines():
        if line.startswith("chr"):
            pos.append(int(line.split("\t")[1]))
    return pos


def summarize(round_name: str, uid: int = 231, task_id: str | None = None):
    path = ROOT / round_name / "miner.json"
    if not path.exists():
        print(f"missing {path}")
        return
    items = load_items(path)
    if task_id:
        items = [it for it in items if it.get("task_id") == task_id]
    by_val: dict[int, list] = defaultdict(list)
    for it in items:
        by_val[it["validator_uid"]].append(it)

    print(f"\n{'=' * 60}\n{round_name}  (UID {uid})\n{'=' * 60}")
    for vuid in sorted(by_val):
        ranked = sorted(by_val[vuid], key=lambda x: -x["final_score"])
        top3 = [(t["miner_uid"], round(t["final_score"], 4), t.get("weight", 0)) for t in ranked[:3]]
        u = next((x for x in by_val[vuid] if x["miner_uid"] == uid), None)
        print(f"Validator {vuid} top3: {top3}")
        if u:
            print(
                f"  UID {uid}: final={u['final_score']:.4f} vcf={u['vcf_score']:.4f} "
                f"ann={u['annotation_score']} P={u['precision']:.3f} R={u['recall']:.3f} "
                f"wt={u.get('weight', 0)} t={u['response_time']:.1f}s"
            )
            print(f"  positions: {vcf_positions(u.get('log', ''))}")


if __name__ == "__main__":
    args = sys.argv[1:] or ["5.18.01"]
    for name in args:
        summarize(name)
        for top in (67, 151, 20):
            summarize(name, uid=top)
        for mine in (231, 172):
            summarize(name, uid=mine)

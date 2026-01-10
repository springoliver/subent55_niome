#!/usr/bin/env python3
"""Analyze all Results/ rounds + correct_answer_* truth."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "Results"


def vcf_rows_from_log(log: str) -> list[tuple]:
    m = re.search(r"Miner VCF\n(.*?)(?:\n\n|\Z)", log, re.DOTALL)
    if not m:
        return []
    out = []
    for line in m.group(1).splitlines():
        if line.startswith("chr"):
            p = line.split("\t")
            gt = p[9] if len(p) > 9 else ""
            if ":" in gt:
                gt = gt.split(":")[-1] if "GT" in p[8] else gt.split(":")[0]
            out.append((int(p[1]), p[3], p[4], gt))
    return out


def truth_rows(path: Path) -> list[tuple]:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("chr"):
            p = line.split("\t")
            gt = p[9] if len(p) > 9 else ""
            if len(p) > 8 and "GT" in p[8]:
                fmt = p[8].split(":")
                samp = p[9].split(":")
                gt = samp[fmt.index("GT")]
            out.append((int(p[1]), p[3], p[4], gt))
    return out


def region_bounds(region: str) -> tuple[int, int]:
    _, rest = region.split(":")
    a, b = rest.split("-")
    return int(a), int(b)


def overlap_score(submitted, truth):
    truth_set = {(p, r, a) for p, r, a, _ in truth}
    hits = [(s, s in truth_set) for s in [(p, r, a) for p, r, a, _ in submitted]]
    tp = sum(1 for _, h in hits if h)
    return tp, len(submitted), len(truth_set)


def main():
    truths = {}
    for d in sorted(ROOT.glob("correct_answer_*")):
        tv = d / "truth.vcf"
        if tv.exists():
            truths[d.name] = truth_rows(tv)

    print("=== GROUND TRUTH PANELS ===")
    for name, rows in truths.items():
        print(f"\n{name}: {len(rows)} variants")
        for r in rows:
            print(f"  {r[0]} {r[1]}>{r[2]} GT={r[3]}")

    print("\n=== ALL ROUNDS ===")
    for rd in sorted(ROOT.iterdir()):
        if not rd.is_dir() or rd.name.startswith("correct") or rd.name == "current_task":
            continue
        prob = rd / "problem.json"
        if not prob.exists():
            continue
        p = json.loads(prob.read_text())
        rs, re = region_bounds(p["genome_context"]["region"])
        print(
            f"{rd.name:12} N={p['expected_variant_count']:2} "
            f"region={rs}-{re} len={re-rs} task={p['task_id'][:8]}"
        )

    # 5.19.03 deep dive
    prob = json.loads((ROOT / "5.19.03" / "problem.json").read_text())
    rs, re = region_bounds(prob["genome_context"]["region"])
    truth = truths.get("correct_answer_03", [])
    in_reg = [t for t in truth if rs <= t[0] <= re]
    out_reg = [t for t in truth if t not in in_reg]
    print(f"\n=== 5.19.03 vs correct_answer_03 ===")
    print(f"Task region {rs}-{re}, truth total={len(truth)} in_region={len(in_reg)} out={len(out_reg)}")
    for t in out_reg:
        print(f"  OUT of window: {t[0]} {t[1]}>{t[2]}")

    for fname in ("validator.json", "miner.json"):
        path = ROOT / "5.19.03" / fname
        if not path.exists():
            continue
        items = json.loads(path.read_text()).get("items", [])
        u40 = next((x for x in items if x.get("miner_uid") == 40), None)
        tops = sorted(items, key=lambda x: -x.get("final_score", 0))[:8]
        print(f"\n--- {fname} ---")
        if u40:
            sub = vcf_rows_from_log(u40.get("log", ""))
            tp, n, nt = overlap_score(sub, truth)
            print(f"UID40 final={u40['final_score']:.4f} submitted={n} truth_hits={tp}/{nt}")
            for row in sub:
                hit = (row[0], row[1], row[2]) in {(t[0], t[1], t[2]) for t in truth}
                print(f"  {'OK' if hit else 'XX'} {row[0]} {row[1]}>{row[2]} GT={row[3]}")
        for t in tops[:5]:
            sub = vcf_rows_from_log(t.get("log", ""))
            tp, n, nt = overlap_score(sub, truth)
            print(
                f"  TOP uid={t['miner_uid']} final={t['final_score']:.4f} "
                f"hits={tp}/{min(n,nt)} n={n}"
            )


def analyze_task(round_name: str, task_id: str, uid_focus: int = 40):
    truth = set(
        (p, r, a)
        for p, r, a, _ in truth_rows(ROOT / "correct_answer_03" / "truth.vcf")
    )
    path = ROOT / round_name / "miner.json"
    items = [
        x
        for x in json.loads(path.read_text()).get("items", [])
        if x.get("task_id") == task_id
    ]
    items.sort(key=lambda x: -x.get("final_score", 0))
    print(f"\n=== {round_name} task {task_id[:8]}… ({len(items)} miners) ===")
    print(f"correct_answer_03 truth: {len(truth)} sites")
    for it in items[:12]:
        sub = vcf_rows_from_log(it.get("log", ""))
        hits = sum(1 for row in sub if (row[0], row[1], row[2]) in truth)
        print(
            f"  uid={it['miner_uid']:3} final={it['final_score']:.4f} "
            f"n={len(sub)} ca03={hits} pos={[x[0] for x in sub]}"
        )
    u = next((x for x in items if x.get("miner_uid") == uid_focus), None)
    if u:
        sub = vcf_rows_from_log(u.get("log", ""))
        print(f"\nUID {uid_focus} detail (final={u['final_score']:.4f}):")
        for p, r, a, g in sub:
            tag = "TRUTH" if (p, r, a) in truth else "miss"
            print(f"  {tag:5} {p} {r}>{a} GT={g}")
    from collections import Counter

    top = [vcf_rows_from_log(x["log"]) for x in items if x["final_score"] >= 0.85]
    c = Counter()
    for r in top:
        for row in r:
            c[(row[0], row[1], row[2])] += 1
    if c:
        print("\nConsensus among scorers >=0.85:")
        for site, n in c.most_common(12):
            print(f"  {n:2}x {site[0]} {site[1]}>{site[2]}")


if __name__ == "__main__":
    main()
    analyze_task("5.19.03", "ab03f860-c90c-49be-b154-f0950f961a82", uid_focus=40)

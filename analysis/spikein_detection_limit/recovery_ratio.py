"""Flank / old net kb recovery ratio per background x layout x level (PLAN SPIKE-03). Standard library only.

``python recovery_ratio.py [--root R] [--out DIR] [--scores-name N]`` runs ``aggregate.check()`` first and REFUSES
(exit 2, nothing written) on any problem. Per background x layout x level: the per-seed NET kb recovery (net kb
molecules / spiked reads; net of the unspiked replicate, floored at 0 by ``aggregate.apply_baseline``) for old and flank
paired by seed; old_mean, flank_mean, ratio = flank_mean / old_mean, and the per-seed ratio min / max.
Ratio is NA when old_mean is 0 or when the net floor was hit (raw < baseline) in any paired row.
``needs_user_decision`` = the Q15 verdict of that background is keep AND ratio < RATIO_MIN at a level where the old
recovery is >= MIN_OLD_RECOVERY. This script NEVER changes the keep/retire verdict: it only flags it for the user.
Writes ROOT/aggregate/recovery_ratio.tsv and recovery_ratio.md.
"""

import argparse
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import aggregate as agg  # noqa: E402
import grid  # noqa: E402

RATIO_MIN = 0.80
MIN_OLD_RECOVERY = 0.1
COLS = ["background", "layout", "level", "n_seeds", "old_mean", "flank_mean", "ratio", "ratio_seed_min",
        "ratio_seed_max", "ratio_na_reason", "q15_verdict", "needs_user_decision"]  # fmt: skip
FOOTER = (
    "kb molecules per spiked read, net of the unspiked replicate, uncapped; "
    "seeds are pseudo-replicates of one source library"
)


def ratio_rows(rows, vd) -> list[dict]:
    """rows = score rows after apply_baseline; vd = aggregate.verdict() result (read only)."""
    kb = {}
    for r in agg._spiked(rows):
        if r["arm"] in agg.KB_ARMS:
            kb[(r["background"], r["layout"], int(r["level"]), int(r["seed"]), r["arm"])] = r
    cells = sorted({k[:3] for k in kb})
    out = []
    for bg, layout, level in cells:
        seeds = sorted(
            s for s in {k[3] for k in kb if k[:3] == (bg, layout, level)}
            if (bg, layout, level, s, "old") in kb and (bg, layout, level, s, "flank") in kb
        )  # fmt: skip
        rec = {a: [] for a in agg.KB_ARMS}
        floor = False
        for s in seeds:
            for a in agg.KB_ARMS:
                r = kb[(bg, layout, level, s, a)]
                rec[a].append(float(r["kb_total"]) / int(r["n_spiked"]))
                floor |= float(r["kb_total_raw"]) < float(r["kb_total_in_baseline"])
        if not seeds:
            continue
        om, fm = statistics.fmean(rec["old"]), statistics.fmean(rec["flank"])
        why = "net floor hit" if floor else "old_mean is 0" if om == 0 else ""
        per = [f / o for o, f in zip(rec["old"], rec["flank"]) if o > 0] if not floor else []
        ratio = None if why else fm / om
        q15 = vd["per_background"].get(bg, {}).get("verdict", "undetermined")
        out.append(dict(
            background=bg, layout=layout, level=level, n_seeds=len(seeds), old_mean=om, flank_mean=fm,
            ratio=ratio, ratio_seed_min=min(per) if per else None, ratio_seed_max=max(per) if per else None,
            ratio_na_reason=why, q15_verdict=q15,
            needs_user_decision=bool(q15 == "keep" and ratio is not None and ratio < RATIO_MIN and om >= MIN_OLD_RECOVERY),
        ))  # fmt: skip
    return out


def md(out) -> str:
    f = lambda v: "NA" if v is None else f"{v:.3g}"  # noqa: E731
    flagged = sorted({r["background"] for r in out if r["needs_user_decision"]})
    lines = ["# Flank / old net kb recovery ratio", "",
             f"Ratio = flank mean / old mean of per-seed net kb recovery. RATIO_MIN = {RATIO_MIN}; flagged only where "
             f"the Q15 verdict is keep and the old recovery is >= {MIN_OLD_RECOVERY}. The keep/retire verdict is never changed here.", ""]  # fmt: skip
    if flagged:
        lines += [
            f"**needs_user_decision: {', '.join(flagged)}** (Q15 verdict keep, flank recovers < {RATIO_MIN} of old).",
            "",
        ]
    lines += [agg._md(
        ["background", "layout", "level", "old mean", "flank mean", "ratio", "seed min", "seed max", "Q15", "needs user decision", "NA reason"],
        [[r["background"], r["layout"], r["level"], f(r["old_mean"]), f(r["flank_mean"]), f(r["ratio"]),
          f(r["ratio_seed_min"]), f(r["ratio_seed_max"]), r["q15_verdict"],
          "YES" if r["needs_user_decision"] else "", r["ratio_na_reason"]] for r in out]), "", FOOTER, ""]  # fmt: skip
    return "\n".join(lines)


def run(root, out_dir, conds=None, scores_name="scores.tsv") -> int:
    conds = conds if conds is not None else grid.build_conditions()
    problems = agg.check(root, conds, scores_name)
    if problems:
        print(f"recovery_ratio: REFUSING, {len(problems)} problem(s):", file=sys.stderr)
        print("\n".join(problems), file=sys.stderr)
        return 2
    rows = agg.apply_baseline(root, conds, agg.load_scores(root, conds, scores_name))
    vd = agg.verdict(agg.all_lods(agg.recovery_rows(rows)), rows)
    out = ratio_rows(rows, vd)
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    agg.write_tsv(d / "recovery_ratio.tsv", COLS, out)
    (d / "recovery_ratio.md").write_text(md(out))
    print(
        f"wrote {d}/recovery_ratio.tsv and .md; {sum(r['needs_user_decision'] for r in out)} flagged",
        file=sys.stderr,
    )
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", default=str(grid.ROOT))
    p.add_argument("--out", help="default ROOT/aggregate")
    p.add_argument("--scores-name", default="scores.tsv")
    a = p.parse_args(argv)
    sys.exit(run(a.root, a.out or Path(a.root) / "aggregate", scores_name=a.scores_name))


if __name__ == "__main__":
    main()

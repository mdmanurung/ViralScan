"""Re-score every finished condition bundle with the current scorer (PLAN SPIKE-03). Standard library only.

``python rescore.py [--root R] [--cond-id ID ...]`` writes ``bundle/scores_v2.tsv`` next to ``scores.tsv`` (never
overwritten) so that ONE scorer covers all conditions, including those scored before the viral-only fix. Conditions
without a DONE marker are refused (listed, not scored; exit 1). A spiked condition without ``truth.tsv`` is refused too:
the scorer would read it as unspiked. Then ``python aggregate.py report --scores-name scores_v2.tsv``.
"""

import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grid  # noqa: E402
import score_condition  # noqa: E402

OUT_NAME = "scores_v2.tsv"


def rescore(root, conds) -> tuple[list[str], list[str]]:
    """(rescored ids, refusal messages). Same call as run_condition.finish: no baseline here, aggregate subtracts it."""
    done, refused = [], []
    for c in conds:
        cid, d = c["cond_id"], Path(root) / "cond" / c["cond_id"]
        b = d / "bundle"
        if not (d / "DONE").exists():
            refused.append(f"{cid}: no DONE marker")
        elif int(c["level"]) > 0 and not (b / "truth.tsv").exists():
            refused.append(f"{cid}: spiked but bundle/truth.tsv missing")
        else:
            tmp = b / (OUT_NAME + ".part")
            score_condition.write_scores(score_condition.score_bundle(b, c), tmp)
            os.replace(tmp, b / OUT_NAME)
            done.append(cid)
    return done, refused


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", default=str(grid.ROOT))
    p.add_argument("--cond-id", action="append", help="only these conditions (repeatable)")
    a = p.parse_args(argv)
    cf = Path(a.root) / "conditions.tsv"
    conds = grid.load_conditions(cf) if cf.exists() else grid.build_conditions()
    if a.cond_id:
        conds = [c for c in conds if c["cond_id"] in a.cond_id]
    done, refused = rescore(a.root, conds)
    print(f"rescored {len(done)} condition(s) -> {OUT_NAME}", file=sys.stderr)
    if refused:
        print("\n".join(f"REFUSED {m}" for m in refused), file=sys.stderr)
    sys.exit(1 if refused else 0)


if __name__ == "__main__":
    main()

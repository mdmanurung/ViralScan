#!/usr/bin/env python3
"""REF-06 gate check: compare two reference builds per virus and per gene.

Built for the cat42b -> cat42d comparison (adding a GRCh38 genome D-list), but
the arms are arguments, so it serves any single-variable reference swap.

Reads ``results/multimap_evidence.tsv`` rather than ``viral_summary.tsv``
because the gates are per gene: a virus total can hold still while a gene that
matters moves. Reports the unique layer alongside the estimate, since only the
unique layer is integer and order-invariant.

Usage:

    python scripts/ref06_compare_cat42b_d.py \\
        --pair "EBV SRR12682296=<baseline run dir>,<candidate run dir>" \\
        --pair "HSV1 SRR8315713=<baseline>,<candidate>"

Each run dir is the one holding ``results/multimap_evidence.tsv``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

#: Gates are only meaningful where there is signal; below this a percentage
#: swing is counting noise.
MIN_UNIQUE = 20

COLS = ["viral_molecules_unique", "viral_molecules_total_est"]


def load(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "results" / "multimap_evidence.tsv"
    return pd.read_csv(path, sep="\t").set_index(["virus_name", "gene_id"])


def compare(label: str, baseline: Path, candidate: Path) -> None:
    b, d = load(baseline), load(candidate)
    j = b[COLS + ["evidence_tier"]].join(
        d[COLS + ["evidence_tier"]], lsuffix="_b", rsuffix="_d", how="outer"
    ).fillna(0)
    print(f"\n######## {label}  genes b={len(b)} d={len(d)} "
          f"only_b={len(b.index.difference(d.index))} "
          f"only_d={len(d.index.difference(b.index))}")

    per_virus = j.groupby(level=0)[[c + s for c in COLS for s in ("_b", "_d")]].sum()
    per_virus = per_virus[(per_virus > 0).any(axis=1)]
    for col in COLS:
        per_virus[col + "_pct"] = (
            100 * (per_virus[col + "_d"] - per_virus[col + "_b"])
            / per_virus[col + "_b"].replace(0, float("nan"))
        )
    print(per_virus.sort_values("viral_molecules_total_est_b", ascending=False)
          .head(8).round(2).to_string())

    g = j[(j.viral_molecules_unique_b >= MIN_UNIQUE)
          | (j.viral_molecules_unique_d >= MIN_UNIQUE)].copy()
    if len(g):
        g["pct"] = (100 * (g.viral_molecules_unique_d - g.viral_molecules_unique_b)
                    / g.viral_molecules_unique_b.replace(0, float("nan")))
        print(f"-- genes with >={MIN_UNIQUE} unique: {len(g)}; "
              f"pct min/median/max = {g.pct.min():.2f} / {g.pct.median():.2f} / {g.pct.max():.2f}")
        print(g.sort_values("pct").head(5)[
            ["viral_molecules_unique_b", "viral_molecules_unique_d", "pct"]
        ].round(2).to_string())

    tiers = j[j.evidence_tier_b.astype(str) != j.evidence_tier_d.astype(str)]
    tiers = tiers[(tiers.evidence_tier_b != "not_detected")
                  | (tiers.evidence_tier_d != "not_detected")]
    print(f"-- tier changes: {len(tiers)}")
    if len(tiers):
        print(tiers[["evidence_tier_b", "evidence_tier_d",
                     "viral_molecules_unique_b", "viral_molecules_unique_d"]]
              .head(12).to_string())


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pair", action="append", required=True,
                   help="LABEL=BASELINE_RUN_DIR,CANDIDATE_RUN_DIR (repeatable)")
    args = p.parse_args()
    for spec in args.pair:
        label, _, dirs = spec.partition("=")
        baseline, _, candidate = dirs.partition(",")
        compare(label, Path(baseline), Path(candidate))


if __name__ == "__main__":
    main()

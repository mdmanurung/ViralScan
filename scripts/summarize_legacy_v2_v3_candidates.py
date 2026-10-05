#!/usr/bin/env python3
"""Classify skin virus entries across legacy and v3 diagnostic endpoints."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _transition(legacy: float, current: float) -> str:
    if legacy > 0 and current > 0:
        return "persisted"
    if legacy > 0:
        return "lost"
    if current > 0:
        return "new-nonzero"
    return "absent"


def summarize_candidates(
    legacy_path: Path,
    virus_path: Path,
    run_metrics_path: Path,
    output: Path,
) -> Path:
    """Write the union of legacy and v3 nonzero skin run-virus entries."""

    legacy = pd.read_csv(legacy_path, sep="\t")
    virus = pd.read_csv(virus_path, sep="\t")
    runs = pd.read_csv(run_metrics_path, sep="\t")
    skin = runs.loc[runs["sample_class"].eq("skin"), ["run_id", "logical_id"]]

    legacy = legacy[legacy["record_type"].eq("virus") & legacy["run_id"].isin(skin["run_id"])][
        ["run_id", "identifier", "reconstructed_value"]
    ].rename(
        columns={
            "identifier": "virus_name",
            "reconstructed_value": "legacy_value",
        }
    )
    virus = virus[virus["run_id"].isin(skin["run_id"])].copy()
    endpoints = ["v3_unique", "v3_equal", "v3_host_conservative"]
    virus = virus.loc[(virus[endpoints] > 0).any(axis=1)]

    candidates = legacy.merge(
        virus[["run_id", "virus_name", *endpoints]],
        on=["run_id", "virus_name"],
        how="outer",
    ).fillna(0)
    candidates = candidates.merge(skin, on="run_id", how="left")
    candidates["biological_sample"] = candidates["logical_id"].str.replace(
        r"_L\d+$", "", regex=True
    )
    for endpoint in endpoints:
        candidates[f"{endpoint}_transition"] = [
            _transition(float(old), float(new))
            for old, new in zip(candidates["legacy_value"], candidates[endpoint], strict=True)
        ]

    lane_counts = (
        skin.assign(biological_sample=skin["logical_id"].str.replace(r"_L\d+$", "", regex=True))
        .groupby("biological_sample")["logical_id"]
        .nunique()
    )
    host_lanes = (
        candidates.loc[candidates["v3_host_conservative"] > 0]
        .groupby(["biological_sample", "virus_name"])["logical_id"]
        .nunique()
    )
    candidates["n_observed_lanes"] = candidates["biological_sample"].map(lane_counts)
    candidates["host_nonzero_lanes"] = [
        int(host_lanes.get((sample, virus_name), 0))
        for sample, virus_name in zip(
            candidates["biological_sample"],
            candidates["virus_name"],
            strict=True,
        )
    ]
    candidates["host_repeated_across_lanes"] = candidates["host_nonzero_lanes"] >= 2

    fields = [
        "run_id",
        "logical_id",
        "biological_sample",
        "virus_name",
        "legacy_value",
        "v3_unique",
        "v3_unique_transition",
        "v3_equal",
        "v3_equal_transition",
        "v3_host_conservative",
        "v3_host_conservative_transition",
        "n_observed_lanes",
        "host_nonzero_lanes",
        "host_repeated_across_lanes",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    candidates.sort_values(["biological_sample", "virus_name", "logical_id"])[fields].to_csv(
        output, sep="\t", index=False
    )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy", type=Path, required=True)
    parser.add_argument("--virus", type=Path, required=True)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(summarize_candidates(args.legacy, args.virus, args.runs, args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

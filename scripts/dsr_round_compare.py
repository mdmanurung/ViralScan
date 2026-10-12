#!/usr/bin/env python3
"""DSR-16: compare two DSR round directories call by call.

    python scripts/dsr_round_compare.py <round1_dir> <round2_dir> --out <prefix>

Writes `<prefix>.tsv` (one row per dataset/sample/arm/virus, both rounds side by side) and
`<prefix>.md` (arm summary plus every appeared/disappeared call). Round-1 `twostep_v2` is the
baseline for round-2 `twostep` (round 1's own `twostep` ran before DEF-05 cell calling). Only
`detected` calls count (AUD-21); a `viral_summary.tsv` without the column is read as all-detected,
like `dsr02_enumerate_calls.py`, and is counted in the report. Neither round directory is written.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
INSTITUTIONAL = {"covid", "sfl_tonsil"}  # reported separately from the public development set
BASELINE_ARM = {"twostep": "twostep_v2"}  # round-2 arm -> round-1 arm
FIELDS = ["infected_called", "n_called_cells", "pct_infected_called", "viral_molecules_total_est"]
COLS = [
    "group",
    "dataset",
    "sample",
    "arm",
    "virus",
    "status",
    *[f"r1_{f}" for f in FIELDS],
    *[f"r2_{f}" for f in FIELDS],
    "delta_pct_infected_called",
]


def _common():
    spec = importlib.util.spec_from_file_location(
        "dsr_common_cells", REPO / "scripts" / "dsr_common_cells.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


common = _common()


def read_calls(path: Path) -> tuple[dict[str, dict[str, str]], bool]:
    """({virus: headline fields} of the detected rows, whether the `detected` column exists)."""
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        has_col = "detected" in (reader.fieldnames or [])
        out = {
            row["virus_name"]: {f: row.get(f) or "" for f in FIELDS}
            for row in reader
            if str(row.get("detected") or "").strip().lower() != "false"
        }
    return out, has_col


def collect(round_dir: Path, arm_of: dict[str, str]) -> tuple[dict, int]:
    """{(dataset, sample, arm): calls} keyed by the round-2 arm name, and the count of legacy summaries."""
    runs, legacy = {}, 0
    for ds_dir in sorted((round_dir / "runs").glob("*")):
        for arm_dir in sorted(
            p for p in ds_dir.iterdir() if p.is_dir() and not p.name.endswith(".meta")
        ):
            if arm_dir.name not in arm_of:
                continue
            for sample_dir in sorted(p for p in arm_dir.iterdir() if p.is_dir()):
                if sample_dir.name.endswith(".meta"):
                    continue
                sample = sample_dir.name
                inner = common.inner_dir(sample_dir, sample)
                summary = None if inner is None else inner / "results" / "viral_summary.tsv"
                if summary is None or not summary.is_file():
                    continue
                if not common.arm_is_final(round_dir, ds_dir.name, arm_dir.name, sample):
                    continue  # a twostep_v2 copy that has not finished still holds v1's results
                calls, has_col = read_calls(summary)
                legacy += not has_col
                runs[(ds_dir.name, sample, arm_of[arm_dir.name])] = calls
    return runs, legacy


def compare(round1: Path, round2: Path):
    arms2 = ["combined_off", "combined_artefact", "twostep"]
    r1, legacy1 = collect(round1, {BASELINE_ARM.get(a, a): a for a in arms2})
    r2, legacy2 = collect(round2, {a: a for a in arms2})
    rows = []
    for key in sorted(set(r1) | set(r2)):
        dataset, sample, arm = key
        calls1, calls2 = r1.get(key), r2.get(key)
        # Keep negative runs visible even when neither round has a detected virus.
        for virus in sorted(set(calls1 or {}) | set(calls2 or {})) or [""]:
            a, b = (calls1 or {}).get(virus), (calls2 or {}).get(virus)
            if calls1 is None:
                status = "run_only_round2"
            elif calls2 is None:
                status = "run_only_round1"
            elif not virus:
                status = "no_calls"
            else:
                status = "both" if a and b else "appeared" if b else "disappeared"
            row = {
                "group": "institutional" if dataset in INSTITUTIONAL else "public",
                "dataset": dataset,
                "sample": sample,
                "arm": arm,
                "virus": virus,
                "status": status,
                "delta_pct_infected_called": "",
            }
            for prefix, call in (("r1_", a), ("r2_", b)):
                row.update({prefix + f: (call or {}).get(f, "") for f in FIELDS})
            if a and b:
                try:
                    delta = float(b["pct_infected_called"]) - float(a["pct_infected_called"])
                    row["delta_pct_infected_called"] = f"{delta:.4f}"
                except ValueError:
                    pass
            rows.append(row)
    return rows, legacy1, legacy2


def _top(rows: list[dict], prefix: str) -> str:
    best = max(
        (r for r in rows if r[prefix + "viral_molecules_total_est"] != ""),
        key=lambda r: float(r[prefix + "viral_molecules_total_est"] or 0),
        default=None,
    )
    if best is None:
        return "-"
    return f"{best['virus']} {best[prefix + 'pct_infected_called']}%"


def _cells(row: dict, prefix: str) -> str:
    if row[prefix + "infected_called"] == "":
        return "-"
    return f"{row[prefix + 'infected_called']}/{row[prefix + 'n_called_cells']}"


def markdown(rows: list[dict], legacy1: int, legacy2: int, round1: Path, round2: Path) -> str:
    out = [
        "# DSR round 2 vs round 1",
        "",
        f"- round 1: `{round1}`  (baseline for `twostep` is `twostep_v2`)",
        f"- round 2: `{round2}`",
        "- Calls are the `detected` rows of `viral_summary.tsv` (AUD-21); `pct_infected_called` is each "
        "run's own called-cell denominator (use `dsr_common_cells.py` to compare arms with each other).",
        "- A blank virus in the TSV records a run with no detected calls; it is not a viral call.",
    ]
    if legacy1 or legacy2:
        out.append(
            f"- `viral_summary.tsv` without a `detected` column (read as all-detected): "
            f"{legacy1} in round 1, {legacy2} in round 2. Those hold candidate-only viruses too, so a "
            "`disappeared` call may only be a candidate that round 2 no longer lists."
        )
    for group in ("public", "institutional"):
        sub = [r for r in rows if r["group"] == group]
        out += ["", f"## {group.capitalize()} datasets", ""]
        if not sub:
            out.append("No runs.")
            continue
        out += [
            "| dataset | sample | arm | round 1 | round 2 | appeared | disappeared | top call r1 | top call r2 |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        by_run: dict[tuple, list[dict]] = defaultdict(list)
        for r in sub:
            by_run[(r["dataset"], r["sample"], r["arm"])].append(r)
        for (ds, sample, arm), rs in by_run.items():
            n1 = sum(
                bool(r["virus"]) and r["status"] in ("both", "disappeared", "run_only_round1")
                for r in rs
            )
            n2 = sum(
                bool(r["virus"]) and r["status"] in ("both", "appeared", "run_only_round2")
                for r in rs
            )
            app = sum(r["status"] == "appeared" for r in rs)
            dis = sum(r["status"] == "disappeared" for r in rs)
            only = {r["status"] for r in rs} & {"run_only_round1", "run_only_round2"}
            note = f" ({', '.join(sorted(only))})" if only else ""
            out.append(
                f"| {ds} | {sample} | {arm}{note} | {n1} | {n2} | {app} | {dis} | {_top(rs, 'r1_')} | {_top(rs, 'r2_')} |"
            )
        out += ["", "### Appeared and disappeared calls", ""]
        flagged = [r for r in sub if r["status"] in ("appeared", "disappeared")]
        if not flagged:
            out.append("None.")
        else:
            out += [
                "| dataset | sample | arm | virus | status | r1 infected/cells | r2 infected/cells | r1 molecules | r2 molecules |",
                "|---|---|---|---|---|---|---|---|---|",
            ]
            for r in flagged:
                out.append(
                    f"| {r['dataset']} | {r['sample']} | {r['arm']} | {r['virus']} | {r['status']} | "
                    f"{_cells(r, 'r1_')} | {_cells(r, 'r2_')} | {r['r1_viral_molecules_total_est'] or '-'} | "
                    f"{r['r2_viral_molecules_total_est'] or '-'} |"
                )
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round1", type=Path)
    ap.add_argument("round2", type=Path)
    ap.add_argument(
        "--out", type=Path, required=True, help="output prefix; writes PREFIX.tsv and PREFIX.md"
    )
    args = ap.parse_args(argv)
    rows, legacy1, legacy2 = compare(args.round1, args.round2)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with Path(f"{args.out}.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    Path(f"{args.out}.md").write_text(markdown(rows, legacy1, legacy2, args.round1, args.round2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

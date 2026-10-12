#!/usr/bin/env python3
"""VERDICT-01 calibration gate: do the verdict rules call known positives and known artefacts correctly?

    python scripts/dsr02_calibration.py <evidence_root> <labels.tsv> [--reference-cells <round>/reference_cells] \\
        [--out calibration.json] > calibration.tsv

`labels.tsv` (tab-separated, written by the user, never inferred) has `evidence_dir` (a directory name under
the evidence root), `expected` (`positive` or `artefact`) and an optional `note`. A positive passes when its
verdict is `viral_best`; an artefact passes when it is anything else. The rules and thresholds are those in
`scripts/dsr02_verdicts.py` (PROVISIONAL until the user freezes them); this script only reports. `--out` writes the
record that `scripts/dsr03_evidence_cards.py` reads: the rule values in force, a hash of the labels, per-row
results and `gate_passed`. The gate passes only when every labelled directory exists, there is at least one
positive and one artefact (a gate without both is vacuous), and every row passes. Changing a rule value
invalidates an old record (`gate_status`).

Choosing the labelled set and freezing the rules is the user's decision; nothing here picks either.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXPECTED = ("positive", "artefact")
OUT_COLS = [
    "evidence_dir",
    "expected",
    "verdict",
    "molecule_scope",
    "virus_best_molecules",
    "host_best_molecules",
    "tie_molecules",
    "passed",
    "reason",
]


def _verdicts():
    spec = importlib.util.spec_from_file_location(
        "dsr02_verdicts", REPO / "scripts" / "dsr02_verdicts.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def current_rule() -> dict[str, float]:
    """The rule values a calibration record is only valid for."""
    v = _verdicts()
    return {
        "min_virus_best_fraction": v.MIN_VIRUS_BEST_FRACTION,
        "min_complex_body_fraction": v.MIN_COMPLEX_BODY_FRACTION,
    }


def read_labels(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        if r.get("expected") not in EXPECTED:
            raise SystemExit(f"{path}: {r.get('evidence_dir')!r} has expected={r.get('expected')!r}; use {EXPECTED}")
    return rows


def evaluate(
    evidence_root: Path, labels: list[dict[str, str]], reference_cells: Path | None = None
) -> list[dict[str, object]]:
    v = _verdicts()
    out: list[dict[str, object]] = []
    for lab in labels:
        d = evidence_root / lab["evidence_dir"]
        row: dict[str, object] = {"evidence_dir": lab["evidence_dir"], "expected": lab["expected"]}
        if not (d / "interpretation_flags.tsv").is_file():
            out.append({**row, "verdict": "missing", "passed": False, "reason": "no evidence directory"})
            continue
        cells = v.frozen_cells(reference_cells, d.name) if reference_cells else None
        got = v.verdict(d, cells)
        if reference_cells and cells is None:
            got["molecule_scope"] = "no_frozen_cells"
        positive = lab["expected"] == "positive"
        passed = (got["verdict"] == "viral_best") == positive
        reason = "" if passed else (
            f"positive called {got['verdict']}" if positive else "artefact called viral_best"
        )
        out.append(
            {**row, **{k: got[k] for k in OUT_COLS if k in got and k not in row}, "passed": passed, "reason": reason}
        )
    return out


def record(rows: list[dict[str, object]], labels_path: Path) -> dict[str, object]:
    n_pos = sum(r["expected"] == "positive" for r in rows)
    n_art = sum(r["expected"] == "artefact" for r in rows)
    all_pass = bool(rows) and all(r["passed"] for r in rows)
    return {
        "rule": current_rule(),
        "labels_sha256": hashlib.sha256(labels_path.read_bytes()).hexdigest(),
        "n_positive": n_pos,
        "n_artefact": n_art,
        "all_rows_pass": all_pass,
        "gate_passed": all_pass and n_pos > 0 and n_art > 0,
        "rows": rows,
    }


def gate_status(record_path: Path | None) -> tuple[bool, str]:
    """(passed, reason) for a calibration record under the rule values in force now."""
    if record_path is None or not record_path.is_file():
        return False, "no calibration record"
    rec = json.loads(record_path.read_text())
    if rec.get("rule") != current_rule():
        return False, "calibration record was made under different rule values"
    if not rec.get("gate_passed"):
        return False, "calibration gate not passed"
    return True, ""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("evidence_root", type=Path)
    ap.add_argument("labels", type=Path)
    ap.add_argument("--reference-cells", type=Path)
    ap.add_argument("--out", type=Path, help="write the calibration record (JSON)")
    args = ap.parse_args(argv)
    rows = evaluate(args.evidence_root, read_labels(args.labels), args.reference_cells)
    w = csv.DictWriter(sys.stdout, fieldnames=OUT_COLS, delimiter="\t", lineterminator="\n", restval="")
    w.writeheader()
    w.writerows(rows)
    rec = record(rows, args.labels)
    if args.out:
        args.out.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
    print(f"gate_passed={rec['gate_passed']} ({rec['n_positive']} positive, {rec['n_artefact']} artefact)", file=sys.stderr)
    return 0 if rec["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

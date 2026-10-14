"""VERDICT-01: scripts/dsr02_calibration.py reports the gate and never picks labels or thresholds."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "dsr02_calibration", Path(__file__).resolve().parents[1] / "scripts" / "dsr02_calibration.py"
)
cal = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cal)

QC = "count_layer\treference\treference_class\treads\tcomplex_body_fraction\n"


def _ev(root, name, virus, host):
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "alignment_qc.tsv").write_text(QC + "deduplicated\tVIRUS|x\tvirus\t5\t0.9\n")
    (d / "interpretation_flags.tsv").write_text(
        "flag\tstatus\nhost_homology\tnot_flagged\nlow_complexity\tnot_flagged\n"
    )
    (d / "viral_reads.fasta").write_text(">r\nACGT\n")
    (d / "molecule_verdict_summary.tsv").write_text(
        "scope\tverdict\tmolecules\n"
        + "".join(
            f"counted\t{k}\t{n}\n"
            for k, n in (("virus_best", virus), ("host_best", host), ("tie", 0), ("unaligned", 0))
        )
    )


def _setup(tmp_path, labels):
    root = tmp_path / "ev"
    _ev(root, "ebv", 40, 1)  # called viral_best
    _ev(root, "telomere", 0, 8)  # called host_best
    _ev(root, "leaky", 30, 2)  # an artefact the rules would pass: the gate must fail
    path = tmp_path / "labels.tsv"
    path.write_text("evidence_dir\texpected\tnote\n" + "".join(f"{a}\t{b}\t\n" for a, b in labels))
    return root, path


def test_gate_passes_only_when_positives_and_artefacts_are_both_right(tmp_path):
    root, labels = _setup(tmp_path, [("ebv", "positive"), ("telomere", "artefact")])
    rows = cal.evaluate(root, cal.read_labels(labels))
    rec = cal.record(rows, labels)
    assert [r["passed"] for r in rows] == [True, True] and rec["gate_passed"]


def test_an_artefact_called_viral_best_fails_the_gate(tmp_path):
    root, labels = _setup(
        tmp_path, [("ebv", "positive"), ("telomere", "artefact"), ("leaky", "artefact")]
    )
    rec = cal.record(cal.evaluate(root, cal.read_labels(labels)), labels)
    assert not rec["gate_passed"] and not rec["all_rows_pass"]
    bad = [r for r in rec["rows"] if not r["passed"]]
    assert [r["evidence_dir"] for r in bad] == ["leaky"]
    assert bad[0]["reason"] == "artefact called viral_best"


def test_a_gate_without_both_classes_or_with_a_missing_directory_is_not_passed(tmp_path):
    root, labels = _setup(tmp_path, [("ebv", "positive")])
    assert not cal.record(cal.evaluate(root, cal.read_labels(labels)), labels)["gate_passed"]
    root, labels = _setup(tmp_path, [("ebv", "positive"), ("absent", "artefact")])
    rows = cal.evaluate(root, cal.read_labels(labels))
    assert rows[1]["verdict"] == "missing" and not cal.record(rows, labels)["gate_passed"]


def test_unknown_expected_value_is_refused(tmp_path):
    _, labels = _setup(tmp_path, [("ebv", "probably")])
    with pytest.raises(SystemExit):
        cal.read_labels(labels)


def test_record_is_only_valid_under_the_rule_values_that_made_it(tmp_path, monkeypatch):
    root, labels = _setup(tmp_path, [("ebv", "positive"), ("telomere", "artefact")])
    out = tmp_path / "cal.json"
    assert cal.main([str(root), str(labels), "--out", str(out)]) == 0
    assert cal.gate_status(out) == (True, "")
    rec = json.loads(out.read_text())
    rec["rule"]["min_virus_best_fraction"] = 0.9  # as if the threshold changed after calibration
    out.write_text(json.dumps(rec))
    ok, why = cal.gate_status(out)
    assert not ok and "different rule" in why
    assert cal.gate_status(None) == (False, "no calibration record")
    assert cal.gate_status(tmp_path / "nope.json")[0] is False

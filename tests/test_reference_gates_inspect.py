"""CAT-40 / CAT-19: record ``kallisto inspect`` on every build; fail on regression."""

from __future__ import annotations

import json
import random
import subprocess
from pathlib import Path

import pytest

from viralscan import reference_gates as gates
from viralscan.validation import tool_path

#: Real output of kallisto 0.52.0 on the 2,343-genome index (stderr+stdout, 2026-10-12).
RECORDED = """\
[index] k-mer length: 31
[index] number of targets: 471,948
[index] number of k-mers: 84,861,133
[index] number of D-list k-mers: 632,265
[inspect] Index version number = 13
[inspect] number of unitigs = 519982
[inspect] minimizer length = 23
[inspect] max EC size = 9921
[inspect] number of ECs discarded = 0
"""


def test_parse_recorded_text() -> None:
    got = gates.parse_inspect(RECORDED)
    assert got["max_ec_size"] == 9921
    assert got["discarded_ec_count"] == 0
    assert got["n_targets"] == 471948 and got["n_kmers"] == 84861133 and got["k"] == 31


def test_parse_rejects_text_without_gated_metrics() -> None:
    with pytest.raises(ValueError, match="discarded_ec_count"):
        gates.parse_inspect("[inspect] max EC size = 3\n")


def test_regression_gate() -> None:
    base = {"max_ec_size": 10, "discarded_ec_count": 0}
    assert gates.check_ec_regression(dict(base), base) == []
    assert gates.check_ec_regression({"max_ec_size": 9, "discarded_ec_count": 0}, base) == []
    both = gates.check_ec_regression({"max_ec_size": 11, "discarded_ec_count": 2}, base)
    assert len(both) == 2 and "10 -> 11" in both[0]
    assert gates.check_ec_regression({"max_ec_size": 10}, base)  # missing metric is a failure


def test_cli_from_text_records_and_gates(tmp_path: Path) -> None:
    text = tmp_path / "inspect.txt"
    text.write_text(RECORDED)
    out, base = tmp_path / "rec.json", tmp_path / "base.json"
    assert gates.main(["inspect", "--from-text", str(text), "--out", str(base)]) == 0
    assert json.loads(base.read_text())["max_ec_size"] == 9921
    assert (
        gates.main(
            ["inspect", "--from-text", str(text), "--out", str(out), "--baseline", str(base)]
        )
        == 0
    )
    worse = tmp_path / "worse.txt"
    worse.write_text(RECORDED.replace("9921", "9922"))
    assert (
        gates.main(
            ["inspect", "--from-text", str(worse), "--out", str(out), "--baseline", str(base)]
        )
        == 1
    )
    assert (
        gates.main(["inspect", "--from-text", str(text), "--out", str(out), "--max-ec-size", "5"])
        == 1
    )


def test_record_index_inspect_is_non_fatal(tmp_path: Path) -> None:
    index = tmp_path / "not_an_index.idx"
    index.write_bytes(b"x")
    assert gates.record_index_inspect(index) is None
    assert not (tmp_path / gates.INSPECT_RECORD).exists()


def _tiny_index(tmp_path: Path, shared: bool) -> Path:
    rng = random.Random(1)
    base = "".join(rng.choice("ACGT") for _ in range(600))
    other = base if shared else "".join(rng.choice("ACGT") for _ in range(600))
    fasta = tmp_path / "t.fa"
    fasta.write_text(f">A\n{base}\n>B\n{base}\n>C\n{other}\n")
    index = tmp_path / "t.idx"
    subprocess.run(
        [tool_path("kallisto"), "index", "-i", str(index), str(fasta)],
        check=True,
        capture_output=True,
    )
    return index


@pytest.mark.integration
@pytest.mark.skipif(tool_path("kallisto") is None, reason="kallisto not installed")
def test_tiny_index_is_recorded_and_a_larger_ec_regresses(tmp_path: Path) -> None:
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    tight = _tiny_index(tmp_path / "a", shared=False)
    record = gates.record_index_inspect(tight)
    assert record is not None and (tmp_path / "a" / gates.INSPECT_RECORD).is_file()
    assert record["max_ec_size"] >= 2  # A and B are identical, so one EC holds both
    loose = _tiny_index(tmp_path / "b", shared=True)
    now = gates.run_inspect(loose)
    assert now["max_ec_size"] > record["max_ec_size"]
    assert gates.check_ec_regression(now, record)
    assert gates.check_ec_regression(record, now) == []


def test_build_gates_record_and_fail_on_regression_or_unresolved_markers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from viralscan.gene_programs import load_catalogue

    t2g = tmp_path / "panel.t2g"
    t2g.write_text(
        "".join(f"t{i}\t{r['gene_id_bundled']}\n" for i, r in enumerate(load_catalogue()))
    )
    index = tmp_path / "panel.idx"
    now = {"max_ec_size": 12, "discarded_ec_count": 0}
    monkeypatch.setattr(gates, "run_inspect", lambda _index: dict(now))
    assert gates.run_build_gates(tmp_path, index, t2g) == []
    assert json.loads((tmp_path / gates.INSPECT_RECORD).read_text())["max_ec_size"] == 12
    assert (tmp_path / "programme_markers.json").is_file()
    baseline = tmp_path / "base.json"
    baseline.write_text(json.dumps({"max_ec_size": 11, "discarded_ec_count": 0}))
    assert any("regressed" in f for f in gates.run_build_gates(tmp_path, index, t2g, baseline))
    t2g.write_text("t0\tENSG1\n")
    assert gates.run_build_gates(tmp_path, index, t2g)  # markers no longer resolve


def test_scale_gate_writes_json(tmp_path: Path) -> None:
    gtf = tmp_path / "v.gtf"
    gtf.write_text('c\tt\texon\t1\t9\t.\t+\t.\tgene_id "G"; transcript_id "T";\n')
    line = gates.run_scale_gate(gtf, None, tmp_path / "scale_check.json")
    assert "1 genes" in line and (tmp_path / "scale_check.json").is_file()

"""DSR-16: explore ledger is append-only, hash-chained and deny-list aware."""

from __future__ import annotations

import pytest

from scripts import explore_ledger as el

REC = {
    "timestamp": "2026-10-13T00:00:00Z",
    "dataset": "ebv",
    "sample": "SRR12682296",
    "arm": "twostep",
    "code_sha": "abc1234",
    "reference_id": "cat42d",
    "label": "explore",
}


def _ledger(tmp_path, n=3):
    path = tmp_path / "explore_ledger.tsv"
    for i in range(n):
        el.append(path, {**REC, "sample": f"S{i}"}, [])
    return path


def test_append_and_check_clean(tmp_path):
    path = _ledger(tmp_path)
    assert el.check(path, ["GSE189536"]) == []
    assert len(path.read_text().splitlines()) == 4


def test_edited_row_is_detected(tmp_path):
    path = _ledger(tmp_path)
    path.write_text(path.read_text().replace("S1", "S9"))
    assert any("chain broken" in e for e in el.check(path, []))


def test_removed_middle_row_is_detected(tmp_path):
    path = _ledger(tmp_path)
    lines = path.read_text().splitlines()
    path.write_text("\n".join(lines[:2] + lines[3:]) + "\n")
    assert el.check(path, [])


def test_removed_tail_needs_expected_head(tmp_path):
    path = _ledger(tmp_path)
    head = el.read_rows(path)[-1][-1]
    path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")
    assert el.check(path, []) == []
    assert el.check(path, [], expect_head=head)


def test_deny_listed_accession_refused_and_flagged(tmp_path):
    with pytest.raises(ValueError, match="GSE256141"):
        el.append(tmp_path / "l.tsv", {**REC, "dataset": "x_gse256141"}, ["GSE256141"])
    path = _ledger(tmp_path, 1)
    prev = el.read_rows(path)[-1][-1]  # a row smuggled in without append(), chain recomputed
    fields = ["2026-10-13T00:00:01Z", "gse189536", "SRR1", "twostep", "s", "r", "explore"]
    with path.open("a") as fh:
        fh.write("\t".join([*fields, prev, el.row_hash(prev, fields)]) + "\n")
    assert any("deny-listed" in e for e in el.check(path, ["GSE189536"]))


def test_cli_exit_codes(tmp_path):
    led, deny = str(tmp_path / "l.tsv"), tmp_path / "deny.txt"
    deny.write_text("# c\nGSE298312\n")
    base = ["--ledger", led, "--deny-list", str(deny)]
    ok = [
        "append",
        "--dataset",
        "ebv",
        "--sample",
        "S",
        "--arm",
        "twostep",
        "--code-sha",
        "abc",
        "--reference-id",
        "cat42d",
    ]
    assert el.main([*base, *ok]) == 0
    assert el.main([*base, "check"]) == 0
    bad = [*ok]
    bad[2] = "gse298312"
    assert el.main([*base, *bad]) == 1


@pytest.mark.parametrize("separator", ["\t", "\r", "\n"])
def test_field_separator_refused_without_modifying_ledger(tmp_path, separator):
    path = _ledger(tmp_path, 1)
    before = path.read_bytes()
    with pytest.raises(ValueError, match="tabs or newlines"):
        el.append(path, {**REC, "sample": "S" + separator + "forged"}, [])
    assert path.read_bytes() == before


def test_missing_ledger_cli_check_fails(tmp_path):
    assert el.main(["--ledger", str(tmp_path / "missing.tsv"), "check"]) == 1


def test_truncated_newline_cannot_be_appended(tmp_path):
    path = _ledger(tmp_path, 1)
    path.write_bytes(path.read_bytes().rstrip(b"\n"))
    before = path.read_bytes()
    with pytest.raises(ValueError, match="final newline"):
        el.append(path, REC, [])
    assert path.read_bytes() == before
    with pytest.raises(ValueError, match="final newline"):
        el.check(path, [])

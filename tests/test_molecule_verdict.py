"""VERDICT-01: per-molecule host-vs-virus verdicts from a competitive SAM."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

from viralscan import molecule_verdict as mv


def _sam(qname: str, flag: int, ref: str, score: int | None) -> str:
    tags = [] if score is None else [f"AS:i:{score}"]
    return "\t".join(
        [qname, str(flag), ref, "1", "60", "50M", "*", "0", "0", "A" * 50, "I" * 50, *tags]
    )


def _verdicts(lines: list[str], counted: set[int] | None = None) -> dict[tuple[str, str], str]:
    return {(r["cb"], r["ub"]): r["verdict"] for r in mv.molecule_rows(lines, counted or set())}


def test_best_side_wins_and_equal_scores_tie() -> None:
    got = _verdicts(
        [
            _sam("AAA_U1_1|r", 0, "VIRUS|x", 90),  # virus beats the host secondary
            _sam("AAA_U1_1|r", 256, "HOST|chr1", 60),
            _sam("AAA_U2_2|r", 0, "HOST|chr1", 80),  # host beats the viral secondary
            _sam("AAA_U2_2|r", 256, "VIRUS|x", 70),
            _sam("AAA_U3_3|r", 0, "HOST|chr1", 75),  # equal: the primary is a coin flip, so a tie
            _sam("AAA_U3_3|r", 256, "VIRUS|x", 75),
            _sam("AAA_U4_4|r", 0, "VIRUS|x", 40),  # the host fell under minimap2's secondary cutoff
        ]
    )
    assert got == {
        ("AAA", "U1"): "virus_best",
        ("AAA", "U2"): "host_best",
        ("AAA", "U3"): "tie",
        ("AAA", "U4"): "virus_best",
    }


def test_unmapped_and_supplementary_alignments_do_not_score() -> None:
    got = _verdicts(
        [
            _sam("AAA_U1_1|r", 4, "*", None),
            _sam("AAA_U2_2|r", 0, "VIRUS|x", 50),
            _sam("AAA_U2_2|r", 2048, "HOST|chr1", 99),  # supplementary: a partial hit, ignored
        ]
    )
    assert got == {("AAA", "U2"): "virus_best"}  # the unmapped read leaves no molecule


def test_a_molecule_takes_its_best_score_over_its_reads() -> None:
    got = _verdicts(
        [
            _sam("AAA_U1_1|r", 0, "HOST|chr1", 60),
            _sam("AAA_U1_2|r", 0, "VIRUS|x", 80),  # a second read of the same UMI
        ]
    )
    assert got == {("AAA", "U1"): "virus_best"}


def test_summary_lists_every_verdict_and_separates_counted(tmp_path: Path) -> None:
    lineage = tmp_path / "read_lineage.tsv.gz"
    with gzip.open(lineage, "wt", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(["read_number", "read_id", "cb", "ub", "assigned_weight"])
        w.writerow([1, "r", "AAA", "U1", "1.0"])
        w.writerow([2, "r", "AAA", "U2", "0.0"])  # host/virus-ambiguous: not counted by the call
        w.writerow([3, "r", "AAA", "U3", ""])
    sam = [
        _sam("AAA_U1_1|r", 0, "VIRUS|x", 90),
        _sam("AAA_U2_2|r", 0, "HOST|chr1", 80),
        _sam("AAA_U3_3|r", 0, "HOST|chr1", 80),
    ]
    summary = mv.write_molecule_verdicts(sam, lineage, tmp_path)
    got = {(r["scope"], r["verdict"]): r["molecules"] for r in summary}
    assert len(summary) == 8
    assert got["all", "virus_best"] == 1 and got["all", "host_best"] == 2
    assert got["counted", "virus_best"] == 1 and got["counted", "host_best"] == 0
    with gzip.open(tmp_path / "molecule_verdicts.tsv.gz", "rt") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 3
    assert (
        (tmp_path / "molecule_verdict_summary.tsv")
        .read_text()
        .startswith("scope\tverdict\tmolecules")
    )

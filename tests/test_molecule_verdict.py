"""VERDICT-01: per-molecule host-vs-virus verdicts from a competitive SAM."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

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


def test_post_alignment_selection_keeps_competition_and_unmapped_molecules() -> None:
    lines = [_sam(f"CB_U1_{i}|r", 0, "VIRUS|x", 90) for i in range(1, 6)]
    lines += [_sam("CB_U1_1|r", 256, "HOST|chr1", 80), _sam("CB_U2_6|r", 4, "*", None)]
    selected = mv.deduplicate_competitive_alignments(lines)
    assert {line.split("\t")[0] for line in selected} == {"CB_U1_1|r", "CB_U2_6|r"}
    assert len(selected) == 3
    assert selected == mv.deduplicate_competitive_alignments(reversed(lines))[::-1]
    with pytest.raises(ValueError, match="Missing corrected CB/UMI"):
        mv.deduplicate_competitive_alignments([_sam("anonymous", 0, "VIRUS|x", 90)])


def test_query_fasta_contains_two_post_alignment_molecules(tmp_path: Path) -> None:
    lines = [_sam(f"CB_U1_{i}|r", 0, "VIRUS|x", 90) for i in range(1, 6)]
    lines += [_sam("CB_U2_6|r", 0, "VIRUS|x", 90)]
    source, output = tmp_path / "raw.fa", tmp_path / "umi.fa"
    source.write_text(
        "".join(f">CB_{'U1' if i < 6 else 'U2'}_{i}|r\n{'A' * 50}\n" for i in range(1, 7))
    )
    mv.write_umi_query_fasta(lines, source, output)
    assert output.read_text().count(">") == 2


def test_representative_is_the_read_with_the_best_viral_alignment() -> None:
    lines = [
        _sam("CB_U1_1|r", 0, "HOST|chr1", 100),  # a host-primary read of the same molecule
        _sam("CB_U1_2|r", 0, "VIRUS|x", 60),
        _sam("CB_U1_3|r", 0, "VIRUS|x", 90),
    ]
    reps, _by_read, n_reads = mv.select_representatives(lines)
    assert reps == {("CB", "U1"): "CB_U1_3|r"} and n_reads[("CB", "U1")] == 3
    assert mv.select_representatives(reversed(lines))[0] == reps  # input order never decides

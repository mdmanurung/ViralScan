"""EVID-02: ``viralscan evidence --splice`` on a synthetic viral genome with one two-exon minus-strand gene."""

import csv
import gzip
import json
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from viralscan import splice

needs_tools = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and samtools are required",
)

# 1-based layout: exon 1 = 101..200, intron = 201..400, exon 2 = 401..500; the gene is on the minus strand,
# so the intron reads CT...AC on the genome strand (GT...AG in transcript orientation).
EXON1, INTRON, EXON2 = (101, 200), (201, 400), (401, 500)


def _genome() -> str:
    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(600)]
    seq[200:202] = "CT"  # intron start, 1-based 201-202
    seq[398:400] = "AC"  # intron end, 1-based 399-400
    return "".join(seq)


GENOME = _genome()
SPLICED = GENOME[EXON1[0] - 1 : EXON1[1]] + GENOME[EXON2[0] - 1 : EXON2[1]]  # genome orientation


def _read(left: int, right: int) -> str:
    """Genome-orientation read: *left* bases of exon 1 then *right* bases of exon 2."""
    return SPLICED[100 - left : 100 + right]


def _write_gz(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _sam_record(
    number: int, cb: str, ub: str, flag: int, seq: str, rname: str = "VIRUS|chrV"
) -> str:
    stored = splice.revcomp(seq) if flag & 0x10 else seq
    cigar = "30M" if seq == "*" else f"{len(seq)}M"
    return "\t".join(
        [f"{cb}_{ub}_{number}|read{number}", str(flag), rname, "1", "60", cigar, "*", "0", "0", stored,
         "*", "AS:i:60"]
    )  # fmt: skip


#: (molecule, counted, verdict, [(read_number, flag, sequence)])
MOLECULES = [
    (("CB1", "U1"), True, "virus_best", [(1, 0, _read(30, 30)), (2, 16, _read(35, 25))]),
    (("CB1", "U2"), True, "virus_best", [(3, 0, _read(30, 30)), (4, 0, GENOME[110:170])]),
    (("CB2", "U3"), True, "host_best", [(5, 0, _read(30, 30))]),  # spliced-looking, but host won
    (("CB2", "U4"), False, "virus_best", [(6, 0, _read(30, 30))]),  # not counted
    (("CB2", "U5"), True, "tie", [(7, 0, _read(30, 30))]),
    (("CB2", "U6"), True, "unaligned", [(8, 4, GENOME[110:170])]),
    (("CB3", "U7"), True, "virus_best", [(9, 0, _read(4, 56))]),  # 4 bp anchor: must not be spliced
]


def _evidence_dir(tmp_path: Path, molecules: list = MOLECULES) -> Path:
    ev = tmp_path / "evidence"
    ev.mkdir()
    _write_gz(
        ev / "molecule_verdicts.tsv.gz",
        ["cb", "ub", "reads", "counted", "virus_as", "host_as", "verdict"],
        [[cb, ub, len(r), counted, "", "", verdict] for (cb, ub), counted, verdict, r in molecules],
    )
    _write_gz(
        ev / "read_lineage.tsv.gz",
        ["read_number", "read_id", "cb", "ub", "assigned_weight"],
        [[n, f"read{n}", cb, ub, 1] for (cb, ub), _, _, r in molecules for n, _, _ in r],
    )
    sam = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:VIRUS|chrV\tLN:600", "@SQ\tSN:HOST|chr1\tLN:1000"]
    for (cb, ub), _, _, reads in molecules:
        sam += [_sam_record(n, cb, ub, flag, seq) for n, flag, seq in reads]
    # a hard-clipped supplementary record of read 1 and a secondary host record of read 2
    sam.append(_sam_record(1, "CB1", "U1", 0x800, _read(10, 10)))
    sam.append(_sam_record(2, "CB1", "U1", 0x100, "*", "HOST|chr1"))
    sam_path = ev / "reads.sam"
    sam_path.write_text("\n".join(sam) + "\n")
    subprocess.run(
        ["samtools", "view", "-b", "-o", str(ev / "competitive_reads.raw.bam"), str(sam_path)],
        check=True,
    )
    sam_path.unlink()
    return ev


def _fasta(tmp_path: Path) -> str:
    path = tmp_path / "viral.fa"
    path.write_text(f">chrV some virus\n{GENOME}\n>chrOther\n{'ACGT' * 50}\n")
    return str(path)


def _gtf(tmp_path: Path, feature: str) -> Path:
    path = tmp_path / f"{feature}.gtf"
    rows = [
        f'chrV\ttest\t{feature}\t{start}\t{end}\t.\t-\t.\tgene_id "G1"; transcript_id "T1";'
        for start, end in (EXON1, EXON2)
    ]
    path.write_text("\n".join(rows) + "\n")
    return path


def _table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@needs_tools
def test_spliced_molecules_junction_and_exclusions(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path)
    row = splice.run_splice(
        ev, viral_fasta=_fasta(tmp_path), annotation_gtfs=[_gtf(tmp_path, "exon")]
    )
    out = ev / "splice"
    summary = _table(out / "splice_summary.tsv")[0]
    assert summary == {k: str(v) for k, v in row.items()}
    # the five categories of docs/evidence_unit.md
    assert [summary[k] for k in splice.CATEGORY_FIELDS] == ["3", "1", "1", "1", "1"]
    assert summary["status"] == "assessed"
    assert (
        summary["reads_realigned"] == "5"
    )  # M1 (2) + M2 (2) + M7 (1); host-best, tie, uncounted never realigned
    assert summary["spliced_reads"] == "3"
    assert summary["spliced_molecules"] == "2"
    assert summary["spliced_molecules_annotated"] == "2"
    assert summary["unique_junctions"] == "1"
    assert summary["annotation_status"] == "annotated_introns_available"
    assert "does not establish" in summary["note"]
    junctions = _table(out / "junctions.tsv")
    assert junctions == [
        {
            "contig": "chrV",
            "start": "201",
            "end": "400",
            "strand": "-",
            "motif": "GT-AG",
            "motif_canonical": "True",
            "supporting_molecules": "2",  # M1 has two spliced reads and counts once
            "supporting_reads": "3",
            "annotated": "True",
            "annotated_transcript": "T1",
            "counted_in_spliced_molecules": "True",
        }
    ]
    assert summary["unique_junctions_excluded_noncanonical"] == "0"
    assert summary["molecules_excluded_noncanonical"] == "0"
    assert (out / "splice.bam.bai").is_file()


@needs_tools
def test_receipt_hashes_inputs_outputs_and_records_the_command(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path)
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path), annotation_gtfs=[_gtf(tmp_path, "exon")])
    manifest = json.loads((ev / "splice" / "splice_manifest.json").read_text())
    assert set(manifest["inputs"]) >= {
        "molecule_verdicts.tsv.gz",
        "read_lineage.tsv.gz",
        "competitive_reads.raw.bam",
    }
    assert {"splice.bam", "junctions.tsv", "splice_summary.tsv", "target_genome.fa"} <= set(
        manifest["outputs"]
    )
    assert manifest["minimap2"]["options"] == ["-ax", "splice:sr", "--secondary=no", "-G", "100000"]
    assert manifest["minimap2"]["version"]
    assert manifest["constants"]["MIN_ANCHOR"] == splice.MIN_ANCHOR
    assert manifest["minimap2"]["target_contigs"] == ["chrV"]  # chrOther was not hit


@needs_tools
def test_cds_only_annotation_is_not_evaluable_and_leaves_annotated_empty(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path)
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path), annotation_gtfs=[_gtf(tmp_path, "CDS")])
    summary = _table(ev / "splice" / "splice_summary.tsv")[0]
    assert summary["annotation_status"] == "intron_matching_not_evaluable"
    assert summary["spliced_molecules"] == "2"
    assert summary["spliced_molecules_annotated"] == ""
    junction = _table(ev / "splice" / "junctions.tsv")[0]
    assert junction["annotated"] == "" and junction["annotated_transcript"] == ""


@needs_tools
def test_no_annotation_file_is_not_evaluable(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path)
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path))
    assert (
        _table(ev / "splice" / "splice_summary.tsv")[0]["annotation_status"]
        == "intron_matching_not_evaluable"
    )


def test_competitor_mode_is_not_assessed_and_computes_nothing(tmp_path: Path) -> None:
    ev = tmp_path / "evidence"
    ev.mkdir()
    splice.run_splice(ev, viral_fasta=None, competitor_mode=True)
    summary = _table(ev / "splice" / "splice_summary.tsv")[0]
    assert summary["status"] == "not_assessed"
    assert summary["reason"] == "molecule_verdict_not_produced_in_competitor_mode"
    assert all(summary[k] == "" for k in splice.CATEGORY_FIELDS)
    assert _table(ev / "splice" / "junctions.tsv") == []
    assert not (ev / "splice" / "splice.bam").exists()
    assert (
        json.loads((ev / "splice" / "splice_manifest.json").read_text())["status"] == "not_assessed"
    )


@needs_tools
def test_no_counted_virus_best_molecule_is_not_evaluable_with_counts(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path, [m for m in MOLECULES if m[2] != "virus_best" or not m[1]])
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path), annotation_gtfs=[_gtf(tmp_path, "exon")])
    summary = _table(ev / "splice" / "splice_summary.tsv")[0]
    assert summary["status"] == "not_evaluable"
    assert summary["reason"] == "no_counted_virus_best_molecules"
    assert [summary[k] for k in splice.CATEGORY_FIELDS] == ["0", "1", "1", "1", "1"]
    assert not (ev / "splice" / "splice.bam").exists()


@needs_tools
def test_empty_verdicts_file(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path, [])
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path))
    summary = _table(ev / "splice" / "splice_summary.tsv")[0]
    assert summary["status"] == "not_evaluable"
    assert [summary[k] for k in splice.CATEGORY_FIELDS] == ["0"] * 5


def test_missing_verdicts_is_not_assessed(tmp_path: Path) -> None:
    ev = tmp_path / "evidence"
    ev.mkdir()
    splice.run_splice(ev, viral_fasta=None)
    assert _table(ev / "splice" / "splice_summary.tsv")[0]["reason"] == "molecule_verdicts_missing"


@needs_tools
def test_member_read_without_a_primary_record_is_refused(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path)
    with gzip.open(ev / "read_lineage.tsv.gz", "at", newline="") as handle:
        handle.write("99\tread99\tCB1\tU1\t1\n")
    with pytest.raises(splice.SpliceError, match="no primary record"):
        splice.run_splice(ev, viral_fasta=_fasta(tmp_path))


# --- version check (no aligner needed) --------------------------------------------------------


@pytest.mark.parametrize(
    "text,ok", [("2.28-r1209", False), ("2.29-r1283", True), ("2.31-r1307", True), ("3.0", True)]
)
def test_version_gate(monkeypatch: pytest.MonkeyPatch, text: str, ok: bool) -> None:
    monkeypatch.setattr(splice, "minimap2_version", lambda: text)
    if ok:
        assert splice.require_splice_minimap2() == text
    else:
        with pytest.raises(splice.SpliceError, match=r"minimap2 >= 2\.29.*no fallback"):
            splice.require_splice_minimap2()


def test_old_minimap2_stops_run_splice_before_any_alignment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    if shutil.which("samtools") is None:
        pytest.skip("samtools is required")
    ev = _evidence_dir(tmp_path)
    monkeypatch.setattr(splice, "minimap2_version", lambda: "2.24-r1122")
    with pytest.raises(splice.SpliceError, match="2.24"):
        splice.run_splice(ev, viral_fasta=_fasta(tmp_path))
    assert not (ev / "splice" / "splice.bam").exists()


def test_unreadable_version_is_an_error() -> None:
    with pytest.raises(splice.SpliceError):
        splice.parse_minimap2_version("not a version")


# --- junction rules on single records (no aligner needed) --------------------------------------


def test_short_anchor_is_refused_on_either_side() -> None:
    assert splice.junctions_of(0, 100, "5M200N55M", 0, 60) == ([], {"anchor": 1})
    assert splice.junctions_of(0, 100, "55M200N7M", 0, 62) == ([], {"anchor": 1})
    accepted, _ = splice.junctions_of(0, 100, "8M200N52M", 0, 60)
    assert accepted == [(108, 307)]


def test_anchor_is_measured_between_neighbouring_introns() -> None:
    accepted, rejected = splice.junctions_of(0, 1, "30M100N5M100N30M", 0, 65)
    assert accepted == [] and rejected == {"anchor": 2}


def test_supplementary_secondary_and_unmapped_records_report_nothing() -> None:
    for flag in (0x800, 0x100, 0x4):
        assert splice.junctions_of(flag, 1, "30M200N30M", 0, 60) == ([], {})


def test_intron_length_bounds() -> None:
    assert splice.junctions_of(0, 1, "30M19N30M", 0, 60)[1] == {"intron_length": 1}
    assert splice.junctions_of(0, 1, "30M20N30M", 0, 60)[0] == [(31, 50)]
    assert splice.junctions_of(0, 1, "30M100001N30M", 0, 60)[1] == {"intron_length": 1}


def test_mismatch_and_clip_caps_refuse_the_whole_record() -> None:
    assert splice.junctions_of(0, 1, "30M200N30M", 5, 60) == ([], {"mismatch_fraction": 1})
    assert splice.junctions_of(0, 1, "30M200N30M", 4, 60)[0] == [(31, 230)]
    assert splice.junctions_of(0, 1, "20S20M200N20M", 0, 60) == ([], {"clip_fraction": 1})


def test_deletion_and_insertion_do_not_shift_the_intron() -> None:
    assert splice.junctions_of(0, 10, "20M2D10M200N30M", 1, 60)[0] == [(42, 241)]
    assert splice.junctions_of(0, 10, "20M2I10M200N28M", 1, 60)[0] == [(40, 239)]


def test_motif_in_transcript_orientation() -> None:
    assert splice.motif_of(GENOME, 201, 400) == ("GT-AG", True, "-")
    plus = "AA" + "GT" + "C" * 30 + "AG" + "TT"  # intron 3..36 on the forward strand
    assert splice.motif_of(plus, 3, 36) == ("GT-AG", True, "+")
    assert splice.motif_of("AAGGCCCCCCTTAA", 3, 12) == ("GG-TT", False, ".")


# --- annotated introns --------------------------------------------------------------------------


def test_annotated_introns_use_exon_gaps_of_one_transcript_only(tmp_path: Path) -> None:
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        "\n".join(
            f'{c}\tt\t{f}\t{s}\t{e}\t.\t+\t.\tgene_id "G"; transcript_id "{t}";'
            for c, f, s, e, t in [
                ("c1", "exon", 1, 10, "T1"),
                ("c1", "exon", 51, 60, "T1"),
                ("c1", "exon", 1, 10, "T2"),  # T1 and T2 exons are never paired with each other
                ("c1", "exon", 51, 60, "T2"),
                ("c1", "exon", 11, 20, "T3"),  # abutting exons: no intron
                ("c1", "exon", 21, 30, "T3"),
                ("c2", "CDS", 1, 10, "T4"),
                ("c2", "CDS", 51, 60, "T4"),  # CDS gap: never an intron
            ]
        )
        + "\n"
    )
    introns, evaluable = splice.annotated_introns([gtf], ["c1", "c2", "c3"])
    assert introns == {"c1": {(11, 50): ["T1", "T2"]}}
    assert evaluable == {"c1"}


# --- wiring into the evidence verb ----------------------------------------------------------------


def test_splice_flag_is_off_by_default() -> None:
    from viralscan.cli.parsers import parse_args

    base = ["evidence", "--run-dir", "r", "--virus", "EBV"]
    assert parse_args(base).splice is False
    assert parse_args([*base, "--splice"]).splice is True


def test_splice_without_viral_fasta_is_refused(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    import argparse

    from viralscan.scripts.evidence_run import run_evidence

    args = argparse.Namespace(
        run_dir=str(tmp_path), virus="EBV", splice=True, viral_fasta=None, review_candidates=False
    )
    with pytest.raises(SystemExit):
        run_evidence(args)
    assert "--splice requires --viral-fasta" in caplog.text


@needs_tools
def test_evidence_manifest_lists_the_splice_outputs(tmp_path: Path) -> None:
    from viralscan.scripts.evidence_run import _write_evidence_manifest

    ev = _evidence_dir(tmp_path)
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path))
    replay = {"strand": None, "technology": "10xv3", "kallisto_technology": "10XV3"}
    path = _write_evidence_manifest(
        ev,
        target_label="T",
        target_genes=["g"],
        method="m",
        run_fingerprint="f",
        references={},
        replay=replay,
    )
    outputs = json.loads(path.read_text())["outputs"]
    for name in ("splice.bam", "junctions.tsv", "splice_summary.tsv", "splice_manifest.json"):
        assert f"splice/{name}" in outputs


@needs_tools
def test_annotation_match_is_exact_and_canonical_unannotated_still_counts(tmp_path: Path) -> None:
    """EVID-DECISIONS-01: a gap shifted by 1 bp is not annotated; a canonical motif still counts."""
    ev = _evidence_dir(tmp_path)
    shifted = tmp_path / "shifted.gtf"
    shifted.write_text(
        "\n".join(
            f'chrV\ttest\texon\t{start}\t{end}\t.\t-\t.\tgene_id "G1"; transcript_id "T1";'
            for start, end in ((EXON1[0], EXON1[1] + 1), EXON2)
        )
        + "\n"
    )
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path), annotation_gtfs=[shifted])
    out = ev / "splice"
    (junction,) = _table(out / "junctions.tsv")
    assert junction["annotated"] == "False"
    summary = _table(out / "splice_summary.tsv")[0]
    assert summary["spliced_molecules_annotated"] == "0"
    assert summary["spliced_molecules"] == "2"  # GT-AG on the minus strand is canonical
    assert junction["counted_in_spliced_molecules"] == "True"


@needs_tools
def test_noncanonical_unannotated_junction_is_excluded_but_counted_separately(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ev = _evidence_dir(tmp_path)
    monkeypatch.setattr(splice, "motif_of", lambda genome, start, end: ("AA-CC", False, "."))
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path), annotation_gtfs=[_gtf(tmp_path, "intron")])
    summary = _table(ev / "splice" / "splice_summary.tsv")[0]
    assert summary["spliced_molecules"] == "0"
    assert summary["unique_junctions"] == "1"
    assert summary["unique_junctions_excluded_noncanonical"] == "1"
    assert summary["molecules_excluded_noncanonical"] == "2"
    assert _table(ev / "splice" / "junctions.tsv")[0]["counted_in_spliced_molecules"] == "False"


@needs_tools
def test_noncanonical_annotated_junction_still_counts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ev = _evidence_dir(tmp_path)
    monkeypatch.setattr(splice, "motif_of", lambda genome, start, end: ("AA-CC", False, "."))
    splice.run_splice(ev, viral_fasta=_fasta(tmp_path), annotation_gtfs=[_gtf(tmp_path, "exon")])
    summary = _table(ev / "splice" / "splice_summary.tsv")[0]
    assert summary["spliced_molecules"] == "2"
    assert summary["spliced_molecules_annotated"] == "2"
    assert summary["unique_junctions_excluded_noncanonical"] == "0"

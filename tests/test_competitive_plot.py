"""ALIGN-PLOT-01: read-by-read host-vs-virus figure and the cross-call locus table."""

from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

from viralscan import competitive_plot as cp
from viralscan.menu import build_parser

pytestmark = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")

HEADER = "@HD\tVN:1.6\n@SQ\tSN:HOST|chr7\tLN:159345973\n@SQ\tSN:HOST|chr15\tLN:101991189\n@SQ\tSN:HOST|chr2\tLN:242193529\n@SQ\tSN:VIRUS|KP1.1\tLN:3800\n"
TSO_READ = (cp.TSO + "GGTAAACTGGGTCTTTGTCATTTGCTTTAAAAAAAAAAAAAAAT" + "A" * 30)[:89]


def _name(i: int, lib: str = "SRR1") -> str:
    return f"CB{i:04d}_UMI{i:04d}_{i}|{lib}.{i}"


def _sam(
    name: str, flag: int, ref: str, pos: int, cigar: str, seq: str, score: int, extra: str = ""
) -> str:
    rec = [name, str(flag), ref, str(pos), "10", cigar, "*", "0", "0", seq, "*", f"AS:i:{score}"]
    return "\t".join(rec + ([extra] if extra else []))


def _evidence_dir(
    root: Path,
    name: str,
    sam_lines: list[str],
    verdicts: list[tuple[str, int, str, str, str]],
    fasta: dict[str, str] | None = None,
) -> Path:
    """verdicts: (read, length, virus_as, host_as, repeat_class); verdict is derived the way the pipeline does."""
    d = root / name
    d.mkdir()
    sam = d / "x.sam"
    sam.write_text(HEADER + "\n".join(sam_lines) + ("\n" if sam_lines else ""))
    subprocess.run(
        ["samtools", "view", "-b", "-o", str(d / "competitive_reads.raw.bam"), str(sam)], check=True
    )
    sam.unlink()
    with gzip.open(d / "read_verdicts.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "read",
                "cb",
                "umi",
                "read_number",
                "length",
                "virus_as",
                "host_as",
                "verdict",
                "repeat_class",
                "counted",
            ]
        )
        for read, length, v, h, rc in verdicts:
            va, ha = (int(v) if v else None), (int(h) if h else None)
            from viralscan.molecule_verdict import verdict_of

            w.writerow([read, "", "", "", length, v, h, verdict_of(va, ha), rc, 1])
    if fasta:
        (d / "viral_reads.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in fasta.items()))
    return d


@pytest.fixture
def evidence(tmp_path: Path) -> Path:
    tt, tso, host, un = _name(1), _name(2), _name(3), _name(4)
    seq = "ACGT" * 22 + "AC"  # 90 nt filler
    tt_seq = seq[:86] + "A" * 12
    host_seq = seq + "A" * 8  # 98 nt
    rc = host_seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    sam = [
        _sam(tt, 0, "VIRUS|KP1.1", 2798, "12S31M55S", tt_seq, 62),
        _sam(tso, 0, "HOST|chr7", 157386233, "29S60M", TSO_READ, 98),
        # reverse-strand primary plus a hard-clipped supplementary on another chromosome
        _sam(
            host,
            16,
            "HOST|chr7",
            157386214,
            "36S62M",
            rc,
            124,
            "SA:Z:HOST|chr15,54926032,-,50M48S,1,0;",
        ),
        _sam(host, 2064, "HOST|chr15", 54926032, "50M48H", rc[:50], 100),
        _sam(un, 4, "*", 0, "*", seq, 0).replace("\tAS:i:0", ""),
    ]
    verdicts = [
        (tt, 98, "62", "", "none"),
        (tso, 89, "", "98", "tso"),
        (host, 98, "", "124", "none"),
        (un, 90, "", "", "none"),
    ]
    return _evidence_dir(
        tmp_path,
        "ev_a",
        sam,
        verdicts,
        {tt: tt_seq, tso: TSO_READ, host: host_seq},
    )


def test_reverse_strand_segment_is_flipped_into_read_orientation() -> None:
    seg = cp._segment("HOST|chr7", 100, True, "50M48H", "supplementary", 1)
    assert seg is not None
    assert (seg.q0, seg.q1, seg.read_length, seg.strand) == (48, 98, 98, "-")
    fwd = cp._segment("VIRUS|v", 5, False, "12S31M55S", "primary", 62)
    assert fwd is not None and (fwd.q0, fwd.q1, fwd.side, fwd.label) == (12, 43, "virus", "v:5")


def test_poly_runs_bridge_one_base_and_need_ten_nt() -> None:
    lane = cp.ReadLane(
        "r",
        60,
        None,
        None,
        "unaligned",
        "none",
        sequence="GG" + "A" * 15 + "T" + "A" * 8 + "CC" + "T" * 9 + "GG" + "T" * 12,
    )
    assert lane.poly_runs() == [(2, 26, "A"), (39, 51, "T")]


def test_lanes_carry_every_alignment_and_the_tso(evidence: Path) -> None:
    lanes = {ln.name: ln for ln in cp.load_lanes(evidence)}
    host = lanes[_name(3)]
    assert sorted((s.ref, s.kind) for s in host.segments) == [
        ("HOST|chr15", "supplementary"),
        ("HOST|chr7", "primary"),
    ]
    assert lanes[_name(2)].tso_span() == (0, len(cp.TSO))
    assert lanes[_name(4)].segments == [] and lanes[_name(4)].verdict == "unaligned"
    assert lanes[_name(1)].poly_runs()[0][2] == "A"


def test_sa_tag_alone_yields_the_supplementary_segment(tmp_path: Path) -> None:
    n = _name(1)
    seq = "ACGT" * 25
    d = _evidence_dir(
        tmp_path,
        "sa",
        [_sam(n, 0, "HOST|chr7", 100, "60M40S", seq, 90, "SA:Z:HOST|chr15,500,-,40M60S,1,0;")],
        [(n, 100, "", "90", "none")],
    )
    (lane,) = cp.load_lanes(d)
    assert [(s.ref, s.kind, s.q0, s.q1) for s in lane.segments] == [
        ("HOST|chr7", "primary", 0, 60),
        ("HOST|chr15", "supplementary", 60, 100),
    ]


def test_figure_written_and_evidence_dir_untouched(evidence: Path, tmp_path: Path) -> None:
    before = sorted(p.name for p in evidence.iterdir())
    out = tmp_path / "plots"
    paths = cp.plot_competitive_alignment(evidence, out)
    assert [Path(p).name for p in paths] == [
        "competitive_alignment.png",
        "competitive_alignment.pdf",
    ]
    assert all(Path(p).stat().st_size > 1000 for p in paths)
    assert sorted(p.name for p in evidence.iterdir()) == before


def test_default_output_is_the_evidence_dir(evidence: Path) -> None:
    cp.plot_competitive_alignment(evidence)
    assert (evidence / "competitive_alignment.png").is_file()


def test_zero_reads_and_unaligned_only(tmp_path: Path) -> None:
    empty = _evidence_dir(tmp_path, "empty", [], [])
    assert cp.load_lanes(empty) == []
    assert (tmp_path / "e").exists() is False
    cp.plot_competitive_alignment(empty, tmp_path / "o1")
    n = _name(1)
    unal = _evidence_dir(
        tmp_path,
        "unal",
        [_sam(n, 4, "*", 0, "*", "ACGT" * 22 + "AC", 0)],
        [(n, 90, "", "", "none")],
        {n: "ACGT" * 22},
    )
    (lane,) = cp.load_lanes(unal)
    assert lane.verdict == "unaligned" and lane.length == 90
    cp.plot_competitive_alignment(unal, tmp_path / "o2")
    assert (tmp_path / "o1" / "competitive_alignment.png").is_file() and (
        tmp_path / "o2" / "competitive_alignment.pdf"
    ).is_file()


def test_verdicts_derived_from_the_bam_when_the_table_is_missing(evidence: Path) -> None:
    (evidence / "read_verdicts.tsv.gz").unlink()
    lanes = {ln.name: ln for ln in cp.load_lanes(evidence)}
    assert (lanes[_name(1)].verdict, lanes[_name(2)].verdict) == ("virus_best", "host_best")
    assert lanes[_name(2)].repeat_class == "tso"


def test_above_the_lane_cap_the_lane_block_is_replaced_by_the_length_histogram(
    tmp_path: Path,
) -> None:
    sam, verdicts = [], []
    for i in range(1, 251):
        n = _name(i)
        sam.append(_sam(n, 0, "VIRUS|KP1.1", 100 + i, "50M", "ACGTA" * 10, 80))
        verdicts.append((n, 50, "80", "", "none"))
    d = _evidence_dir(tmp_path, "big", sam, verdicts)
    lanes = cp.load_lanes(d)
    assert len(lanes) == 250 and {ln.viral_aligned() for ln in lanes} == {50}
    cp.plot_competitive_alignment(
        d, tmp_path / "o"
    )  # 250 molecules > 60: no lane block, histogram panel D
    cp.plot_competitive_alignment(d, tmp_path / "raw", layer="raw")
    assert (tmp_path / "o" / "competitive_alignment.png").is_file()
    assert (tmp_path / "raw" / "competitive_alignment.png").is_file()


def test_cross_call_ranks_the_recurrent_locus_first(tmp_path: Path) -> None:
    dirs = []
    for k, lib in enumerate(("SRR1", "SRR2", "SRR3")):
        n = _name(k + 1, lib)
        dirs.append(
            _evidence_dir(
                tmp_path,
                f"call{k}",
                [_sam(n, 0, "HOST|chr7", 157386200 + k * 9, "60M", "ACGTA" * 12, 90)],
                [(n, 60, "", "90", "none")],
            )
        )
    n = _name(9, "SRR9")
    dirs.append(
        _evidence_dir(
            tmp_path,
            "call3",
            [_sam(n, 0, "HOST|chr2", 5000, "60M", "ACGTA" * 12, 90)],
            [(n, 60, "", "90", "none")],
        )
    )
    rows = cp.cross_call_loci(dirs, tmp_path / "out" / "cross_call")
    assert (
        rows[0]["reference"],
        rows[0]["bin_start"],
        rows[0]["n_libraries"],
        rows[0]["n_calls"],
    ) == ("chr7", 157386001, 3, 3)
    with (tmp_path / "out" / "cross_call.tsv").open() as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    assert list(table[0]) == list(cp.CROSS_CALL_FIELDS) and len(table) == 2
    assert (tmp_path / "out" / "cross_call.png").is_file()


def test_cli_subcommand_writes_figures_and_cross_call(evidence: Path, tmp_path: Path) -> None:
    args = build_parser().parse_args(
        [
            "evidence-plot",
            "--evidence-dir",
            str(evidence),
            "-o",
            str(tmp_path / "p"),
            "--cross-call",
            str(tmp_path / "p" / "x"),
        ]
    )
    assert args._subcommand == "evidence-plot"
    cp.run_evidence_plot(args)
    assert (tmp_path / "p" / "ev_a" / "competitive_alignment.png").is_file() and (
        tmp_path / "p" / "x.tsv"
    ).is_file()


def test_missing_inputs_raise(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        cp.load_lanes(tmp_path)


def test_duplicate_reads_count_as_two_molecules_and_retain_host_competition(tmp_path: Path) -> None:
    names = [f"CB_U1_{i}|SRR1.{i}" for i in range(1, 6)] + ["CB_U2_6|SRR1.6"]
    seq = "ACGTA" * 12
    sam = [_sam(n, 0, "VIRUS|KP1.1", 100, "60M", seq, 90) for n in names]
    sam += [_sam(n, 256, "HOST|chr7", 1000, "60M", "*", 80) for n in names]
    d = _evidence_dir(tmp_path, "duplicates", sam, [(n, 60, "90", "80", "none") for n in names])
    lanes = cp.load_lanes(d)
    assert len(lanes) == 2 and len(cp.load_lanes(d, "raw")) == 6
    assert all({s.side for s in lane.segments} == {"host", "virus"} for lane in lanes)
    rows = cp.cross_call_loci([d], tmp_path / "dedup")
    assert all(r["count_layer"] == "umi" and r["n_molecules"] == 2 for r in rows)


def test_cross_call_deduplicates_molecules_even_when_representative_reads_differ(
    tmp_path: Path,
) -> None:
    dirs = []
    for i in (1, 2):
        n = f"CB_U1_{i}|SRR1.{i}"
        dirs.append(
            _evidence_dir(
                tmp_path,
                f"call{i}",
                [_sam(n, 0, "HOST|chr7", 100, "60M", "ACGTA" * 12, 90)],
                [(n, 60, "", "90", "none")],
            )
        )
    (row,) = cp.cross_call_loci(dirs, tmp_path / "shared")
    assert row["n_molecules"] == 1 and row["n_calls"] == 2

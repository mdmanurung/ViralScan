"""Live samtools coordinate parity for opt-in coverage maps."""

import json
import shutil
import subprocess

import pandas as pd
import pytest

from viralscan.evidence import CoverageMapError, coverage_map_figure, resolve_reference_annotation


@pytest.fixture
def evidence_dir(tmp_path):
    if not shutil.which("samtools"):
        pytest.skip("samtools required for live BAM fixture")
    d = tmp_path / "evidence"
    d.mkdir()
    sam = d / "reads.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:TARGET|X\tLN:1000\n@SQ\tSN:chr1\tLN:100000000\n"
        + "".join(
            f"r{i}\t0\tTARGET|X\t{p}\t60\t1M\t*\t0\t0\tA\tI\n" for i, p in enumerate([5, 5, 900])
        )
    )
    for name in ["raw", "umi_dedup"]:
        subprocess.run(
            ["samtools", "view", "-b", "-o", str(d / f"competitive_reads.{name}.bam"), str(sam)],
            check=True,
        )
    return d


def table(d, target="TARGET|X"):
    return pd.read_csv(d / f"coverage_map.{target}.tsv.gz", sep="\t")


def test_gap_coordinate_parity_missing_annotation(evidence_dir):
    d = evidence_dir
    meta = coverage_map_figure(d)
    assert meta["status"] == "annotation_missing"
    data = table(d)
    assert len(data) == 1000 and data.raw_depth.iloc[4] == 2 and data.raw_depth.iloc[899] == 1
    assert data.raw_depth.iloc[5:899].eq(0).all()
    depth = subprocess.check_output(
        ["samtools", "depth", str(d / "competitive_reads.raw.bam")], text=True
    )
    for line in depth.splitlines():
        ref, position, value = line.split("\t")
        assert data.raw_depth.iloc[int(position) - 1] == int(value)
    assert {p.suffix for p in d.glob("coverage_map.*")} >= {".png", ".pdf", ".gz", ".json"}


def test_annotation_mismatch_and_missing_target(evidence_dir):
    d = evidence_dir
    gtf = d / "genes.gtf"
    gtf.write_text('TARGET|X\tt\tgene\t1\t1001\t.\t+\t.\tgene_id "G";\n')
    with pytest.raises(CoverageMapError, match="1001.*1000") as error:
        coverage_map_figure(d, annotation_path=gtf)
    assert error.value.status == "annotation_mismatch"
    assert (
        json.loads((d / "coverage_map.TARGET|X.json").read_text())["status"]
        == "annotation_mismatch"
    )
    with pytest.raises(CoverageMapError, match="available references") as error:
        coverage_map_figure(d, "absent")
    assert error.value.status == "target_not_in_bam"


def test_multiple_targets_require_selection(evidence_dir, monkeypatch):
    import viralscan.evidence as ev

    monkeypatch.setattr(ev, "reference_lengths", lambda _: {"TARGET|A": 100, "TARGET|B": 100})
    with pytest.raises(CoverageMapError, match="--target"):
        coverage_map_figure(evidence_dir)


def test_zero_coverage(evidence_dir):
    d = evidence_dir
    meta = coverage_map_figure(d, "chr1")
    assert meta["status"] == "no_coverage" and meta["bin_size"] == 2000
    data = table(d, "chr1")
    assert len(data) == 50000 and data.raw_depth.eq(0).all()


def test_transcript_and_circular_annotation(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('NM_123.1\tt\tgene\t1\t9000\t.\t+\t.\tgene_id "G";\n')
    assert resolve_reference_annotation(gtf, "NM_123.1", 1000) == ([], "transcript_reference")
    gtf.write_text('V\tt\texon\t1\t1000\t.\t+\t.\tgene_id "G";\n')
    assert resolve_reference_annotation(gtf, "V", 1000) == ([], "transcript_reference")
    gtf.write_text('V\tt\tCDS\t900\t100\t.\t-\t.\tgene_id "G";\n')
    features, status = resolve_reference_annotation(gtf, "V", 1000)
    assert status == "ok" and features[0][:3] == (900, 100, "-")


def test_gene_track_and_cli(evidence_dir):
    gtf = evidence_dir / "genes.gtf"
    gtf.write_text('TARGET|X\tt\tgene\t5\t900\t.\t+\t.\tgene_id "G";\n')
    assert coverage_map_figure(evidence_dir, annotation_path=gtf)["status"] == "ok"


def test_coverage_map_cli_dispatch(evidence_dir):
    from viralscan.competitive_plot import run_evidence_plot
    from viralscan.menu import build_parser

    args = build_parser().parse_args(
        ["evidence-plot", "--coverage-map", "--evidence-dir", str(evidence_dir)]
    )
    run_evidence_plot(args)
    assert table(evidence_dir).shape[0] == 1000


def test_transcript_map_suppresses_genome_overlay(evidence_dir, monkeypatch):
    import viralscan.evidence as ev

    monkeypatch.setattr(ev, "reference_lengths", lambda _: {"NM_123.1": 1000})
    gtf = evidence_dir / "genome.gtf"
    gtf.write_text('NM_123.1\tt\tgene\t1\t9000\t.\t+\t.\tgene_id "G";\n')
    original_depth = ev.coverage_depth_points
    monkeypatch.setattr(
        ev,
        "coverage_depth_points",
        lambda bam: [dict(r, reference="NM_123.1") for r in original_depth(bam)],
    )
    metadata = coverage_map_figure(evidence_dir, "NM_123.1", gtf)
    assert metadata["status"] == metadata["annotation_status"] == "transcript_reference"


def _rename_target(d, new):
    # Rewrite the BAM reference name so the fixture can mimic other reference classes.
    for name in ["raw", "umi_dedup"]:
        bam = d / f"competitive_reads.{name}.bam"
        sam = subprocess.check_output(["samtools", "view", "-h", str(bam)], text=True)
        sam = sam.replace("TARGET|X", new)
        subprocess.run(
            ["samtools", "view", "-b", "-o", str(bam), "-"], input=sam, text=True, check=True
        )


def test_viral_prefixed_reference_accepts_bare_gtf_seqid(evidence_dir):
    d = evidence_dir
    gtf = d / "genes.gtf"
    gtf.write_text('X\tt\tgene\t5\t900\t.\t+\t.\tgene_id "G";\n')
    meta = coverage_map_figure(d, annotation_path=gtf)
    assert meta["status"] == "ok"
    saved = json.loads((d / "coverage_map.TARGET|X.json").read_text())
    assert saved["annotation_seqid"] == "X" and saved["annotation_seqid_normalised"] is True


def test_exact_seqid_not_marked_normalised(evidence_dir):
    d = evidence_dir
    gtf = d / "genes.gtf"
    gtf.write_text('TARGET|X\tt\tgene\t5\t900\t.\t+\t.\tgene_id "G";\n')
    coverage_map_figure(d, annotation_path=gtf)
    saved = json.loads((d / "coverage_map.TARGET|X.json").read_text())
    assert saved["annotation_seqid"] == "TARGET|X" and saved["annotation_seqid_normalised"] is False


def test_host_prefixed_reference_is_not_prefix_stripped(evidence_dir):
    d = evidence_dir
    _rename_target(d, "HOST|X")
    gtf = d / "genes.gtf"
    gtf.write_text('X\tt\tgene\t5\t900\t.\t+\t.\tgene_id "G";\n')
    with pytest.raises(CoverageMapError) as error:
        coverage_map_figure(d, "HOST|X", annotation_path=gtf)
    assert error.value.status == "annotation_mismatch"

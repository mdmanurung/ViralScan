"""Native Cutadapt checks for the exploratory prototype, not package-default qualification."""

import gzip
import importlib.util
import json
from pathlib import Path

import pytest

from tests._fastq import write_fastq
from viralscan.scripts import read_filter as rf

spec = importlib.util.spec_from_file_location(
    "cleanup_pilot", Path(__file__).parents[1] / "analysis/viral_read_cleanup/pilot.py"
)
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)

BODY = "CGACGTGCTACGATCGTAGCTAGCATCGATGCTAGCTGACTACGATCGAGCTACGATCG"
R1 = "ACGTACGTACGTACGTAAACCCGGGTTT"


@pytest.mark.parametrize("orientation", ["sense", "antisense"])
def test_native_cleanup_preserves_body_barcode_and_qualities(tmp_path, orientation):
    pytest.importorskip("cutadapt")
    body = BODY if orientation == "sense" else pilot._revcomp(BODY)
    tail_read = body + "A" * 25 if orientation == "sense" else "T" * 25 + body
    r1, r2 = tmp_path / "input1.fastq.gz", tmp_path / "input2.fastq.gz"
    write_fastq(r1, [("tail/1", R1), ("internal/1", R1), ("short/1", R1)])
    internal = BODY[:30] + "A" * 20 + BODY[30:]
    short = "A" * 90 if orientation == "sense" else "T" * 90
    write_fastq(
        r2, [("tail/2", rf.TSO_FULL + tail_read), ("internal/2", internal), ("short/2", short)]
    )
    out = tmp_path / "out"
    pilot.transform(r1, r2, out, "cleanup", orientation)
    result = list(pilot.pairs(out / "R1.fastq.gz", out / "R2.fastq.gz"))
    assert [a[1].strip() for a, _ in result] == [R1, R1]
    assert [b[1].strip() for _, b in result] == [body, internal]
    assert json.loads((out / "audit.json").read_text())["removed"] == 1
    with gzip.open(out / "fragment_lineage.tsv.gz", "rt") as f:
        assert "short\tremoved\tshort_after_trim" in f.read()


def test_off_preserves_short_reads_and_refuses_overwrite(tmp_path):
    r1, r2 = tmp_path / "a.fastq", tmp_path / "b.fastq"
    write_fastq(r1, [("short", R1)])
    write_fastq(r2, [("short", "ACG")])
    out = tmp_path / "out"
    pilot.transform(r1, r2, out, "off")
    assert list(pilot.pairs(out / "R1.fastq.gz", out / "R2.fastq.gz"))[0][1][1].strip() == "ACG"
    with pytest.raises(FileExistsError):
        pilot.transform(r1, r2, out, "off")


@pytest.mark.parametrize("fault", ["mismatch", "duplicate", "truncated", "quality"])
def test_invalid_fastq_has_no_success_receipt(tmp_path, fault):
    r1, r2 = tmp_path / "a.fastq", tmp_path / "b.fastq"
    a = [("one", R1), ("one" if fault == "duplicate" else "two", R1)]
    b = [("wrong" if fault == "mismatch" else "one", BODY), (a[1][0], BODY)]
    write_fastq(r1, a)
    write_fastq(r2, b)
    if fault == "truncated":
        r2.write_text(r2.read_text().rsplit("\n", 2)[0] + "\n")
    if fault == "quality":
        lines = r2.read_text().splitlines()
        lines[3] = "I"
        r2.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        pilot.transform(r1, r2, tmp_path / "out", "off")
    assert not (tmp_path / "out/audit.json").exists()


def test_sample_rejects_unsupported_barcode_geometry(tmp_path):
    r1, r2 = tmp_path / "a.fastq", tmp_path / "b.fastq"
    write_fastq(r1, [("one", R1)])
    write_fastq(r2, [("one", BODY)])
    with pytest.raises(ValueError, match="geometry"):
        pilot.sample(r1, r2, tmp_path / "out", 1, 26)


def test_candidate_adjudication_keeps_host_and_sibling_ambiguity(tmp_path):
    out = tmp_path / "sample/off"
    quant = out / "quant"
    quant.mkdir(parents=True)
    write_fastq(out / "R1.fastq.gz", [(str(i), R1) for i in range(3)])
    write_fastq(out / "R2.fastq.gz", [(str(i), BODY) for i in range(3)])
    (quant / "transcripts.txt").write_text("target\nhost\nsibling\n")
    (quant / "matrix.ec").write_text("0\t0\n1\t0,1\n2\t0,2\n")
    (quant / "bus.tsv").write_text(
        "ACGT\tAAAA\t0\t1\t0\nACGT\tAAAC\t1\t1\t1\nACGT\tAAAG\t2\t1\t2\n"
    )
    t2g, viral = tmp_path / "t2g.txt", tmp_path / "viral.txt"
    t2g.write_text("target\tgene\nhost\thostgene\nsibling\tothergene\n")
    viral.write_text("target\tgene\t\t\tNC_target\nsibling\tothergene\t\t\tNC_sibling\n")
    cfg = tmp_path / "config.json"
    cfg.write_text(
        json.dumps(
            {
                "output": str(tmp_path),
                "samples": {"sample": {"accessions": ["NC_target"]}},
                "t2g": str(t2g),
                "viral_t2g": str(viral),
            }
        )
    )
    pilot.assess(cfg, "sample", "off")
    summary = json.loads((out / "candidate_summary.json").read_text())
    assert summary["read_classes"] == {
        "unique_target": 1,
        "host_virus_ambiguous": 1,
        "virus_virus_ambiguous": 1,
    }
    assert summary["unique_target_raw_cb_umi_keys"] == 1
    assert len(list(pilot.records(out / "candidate_reads.fastq"))) == 1

"""scripts/compare_reference_reproducibility.py (PLAN REF-04): identical vs differing audits."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import compare_reference_reproducibility as cmp  # noqa: E402

from viralscan.scripts import build_reference as builder  # noqa: E402


def _audit(tmp_path, name, fasta_text):
    out = tmp_path / name
    out.mkdir()
    fasta = out / "input.fa"
    fasta.write_text(fasta_text)
    gtf = out / "input.gtf"
    gtf.write_text('V1\ts\texon\t1\t10\t.\t+\t.\tgene_id "V1_g"; transcript_id "V1_tx";\n')
    manifest = builder.write_reference_manifest(
        fasta,
        out / "reference_manifest.json",
        profile="curated",
        host_species="none",
        viral_identifiers={"V1"},
    )
    builder.record_reference_build(manifest, fasta=fasta, gtf=gtf)
    return out / "reference_reproducibility.json"


def test_identical_rebuilds_exit_zero_and_a_changed_fasta_is_reported_per_field(tmp_path, capsys):
    a = _audit(tmp_path, "a", ">V1\nACGTCAGTAC\n")
    b = _audit(tmp_path, "b", ">V1\nACGTCAGTAC\n")  # different paths and dates, same content
    assert cmp.main([str(a), str(b)]) == 0

    c = _audit(tmp_path, "c", ">V1\nACGTCAGTAA\n")
    assert cmp.main([str(a), str(c)]) == 1
    printed = capsys.readouterr().out
    assert "content_sha256:" in printed
    assert "content_files.fasta:" in printed
    assert "content_files.gtf" not in printed


def test_a_missing_file_entry_is_a_difference():
    base = {"content_sha256": "x", "content_files": {"fasta": "1", "gtf": "2"}}
    other = {"content_sha256": "x", "content_files": {"fasta": "1"}}
    assert cmp.compare_reproducibility(base, other) == ["content_files.gtf: 2 != None"]
    assert cmp.compare_reproducibility(base, json.loads(json.dumps(base))) == []

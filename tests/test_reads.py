"""viralscan.reads: FASTA loading, lineage loading and the extract_virus_reads dispatch."""

import gzip
import sys

import pytest

from viralscan import reads


def _fake_run_evidence(args):
    from pathlib import Path

    out = Path(args.output)
    out.mkdir(parents=True)
    (out / "viral_reads.fasta").write_text(">AAAC_TTTT_1\nACGT\nAC\n>AAAG_TTTA_2 extra\nGGGG\n")
    with gzip.open(out / "read_lineage.tsv.gz", "wt") as fh:
        fh.write("read\tcell\n1\tAAAC\n")


def test_extract_and_load(tmp_path, monkeypatch):
    import viralscan.scripts.evidence_run as er

    seen = []
    monkeypatch.setattr(er, "run_evidence", lambda a: (seen.append(a), _fake_run_evidence(a)))
    got = reads.extract_virus_reads("run", "Human papillomavirus 16", tmp_path, cores=2)
    vr = got["Human papillomavirus 16"]
    assert (
        vr.out_dir.name == "Human_papillomavirus_16" and seen[0].virus == "Human papillomavirus 16"
    )
    assert (
        seen[0].viral_fasta is None and seen[0].cores == 2
    )  # extraction only: no alignment, no BLAST
    assert vr.sequences() == {"AAAC_TTTT_1": "ACGTAC", "AAAG_TTTA_2": "GGGG"}
    assert vr.lineage() == [{"read": "1", "cell": "AAAC"}]


def test_failure_raises_not_exits(tmp_path, monkeypatch):
    import viralscan.scripts.evidence_run as er

    monkeypatch.setattr(er, "run_evidence", lambda a: sys.exit(1))
    with pytest.raises(RuntimeError, match="read extraction failed"):
        reads.extract_virus_reads("run", ["x"], tmp_path)

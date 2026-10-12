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


def test_colliding_selectors_use_distinct_stable_outputs(tmp_path, monkeypatch):
    import viralscan.scripts.evidence_run as er

    seen = []
    monkeypatch.setattr(er, "run_evidence", lambda a: (seen.append(a), _fake_run_evidence(a)))
    names = ["X-1", "X 1", "Unrelated"]
    forward = reads.extract_virus_reads("run", iter(names), tmp_path / "forward")
    backward = reads.extract_virus_reads("run", reversed(names), tmp_path / "backward")
    assert len({vr.out_dir for vr in forward.values()}) == 3
    assert forward["Unrelated"].out_dir.name == "Unrelated"
    assert {n: vr.out_dir.name for n, vr in forward.items()} == {
        n: vr.out_dir.name for n, vr in backward.items()
    }
    assert all(vr.fasta.exists() for vr in forward.values())
    assert len(seen) == 6


def test_duplicate_selector_extracts_once(tmp_path, monkeypatch):
    import viralscan.scripts.evidence_run as er

    seen = []
    monkeypatch.setattr(er, "run_evidence", lambda a: (seen.append(a), _fake_run_evidence(a)))
    reads.extract_virus_reads("run", ["X", "X"], tmp_path)
    assert len(seen) == 1


def test_empty_output_slug_fails_before_extraction(tmp_path, monkeypatch):
    import viralscan.scripts.evidence_run as er

    seen = []
    monkeypatch.setattr(er, "run_evidence", seen.append)
    with pytest.raises(ValueError, match="filesystem-safe name"):
        reads.extract_virus_reads("run", ["X", "--"], tmp_path)
    assert not seen


def test_disambiguated_name_collision_fails_before_extraction(tmp_path, monkeypatch):
    import hashlib

    import viralscan.scripts.evidence_run as er

    seen = []
    monkeypatch.setattr(er, "run_evidence", seen.append)
    reserved = "X_1_" + hashlib.sha256(b"X-1").hexdigest()[:8]
    with pytest.raises(ValueError, match="colliding output names"):
        reads.extract_virus_reads("run", ["X-1", "X 1", reserved], tmp_path)
    assert not seen

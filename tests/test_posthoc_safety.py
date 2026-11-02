"""REVIEW-FIX-E: post-hoc verbs never leave a dishonest or destroyed tree."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from viralscan import menu
from viralscan.cli.parsers import parse_args
from viralscan.run_safety import (
    RunSafetyError,
    _fingerprint,
    _write_manifest_atomic,
    require_fresh_output,
    write_run_complete,
)


def make_run(root: Path, sample: str = "S1") -> Path:
    """A finished synthetic v3 run; returns the sample dir."""
    root.mkdir(parents=True)
    manifest = {
        "schema_version": "3.0.0",
        "viralscan_version": "x",
        "allocation_method": "uniform",
        "options": {"multimapping": True, "umap": True, "gene_programs": False},
        "input_fingerprints": {},
        "completion_marker": True,
        "layout_version": 3,
    }
    manifest["run_fingerprint"] = _fingerprint(manifest)
    _write_manifest_atomic(root, manifest)
    out = root / sample
    argv = ["-i", "i", "-t", "t", "--r1", "a_R1.fastq", "--r2", "b_R2.fastq", "-o", "o"]
    cfg = menu._build_run_config(parse_args(argv), str(out) + "/", "i", "t", None, "a", "b")
    menu._write_run_config(cfg)
    files = {
        "results/viral_summary.tsv": "virus\tn\nEBV\t3\n",
        "results/virus_identity.tsv": "a\tb\n",
        "results/analysis_genes.txt": "g1\n",
        "results/found_genes.txt": "g1;3\n",
        "summary.txt": "Headline\n\nRuntime: 1.0 seconds.\n\nAmount of transcripts in data: 5\n",
        "report.html": "<html/>",
        "log/multimap.done": "",
        "log/detection.done": "",
        "log/umap.done": "",
        "kb-python/counts_unfiltered/adata_multimap.h5ad": "x",
        "kb-python/counts_unfiltered/adata.h5ad": "x",
    }
    for rel, text in files.items():
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        (out / rel).write_text(text)
    write_run_complete(root)
    return out


def rerun_args(root: Path, output: Path, **over):
    base = dict(
        run_dir=str(root), output=str(output), multimap_method="em-global", cores=1,
        multimap_em_max_iter=None, multimap_em_tol=None, verbose=False, quiet=True, force=False,
    )  # fmt: skip
    return SimpleNamespace(**{**base, **over})


@pytest.fixture
def stubbed(monkeypatch):
    monkeypatch.setattr(menu, "_check_required_tools", lambda: None)
    monkeypatch.setattr(menu, "_check_cell_caller_tools", lambda *a, **k: None)


def test_guard_refuses_primary_tree_even_with_force(tmp_path):
    root = tmp_path / "run"
    sample = make_run(root)
    for bad in (sample / "results" / "newdir", sample / "newdir", sample / "results", root / "x"):
        with pytest.raises(RunSafetyError, match="inside the ViralScan run"):
            require_fresh_output(bad, force=True)
    require_fresh_output(root / "posthoc" / "x")
    require_fresh_output(sample / "posthoc" / "x")
    require_fresh_output(tmp_path / "elsewhere")


def test_em_tol_zero_rejected_by_parser():
    argv = ["rerun-multimap", "--run-dir", "r", "--multimap-method", "em-global"]
    with pytest.raises(SystemExit):
        parse_args([*argv, "--multimap-em-tol", "0"])


def test_bad_override_dies_before_copy(tmp_path, stubbed):
    root = tmp_path / "run"
    make_run(root)
    out = tmp_path / "out"
    with pytest.raises(SystemExit) as exc:
        menu._run_rerun_multimap(rerun_args(root, out, multimap_em_tol=0.0))
    assert exc.value.code == 2 and not out.exists()


def test_failed_rerun_leaves_unfinished_tree(tmp_path, stubbed, monkeypatch):
    root = tmp_path / "run"
    make_run(root)
    out = tmp_path / "out"

    def boom(*a, **k):
        raise SystemExit(1)

    monkeypatch.setattr(menu, "_run_snakemake", boom)
    with pytest.raises(SystemExit):
        menu._run_rerun_multimap(rerun_args(root, out))
    assert "NOT FINISHED" in (out / "S1" / "README.txt").read_text()
    assert "NOT FINISHED" in (out / "README.txt").read_text()
    assert not (out / "S1" / "results" / "artifact_inventory.json").exists()
    assert not (out / "combined_summary.tsv").exists()
    assert not (out / "run_complete.json").exists()


def test_v2_source_loses_legacy_files(tmp_path, stubbed, monkeypatch):
    root = tmp_path / "run"
    s = make_run(root)
    (s / "results" / "analysis_genes.txt").rename(s / "log" / "analysis.txt")
    (s / "results" / "found_genes.txt").rename(s / "log" / "found_genes.txt")
    (s / "plots").mkdir()
    (s / "plots" / "SuperExpressor_EBV.tsv").write_text("old\n")
    monkeypatch.setattr(menu, "_run_snakemake", lambda *a, **k: None)
    out = tmp_path / "out"
    menu._run_rerun_multimap(rerun_args(root, out))
    o = out / "S1"
    assert (o / "results" / "analysis_genes.txt").exists()
    assert not (o / "log" / "analysis.txt").exists()
    assert not (o / "log" / "found_genes.txt").exists()
    assert not list((o / "plots").glob("SuperExpressor_*"))


def _hostresponse_args(root: Path, out: Path, force: bool):
    argv = ["hostresponse", "--run-dir", str(root / "S1"), "--host-h5ad", "h.h5ad"]
    args = parse_args([*argv, "--output", str(out)])
    args.force, args.quiet = force, True
    return args


def test_hostresponse_force_failure_restores_old_outputs(tmp_path, monkeypatch):
    root = tmp_path / "run"
    make_run(root)
    out = tmp_path / "hr"
    out.mkdir()
    (out / "EBV_gene_weights.csv").write_text("precious\n")
    (out / "notes_mine.csv").write_text("mine\n")
    from viralscan.scripts import hostresponse as hr

    def clear_then_fail(**kw):
        hr.clear_stale_virus_outputs(kw["out_dir"])
        raise RuntimeError("model blew up")

    monkeypatch.setattr(hr, "run_hostresponse", clear_then_fail)
    with pytest.raises(RuntimeError):
        menu._run_hostresponse_subcommand(_hostresponse_args(root, out, force=True))
    assert (out / "EBV_gene_weights.csv").read_text() == "precious\n"
    assert (out / "notes_mine.csv").exists()
    assert not (tmp_path / ".hr.prev").exists()


def test_hostresponse_too_few_barcodes_exits_nonzero(tmp_path, monkeypatch):
    root = tmp_path / "run"
    make_run(root)
    from viralscan.scripts import hostresponse as hr

    monkeypatch.setattr(hr, "run_hostresponse", lambda **kw: {"status": "too_few_shared_barcodes"})
    with pytest.raises(SystemExit) as exc:
        menu._run_hostresponse_subcommand(_hostresponse_args(root, tmp_path / "hr", force=False))
    assert exc.value.code == 1


def test_auto_evidence_drops_demoted_genus(tmp_path, monkeypatch):
    from viralscan.scripts import auto_evidence as ae

    ev = tmp_path / "results" / "evidence"
    for genus in ("Alpha", "Beta"):
        (ev / genus).mkdir(parents=True)
        (ev / genus / "evidence_manifest.json").write_text(json.dumps({}))
    (ev / "review_candidates").mkdir()
    (ev / "review_candidates" / "evidence_review.tsv").write_text("x")
    monkeypatch.setattr(ae, "load_run_identity", lambda d: object())
    monkeypatch.setattr(ae, "anello_genera", lambda s, i: ["Alpha"])
    monkeypatch.setattr(
        ae,
        "extract_virus_reads",
        lambda run_dir, genera, out, **kw: {
            g: SimpleNamespace(out_dir=Path(out) / g) for g in genera
        },
    )
    cfg = SimpleNamespace(
        output=str(tmp_path), review_candidates=False, cores=1, viral_fasta=None, host_fasta=None
    )
    ae.run(SimpleNamespace(config=cfg), str(tmp_path / "done"))
    assert (ev / "Alpha").is_dir() and not (ev / "Beta").exists()
    assert (ev / "review_candidates").is_dir()

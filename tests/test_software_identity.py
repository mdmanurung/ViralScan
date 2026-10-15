"""Runtime provenance describes code without changing resume compatibility."""

import argparse
import json

import pytest

from viralscan import run_safety
from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.scripts.analysis import _chemistry_sanity


def test_identity_change_does_not_change_fingerprint_or_resume(tmp_path, monkeypatch):
    args = argparse.Namespace(multimap_method="equal", cell_calling="none")
    monkeypatch.setattr(run_safety, "software_identity", lambda: {"git_commit": "a" * 40})
    previous = run_safety.build_run_manifest(args)
    run_safety.prepare_output_directory(tmp_path, previous, resume=False, overwrite=False, yes=True)
    monkeypatch.setattr(run_safety, "software_identity", lambda: {"git_commit": "b" * 40})
    current = run_safety.build_run_manifest(args)
    assert current["software_identity"] != previous["software_identity"]
    assert current["run_fingerprint"] == previous["run_fingerprint"]
    assert (
        run_safety.prepare_output_directory(
            tmp_path, current, resume=True, overwrite=False, yes=True
        )
        == "resume"
    )


def test_installed_identity_hashes_relative_runtime_bytes(tmp_path, monkeypatch):
    package = tmp_path / "site-packages" / "viralscan"
    package.mkdir(parents=True)
    source = package / "run_safety.py"
    source.write_text("first")
    monkeypatch.setattr(run_safety, "__file__", str(source))
    first = run_safety.software_identity()
    assert first["kind"] == "installed" and "git_commit" not in first
    (package / "local.gtf").write_text("ignored reference input")
    assert run_safety.software_identity() == first
    source.write_text("second")
    assert run_safety.software_identity()["build_sha256"] != first["build_sha256"]


@pytest.mark.parametrize("p,failed", [(8.2, True), (25.0, False), (73.4, False)])
def test_chemistry_gate_retains_receipt_and_manifest_before_failure(tmp_path, p, failed):
    sample = tmp_path / "sample"
    kb = sample / "kb-python"
    kb.mkdir(parents=True)
    (kb / "run_info.json").write_text(json.dumps({"p_pseudoaligned": p}))
    (tmp_path / "run_manifest.json").write_text('{"run_fingerprint": "unchanged"}')
    ctx = RunContext.from_config(RunConfig(output=str(sample), require_chemistry_sanity=True))
    if failed:
        with pytest.raises(RuntimeError, match="chemistry_sanity.json"):
            _chemistry_sanity(ctx)
    else:
        _chemistry_sanity(ctx)
    report = json.loads((sample / "results" / "chemistry_sanity.json").read_text())
    manifest = json.loads((tmp_path / "run_manifest.json").read_text())
    assert report["p_pseudoaligned"] == p
    assert manifest["chemistry_sanity"]["sample"] == report
    assert manifest["run_fingerprint"] == "unchanged"


@pytest.mark.parametrize("raw", [None, "{}", "bad json", '{"p_pseudoaligned": -1}'])
def test_required_chemistry_gate_fails_closed_on_unavailable_or_bad_info(tmp_path, raw):
    kb = tmp_path / "kb-python"
    kb.mkdir()
    if raw is not None:
        (kb / "run_info.json").write_text(raw)
    ctx = RunContext.from_config(RunConfig(output=str(tmp_path), require_chemistry_sanity=True))
    with pytest.raises(RuntimeError, match="chemistry_sanity.json"):
        _chemistry_sanity(ctx)
    report = json.loads((tmp_path / "results" / "chemistry_sanity.json").read_text())
    assert report["status"] in {"unavailable", "failed"}


def test_host_subtracted_low_rate_does_not_fail_gate(tmp_path):
    kb = tmp_path / "kb-python"
    kb.mkdir()
    (kb / "run_info.json").write_text('{"p_pseudoaligned": 2.0}')
    _chemistry_sanity(
        RunContext.from_config(
            RunConfig(
                output=str(tmp_path), require_chemistry_sanity=True, host_filter_aligner="star"
            )
        )
    )


_HOST_ROW = "ENST1.1\tENSG1.1\t\t\tENST1.1\t1\t400\t+"
_VIRAL_ROW = "tx1\tEPSTEIN_HHV4_BZLF1\tBZLF1\t\tNC_007605.1\t1\t100\t+"


def _run_analysis(tmp_path, monkeypatch, t2g_lines, **config):
    """Run the analysis step on a low-rate kb run (p=0.3) and return the chemistry report."""
    from viralscan.scripts import analysis

    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "panel.gtf").write_text(
        'NC_007605.1\t.\tgene\t1\t100\t.\t+\t.\tgene_id "EPSTEIN_HHV4_BZLF1";\n'
    )
    monkeypatch.setattr(analysis, "ensure_viral_data", lambda cache_dir_arg=None: cache)
    t2g = tmp_path / "t2g.txt"
    t2g.write_text("".join(line + "\n" for line in t2g_lines))
    out = tmp_path / "out"
    (out / "log").mkdir(parents=True)
    (out / "kb-python").mkdir()
    (out / "kb-python" / "run_info.json").write_text('{"p_pseudoaligned": 0.3}')
    cfg = RunConfig(
        output=f"{out}/",
        transcripts=str(t2g),
        anellovirus_gene_ids=False,
        require_chemistry_sanity=True,
        **config,
    )
    analysis.run(RunContext.from_config(cfg))
    return json.loads((out / "results" / "chemistry_sanity.json").read_text())


def test_virus_only_index_without_host_filter_is_not_a_chemistry_error(tmp_path, monkeypatch):
    # FN3: a virus-only index has no host transcriptome, so a low rate is expected, not wrong.
    report = _run_analysis(tmp_path, monkeypatch, [_VIRAL_ROW])
    assert report["host_in_index"] is False and report["findings"] == []


def test_combined_index_without_host_filter_still_errors(tmp_path, monkeypatch):
    with pytest.raises(RuntimeError, match="chemistry_sanity.json"):
        _run_analysis(tmp_path, monkeypatch, [_HOST_ROW, _VIRAL_ROW])


def test_two_step_forward_strand_is_a_warning_not_a_failure(tmp_path, monkeypatch):
    report = _run_analysis(
        tmp_path,
        monkeypatch,
        [_HOST_ROW, _VIRAL_ROW],
        host_filter_aligner="starsolo",
        strand="forward",
    )
    assert [(f["level"], f["check"]) for f in report["findings"]] == [("warning", "strand")]

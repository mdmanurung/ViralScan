"""RUNDIR-01: post-hoc verbs find their run by ``--run-dir <run root>`` and an optional ``--sample``."""

from __future__ import annotations

from pathlib import Path

import pytest

from viralscan import menu, posthoc_record, rundir
from viralscan.cli.parsers import parse_args
from viralscan.discovery import DiscoveryError
from viralscan.run_safety import RunSafetyError, require_fresh_output

SAMPLE_VERB_ARGV = {
    "evidence": ["evidence", "--virus", "EBV"],
    "rerun-programs": ["rerun-programs"],
    "hostresponse": ["hostresponse", "--host-h5ad", "h.h5ad"],
    "host-loss": ["host-loss", "--truth", "t.tsv"],
}


def make_root(tmp_path: Path, samples=("S1",), complete: bool = True) -> Path:
    root = tmp_path / "run"
    for sample in samples:
        (root / sample / "kb-python").mkdir(parents=True)
        (root / sample / "config.yaml").write_text("output: x\n")
    (root / "run_manifest.json").write_text("{}\n")
    if complete:
        (root / "run_complete.json").write_text("{}\n")
    return root


def test_root_with_one_sample_resolves_it(tmp_path):
    root = make_root(tmp_path)
    assert rundir.resolve_sample(root) == (root / "S1").resolve()


def test_root_with_several_samples_needs_sample(tmp_path):
    root = make_root(tmp_path, samples=("S1", "S2"))
    with pytest.raises(DiscoveryError, match="--sample"):
        rundir.resolve_sample(root)
    assert rundir.resolve_sample(root, "S2") == (root / "S2").resolve()
    with pytest.raises(DiscoveryError, match="not in run"):
        rundir.resolve_sample(root, "S9")


def test_unfinished_run_is_refused(tmp_path):
    root = make_root(tmp_path, complete=False)
    with pytest.raises(DiscoveryError, match="not a completed ViralScan run"):
        rundir.resolve_sample(root)


def test_sample_dir_is_refused_unless_legacy_form_is_allowed(tmp_path):
    sample = make_root(tmp_path) / "S1"
    with pytest.raises(DiscoveryError):
        rundir.resolve_sample(sample)


def test_legacy_sample_dir_works_with_a_deprecation_warning(tmp_path, capsys):
    root = make_root(tmp_path)
    sample = (root / "S1").resolve()
    assert rundir.resolve_sample(sample, sample_dir_ok=True) == sample
    err = capsys.readouterr().err
    assert "deprecated" in err and "3.1.0" in err
    assert f"--run-dir {root.resolve()} --sample S1" in err


def test_legacy_sample_dir_with_a_wrong_sample_is_refused(tmp_path):
    sample = make_root(tmp_path) / "S1"
    with pytest.raises(DiscoveryError, match="does not match"):
        rundir.resolve_sample(sample, "S2", sample_dir_ok=True)


def test_root_form_prints_no_warning(tmp_path, capsys):
    rundir.resolve_sample(make_root(tmp_path), sample_dir_ok=True)
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("verb", sorted(SAMPLE_VERB_ARGV))
def test_every_sample_verb_takes_sample(verb):
    args = parse_args([*SAMPLE_VERB_ARGV[verb], "--run-dir", "run", "--sample", "S2"])
    assert args.sample == "S2"
    assert parse_args([*SAMPLE_VERB_ARGV[verb], "--run-dir", "run"]).sample is None


def _capture_handler(monkeypatch, verb: str, seen: list) -> None:
    def grab(args):
        seen.append(args)

    if verb == "evidence":
        from viralscan.scripts import evidence_run

        monkeypatch.setattr(evidence_run, "run_evidence", grab)
    else:
        monkeypatch.setattr(
            menu,
            {
                "rerun-programs": "_run_rerun_programs",
                "hostresponse": "_run_hostresponse_subcommand",
                "host-loss": "_run_host_loss",
            }[verb],
            grab,
        )


@pytest.mark.parametrize("verb", sorted(SAMPLE_VERB_ARGV))
def test_dispatch_hands_the_sample_directory_to_the_verb(tmp_path, monkeypatch, verb, capsys):
    root = make_root(tmp_path, samples=("S1", "S2"))
    seen: list = []
    _capture_handler(monkeypatch, verb, seen)
    argv = [*SAMPLE_VERB_ARGV[verb], "--run-dir", str(root), "--sample", "S2"]
    menu._run_subcommand(parse_args(argv))
    assert [a.run_dir for a in seen] == [str((root / "S2").resolve())]
    assert seen[0]._posthoc_sample == "S2"
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("verb", sorted(SAMPLE_VERB_ARGV))
def test_dispatch_still_takes_a_sample_directory(tmp_path, monkeypatch, verb, capsys):
    root = make_root(tmp_path)
    seen: list = []
    _capture_handler(monkeypatch, verb, seen)
    menu._run_subcommand(parse_args([*SAMPLE_VERB_ARGV[verb], "--run-dir", str(root / "S1")]))
    assert [a.run_dir for a in seen] == [str((root / "S1").resolve())]
    assert "--sample S1" in capsys.readouterr().err


def test_dispatch_of_an_ambiguous_root_is_one_error_line(tmp_path, caplog):
    root = make_root(tmp_path, samples=("S1", "S2"))
    with pytest.raises(SystemExit) as exc:
        menu._dispatch(parse_args(["rerun-programs", "--run-dir", str(root)]))
    assert exc.value.code == 1
    assert "choose one with --sample" in caplog.text


def test_status_warns_for_a_sample_directory(tmp_path, capsys):
    root = make_root(tmp_path)
    menu._run_subcommand(parse_args(["status", "--run-dir", str(root / "S1")]))
    assert "deprecated" in capsys.readouterr().err
    menu._run_subcommand(parse_args(["status", "--run-dir", str(root)]))
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("force", [False, True])
def test_guard_refuses_output_in_the_primary_tree(tmp_path, force):
    root = make_root(tmp_path)
    for inside in (root / "S1" / "results" / "x", root / "S1", root / "other"):
        with pytest.raises(RunSafetyError, match="inside the ViralScan run"):
            require_fresh_output(inside, force=force)
    require_fresh_output(root / "S1" / "posthoc" / "evidence" / "EBV", force=force)


def test_discover_output_uses_the_same_guard(tmp_path):
    from viralscan import discovery

    root = make_root(tmp_path)
    with pytest.raises(DiscoveryError, match="inside the ViralScan run"):
        discovery.check_output_dir(root / "S1" / "results" / "discovery", root)
    discovery.check_output_dir(root / "S1" / "posthoc" / "discovery", root)


def test_record_of_a_root_form_verb_names_the_sample_and_hashes_its_config(tmp_path):
    root = make_root(tmp_path)
    args = parse_args(["rerun-programs", "--run-dir", str(root)])
    current = posthoc_record.Attempt("rerun-programs", args, [])
    current.start()
    rundir.use_sample_dir(args)  # what _run_subcommand does before the verb
    current.finish("succeeded", 0, None)
    record = current.record
    assert record["sample"] == "S1"
    assert record["source_root"] == str(root.resolve())
    assert {Path(i["path"]).name for i in record["inputs"]} == {"run_manifest.json", "config.yaml"}
    assert record["settings"]["run_dir"] == str(root)  # what the user typed

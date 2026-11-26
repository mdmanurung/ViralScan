"""`run` fails before kb count when the bundled annotations are missing (PLAN PREFLIGHT-GTF-01)."""

import pytest

from viralscan import menu
from viralscan.cli.parsers import build_parser
from viralscan.data_fetch import ViralScanDataError
from viralscan.validation import doctor_report


@pytest.fixture
def empty_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("VIRALSCAN_CACHE", str(tmp_path / "cache"))


@pytest.fixture
def run_argv(tmp_path):
    paths = {n: tmp_path / n for n in ("idx", "t2g", "a_R1.fastq", "a_R2.fastq", "my.gtf")}
    for p in paths.values():
        p.write_text("")
    return [
        "run", "-i", str(paths["idx"]), "-t", str(paths["t2g"]),
        "--r1", str(paths["a_R1.fastq"]), "--r2", str(paths["a_R2.fastq"]),
        "-o", str(tmp_path / "out"), "--yes",
    ]  # fmt: skip


def _boom(*_a, **_k):
    raise AssertionError("work started after a failed preflight")


def test_no_gtf_and_empty_cache_fails_before_any_work(run_argv, empty_cache, monkeypatch):
    monkeypatch.setattr(menu, "_check_required_tools", lambda: None)
    monkeypatch.setattr(menu, "validate_run_request", lambda a: None)
    for name in ("_run_snakemake", "acquire_run_lock", "_resolve_chemistry", "_check_sample_ids"):
        monkeypatch.setattr(menu, name, _boom)
    monkeypatch.setattr("subprocess.run", _boom)
    monkeypatch.setattr("subprocess.Popen", _boom)
    with pytest.raises(ViralScanDataError, match="viralscan data fetch"):
        menu._run_subcommand(build_parser().parse_args(run_argv))
    with pytest.raises(SystemExit) as exc:
        menu._dispatch(build_parser().parse_args(run_argv))
    assert exc.value.code == 1


def test_gtf_given_passes_preflight(run_argv, tmp_path, empty_cache):
    args = build_parser().parse_args(run_argv + ["--gtf", str(tmp_path / "my.gtf")])
    menu._preflight_bundled_annotations(args)


def test_index_build_manifest_passes_preflight(run_argv, tmp_path, empty_cache):
    args = build_parser().parse_args(run_argv)
    with pytest.raises(ViralScanDataError):
        menu._preflight_bundled_annotations(args)
    (tmp_path / "idx.build_manifest.json").write_text("{}")
    menu._preflight_bundled_annotations(args)


def test_doctor_warns_but_stays_ok_without_annotations(empty_cache):
    status = doctor_report("pip")["viral_annotations"]
    assert status["present"] is False
    assert "viralscan data fetch" in status["warning"]
    assert "--gtf" in status["warning"]
    assert doctor_report("pip")["ok"] is True


def test_check_env_prints_the_warning(empty_cache, monkeypatch, capsys):
    args = build_parser().parse_args(["check", "env", "--profile", "pip"])
    monkeypatch.setattr(menu, "create_help", lambda: args)
    menu._run_subcommand(args)
    assert "warning: Viral reference annotations were not found" in capsys.readouterr().err

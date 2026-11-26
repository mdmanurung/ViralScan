"""VALIDATE-ONCE-01: one validation pass on the typed RunConfig, before the output directory."""

from pathlib import Path

import pytest

from viralscan import menu
from viralscan.cli.parsers import build_parser
from viralscan.runconfig import RunConfig


@pytest.fixture
def args(tmp_path):
    files = [tmp_path / n for n in ("idx", "t2g", "a_R1.fastq", "a_R2.fastq")]
    for f in files:
        f.write_text("")
    return build_parser().parse_args(
        ["run", "-i", str(files[0]), "-t", str(files[1]), "--r1", str(files[2])]
        + ["--r2", str(files[3]), "-o", str(tmp_path / "out")]
    )


def test_validation_returns_the_checked_options_and_leaves_the_output_alone(args, monkeypatch):
    checks = []
    real = RunConfig.check_options.__func__
    monkeypatch.setattr(
        RunConfig,
        "check_options",
        classmethod(lambda cls, cfg: checks.append(cfg) or real(cls, cfg)),
    )
    checked = menu.validate_run_request(args)
    assert len(checks) == 1
    assert checked["cores"] == args.cores
    assert not Path(args.output).exists()


def test_a_multi_sample_run_checks_the_options_once(tmp_path, monkeypatch):
    from unittest.mock import patch

    names = ("a_R1.fastq", "a_R2.fastq", "b_R1.fastq", "b_R2.fastq", "index.idx", "t2g.txt")
    (tmp_path / "index.idx.build_manifest.json").write_text("{}")  # no bundled panel needed
    for name in names:
        (tmp_path / name).write_text("@r\nACGT\n+\nIIII\n" if "fastq" in name else "t\tg\n")
    p = lambda n: str(tmp_path / n)  # noqa: E731
    argv = ["viralscan", "-i", p("index.idx"), "-t", p("t2g.txt"), "-o", p("out")]
    argv += ["--r1", f"{p('a_R1.fastq')},{p('b_R1.fastq')}", "-x", "10xv3", "--yes"]
    argv += ["--r2", f"{p('a_R2.fastq')},{p('b_R2.fastq')}", "--strand", "unstranded"]
    monkeypatch.setattr("sys.argv", argv)
    checks, builds = [], []
    real_check = RunConfig.check_options.__func__
    real_build = RunConfig.from_snakemake_config.__func__
    monkeypatch.setattr(
        RunConfig,
        "check_options",
        classmethod(lambda cls, cfg: checks.append(cfg) or real_check(cls, cfg)),
    )
    monkeypatch.setattr(
        RunConfig,
        "from_snakemake_config",
        classmethod(
            lambda cls, cfg, checked=None: (
                builds.append(checked) or real_build(cls, cfg, checked=checked)
            )
        ),
    )
    with (
        patch.object(menu, "_check_required_tools"),
        patch.object(menu, "_check_cell_caller_tools"),
        patch.object(menu, "_resolve_chemistry", return_value={}),
        patch.object(menu, "_run_snakemake"),
    ):
        menu._main()
    assert len(checks) == 1  # the probe; neither sample re-runs the flag rules
    assert len(builds) == 2
    assert all(b is not None for b in builds)


def test_a_build_without_checked_options_still_runs_the_rules(args):
    cfg = menu._run_config_options(args, "o/", "i", "t", "r1", "r2")
    cfg["cores"] = 0
    with pytest.raises(ValueError, match="cores must be >= 1"):
        RunConfig.from_snakemake_config(cfg)


@pytest.mark.parametrize(
    ("flag", "value", "message"),
    [
        ("cores", 0, "cores must be >= 1, got 0."),
        ("emptydrops_fdr", 1.5, "emptydrops_fdr must be in (0, 1), got 1.5."),
        ("hvg_max_mean", 0.0, "hvg_min_mean (0.0125) must be < hvg_max_mean (0.0)."),
    ],
)
def test_numeric_range_is_refused_once_by_runconfig(args, caplog, flag, value, message):
    setattr(args, flag, value)
    menu.errorhandler(args)  # errorhandler no longer repeats the RunConfig rule
    with pytest.raises(SystemExit) as exc:
        menu.validate_run_request(args)
    assert exc.value.code == 2
    assert message in caplog.text


def test_cli_only_message_comes_before_the_runconfig_message(args, caplog):
    args.cores = 0
    args.cell_calling = "external"  # errorhandler: --cell-calling external requires a file
    with pytest.raises(SystemExit) as exc:
        menu.validate_run_request(args)
    assert exc.value.code == 2
    assert "--cell-calling external requires --called-cells-file." in caplog.text
    assert "cores must be" not in caplog.text


def test_main_calls_the_single_entry(monkeypatch):
    import inspect

    src = inspect.getsource(menu._run_subcommand)
    assert "validate_run_request(args)" in src
    assert "errorhandler(args)" not in src
    assert "_precheck_run_config(args)" not in src

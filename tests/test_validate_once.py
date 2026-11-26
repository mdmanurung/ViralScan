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


def test_one_runconfig_is_built_and_the_output_directory_is_untouched(args, monkeypatch):
    built = []
    real = RunConfig.from_snakemake_config.__func__
    monkeypatch.setattr(
        RunConfig,
        "from_snakemake_config",
        classmethod(lambda cls, cfg: built.append(cfg) or real(cls, cfg)),
    )
    menu.validate_run_request(args)
    assert len(built) == 1
    assert not Path(args.output).exists()


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

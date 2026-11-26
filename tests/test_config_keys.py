"""CONFIG-KEYS-01: the four cell-calling config keys are `run` flags at their current defaults."""

import argparse
import dataclasses
import logging
from pathlib import Path

import pytest

import viralscan.run_safety as rs
from viralscan import menu
from viralscan.cli.parsers import _build_run_parser, parse_args
from viralscan.defaults import DEFAULTS
from viralscan.runconfig import RunConfig

FLAGS = {
    "emptydrops_fdr": ("--emptydrops-fdr", 0.05),
    "emptydrops_lower": ("--emptydrops-lower", 200),
    "knee_min_umi": ("--knee-min-umi", 20.0),
    "cell_caller_rscript": ("--cell-caller-rscript", "/opt/R/bin/Rscript"),
}


def _run_parser() -> argparse.ArgumentParser:
    sub = argparse.ArgumentParser().add_subparsers()
    _build_run_parser(sub)
    return sub.choices["run"]


def _defaults() -> dict:
    return {a.dest: a.default for a in _run_parser()._actions if a.dest != "help"}


def _args(tmp_path: Path, **extra) -> argparse.Namespace:
    for name, text in [("R1.fastq", "r1\n"), ("R2.fastq", "r2\n"), ("idx.idx", "x\n")]:
        (tmp_path / name).write_text(text)
    (tmp_path / "t2g.txt").write_text("a\tb\n")
    base = _defaults()
    base.update(
        sample1=str(tmp_path / "R1.fastq"),
        sample2=str(tmp_path / "R2.fastq"),
        index=str(tmp_path / "idx.idx"),
        transcripts=str(tmp_path / "t2g.txt"),
        output=str(tmp_path / "out"),
        strand="forward",
    )
    base.update(extra)
    return argparse.Namespace(**base)


@pytest.fixture(autouse=True)
def _fixed_version(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(rs, "__version__", "0.0.test")


@pytest.mark.parametrize("key", sorted(FLAGS))
def test_flag_default_is_the_defaults_table_value(key: str) -> None:
    assert _defaults()[key] == DEFAULTS[key]
    assert RunConfig.__dataclass_fields__[key].default == DEFAULTS[key]


def test_defaults_are_unchanged() -> None:
    """Impactful choices (decision row D-CC): a flag must not move a value."""
    assert {k: DEFAULTS[k] for k in FLAGS} == {
        "emptydrops_fdr": 0.01,
        "emptydrops_lower": 100,
        "knee_min_umi": 10.0,
        "cell_caller_rscript": "Rscript",
    }


@pytest.mark.parametrize("key", sorted(FLAGS))
def test_flag_reaches_the_run_config(key: str) -> None:
    flag, value = FLAGS[key]
    args = parse_args(
        [
            "-i",
            "i",
            "-t",
            "t",
            "--r1",
            "a_R1.fastq",
            "--r2",
            "b_R2.fastq",
            "-o",
            "o",
            flag,
            str(value),
        ]
    )
    cfg = menu._build_run_config(args, "out/", "i", "t", "R1.fq", "R2.fq")
    assert getattr(cfg, key) == value
    plain = parse_args(
        ["run", "-i", "i", "-t", "t", "--r1", "a_R1.fastq", "--r2", "b_R2.fastq", "-o", "o"]
    )
    default_cfg = menu._build_run_config(plain, "out/", "i", "t", "R1.fq", "R2.fq")
    assert getattr(default_cfg, key) == DEFAULTS[key]


@pytest.mark.parametrize(("flag", "bad"), [("--emptydrops-fdr", "1.5"), ("--knee-min-umi", "0")])
def test_out_of_range_value_is_refused(flag: str, bad: str) -> None:
    args = parse_args(
        [
            "run",
            "-i",
            "i",
            "-t",
            "t",
            "--r1",
            "a_R1.fastq",
            "--r2",
            "b_R2.fastq",
            "-o",
            "o",
            flag,
            bad,
        ]
    )
    with pytest.raises(ValueError, match="must be"):
        menu._build_run_config(args, "out/", "i", "t", "R1.fq", "R2.fq")


@pytest.mark.parametrize("key", sorted(FLAGS))
def test_a_default_flag_hashes_like_an_unset_one(tmp_path: Path, key: str) -> None:
    unset = _args(tmp_path, **{key: None})  # a caller that leaves the option unset
    explicit = _args(tmp_path, **{key: DEFAULTS[key]})
    a = rs.build_run_manifest(unset)
    b = rs.build_run_manifest(explicit)
    assert a["options"][key] == DEFAULTS[key]
    assert a["run_fingerprint"] == b["run_fingerprint"]


@pytest.mark.parametrize("key", sorted(FLAGS))
def test_a_non_default_flag_changes_the_fingerprint(tmp_path: Path, key: str) -> None:
    value = FLAGS[key][1]
    assert value != DEFAULTS[key]
    plain = rs.build_run_manifest(_args(tmp_path))
    changed = rs.build_run_manifest(_args(tmp_path, **{key: value}))
    assert changed["run_fingerprint"] != plain["run_fingerprint"]
    assert changed["options"][key] == value


def test_other_options_keep_their_fingerprint_inputs(tmp_path: Path) -> None:
    """The four flags add four option keys at their defaults and nothing else."""
    options = rs.build_run_manifest(_args(tmp_path))["options"]
    assert {k: options[k] for k in FLAGS} == {k: DEFAULTS[k] for k in FLAGS}


@pytest.mark.parametrize("key", ["overwrite", "f1", "reference"])
def test_unreachable_run_config_fields_are_gone(key: str) -> None:
    assert key not in {f.name for f in dataclasses.fields(RunConfig)}
    assert key not in RunConfig().to_dict()


def test_an_old_config_yaml_with_the_removed_keys_still_loads(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text("output: out/\noverwrite: 'yes'\nf1: cdna.fa\nreference: false\n")
    assert RunConfig.from_yaml(path).output == "out/"


def test_ignored_flag_warning_names_the_new_flags(caplog: pytest.LogCaptureFixture) -> None:
    args = parse_args(
        [
            "-i",
            "i",
            "-t",
            "t",
            "--r1",
            "a_R1.fastq",
            "--r2",
            "b_R2.fastq",
            "-o",
            "o",
            "--cell-calling",
            "none",
        ]
        + ["--emptydrops-fdr", "0.05", "--emptydrops-lower", "200", "--knee-min-umi", "20"]
        + ["--cell-caller-rscript", "R2script"]
    )
    with caplog.at_level(logging.WARNING):
        menu._warn_ignored_flags(args)
    text = caplog.text
    for flag in ("--emptydrops-fdr", "--emptydrops-lower", "--cell-caller-rscript"):
        assert flag in text
    assert "--knee-min-umi ignored: cell calling is none, not knee" in text

"""FLAGS-01: removed no-effect flags, deprecated renames and the essential ``run --help``."""

import argparse
import re
from pathlib import Path

import pytest

import viralscan.run_safety as rs
from viralscan.cli.parsers import (
    DEPRECATED_OPTIONS,
    RUN_ESSENTIAL_DESTS,
    build_parser,
    parse_args,
)

RUN = ["run", "-o", "o", "--r1", "a", "--r2", "b", "-t", "t", "-i", "i"]


@pytest.mark.parametrize(
    ("extra", "needle"),
    [
        (["--f1", "x.fa"], "--f1 was removed"),
        (["-f1", "x.fa"], "-f1 was removed"),
        (["--se-threshold", "3"], "--se-threshold was removed"),
        (["--multimap-primary-call", "selected-method"], "--multimap-primary-call was removed"),
        (["--positive-control-scope", "virus_key"], "virus_key was removed"),
        (["--positive-control-scope=virus_key"], "virus_key was removed"),
    ],
)
def test_removed_flag_gives_a_clear_error(extra, needle, capsys):
    with pytest.raises(SystemExit) as exc:
        parse_args(RUN + extra)
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert needle in err
    assert "unrecognized arguments" not in err and "invalid choice" not in err


def test_virus_key_error_names_the_replacements(capsys):
    with pytest.raises(SystemExit):
        parse_args(RUN + ["--positive-control-scope", "virus_key"])
    assert "exact_sequence" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("argv", "dest", "value"),
    [
        (["run", "--programme-min-umi", "2"], "programme_min_umi", 2.0),
        (["run", "--hostresponse-min-umi", "7"], "hostresponse_min_umi", 7),
        (["rerun", "programs", "--run-dir", "r", "--programme-min-umi", "2"], "programme_min_umi", 2.0),
        (["rerun", "hostresponse", "--run-dir", "r", "--host-h5ad", "h", "--min-umi", "5"], "hostresponse_min_umi", 5),
        (["discover", "run", "--run-dir", "r", "--kraken2-db", "k", "--host-subtract", "starsolo"], "host_subtract", "starsolo"),
        (["discover", "run", "--run-dir", "r", "--kraken2-db", "k", "--called-cells=c.tsv"], "called_cells", "c.tsv"),
        (["translated", "--run-dir", "r", "--index-dir", "d", "--cells", "c.tsv"], "cells", "c.tsv"),
    ],
)  # fmt: skip
def test_old_spelling_still_parses_with_a_deprecation_line(argv, dest, value, capsys):
    if argv[0] == "run":
        argv = RUN + argv[1:]
    args = parse_args(argv)
    assert getattr(args, dest) == value
    err = capsys.readouterr().err
    assert "is deprecated; use" in err and "3.1.0" in err


def test_every_new_spelling_parses_without_a_warning(capsys):
    args = parse_args(RUN + ["--programme-min-molecules", "2", "--hostresponse-min-molecules", "7"])
    assert (args.programme_min_umi, args.hostresponse_min_umi) == (2.0, 7)
    assert capsys.readouterr().err == ""


def test_each_alias_points_at_a_real_option():
    sub = next(
        a for a in build_parser()._actions if isinstance(a, argparse._SubParsersAction)
    ).choices
    canonical = {"rerun-programs": "rerun programs", "hostresponse": "rerun hostresponse"}
    for verb, renames in DEPRECATED_OPTIONS.items():
        parser = sub[verb] if verb in sub else sub[canonical[verb].split()[0]]
        if verb in canonical:
            group = sub[canonical[verb].split()[0]]
            parser = next(
                a for a in group._actions if isinstance(a, argparse._SubParsersAction)
            ).choices[canonical[verb].split()[1]]
        if verb == "discover":
            group = sub["discover"]
            parser = next(
                a for a in group._actions if isinstance(a, argparse._SubParsersAction)
            ).choices["run"]
        options = {o for a in parser._actions for o in a.option_strings}
        for old, new in renames.items():
            assert new in options and old not in options, (verb, old, new)


def _run_parser(essential: bool) -> argparse.ArgumentParser:
    root = build_parser(essential_help=essential)
    return next(a for a in root._actions if isinstance(a, argparse._SubParsersAction)).choices[
        "run"
    ]


def test_run_help_lists_only_the_essential_flags_and_help_all_the_rest():
    essential_dests = {
        a.dest
        for a in _run_parser(True)._actions
        if a.help != argparse.SUPPRESS and a.option_strings
    }
    assert essential_dests == RUN_ESSENTIAL_DESTS | {"help", "help_all"}
    assert 14 <= len(RUN_ESSENTIAL_DESTS) <= 17
    shown = _run_parser(True).format_help()
    full = _run_parser(False).format_help()
    assert "--gene-programs" not in shown and "--gene-programs" in full
    assert "--help-all" in shown


def test_hidden_flags_still_parse():
    args = parse_args(RUN + ["--gene-programs", "--emptydrops-seed", "5"])
    assert args.gene_programs is True and args.emptydrops_seed == 5


def test_help_all_prints_every_option(capsys):
    with pytest.raises(SystemExit) as exc:
        parse_args(["run", "--help-all"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    for flag in ("--gene-programs", "--programme-min-molecules", "--emptydrops-seed"):
        assert flag in out
    for gone in ("--f1", "--se-threshold", "--multimap-primary-call"):
        assert gone not in out


def test_fingerprint_of_kept_flags_is_unchanged(tmp_path, monkeypatch):
    """The digest was taken on origin/main before FLAGS-01 with the same argv and inputs."""
    monkeypatch.setattr(rs, "__version__", "0.0.test")
    monkeypatch.chdir(tmp_path)
    for name in ("R1.fastq", "R2.fastq", "idx.idx"):
        Path(name).write_text("x\n")
    Path("t2g.txt").write_text("a\tb\n")
    argv = [
        "run", "-t", "t2g.txt", "-i", "idx.idx", "-o", "o", "--r1", "R1.fastq",
        "--r2", "R2.fastq", "--technology", "10xv3", "--strand", "forward",
        "--detection-threshold", "3", "--programme-min-breadth", "3",
        "--multimap-method", "equal", "--cell-calling", "none", "--umap",
    ]  # fmt: skip
    manifest = rs.build_run_manifest(parse_args(argv))
    assert re.fullmatch(r"[0-9a-f]{64}", manifest["run_fingerprint"])
    assert manifest["run_fingerprint"] == (
        "5281c53c4d30f81811aac9f1fef71ee564c97c98baf5d8ef8e48820fb5bb75ab"
    )
    # The removed flags stay hashed at their fixed values, which is what keeps the digest.
    assert manifest["options"]["se_threshold"] == 10
    assert manifest["options"]["multimap_primary_call"] == "selected-method"
    assert not any(key.startswith("f1:") for key in manifest["input_fingerprints"])


def test_accession_lists_take_commas_or_spaces():
    base = ["build", "ref", "--host", "human"]
    for tokens in (["A", "B"], ["A,B"], ["A,", "B"]):
        assert parse_args(base + ["--virus-accessions", *tokens]).virus_accessions == ["A", "B"]

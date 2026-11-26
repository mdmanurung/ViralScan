"""CMDGROUP-01: grouped verbs, and the old top-level names kept as warning aliases."""

from __future__ import annotations

import pytest

from viralscan import menu, posthoc_record, rundir
from viralscan.cli.parsers import COMMAND_GROUPS, parse_args

# group path -> (canonical verb, minimal argv after the path that parses)
CASES = {
    "check chemistry": ("check-chemistry", ["--r1", "a", "--r2", "b", "-i", "i", "-t", "t"]),
    "check whitelist": ("check-whitelist", ["--r1", "a", "-w", "w", "-x", "10xv3"]),
    "check run": ("validate-run", ["--run-dir", "r"]),
    "check env": ("doctor", []),
    "discover run": ("discover", ["--run-dir", "r", "--kraken2-db", "db"]),
    "discover collapse": ("discover-collapse", ["d", "--rank", "genus"]),
    "discover states": ("discover-states", ["d", "--run-dir", "r"]),
    "rerun multimap": ("rerun-multimap", ["--run-dir", "r", "--multimap-method", "em"]),
    "rerun programs": ("rerun-programs", ["--run-dir", "r"]),
    "rerun hostresponse": ("hostresponse", ["--run-dir", "r", "--host-h5ad", "h"]),
    "build ref": ("build-ref", ["--host", "human", "-o", "o"]),
    "build translated-index": (
        "translated-index",
        ["--proteins-dir", "p", "--host-mask", "m", "-o", "o"],
    ),
    "build star-index": (
        "star-index",
        [f"--{n}=x" for n in ("host-fasta", "host-gtf", "viral-fasta", "viral-gtf", "out-dir")],
    ),
    "dev host-loss": ("host-loss", ["--run-dir", "r", "--truth", "t"]),
}


def test_every_group_word_is_covered():
    listed = {f"{g} {w}" for g, words in COMMAND_GROUPS.items() for w in words}
    assert listed == set(CASES)


@pytest.mark.parametrize("path", CASES)
def test_group_path_dispatches_to_the_canonical_verb(path):
    verb, rest = CASES[path]
    args = parse_args([*path.split(), *rest])
    assert args._subcommand == verb
    assert getattr(args, "_alias_for", None) is None


@pytest.mark.parametrize("path", CASES)
def test_old_name_is_an_alias_of_the_same_verb(path):
    verb, rest = CASES[path]
    args = parse_args([verb, *rest])
    assert args._subcommand == verb
    assert args._alias_for == path


def test_old_and_new_spelling_parse_to_the_same_options():
    for path, (verb, rest) in CASES.items():
        new, old = vars(parse_args([*path.split(), *rest])), vars(parse_args([verb, *rest]))
        old.pop("_alias_for")
        new.pop("_group_command", None)
        new.pop("_missing_group", None)
        old.pop("_legacy", None)
        assert new == old, path


def test_keyed_on_the_canonical_verb():
    assert "rerun-programs" in rundir.SAMPLE_VERBS
    for path in ("rerun programs", "rerun hostresponse", "dev host-loss"):
        verb, _ = CASES[path]
        assert verb in rundir.SAMPLE_VERBS
        assert verb in posthoc_record.VERBS


def test_group_without_command_is_a_usage_error():
    with pytest.raises(SystemExit) as exc:
        menu._run_subcommand(parse_args(["rerun"]))
    assert exc.value.code == 2


def test_old_name_warns_and_new_name_does_not(monkeypatch, capsys):
    monkeypatch.setattr(menu, "_dispatch", lambda args: None)
    monkeypatch.setattr("sys.argv", ["viralscan", "validate-run", "--run-dir", "r"])
    menu._main()
    err = capsys.readouterr().err
    assert "viralscan validate-run is deprecated; use viralscan check run." in err
    assert "3.1.0" in err
    monkeypatch.setattr("sys.argv", ["viralscan", "check", "run", "--run-dir", "r"])
    menu._main()
    assert "deprecated" not in capsys.readouterr().err


def test_hidden_aliases_stay_out_of_the_help(capsys):
    with pytest.raises(SystemExit):
        parse_args(["--help"])
    out = capsys.readouterr().out
    assert "check chemistry" in out
    assert "  validate-run" not in out
    assert "dev" not in out.split("Subcommands:")[1].split("The old names")[0]


def test_removed_option_is_named_under_the_group_path(capsys):
    with pytest.raises(SystemExit) as exc:
        parse_args(["rerun", "hostresponse", "--run-dir", "r", "--host-h5ad", "h", "-o", "x"])
    assert exc.value.code == 2
    assert "viralscan rerun hostresponse: error: -o was removed; use --run-dir" in (
        capsys.readouterr().err
    )

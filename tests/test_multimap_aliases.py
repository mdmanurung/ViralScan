"""--multimap-method shows 3 options; short names are CLI aliases only (PLAN Q72-Q74)."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from viralscan.defaults import DEFAULTS, MULTIMAP_METHOD_ALIASES, MULTIMAP_METHODS
from viralscan.menu import build_parser
from viralscan.run_safety import build_run_manifest

EM_TARGET = MULTIMAP_METHOD_ALIASES["em"]
HIDDEN = ("unique-weighted", "sibling-weighted", "em-global")


def _main(*extra: str) -> argparse.Namespace:
    return build_parser().parse_args(["-s1", "a", "-s2", "b", "-o", "o", *extra])


def _rerun(method: str) -> argparse.Namespace:
    argv = ["rerun-multimap", "--run-dir", "x", "-o", "y", "--multimap-method", method]
    return build_parser().parse_args(argv)


def _action(parser: argparse.ArgumentParser) -> argparse.Action:
    return next(a for a in parser._actions if "--multimap-method" in a.option_strings)


def test_em_alias_targets_an_em_engine() -> None:
    assert EM_TARGET in {"em-cell", "em-global"}


@pytest.mark.parametrize(
    ("given", "canonical"),
    [("conservative", "host-conservative"), ("equal", "equal"), ("em", EM_TARGET)],
)
def test_aliases_resolve_to_canonical_in_both_parsers(given: str, canonical: str) -> None:
    assert _main("--multimap-method", given).multimap_method == canonical
    assert _rerun(given).multimap_method == canonical


def test_default_is_conservative_and_canonical() -> None:
    assert _main().multimap_method == "host-conservative" == DEFAULTS["multimap_method"]


@pytest.mark.parametrize("method", MULTIMAP_METHODS)
def test_every_canonical_name_still_accepted_silently(
    method: str, capsys: pytest.CaptureFixture[str]
) -> None:
    assert _main("--multimap-method", method).multimap_method == method
    assert _rerun(method).multimap_method == method
    captured = capsys.readouterr()
    assert captured.out == captured.err == ""


def test_help_lists_exactly_three_options() -> None:
    for parser in (
        build_parser(),
        build_parser()._subparsers._group_actions[0].choices["rerun-multimap"],
    ):
        action = _action(parser)
        assert action.choices is None
        assert action.metavar == "{conservative,em,equal}"
        text = parser.format_help()
        assert "{conservative,em,equal}" in text
    main_action = _action(build_parser())
    for hidden in (*HIDDEN, "em-cell", "host-conservative"):
        assert hidden not in (main_action.help or "")


def test_unknown_method_rejected_and_message_names_the_three(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit):
        _main("--multimap-method", "bogus")
    err = capsys.readouterr().err
    assert "conservative, em, equal" in err
    assert "unique-weighted" not in err


@pytest.mark.parametrize(
    ("given", "canonical"), [("conservative", "host-conservative"), ("em", EM_TARGET)]
)
def test_manifest_records_canonical_name(tmp_path: Path, given: str, canonical: str) -> None:
    r1, r2 = tmp_path / "R1.fastq", tmp_path / "R2.fastq"
    r1.write_text("r1\n")
    r2.write_text("r2\n")
    args = build_parser().parse_args(
        ["-s1", str(r1), "-s2", str(r2), "-o", str(tmp_path / "out"), "--multimap-method", given]
    )
    assert build_run_manifest(args)["allocation_method"] == canonical

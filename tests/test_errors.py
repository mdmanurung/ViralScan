"""ERRORS-01: one ``ViralScanError`` base; ``menu._dispatch`` turns it into ``ERROR`` and exit 1."""

import argparse
import ast
import logging
from pathlib import Path

import pytest

from viralscan import chemistry, data_fetch, menu
from viralscan.errors import ViralScanError

SRC = Path(__file__).resolve().parents[1] / "src" / "viralscan"
_BUILTIN_BASES = {"Exception", "ValueError", "RuntimeError", "OSError", "KeyError", "TypeError"}


def _exception_classes() -> dict[str, list[str]]:
    """Every exception class defined in the package: name -> base names (by AST)."""
    classes: dict[str, list[str]] = {}
    for path in sorted(SRC.rglob("*.py")):
        if "vendor" in path.parts:
            continue
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.ClassDef):
                bases = [b.id for b in node.bases if isinstance(b, ast.Name)]
                bases += [b.attr for b in node.bases if isinstance(b, ast.Attribute)]
                classes[node.name] = bases
    found: set[str] = set()
    changed = True
    while changed:
        changed = False
        for name, bases in classes.items():
            if name not in found and any(b in _BUILTIN_BASES or b in found for b in bases):
                found.add(name)
                changed = True
    return {name: classes[name] for name in found}


def test_every_package_exception_inherits_the_base() -> None:
    classes = _exception_classes()
    assert "ChemistryError" in classes and "NCBIFetchError" in classes
    rooted = {n for n, bases in classes.items() if "ViralScanError" in bases}
    rooted |= {n for n, bases in classes.items() if any(b in rooted for b in bases)}
    assert set(classes) - {"ViralScanError"} <= rooted, sorted(set(classes) - rooted)


def test_package_exceptions_keep_their_older_bases() -> None:
    assert issubclass(chemistry.ChemistryError, ValueError)
    assert issubclass(data_fetch.ViralScanDataError, RuntimeError)
    assert issubclass(chemistry.ChemistryError, ViralScanError)
    assert issubclass(data_fetch.ViralScanDataError, ViralScanError)
    assert issubclass(ViralScanError, Exception) and not issubclass(ViralScanError, ValueError)


def test_dispatch_turns_the_base_into_one_error_line_and_exit_1(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    def boom(args: argparse.Namespace) -> None:
        raise chemistry.ChemistryError("no such chemistry")

    monkeypatch.setattr(menu, "_run_subcommand", boom)
    with caplog.at_level(logging.ERROR), pytest.raises(SystemExit) as info:
        menu._dispatch(argparse.Namespace())
    assert info.value.code == 1
    assert "no such chemistry" in caplog.text


def test_dispatch_leaves_other_errors_alone(monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(args: argparse.Namespace) -> None:
        raise ValueError("plain")

    monkeypatch.setattr(menu, "_run_subcommand", boom)
    with pytest.raises(ValueError, match="plain"):
        menu._dispatch(argparse.Namespace())


def test_data_fetch_failure_still_exits_1(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(**kwargs: object) -> None:
        raise data_fetch.ViralScanDataError("checksum mismatch")

    monkeypatch.setattr(data_fetch, "fetch_viral_data", fail)
    args = argparse.Namespace(
        _subcommand="data-fetch",
        verbose=False,
        quiet=False,
        cache_dir=None,
        url=None,
        sha256=None,
        force=False,
    )
    with pytest.raises(SystemExit) as info:
        menu._dispatch(args)
    assert info.value.code == 1


def test_evidence_plot_usage_error_still_exits_2(monkeypatch: pytest.MonkeyPatch) -> None:
    from viralscan import competitive_plot

    def bad(args: argparse.Namespace) -> None:
        raise chemistry.ChemistryError("bad flag")  # a ViralScanError that is also a ValueError

    monkeypatch.setattr(competitive_plot, "run_evidence_plot", bad)
    args = argparse.Namespace(_subcommand="evidence-plot", verbose=False, quiet=False)
    with pytest.raises(SystemExit) as info:
        menu._dispatch(args)
    assert info.value.code == 2

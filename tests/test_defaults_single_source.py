"""One defaults table: parser, RunConfig and script fallbacks must agree with DEFAULTS.

PLAN CONTRACT-DEFAULTS-01. A default copied into a second place drifts; this guard fails when it does.
"""

import argparse
import dataclasses
import inspect

import pytest

from viralscan.cli.parsers import _build_run_parser
from viralscan.defaults import DEFAULTS
from viralscan.runconfig import RunConfig

# Parser flags whose default is None on purpose: RunConfig resolves them from DEFAULTS.
_RESOLVED_LATER = {
    "detection_min_cells",
    "programme_mixed_veto",
    "programme_mixed_veto_alpha",
    "auto_evidence",
    *(key for key in DEFAULTS if key.startswith("hostresponse_")),
}


def _run_parser() -> argparse.ArgumentParser:
    subparsers = argparse.ArgumentParser().add_subparsers()
    _build_run_parser(subparsers)
    return subparsers.choices["run"]


def test_run_parser_defaults_equal_defaults_table():
    parser = _run_parser()
    wrong = {}
    for key, value in DEFAULTS.items():
        if key in _RESOLVED_LATER:
            continue
        got = parser.get_default(key)
        if got is None and not any(a.dest == key for a in parser._actions):
            continue  # config-only key, no CLI flag
        if got != value:
            wrong[key] = (got, value)
    assert not wrong, f"parser default differs from DEFAULTS: {wrong}"


def test_runconfig_field_defaults_equal_defaults_table():
    fields = {f.name: f for f in dataclasses.fields(RunConfig)}
    wrong = {}
    for key, value in DEFAULTS.items():
        field = fields.get(key)
        if field is None or field.default is dataclasses.MISSING:
            continue
        if field.default != value:
            wrong[key] = (field.default, value)
    assert not wrong, f"RunConfig default differs from DEFAULTS: {wrong}"


def test_multimap_method_default_is_canonical():
    assert _run_parser().get_default("multimap_method") == DEFAULTS["multimap_method"]
    assert DEFAULTS["multimap_method"] == "host-conservative"


@pytest.mark.parametrize(
    "module_name, banned",
    [
        ("viralscan.scripts.gene_programs", ('"catalogue")', ", 0.01)")),
        ("viralscan.scripts.cellcalling", ('"knee_min_umi", 10.0',)),
    ],
)
def test_scripts_do_not_hold_a_second_copy(module_name, banned):
    import importlib

    source = inspect.getsource(importlib.import_module(module_name))
    for text in banned:
        assert text not in source, f"{module_name} still hard-codes {text!r}"

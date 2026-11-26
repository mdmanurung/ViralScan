"""Every Snakemake worker exposes one ``run`` and keeps no module-level run state.

``WORKER_SHAPE_DEFERRED`` lists the four workers that still read module globals;
shrink it as each one is converted (PLAN ``WORKER-SHAPE-01``).
"""

import ast
import importlib
import inspect

import pytest

from viralscan.run_context import RunContext

CONVERTED = [
    "analysis",
    "anello_align",
    "auto_evidence",
    "host_filter",
    "hostresponse",
    "read_filter",
    "scautoqc",
    "star_confirm",
]
WORKER_SHAPE_DEFERRED = ["detection", "gene_programs", "multimap", "umap"]


@pytest.mark.parametrize("name", CONVERTED + WORKER_SHAPE_DEFERRED)
def test_worker_has_run_taking_a_context(name: str) -> None:
    module = importlib.import_module(f"viralscan.scripts.{name}")
    first = next(iter(inspect.signature(module.run).parameters.values()))
    assert first.name == "ctx"
    assert RunContext  # the context type every worker receives


@pytest.mark.parametrize("name", CONVERTED)
def test_converted_worker_has_no_global_statement(name: str) -> None:
    module = importlib.import_module(f"viralscan.scripts.{name}")
    tree = ast.parse(inspect.getsource(module))
    assert not [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Global)]

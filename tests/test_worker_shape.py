"""Every Snakemake worker exposes one ``run`` and keeps no module-level run state.

Run state (config, output, kb, identity, file) flows through ``ctx`` and function
arguments (PLAN ``WORKER-SHAPE-01``).
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
    "detection",
    "gene_programs",
    "host_filter",
    "hostresponse",
    "multimap",
    "read_filter",
    "scautoqc",
    "star_confirm",
    "umap",
]
RUN_STATE_NAMES = {"config", "output", "kb", "identity", "file"}


@pytest.mark.parametrize("name", CONVERTED)
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


@pytest.mark.parametrize("name", CONVERTED)
def test_converted_worker_assigns_no_run_state_at_module_level(name: str) -> None:
    """No module-level ``config = ...``; the Snakemake entry block only calls ``run``."""
    module = importlib.import_module(f"viralscan.scripts.{name}")
    tree = ast.parse(inspect.getsource(module))
    assigned = set()
    for node in tree.body:
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign)) and getattr(node, "value", None):
            targets = [node.target]
        else:
            continue
        assigned |= {t.id for t in targets if isinstance(t, ast.Name)}
    assert not assigned & RUN_STATE_NAMES
    assert not [n for n in tree.body if isinstance(n, ast.If) and _assigns(n, RUN_STATE_NAMES)]


def _assigns(node: ast.AST, names: set[str]) -> bool:
    return any(
        isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store) and n.id in names
        for n in ast.walk(node)
    )

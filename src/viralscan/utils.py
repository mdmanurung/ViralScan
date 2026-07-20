"""Small shared helpers used across the Snakemake worker scripts."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Union, cast

import numpy as np
import yaml

# Module-level logger for scripts that import this module.
logger = logging.getLogger("viralscan")


def configure_logging(verbose: bool = False, quiet: bool = False) -> None:
    """Configure the root *viralscan* logger.

    Called once from ``menu.py`` after parsing CLI args.  Snakemake worker
    scripts run in their own processes so they call this themselves from
    the ``setup_script_logging()`` helper below.

    Priority: quiet > verbose > default (INFO).
    """
    level = logging.WARNING if quiet else (logging.DEBUG if verbose else logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(levelname)s  %(message)s"))
    root = logging.getLogger("viralscan")
    root.setLevel(level)
    root.handlers.clear()
    root.addHandler(handler)


def setup_script_logging() -> logging.Logger:
    """Minimal logging setup for Snakemake worker scripts.

    Each script runs in its own interpreter so it must configure its own
    handler.  Returns the ``viralscan`` logger ready to use.
    """
    configure_logging()
    return logger


def split_comma_paths(value: str | None) -> list[str]:
    """Split a comma-separated CLI/config value into non-empty trimmed paths."""
    if value is None:
        return []
    return [part.strip() for part in value.split(",") if part.strip()]


def matrix_for_genes(adata: Any, matrix: Any, gene_names: list[str]) -> Any:
    """Return columns from *matrix* matching AnnData ``var_names`` gene names.

    AnnData slicing couples the selected matrix to ``adata.X``. Detection needs a
    different interface: callers may pass an alternate primary-call matrix
    (for example ``counts_unique_viral``) while still using the same AnnData
    metadata. This helper centralises the name-to-column mapping so statistics,
    enrichment, and plots cannot drift.
    """
    if not gene_names:
        return matrix[:, []]
    indices = adata.var_names.get_indexer(gene_names)
    missing = [gene for gene, idx in zip(gene_names, indices) if idx < 0]
    if missing:
        raise KeyError(f"Genes not present in AnnData var_names: {missing}")
    return matrix[:, np.asarray(indices, dtype=int)]


def resolve_count_matrix(viral_count_matrix: Any, adata: Any) -> Any:
    """Return the primary-call matrix, defaulting to ``adata.X``.

    Centralises the ``viral_count_matrix if ... is not None else adata.X`` guard
    shared by Detection, enrichment, and plotting so a new primary-call mode is
    wired in one place instead of five.
    """
    return viral_count_matrix if viral_count_matrix is not None else adata.X


def load_config(path: Union[str, Path]) -> dict[str, Any]:
    """Read a YAML config file and return it as a plain ``dict``.

    Centralised so every script can use the same loader (and so we have
    one place to evolve the boolean-normalisation work tracked in
    PLAN §1.6).
    """
    with open(path, encoding="utf-8") as f:
        loaded = cast(Any, yaml.safe_load(f))
    if not isinstance(loaded, dict):
        raise ValueError(f"Config file {path} did not contain a YAML mapping.")
    return cast(dict[str, Any], loaded)

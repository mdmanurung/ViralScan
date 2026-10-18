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
    if not adata.var_names.is_unique:
        raise ValueError(
            "matrix_for_genes requires unique adata.var_names — the reference has "
            "duplicate gene IDs. Rebuild the reference with unique gene_ids."
        )
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


#: Relative slack on the infected cut-off. An allocation that splits one molecule
#: over several genes or cells (1/3 + 1/3 + 1/3 in float32 is 0.99999994) must not
#: fall just under a cut-off of 1 and be demoted to "sub-threshold".
INFECTED_RELATIVE_TOLERANCE = 1e-6

#: Private byte budget for the dense row-block reduction in
#: :func:`virus_molecules_per_cell`. Every selected column is kept in its original
#: order and each row keeps its full width, so the per-row float64 reduction is
#: unchanged while a whole-matrix dense copy is avoided. Only whole rows are
#: densified, so a single row wider than the budget is still densified (and
#: reduced) on its own; a row's reduction is never split across blocks. Read at
#: call time so the qualification harness can pin reduced budgets.
_ROW_BLOCK_BYTES = 32 * 1024 * 1024


def _sparse_row_block(n_rows: int, n_cols: int) -> int | None:
    """Rows per dense block, or ``None`` when the whole matrix fits the budget.

    Returns ``None`` when there are no columns, when nothing needs blocking, or
    when the whole matrix already fits. A lone oversized row still resolves to
    ``1`` so it is densified alone rather than split.
    """
    if n_cols <= 0:
        return None
    rows = max(1, _ROW_BLOCK_BYTES // (n_cols * np.dtype(float).itemsize))
    if rows >= n_rows:
        return None
    return rows


def virus_molecules_per_cell(adata: Any, matrix: Any, gene_ids: list[str]) -> Any:
    """Per-cell molecules summed over *gene_ids*: a dense 1-D float array.

    *matrix* is the primary-call (selected-method) matrix; fractional weights are kept.
    """
    sub = matrix_for_genes(adata, matrix, [g for g in gene_ids if g in adata.var_names])
    if hasattr(sub, "toarray"):
        # Only a CSR densification is provably C-order, so each blocked row is the
        # same contiguous width-N float64 slice the unblocked ``sum(axis=1)`` would
        # reduce, in the same element order. Dense views, CSC and every other layout
        # keep the original single-shot reduction rather than risk a
        # layout-sensitive bit change.
        # Only real dtypes whose per-block ``toarray`` + float64 cast can never warn
        # (bool/int/uint/float32/64; scipy.sparse has no float16: exact or correctly-rounded, no cast warning,
        # no overflow). Complex (ComplexWarning on every cast, so once per block),
        # longdouble (cast can overflow/warn per block), object and anything else
        # take the original single-shot path so warning multiplicity is baseline's.
        if (
            getattr(sub, "format", None) == "csr"
            and sub.dtype.kind in "biuf"
            and sub.dtype.itemsize <= 8
        ):
            n_rows, n_cols = int(sub.shape[0]), int(sub.shape[1])
            block = _sparse_row_block(n_rows, n_cols)
            if block is not None:
                out = np.empty(n_rows, dtype=float)
                # Block-level floating-point warnings/errors are suppressed: if any
                # row sum is nonfinite (NaN/Inf input or overflow) the original
                # single-shot path below reruns, so warning category, message,
                # multiplicity and raised errors are exactly the original's.
                with np.errstate(over="ignore", invalid="ignore"):
                    for start in range(0, n_rows, block):
                        stop = min(start + block, n_rows)
                        dense = np.asarray(sub[start:stop].toarray(), dtype=float)
                        out[start:stop] = dense.sum(axis=1)
                if np.isfinite(out).all():
                    return out.reshape(-1)
        sub = sub.toarray()
    return np.asarray(sub, dtype=float).sum(axis=1).reshape(-1)


def infected_masks(per_cell: Any, min_molecules: float) -> tuple[Any, Any]:
    """The one "infected" rule: ``(infected, subthreshold)`` boolean masks over cells.

    A cell is infected for a virus when its summed allocated molecules reach
    *min_molecules* (``--infected-min-molecules``). A cell with ``0 < weight <
    min_molecules`` is sub-threshold: still reported (``infected_subthreshold_*``),
    never silently dropped. Detection, cell-type enrichment and the UMAP all use this.
    """
    weight = np.asarray(per_cell, dtype=float)
    infected = weight >= float(min_molecules) * (1.0 - INFECTED_RELATIVE_TOLERANCE)
    return infected, (weight > 0) & ~infected


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

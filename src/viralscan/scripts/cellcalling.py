"""Cell-calling for ViralScan — label real (non-empty-droplet) barcodes.

Reporting viral detection rates over *all* barcodes (empty droplets included)
produces meaningless "% infected" denominators. This is the recurring bug behind
the HSV-1 apparent 25x discrepancy (P22.5) and the covid empty-droplet artifact
(finding F-005): e.g. Torque teno virus reads as 37.6% of all barcodes but ~90%
of *called cells* at >=5 UMI. This module labels which barcodes are real cells so
downstream stats can be reported over them (while still also reporting the
all-barcode denominator, so the choice is never hidden).

Methods (config ``cell_calling``):
  - ``auto``      : use an external list when supplied, otherwise emptyDrops.
  - ``external``  : use a provided barcode list (e.g. CellRanger / STARsolo called
                    cells). PREFERRED when a matched run exists — it handles the
                    chemistry/barcode space correctly and calls cells with full
                    annotation coverage.
  - ``emptydrops``: DropletUtils::emptyDrops via ``emptydrops.R`` (gold standard;
                    needs R + DropletUtils on ``cell_caller_rscript``).
  - ``knee``      : explicit pure-Python barcode-rank knee approximation.
  - ``none``      : every barcode treated as a cell (legacy behaviour).

All callers return a boolean mask aligned to ``obs_names``.
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path

import numpy as np

from viralscan.defaults import DEFAULTS

log = logging.getLogger("viralscan")


def _strip_suffix(bc: str) -> str:
    """Drop a trailing 10x ``-1`` (or ``-N``) gem-group suffix if present."""
    i = bc.rfind("-")
    if i != -1 and bc[i + 1 :].isdigit():
        return bc[:i]
    return bc


class CellCallingError(RuntimeError):
    """Raised when cell calling cannot produce a trustworthy mask.

    Cell calling sets the denominator for every reported viral rate. A failure
    that silently falls back to "every barcode is a cell" does not lose the
    result, it changes what the result means, because barcodes are mostly empty
    droplets. So every failure here is fatal, and the only way to report over all
    barcodes is to ask for it with cell_calling=none.
    """


def external_cells(obs_names, barcode_file: os.PathLike[str] | str) -> np.ndarray:
    """Mask of *obs_names* present in an external called-cell list.

    The list may be plain text or ``.gz``, one barcode per line, with or without
    a ``-1`` suffix (both the list and obs_names are suffix-stripped before match).
    """
    import gzip

    path = Path(barcode_file)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        wanted = {_strip_suffix(line.strip()) for line in fh if line.strip()}
    if not wanted:
        raise CellCallingError(f"cell_calling=external: the called-cell list {path} is empty")

    canonical = [_strip_suffix(str(b)) for b in obs_names]
    duplicates = len(canonical) - len(set(canonical))
    if duplicates:
        # Two raw barcodes collapsing to one canonical form makes membership
        # ambiguous, and the ambiguity is invisible in the resulting mask.
        raise CellCallingError(
            f"cell_calling=external: {duplicates} barcodes collide after suffix "
            "stripping, so external membership is ambiguous"
        )

    mask = np.array([bc in wanted for bc in canonical], dtype=bool)
    log.info(
        "cell_calling=external: %d/%d barcodes matched the called-cell list",
        int(mask.sum()),
        len(mask),
    )
    if mask.sum() == 0:
        raise CellCallingError(
            f"cell_calling=external: none of {len(mask)} barcodes matched the "
            f"{len(wanted)} in {path}. The barcode spaces almost certainly differ "
            "in orientation or translation. Continuing would report viral rates "
            "over all barcodes while labelling them called-cell rates."
        )
    return mask


def knee_cells(total_umi, min_umi: float = 10.0) -> np.ndarray:
    """Barcode-rank knee via a kneedle-style max-distance-from-chord heuristic.

    On the log10(rank) vs log10(total-UMI) curve (barcodes with total >= *min_umi*,
    sorted descending), draw the chord from the first to the last point and take the
    knee as the point of maximum perpendicular distance below that chord. Barcodes
    with total UMI at or above the knee value are cells.

    This is an approximation for when no external list / emptyDrops is available;
    prefer ``external`` or ``emptydrops`` for publication-grade calls.
    """
    total = np.asarray(total_umi, dtype=float)
    order = np.argsort(total)[::-1]
    kept = total[order] >= min_umi
    if kept.sum() < 3:
        # too few barcodes to find a knee — treat everything above min_umi as cells
        return total >= max(min_umi, 1.0)

    y = np.log10(total[order][kept])
    x = np.log10(np.arange(1, kept.sum() + 1))
    # perpendicular distance of each point from the chord (x0,y0)->(x1,y1)
    x0, y0, x1, y1 = x[0], y[0], x[-1], y[-1]
    dx, dy = x1 - x0, y1 - y0
    denom = np.hypot(dx, dy) or 1.0
    dist = ((y - y0) * dx - (x - x0) * dy) / denom  # signed; below chord is negative
    knee_i = int(np.argmin(dist))  # most-below-chord point
    knee_val = float(10 ** y[knee_i])
    mask = total >= knee_val
    log.info(
        "cell_calling=knee: knee at total>=%.0f -> %d/%d cells",
        knee_val,
        int(mask.sum()),
        len(mask),
    )
    return mask


def emptydrops_cells(
    obs_names, matrix_dir, *, rscript, fdr, lower, niters, seed
) -> np.ndarray:
    """Mask from DropletUtils::emptyDrops via the bundled ``emptydrops.R``.

    Every parameter is required and keyword-only. emptyDrops is a Monte-Carlo
    test, so ``seed`` and ``niters`` change which barcodes are called; a default
    here would let a run silently use a value that no configuration declared,
    which is how the protocol-frozen ``seeds.cell_calling`` came to be ignored.
    """
    script = Path(__file__).with_name("emptydrops.R")
    out_tsv = Path(matrix_dir) / "emptydrops_cells.tsv"
    cmd = [
        rscript,
        str(script),
        str(matrix_dir),
        str(out_tsv),
        str(fdr),
        str(lower),
        str(niters),
        str(seed),
    ]
    log.info("cell_calling=emptydrops: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)  # list form, no shell (CLAUDE.md §1.2)

    called: set[str] = set()
    with open(out_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        bc_i, cell_i = header.index("barcode"), header.index("is_cell")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[cell_i].strip() in ("TRUE", "True", "1"):
                called.add(parts[bc_i])
    mask = np.array([str(b) in called for b in obs_names], dtype=bool)
    if not mask.any():
        # external_cells raises on zero matches; this path did not, so the SW-11
        # promise that every caller fails closed was not quite true. A zero-cell
        # mask makes every called-cell rate a 0/0, reported as 0.0 rather than as
        # a failure.
        raise CellCallingError(
            f"cell_calling=emptydrops called zero cells from {out_tsv}. Check the "
            "matrix depth, lower, and FDR, or rerun with --cell-calling none to "
            "report over all barcodes deliberately."
        )
    log.info("cell_calling=emptydrops: %d/%d cells", int(mask.sum()), len(mask))
    return mask


def call_cells(adata, config, matrix_dir=None) -> np.ndarray:
    """Return a boolean cell mask over ``adata.obs_names`` per ``config.cell_calling``.

    Parameters
    ----------
    matrix_dir : path | None
        kb ``counts_unfiltered`` directory — required for ``emptydrops`` (its ``.mtx``
        is what DropletUtils reads). Ignored by the other methods.

    Recognised config attributes (all optional, falling back to ``DEFAULTS``):
      cell_calling        : auto|external|emptydrops|knee|none   (default: auto)
      called_cells_file   : path (required for external)
      cell_caller_rscript : Rscript path (default: "Rscript")
      emptydrops_fdr/emptydrops_lower/emptydrops_niters/emptydrops_seed
      knee_min_umi

    ``emptydrops_seed`` and ``emptydrops_niters`` govern a Monte-Carlo test, so
    they change which barcodes are called. Both come from the configuration; a
    run under a frozen protocol sets ``emptydrops_seed`` from that protocol's
    ``seeds.cell_calling``.
    """
    method = str(getattr(config, "cell_calling", "auto") or "auto").lower()
    obs = adata.obs_names

    if len(obs) == 0:
        raise CellCallingError("cell calling requires at least one barcode")

    if method == "auto":
        method = "external" if getattr(config, "called_cells_file", None) else "emptydrops"
        log.info("cell_calling=auto selected %s", method)

    if method == "none":
        log.info("cell_calling=none: all %d barcodes treated as cells", adata.n_obs)
        return np.ones(adata.n_obs, dtype=bool)

    if method == "external":
        f = getattr(config, "called_cells_file", None)
        if not f:
            raise CellCallingError("cell_calling=external requires config.called_cells_file")
        return external_cells(obs, f)

    if method == "emptydrops":
        mdir = matrix_dir or getattr(config, "cell_caller_matrix_dir", None)
        if not mdir:
            raise CellCallingError(
                "cell_calling=emptydrops requires the kb counts_unfiltered "
                "directory (pass matrix_dir=...)"
            )
        return emptydrops_cells(
            obs,
            mdir,
            rscript=getattr(config, "cell_caller_rscript", DEFAULTS["cell_caller_rscript"]),
            fdr=float(getattr(config, "emptydrops_fdr", DEFAULTS["emptydrops_fdr"])),
            lower=float(getattr(config, "emptydrops_lower", DEFAULTS["emptydrops_lower"])),
            niters=int(getattr(config, "emptydrops_niters", DEFAULTS["emptydrops_niters"])),
            seed=int(getattr(config, "emptydrops_seed", DEFAULTS["emptydrops_seed"])),
        )

    if method == "knee" and hasattr(adata.X, "sum"):
        import scipy.sparse as sp

        total = (
            np.asarray(adata.X.sum(axis=1)).ravel() if sp.issparse(adata.X) else adata.X.sum(axis=1)
        )
    elif method == "knee":
        total = np.asarray(adata.X).sum(axis=1)
    else:
        raise CellCallingError(f"Unknown cell_calling method: {method!r}")
    return knee_cells(total, min_umi=float(getattr(config, "knee_min_umi", 10.0)))

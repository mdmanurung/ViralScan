"""Cell-calling for ViralScan — label real (non-empty-droplet) barcodes.

Reporting viral detection rates over *all* barcodes (empty droplets included)
produces meaningless "% infected" denominators. This is the recurring bug behind
the HSV-1 apparent 25x discrepancy (P22.5) and the covid empty-droplet artifact
(finding F-005): e.g. Torque teno virus reads as 37.6% of all barcodes but ~90%
of *called cells* at >=5 UMI. This module labels which barcodes are real cells so
downstream stats can be reported over them (while still also reporting the
all-barcode denominator, so the choice is never hidden).

Methods (config ``cell_calling``):
  - ``external``  : use a provided barcode list (e.g. CellRanger / STARsolo called
                    cells). PREFERRED when a matched run exists — it handles the
                    chemistry/barcode space correctly and calls cells with full
                    annotation coverage.
  - ``emptydrops``: DropletUtils::emptyDrops via ``emptydrops.R`` (gold standard;
                    needs R + DropletUtils on ``cell_caller_rscript``).
  - ``knee``      : pure-Python barcode-rank knee (dependency-free default).
  - ``none``      : every barcode treated as a cell (legacy behaviour).

All callers return a boolean mask aligned to ``obs_names``.
"""
from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path

import numpy as np

log = logging.getLogger("viralscan")


def _strip_suffix(bc: str) -> str:
    """Drop a trailing 10x ``-1`` (or ``-N``) gem-group suffix if present."""
    i = bc.rfind("-")
    if i != -1 and bc[i + 1:].isdigit():
        return bc[:i]
    return bc


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
    mask = np.array([_strip_suffix(str(b)) in wanted for b in obs_names], dtype=bool)
    log.info("cell_calling=external: %d/%d barcodes matched the called-cell list",
             int(mask.sum()), len(mask))
    if mask.sum() == 0:
        log.warning("cell_calling=external matched 0 cells — barcode spaces may differ "
                    "(orientation/translation). Falling back is the caller's decision.")
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
    knee_i = int(np.argmin(dist))                   # most-below-chord point
    knee_val = float(10 ** y[knee_i])
    mask = total >= knee_val
    log.info("cell_calling=knee: knee at total>=%.0f -> %d/%d cells",
             knee_val, int(mask.sum()), len(mask))
    return mask


def emptydrops_cells(obs_names, matrix_dir, rscript="Rscript", fdr=0.01,
                     lower=100, niters=10000, seed=100) -> np.ndarray:
    """Mask from DropletUtils::emptyDrops via the bundled ``emptydrops.R``."""
    script = Path(__file__).with_name("emptydrops.R")
    out_tsv = Path(matrix_dir) / "emptydrops_cells.tsv"
    cmd = [rscript, str(script), str(matrix_dir), str(out_tsv),
           str(fdr), str(lower), str(niters), str(seed)]
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
    log.info("cell_calling=emptydrops: %d/%d cells", int(mask.sum()), len(mask))
    return mask


def call_cells(adata, config) -> np.ndarray:
    """Return a boolean cell mask over ``adata.obs_names`` per ``config.cell_calling``.

    Recognised config attributes (all optional, with sensible defaults):
      cell_calling            : external|emptydrops|knee|none   (default: knee)
      called_cells_file       : path (required for external)
      cell_caller_rscript     : Rscript path (default: "Rscript")
      emptydrops_fdr/lower/niters, cell_caller_matrix_dir, knee_min_umi
    """
    method = str(getattr(config, "cell_calling", "knee") or "knee").lower()
    obs = adata.obs_names

    if method == "none":
        log.info("cell_calling=none: all %d barcodes treated as cells", adata.n_obs)
        return np.ones(adata.n_obs, dtype=bool)

    if method == "external":
        f = getattr(config, "called_cells_file", None)
        if not f:
            raise ValueError("cell_calling=external requires config.called_cells_file")
        return external_cells(obs, f)

    if method == "emptydrops":
        mdir = getattr(config, "cell_caller_matrix_dir", None)
        if not mdir:
            raise ValueError("cell_calling=emptydrops requires config.cell_caller_matrix_dir "
                             "(the kb counts_unfiltered directory)")
        return emptydrops_cells(
            obs, mdir,
            rscript=getattr(config, "cell_caller_rscript", "Rscript"),
            fdr=float(getattr(config, "emptydrops_fdr", 0.01)),
            lower=float(getattr(config, "emptydrops_lower", 100)),
            niters=int(getattr(config, "emptydrops_niters", 10000)),
        )

    # default: knee
    if hasattr(adata.X, "sum"):
        import scipy.sparse as sp
        total = (np.asarray(adata.X.sum(axis=1)).ravel()
                 if sp.issparse(adata.X) else adata.X.sum(axis=1))
    else:
        total = np.asarray(adata.X).sum(axis=1)
    return knee_cells(total, min_umi=float(getattr(config, "knee_min_umi", 10.0)))

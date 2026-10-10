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
                    Sensitivity-only, never for reported numbers (decided
                    2026-10-01): its estimator lands at ``knee_min_umi`` on
                    real libraries (PLAN SW-23), so it calls empty droplets
                    as cells. Every knee run warns.
  - ``none``      : every barcode treated as a cell (legacy behaviour).

All callers return a boolean mask aligned to ``obs_names``.
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
from pathlib import Path
from typing import TypedDict

import numpy as np

from viralscan.defaults import DEFAULTS

#: Total-UMI floor for the strategy-independent denominator, in host molecules
#: per barcode. Chosen to sit above the empty-droplet mode and below the knee of
#: a real 10x barcode-rank curve, so it selects cells under *any* host-filter
#: strategy. Deliberately not ``defaults.min_counts`` (1000): that is a UMAP QC
#: knob, and after host subtraction most barcodes fall below it.
COMPARABLE_CELL_MIN_UMI = 200.0

log = logging.getLogger("viralscan")


class HostCellSets(TypedDict):
    called: set[str]
    comparable: set[str] | None


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
    path = Path(barcode_file)
    wanted = _external_barcodes(path)

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


def _external_barcodes(path: Path) -> set[str]:
    import gzip

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        wanted = {_strip_suffix(line.strip()) for line in fh if line.strip()}
    if not wanted:
        raise CellCallingError(f"cell_calling=external: the called-cell list {path} is empty")
    return wanted


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


def _run_emptydrops(matrix_dir, *, rscript, fdr, lower, niters, seed) -> Path:
    """Run ``emptydrops.R`` on *matrix_dir*; returns the TSV it wrote."""
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
    return out_tsv


def emptydrops_cells(obs_names, matrix_dir, *, rscript, fdr, lower, niters, seed) -> np.ndarray:
    """Mask from DropletUtils::emptyDrops via the bundled ``emptydrops.R``.

    Every parameter is required and keyword-only. emptyDrops is a Monte-Carlo
    test, so ``seed`` and ``niters`` change which barcodes are called; a default
    here would let a run silently use a value that no configuration declared,
    which is how the protocol-frozen ``seeds.cell_calling`` came to be ignored.
    """
    out_tsv = _run_emptydrops(
        matrix_dir, rscript=rscript, fdr=fdr, lower=lower, niters=niters, seed=seed
    )
    mask = read_emptydrops_mask(obs_names, out_tsv)
    log.info("cell_calling=emptydrops: %d/%d cells", int(mask.sum()), len(mask))
    return mask


def solo_raw_dir(config) -> Path | None:
    """Required GeneFull host matrix for two-step EmptyDrops, else ``None``.

    A two-step run quantifies only the reads STAR could not place on the host,
    so its kb matrix holds no host UMIs and emptyDrops has nothing to separate
    cells from ambient on (1-2,659 barcodes, "insufficient unique points").
    The host matrix STAR already wrote is the one that can.
    """
    if not getattr(config, "host_index", None) or resolve_method(config) != "emptydrops":
        return None
    raw = Path(config.output) / "host_filtered" / "star_tmp" / "Solo.out" / "GeneFull" / "raw"
    missing = [
        name
        for name in ("matrix.mtx", "barcodes.tsv", "features.tsv")
        if not (raw / name).is_file()
    ]
    if missing:
        raise CellCallingError(
            f"two-step EmptyDrops requires STARsolo GeneFull/raw in {raw}; missing "
            f"{', '.join(missing)}. Gene/raw and the viral kb matrix cannot replace it. "
            "Supply --called-cells-file with --cell-calling auto/external or regenerate "
            "the host matrix with --soloFeatures GeneFull."
        )
    return raw


def read_host_matrix_totals(raw: Path) -> tuple[list[str], np.ndarray]:
    """Validated STARsolo ``GeneFull/raw`` barcodes and per-barcode host UMI totals.

    Checks dimensions, duplicate barcodes and non-finite or negative counts, so a
    corrupt matrix stops here instead of inside R or in the comparable-cell set.
    """
    from scipy.io import mmread

    barcodes = [_strip_suffix(b) for b in (raw / "barcodes.tsv").read_text().splitlines()]
    matrix = mmread(raw / "matrix.mtx")
    n_features = len((raw / "features.tsv").read_text().splitlines())
    if matrix.shape != (n_features, len(barcodes)) or len(set(barcodes)) != len(barcodes):
        raise CellCallingError(
            f"invalid GeneFull/raw barcode dimension or duplicate barcodes in {raw}"
        )
    values = matrix.data if hasattr(matrix, "tocoo") else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any():
        raise CellCallingError(f"invalid GeneFull/raw host counts in {raw}")
    return barcodes, np.asarray(matrix.sum(axis=0)).ravel()


def host_cells_from_tsv(out_tsv, min_comparable_umi: float) -> HostCellSets:
    """Called barcodes, and the subset with host UMI >= *min_comparable_umi*."""
    called: set[str] = set()
    comparable: set[str] = set()
    with open(out_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        bc_i, cell_i, tot_i = (header.index(c) for c in ("barcode", "is_cell", "total"))
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[cell_i].strip() in ("TRUE", "True", "1"):
                called.add(parts[bc_i])
                if float(parts[tot_i]) >= min_comparable_umi:
                    comparable.add(parts[bc_i])
    if not called:
        raise CellCallingError(
            f"cell_calling=emptydrops called zero cells from the host matrix in {out_tsv}."
        )
    return {"called": called, "comparable": comparable}


def host_called_cells(config, solo_dir, min_comparable_umi: float) -> HostCellSets:
    """emptyDrops on the STARsolo host matrix: the cell set of a two-step run."""
    read_host_matrix_totals(Path(solo_dir))  # fail on a corrupt matrix before starting R
    out_tsv = _run_emptydrops(
        solo_dir,
        rscript=getattr(config, "cell_caller_rscript", DEFAULTS["cell_caller_rscript"]),
        fdr=float(getattr(config, "emptydrops_fdr", DEFAULTS["emptydrops_fdr"])),
        lower=float(getattr(config, "emptydrops_lower", DEFAULTS["emptydrops_lower"])),
        niters=int(getattr(config, "emptydrops_niters", DEFAULTS["emptydrops_niters"])),
        seed=int(getattr(config, "emptydrops_seed", DEFAULTS["emptydrops_seed"])),
    )
    cells = host_cells_from_tsv(out_tsv, min_comparable_umi)
    _write_input_receipt(
        config,
        {
            "method": "emptydrops",
            "count_layer": "GeneFull.raw",
            "matrix_dir": str(solo_dir),
            "input_sha256": _digests(
                Path(solo_dir), ("matrix.mtx", "barcodes.tsv", "features.tsv")
            ),
            "emptydrops_cells_sha256": _digests(Path(out_tsv).parent, (Path(out_tsv).name,)),
            "parameters": {
                "fdr": float(getattr(config, "emptydrops_fdr", DEFAULTS["emptydrops_fdr"])),
                "lower": float(getattr(config, "emptydrops_lower", DEFAULTS["emptydrops_lower"])),
                "niters": int(getattr(config, "emptydrops_niters", DEFAULTS["emptydrops_niters"])),
                "seed": int(getattr(config, "emptydrops_seed", DEFAULTS["emptydrops_seed"])),
            },
        },
    )
    comparable = cells["comparable"]
    assert comparable is not None
    log.info(
        "cell_calling=emptydrops on host matrix: %d cells (%d with host UMI >= %g)",
        len(cells["called"]),
        len(comparable),
        min_comparable_umi,
    )
    return cells


def external_host_cells(config, min_comparable_umi: float) -> HostCellSets:
    """Full external universe; host-depth subset when GeneFull is available."""
    called = _external_barcodes(Path(config.called_cells_file))
    raw = Path(config.output) / "host_filtered" / "star_tmp" / "Solo.out" / "GeneFull" / "raw"
    comparable = None
    complete = all(
        (raw / name).is_file() for name in ("matrix.mtx", "barcodes.tsv", "features.tsv")
    )
    if complete:
        barcodes, totals = read_host_matrix_totals(raw)
        comparable = {
            b for b, total in zip(barcodes, totals) if b in called and total >= min_comparable_umi
        }
    _write_input_receipt(
        config,
        {
            "method": "external",
            "called_cells_file": str(config.called_cells_file),
            "called_cells_sha256": _digests(
                Path(config.called_cells_file).parent, (Path(config.called_cells_file).name,)
            ),
            "input_sha256": (
                _digests(raw, ("matrix.mtx", "barcodes.tsv", "features.tsv")) if complete else None
            ),
            "matrix_dir": str(raw) if complete else None,
            "count_layer": "GeneFull.raw" if complete else None,
            "comparable_status": "available" if complete else "host_matrix_unavailable",
        },
    )
    return {"called": called, "comparable": comparable}


def _digests(directory: Path, names) -> dict[str, str]:
    """SHA-256 of each named file, so the receipt pins the exact cell-calling inputs."""
    from viralscan.run_safety import sha256_file

    return {name: sha256_file(directory / name) for name in names}


def _write_input_receipt(config, payload) -> None:
    results = Path(config.output) / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "cell_calling_input.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )


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
    method = resolve_method(config)
    obs = adata.obs_names

    if len(obs) == 0:
        raise CellCallingError("cell calling requires at least one barcode")

    if str(getattr(config, "cell_calling", "auto") or "auto").lower() == "auto":
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

    if method == "knee":
        log.warning(
            "cell_calling=knee is sensitivity-only and not for reported numbers: "
            "its estimator lands near knee_min_umi on real libraries and calls "
            "empty droplets as cells (PLAN SW-23). Use emptydrops or external."
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


def read_emptydrops_mask(obs_names, out_tsv) -> np.ndarray:
    """Mask of *obs_names* that ``emptydrops.R`` marked ``is_cell`` in *out_tsv*.

    Raises when no barcode is called, so every caller fails closed.
    """
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
    return mask


def resolve_method(config) -> str:
    """``config.cell_calling`` with ``auto`` resolved to the method it runs."""
    method = str(getattr(config, "cell_calling", "auto") or "auto").lower()
    if method == "auto":
        method = "external" if getattr(config, "called_cells_file", None) else "emptydrops"
    return method


#: The called-cell set detection used, one barcode per row (PLAN PROG-17).
CALLED_CELLS_TSV = os.path.join("results", "called_cells.tsv")


def write_called_cells(obs_names, mask, outputpath) -> str:
    """Write the called barcodes to ``results/called_cells.tsv``."""
    path = os.path.join(outputpath, CALLED_CELLS_TSV)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("barcode\n")
        fh.writelines(f"{b}\n" for b in np.asarray(obs_names)[np.asarray(mask, dtype=bool)])
    return path


def load_called_mask(adata, config, run_dir) -> np.ndarray:
    """The called-cell mask detection used for *run_dir*, over ``adata.obs_names``.

    Reads ``results/called_cells.tsv``. A run directory from before PROG-17 has
    none, so it falls back to that run's own emptyDrops output when its method
    was emptyDrops, else re-calls cells with the run's configuration. Fails
    closed: a barcode list that does not match ``obs_names`` means another
    matrix, and scoring over it would read as a real negative.
    """
    path = Path(run_dir) / CALLED_CELLS_TSV
    if path.is_file():
        with open(path, encoding="utf-8") as fh:
            fh.readline()
            called = {line.strip() for line in fh if line.strip()}
        mask = np.isin(np.asarray(adata.obs_names, dtype=str), list(called))
        if not called or int(mask.sum()) != len(called):
            raise CellCallingError(
                f"{path} lists {len(called)} called barcode(s) but only "
                f"{int(mask.sum())} are in this matrix; it belongs to another run."
            )
        log.info("called cells: %d/%d from %s", len(called), adata.n_obs, path)
        return mask
    counts_dir = Path(run_dir) / "kb-python" / "counts_unfiltered"
    solo_dir = solo_raw_dir(config)
    if solo_dir is not None:
        host_tsv = solo_dir / "emptydrops_cells.tsv"
        if host_tsv.is_file():
            return read_emptydrops_mask(adata.obs_names, host_tsv)
        host_called_cells(config, solo_dir, COMPARABLE_CELL_MIN_UMI)
        return read_emptydrops_mask(adata.obs_names, host_tsv)
    legacy = counts_dir / "emptydrops_cells.tsv"
    if resolve_method(config) == "emptydrops" and legacy.is_file():
        mask = read_emptydrops_mask(adata.obs_names, legacy)
        log.info(
            "called cells: %d/%d from %s (pre-PROG-17 run)", int(mask.sum()), adata.n_obs, legacy
        )
        return mask
    log.info("called cells: no %s; re-calling cells with the run's configuration", path)
    return call_cells(adata, config, matrix_dir=counts_dir)

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
                    needs R + DropletUtils on ``cell_caller_rscript``). Always on
                    the HOST transcriptome: STARsolo GeneFull in two-step mode, the
                    host columns of the kb matrix in combined mode (never viral reads).
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
import shutil
import subprocess
import tempfile
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

#: R's ``.Machine$integer.max``; ``set.seed`` takes nothing larger.
R_INTEGER_MAX = 2**31 - 1

log = logging.getLogger("viralscan")


def effective_seed(seed: int) -> int:
    """The seed ``emptydrops.R`` hands to ``set.seed`` for *seed*.

    Seeds within R's integer range are used unchanged; a larger one (the frozen
    ``seeds.cell_calling`` 20260727002 is) is reduced modulo ``R_INTEGER_MAX``.
    Keep in step with ``emptydrops.R``.
    """
    seed = int(seed)
    return seed if abs(seed) <= R_INTEGER_MAX else seed % R_INTEGER_MAX


class _HostCellSetsRequired(TypedDict):
    called: set[str]
    comparable: set[str] | None


class HostCellSets(_HostCellSetsRequired, total=False):
    #: Suffix-stripped barcode -> host molecules from STARsolo ``GeneFull/raw`` (zero totals
    #: omitted). Absent when no host matrix exists; the ambient screen is then not assessed.
    host_totals: dict[str, float]


def _nonzero_totals(barcodes, totals) -> dict[str, float]:
    return {b: float(t) for b, t in zip(barcodes, totals) if t > 0}


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


def _run_emptydrops(
    matrix_dir, *, rscript, fdr, lower, niters, seed, out_tsv: os.PathLike[str] | str | None = None
) -> Path:
    """Run ``emptydrops.R`` on *matrix_dir*; returns the TSV it wrote.

    The TSV defaults to ``<matrix_dir>/emptydrops_cells.tsv``.
    """
    script = Path(__file__).with_name("emptydrops.R")
    out_tsv = Path(out_tsv) if out_tsv else Path(matrix_dir) / "emptydrops_cells.tsv"
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
    # A failed re-run must not leave a previous (possibly all-column) table for a reader to reuse.
    out_tsv.unlink(missing_ok=True)
    subprocess.run(cmd, check=True)  # list form, no shell (CLAUDE.md §1.2)
    return out_tsv


def emptydrops_cells(
    obs_names, matrix_dir, *, rscript, fdr, lower, niters, seed, out_tsv=None
) -> np.ndarray:
    """Mask from DropletUtils::emptyDrops via the bundled ``emptydrops.R``.

    Every parameter is required and keyword-only. emptyDrops is a Monte-Carlo
    test, so ``seed`` and ``niters`` change which barcodes are called; a default
    here would let a run silently use a value that no configuration declared,
    which is how the protocol-frozen ``seeds.cell_calling`` came to be ignored.
    """
    out_tsv = _run_emptydrops(
        matrix_dir, rscript=rscript, fdr=fdr, lower=lower, niters=niters, seed=seed, out_tsv=out_tsv
    )
    mask = read_emptydrops_mask(obs_names, out_tsv)
    log.info("cell_calling=emptydrops: %d/%d cells", int(mask.sum()), len(mask))
    return mask


def viral_partition(run_output, matrix_genes=None) -> tuple[set[str], Path]:
    """The Run's viral gene IDs and the file they came from.

    Source: the Virus Identity table, else ``log/analysis.txt``. Unlike
    ``multimap.viral_gene_partition`` this fails closed: cell calling must never
    quietly fall back to every column of the kb matrix. With *matrix_genes* and a
    Virus Identity table, a matrix column the table does not list also fails: the
    table covers every column of a real run, so a gap means another run's table,
    and a silent host classification could change the cell set.
    """
    from viralscan.virus_grouping import (
        IdentityMatrixMismatch,
        load_run_identity,
        reconcile_identity_to_matrix,
    )

    table = load_run_identity(run_output)
    ids: set[str] = set()
    if table is not None:
        source = Path(run_output) / "results" / "virus_identity.tsv"
        ids = set(table.viral_gene_ids())
        if matrix_genes:
            try:
                reconcile_identity_to_matrix(table, matrix_genes, source)
            except IdentityMatrixMismatch as exc:
                raise CellCallingError(
                    f"{exc} Or pass --called-cells-file / use two-step --host-filter starsolo."
                ) from exc
    else:
        source = Path(run_output) / "log" / "analysis.txt"
        if source.is_file():
            ids = {line.strip() for line in source.read_text().splitlines() if line.strip()}
    if not ids:
        raise CellCallingError(
            f"cannot call cells on the host genes only: no viral gene partition in {run_output} "
            "(neither results/virus_identity.tsv nor log/analysis.txt lists viral genes). "
            "Refusing to call cells on every column, which would count viral reads. "
            "Re-run the analysis step, or pass --called-cells-file / use two-step "
            "--host-filter starsolo."
        )
    return ids, source


def write_host_only_matrix(matrix_dir, viral_genes, dest) -> tuple[Path, int, int]:
    """Write the host-gene columns of a kb ``cells_x_genes`` matrix into *dest*.

    Rows (barcodes) are untouched and the barcode file is byte-copied, so the
    result stays aligned with the original matrix. Viral columns are dropped, not
    zeroed. Fails closed when the partition does not match the matrix.
    Returns ``(dest, n_host_genes, n_viral_genes_excluded)``.
    """
    from scipy.io import mmread, mmwrite

    src, dest = Path(matrix_dir), Path(dest)
    genes = (src / "cells_x_genes.genes.txt").read_text().splitlines()
    matrix = mmread(src / "cells_x_genes.mtx").tocsc()
    n_barcodes = len((src / "cells_x_genes.barcodes.txt").read_text().splitlines())
    if matrix.shape != (n_barcodes, len(genes)):
        raise CellCallingError(
            f"kb matrix in {src} is {matrix.shape}, but it lists {n_barcodes} barcodes "
            f"and {len(genes)} genes."
        )
    is_viral = np.array([g in viral_genes for g in genes], dtype=bool)
    if not is_viral.any():
        raise CellCallingError(
            f"none of the {len(viral_genes)} viral gene IDs of this run is a column of {src}; "
            "the partition belongs to another reference. Re-run the analysis step, or pass "
            "--called-cells-file / use two-step --host-filter starsolo."
        )
    if is_viral.all():
        raise CellCallingError(
            f"every column of {src} is viral (a virus-only index): there is no host "
            "transcriptome to call cells on. Use two-step mode (--host-filter starsolo, which "
            "calls cells on the STARsolo host matrix) or pass --called-cells-file."
        )
    host = np.flatnonzero(~is_viral)
    mmwrite(dest / "cells_x_genes.mtx", matrix[:, host])
    shutil.copyfile(src / "cells_x_genes.barcodes.txt", dest / "cells_x_genes.barcodes.txt")
    (dest / "cells_x_genes.genes.txt").write_text("".join(f"{genes[i]}\n" for i in host))
    return dest, len(host), int(is_viral.sum())


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
    barcodes, totals = read_host_matrix_totals(Path(solo_dir))  # fail on a corrupt matrix first
    _clear_input_receipt(config)
    rscript = getattr(config, "cell_caller_rscript", DEFAULTS["cell_caller_rscript"])
    params = _emptydrops_params(config)
    out_tsv = _run_emptydrops(solo_dir, rscript=rscript, **params)
    cells = host_cells_from_tsv(out_tsv, min_comparable_umi)
    cells["host_totals"] = _nonzero_totals(barcodes, totals)
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
            "parameters": {**params, "effective_seed": effective_seed(params["seed"])},
            **_tool_provenance(rscript),
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
    _clear_input_receipt(config)
    raw = Path(config.output) / "host_filtered" / "star_tmp" / "Solo.out" / "GeneFull" / "raw"
    comparable = None
    host_totals: dict[str, float] | None = None
    complete = all(
        (raw / name).is_file() for name in ("matrix.mtx", "barcodes.tsv", "features.tsv")
    )
    if complete:
        barcodes, totals = read_host_matrix_totals(raw)
        comparable = {
            b for b, total in zip(barcodes, totals) if b in called and total >= min_comparable_umi
        }
        host_totals = _nonzero_totals(barcodes, totals)
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
    cells: HostCellSets = {"called": called, "comparable": comparable}
    if host_totals is not None:
        cells["host_totals"] = host_totals
    return cells


def _digests(directory: Path, names) -> dict[str, str]:
    """SHA-256 of each named file, so the receipt pins the exact cell-calling inputs."""
    from viralscan.run_safety import sha256_file

    return {name: sha256_file(directory / name) for name in names}


def _emptydrops_params(config) -> dict:
    return {
        "fdr": float(getattr(config, "emptydrops_fdr", DEFAULTS["emptydrops_fdr"])),
        "lower": float(getattr(config, "emptydrops_lower", DEFAULTS["emptydrops_lower"])),
        "niters": int(getattr(config, "emptydrops_niters", DEFAULTS["emptydrops_niters"])),
        "seed": int(getattr(config, "emptydrops_seed", DEFAULTS["emptydrops_seed"])),
    }


def _tool_provenance(rscript) -> dict:
    """R, DropletUtils and ViralScan code identity for the receipt; ``None`` where unavailable.

    The R side honours ``VIRALSCAN_R_LIBS`` as ``emptydrops.R`` does. Never raises:
    provenance must not stop a cell call.
    """
    code = (
        'l <- Sys.getenv("VIRALSCAN_R_LIBS"); if (nzchar(l) && dir.exists(l)) .libPaths(c(l, .libPaths()));'
        ' cat(R.version.string, as.character(packageVersion("DropletUtils")), sep = "\n")'
    )
    r_version = dropletutils_version = None
    try:
        out = subprocess.run(
            [rscript, "-e", code], check=True, capture_output=True, text=True, timeout=120
        ).stdout.splitlines()
        if len(out) == 2:
            r_version, dropletutils_version = out
    except Exception as exc:  # provenance never stops a cell call
        log.warning("cell calling: could not query the R and DropletUtils versions (%s)", exc)
    try:
        from viralscan.run_safety import software_identity

        identity = software_identity()
    except Exception as exc:
        log.warning("cell calling: could not read the ViralScan code identity (%s)", exc)
        identity = {}
    return {
        "r_version": r_version,
        "dropletutils_version": dropletutils_version,
        "code_sha": identity.get("git_commit"),  # None outside a git checkout
        "build_sha256": identity.get("build_sha256"),
    }


def _clear_input_receipt(config) -> None:
    """Remove a previous run's receipt, so a run that writes none cannot be mistaken for it."""
    output = getattr(config, "output", None)
    if output:
        (Path(output) / "results" / "cell_calling_input.json").unlink(missing_ok=True)


def _write_input_receipt(config, payload) -> None:
    results = Path(config.output) / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "cell_calling_input.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )


#: Prefix of the scratch host-only matrix directory inside ``counts_unfiltered``.
_SCRATCH_PREFIX = ".host_only_"


def _host_only_emptydrops(obs, matrix_dir, config) -> np.ndarray:
    """emptyDrops on the host gene columns of the combined kb matrix (PLAN CELLS-02).

    Cells are called on the host transcriptome only, never on viral reads. The
    called set is written to ``<matrix_dir>/emptydrops_cells.tsv`` as before
    (its ``total`` column is now host UMIs), so every reader keeps working.
    """
    output = getattr(config, "output", None)
    if not output:
        raise CellCallingError("host-only cell calling requires config.output (the run directory)")
    matrix_dir = Path(matrix_dir)
    _clear_input_receipt(config)
    matrix_genes = (matrix_dir / "cells_x_genes.genes.txt").read_text().splitlines()
    viral, partition = viral_partition(output, matrix_genes)
    params = _emptydrops_params(config)
    rscript = getattr(config, "cell_caller_rscript", DEFAULTS["cell_caller_rscript"])
    out_tsv = matrix_dir / "emptydrops_cells.tsv"
    # A SIGKILL skips the TemporaryDirectory cleanup; sweep what an earlier run left behind.
    for stale in matrix_dir.glob(f"{_SCRATCH_PREFIX}*"):
        shutil.rmtree(stale, ignore_errors=True)
    with tempfile.TemporaryDirectory(dir=matrix_dir, prefix=_SCRATCH_PREFIX) as tmp:
        _, n_host, n_viral = write_host_only_matrix(matrix_dir, viral, tmp)
        log.info(
            "cell_calling=emptydrops on %d host genes (%d viral columns excluded)", n_host, n_viral
        )
        mask = emptydrops_cells(obs, tmp, rscript=rscript, out_tsv=out_tsv, **params)
    names = ("cells_x_genes.mtx", "cells_x_genes.barcodes.txt", "cells_x_genes.genes.txt")
    _write_input_receipt(
        config,
        {
            "method": "emptydrops",
            "count_layer": "kb.host_genes",
            "matrix_dir": str(matrix_dir),
            "host_genes": n_host,
            "viral_genes_excluded": n_viral,
            "partition_source": partition.relative_to(output).as_posix(),
            "partition_sha256": _digests(partition.parent, (partition.name,))[partition.name],
            "input_sha256": _digests(matrix_dir, names),
            "emptydrops_cells_sha256": _digests(matrix_dir, (out_tsv.name,)),
            "parameters": {**params, "effective_seed": effective_seed(params["seed"])},
            **_tool_provenance(rscript),
        },
    )
    return mask


def call_cells(adata, config, matrix_dir=None) -> np.ndarray:
    """Return a boolean cell mask over ``adata.obs_names`` per ``config.cell_calling``.

    Parameters
    ----------
    matrix_dir : path | None
        kb ``counts_unfiltered`` directory — required for ``emptydrops`` (its ``.mtx``
        is what DropletUtils reads). Ignored by the other methods. emptyDrops sees the
        HOST gene columns only (viral columns are dropped via ``viral_partition``).

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

    if method != "emptydrops":
        _clear_input_receipt(config)  # only emptyDrops (and two-step) write one

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
        return _host_only_emptydrops(obs, mdir, config)

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

#!/usr/bin/env python3
"""Targeted fast profiling of ViralScan's 4 multimapping methods.

Rationale for this script vs profile_multimap.py
-------------------------------------------------
profile_multimap.py timed all 4 methods on the full 103M-row BUS file.
The 'equal' method took 19,964s (5.55 hours) on the OLD code (itertuples).
Two performance refactors landed before this run:
  - b7e9635: hoist EC invariants + zip over numpy arrays (~15% faster)
  - 2c2e6f0: vectorised em_gene_abundances (sparse matvec)
This script profiles the CURRENT code on subsamples and extrapolates.

Strategy
--------
P0: Load all inputs once. Assert shapes. Log key counts.
P3: cProfile on 1M-row subsample — equal and em methods. Save text immediately.
P2: EM convergence on 5M-row subsample — get iteration count, L1-delta trajectory.
    Extrapolate wall-time from subsample wall-time × (103M / 5M).
P5: Register headline values. Write timing_table.tsv and em_convergence.json.

Subsampling rationale: the zip loop is O(n_bus_records); the 1M/5M slice is
deterministic (first-N rows) and the relative cost attribution is invariant to
subsample size. Wall-time extrapolation is valid because the loop is linear in
row count. EM convergence iteration count is NOT affected by subsample size
(it depends on the ec_counts distribution, which is subsampled here — we note
this and recommend confirming on full data).

ANALYSIS_OK waivers
-------------------
[subsample-extrapolation] Wall-time extrapolated as t_sub × (n_full / n_sub).
    Loop is linear in n_bus_records; relative cost attribution (hot functions)
    is invariant to subsample size. Anchor: 19,964s for OLD itertuples code on
    full data. New code is ~15% faster per b7e9635 commit message.
[runtime-assert]    All asserts are developer tripwires in a non-production
    analysis script; this script is never run with -O.
[optional-input]    gene_names file fallback (same as profile_multimap.py).
[best-effort-fan-out] cProfile text parse failures logged; do not affect numbers.
[layer-choice]      adata.X = raw counts from kb-python (no .raw, no transform).
[file-selection]    All paths pinned to absolute paths.
[sample-filter]     bus_df.dropna: NaN rows from malformed BUS lines; < 1% expected.
[threshold]         EM_MAX_ITER=100, EM_TOL=1e-6, PSEUDOCOUNT=1.0 match viralscan
    defaults (confirmed: multimapping.py). CPROFILE_N_ROWS=1_000_000 and
    WALLTIME_N_ROWS=5_000_000 are analysis-only slices; not viralscan defaults.
"""

from __future__ import annotations

import cProfile
import io
import json
import os
import pstats
import re
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

# ---------------------------------------------------------------------------
# Paths — all absolute, pinned to specific run  ANALYSIS_OK[file-selection]
# ---------------------------------------------------------------------------
KB = Path(
    "/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_runs"
    "/reference_strategy_2026-06-28_fresh12b/runs/ebv__viralscan__combined"
    "/SRR12682296/kb-python"
)
BUS_TXT          = KB / "output.bus.txt"
EC_FILE          = KB / "matrix.ec"
TRANSCRIPTS_FILE = KB / "transcripts.txt"
BARCODES_FILE    = KB / "counts_unfiltered" / "cells_x_genes.barcodes.txt"
GENES_FILE       = KB / "counts_unfiltered" / "cells_x_genes.genes.txt"
# ANALYSIS_OK[file-selection]: pinned to absolute path; t2g line count = 230901
# matches transcripts.txt (confirmed by wc -l).
T2G_FILE   = Path(
    "/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_showcase"
    "/viralscan_showcase/fullrun/refs/merged/t2g_plus_anellovirus.txt"
)
ADATA_FILE = KB / "counts_unfiltered" / "adata.h5ad"

REPO_ROOT            = Path("/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan")
REGISTER_VALUE_PATH  = Path(
    "/home/mdmanurung/.claude/plugins/marketplaces/mycelium/skills/core/scripts"
)

sys.path.insert(0, str(REGISTER_VALUE_PATH))
sys.path.insert(0, str(REPO_ROOT / "src"))

from register_value import register_value  # noqa: E402
from viralscan.scripts.multimap import (  # noqa: E402
    load_barcodes,
    read_ec,
    strip_10x_suffix,
)
from viralscan.multimapping import build_multimap_layers  # noqa: E402

# ---------------------------------------------------------------------------
# Constants  ANALYSIS_OK[threshold]
# ---------------------------------------------------------------------------
EM_MAX_ITER       = 100
EM_TOL            = 1e-6
PSEUDOCOUNT       = 1.0
CPROFILE_N_ROWS   = 1_000_000     # ANALYSIS_OK[duplicate-config-source]: intentionally different from profile_multimap.py (10M); 1M is faster for the current zip-loop code
WALLTIME_N_ROWS   = 5_000_000     # 5M rows for wall-time subsample
FULL_DATA_N_ROWS  = 103_145_071   # full BUS file row count (Phase 0 confirmed)
METHODS           = ["equal", "host-conservative", "unique-weighted", "em"]
OLD_CODE_EQUAL_S  = 19_964.2      # anchor: 'equal' on OLD itertuples code, full data


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[fast_profile] {msg}", flush=True)


def check_nonzero(name: str, value: int) -> None:
    """Raise ValueError if value <= 0.  ANALYSIS_OK[runtime-assert]"""
    if value <= 0:
        raise ValueError(f"Expected {name} > 0, got {value}")


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")
    log(f"  written: {path}")


# ---------------------------------------------------------------------------
# Phase 0: Load inputs
# ---------------------------------------------------------------------------
def load_all_inputs() -> dict:
    """Load every input once; check shapes; log row counts."""
    log("Phase 0: Loading inputs...")
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    barcode_to_idx, n_cells = load_barcodes(str(BARCODES_FILE))
    check_nonzero("n_cells", n_cells)
    log(f"  barcodes: n_cells={n_cells:,}")

    import anndata as ad  # noqa: PLC0415 — import here to avoid top-level anndata dep check
    adata_orig = ad.read_h5ad(str(ADATA_FILE))
    gene_ids = list(adata_orig.var_names)
    n_genes = len(gene_ids)
    check_nonzero("n_genes", n_genes)
    log(f"  genes: n_genes={n_genes:,}")

    with open(GENES_FILE) as f:
        gene_ids_file = [line.strip() for line in f]
    if len(gene_ids_file) != n_genes:
        raise ValueError(
            f"Gene count mismatch: adata has {n_genes}, genes.txt has {len(gene_ids_file)}"
        )
    log(f"  gene list cross-check: OK (n={n_genes:,})")

    # ANALYSIS_OK[optional-input]: gene_names file may be absent
    gene_names_file = KB / "counts_unfiltered" / "cells_x_genes.genes.names.txt"
    if gene_names_file.exists():
        with open(gene_names_file) as f:
            gene_names = [line.strip() for line in f]
        if len(gene_names) != n_genes:
            raise ValueError(f"Gene names count mismatch: {len(gene_names)} vs {n_genes}")
    else:
        gene_names = gene_ids.copy()
        log("  WARNING: no gene_names file found, using gene_ids as gene_names")

    with open(TRANSCRIPTS_FILE) as f:
        transcripts = [line.strip() for line in f]
    check_nonzero("n_transcripts", len(transcripts))
    log(f"  transcripts: n={len(transcripts):,}")

    t2g = pd.read_csv(T2G_FILE, sep=r"\s+", header=None, usecols=[0, 1],
                      names=["transcript", "gene"])
    t2g_map = dict(zip(t2g["transcript"], t2g["gene"]))
    check_nonzero("n_t2g_entries", len(t2g_map))
    log(f"  t2g: n_mappings={len(t2g_map):,}")

    log("  reading EC file...")
    ec_map = read_ec(str(EC_FILE), transcripts, t2g_map, gene_ids)
    check_nonzero("n_ec_map", len(ec_map))
    n_multi_ec = sum(1 for genes in ec_map.values() if len(genes) > 1)
    log(f"  ec_map: n_ECs={len(ec_map):,}, n_multi-gene ECs={n_multi_ec:,}")

    log("  reading output.bus.txt (may take a moment)...")
    bus_df = pd.read_csv(
        str(BUS_TXT),
        sep="\t",
        header=None,
        names=["barcode", "umi", "ec", "count"],
        usecols=["barcode", "ec", "count"],
        dtype={"barcode": "category", "ec": "int32", "count": "int32"},
    )
    n_bus_before = len(bus_df)
    # ANALYSIS_OK[sample-filter]: NaN rows from partial/corrupt BUS text lines.
    bus_df = bus_df.dropna()
    n_bus_after = len(bus_df)
    n_dropped = n_bus_before - n_bus_after
    log(f"  bus records: loaded={n_bus_before:,}, after dropna={n_bus_after:,}, "
        f"dropped={n_dropped:,}")
    if n_dropped > n_bus_before * 0.01:
        raise ValueError(
            f"Too many dropped bus records: {n_dropped:,} / {n_bus_before:,}"
        )
    check_nonzero("n_bus_records", n_bus_after)

    # Strip trailing '-1' lane suffix from barcodes
    barcode_col = bus_df["barcode"]
    if isinstance(barcode_col.dtype, pd.CategoricalDtype):
        stripped = barcode_col.cat.categories.map(strip_10x_suffix)
        if stripped.is_unique:
            bus_df = bus_df.copy()
            bus_df["barcode"] = barcode_col.cat.rename_categories(stripped)
        else:
            bus_df = bus_df.copy()
            bus_df["barcode"] = barcode_col.astype("string").map(strip_10x_suffix)
    else:
        bus_df = bus_df.copy()
        bus_df["barcode"] = barcode_col.map(strip_10x_suffix)

    # Viral gene indices: non-ENSG gene IDs (heuristic, no analysis.txt available)
    viral_gene_indices: set[int] = {
        i for i, gid in enumerate(gene_ids) if not str(gid).startswith("ENSG")
    }
    n_viral = len(viral_gene_indices)
    log(f"  viral_gene_indices: n_viral_genes={n_viral:,} (non-ENSG heuristic)")
    check_nonzero("n_viral_genes", n_viral)

    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_load_kb = rss_after - rss_before
    log(f"  RSS after loading: watermark_delta={rss_load_kb:,} KB "
        f"({rss_load_kb / 1024:.0f} MB), abs_watermark={rss_after:,} KB")

    return dict(
        barcode_to_idx=barcode_to_idx,
        n_cells=n_cells,
        gene_ids=gene_ids,
        gene_names=gene_names,
        n_genes=n_genes,
        ec_map=ec_map,
        n_multi_ec=n_multi_ec,
        viral_gene_indices=viral_gene_indices,
        n_viral=n_viral,
        bus_df=bus_df,
        n_bus=n_bus_after,
        adata_orig=adata_orig,
        original_counts=adata_orig.X,  # ANALYSIS_OK[layer-choice]: raw counts
        rss_load_kb=rss_load_kb,
    )


# ---------------------------------------------------------------------------
# Phase 3 (fast): cProfile on 1M-row subsample
# ---------------------------------------------------------------------------
def run_cprofile_subsample(inputs: dict, out_dir: Path) -> tuple[str, str]:
    """Run cProfile on a fixed 1M-row slice; write text immediately."""
    log(f"\nPhase 3 (cProfile on {CPROFILE_N_ROWS:,}-row subsample)...")

    bus_full = inputs["bus_df"]
    # ANALYSIS_OK[sample-filter]: subsample is for cProfile only — timing and EM
    # use WALLTIME_N_ROWS or full data. Slice is deterministic (first-N rows).
    bus_sub = bus_full.iloc[:CPROFILE_N_ROWS].copy()
    log(f"  subsampled to first {CPROFILE_N_ROWS:,} rows (of {len(bus_full):,})")

    def _profile_one(method_name: str) -> str:
        profiler = cProfile.Profile()
        profiler.enable()
        _ = build_multimap_layers(
            bus_df=bus_sub,
            barcode_to_idx=inputs["barcode_to_idx"],
            ec_map=inputs["ec_map"],
            n_cells=inputs["n_cells"],
            n_genes=inputs["n_genes"],
            viral_gene_indices=inputs["viral_gene_indices"],
            original_counts=inputs["original_counts"],
            method=method_name,
            pseudocount=PSEUDOCOUNT,
            em_max_iter=EM_MAX_ITER,
            em_tol=EM_TOL,
        )
        profiler.disable()
        sio = io.StringIO()
        pstats.Stats(profiler, stream=sio).sort_stats("cumulative").print_stats(40)
        return sio.getvalue()

    t0 = time.perf_counter()
    profile_equal_text = _profile_one("equal")
    t1 = time.perf_counter()
    log(f"  cProfile equal: {t1-t0:.1f}s (with profiling overhead)")

    # Write immediately so data survives if EM cProfile hangs
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "cprofile_equal.txt").write_text(profile_equal_text)
    log(f"  written: {out_dir}/cprofile_equal.txt")

    log(f"  cProfile equal (top 12 lines):")
    for line in profile_equal_text.split("\n")[3:15]:
        log(f"    {line}")

    t2 = time.perf_counter()
    profile_em_text = _profile_one("em")
    t3 = time.perf_counter()
    log(f"  cProfile em: {t3-t2:.1f}s (with profiling overhead)")

    (out_dir / "cprofile_em.txt").write_text(profile_em_text)
    log(f"  written: {out_dir}/cprofile_em.txt")

    log(f"  cProfile em (top 12 lines):")
    for line in profile_em_text.split("\n")[3:15]:
        log(f"    {line}")

    return profile_equal_text, profile_em_text


# ---------------------------------------------------------------------------
# Phase 2 (fast): Wall-time on 5M subsample + EM convergence
# ---------------------------------------------------------------------------
def run_walltime_subsample(inputs: dict, out_dir: Path) -> dict:
    """Time all 4 methods on 5M-row subsample; extrapolate to full-data."""
    log(f"\nPhase 2 (wall-time on {WALLTIME_N_ROWS:,}-row subsample)...")
    log("  Note: extrapolation assumes loop is linear in n_bus_records.")
    log(f"  Anchor for validation: equal on OLD itertuples code = {OLD_CODE_EQUAL_S:.1f}s (full data).")

    bus_full = inputs["bus_df"]
    # ANALYSIS_OK[sample-filter]: subsample for wall-time measurement only.
    bus_sub = bus_full.iloc[:WALLTIME_N_ROWS].copy()
    log(f"  subsampled to first {WALLTIME_N_ROWS:,} rows (of {len(bus_full):,})")

    n_full  = FULL_DATA_N_ROWS
    n_sub   = WALLTIME_N_ROWS
    scale   = n_full / n_sub

    wall_sub: dict[str, float] = {}
    wall_extrap: dict[str, float] = {}
    em_convergence: dict = {}

    for method in METHODS:
        log(f"  timing method={method!r} on {n_sub:,} rows...")
        t0 = time.perf_counter()
        layers = build_multimap_layers(
            bus_df=bus_sub,
            barcode_to_idx=inputs["barcode_to_idx"],
            ec_map=inputs["ec_map"],
            n_cells=inputs["n_cells"],
            n_genes=inputs["n_genes"],
            viral_gene_indices=inputs["viral_gene_indices"],
            original_counts=inputs["original_counts"],
            method=method,
            pseudocount=PSEUDOCOUNT,
            em_max_iter=EM_MAX_ITER,
            em_tol=EM_TOL,
        )
        t1 = time.perf_counter()
        elapsed = t1 - t0

        # ANALYSIS_OK[runtime-assert]: validating own output shape
        expected = (inputs["n_cells"], inputs["n_genes"])
        for lname, lmat in [
            ("equal", layers.equal),
            ("host_conservative", layers.host_conservative),
            ("unique_weighted", layers.unique_weighted),
            ("corrected", layers.corrected),
        ]:
            if lmat.shape != expected:
                raise ValueError(
                    f"Layer {lname} shape mismatch: {lmat.shape} != {expected}"
                )

        wall_sub[method] = elapsed
        wall_extrap[method] = elapsed * scale
        log(f"    subsample={elapsed:.2f}s  →  extrap_full={wall_extrap[method]:.0f}s "
            f"  ({wall_extrap[method]/3600:.2f}h)")

    # Write intermediate result immediately
    sub_result = {
        "n_sub_rows": n_sub,
        "n_full_rows": n_full,
        "scale_factor": scale,
        "wall_sub_s": {k: round(v, 3) for k, v in wall_sub.items()},
        "wall_extrap_s": {k: round(v, 0) for k, v in wall_extrap.items()},
        "anchor_old_code_equal_s": OLD_CODE_EQUAL_S,
        "note": (
            "equal and host-conservative and unique-weighted share ONE bus pass; "
            "EM adds em_records allocation loop after the shared pass. "
            "em_gene_abundances is now vectorised (sparse matvec, commit 2c2e6f0). "
            "Extrapolation valid: loop is O(n_bus_records)."
        ),
    }
    write_json(out_dir / "walltime_subsample.json", sub_result)

    return dict(
        wall_sub=wall_sub,
        wall_extrap=wall_extrap,
        em_convergence=em_convergence,
    )


# ---------------------------------------------------------------------------
# Phase 2b: EM convergence trajectory
# ---------------------------------------------------------------------------
def run_em_convergence(inputs: dict, out_dir: Path) -> dict:
    """Run instrumented EM on a subsample to get convergence iteration count."""
    log(f"\nPhase 2b (EM convergence on {WALLTIME_N_ROWS:,}-row subsample)...")
    log("  NOTE: em_gene_abundances is now vectorised (sparse matvec, commit 2c2e6f0).")
    log("  The instrumented loop below replicates the NEW vectorised version.")

    bus_full = inputs["bus_df"]
    bus_sub  = bus_full.iloc[:WALLTIME_N_ROWS].copy()

    barcode_to_idx = inputs["barcode_to_idx"]
    ec_map         = inputs["ec_map"]

    # Collect em_ec_counts from subsample
    log("  collecting em_ec_counts from bus subsample (one pass)...")
    t_pass_start = time.perf_counter()
    em_ec_counts: dict[tuple[int, ...], float] = {}
    n_multi_bus = 0

    bc_arr  = bus_sub["barcode"].to_numpy()
    ec_arr  = bus_sub["ec"].to_numpy()
    cnt_arr = bus_sub["count"].to_numpy()

    for bc, ec_raw, count_raw in zip(bc_arr, ec_arr, cnt_arr):
        if pd.isna(ec_raw):
            continue
        cell_idx = barcode_to_idx.get(bc)
        if cell_idx is None:
            continue
        info = ec_map.get(int(ec_raw))
        if not info:
            continue
        genes_in_ec = list(info)
        if len(genes_in_ec) <= 1:
            continue
        distinct_key = tuple(dict.fromkeys(genes_in_ec))
        em_ec_counts[distinct_key] = em_ec_counts.get(distinct_key, 0.0) + float(count_raw)
        n_multi_bus += 1

    t_pass_end = time.perf_counter()
    log(f"  bus pass: {t_pass_end - t_pass_start:.1f}s, "
        f"n_multi_bus_rows={n_multi_bus:,}, n_distinct_multi_ec_keys={len(em_ec_counts):,}")

    original_counts = inputs["original_counts"]
    unique_per_gene = np.asarray(
        original_counts.sum(axis=0)
        if sparse.issparse(original_counts)
        else original_counts.sum(axis=0)
    ).reshape(-1)

    if len(unique_per_gene) != inputs["n_genes"]:
        raise ValueError(  # ANALYSIS_OK[runtime-assert]
            f"unique_per_gene shape mismatch: {len(unique_per_gene)} vs {inputs['n_genes']}"
        )
    log(f"  unique_per_gene: shape={unique_per_gene.shape}, sum={unique_per_gene.sum():.0f}")

    # Instrumented EM loop — replicates vectorised em_gene_abundances
    # (src/viralscan/multimapping.py, commit 2c2e6f0 / b7e9635)
    n_ec = len(em_ec_counts)
    log(f"  building sparse incidence matrix ({n_ec:,} ECs × {inputs['n_genes']:,} genes)...")
    t_build_start = time.perf_counter()
    ec_rows_list: list[np.ndarray] = []
    ec_cols_list: list[np.ndarray] = []
    counts_vec  = np.empty(n_ec, dtype=float)
    genes_per_ec = np.empty(n_ec, dtype=float)
    for e, (genes, count) in enumerate(em_ec_counts.items()):
        g = np.asarray(genes, dtype=int)
        ec_rows_list.append(np.full(g.shape[0], e, dtype=int))
        ec_cols_list.append(g)
        counts_vec[e]   = count
        genes_per_ec[e] = g.shape[0]
    row_idx    = np.concatenate(ec_rows_list)
    col_idx    = np.concatenate(ec_cols_list)
    from scipy.sparse import csr_matrix  # noqa: PLC0415
    incidence  = csr_matrix(
        (np.ones(row_idx.shape[0], dtype=float), (row_idx, col_idx)),
        shape=(n_ec, inputs["n_genes"]),
    )
    incidence_t = incidence.T.tocsr()
    t_build_end = time.perf_counter()
    log(f"  incidence matrix built: {t_build_end - t_build_start:.2f}s, "
        f"nnz={incidence.nnz:,}")

    theta = unique_per_gene.astype(float) + float(PSEUDOCOUNT)
    iter_deltas: list[float] = []
    converged_iter: int | None = None

    log(f"  starting EM (max_iter={EM_MAX_ITER}, tol={EM_TOL})...")
    t_em_start = time.perf_counter()
    for i_iter in range(int(EM_MAX_ITER)):
        s    = incidence @ theta
        good = s > 1e-12
        weighted = theta * (incidence_t @ np.where(good, counts_vec / np.where(good, s, 1.0), 0.0))
        equal_mass = incidence_t @ np.where(good, 0.0, counts_vec / genes_per_ec)
        new   = unique_per_gene + weighted + equal_mass
        denom = float(theta.sum()) or 1.0
        l1_delta = float(np.abs(new - theta).sum()) / denom
        iter_deltas.append(l1_delta)
        theta = new
        if l1_delta < EM_TOL:
            converged_iter = i_iter + 1
            log(f"  EM converged at iteration {converged_iter}, L1_delta={l1_delta:.2e}")
            break
    else:
        converged_iter = EM_MAX_ITER
        log(f"  EM did NOT converge in {EM_MAX_ITER} iters; final L1_delta={iter_deltas[-1]:.2e}")
    t_em_end = time.perf_counter()
    em_loop_s = t_em_end - t_em_start
    log(f"  EM loop wall-time (vectorised, {converged_iter} iters): {em_loop_s:.3f}s")
    log(f"  First 5 L1_deltas: {[f'{d:.3e}' for d in iter_deltas[:5]]}")
    log(f"  Last 5 L1_deltas:  {[f'{d:.3e}' for d in iter_deltas[-5:]]}")

    conv_data = {
        "subsample_n_rows": WALLTIME_N_ROWS,
        "converged_iter": converged_iter,
        "n_multi_bus_rows": n_multi_bus,
        "n_em_ec_keys": n_ec,
        "incidence_nnz": int(incidence.nnz),
        "em_loop_time_s": round(em_loop_s, 4),
        "bus_pass_time_s": round(t_pass_end - t_pass_start, 2),
        "incidence_build_time_s": round(t_build_end - t_build_start, 4),
        "l1_delta_trajectory": [round(d, 8) for d in iter_deltas],
        "note": (
            "em_gene_abundances is vectorised (sparse matvec per iteration). "
            "Convergence iteration count from subsample may differ slightly from "
            "full-data run (different ec_counts distribution). "
            "EM loop time is for the vectorised version only."
        ),
    }
    write_json(out_dir / "em_convergence.json", conv_data)

    return conv_data


# ---------------------------------------------------------------------------
# Phase 5: Summary table and register_value
# ---------------------------------------------------------------------------
def build_summary(
    inputs: dict,
    wall_result: dict,
    em_conv: dict,
    profile_equal_text: str,
    profile_em_text: str,
    out_dir: Path,
) -> None:
    log("\nPhase 5: Summary and value registration...")

    wall_sub    = wall_result["wall_sub"]
    wall_extrap = wall_result["wall_extrap"]

    rows = []
    for method in METHODS:
        rows.append({
            "method": method,
            "wall_sub_s": round(wall_sub[method], 3),
            "wall_extrap_s": round(wall_extrap[method], 0),
            "wall_extrap_h": round(wall_extrap[method] / 3600, 2),
            "anchor_notes": (
                "NEW zip-loop code; extrap from 5M-row slice; "
                f"OLD itertuples anchor={OLD_CODE_EQUAL_S}s"
                if method == "equal" else ""
            ),
        })

    df = pd.DataFrame(rows)
    # ANALYSIS_OK[runtime-assert]: validating own output
    if df.shape[0] != len(METHODS):
        raise ValueError(f"Summary row count mismatch: {df.shape[0]} != {len(METHODS)}")

    log("\n=== TIMING TABLE (extrapolated to full 103M-row BUS) ===")
    log(df.to_string(index=False))

    timing_path = out_dir / "timing_table.tsv"
    df.to_csv(timing_path, sep="\t", index=False)
    log(f"\n  timing table written to: {timing_path}")

    # Register headline values
    def reg(key: str, val: object, prov: str) -> None:
        register_value(key, val, provenance=prov)

    reg("n_bus_records",    inputs["n_bus"],       "KB/output.bus.txt:wc-l")
    reg("n_cells",          inputs["n_cells"],     "KB/counts_unfiltered/cells_x_genes.barcodes.txt:wc-l")
    reg("n_genes",          inputs["n_genes"],     "KB/counts_unfiltered/adata.h5ad:adata.n_vars")
    reg("n_ec_map",         len(inputs["ec_map"]), "KB/matrix.ec:parsed entries")
    reg("n_multi_gene_ec",  inputs["n_multi_ec"],  "KB/matrix.ec:entries with >1 gene")
    reg("n_viral_genes",    inputs["n_viral"],     "gene_ids:non-ENSG count")
    reg("rss_load_mb",      round(inputs["rss_load_kb"] / 1024), "resource.getrusage:after Phase0 load")

    for method in METHODS:
        key_sub    = f"wall_sub_s_{method.replace('-', '_')}"
        key_extrap = f"wall_extrap_s_{method.replace('-', '_')}"
        reg(key_sub,    round(wall_sub[method], 3),    f"time.perf_counter:{WALLTIME_N_ROWS}row-subsample")
        reg(key_extrap, round(wall_extrap[method], 0), f"extrapolated×{FULL_DATA_N_ROWS/WALLTIME_N_ROWS:.1f}")

    reg("em_converged_iter",     em_conv["converged_iter"],
        f"instrumented EM vectorised:{WALLTIME_N_ROWS}row-subsample")
    reg("em_loop_time_s",        em_conv["em_loop_time_s"],
        "instrumented vectorised EM loop wall-time (subsample)")
    reg("em_incidence_build_s",  em_conv["incidence_build_time_s"],
        "sparse incidence matrix construction time")
    reg("em_n_ec_keys",          em_conv["n_em_ec_keys"],
        f"distinct multi-gene EC keys in {WALLTIME_N_ROWS}row-subsample")

    reg("anchor_old_code_equal_s", OLD_CODE_EQUAL_S,
        "profile_multimap.py Phase1: equal on old itertuples code, full data")

    log("  values registered to analysis/multimap_profiling/outputs/numbers.json")

    # Print recommendation
    _print_recommendation(inputs, wall_extrap, em_conv, profile_equal_text, profile_em_text)


# ---------------------------------------------------------------------------
# Recommendation
# ---------------------------------------------------------------------------
def _parse_func(profile_text: str, func_pattern: str) -> tuple[float | None, float | None]:
    """Return (tottime, cumtime) for first matching function, or (None, None)."""
    for line in profile_text.splitlines():
        if func_pattern in line:
            m = re.match(r"^\s*(\d+/?\d*)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(.+)$", line)
            if m:
                try:
                    return float(m.group(2)), float(m.group(4))
                except ValueError:
                    pass  # ANALYSIS_OK[best-effort-fan-out]: non-matching cProfile header lines
    return None, None


def _extract_total_time(profile_text: str) -> float | None:
    """Extract total profiled time from cProfile header."""
    for line in profile_text.splitlines():
        m = re.search(r"([\d.]+) seconds", line)
        if m:
            return float(m.group(1))
    return None


def _print_recommendation(
    inputs: dict,
    wall_extrap: dict[str, float],
    em_conv: dict,
    profile_equal_text: str,
    profile_em_text: str,
) -> None:
    log("\n=== SPEEDUP RECOMMENDATION ===")

    t_equal = wall_extrap["equal"]
    t_em    = wall_extrap["em"]
    conv_iter = em_conv["converged_iter"]
    em_loop_s = em_conv["em_loop_time_s"]
    n_em_keys = em_conv["n_em_ec_keys"]

    # --- Hypothesis (a): _matrix_value dominates non-EM methods ---
    total_equal = _extract_total_time(profile_equal_text)
    bm_tot, bm_cum = _parse_func(profile_equal_text, "build_multimap_layers")
    matv_tot, matv_cum = _parse_func(profile_equal_text, "_matrix_value")

    log(f"\n(a) Hot-function attribution (equal, {CPROFILE_N_ROWS:,}-row cProfile):")
    if total_equal:
        log(f"  Total profiled time: {total_equal:.1f}s")
    if bm_cum and total_equal:
        log(f"  build_multimap_layers: tottime={bm_tot:.1f}s, cumtime={bm_cum:.1f}s "
            f"({bm_cum/total_equal:.1%} of total)")
    else:
        log("  build_multimap_layers: not found in top-40 profile lines")
    if matv_cum and total_equal:
        log(f"  _matrix_value: tottime={matv_tot:.1f}s, cumtime={matv_cum:.1f}s "
            f"({matv_cum/total_equal:.1%} of total)")
        if matv_cum > 0.5 * (total_equal or 1):
            log("  CONFIRMED: _matrix_value dominates non-EM build time (>50%).")
        elif matv_cum > 0.3 * (total_equal or 1):
            log("  LIKELY: _matrix_value is major contributor (>30%).")
        else:
            log(f"  _matrix_value = {matv_cum/total_equal:.1%}; check cprofile_equal.txt for full picture.")
    else:
        log("  _matrix_value: not found in top-40 profile lines")

    # --- Hypothesis (b): EM em_gene_abundances cost ---
    total_em = _extract_total_time(profile_em_text)
    em_bm_tot, em_bm_cum = _parse_func(profile_em_text, "build_multimap_layers")
    emga_tot, emga_cum   = _parse_func(profile_em_text, "em_gene_abundances")

    log(f"\n(b) EM method hot-function attribution (em, {CPROFILE_N_ROWS:,}-row cProfile):")
    if total_em:
        log(f"  Total profiled time (em): {total_em:.1f}s")
    if emga_cum and total_em:
        log(f"  em_gene_abundances: tottime={emga_tot:.1f}s, cumtime={emga_cum:.1f}s "
            f"({emga_cum/total_em:.1%} of EM total)")
    else:
        log("  em_gene_abundances: not found in top-40 profile lines")
    if em_bm_tot and bm_tot and total_equal:
        em_records_extra = (em_bm_tot or 0) - (bm_tot or 0)
        log(f"  em_records allocation loop (build tottime delta): "
            f"{em_records_extra:.1f}s ({em_records_extra/total_equal:.1%} of equal total)")

    # --- Convergence ---
    log(f"\n(iii) EM convergence:")
    log(f"  EM converged at iteration {conv_iter} / {EM_MAX_ITER} "
        f"(vectorised sparse-matvec, {WALLTIME_N_ROWS:,}-row subsample)")
    log(f"  EM loop time (vectorised): {em_loop_s:.3f}s for {conv_iter} iterations")
    log(f"  EM is now fast (sparse matvec); the bus-pass (zip loop) dominates all methods.")
    if conv_iter < EM_MAX_ITER:
        log(f"  Recommended: em_max_iter={conv_iter + 5} (add small buffer). "
            f"Saves {(1 - conv_iter/EM_MAX_ITER)*100:.0f}% of EM loop overhead.")

    # --- Overall priority ranking ---
    log(f"\n=== PRIORITY ORDER ===")
    log(f"  Context: current code uses zip(bc_arr, ec_arr, cnt_arr) over numpy arrays.")
    log(f"  ALL 3 non-EM methods share ONE pass; EM adds em_records allocation loop after.")
    log(f"  em_gene_abundances is NOW VECTORISED (sparse matvec, commit 2c2e6f0).")
    log(f"")
    log(f"  1. Vectorize the main bus-pass (ii) — affects ALL methods, highest total gain.")
    log(f"     Current (zip): extrap ~{t_equal:.0f}s for 'equal' method.")
    log(f"     Strategy: groupby (barcode, ec) → merge ec_info → vectorized per-EC ops.")
    log(f"     Expected: <10s (100–1000× speedup).")
    log(f"")
    if matv_cum and total_equal and matv_cum > 0.3 * total_equal:
        log(f"  1a. Within the bus-pass: eliminate _matrix_value hot-spot.")
        log(f"     _matrix_value = {matv_cum/total_equal:.1%} of build time in cProfile.")
        log(f"     Fix: precompute original_counts as dense row-slices or CSR row view.")
        log(f"     This alone gives a meaningful speedup without full vectorization.")
        log(f"")
    log(f"  2. Lower em_max_iter (iii) — trivial 1-line change.")
    if conv_iter < EM_MAX_ITER:
        log(f"     Set em_max_iter={conv_iter + 5}. EM loop is already fast (vectorised).")
    log(f"")
    log(f"  3. Vectorize main-pass fully (pandas groupby + merge) to close the gap entirely.")
    log(f"")
    log(f"  BUILD vs SWAP distinction (important for user recommendation):")
    log(f"  - Non-EM method selection = O(1) layer pointer in h5ad. Re-running multimap")
    log(f"    to switch equal→host-conservative = zero reprocessing in ViralScan.")
    log(f"  - EM requires a full build pass (~{t_em:.0f}s). Switching to/from EM = full reprocess.")
    log(f"  - Recommendation: default to 'equal' (fastest, transparent). Use 'em' only")
    log(f"    when per-gene EM abundance estimates are explicitly needed.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    log("=== ViralScan multimap FAST profiling ===")
    log(f"  Code version: current HEAD (zip loop + vectorised EM)")
    log(f"  Anchor: equal on OLD itertuples code = {OLD_CODE_EQUAL_S:.1f}s (full data)")
    log(f"  cProfile slice: {CPROFILE_N_ROWS:,} rows")
    log(f"  Wall-time slice: {WALLTIME_N_ROWS:,} rows (extrapolated × "
        f"{FULL_DATA_N_ROWS/WALLTIME_N_ROWS:.1f})")

    out_dir = Path(__file__).parent.parent / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Phase 0
    inputs = load_all_inputs()

    # Phase 3 — cProfile (1M rows); writes outputs immediately
    profile_equal_text, profile_em_text = run_cprofile_subsample(inputs, out_dir)

    # Phase 2 — wall-time (5M rows); extrapolate
    wall_result = run_walltime_subsample(inputs, out_dir)

    # Phase 2b — EM convergence
    em_conv = run_em_convergence(inputs, out_dir)

    # Phase 5 — summary + register_value
    build_summary(inputs, wall_result, em_conv, profile_equal_text, profile_em_text, out_dir)

    log("\n=== fast profiling complete ===")


if __name__ == "__main__":
    main()

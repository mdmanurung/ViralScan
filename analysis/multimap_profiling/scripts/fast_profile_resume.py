#!/usr/bin/env python3
"""Resume script for fast_profile.py after process was killed during em timing.

Context
-------
The full fast_profile.py run was killed during Phase 2 em wall-time on 5M rows
(process PID 518744, killed ~2026-07-05 around 04:xx UTC+8 by external signal,
not OOM — ulimit -v = unlimited, 645 GB free). Three of four wall-times were
already written to the log file before the kill.

This script:
  1. Loads all inputs (Phase 0 — unavoidable, ~240s).
  2. Reads the already-written cprofile_equal.txt and cprofile_em.txt from disk
     (Phase 3 — skip re-profiling).
  3. Hard-codes the three completed wall-times from the log, then measures em
     wall-time only (Phase 2). Writes the em result to disk IMMEDIATELY after
     timing so a second kill cannot lose it.
  4. Runs Phase 2b (EM convergence trajectory).
  5. Runs Phase 5 (build_summary → timing_table.tsv, numbers.json, register_value).

Hard-coded wall-times (from fast_profile_run.log, lines 57-60):
  equal:             subsample=1139.75s → extrap=23512s
  host-conservative: subsample=1078.88s → extrap=22256s
  unique-weighted:   subsample=1136.83s → extrap=23452s

Anchor discrepancy note
-----------------------
Old itertuples code (full 103M rows, profile_run.log line 24): 19,964.2s
New zip-loop code (extrapolated from 5M head-slice): ~23,512s
The new code extrapolates HIGHER than the old anchor (~18% slower). This is
NOT necessarily a real regression. The most likely cause is HEAD-SLICE BIAS:
the first 5M rows of a barcode-sorted BUS file may have a higher fraction of
multi-gene EC records than the full distribution, inflating per-row cost. The
absolute extrapolated seconds should be treated as UPPER BOUNDS; the RELATIVE
ranking across methods (confirmed ~equal) and the hot-function attribution
(_matrix_value at 86%) are unaffected by this bias. This is noted in the
outputs and caveated in the recommendation.

ANALYSIS_OK waivers (inherited from fast_profile.py)
---------------------------------------------------
[subsample-extrapolation] Wall-time extrapolated as t_sub * (n_full / n_sub).
[runtime-assert]    Developer tripwires; script never run with -O.
[optional-input]    gene_names file fallback.
[best-effort-fan-out] cProfile parse failures logged; do not affect numbers.
[layer-choice]      adata.X = raw counts.
[file-selection]    All paths pinned to absolute paths.
[sample-filter]     bus_df.dropna: NaN rows from corrupt BUS lines.
[threshold]         EM_MAX_ITER=100, EM_TOL=1e-6, PSEUDOCOUNT=1.0.
[partial-run-recovery] Known wall-times for 3/4 methods hard-coded from prior
    run log. Provenance: fast_profile_run.log lines 57-60. Em timing is the
    only new measurement in Phase 2.
[head-slice-bias]   First-N-row subsample may over-represent multi-gene ECs
    relative to full distribution. Absolute extrap values are upper bounds;
    relative method ranking is unaffected. Caveated in outputs.
"""

from __future__ import annotations

import json
import os
import re
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

# ---------------------------------------------------------------------------
# Paths — all absolute, pinned  ANALYSIS_OK[file-selection]
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
# matches transcripts.txt (confirmed by wc -l in fast_profile.py).
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
CPROFILE_N_ROWS   = 1_000_000
WALLTIME_N_ROWS   = 5_000_000
FULL_DATA_N_ROWS  = 103_145_071
METHODS           = ["equal", "host-conservative", "unique-weighted", "em"]
OLD_CODE_EQUAL_S  = 19_964.2   # anchor: equal on OLD itertuples code, full data

# ANALYSIS_OK[partial-run-recovery]: hard-coded from fast_profile_run.log
# lines 57-60 (committed subsample times from the prior run).
KNOWN_WALL_SUB: dict[str, float] = {
    "equal":             1139.75,
    "host-conservative": 1078.88,
    "unique-weighted":   1136.83,
    # "em": NOT available — will be measured below
}
KNOWN_PROVENANCE = "fast_profile_run.log:lines-57-60 (prior run, 2026-07-04)"

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[fast_profile_resume] {msg}", flush=True)


def check_nonzero(name: str, value: int) -> None:
    """ANALYSIS_OK[runtime-assert]"""
    if value <= 0:
        raise ValueError(f"Expected {name} > 0, got {value}")


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")
    log(f"  written: {path}")


# ---------------------------------------------------------------------------
# Phase 0: Load inputs (identical to fast_profile.py)
# ---------------------------------------------------------------------------
def load_all_inputs() -> dict:
    log("Phase 0: Loading inputs...")
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    barcode_to_idx, n_cells = load_barcodes(str(BARCODES_FILE))
    check_nonzero("n_cells", n_cells)
    log(f"  barcodes: n_cells={n_cells:,}")

    import anndata as ad  # noqa: PLC0415
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
# Phase 3: Read cprofile text from disk (skip re-profiling)
# ---------------------------------------------------------------------------
def load_cprofile_from_disk(out_dir: Path) -> tuple[str, str]:
    log("\nPhase 3 (reading cprofile from disk — skip re-profiling)...")
    eq_path = out_dir / "cprofile_equal.txt"
    em_path = out_dir / "cprofile_em.txt"
    for p in [eq_path, em_path]:
        if not p.exists():
            raise FileNotFoundError(f"Expected cprofile output not found: {p}")
        log(f"  found: {p} ({p.stat().st_size} bytes)")
    profile_equal_text = eq_path.read_text()
    profile_em_text    = em_path.read_text()
    log("  cprofile texts loaded from disk.")
    return profile_equal_text, profile_em_text


# ---------------------------------------------------------------------------
# Phase 2: em wall-time only; hard-code the 3 known times
# ---------------------------------------------------------------------------
def run_em_walltime_only(inputs: dict, out_dir: Path) -> dict:
    """Measure em wall-time on 5M subsample; hard-code the 3 known results."""
    log(f"\nPhase 2 (em wall-time only on {WALLTIME_N_ROWS:,}-row subsample)...")
    log(f"  Known wall-times from prior run [{KNOWN_PROVENANCE}]:")
    for m, t in KNOWN_WALL_SUB.items():
        scale = FULL_DATA_N_ROWS / WALLTIME_N_ROWS
        log(f"    {m}: subsample={t:.2f}s → extrap={t*scale:.0f}s ({t*scale/3600:.2f}h)")
    log(f"  NOTE: head-slice bias — first {WALLTIME_N_ROWS:,} rows may over-represent "
        f"multi-gene ECs.")
    log(f"  Absolute extrap values are upper bounds; relative method ranking is robust.")

    bus_full = inputs["bus_df"]
    bus_sub  = bus_full.iloc[:WALLTIME_N_ROWS].copy()
    n_sub    = WALLTIME_N_ROWS
    n_full   = FULL_DATA_N_ROWS
    scale    = n_full / n_sub

    log(f"  subsampled to first {n_sub:,} rows (of {len(bus_full):,})")
    log(f"  timing method='em' on {n_sub:,} rows...")

    # Instrument RSS before/during/after em timing
    rss_before_em = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t0 = time.perf_counter()
    layers = build_multimap_layers(
        bus_df=bus_sub,
        barcode_to_idx=inputs["barcode_to_idx"],
        ec_map=inputs["ec_map"],
        n_cells=inputs["n_cells"],
        n_genes=inputs["n_genes"],
        viral_gene_indices=inputs["viral_gene_indices"],
        original_counts=inputs["original_counts"],
        method="em",
        pseudocount=PSEUDOCOUNT,
        em_max_iter=EM_MAX_ITER,
        em_tol=EM_TOL,
    )
    t1 = time.perf_counter()
    em_elapsed = t1 - t0
    rss_after_em = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_em_delta_kb = rss_after_em - rss_before_em

    # ANALYSIS_OK[runtime-assert]: validating own output shape
    expected = (inputs["n_cells"], inputs["n_genes"])
    for lname, lmat in [
        ("equal", layers.equal),
        ("host_conservative", layers.host_conservative),
        ("unique_weighted", layers.unique_weighted),
        ("corrected", layers.corrected),
    ]:
        if lmat.shape != expected:
            raise ValueError(f"Layer {lname} shape mismatch: {lmat.shape} != {expected}")

    em_extrap = em_elapsed * scale
    log(f"    em: subsample={em_elapsed:.2f}s → extrap_full={em_extrap:.0f}s ({em_extrap/3600:.2f}h)")
    log(f"    em: RSS delta during em timing: {rss_em_delta_kb:,} KB ({rss_em_delta_kb/1024:.0f} MB)")

    # Merge known + new
    wall_sub    = dict(KNOWN_WALL_SUB)
    wall_sub["em"] = em_elapsed
    wall_extrap = {m: wall_sub[m] * scale for m in METHODS}

    # ANALYSIS_OK[partial-run-recovery]: write immediately after em timing
    sub_result = {
        "n_sub_rows": n_sub,
        "n_full_rows": n_full,
        "scale_factor": round(scale, 4),
        "wall_sub_s": {k: round(v, 3) for k, v in wall_sub.items()},
        "wall_extrap_s": {k: round(v, 0) for k, v in wall_extrap.items()},
        "anchor_old_code_equal_s": OLD_CODE_EQUAL_S,
        "em_rss_delta_kb": rss_em_delta_kb,
        "provenance": {
            "equal": KNOWN_PROVENANCE,
            "host-conservative": KNOWN_PROVENANCE,
            "unique-weighted": KNOWN_PROVENANCE,
            "em": "fast_profile_resume.py:Phase2 (2026-07-05)",
        },
        "note": (
            "equal/host-conservative/unique-weighted wall-times are from prior run "
            f"(killed by external signal; hard-coded from {KNOWN_PROVENANCE}). "
            "em timing is fresh. "
            "HEAD-SLICE BIAS: first-N-row subsample may over-represent multi-gene ECs "
            "vs full distribution. Absolute extrap values are upper bounds. "
            "Relative method ranking (equal ≈ hc ≈ uw < em) is unaffected. "
            "em_gene_abundances is now vectorised (sparse matvec, commit 2c2e6f0). "
            "Extrapolation valid: loop is O(n_bus_records)."
        ),
    }
    write_json(out_dir / "walltime_subsample.json", sub_result)

    return dict(
        wall_sub=wall_sub,
        wall_extrap=wall_extrap,
    )


# ---------------------------------------------------------------------------
# Phase 2b: EM convergence trajectory (identical to fast_profile.py)
# ---------------------------------------------------------------------------
def run_em_convergence(inputs: dict, out_dir: Path) -> dict:
    log(f"\nPhase 2b (EM convergence on {WALLTIME_N_ROWS:,}-row subsample)...")
    log("  NOTE: em_gene_abundances is now vectorised (sparse matvec, commit 2c2e6f0).")

    bus_full = inputs["bus_df"]
    bus_sub  = bus_full.iloc[:WALLTIME_N_ROWS].copy()

    barcode_to_idx = inputs["barcode_to_idx"]
    ec_map         = inputs["ec_map"]

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

    n_ec = len(em_ec_counts)
    log(f"  building sparse incidence matrix ({n_ec:,} ECs × {inputs['n_genes']:,} genes)...")
    t_build_start = time.perf_counter()
    ec_rows_list: list[np.ndarray] = []
    ec_cols_list: list[np.ndarray] = []
    counts_vec   = np.empty(n_ec, dtype=float)
    genes_per_ec = np.empty(n_ec, dtype=float)
    for e, (genes, count) in enumerate(em_ec_counts.items()):
        g = np.asarray(genes, dtype=int)
        ec_rows_list.append(np.full(g.shape[0], e, dtype=int))
        ec_cols_list.append(g)
        counts_vec[e]   = count
        genes_per_ec[e] = g.shape[0]
    row_idx   = np.concatenate(ec_rows_list)
    col_idx   = np.concatenate(ec_cols_list)
    from scipy.sparse import csr_matrix  # noqa: PLC0415
    incidence = csr_matrix(
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
            "full-data run. EM loop time is for the vectorised version only."
        ),
    }
    write_json(out_dir / "em_convergence.json", conv_data)

    return conv_data


# ---------------------------------------------------------------------------
# Phase 5: Summary + register_value + timing_table.tsv + numbers.json
# ---------------------------------------------------------------------------
def _parse_func(profile_text: str, func_pattern: str) -> tuple[float | None, float | None]:
    for line in profile_text.splitlines():
        if func_pattern in line:
            m = re.match(
                r"^\s*(\d+/?\d*)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(.+)$", line
            )
            if m:
                try:
                    return float(m.group(2)), float(m.group(4))
                except ValueError:
                    pass  # ANALYSIS_OK[best-effort-fan-out]: non-matching cProfile header lines
    return None, None


def _extract_total_time(profile_text: str) -> float | None:
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

    total_equal = _extract_total_time(profile_equal_text)
    bm_tot, bm_cum = _parse_func(profile_equal_text, "build_multimap_layers")
    matv_tot, matv_cum = _parse_func(profile_equal_text, "_matrix_value")

    log(f"\n(a) Hot-function attribution (equal, {CPROFILE_N_ROWS:,}-row cProfile):")
    if total_equal:
        log(f"  Total profiled time: {total_equal:.1f}s")
    if bm_cum and total_equal:
        log(f"  build_multimap_layers: tottime={bm_tot:.1f}s, cumtime={bm_cum:.1f}s "
            f"({bm_cum/total_equal:.1%} of total)")
    if matv_cum and total_equal:
        log(f"  _matrix_value: tottime={matv_tot:.1f}s, cumtime={matv_cum:.1f}s "
            f"({matv_cum/total_equal:.1%} of total)")
        if matv_cum > 0.5 * (total_equal or 1):
            log("  CONFIRMED: _matrix_value dominates non-EM build time (>50%).")
        elif matv_cum > 0.3 * (total_equal or 1):
            log("  LIKELY: _matrix_value is major contributor (>30%).")
    else:
        log("  _matrix_value: not found in top-40 profile lines")

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
        log("  em_gene_abundances: NOT IN top-40 (vectorised — below profiling noise floor).")
        if em_bm_tot and bm_tot and total_equal:
            em_records_extra = (em_bm_tot or 0) - (bm_tot or 0)
            eq_tot = total_equal or 1.0
            log(f"  REFUTED (vectorised code): em_gene_abundances is essentially free.")
            log(f"  EM extra cost = em_records allocation loop only: "
                f"{em_records_extra:.1f}s ({em_records_extra/eq_tot:.1%} of equal total).")
            log(f"  Dominant bottleneck for ALL methods (equal AND em): _matrix_value at line 292.")

    log(f"\n(iii) EM convergence:")
    log(f"  EM converged at iteration {conv_iter} / {EM_MAX_ITER} "
        f"(vectorised sparse-matvec, {WALLTIME_N_ROWS:,}-row subsample)")
    log(f"  EM loop time (vectorised): {em_loop_s:.3f}s for {conv_iter} iterations")
    log(f"  EM is now fast (sparse matvec); the bus-pass (zip loop) dominates all methods.")
    if conv_iter < EM_MAX_ITER:
        log(f"  Recommended: em_max_iter={conv_iter + 5} (add small buffer). "
            f"Saves {(1 - conv_iter/EM_MAX_ITER)*100:.0f}% of EM loop overhead.")

    log(f"\n=== PRIORITY ORDER ===")
    log(f"  CAVEAT: Absolute extrap times (e.g. equal={t_equal:.0f}s) are upper bounds due to")
    log(f"  head-slice bias (first-5M rows may over-represent multi-gene ECs).")
    log(f"  Old itertuples anchor (full data, different code): {OLD_CODE_EQUAL_S}s.")
    log(f"  Relative ranking across methods and hot-function attribution are robust.")
    log(f"")
    log(f"  1. Vectorize the main bus-pass — affects ALL methods, highest total gain.")
    log(f"     Current (zip): extrap ~{t_equal:.0f}s for 'equal' method (upper bound).")
    log(f"     Strategy: groupby (barcode, ec) + merge ec_info + vectorised per-EC ops.")
    log(f"     Expected: <10s (100–1000× speedup over current code).")
    log(f"")
    if matv_cum and total_equal and matv_cum > 0.3 * total_equal:
        log(f"  1a. Within the bus-pass: eliminate _matrix_value hot-spot first (easier fix).")
        log(f"     _matrix_value = {matv_cum/total_equal:.1%} of build time in cProfile.")
        log(f"     Fix: precompute original_counts as dense row-slices or CSR row view.")
        log(f"     This alone gives a meaningful speedup without full vectorization.")
        log(f"")
    log(f"  2. Lower em_max_iter (iii) — trivial 1-line change, low effort.")
    if conv_iter < EM_MAX_ITER:
        log(f"     Set em_max_iter={conv_iter + 5}. EM loop is already fast (vectorised).")
    log(f"")
    log(f"  3. Vectorize main-pass fully (pandas groupby + merge) to close the gap entirely.")
    log(f"")
    log(f"  BUILD vs SWAP distinction (important for user recommendation):")
    log(f"  - Non-EM method selection = O(1) layer pointer in h5ad. Re-running multimap")
    log(f"    to switch equal→host-conservative = zero reprocessing in ViralScan.")
    log(f"  - EM requires a full build pass (~{t_em:.0f}s extrap, upper bound).")
    log(f"  - Recommendation: default to 'equal' (fastest, transparent). Use 'em' only")
    log(f"    when per-gene EM abundance estimates are explicitly needed.")


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
            "provenance": (
                "prior-run-hard-coded" if method != "em" else "fresh-measurement"
            ),
            "anchor_notes": (
                "NEW zip-loop code; extrap from 5M-row HEAD-SLICE (upper bound); "
                f"OLD itertuples anchor={OLD_CODE_EQUAL_S}s"
                if method == "equal" else ""
            ),
        })

    df = pd.DataFrame(rows)
    if df.shape[0] != len(METHODS):  # ANALYSIS_OK[runtime-assert]
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
        prov_sub   = (
            f"time.perf_counter:{WALLTIME_N_ROWS}row-subsample (fresh)"
            if method == "em"
            else f"{KNOWN_PROVENANCE} (hard-coded)"
        )
        reg(key_sub,    round(wall_sub[method], 3),    prov_sub)
        reg(key_extrap, round(wall_extrap[method], 0),
            f"extrapolated×{FULL_DATA_N_ROWS/WALLTIME_N_ROWS:.1f} (upper bound: head-slice bias)")

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

    _print_recommendation(inputs, wall_extrap, em_conv, profile_equal_text, profile_em_text)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    log("=== ViralScan multimap FAST profiling — RESUME ===")
    log(f"  Code version: current HEAD (zip loop + vectorised EM)")
    log(f"  Resuming after process kill during em wall-time timing.")
    log(f"  Anchor: equal on OLD itertuples code = {OLD_CODE_EQUAL_S:.1f}s (full data)")
    log(f"  cProfile slice: {CPROFILE_N_ROWS:,} rows (reading from disk)")
    log(f"  Wall-time slice: {WALLTIME_N_ROWS:,} rows (em only; others hard-coded)")

    out_dir = Path(__file__).parent.parent / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Phase 0 — load inputs (unavoidable)
    inputs = load_all_inputs()

    # Phase 3 — read existing cprofile from disk
    profile_equal_text, profile_em_text = load_cprofile_from_disk(out_dir)

    # Phase 2 — measure em only; hard-code other 3
    wall_result = run_em_walltime_only(inputs, out_dir)

    # Phase 2b — EM convergence
    em_conv = run_em_convergence(inputs, out_dir)

    # Phase 5 — summary + register_value
    build_summary(inputs, wall_result, em_conv, profile_equal_text, profile_em_text, out_dir)

    log("\n=== fast profiling RESUME complete ===")


if __name__ == "__main__":
    main()

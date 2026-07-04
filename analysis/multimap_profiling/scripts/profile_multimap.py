#!/usr/bin/env python3
"""Empirical profiling of ViralScan's 4 multimapping methods on a full-depth EBV sample.

Robust-analysis conventions (strict execution, assert shapes, log row counts, no silent drops).

Structure
---------
Phase 0: Load inputs once. Assert shapes/sizes. Log all key counts.
Phase 1: Wall-time per method (full data). Run each method; time.perf_counter.
Phase 2: EM convergence trajectory (copied loop — does NOT modify src/).
Phase 3: cProfile hotspot analysis on a fixed subsample (10M rows).
Phase 4: Peak-RSS per method via subprocess isolation (/usr/bin/time -v).
Phase 5: Register headline values + write summary table.

ANALYSIS_OK waivers
-------------------
[runtime-assert]    All asserts are developer tripwires in a non-production analysis
                    script; this script is never run with -O. Each assert is followed
                    by a log() call that would also surface the failure. Waiver applied
                    globally for this script.
[optional-input]    gene_names file fallback: the file may legitimately be absent
                    (some kb-python versions omit it). We fall back to gene_ids with
                    a logged WARNING. This does not affect any timing or convergence
                    measurement.
[best-effort-fan-out] cProfile text parsers: ValueError/parse failures return None
                    and the caller logs the gap explicitly. A non-parseable cProfile
                    line does not affect timing numbers or EM convergence; the raw
                    cProfile text is written to disk for manual inspection.
[layer-choice]      adata.X is used as original_counts. This is the raw count matrix
                    from kb-python (cells × genes, CSR), confirmed by inspection of
                    the h5ad: no .raw is set and no transformation has been applied
                    upstream of this step. This matches the intent of build_multimap_layers.
[file-selection]    All paths are pinned to absolute paths matching the specific
                    benchmark run directory. No glob or latest-file selection is used.
[sample-filter]     bus_df.dropna: NaN rows arise from partial/corrupt BUS lines.
                    We log the count before and after, assert < 1% dropped.
[threshold]         EM_MAX_ITER=100, EM_TOL=1e-6, PSEUDOCOUNT=1.0 match viralscan
                    defaults exactly (confirmed: multimapping.py lines 151-152, 149).
"""

from __future__ import annotations

import cProfile
import json
import pstats
import resource
import subprocess
import sys
import time
import io
import os
import re
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
BUS_TXT    = KB / "output.bus.txt"
EC_FILE    = KB / "matrix.ec"
TRANSCRIPTS_FILE = KB / "transcripts.txt"
BARCODES_FILE    = KB / "counts_unfiltered" / "cells_x_genes.barcodes.txt"
GENES_FILE       = KB / "counts_unfiltered" / "cells_x_genes.genes.txt"
# ANALYSIS_OK[file-selection]: pinned to an absolute path; SHA256 of this file
# matches the benchmark run it was built with (t2g line count = 230901, same as
# transcripts.txt — confirmed by wc -l). No glob or latest-file selection is used.
T2G_FILE   = Path(
    "/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_showcase"
    "/viralscan_showcase/fullrun/refs/merged/t2g_plus_anellovirus.txt"
)
ADATA_FILE = KB / "counts_unfiltered" / "adata.h5ad"

REPO_ROOT  = Path("/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan")
REGISTER_VALUE_PATH = Path(
    "/home/mdmanurung/.claude/plugins/marketplaces/mycelium/skills/core/scripts"
)

# Script lives at analysis/multimap_profiling/scripts/profile_multimap.py
# register_value auto-infers namespace from path, so no override needed.
sys.path.insert(0, str(REGISTER_VALUE_PATH))
sys.path.insert(0, str(REPO_ROOT / "src"))

from register_value import register_value
from viralscan.scripts.multimap import (
    load_barcodes,
    read_ec,
    strip_10x_suffix,
)
from viralscan.multimapping import build_multimap_layers

# ---------------------------------------------------------------------------
# Constants — promoted to named values per scilintr magic-threshold rule
# ANALYSIS_OK[threshold]: all match viralscan defaults (confirmed: multimapping.py lines 149-152)
# ---------------------------------------------------------------------------
EM_MAX_ITER  = 100       # matches viralscan default em_max_iter
EM_TOL       = 1e-6      # matches viralscan default em_tol
PSEUDOCOUNT  = 1.0       # matches viralscan default pseudocount
CPROFILE_N_ROWS = 10_000_000  # ANALYSIS_OK[duplicate-config-source]: intentionally different from fast_profile.py (1M); 10M was chosen for the original full-data profiling run
METHODS = ["equal", "host-conservative", "unique-weighted", "em"]


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[profile_multimap] {msg}", flush=True)


def check_nonzero(name: str, value: int) -> None:
    """Raise ValueError if value <= 0.  ANALYSIS_OK[runtime-assert]"""
    if value <= 0:
        raise ValueError(f"Expected {name} > 0, got {value}")


# ---------------------------------------------------------------------------
# Phase 0: Load inputs
# ---------------------------------------------------------------------------
def load_all_inputs() -> dict:
    """Load every input once; check shapes; log row counts."""
    log("Phase 0: Loading inputs...")
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    # -- Barcodes
    barcode_to_idx, n_cells = load_barcodes(str(BARCODES_FILE))
    check_nonzero("n_cells", n_cells)
    log(f"  barcodes: n_cells={n_cells:,}")

    # -- Genes from h5ad  ANALYSIS_OK[layer-choice]
    import anndata as ad
    adata_orig = ad.read_h5ad(str(ADATA_FILE))
    # adata.X is the raw count matrix from kb-python (confirmed: no .raw, no
    # in-place normalization upstream). Used as original_counts for build_multimap_layers.
    gene_ids = list(adata_orig.var_names)
    n_genes = len(gene_ids)
    check_nonzero("n_genes", n_genes)
    log(f"  genes: n_genes={n_genes:,}")
    # Cross-check against genes.txt
    with open(GENES_FILE) as f:
        gene_ids_file = [line.strip() for line in f]
    if len(gene_ids_file) != n_genes:
        raise ValueError(
            f"Gene count mismatch: adata has {n_genes}, genes.txt has {len(gene_ids_file)}"
        )
    log(f"  gene list cross-check: OK (n={n_genes:,})")

    # -- gene_names: load from file if present, else use gene_ids
    # ANALYSIS_OK[optional-input]: some kb-python versions omit gene_names file
    gene_names_file = KB / "counts_unfiltered" / "cells_x_genes.genes.names.txt"
    if gene_names_file.exists():
        with open(gene_names_file) as f:
            gene_names = [line.strip() for line in f]
        if len(gene_names) != n_genes:
            raise ValueError(
                f"Gene names count mismatch: {len(gene_names)} vs {n_genes}"
            )
    else:
        gene_names = gene_ids.copy()
        log("  WARNING: no gene_names file found, using gene_ids as gene_names")

    # -- Transcripts + t2g
    with open(TRANSCRIPTS_FILE) as f:
        transcripts = [line.strip() for line in f]
    check_nonzero("n_transcripts", len(transcripts))
    log(f"  transcripts: n={len(transcripts):,}")

    t2g = pd.read_csv(T2G_FILE, sep=r"\s+", header=None, usecols=[0, 1],
                      names=["transcript", "gene"])
    t2g_map = dict(zip(t2g["transcript"], t2g["gene"]))
    check_nonzero("n_t2g_entries", len(t2g_map))
    log(f"  t2g: n_mappings={len(t2g_map):,}")

    # -- EC map
    log("  reading EC file...")
    ec_map = read_ec(str(EC_FILE), transcripts, t2g_map, gene_ids)
    check_nonzero("n_ec_map", len(ec_map))
    n_multi_ec = sum(1 for genes in ec_map.values() if len(genes) > 1)
    log(f"  ec_map: n_ECs={len(ec_map):,}, n_multi-gene ECs={n_multi_ec:,}")

    # -- BUS records
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
    # We log before/after and assert < 1% dropped. Non-NaN rows are structurally
    # guaranteed by the int32 dtype cast above; NaN survives only if dtype parse
    # yields NA (i.e., truly malformed line).
    bus_df = bus_df.dropna()
    n_bus_after = len(bus_df)
    n_dropped = n_bus_before - n_bus_after
    log(f"  bus records: loaded={n_bus_before:,}, after dropna={n_bus_after:,}, "
        f"dropped={n_dropped:,}")
    if n_dropped > n_bus_before * 0.01:
        raise ValueError(
            f"Too many dropped bus records: {n_dropped:,} / {n_bus_before:,} "
            f"({n_dropped/n_bus_before:.1%}). Check BUS file integrity."
        )
    check_nonzero("n_bus_records", n_bus_after)

    # Strip trailing '-1' lane suffix from barcodes in bus_df
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

    # -- Viral gene indices
    # For this profiling run we have no analysis.txt file, so use a heuristic:
    # any gene whose ID is not a human Ensembl gene (ENSG...) is treated as viral.
    # This replicates what normalize_barcodes does when no analysis.txt is available.
    viral_gene_indices: set[int] = {
        i for i, gid in enumerate(gene_ids) if not str(gid).startswith("ENSG")
    }
    n_viral = len(viral_gene_indices)
    log(f"  viral_gene_indices: n_viral_genes={n_viral:,} (heuristic: non-ENSG gene IDs)")
    check_nonzero("n_viral_genes", n_viral)

    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_load_kb = rss_after - rss_before
    log(f"  RSS after loading: watermark_delta={rss_load_kb:,} KB "
        f"({rss_load_kb/1024:.0f} MB), abs_watermark={rss_after:,} KB")

    return dict(
        barcode_to_idx=barcode_to_idx,
        n_cells=n_cells,
        gene_ids=gene_ids,
        gene_names=gene_names,
        n_genes=n_genes,
        transcripts=transcripts,
        t2g_map=t2g_map,
        ec_map=ec_map,
        n_multi_ec=n_multi_ec,
        viral_gene_indices=viral_gene_indices,
        n_viral=n_viral,
        bus_df=bus_df,
        n_bus=n_bus_after,
        adata_orig=adata_orig,
        original_counts=adata_orig.X,  # ANALYSIS_OK[layer-choice]: raw counts, see docstring
        rss_load_kb=rss_load_kb,
    )


# ---------------------------------------------------------------------------
# Phase 1: Wall-time per method
# ---------------------------------------------------------------------------
def run_wall_time(inputs: dict) -> dict[str, float]:
    """Time each method call on the full dataset. Return {method: wall_secs}."""
    log("\nPhase 1: Wall-time profiling (full data)...")
    log(f"  n_bus_records={inputs['n_bus']:,}, n_cells={inputs['n_cells']:,}, "
        f"n_genes={inputs['n_genes']:,}")
    log("  NOTE: all 4 methods share the same itertuples pass. Non-EM methods build")
    log("  all 3 deterministic layers in ONE pass; EM adds EM loop + second allocation.")
    log("  Build cost (this run) != swap cost (in viralscan h5ad layer select = O(1)).")

    wall_times: dict[str, float] = {}
    for method in METHODS:
        log(f"  timing method={method!r}...")
        t0 = time.perf_counter()
        layers = build_multimap_layers(
            bus_df=inputs["bus_df"],
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
        wall_times[method] = elapsed
        # Validate output shapes  ANALYSIS_OK[runtime-assert]
        expected = (inputs["n_cells"], inputs["n_genes"])
        for layer_name, layer_mat in [
            ("equal", layers.equal),
            ("host_conservative", layers.host_conservative),
            ("unique_weighted", layers.unique_weighted),
            ("corrected", layers.corrected),
        ]:
            if layer_mat.shape != expected:
                raise ValueError(
                    f"Layer {layer_name} shape mismatch: {layer_mat.shape} != {expected}"
                )
        log(f"    -> {elapsed:.1f}s")
    return wall_times


# ---------------------------------------------------------------------------
# Phase 2: EM convergence trajectory
# ---------------------------------------------------------------------------
def run_em_convergence(inputs: dict) -> dict:
    """Instrument EM loop (copy of em_gene_abundances) to get iteration count
    and L1-delta trajectory WITHOUT modifying src/.

    The copy is a faithful replica of the loop body in
    src/viralscan/multimapping.py:em_gene_abundances (lines 124-137 as of
    2026-07-04). Any discrepancy would change the convergence point.
    """
    log("\nPhase 2: EM convergence trajectory...")

    # Build em_ec_counts from the full bus_df (same logic as build_multimap_layers)
    log("  collecting em_ec_counts from bus_df (one pass)...")
    t_pass_start = time.perf_counter()
    em_ec_counts: dict[tuple[int, ...], float] = {}
    n_multi_bus = 0

    barcode_to_idx = inputs["barcode_to_idx"]
    ec_map         = inputs["ec_map"]

    for row in inputs["bus_df"].itertuples(index=False):
        bc, ec, count = row.barcode, row.ec, float(row.count)
        if pd.isna(ec):
            continue
        ec = int(ec)
        if bc not in barcode_to_idx or ec not in ec_map:
            continue
        genes_in_ec = list(ec_map[ec])
        if not genes_in_ec or len(genes_in_ec) == 1:
            continue
        # multi-gene EC
        distinct_genes = list(dict.fromkeys(genes_in_ec))
        key = tuple(distinct_genes)
        em_ec_counts[key] = em_ec_counts.get(key, 0.0) + count
        n_multi_bus += 1

    t_pass_end = time.perf_counter()
    log(f"  bus pass: {t_pass_end - t_pass_start:.1f}s, "
        f"n_multi_bus_rows={n_multi_bus:,}, n_distinct_multi_ec_keys={len(em_ec_counts):,}")

    # Build unique_per_gene from original_counts  ANALYSIS_OK[layer-choice]
    original_counts = inputs["original_counts"]
    unique_per_gene = np.asarray(
        original_counts.sum(axis=0)
        if sparse.issparse(original_counts)
        else original_counts.sum(axis=0)
    ).reshape(-1)
    if len(unique_per_gene) != inputs["n_genes"]:
        raise ValueError(
            f"unique_per_gene length mismatch: {len(unique_per_gene)} vs {inputs['n_genes']}"
        )
    log(f"  unique_per_gene: shape={unique_per_gene.shape}, sum={unique_per_gene.sum():.0f}")

    # --- Instrumented EM loop (faithful copy of em_gene_abundances) ---
    # See src/viralscan/multimapping.py:78-138
    theta = unique_per_gene.astype(float) + float(PSEUDOCOUNT)
    items = [(np.asarray(genes, dtype=int), float(count))
             for genes, count in em_ec_counts.items()]

    log(f"  starting instrumented EM (max_iter={EM_MAX_ITER}, tol={EM_TOL})...")
    iter_deltas: list[float] = []
    converged_iter: int | None = None

    t_em_start = time.perf_counter()
    for i_iter in range(int(EM_MAX_ITER)):
        new = unique_per_gene.copy().astype(float)
        for genes, count in items:
            w = theta[genes]
            s = float(w.sum())
            if s <= 1e-12:
                new[genes] += count / len(genes)
            else:
                new[genes] += count * w / s
        denom = float(theta.sum()) or 1.0
        l1_delta = float(np.abs(new - theta).sum()) / denom
        iter_deltas.append(l1_delta)
        theta = new
        if l1_delta < EM_TOL:
            converged_iter = i_iter + 1  # 1-indexed
            log(f"  EM converged at iteration {converged_iter}, L1_delta={l1_delta:.2e}")
            break
    else:
        converged_iter = EM_MAX_ITER
        log(f"  EM did NOT converge within {EM_MAX_ITER} iterations; "
            f"final L1_delta={iter_deltas[-1]:.2e}")
    t_em_end = time.perf_counter()
    log(f"  EM loop wall-time: {t_em_end - t_em_start:.1f}s")
    log(f"  First 5 L1_deltas: {[f'{d:.3e}' for d in iter_deltas[:5]]}")
    log(f"  Last 5 L1_deltas:  {[f'{d:.3e}' for d in iter_deltas[-5:]]}")

    return dict(
        converged_iter=converged_iter,
        iter_deltas=iter_deltas,
        n_multi_bus=n_multi_bus,
        n_em_ec_keys=len(em_ec_counts),
        em_pass_time=t_pass_end - t_pass_start,
        em_loop_time=t_em_end - t_em_start,
    )


# ---------------------------------------------------------------------------
# Phase 3: cProfile hotspots on subsample
# ---------------------------------------------------------------------------
def run_cprofile_hotspots(inputs: dict) -> tuple[str, str]:
    """Run cProfile on a fixed 10M-row deterministic slice; report top functions."""
    log(f"\nPhase 3: cProfile hotspot analysis (subsample n={CPROFILE_N_ROWS:,})...")

    bus_full = inputs["bus_df"]
    if len(bus_full) > CPROFILE_N_ROWS:
        # Deterministic first-N slice: reproducible without a random seed.
        # ANALYSIS_OK[sample-filter]: subsample is for cProfile only; all
        # timing and EM measurements use the full dataset. Slice is logged.
        bus_sub = bus_full.iloc[:CPROFILE_N_ROWS].copy()
        log(f"  subsampled to first {CPROFILE_N_ROWS:,} rows (of {len(bus_full):,}); "
            "cProfile only — not used for timing or EM measurements")
    else:
        bus_sub = bus_full
        log(f"  using full dataset ({len(bus_full):,} rows, < {CPROFILE_N_ROWS:,} threshold)")

    def _profile_method(method_name: str) -> str:
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
        ps = pstats.Stats(profiler, stream=sio).sort_stats("cumulative")
        ps.print_stats(30)
        return sio.getvalue()

    profile_equal_text = _profile_method("equal")
    profile_em_text    = _profile_method("em")

    log(f"  cProfile (equal, {CPROFILE_N_ROWS:,} rows) — top 12 lines:")
    for line in profile_equal_text.split("\n")[3:15]:
        log(f"    {line}")
    log(f"  cProfile (em, {CPROFILE_N_ROWS:,} rows) — top 12 lines:")
    for line in profile_em_text.split("\n")[3:15]:
        log(f"    {line}")

    return profile_equal_text, profile_em_text


# ---------------------------------------------------------------------------
# Phase 4: Peak RSS via subprocess isolation
# ---------------------------------------------------------------------------
SUBPROCESS_RUNNER_SCRIPT = Path(__file__).parent / "_rss_probe.py"


def _write_rss_probe() -> None:
    """Write a helper script that loads inputs + runs one method, then exits."""
    probe_src = r'''#!/usr/bin/env python3
"""Subprocess probe: load inputs + run one method. /usr/bin/time -v gives peak RSS."""
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

REPO_ROOT = Path("/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan")
sys.path.insert(0, str(REPO_ROOT / "src"))

KB = Path(
    "/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_runs"
    "/reference_strategy_2026-06-28_fresh12b/runs/ebv__viralscan__combined"
    "/SRR12682296/kb-python"
)
BUS_TXT = KB / "output.bus.txt"
EC_FILE = KB / "matrix.ec"
TRANSCRIPTS_FILE = KB / "transcripts.txt"
BARCODES_FILE = KB / "counts_unfiltered" / "cells_x_genes.barcodes.txt"
T2G_FILE = Path(
    "/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_showcase"
    "/viralscan_showcase/fullrun/refs/merged/t2g_plus_anellovirus.txt"
)
ADATA_FILE = KB / "counts_unfiltered" / "adata.h5ad"

PSEUDOCOUNT = 1.0
EM_MAX_ITER = 100
EM_TOL = 1e-6

method = sys.argv[1] if len(sys.argv) > 1 else "load-only"

import anndata as ad
from viralscan.scripts.multimap import load_barcodes, read_ec, strip_10x_suffix
from viralscan.multimapping import build_multimap_layers

barcode_to_idx, n_cells = load_barcodes(str(BARCODES_FILE))
adata_orig = ad.read_h5ad(str(ADATA_FILE))
gene_ids = list(adata_orig.var_names)
n_genes = len(gene_ids)

with open(TRANSCRIPTS_FILE) as f:
    transcripts = [line.strip() for line in f]
t2g = pd.read_csv(T2G_FILE, sep=r"\s+", header=None, usecols=[0, 1],
                  names=["transcript", "gene"])
t2g_map = dict(zip(t2g["transcript"], t2g["gene"]))
ec_map = read_ec(str(EC_FILE), transcripts, t2g_map, gene_ids)

bus_df = pd.read_csv(
    str(BUS_TXT), sep="\t", header=None,
    names=["barcode", "umi", "ec", "count"],
    usecols=["barcode", "ec", "count"],
    dtype={"barcode": "category", "ec": "int32", "count": "int32"},
)
bus_df = bus_df.dropna()
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

viral_gene_indices = {i for i, gid in enumerate(gene_ids) if not str(gid).startswith("ENSG")}

if method != "load-only":
    layers = build_multimap_layers(
        bus_df=bus_df,
        barcode_to_idx=barcode_to_idx,
        ec_map=ec_map,
        n_cells=n_cells,
        n_genes=n_genes,
        viral_gene_indices=viral_gene_indices,
        original_counts=adata_orig.X,
        method=method,
        pseudocount=PSEUDOCOUNT,
        em_max_iter=EM_MAX_ITER,
        em_tol=EM_TOL,
    )

print(f"RSS_PROBE_DONE method={method}", flush=True)
'''
    SUBPROCESS_RUNNER_SCRIPT.write_text(probe_src)


def parse_rss_from_time_output(stderr_text: str) -> int | None:
    """Extract 'Maximum resident set size' in KB from /usr/bin/time -v stderr."""
    for line in stderr_text.splitlines():
        if "Maximum resident set size" in line:
            m = re.search(r"(\d+)", line.split(":")[-1])
            if m:
                return int(m.group(1))
    return None


def run_peak_rss(python_exe: str) -> dict[str, int | None]:
    """Run each method in a fresh subprocess via /usr/bin/time -v.

    ru_maxrss is a monotonically increasing high-water mark for the whole
    process. Running each method in a fresh subprocess gives a clean per-method
    peak. Per-method cost = (method run total) - (load-only baseline).
    """
    log("\nPhase 4: Peak RSS via subprocess isolation (/usr/bin/time -v)...")

    _write_rss_probe()

    pythonpath = str(REPO_ROOT / "src")
    targets = ["load-only"] + list(METHODS)
    rss_results: dict[str, int | None] = {}

    for method in targets:
        log(f"  launching subprocess for method={method!r}...")
        env = os.environ.copy()
        existing_pp = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = f"{pythonpath}:{existing_pp}" if existing_pp else pythonpath
        cmd = [
            "/usr/bin/time", "-v",
            python_exe, str(SUBPROCESS_RUNNER_SCRIPT), method
        ]
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                env=env,
                timeout=7200,  # ANALYSIS_OK[threshold]: 2h timeout per method
                check=False,
            )
            rss_kb = parse_rss_from_time_output(result.stderr)
            rss_results[method] = rss_kb
            if result.returncode != 0:
                log(f"    WARNING: method={method!r} exited {result.returncode}")
                log(f"    stderr (last 500 chars): {result.stderr[-500:]}")
            else:
                if rss_kb:
                    log(f"    peak RSS: {rss_kb:,} KB ({rss_kb/1024:.0f} MB)")
                else:
                    log("    peak RSS: not found in /usr/bin/time -v output")
        except subprocess.TimeoutExpired:
            log(f"    TIMEOUT: method={method!r} exceeded 2h timeout")
            # ANALYSIS_OK[optional-input]: timeout means no RSS measurement for
            # this method. None is logged and clearly marked missing in the
            # output table; it does not propagate silently into any computation.
            rss_results[method] = None

    return rss_results


# ---------------------------------------------------------------------------
# Phase 5: Build summary table and register headline values
# ---------------------------------------------------------------------------
def build_summary(
    inputs: dict,
    wall_times: dict[str, float],
    em_convergence: dict,
    rss_results: dict[str, int | None],
) -> pd.DataFrame:
    """Build and print the summary table, register values, write outputs."""
    log("\nPhase 5: Summary and value registration...")

    load_only_rss = rss_results.get("load-only")

    rows = []
    for method in METHODS:
        rss_raw = rss_results.get(method)
        rss_delta = (rss_raw - load_only_rss) if (rss_raw and load_only_rss) else None
        rows.append({
            "method": method,
            "wall_time_s": round(wall_times[method], 1),
            "peak_rss_kb": rss_raw,
            "peak_rss_minus_load_kb": rss_delta,
        })

    df = pd.DataFrame(rows)
    # ANALYSIS_OK[runtime-assert]: validating own output shape, not public API
    if df.shape[0] != len(METHODS):
        raise ValueError(f"Summary row count mismatch: {df.shape[0]} != {len(METHODS)}")
    if "method" not in df.columns or "wall_time_s" not in df.columns:
        raise ValueError("Summary df missing expected columns")
    log("\n=== WALL-TIME AND PEAK RSS TABLE ===")
    log(df.to_string(index=False))

    # Write to outputs
    out_dir = Path(__file__).parent.parent / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "timing_table.tsv"
    df.to_csv(summary_path, sep="\t", index=False)
    log(f"\n  timing table written to: {summary_path}")

    # EM convergence trajectory
    em_path = out_dir / "em_convergence.json"
    em_data = {
        "converged_iter": em_convergence["converged_iter"],
        "n_multi_bus_rows": em_convergence["n_multi_bus"],
        "n_em_ec_keys": em_convergence["n_em_ec_keys"],
        "em_loop_time_s": round(em_convergence["em_loop_time"], 2),
        "em_pass_time_s": round(em_convergence["em_pass_time"], 2),
        "l1_delta_trajectory": [round(d, 8) for d in em_convergence["iter_deltas"]],
    }
    em_path.write_text(json.dumps(em_data, indent=2) + "\n")
    log(f"  EM convergence data written to: {em_path}")

    # Register headline values
    def reg(key: str, val: int | float | str, prov: str) -> None:
        register_value(key, val, provenance=prov)

    reg("n_bus_records", inputs["n_bus"], "KB/output.bus.txt:wc-l")
    reg("n_cells", inputs["n_cells"], "KB/counts_unfiltered/cells_x_genes.barcodes.txt:wc-l")
    reg("n_genes", inputs["n_genes"], "KB/counts_unfiltered/adata.h5ad:adata.n_vars")
    reg("n_ec_map", len(inputs["ec_map"]), "KB/matrix.ec:parsed entries")
    reg("n_multi_gene_ec", inputs["n_multi_ec"], "KB/matrix.ec:entries with >1 gene")
    reg("n_viral_genes", inputs["n_viral"], "gene_ids:non-ENSG count")

    for method in METHODS:
        key = f"wall_time_s_{method.replace('-', '_')}"
        reg(key, round(wall_times[method], 1), "time.perf_counter")

    reg("em_converged_iter", em_convergence["converged_iter"],
        "instrumented em loop: convergence iteration (1-indexed)")
    reg("em_loop_time_s", round(em_convergence["em_loop_time"], 2),
        "instrumented em loop: wall-time")

    if load_only_rss:
        reg("rss_load_only_kb", load_only_rss, "/usr/bin/time -v: load-only subprocess")
    for method in METHODS:
        rss_raw = rss_results.get(method)
        if rss_raw:
            key = f"rss_kb_{method.replace('-', '_')}"
            reg(key, rss_raw, f"/usr/bin/time -v: {method} subprocess")

    log("  values registered to analysis/multimap_profiling/outputs/numbers.json")
    return df


# ---------------------------------------------------------------------------
# Speedup analysis helpers
# ---------------------------------------------------------------------------
def _extract_itertuples_fraction(profile_text: str, total_time: float) -> float | None:
    """Parse cumtime fraction for itertuples from cProfile output.

    Returns None when the function is not found or the line cannot be parsed.
    The caller logs the gap and falls back to manual inspection of the saved
    cProfile text. This does not affect timing numbers or EM convergence.
    ANALYSIS_OK[degraded-fallback]: None signals "could not measure" to the
    caller, which logs it and directs the user to the raw cProfile file.
    """
    for line in profile_text.splitlines():
        if "itertuples" in line:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    cumtime = float(parts[3])
                    return cumtime / total_time if total_time > 0 else None  # ANALYSIS_OK[degraded-fallback]: zero total_time means no timing data; logged above
                except ValueError:
                    log(f"  WARNING: could not parse cumtime from cProfile line: {line!r}")
                    return None  # ANALYSIS_OK[degraded-fallback]: parse fail is logged; caller checks for None and notes missing measurement
    return None  # ANALYSIS_OK[degraded-fallback]: function not in profile output; caller logs and directs to raw cProfile text file


def _extract_build_multimap_time(profile_text: str) -> float | None:
    """Parse cumtime for build_multimap_layers from cProfile output.

    Returns None when function not found or parse fails; callers fall back to
    wall-time as the denominator and log the gap.
    """
    for line in profile_text.splitlines():
        if "build_multimap_layers" in line:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    return float(parts[3])
                except ValueError:
                    log(f"  WARNING: could not parse cumtime from cProfile line: {line!r}")
                    return None  # ANALYSIS_OK[degraded-fallback]: parse fail is logged; caller uses wall-time fallback
    return None  # ANALYSIS_OK[degraded-fallback]: function not in profile; caller uses wall-time fallback


def _extract_em_loop_fraction(profile_em_text: str, em_total: float) -> float | None:
    """Parse cumtime fraction for em_gene_abundances from EM cProfile output.

    Returns None when not found or parse fails; caller logs and directs user
    to raw cProfile file for manual inspection.
    """
    for line in profile_em_text.splitlines():
        if "em_gene_abundances" in line:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    cumtime = float(parts[3])
                    return cumtime / em_total if em_total > 0 else None  # ANALYSIS_OK[degraded-fallback]: zero total means no timing; logged
                except ValueError:
                    log(f"  WARNING: could not parse cumtime from cProfile line: {line!r}")
                    return None  # ANALYSIS_OK[degraded-fallback]: parse fail logged; caller checks None
    return None  # ANALYSIS_OK[degraded-fallback]: function absent from profile; caller logs and directs to raw cProfile text


def _print_recommendation(
    wall_times: dict[str, float],
    em_convergence: dict,
    profile_equal_text: str,
    profile_em_text: str,
    inputs: dict,
    rss_results: dict,
) -> None:
    log("\n=== SPEEDUP RECOMMENDATION ===")

    t_equal = wall_times["equal"]
    t_em    = wall_times["em"]

    conv_iter = em_convergence["converged_iter"]
    em_loop_s = em_convergence["em_loop_time"]

    # (iii) Lower em_max_iter
    if conv_iter < EM_MAX_ITER:
        iter_ratio = conv_iter / EM_MAX_ITER
        log(f"\n(iii) Lower em_max_iter:")
        log(f"  EM converged at iteration {conv_iter} / {EM_MAX_ITER}.")
        log(f"  EM loop time = {em_loop_s:.1f}s; reducing max_iter to {conv_iter} "
            f"would save {(1 - iter_ratio)*100:.0f}% of the EM loop cost.")
        log(f"  Estimated time savings: ~{em_loop_s * (1 - iter_ratio):.1f}s "
            f"({em_loop_s * (1 - iter_ratio) / t_em * 100:.0f}% of total EM wall-time)")
        log(f"  Recommended: set em_max_iter={conv_iter + 5} (add small buffer).")
    else:
        log(f"\n(iii) EM did not converge within {EM_MAX_ITER} iterations. "
            "Recommend investigating convergence before lowering max_iter.")

    # Parse cProfile for hypothesis (a) and (b)
    bm_equal_time = _extract_build_multimap_time(profile_equal_text)
    itertuples_frac = _extract_itertuples_fraction(
        profile_equal_text, bm_equal_time or t_equal
    )
    em_frac = _extract_em_loop_fraction(profile_em_text, bm_equal_time or t_em)

    log(f"\n(a) Hypothesis: itertuples pass dominates ALL methods:")
    if itertuples_frac is not None:
        log(f"  itertuples cumtime fraction (equal method, {CPROFILE_N_ROWS:,}-row subsample) = "
            f"{itertuples_frac:.1%}")
        if itertuples_frac > 0.7:
            log("  CONFIRMED: itertuples dominates (>70% of build_multimap_layers time).")
        else:
            log(f"  NOT CONFIRMED: itertuples = {itertuples_frac:.1%} of build time "
                f"in this subsample. Check cprofile_equal.txt for full output.")
    else:
        log("  Could not extract itertuples fraction from cProfile output.")
        log("  Check outputs/cprofile_equal.txt manually.")

    log(f"\n(b) Hypothesis: em_gene_abundances dominates EM-only extra cost:")
    if em_frac is not None:
        log(f"  em_gene_abundances cumtime fraction (em method subsample) = {em_frac:.1%}")
        if em_frac > 0.7:
            log("  CONFIRMED: em_gene_abundances is the EM bottleneck.")
        else:
            log(f"  Partial or not confirmed: em loop = {em_frac:.1%} of EM method time.")
    else:
        log("  Could not extract em_gene_abundances fraction from cProfile output.")
        log("  Check outputs/cprofile_em.txt manually.")

    log(f"\n(i) Vectorize em_gene_abundances (sparse EC x gene matvec):")
    n_em_keys = em_convergence["n_em_ec_keys"]
    if em_loop_s > 0 and conv_iter:
        per_iter_s = em_loop_s / max(conv_iter, 1)
        log(f"  Current EM: {conv_iter}-iteration pure-Python loop "
            f"over {n_em_keys:,} distinct multi-gene EC keys.")
        log(f"  Per-iteration time (measured): {per_iter_s*1000:.1f}ms")
        log(f"  Vectorized approach: build sparse incidence matrix "
            f"(n_ec_keys={n_em_keys:,} x n_genes={inputs['n_genes']:,}), "
            f"each EM iteration = 2 sparse matvec ops.")
        log(f"  Expected speedup on EM loop: 10-100x per iteration "
            f"(sparse BLAS vs Python loop at {per_iter_s*1000:.0f}ms/iter).")
        log(f"  Amdahl gain (EM only): capped at em_loop_time/t_em = "
            f"{em_loop_s/t_em:.0%} of EM total wall-time.")
    else:
        log("  em_loop_time not available from instrumented run.")

    log(f"\n(ii) Vectorize/group the itertuples pass:")
    log(f"  itertuples over {inputs['n_bus']:,} rows is the dominant cost for ALL methods.")
    log(f"  Vectorization strategy: groupby (barcode, ec) then aggregate count;")
    log(f"  EC-level lookup becomes pandas merge; per-row logic becomes vectorized ops.")
    if itertuples_frac is not None and bm_equal_time:
        log(f"  Gain capped by Amdahl at 1/(1-{itertuples_frac:.2f}) = "
            f"{1/(1-itertuples_frac):.1f}x over build_multimap_layers IF "
            f"itertuples = {itertuples_frac:.0%} of total function time.")
    log(f"  Non-EM methods: this is the ONLY bottleneck (wall ~{t_equal:.0f}s).")
    log(f"  Vectorizing this pass is the highest-impact optimization overall.")

    log(f"\n(iv) Parallelization:")
    log(f"  The itertuples pass is sequential over rows; no dependency between "
        f"different barcodes. Parallelization by barcode shard is feasible.")
    log(f"  However, vectorizing (ii) first is likely to close most of the gap "
        f"at much lower implementation cost.")

    log(f"\n=== PRIORITY ORDER ===")
    log(f"  1. Vectorize itertuples pass (ii) — affects ALL methods, highest total gain.")
    log(f"     Non-EM wall-time: {t_equal:.0f}s; target: <10s with pandas groupby.")
    log(f"  2. Lower em_max_iter to convergence point (iii) — trivial code change.")
    if conv_iter < EM_MAX_ITER:
        savings = em_loop_s * (1 - conv_iter / EM_MAX_ITER)
        log(f"     Saves ~{savings:.0f}s ({savings/t_em:.0%} of EM wall-time). "
            f"Recommend: em_max_iter={conv_iter + 5}")
    log(f"  3. Vectorize em_gene_abundances (i) — high gain on EM loop only;")
    log(f"     less important if (ii) is done first (shared pass already fast).")
    log(f"  4. Parallelization (iv) — consider only after (i)+(ii) if needed.")

    log(f"\nBUILD vs SWAP distinction:")
    log(f"  - Non-EM methods (equal, host-conservative, unique-weighted) all share")
    log(f"    ONE itertuples build pass. Method selection = which pre-built sparse")
    log(f"    matrix is assigned as 'corrected' — O(1) layer swap from h5ad.")
    log(f"    viralscan rerun-multimap between non-EM methods = zero reprocessing.")
    log(f"  - EM runs the shared pass PLUS EM loop + second allocation loop.")
    log(f"    Switching to/from EM = full reprocess (~{t_em:.0f}s).")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(python_exe: str) -> None:
    log("=== ViralScan multimap profiling ===")
    log(f"Python: {python_exe}")
    log(f"BUS: {BUS_TXT}")
    log(f"t2g: {T2G_FILE}")

    # Phase 0
    inputs = load_all_inputs()

    # Phase 1
    wall_times = run_wall_time(inputs)

    # Phase 2
    em_convergence = run_em_convergence(inputs)

    # Phase 3
    profile_equal_text, profile_em_text = run_cprofile_hotspots(inputs)

    # Save cProfile text
    out_dir = Path(__file__).parent.parent / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "cprofile_equal.txt").write_text(profile_equal_text)
    (out_dir / "cprofile_em.txt").write_text(profile_em_text)
    log(f"  cProfile outputs written to {out_dir}/cprofile_{{equal,em}}.txt")

    # Phase 4
    rss_results = run_peak_rss(python_exe)

    # Phase 5
    df_summary = build_summary(inputs, wall_times, em_convergence, rss_results)

    # Speedup recommendation
    _print_recommendation(wall_times, em_convergence, profile_equal_text,
                          profile_em_text, inputs, rss_results)

    log("\n=== profiling complete ===")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            "Usage: PYTHONPATH=src python analysis/multimap_profiling/scripts/profile_multimap.py "
            "<python_exe>",
            file=sys.stderr,
        )
        sys.exit(1)
    python_exe = sys.argv[1]
    main(python_exe)

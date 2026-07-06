#!/usr/bin/env python3
"""Register headline numbers from the multimap_profiling analysis.

Writes analysis/multimap_profiling/outputs/numbers.json in the mycelium
register_value format: a JSON object with "namespace" and "values" list,
each entry having "key", "value", "provenance", and "computed_at".

All numbers are from the BEFORE-FIX code state (b7e9635 + 2c2e6f0 applied;
3c53ad7 NOT yet applied), profiled on SRR12682296 (EBV LCL, 103M BUS rows).

Run from the repo root:
    PYTHONPATH=src python analysis/multimap_profiling/scripts/register_values.py

ANALYSIS_OK[file-selection]: reads hardcoded analysis output paths; intended
    to be run after fast_profile.py has written cprofile_equal.txt.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
OUTPUTS = REPO / "analysis" / "multimap_profiling" / "outputs"
OUT_JSON = OUTPUTS / "numbers.json"
NS = "multimap_profiling"
COMPUTED_AT = "analysis/multimap_profiling/scripts/register_values.py"

# ── Provenance strings ─────────────────────────────────────────────────────────

PROV_LOAD = "analysis/multimap_profiling/outputs/fast_profile_run.log:Phase-0"
PROV_CPROFILE = "analysis/multimap_profiling/outputs/cprofile_equal.txt:pre-3c53ad7-code"
PROV_CPROFILE_EM = "analysis/multimap_profiling/outputs/cprofile_em.txt:pre-3c53ad7-code"
PROV_WALLTIME = "analysis/multimap_profiling/outputs/fast_profile_run.log:lines-57-60:pre-3c53ad7-code"
PROV_ANCHOR = "analysis/multimap_profiling/outputs/profile_run.log:line-24:OLD-itertuples-code"
PROV_FIX = "git:3c53ad7:direct-CSR-buffer-access"

# ── Hardcoded from logged outputs ─────────────────────────────────────────────
# Code version: b7e9635 (zip+hoist) + 2c2e6f0 (vectorise EM) applied;
# 3c53ad7 (direct CSR buffer access, kills _matrix_value __getitem__) NOT yet applied.
# Profiling run started 2026-07-05 00:31 UTC+8.

# Structural facts (code-version independent)
N_BUS_RECORDS = 103_145_071
N_CELLS = 848_191
N_GENES = 43_451
N_ECS = 388_677
N_MULTI_GENE_ECS = 317_685
MULTIMAPPING_RATE = 317_685 / 388_677      # = 0.8172...
N_VIRAL_GENES = 4_845
RSS_AFTER_LOAD_MB = 4_984

# cProfile (1M-row subsample) — BEFORE-FIX code (pre-3c53ad7)
CPROFILE_EQUAL_TOTAL_S = 472.231          # build_multimap_layers cumtime
CPROFILE_MATRIX_VALUE_CUMTIME_S = 408.188 # _matrix_value cumtime
CPROFILE_LISTCOMP_CUMTIME_S = 415.033     # <listcomp> cumtime
CPROFILE_MATRIX_VALUE_NCALLS = 5_873_997
CPROFILE_ISINTLIKE_NCALLS = 23_496_002
CPROFILE_GET_CSR_SUBMATRIX_NCALLS = 5_873_997
MATRIX_VALUE_FRACTION = CPROFILE_MATRIX_VALUE_CUMTIME_S / CPROFILE_EQUAL_TOTAL_S

# Wall-time (5M-row subsample, extrapolated ×20.629) — BEFORE-FIX code (pre-3c53ad7)
# NOTE: first-N-row head-slice bias — multi-gene EC rate may be higher in head slice
# than full BUS file. Absolute extrapolations are upper bounds; relative ranking is robust.
SCALE_FACTOR = N_BUS_RECORDS / 5_000_000   # = 20.629×
WALLTIME_EQUAL_SUB_S = 1139.75
WALLTIME_HC_SUB_S = 1078.88
WALLTIME_UW_SUB_S = 1136.83
# em wall-time on 5M subsample was NEVER measured: fast_profile.py was killed during
# em timing, and the resume script (fast_profile_resume.py) ran em on the CURRENT code
# (post-3c53ad7) — mixing code versions in one table would be misleading. em = NOT AVAILABLE.

WALLTIME_EQUAL_EXTRAP_S = WALLTIME_EQUAL_SUB_S * SCALE_FACTOR
WALLTIME_HC_EXTRAP_S = WALLTIME_HC_SUB_S * SCALE_FACTOR
WALLTIME_UW_EXTRAP_S = WALLTIME_UW_SUB_S * SCALE_FACTOR

# OLD itertuples anchor (profile_multimap.py, equal full data) — entirely different code
# version (no zip/hoist, no vectorised EM, no direct CSR buffer access)
OLD_ITERTUPLES_EQUAL_FULL_S = 19_964.2    # actual full-data measurement

# Fix outcome (3c53ad7 — direct CSR buffer access, replaces _matrix_value __getitem__)
# Speedup measured by prior session (.living/last-session.md, commit bc35036):
# equal full-data: 19,964s (itertuples, pre-all-fixes) → ~5.5h estimate → ~1.4h post-fix
# ~4× overall across itertuples→zip + vectorise-EM + direct-CSR-buffer-access commits
FIX_SPEEDUP_APPROX_X = 4.0   # rough ~4× overall; direct CSR alone ~2.9–3.1×


def _reg(values: list, key: str, value, provenance: str) -> None:
    values.append({
        "key": key,
        "value": value,
        "provenance": provenance,
        "computed_at": COMPUTED_AT,
    })


def main() -> int:
    values: list = []

    # ── Structural facts ─────────────────────────────────────────────────────
    _reg(values, "n_bus_records", N_BUS_RECORDS, PROV_LOAD)
    _reg(values, "n_cells", N_CELLS, PROV_LOAD)
    _reg(values, "n_genes", N_GENES, PROV_LOAD)
    _reg(values, "n_ECs", N_ECS, PROV_LOAD)
    _reg(values, "n_multi_gene_ECs", N_MULTI_GENE_ECS, PROV_LOAD)
    _reg(values, "multimapping_rate", round(MULTIMAPPING_RATE, 4), PROV_LOAD)
    _reg(values, "n_viral_genes", N_VIRAL_GENES, PROV_LOAD)
    _reg(values, "rss_after_load_mb", RSS_AFTER_LOAD_MB, PROV_LOAD)

    # ── cProfile findings (pre-3c53ad7 code) ────────────────────────────────
    _reg(values, "cprofile_equal_total_s", CPROFILE_EQUAL_TOTAL_S, PROV_CPROFILE)
    _reg(values, "cprofile_matrix_value_cumtime_s", CPROFILE_MATRIX_VALUE_CUMTIME_S, PROV_CPROFILE)
    _reg(values, "cprofile_matrix_value_fraction",
         round(MATRIX_VALUE_FRACTION, 4), PROV_CPROFILE)
    _reg(values, "cprofile_matrix_value_ncalls", CPROFILE_MATRIX_VALUE_NCALLS, PROV_CPROFILE)
    _reg(values, "cprofile_isintlike_ncalls", CPROFILE_ISINTLIKE_NCALLS, PROV_CPROFILE)
    _reg(values, "cprofile_get_csr_submatrix_ncalls",
         CPROFILE_GET_CSR_SUBMATRIX_NCALLS, PROV_CPROFILE)
    _reg(values, "cprofile_listcomp_cumtime_s", CPROFILE_LISTCOMP_CUMTIME_S, PROV_CPROFILE)

    # ── Wall-time measurements (pre-3c53ad7 code, 5M-row subsample) ─────────
    _reg(values, "walltime_equal_sub_s", WALLTIME_EQUAL_SUB_S, PROV_WALLTIME)
    _reg(values, "walltime_hc_sub_s", WALLTIME_HC_SUB_S, PROV_WALLTIME)
    _reg(values, "walltime_uw_sub_s", WALLTIME_UW_SUB_S, PROV_WALLTIME)
    _reg(values, "walltime_scale_factor", round(SCALE_FACTOR, 3), PROV_WALLTIME)
    _reg(values, "walltime_equal_extrap_s", round(WALLTIME_EQUAL_EXTRAP_S, 0), PROV_WALLTIME)
    _reg(values, "walltime_hc_extrap_s", round(WALLTIME_HC_EXTRAP_S, 0), PROV_WALLTIME)
    _reg(values, "walltime_uw_extrap_s", round(WALLTIME_UW_EXTRAP_S, 0), PROV_WALLTIME)
    _reg(values, "walltime_em_sub_s", None,
         "em wall-time NOT AVAILABLE: fast_profile.py killed during em timing; "
         "resume (fast_profile_resume.py) ran on post-3c53ad7 code — "
         "mixing code versions in one table would be misleading; em not reported")

    # ── Old-code anchor (itertuples) ─────────────────────────────────────────
    _reg(values, "old_itertuples_equal_full_s", OLD_ITERTUPLES_EQUAL_FULL_S, PROV_ANCHOR)

    # ── Fix outcome ──────────────────────────────────────────────────────────
    _reg(values, "fix_commit", "3c53ad7", PROV_FIX)
    _reg(values, "fix_description",
         "direct CSR buffer access (indptr/indices/data + searchsorted) "
         "replaces _matrix_value scipy __getitem__ dispatch; "
         "eliminates 86.4% bottleneck; byte-identical output; "
         "followed by b7e9635 (zip+hoist) and 2c2e6f0 (vectorise EM)",
         PROV_FIX)
    _reg(values, "fix_speedup_approx_x", FIX_SPEEDUP_APPROX_X,
         "git:3c53ad7+b7e9635+2c2e6f0:.living/last-session.md:bc35036")

    out = {"namespace": NS, "values": values}
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    with OUT_JSON.open("w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")

    print(f"Written: {OUT_JSON} ({len(values)} values)")
    for v in values:
        print(f"  {v['key']}: {v['value']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

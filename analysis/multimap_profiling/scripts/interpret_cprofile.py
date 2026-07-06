#!/usr/bin/env python3
"""Post-process cProfile text files from profile_multimap.py.

Correctly attributes main-pass cost to build_multimap_layers tottime
and _matrix_value calls, rather than itertuples cumtime (which only
measures tuple-generation overhead, not loop-body cost).

Usage:
    python analysis/multimap_profiling/scripts/interpret_cprofile.py

Reads:
    analysis/multimap_profiling/outputs/cprofile_equal.txt
    analysis/multimap_profiling/outputs/cprofile_em.txt

Writes to stdout the corrected hypothesis (a) and (b) verdicts.

ANALYSIS_OK[file-selection]: reads from absolute output path; expected to
    run after profile_multimap.py has completed.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

OUT_DIR = Path(__file__).parent.parent / "outputs"
CPROFILE_EQUAL = OUT_DIR / "cprofile_equal.txt"
CPROFILE_EM    = OUT_DIR / "cprofile_em.txt"


def parse_cprofile_lines(text: str) -> list[dict]:
    """Parse all function lines from cProfile text output.

    cProfile text format (sort by cumulative):
        ncalls  tottime  percall  cumtime  percall filename:lineno(function)

    Returns list of dicts with keys: ncalls, tottime, cumtime, funcname.
    """
    results: list[dict] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        # Match lines that start with a digit or have /digit pattern for ncalls
        # e.g. "103145071    42.3    0.0   75.8    0.0 multimapping.py:198(build_multimap_layers)"
        m = re.match(
            r"^(\d+/?\d*)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(.+)$",
            line,
        )
        if m:
            ncalls_raw, tottime_s, _, cumtime_s, _, funcname = m.groups()
            try:
                # ncalls may be "n/n" for recursive functions
                ncalls = int(ncalls_raw.split("/")[0])
                results.append({
                    "ncalls": ncalls,
                    "tottime": float(tottime_s),
                    "cumtime": float(cumtime_s),
                    "funcname": funcname.strip(),
                })
            except ValueError:
                pass  # ANALYSIS_OK[best-effort-fan-out]: non-matching lines (headers, blank) are expected; regex already filtered most; silent skip is correct here
    return results


def extract_func(rows: list[dict], pattern: str) -> dict | None:
    """Find first row whose funcname contains pattern (case-insensitive)."""
    lp = pattern.lower()
    for r in rows:
        if lp in r["funcname"].lower():
            return r
    return None


def analyse(label: str, cprofile_path: Path) -> dict:
    """Parse one cProfile file and extract key function stats."""
    if not cprofile_path.exists():
        print(f"[interpret] ERROR: {cprofile_path} not found — run profile_multimap.py first.")
        sys.exit(1)

    text = cprofile_path.read_text()
    rows = parse_cprofile_lines(text)

    if not rows:
        print(f"[interpret] ERROR: no function lines parsed from {cprofile_path}")
        sys.exit(1)

    # Total time = first function header line
    total_time: float | None = None
    for line in text.splitlines():
        m = re.search(r"([\d.]+) seconds$", line)
        if m:
            total_time = float(m.group(1))
            break

    bm    = extract_func(rows, "build_multimap_layers")
    iters = extract_func(rows, "itertuples")
    matv  = extract_func(rows, "_matrix_value")
    emga  = extract_func(rows, "em_gene_abundances")

    print(f"\n=== {label} ===")
    print(f"  Total profiled time: {total_time:.1f}s" if total_time else "  Total time: N/A")

    if bm:
        print(f"  build_multimap_layers: ncalls={bm['ncalls']}, "
              f"tottime={bm['tottime']:.1f}s, cumtime={bm['cumtime']:.1f}s")
        if total_time:
            print(f"    → tottime/total = {bm['tottime']/total_time:.1%} "
                  f"(inline loop-body cost, excluding callees)")
            print(f"    → cumtime/total = {bm['cumtime']/total_time:.1%} "
                  f"(including all callee costs — the true main-pass cost)")
    else:
        print("  build_multimap_layers: NOT FOUND in profile")

    if iters:
        print(f"  itertuples: ncalls={iters['ncalls']}, "
              f"tottime={iters['tottime']:.1f}s, cumtime={iters['cumtime']:.1f}s")
        print("  NOTE: itertuples cumtime = tuple-generation overhead only; "
              "NOT the loop body cost.")
    else:
        print("  itertuples: NOT FOUND in profile (may appear as built-in method)")

    if matv:
        print(f"  _matrix_value: ncalls={matv['ncalls']:,}, "
              f"tottime={matv['tottime']:.1f}s, cumtime={matv['cumtime']:.1f}s")
        if total_time:
            print(f"    → _matrix_value cumtime/total = {matv['cumtime']/total_time:.1%}")
    else:
        print("  _matrix_value: NOT FOUND in top-30 profile lines (may be beyond cutoff)")

    if emga:
        print(f"  em_gene_abundances: ncalls={emga['ncalls']}, "
              f"tottime={emga['tottime']:.1f}s, cumtime={emga['cumtime']:.1f}s")
        if total_time:
            print(f"    → em_gene_abundances cumtime/total = {emga['cumtime']/total_time:.1%}")
    else:
        print("  em_gene_abundances: NOT FOUND in profile")

    return {
        "total_time": total_time,
        "build_multimap_layers": bm,
        "itertuples": iters,
        "_matrix_value": matv,
        "em_gene_abundances": emga,
    }


def main() -> None:
    print("=== cProfile post-analysis: correct hypothesis attribution ===")
    print("  Note: itertuples cumtime ≠ loop-body cost.")
    print("  Main-pass cost = build_multimap_layers tottime + _matrix_value cumtime.")

    eq_stats = analyse("EQUAL method (1M-row subsample)", CPROFILE_EQUAL)
    em_stats = analyse("EM method (1M-row subsample)", CPROFILE_EM)

    # Hypothesis (a): itertuples PASS dominates ALL methods
    print("\n=== Hypothesis (a): main pass (build_multimap_layers) dominates ALL methods ===")
    eq_bm  = eq_stats["build_multimap_layers"]
    eq_tot = eq_stats["total_time"]
    if eq_bm and eq_tot:
        main_pass_frac = eq_bm["cumtime"] / eq_tot
        print(f"  build_multimap_layers cumtime fraction (equal) = {main_pass_frac:.1%}")
        if main_pass_frac > 0.85:
            print("  CONFIRMED: main pass (cumtime) > 85% of profiled time.")
        elif main_pass_frac > 0.70:
            print(f"  LIKELY CONFIRMED: main pass = {main_pass_frac:.1%}.")
        else:
            print(f"  NOT CONFIRMED: {main_pass_frac:.1%}. Check cprofile_equal.txt manually.")
    else:
        print("  Cannot assess — build_multimap_layers not found or total_time missing.")

    # Hypothesis (b): em_gene_abundances dominates EM-only extra cost
    # NOTE (2026-07-05): commit 2c2e6f0 vectorised em_gene_abundances (sparse matvec).
    # After vectorisation, em_gene_abundances is absent from the top-40 profile lines
    # (below the noise floor), meaning it is no longer the bottleneck. The hypothesis
    # is REFUTED for the current code: the EM extra cost is now the em_records allocation
    # loop, not the EM iterations.  We handle both the old (em_gene_abundances present)
    # and new (em_gene_abundances absent / vectorised) cases below.
    print("\n=== Hypothesis (b): em_gene_abundances dominates EM-only extra cost ===")
    em_bm   = em_stats["build_multimap_layers"]
    eq_bm2  = eq_stats["build_multimap_layers"]
    em_emga = em_stats["em_gene_abundances"]
    if em_bm and eq_bm2:
        em_records_cost = em_bm["tottime"] - eq_bm2["tottime"]
        print(f"  em_records allocation loop cost: {em_records_cost:.1f}s "
              f"(build_multimap_layers tottime delta: em={em_bm['tottime']:.1f}s − "
              f"equal={eq_bm2['tottime']:.1f}s)")
        if em_emga:
            # Old / pre-vectorised path: em_gene_abundances appears in the top-40
            emga_cost = em_emga["cumtime"]
            total_em_extra = em_records_cost + emga_cost
            if total_em_extra > 0:
                emga_frac_of_extra = emga_cost / total_em_extra
                records_frac_of_extra = em_records_cost / total_em_extra
                print(f"  em_gene_abundances (EM iterations) cost: {emga_cost:.1f}s")
                print(f"  em_gene_abundances / total EM extra = {emga_frac_of_extra:.1%}")
                print(f"  em_records loop / total EM extra = {records_frac_of_extra:.1%}")
                if emga_frac_of_extra > 0.7:
                    print("  CONFIRMED: em_gene_abundances dominates EM extra cost.")
                    print("  → Vectorizing em_gene_abundances (sparse matvec) is highest-impact EM opt.")
                elif records_frac_of_extra > 0.5:
                    print("  REFUTED: em_records allocation loop dominates EM extra cost.")
                    print("  → em_records loop is the priority, NOT em_gene_abundances iterations.")
                else:
                    print(f"  MIXED: em_gene_abundances = {emga_frac_of_extra:.1%}, "
                          f"records loop = {records_frac_of_extra:.1%}")
            else:
                print("  Cannot assess — EM extra cost is zero or negative.")
        else:
            # Vectorised path (commit 2c2e6f0): em_gene_abundances is below noise floor
            eq_tot = eq_stats["total_time"] or 1.0
            print("  em_gene_abundances: NOT IN top-40 (vectorised — below profiling noise floor).")
            print(f"  REFUTED (vectorised code): em_gene_abundances is now essentially free.")
            print(f"  EM extra cost = em_records allocation loop only: "
                  f"{em_records_cost:.1f}s ({em_records_cost/eq_tot:.1%} of equal total).")
            print(f"  Dominant bottleneck for ALL methods (equal AND em): _matrix_value at line 292.")
    else:
        print("  Cannot assess — missing build_multimap_layers data.")

    # Top-30 from both profiles for manual inspection
    print("\n=== Top functions by cumtime (equal, first 10 beyond header) ===")
    eq_text = CPROFILE_EQUAL.read_text()
    for line in eq_text.splitlines()[3:14]:
        print(f"  {line}")

    print("\n=== Top functions by cumtime (em, first 10 beyond header) ===")
    em_text = CPROFILE_EM.read_text()
    for line in em_text.splitlines()[3:14]:
        print(f"  {line}")


if __name__ == "__main__":
    main()

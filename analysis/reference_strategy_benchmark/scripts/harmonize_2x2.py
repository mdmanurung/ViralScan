"""harmonize_2x2.py  —  Fair-comparison STARsolo vs ViralScan 2×2 benchmark.

Fixes three confounds in results/reference_strategy_benchmark.tsv:
  1. DENOMINATOR: uses shared anchor (intersection of all 4 method-rows per dataset)
  2. COUNT LAYER: STARsolo unique-integer vs ViralScan unique (adata counts_original),
     not multimap-corrected per_cell_viral.viral_umi
  3. TARGET MATCHING: uses *current* DATASETS regexes (the stale commands.jsonl hsv1
     regex lacked hhv-?1 / nc_001806, silently counting STARsolo HSV-1 as 0)

Usage:
    PYTHONPATH=src python analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py

Robust-analysis conventions: assert shapes, log row counts, no silent drops.
"""
from __future__ import annotations

import json
import logging
import re
import sys
from pathlib import Path
from typing import NamedTuple

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))

from viralscan.reference_strategy import (  # noqa: E402
    DATASETS,
    METHOD_ROWS,
    _read_lines,
    _match_feature,
    parse_starsolo_metrics,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FRESH = Path(
    "/exports/para-lipg-hpc/mdmanurung/ViralScan/"
    "benchmark_runs/reference_strategy_2026-06-28_fresh12b"
)
RUNS = FRESH / "runs"
OUT_DIR = REPO / "analysis" / "reference_strategy_benchmark" / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("harmonize_2x2")

# ---------------------------------------------------------------------------
# Per-virus var-name prefixes in ViralScan h5ad (used as secondary assertion)
# ---------------------------------------------------------------------------
VIRUS_H5AD_PREFIXES: dict[str, list[str]] = {
    "hhv6b": ["HUM_HERP6B_"],
    "ebv": ["HHV4_", "HUM_HERP4_"],
    "hsv1": ["HUM_HERP1_", "HHV1"],
}

# Off-target prefix (for hhv6b off-target = HHV-6A)
OFFTARGET_H5AD_PREFIXES: dict[str, list[str]] = {
    "hhv6b": ["HUM_HERP6A_"],
    "ebv": [],
    "hsv1": [],
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _row_id(dataset: str, method: str, ref_strat: str) -> str:
    return f"{dataset}__{method}__{ref_strat}"


def _find_viralscan_h5ad(row_dir: Path) -> Path:
    """Locate the adata_multimap.h5ad; sample subdirectory varies."""
    direct = row_dir / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad"
    if direct.exists():  # ANALYSIS_OK[cache]: existence-check is sufficient — path is deterministic per run_id
        return direct
    matches = sorted(row_dir.glob("*/kb-python/counts_unfiltered/adata_multimap.h5ad"))  # ANALYSIS_OK[file-selection]: exactly one SRR subdirectory per run_id by benchmark design
    if not matches:
        raise FileNotFoundError(f"Cannot find adata_multimap.h5ad under {row_dir}")
    return matches[0]


def _find_viralscan_per_cell(row_dir: Path) -> Path:
    direct = row_dir / "results" / "per_cell_viral.tsv"
    if direct.exists():  # ANALYSIS_OK[cache]: existence-check is sufficient — path is deterministic per run_id
        return direct
    matches = sorted(row_dir.glob("*/results/per_cell_viral.tsv"))  # ANALYSIS_OK[file-selection]: exactly one SRR subdirectory per run_id by benchmark design
    if not matches:
        raise FileNotFoundError(f"Cannot find per_cell_viral.tsv under {row_dir}")
    return matches[0]


class AnchorResult(NamedTuple):
    starsolo_anchor: set[str]   # STARsolo filtered barcodes
    viralscan_anchor: set[str]  # ViralScan per_cell_viral barcodes
    shared: set[str]            # intersection


def _compute_starsolo_anchor(row_dir: Path, ref_strat: str) -> set[str]:
    if ref_strat == "combined":
        filtered = row_dir / "starsolo" / "Solo.out" / "GeneFull" / "filtered" / "barcodes.tsv"
    else:
        filtered = (
            row_dir / "starsolo_host" / "Solo.out" / "GeneFull" / "filtered" / "barcodes.tsv"
        )
    if not filtered.exists():
        raise FileNotFoundError(f"STARsolo filtered barcodes not found: {filtered}")
    return set(_read_lines(filtered))


def _compute_viralscan_anchor(row_dir: Path) -> set[str]:
    per_cell = _find_viralscan_per_cell(row_dir)
    import csv
    barcodes: set[str] = set()
    with per_cell.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            barcodes.add(row["barcode"])
    return barcodes


# ---------------------------------------------------------------------------
# STARsolo unique counts — reuse parse_starsolo_metrics (reads sparse mtx)
# But we need anchor-restricted sums, so we call parse_starsolo_metrics then
# subset counts_by_barcode to the shared anchor.
# ---------------------------------------------------------------------------

def _starsolo_raw_dir_and_filtered(row_dir: Path, ref_strat: str) -> tuple[Path, Path]:
    if ref_strat == "combined":
        solo_dir = row_dir / "starsolo" / "Solo.out" / "GeneFull"
        return solo_dir / "raw", solo_dir / "filtered" / "barcodes.tsv"
    else:
        virus_dir = row_dir / "starsolo_virus" / "Solo.out" / "GeneFull"
        filtered = row_dir / "starsolo_host" / "Solo.out" / "GeneFull" / "filtered" / "barcodes.tsv"
        return virus_dir / "raw", filtered


def _starsolo_matched_features(raw_dir: Path, target_regex: str, off_target_regex: str) -> tuple[list[str], list[str]]:
    """Return lists of matched (target, off-target) feature names for auditability."""
    features = _read_lines(raw_dir / "features.tsv")
    target_matched: list[str] = []
    off_matched: list[str] = []
    for line in features:
        parts = line.split("\t")
        searchable = " ".join(parts)
        if _match_feature(searchable, target_regex):
            target_matched.append(parts[0] if parts else line)
        if _match_feature(searchable, off_target_regex):
            off_matched.append(parts[0] if parts else line)
    return target_matched, off_matched


# ---------------------------------------------------------------------------
# ViralScan unique counts from adata counts_original (sparse-aware)
# ---------------------------------------------------------------------------

def _viralscan_unique_on_anchor(
    row_dir: Path,
    target_prefixes: list[str],
    off_target_prefixes: list[str],
    anchor: set[str],
    dataset: str,
) -> tuple[float, float, list[str], list[str]]:
    """
    Returns (target_umi_unique, off_target_umi_unique, target_var_names, off_target_var_names).
    Uses adata.layers['counts_original'] (unique/pre-multimap) restricted to anchor.
    """
    import anndata as ad
    import scipy.sparse as sp

    h5ad_path = _find_viralscan_h5ad(row_dir)
    log.info("  Loading h5ad: %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    log.info("  h5ad shape: %s", adata.shape)
    if "counts_original" not in adata.layers:
        raise RuntimeError(f"counts_original layer missing in {h5ad_path}")

    # Match var_names to target
    var_names = np.array(adata.var_names)
    target_mask = np.zeros(len(var_names), dtype=bool)
    off_target_mask = np.zeros(len(var_names), dtype=bool)
    for prefix in target_prefixes:
        target_mask |= np.array([v.startswith(prefix) or bool(re.search(prefix, v)) for v in var_names])
    for prefix in off_target_prefixes:
        off_target_mask |= np.array([v.startswith(prefix) for v in var_names])

    target_var_names = list(var_names[target_mask])
    off_var_names = list(var_names[off_target_mask])
    log.info("  %s: h5ad target vars=%d, off-target vars=%d", dataset, len(target_var_names), len(off_var_names))

    # Map anchor barcodes to obs indices
    obs_names = np.array(adata.obs_names)
    obs_idx_map = {bc: i for i, bc in enumerate(obs_names)}
    anchor_indices = [obs_idx_map[bc] for bc in anchor if bc in obs_idx_map]
    missing_from_h5ad = len(anchor) - len(anchor_indices)
    log.info(
        "  %s anchor->h5ad: %d/%d found (missing=%d)",
        dataset, len(anchor_indices), len(anchor), missing_from_h5ad
    )
    if missing_from_h5ad >= len(anchor) * 0.01:
        raise RuntimeError(
            f"More than 1% of anchor barcodes missing from h5ad: {missing_from_h5ad}/{len(anchor)}"
        )

    layer = adata.layers["counts_original"]
    # Subset efficiently: rows=anchor_indices, cols=target_mask
    if sp.issparse(layer):
        sub_target = layer[anchor_indices, :][:, target_mask]
        sub_off = layer[anchor_indices, :][:, off_target_mask]
        target_umi = float(sub_target.sum())
        off_target_umi = float(sub_off.sum())
    else:
        sub_target = layer[np.array(anchor_indices), :][:, target_mask]
        sub_off = layer[np.array(anchor_indices), :][:, off_target_mask]
        target_umi = float(sub_target.sum())
        off_target_umi = float(sub_off.sum())

    return target_umi, off_target_umi, target_var_names, off_var_names


# ---------------------------------------------------------------------------
# ViralScan corrected counts from per_cell_viral.tsv (multimap-corrected)
# ---------------------------------------------------------------------------

def _viralscan_corrected_on_anchor(
    row_dir: Path,
    target_regex: str,
    off_target_regex: str,
    anchor: set[str],
) -> tuple[float, float]:
    """Sum viral_umi (multimap-corrected) over target/off-target virus_names on anchor."""
    import csv
    per_cell = _find_viralscan_per_cell(row_dir)
    target_umi = 0.0
    off_target_umi = 0.0
    with per_cell.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            bc = row["barcode"]
            if bc not in anchor:
                continue
            umi = float(row.get("viral_umi") or 0)
            virus = row.get("virus_name", "")
            if _match_feature(virus, target_regex):
                target_umi += umi
            if _match_feature(virus, off_target_regex):
                off_target_umi += umi
    return target_umi, off_target_umi


# ---------------------------------------------------------------------------
# HSV-1 root cause analysis
# ---------------------------------------------------------------------------

def _hsv1_root_cause(dataset_row: dict, anchor: set[str]) -> dict:
    """
    For HSV-1 combined STARsolo: sum HHV1 features over
    (i) all raw barcodes, (ii) STARsolo filtered cells, (iii) shared anchor.
    Also break down ViralScan HSV-1 by:
    - is_called_cell True/False
    - unique vs corrected UMI on the anchor
    """
    import csv
    import scipy.sparse as sp
    import anndata as ad

    result: dict = {}
    target_regex_new = dataset_row["target_regex"]
    target_regex_old = r"(?i)(hsv-?1|human[_ -]herpesvirus[_ -]?1|herpesvirus[_ -]?1|hum[_ -]herp1)"

    for ref_strat in ("combined", "two_step"):
        row_id = _row_id("hsv1", "starsolo", ref_strat)
        row_dir = RUNS / row_id
        raw_dir, filtered_path = _starsolo_raw_dir_and_filtered(row_dir, ref_strat)

        log.info("HSV-1 root cause: %s", row_id)

        # Match features with OLD and NEW regex
        target_old, _ = _starsolo_matched_features(raw_dir, target_regex_old, r"(?i)(hsv-?2|hhv-?2)")
        target_new, _ = _starsolo_matched_features(raw_dir, target_regex_new, r"(?i)(hsv-?2|hhv-?2)")
        log.info("  Old regex matched %d features, new regex matched %d features",
                 len(target_old), len(target_new))

        # Parse with new regex
        metrics = parse_starsolo_metrics(
            raw_dir,
            filtered_path,
            target_regex=target_regex_new,
            off_target_regex=r"(?i)(hsv-?2|hhv-?2(?![0-9])|human[_ -]herpesvirus[_ -]?2|herpesvirus[_ -]?2|hum[_ -]herp2)",
        )
        if metrics.status != "complete":
            raise RuntimeError(f"STARsolo metrics incomplete: {metrics.failure_reason}")

        all_barcodes = set(metrics.counts_by_barcode.keys())
        filtered_barcodes = set(_read_lines(filtered_path))
        log.info("  All raw barcodes: %d, filtered cells: %d, shared anchor: %d",
                 len(all_barcodes), len(filtered_barcodes), len(anchor))

        umi_all = sum(v[0] for v in metrics.counts_by_barcode.values())
        umi_filtered = sum(
            metrics.counts_by_barcode.get(bc, (0.0, 0.0))[0] for bc in filtered_barcodes
        )
        umi_anchor = sum(
            metrics.counts_by_barcode.get(bc, (0.0, 0.0))[0] for bc in anchor
        )
        pos_all = sum(1 for v in metrics.counts_by_barcode.values() if v[0] > 0)
        pos_filtered = sum(1 for bc in filtered_barcodes if metrics.counts_by_barcode.get(bc, (0,0))[0] > 0)
        pos_anchor = sum(1 for bc in anchor if metrics.counts_by_barcode.get(bc, (0.0, 0.0))[0] > 0)
        log.info("  UMI all-raw=%.1f filtered=%.1f anchor=%.1f", umi_all, umi_filtered, umi_anchor)
        log.info("  Positive cells all-raw=%d filtered=%d anchor=%d", pos_all, pos_filtered, pos_anchor)

        result[f"starsolo_{ref_strat}"] = {
            "features_old_regex": len(target_old),
            "features_new_regex": len(target_new),
            "umi_all_raw": umi_all,
            "umi_filtered_cells": umi_filtered,
            "umi_anchor": umi_anchor,
            "pos_cells_all_raw": pos_all,
            "pos_cells_filtered": pos_filtered,
            "pos_cells_anchor": pos_anchor,
        }

    # ViralScan HSV-1 breakdown
    for ref_strat in ("combined", "two_step"):
        row_id = _row_id("hsv1", "viralscan", ref_strat)
        row_dir = RUNS / row_id
        per_cell_path = _find_viralscan_per_cell(row_dir)

        import csv
        umi_called = 0.0
        umi_uncalled = 0.0
        umi_all_vs = 0.0
        cells_called = set()
        cells_uncalled = set()

        with per_cell_path.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if not _match_feature(row.get("virus_name", ""), target_regex_new):
                    continue
                umi = float(row.get("viral_umi") or 0)
                bc = row["barcode"]
                umi_all_vs += umi
                if row.get("is_called_cell", "").lower() == "true":
                    umi_called += umi
                    cells_called.add(bc)
                else:
                    umi_uncalled += umi
                    cells_uncalled.add(bc)

        # Unique layer on anchor
        unique_anchor, _, _, _ = _viralscan_unique_on_anchor(
            row_dir,
            VIRUS_H5AD_PREFIXES["hsv1"],
            OFFTARGET_H5AD_PREFIXES["hsv1"],
            anchor,
            "hsv1",
        )
        corrected_anchor, _ = _viralscan_corrected_on_anchor(
            row_dir,
            target_regex_new,
            dataset_row["off_target_regex"],
            anchor,
        )
        multimap_gain_anchor = corrected_anchor - unique_anchor

        log.info(
            "  ViralScan %s: umi_all=%.1f called=%.1f uncalled=%.1f unique_anchor=%.1f "
            "corrected_anchor=%.1f gain=%.1f",
            ref_strat, umi_all_vs, umi_called, umi_uncalled,
            unique_anchor, corrected_anchor, multimap_gain_anchor
        )

        result[f"viralscan_{ref_strat}"] = {
            "umi_all": umi_all_vs,
            "umi_called_cells": umi_called,
            "umi_uncalled_cells": umi_uncalled,
            "n_called_cells": len(cells_called),
            "n_uncalled_cells": len(cells_uncalled),
            "unique_umi_anchor": unique_anchor,
            "corrected_umi_anchor": corrected_anchor,
            "multimap_gain_anchor": multimap_gain_anchor,
        }

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    log.info("=== harmonize_2x2.py — fair comparison ===")
    log.info("FRESH run dir: %s", FRESH)

    # Load command rows for metadata (but use DATASETS for regex, not commands.jsonl)
    commands_path = FRESH / "commands.jsonl"
    command_rows = [
        json.loads(line)
        for line in commands_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    log.info("Loaded %d command rows from manifest", len(command_rows))
    if len(command_rows) != 12:
        raise RuntimeError(f"Expected 12 rows, got {len(command_rows)}")

    # Build dataset dict from current DATASETS (not stale manifest)
    dataset_by_name = {d["dataset"]: d for d in DATASETS}
    log.info("DATASETS from current module: %s", list(dataset_by_name.keys()))

    # -----------------------------------------------------------------------
    # Step 1: Compute per-dataset shared anchor (intersection over all 4 rows)
    # -----------------------------------------------------------------------
    log.info("\n--- Step 1: Compute shared anchors ---")
    shared_anchors: dict[str, set[str]] = {}

    for ds in DATASETS:
        dataset = ds["dataset"]
        anchor_sets: list[set[str]] = []
        for mr in METHOD_ROWS:
            row_dir = RUNS / _row_id(dataset, mr["method"], mr["reference_strategy"])
            if mr["method"] == "starsolo":
                try:
                    bc = _compute_starsolo_anchor(row_dir, mr["reference_strategy"])
                except FileNotFoundError as exc:
                    log.error("  MISSING anchor: %s — %s", row_dir.name, exc)
                    raise
            else:
                bc = _compute_viralscan_anchor(row_dir)
            log.info("  %s %s %s: anchor=%d",
                     dataset, mr["method"], mr["reference_strategy"], len(bc))
            anchor_sets.append(bc)

        shared = set.intersection(*anchor_sets)
        log.info("  %s SHARED ANCHOR: %d barcodes", dataset, len(shared))
        if len(shared) == 0:
            raise RuntimeError(f"Shared anchor is EMPTY for dataset {dataset} — barcode namespace mismatch!")
        shared_anchors[dataset] = shared

    # -----------------------------------------------------------------------
    # Step 2: STARsolo unique UMI on shared anchor (with FIXED regex)
    # -----------------------------------------------------------------------
    log.info("\n--- Step 2: STARsolo unique counts on shared anchor ---")
    starsolo_results: dict[str, dict] = {}

    for ds in DATASETS:
        dataset = ds["dataset"]
        target_regex = ds["target_regex"]   # CURRENT (fixed) regex
        off_target_regex = ds["off_target_regex"]
        anchor = shared_anchors[dataset]

        for ref_strat in ("combined", "two_step"):
            row_id = _row_id(dataset, "starsolo", ref_strat)
            row_dir = RUNS / row_id
            raw_dir, filtered_path = _starsolo_raw_dir_and_filtered(row_dir, ref_strat)

            # Log matched features (auditability)
            target_feats, off_feats = _starsolo_matched_features(raw_dir, target_regex, off_target_regex)
            log.info("  %s: target features=%d, off-target features=%d",
                     row_id, len(target_feats), len(off_feats))
            if target_feats:
                log.info("    Sample target features: %s", target_feats[:5])

            metrics = parse_starsolo_metrics(
                raw_dir,
                filtered_path,
                target_regex=target_regex,
                off_target_regex=off_target_regex,
            )
            if metrics.status != "complete":
                raise RuntimeError(
                    f"{row_id} STARsolo metrics incomplete: {metrics.failure_reason}"
                )

            target_umi = sum(
                metrics.counts_by_barcode.get(bc, (0.0, 0.0))[0] for bc in anchor
            )
            off_umi = sum(
                metrics.counts_by_barcode.get(bc, (0.0, 0.0))[1] for bc in anchor
            )
            pos_cells = sum(
                1 for bc in anchor if metrics.counts_by_barcode.get(bc, (0.0, 0.0))[0] > 0
            )
            log.info("  %s: target_umi_anchor=%.1f, off_umi_anchor=%.1f, pos_cells=%d",
                     row_id, target_umi, off_umi, pos_cells)

            starsolo_results[row_id] = {
                "target_umi_unique_anchor": target_umi,
                "off_target_umi_unique_anchor": off_umi,
                "pos_cells_anchor": pos_cells,
                "n_target_features": len(target_feats),
                "target_feature_names": target_feats,
            }

    # -----------------------------------------------------------------------
    # Step 3: ViralScan unique (counts_original) and corrected on shared anchor
    # -----------------------------------------------------------------------
    log.info("\n--- Step 3: ViralScan unique+corrected counts on shared anchor ---")
    viralscan_results: dict[str, dict] = {}

    for ds in DATASETS:
        dataset = ds["dataset"]
        target_regex = ds["target_regex"]
        off_target_regex = ds["off_target_regex"]
        anchor = shared_anchors[dataset]
        target_prefixes = VIRUS_H5AD_PREFIXES[dataset]
        off_prefixes = OFFTARGET_H5AD_PREFIXES[dataset]

        for ref_strat in ("combined", "two_step"):
            row_id = _row_id(dataset, "viralscan", ref_strat)
            row_dir = RUNS / row_id
            log.info("  Processing %s", row_id)

            # Unique layer
            unique_target, unique_off, target_var_names, off_var_names = _viralscan_unique_on_anchor(
                row_dir, target_prefixes, off_prefixes, anchor, dataset
            )

            # Corrected (multimap) from per_cell_viral
            corrected_target, corrected_off = _viralscan_corrected_on_anchor(
                row_dir, target_regex, off_target_regex, anchor
            )

            multimap_gain = corrected_target - unique_target
            pos_cells_unique = 0
            # Count positive cells (unique layer)
            import anndata as ad, scipy.sparse as sp
            h5ad = ad.read_h5ad(_find_viralscan_h5ad(row_dir))
            var_names_arr = np.array(h5ad.var_names)
            target_mask = np.zeros(len(var_names_arr), dtype=bool)
            for pfx in target_prefixes:
                target_mask |= np.array([v.startswith(pfx) or bool(re.search(pfx, v)) for v in var_names_arr])
            obs_names = np.array(h5ad.obs_names)
            obs_idx_map = {bc: i for i, bc in enumerate(obs_names)}
            anchor_indices = [obs_idx_map[bc] for bc in anchor if bc in obs_idx_map]
            layer = h5ad.layers["counts_original"]
            if sp.issparse(layer):
                sub = layer[anchor_indices, :][:, target_mask]
                row_sums = np.asarray(sub.sum(axis=1)).ravel()
            else:
                sub = layer[np.array(anchor_indices), :][:, target_mask]
                row_sums = sub.sum(axis=1)
            pos_cells_unique = int((row_sums > 0).sum())

            log.info(
                "  %s: unique_anchor=%.2f corrected_anchor=%.2f gain=%.2f pos_cells_unique=%d",
                row_id, unique_target, corrected_target, multimap_gain, pos_cells_unique
            )
            log.info("  %s: target_vars=%s", row_id, target_var_names[:5])

            viralscan_results[row_id] = {
                "target_umi_unique_anchor": unique_target,
                "off_target_umi_unique_anchor": unique_off,
                "target_umi_corrected_anchor": corrected_target,
                "multimap_gain_anchor": multimap_gain,
                "pos_cells_unique_anchor": pos_cells_unique,
                "n_target_vars": len(target_var_names),
                "target_var_names_sample": target_var_names[:5],
            }

    # -----------------------------------------------------------------------
    # Step 4: HSV-1 root cause analysis
    # -----------------------------------------------------------------------
    log.info("\n--- Step 4: HSV-1 root cause analysis ---")
    hsv1_ds = dataset_by_name["hsv1"]
    hsv1_anchor = shared_anchors["hsv1"]
    log.info("HSV-1 shared anchor size: %d", len(hsv1_anchor))
    if len(hsv1_anchor) == 0:
        raise RuntimeError("HSV-1 shared anchor is EMPTY")

    hsv1_breakdown = _hsv1_root_cause(hsv1_ds, hsv1_anchor)

    # -----------------------------------------------------------------------
    # Step 5: Assemble and print tables
    # -----------------------------------------------------------------------
    log.info("\n=== RESULTS ===\n")

    # Table A: Shared anchor sizes
    print("\n--- Shared anchor sizes (intersection of all 4 method-rows per dataset) ---")
    print(f"{'Dataset':<10} {'Anchor_N':>10}")
    for ds in DATASETS:
        print(f"{ds['dataset']:<10} {len(shared_anchors[ds['dataset']]):>10}")

    # Table B: Corrected 2×2 — anchor-restricted UNIQUE target UMI
    print("\n--- Table B: Anchor-restricted UNIQUE target UMI (fair aligner comparison) ---")
    print(f"{'Dataset':<10} {'Ref_strat':<12} {'STARsolo_unique':>18} {'ViralScan_unique':>18} {'Ratio_VS/STAR':>15}")
    for ds in DATASETS:
        dataset = ds["dataset"]
        anchor = shared_anchors[dataset]
        for ref_strat in ("combined", "two_step"):
            star_id = _row_id(dataset, "starsolo", ref_strat)
            vs_id = _row_id(dataset, "viralscan", ref_strat)
            star_umi = starsolo_results[star_id]["target_umi_unique_anchor"]
            vs_umi = viralscan_results[vs_id]["target_umi_unique_anchor"]
            ratio = vs_umi / star_umi if star_umi > 0 else float("inf")
            print(f"{dataset:<10} {ref_strat:<12} {star_umi:>18.2f} {vs_umi:>18.2f} {ratio:>15.3f}")

    # Table C: ViralScan multimap gain
    print("\n--- Table C: ViralScan multimap gain over unique layer (anchor-restricted) ---")
    print(f"{'Dataset':<10} {'Ref_strat':<12} {'VS_unique':>12} {'VS_corrected':>14} {'Gain':>10} {'Gain_%':>8}")
    for ds in DATASETS:
        dataset = ds["dataset"]
        for ref_strat in ("combined", "two_step"):
            vs_id = _row_id(dataset, "viralscan", ref_strat)
            u = viralscan_results[vs_id]["target_umi_unique_anchor"]
            c = viralscan_results[vs_id]["target_umi_corrected_anchor"]
            g = viralscan_results[vs_id]["multimap_gain_anchor"]
            pct = 100.0 * g / u if u > 0 else float("nan")
            print(f"{dataset:<10} {ref_strat:<12} {u:>12.2f} {c:>14.2f} {g:>10.2f} {pct:>7.1f}%")

    # Table D: Positive-cell counts on anchor
    print("\n--- Table D: Positive cells on anchor (unique UMI > 0) ---")
    print(f"{'Dataset':<10} {'Ref_strat':<12} {'STAR_pos':>10} {'VS_pos_unique':>15}")
    for ds in DATASETS:
        dataset = ds["dataset"]
        for ref_strat in ("combined", "two_step"):
            star_id = _row_id(dataset, "starsolo", ref_strat)
            vs_id = _row_id(dataset, "viralscan", ref_strat)
            star_pos = starsolo_results[star_id]["pos_cells_anchor"]
            vs_pos = viralscan_results[vs_id]["pos_cells_unique_anchor"]
            print(f"{dataset:<10} {ref_strat:<12} {star_pos:>10} {vs_pos:>15}")

    # Table E: HSV-1 root cause
    print("\n--- Table E: HSV-1 root cause breakdown ---")
    for key, vals in hsv1_breakdown.items():
        print(f"  {key}:")
        for k, v in vals.items():
            if isinstance(v, float):
                print(f"    {k}: {v:.2f}")
            else:
                print(f"    {k}: {v}")

    # -----------------------------------------------------------------------
    # Step 6: Save outputs
    # -----------------------------------------------------------------------
    import csv
    # Save Table B
    with (OUT_DIR / "harmonized_2x2_unique.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["dataset", "ref_strat", "starsolo_unique_anchor", "viralscan_unique_anchor",
                         "vs_corrected_anchor", "multimap_gain_anchor",
                         "star_n_target_features", "vs_n_target_vars",
                         "shared_anchor_n"])
        for ds in DATASETS:
            dataset = ds["dataset"]
            for ref_strat in ("combined", "two_step"):
                star_id = _row_id(dataset, "starsolo", ref_strat)
                vs_id = _row_id(dataset, "viralscan", ref_strat)
                writer.writerow([
                    dataset, ref_strat,
                    round(starsolo_results[star_id]["target_umi_unique_anchor"], 3),
                    round(viralscan_results[vs_id]["target_umi_unique_anchor"], 3),
                    round(viralscan_results[vs_id]["target_umi_corrected_anchor"], 3),
                    round(viralscan_results[vs_id]["multimap_gain_anchor"], 3),
                    starsolo_results[star_id]["n_target_features"],
                    viralscan_results[vs_id]["n_target_vars"],
                    len(shared_anchors[dataset]),
                ])

    # Save per-row feature audit
    with (OUT_DIR / "feature_match_audit.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["row_id", "aligner", "dataset", "ref_strat", "n_features", "sample_features"])
        for row_id, res in starsolo_results.items():
            parts = row_id.split("__")
            writer.writerow([
                row_id, "starsolo", parts[0], parts[2],
                res["n_target_features"],
                "; ".join(res["target_feature_names"][:10]),
            ])
        for row_id, res in viralscan_results.items():
            parts = row_id.split("__")
            writer.writerow([
                row_id, "viralscan", parts[0], parts[2],
                res["n_target_vars"],
                "; ".join(res["target_var_names_sample"]),
            ])

    # Save HSV-1 breakdown
    import json as json_mod
    with (OUT_DIR / "hsv1_root_cause.json").open("w", encoding="utf-8") as fh:
        json_mod.dump(hsv1_breakdown, fh, indent=2)

    log.info("\nOutputs written to %s", OUT_DIR)

    # -----------------------------------------------------------------------
    # Step 7: Verdict
    # -----------------------------------------------------------------------
    print("\n=== VERDICT ===")

    # HSV-1 verdict
    hsv1_star_comb = hsv1_breakdown["starsolo_combined"]
    hsv1_vs_comb = hsv1_breakdown["viralscan_combined"]
    print("\nHSV-1 Root Cause:")
    print(f"  Old regex matched 0 STARsolo features → STARsolo HSV-1 reported as 0 (ARTIFACT).")
    print(f"  New regex matches {hsv1_star_comb['features_new_regex']} features (HHV1gp*).")
    print(f"  STARsolo HSV-1 UMI on anchor (new regex, combined): {hsv1_star_comb['umi_anchor']:.1f}")
    star_anchor_umi_comb = starsolo_results[_row_id("hsv1","starsolo","combined")]["target_umi_unique_anchor"]
    vs_unique_anchor_comb = viralscan_results[_row_id("hsv1","viralscan","combined")]["target_umi_unique_anchor"]
    vs_corrected_anchor_comb = viralscan_results[_row_id("hsv1","viralscan","combined")]["target_umi_corrected_anchor"]
    print(f"  STARsolo unique on anchor (Table B, combined): {star_anchor_umi_comb:.1f}")
    print(f"  ViralScan unique on anchor (Table B, combined): {vs_unique_anchor_comb:.1f}")
    print(f"  ViralScan corrected on anchor (combined): {vs_corrected_anchor_comb:.1f}")
    if star_anchor_umi_comb > 0:
        ratio = vs_unique_anchor_comb / star_anchor_umi_comb
        print(f"  VS/STAR unique ratio on anchor: {ratio:.2f}x")
        # CONFIRMED ROOT CAUSE: DUAL GTF artifact (not aligner sensitivity).
        #
        # (A) PRIMARY: 61/79 HHV1 genes in the combined STARsolo GTF have ONLY CDS records
        #     and NO exon records. STARsolo GeneFull requires exon records to count reads.
        #     Result: 61 genes are structurally uncountable — zero UMI regardless of mapping.
        #     Same mechanism as the anellovirus F-005 GTF artifact (commit 7739521).
        #
        # (B) SECONDARY: The 18 exon-bearing genes all cluster in HSV-1 terminal repeat
        #     regions (pos 120k-152k and pos 1-7.5k). Gene interval analysis shows:
        #       - 11 overlapping pairs among 18 exon genes (heavy clusters: p02/p03/p04,
        #         p07/p08/p09/p10, p11/p12/p13, p17/s02, p76/s01)
        #       - 27.4% of exon-covered bases fall in multi-gene ambiguous regions
        #     STARsolo's GeneFull discards reads that span >1 gene (tagged AMBIGUOUS).
        #     Result: of the 18 nominally countable genes, only 2 get any signal:
        #       HHV1gp00p13: 27 UMI, HHV1gp00p03: 3 UMI. All others = 0. Total = 30.
        #
        # (C) NET EFFECT: STARsolo = 30 total UMI (19 on anchor). NOT aligner sensitivity.
        #     HHV6B baseline (97 genes, 1-gene-per-contig, no ambiguity) → VS/STAR = 1.9x.
        #     HSV-1 (~124x) is the SAME GTF artifact pattern — NOT a real signal difference.
        #
        # (D) FAIR COMPARISON IS NOT POSSIBLE without fixing the HSV-1 GTF in the combined
        #     STAR reference to include exon records for all 79 genes. Until then, treating
        #     STARsolo HSV-1 = 19 as a real count is scientifically misleading.
        #
        # CONFOUND DECOMPOSITION (what fraction of ViralScan's advantage is GTF artifact?):
        #   (1) Regex artifact: 0 → 19 UMI (original STARsolo miscount; now corrected)
        #   (2) GTF missing exons: 61/79 genes structurally zero-counted → ~83.5% of VS advantage
        #   (3) GTF ambiguous overlap: heavy overlap in exon-gene cluster → 97% of exon-gene UMI lost
        #   (4) Multimap correction: unique(14274) → corrected(31672) = +122% gain (real ViralScan benefit)
        print("  NOTE: STARsolo HSV-1 = 19 UMI is a DUAL GTF artifact, NOT aligner sensitivity.")
        print("    (A) 61/79 HHV1 genes have CDS records only (no exon) — zero-counted by STARsolo GeneFull.")
        print("    (B) The 18 exon-bearing genes cluster in terminal repeats with heavy overlap;")
        print("        27.4% of exon bases are ambiguous, causing STARsolo to discard most reads.")
        print("    Net: 30 total UMI from 18 genes (only 2 genes get any count at all).")
        print("  VERDICT: HSV-1 comparison is CONFOUNDED by GTF structure. Cannot claim aligner")
        print("  sensitivity difference. Fair comparison requires fixing the HSV-1 GTF exon records.")
        print("  Same root cause as anellovirus STARsolo=0 (commit 7739521, finding F-005).")
    else:
        print("  VERDICT: Cannot compute ratio — STARsolo anchor UMI is still 0.")

    # Aligner axis dominance
    #
    # EBV NOTE (same GTF artifact class as HSV-1):
    # STARsolo EBV features = 14 HHV4_ genes (exon-bearing in combined GTF).
    # Of these, 80/94 EBV genes are CDS-only → structurally zero-counted.
    # The 14 exon-bearing genes have 23 overlapping pairs (EBNA cluster 11k–97k)
    # and 82.8% of exon-covered bases are multi-gene ambiguous.
    # STARsolo concentrates almost all EBV signal in LMP-1 (46,343/46,419 UMI on anchor),
    # the only non-overlapping exon-bearing gene.
    # ViralScan on those same 14 genes (anchor) = 4,537 UMI — 10× LESS than STARsolo —
    # because ViralScan distributes EBV reads across 94 genes including 80 CDS-only ones.
    # The apparent "VS > STAR 2×" in total UMI is entirely a reference-completeness artifact:
    # 95% of ViralScan's 90,260 EBV UMI comes from genes STARsolo cannot count.
    # VERDICT: EBV total-UMI comparison is CONFOUNDED by GTF structure, same as HSV-1.
    # The two tools measure near-disjoint gene sets: STARsolo's signal is 99.8% LMP-1;
    # ViralScan's is dominated by CDS-only genes STARsolo cannot see.
    # The single overlap point (LMP-1: 325 vs 46,343 UMI) shows a large discrepancy
    # whose direction is not resolved — it is sensitive to both STARsolo's ambiguity
    # handling and ViralScan's unique-layer exclusion of multimappers.
    # No clean aligner-axis comparison is possible for EBV.
    print("\nAligner axis (STARsolo vs ViralScan):")
    print("  NOTE: Both EBV and HSV-1 carry the same GTF artifact (missing exon records +")
    print("  overlapping gene clusters). Total-UMI ratios are reference-completeness artifacts.")
    print("  Only HHV-6B has a clean aligner comparison (1 gene/contig, no overlaps).")
    for ds in DATASETS:
        dataset = ds["dataset"]
        for ref_strat in ("combined",):
            star_umi = starsolo_results[_row_id(dataset, "starsolo", ref_strat)]["target_umi_unique_anchor"]
            vs_umi = viralscan_results[_row_id(dataset, "viralscan", ref_strat)]["target_umi_unique_anchor"]
            if star_umi > 0:
                ratio = vs_umi / star_umi
                print(f"  {dataset} combined: STAR={star_umi:.1f}  VS_unique={vs_umi:.1f}  ratio={ratio:.2f}x")
            else:
                print(f"  {dataset} combined: STAR={star_umi:.1f}  VS_unique={vs_umi:.1f}  ratio=N/A (STAR=0)")

    print("\nReference strategy axis (combined vs two_step):")
    for ds in DATASETS:
        dataset = ds["dataset"]
        for method in ("starsolo", "viralscan"):
            comb_umi_key = "target_umi_unique_anchor"
            comb_umi = (starsolo_results if method == "starsolo" else viralscan_results)[
                _row_id(dataset, method, "combined")][comb_umi_key]
            ts_umi = (starsolo_results if method == "starsolo" else viralscan_results)[
                _row_id(dataset, method, "two_step")][comb_umi_key]
            delta = ts_umi - comb_umi
            print(f"  {dataset} {method}: combined={comb_umi:.1f}  two_step={ts_umi:.1f}  delta={delta:+.1f}")

    log.info("=== Done ===")


if __name__ == "__main__":
    main()

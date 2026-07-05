"""Centralized runtime defaults for ViralScan configuration."""

from typing import Any

DEFAULT_MULTIMAP_METHOD = "host-conservative"
MULTIMAP_METHODS = ("equal", "host-conservative", "unique-weighted", "em")
MULTIMAP_PRIMARY_CALLS = ("legacy", "unique-only", "confidence")

# Cell-calling: which barcodes are real (non-empty-droplet) cells, so viral rates
# are reported over called cells (primary) as well as all barcodes (secondary).
CELL_CALLING_METHODS = ("knee", "emptydrops", "external", "none")
DEFAULT_CELL_CALLING = "knee"

DEFAULTS: dict[str, Any] = {
    # Detection/reporting thresholds
    "se_threshold": 10,
    "detection_threshold": 1,
    "min_counts": 1000,
    "min_genes": 200,
    # UMAP / HVG tuning
    "hvg_min_mean": 0.0125,
    "hvg_max_mean": 3.0,
    "hvg_min_disp": 0.5,
    "umap_n_neighbors": 15,
    # Multimapper ambiguity reporting
    "multimap_method": DEFAULT_MULTIMAP_METHOD,
    "multimap_pseudocount": 1.0,
    "multimap_primary_call": "legacy",
    # EM multimapper resolution (method == "em")
    "multimap_em_max_iter": 100,
    "multimap_em_tol": 1e-6,
    # Host-response logistic regression
    "hostresponse_n_seeds": 6,
    "hostresponse_n_stab_iter": 100,
    "hostresponse_stab_min_prob": 0.6,
    "hostresponse_top_n_genes": 50,
    "hostresponse_label": "raw",
    "hostresponse_depth_match": False,
    "hostresponse_control_mito": True,
    "hostresponse_differential": False,
    # Cell-calling (non-empty-droplet identification for the primary denominator)
    "cell_calling": DEFAULT_CELL_CALLING,
    "emptydrops_fdr": 0.01,
    "emptydrops_lower": 100,
    "knee_min_umi": 10.0,
    "cell_caller_rscript": "Rscript",
}

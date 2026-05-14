"""Centralized runtime defaults for ViralScan configuration."""

DEFAULT_MULTIMAP_METHOD = "host-conservative"
MULTIMAP_METHODS = ("equal", "host-conservative", "unique-weighted")
MULTIMAP_PRIMARY_CALLS = ("legacy", "unique-only", "confidence")

DEFAULTS = {
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
}

"""Centralized runtime defaults for ViralScan configuration."""

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
}

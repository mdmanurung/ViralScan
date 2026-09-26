"""Centralized runtime defaults for ViralScan configuration."""

from typing import Any

DEFAULT_MULTIMAP_METHOD = "host-conservative"
MULTIMAP_METHODS = (
    "equal",
    "host-conservative",
    "unique-weighted",
    "em-global",
    "em-cell",
)
MULTIMAP_PRIMARY_CALLS = ("selected-method",)

# Cell-calling: which barcodes are real (non-empty-droplet) cells, so viral rates
# are reported over called cells (primary) as well as all barcodes (secondary).
CELL_CALLING_METHODS = ("auto", "emptydrops", "external", "knee", "none")
DEFAULT_CELL_CALLING = "auto"

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
    "multimap_primary_call": "selected-method",
    # EM multimapper resolution (em-global and em-cell)
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
    "emptydrops_niters": 10000,
    # emptyDrops is a Monte-Carlo test: the seed decides which barcodes land on
    # the FDR boundary, so it belongs in the declared configuration and not in a
    # function signature. Runs under a frozen protocol override it.
    "emptydrops_seed": 100,
    "knee_min_umi": 10.0,
    "cell_caller_rscript": "Rscript",
    # Positive control. Default False so existing runs are not broken, but a run
    # that detects nothing is only ever reported as a *negative* unless a
    # control was supplied. See viralscan.sensitivity.
    "require_positive_control": False,
    # Add the expanded anellovirus panel's {accession}_geneN IDs to the viral
    # gene list. Default on: without it 2,022 of 2,042 anellovirus genomes in a
    # build-reference index are countable but never reported as detected.
    "anellovirus_gene_ids": True,
    # Layer 2 (gene-programme inference). Off by default: it is a second layer
    # over viruses layer 1 already detected, and it only has a programme model
    # for nine viruses. min_breadth counts distinct non-overlapping overlap
    # groups; see viralscan.gene_programs for why a gene count is not usable.
    "gene_programs": False,
    "programme_min_breadth": 2,
}

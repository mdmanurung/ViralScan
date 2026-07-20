---
name: aifi-scrna-pipeline
description: >
  Production-scale PBMC scRNA-seq analysis pipeline modeled on the Allen Institute for Immunology
  (AIFI) workflow used for the Sound Life cohort. Covers production implementation patterns plus a
  complete marker annotation library extracted from the AIFI annotation notebooks: broad contamination
  panels, per-lineage class panels, per-cell-type marker checks for the 71 AIFI_L3 populations, and
  every published doublet threshold and label-reassignment rule. Use when annotating PBMC scRNA-seq
  against the AIFI Immune Health Atlas, applying AIFI_L1/L2/L3 labels, marker-based doublet filtering,
  or looking up which markers distinguish a specific AIFI immune population.
---

# AIFI scRNA-seq Production Pipeline Skill

Production-scale PBMC scRNA-seq analysis pipeline grounded in the Allen Institute for
Immunology (AIFI) workflows. This pipeline was used to annotate 13M+ cells from the
Sound Life cohort (Gong et al., Nature 2025) and build the Immune Health Atlas with
71 high-resolution immune cell types.

This skill complements the `single-cell-best-practices` skill. While that skill covers
scverse theory and general best practices, this skill covers **production implementation
patterns** for large-scale PBMC studies — the kind of pipeline you'd build when processing
hundreds of samples with millions of cells across multiple cohorts.

## How to Use This Skill

1. **Identify your analysis stage** from the pipeline overview below
2. **Read the corresponding reference file** for detailed methods, code, and parameters
3. **Use scripts/** for automated pipeline steps where available
4. **Chain stages sequentially** — each stage depends on outputs from previous stages

## When to Use This Skill

Use for production-scale scRNA-seq tasks including:
- Processing >100K cells across multiple samples/batches/cohorts
- Building or applying multi-resolution cell type hierarchies (L1 broad → L2 intermediate → L3 fine)
- CellTypist-based automated labeling with custom PBMC models
- Two-stage doublet removal: scrublet (computational) + marker-based (biological)
- Iterative subclustering within major cell classes for expert annotation
- Harmony-based cohort/batch integration for large datasets
- Pseudobulk aggregation and DESeq2 multi-factor differential expression
- CLR compositional analysis of cell type frequencies
- Normalizing cell proportions to absolute lymphocyte counts (ALC)
- Multi-omics integration: scRNA-seq with flow cytometry, Olink proteomics, HAI serology
- Longitudinal immune profiling (aging, vaccination responses)
- Partitioning large cell classes by metadata (cohort/sex/CMV) for tractable processing

## Pipeline Overview

```
Stage 0: Sample Selection & Metadata Assembly
│  └── references/sample_selection.md
│
Stage 1: CellTypist Multi-Resolution Labeling
│  ├── L1 (9 broad classes) + L2 (29 intermediate types)
│  ├── L3 (71 high-resolution types)
│  └── references/celltypist_labeling.md
│
Stage 2: QC Filtering & Doublet Detection
│  ├── Scrublet computational doublet detection
│  ├── Threshold-based QC (MT%, gene counts)
│  ├── Marker-based doublet filtering within L2 classes
│  └── references/qc_doublet_filtering.md
│
Stage 3: Iterative Subclustering & Annotation Refinement
│  ├── Cluster within L2 classes (Harmony + Leiden)
│  ├── Partition large classes by cohort/sex/status
│  ├── Expert annotation with marker gene review
│  ├── L3-level subclustering, review, and filtering
│  ├── Label reassignment for misclassified cells
│  └── references/subclustering_annotation.md
│
Stage 4: Data Assembly & Frequency Analysis
│  ├── Final label assembly across all cell classes
│  ├── Cell type frequency tabulation per sample
│  ├── CLR transformation of proportions
│  ├── ALC normalization for abundance estimation
│  └── references/data_assembly.md
│
Stage 5: Pseudobulk DE & Downstream Analysis
│  ├── Pseudobulk aggregation (sum UMI, mean normalized)
│  ├── DESeq2 with multi-factor designs
│  ├── Age/Sex/CMV contrasts
│  ├── Longitudinal paired comparisons
│  └── references/pseudobulk_de.md
│
Stage 6: Multi-Omics Integration & Visualization
│  ├── scRNA-seq + flow cytometry cross-validation
│  ├── Olink proteomics integration
│  ├── HAI/serology correlation
│  ├── Composite aging scores
│  ├── Standardized color palettes for 71 cell types
│  └── references/multiomics_visualization.md
```

## Quick-Start: Standard AIFI PBMC Pipeline

This is the core workflow. Each step links to its reference file.

### Stage 1: CellTypist Labeling (references/celltypist_labeling.md)

```python
import celltypist
import scanpy as sc

# Load and normalize for CellTypist (requires log1p-normalized to 10k)
adata = sc.read_10x_h5("sample.h5")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Apply AIFI models at 3 resolutions
model_files = {
    'AIFI_L1': 'ref_pbmc_clean_celltypist_model_AIFI_L1.pkl',  # 9 classes
    'AIFI_L2': 'ref_pbmc_clean_celltypist_model_AIFI_L2.pkl',  # 29 types
    'AIFI_L3': 'ref_pbmc_clean_celltypist_model_AIFI_L3.pkl',  # 71 types
}
for name, model_path in model_files.items():
    predictions = celltypist.annotate(adata, model=model_path)
    adata.obs[f'{name}_prediction'] = predictions.predicted_labels['predicted_labels']
    # Store prediction scores
    prob = predictions.probability_matrix
    scores = [prob.loc[idx, adata.obs[f'{name}_prediction'][idx]]
              for idx in adata.obs.index]
    adata.obs[f'{name}_score'] = scores
```

### Stage 2: QC & Doublet Detection (references/qc_doublet_filtering.md)

```python
# Scrublet doublet detection
import scanpy.external as sce
sce.pp.scrublet(adata)

# QC metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Serial QC filtering (AIFI thresholds)
adata = adata[~adata.obs['predicted_doublet']].copy()          # scrublet
adata = adata[adata.obs['pct_counts_mt'] <= 10].copy()         # MT% ≤ 10
adata = adata[adata.obs['n_genes_by_counts'] >= 200].copy()    # ≥ 200 genes
adata = adata[adata.obs['n_genes_by_counts'] <= 5000].copy()   # ≤ 5000 genes
```

### Stage 3: Clustering within Cell Classes (references/subclustering_annotation.md)

```python
# Standard AIFI scanpy workflow for subclustering
def aifi_cluster_pipeline(adata, harmony_key='cohort', n_neighbors=50,
                          n_pcs=30, resolutions=[1.0, 1.5, 2.0]):
    """AIFI-standard processing: normalize → HVG → PCA → Harmony → cluster."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg)
    sc.tl.pca(adata_hvg, svd_solver='arpack')
    # Transfer PCA to full object
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

    # Harmony integration
    sce.pp.harmony_integrate(adata, key=harmony_key)

    # Neighbors on harmonized space
    sc.pp.neighbors(adata, n_neighbors=n_neighbors,
                    use_rep='X_pca_harmony', n_pcs=n_pcs)
    sc.tl.umap(adata)

    # Leiden at multiple resolutions
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    return adata
```

### Stage 4: Pseudobulk DE (references/pseudobulk_de.md)

```python
# DESeq2 with multi-factor design: ~ Age + CMV + Sex
# Filter genes: ≥ 10% detection across single cells
# Use summed UMI counts per cell type per sample
# Contrasts for each factor independently
```

## AIFI Cell Type Hierarchy

The AIFI Immune Health Atlas defines 3 levels of immune cell classification for PBMCs:

| Level | Resolution | N Types | Description |
|-------|-----------|---------|-------------|
| AIFI_L1 | Broad | 9 | T cell, B cell, NK cell, Monocyte, DC, ILC, Progenitor, Erythrocyte, Platelet |
| AIFI_L2 | Intermediate | 29 | Naive/Memory CD4/CD8 T, Treg, MAIT, gdT, DN T, Naive/Memory/Effector B, Plasma, CD14/CD16/Intermediate Mono, cDC1/cDC2/ASDC/pDC, CD56bright/CD56dim NK, ILC, Progenitor, Erythrocyte, Platelet, Proliferating T/NK |
| AIFI_L3 | High-resolution | 71 | Fine subtypes with functional markers (e.g., GZMK+ CD27+ EM CD8, ISG+ CD14 Mono, Adaptive NK, SOX4+ naive CD4 T) |

For the complete L3 roster, the marker panel used to resolve each lineage, and the published
threshold for every filter and reassignment, see `references/marker_annotation_library.md`,
`references/l2_doublet_filter_rules.md`, and `references/l3_refinement_rules.md`.

## Key Methodological Decisions

1. **CellTypist over scANVI for labeling**: CellTypist uses logistic regression (OvR for L1/L2, multinomial for L3) — fast, interpretable, and works sample-by-sample without requiring joint embedding. Good when you have a well-curated reference atlas.

2. **Two-stage doublet removal**: Computational (scrublet) catches most doublets; marker-based filtering within cell classes catches the rest using biological knowledge (e.g., MS4A1+ clusters within T cells = B cell doublets).

3. **Harmony over scVI for integration**: When the batch effect is primarily cohort-level (not technology-level), Harmony is fast enough for millions of cells and preserves biological variation well.

4. **Partitioning large cell classes**: For populations >1M cells (e.g., CD4 naive T, CD14 Mono), process separately by cohort × sex × status to keep file sizes and computation tractable.

5. **Serial QC filtering**: Apply filters in sequence and track cell counts removed at each step — this creates an audit trail and prevents compounding filter effects.

6. **CLR over proportions**: Centered log-ratio transformation handles the compositional nature of cell frequency data. Always use CLR for statistical comparisons of cell type abundance.

7. **ALC normalization**: Multiply scRNA-seq proportions by absolute lymphocyte counts from clinical blood counts to estimate true cell abundance — critical for immune monitoring studies.

## Reference Files

| Reference | File | Covers |
|-----------|------|--------|
| Sample Selection | `references/sample_selection.md` | HISE retrieval, cohort filtering, metadata assembly, CMV/BMI annotation |
| CellTypist Labeling | `references/celltypist_labeling.md` | Multi-resolution labeling, custom model training, score interpretation, batch processing |
| QC & Doublet Filtering | `references/qc_doublet_filtering.md` | Scrublet, MT/gene thresholds, marker-based doublet tables, low-quality cluster removal |
| Subclustering & Annotation | `references/subclustering_annotation.md` | AIFI scanpy workflow, cell class partitioning, iterative subclustering, expert annotation, label reassignment |
| Data Assembly | `references/data_assembly.md` | Label assembly, frequency tabulation, CLR transformation, ALC normalization |
| Pseudobulk DE | `references/pseudobulk_de.md` | Pseudobulk aggregation, DESeq2 designs, age/sex/CMV contrasts, longitudinal paired tests |
| Multi-Omics & Visualization | `references/multiomics_visualization.md` | Flow cytometry integration, Olink proteomics, HAI serology, color palettes, helper functions |
| **Marker Annotation Library** | `references/marker_annotation_library.md` | **Full marker reference: broad contamination panels, 7 lineage class panels, per-cell-type marker checks for 66 of 71 L3 types, L3 roster by partition** |
| **L2 Doublet Filter Rules** | `references/l2_doublet_filter_rules.md` | **All 113 published thresholds across 29 L2 classes, mean-expression HBB rules, duplicate-key caveat** |
| **L3 Refinement Rules** | `references/l3_refinement_rules.md` | **84 remove/reassign rules across 49 L3 types, scrublet-score filters, CD4/CD8 confusion patterns** |

## Scripts

| Script | Purpose |
|--------|---------|
| `scripts/celltypist_batch_labeling.py` | Batch CellTypist labeling across samples at L1/L2/L3 |
| `scripts/aifi_qc_pipeline.py` | Serial QC filtering with audit trail |
| `scripts/aifi_cluster_pipeline.py` | Standard AIFI: normalize → HVG → PCA → Harmony → Leiden |
| `scripts/marker_doublet_filter.py` | Marker-based doublet filtering with configurable thresholds |
| `scripts/pseudobulk_deseq2.py` | Pseudobulk aggregation and DESeq2 via rpy2 |
| `scripts/clr_frequency_analysis.py` | Cell type frequency tabulation with CLR and ALC normalization |
| `scripts/aifi_markers.py` | **Importable marker library: panels, published thresholds, and AIFI's `remove_cl`/`extract_cl` primitives** |

## Assets (machine-readable)

| Asset | Contents |
|-------|----------|
| `assets/aifi_marker_library.json` | Lineage panels with inline annotations, per-cell-type marker checks, L3 roster |
| `assets/aifi_l2_doublet_thresholds.json` | 113 fraction rules across 29 L2 classes + 4 mean-expression rules, source order preserved |
| `assets/aifi_l3_refinement_rules.csv` | 84 L3 rules: lineage, cell_type, gene, direction, cutoff, update_type, change_to |

All extracted verbatim from `aifimmunology/sound-life-scrna-analysis` (branch `development`,
notebooks `02-reference_labeling/07`, `10a`-`13g`) and cross-checked against the published
methods tables. Counts stated here were verified programmatically at extraction time: the
L3 roster reproduces exactly 71 populations and the L2 dict exactly 29 classes, matching
AIFI's stated hierarchy.

## Key Resources

- [Gong et al., Nature 2025](https://doi.org/10.1038/s41586-025-09686-5) — Multi-omic profiling of age-related immune dynamics
- [AIFI Sound Life Analysis](https://github.com/aifimmunology/sound-life-scrna-analysis) — Full analysis notebooks
- [AIFI Immune Health Atlas](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/) — Reference atlas and models
- [IHA Figure Code](https://github.com/aifimmunology/IHA-Figure) — Manuscript figure generation
- [CellTypist](https://www.celltypist.org/) — Automated cell type annotation
- [PALMO](https://github.com/aifimmunology/PALMO) — Longitudinal single-cell analysis

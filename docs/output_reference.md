# Output Reference

ViralScan writes all results to the directory specified by `--output / -o`.

---

## Directory layout

```
output/
├── config.yaml                  # Per-run configuration written by ViralScan
├── kb-python/                   # Raw kb count output
│   ├── output.bus               # BUS file
│   ├── run_info.json
│   └── ...
├── results/
│   ├── viral_summary.tsv        # Per-virus summary table
│   ├── per_cell_viral.tsv       # Per-barcode viral UMI table
│   ├── report.html              # Self-contained HTML report
│   ├── adata_original.h5ad      # AnnData before multimapping correction
│   ├── adata_multimap.h5ad      # AnnData after multimapping correction
│   └── plots/                   # PNG visualisations
│       ├── histogram_<virus>.png
│       └── ...
└── logs/
    └── snakemake.log
```

---

## `viral_summary.tsv`

Tab-separated, one row per detected virus.

| Column | Description |
|--------|-------------|
| `virus_name` | Human-readable virus name |
| `total_umi` | Total viral UMI across all cells |
| `infected_cells` | Number of cells with ≥ `--detection-threshold` viral UMI |
| `total_cells` | Total cells passing QC |
| `pct_infected` | `infected_cells / total_cells × 100` |
| `umi_per_10k` | `total_umi / total_umi_all × 10 000` |
| `cluster_pvalue` | Fisher's exact test p-value from `viral_neighbor_enrichment` |

---

## `per_cell_viral.tsv`

Tab-separated, one row per cell × detected virus combination.

| Column | Description |
|--------|-------------|
| `barcode` | Cell barcode |
| `virus_name` | Virus name |
| `viral_umi` | Viral UMI count for this cell |
| `total_umi` | Total UMI count for this cell |
| `viral_fraction` | `viral_umi / total_umi` |

---

## `report.html`

A self-contained HTML file with:

- Run metadata (date, parameters, sample names)
- QC summary table
- Per-virus detection table with normalised metrics
- Embedded histogram plots (base64 PNG)
- Interpretation guidance
- `viral_neighbor_enrichment` p-values

Open in any modern browser — no internet connection required.

---

## AnnData files (`.h5ad`)

The AnnData objects can be loaded with [scanpy](https://scanpy.readthedocs.io/):

```python
import scanpy as sc

adata = sc.read_h5ad("output/results/adata_multimap.h5ad")
print(adata)
# Layers: counts_original, counts_corrected
```

Key layers:

| Layer | Description |
|-------|-------------|
| `counts_original` | Raw kb count matrix (unique-mapping reads) |
| `counts_corrected` | Extra multi-mapped read share (additive correction) |

`adata.X` = `counts_original + counts_corrected` (combined count matrix).

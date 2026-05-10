# Output Reference

ViralScan writes one sample directory under the path passed to `--output / -o`.
The sample directory is inferred from the R1 FASTQ filename before the first
underscore. For `sample_R1.fastq.gz`, the run directory is `output/sample/`.

---

## Directory layout

```
output/
└── sample/
    ├── config.yaml
    ├── summary.txt
    ├── report.html
    ├── log/
    │   ├── analysis.txt
    │   ├── detection.done
    │   ├── kb.done
    │   ├── multimap.done
    │   └── umap.done
    ├── kb-python/
    │   ├── counts_unfiltered/
    │   │   ├── adata.h5ad
    │   │   └── adata_multimap.h5ad
    │   ├── output.bus
    │   ├── run_info.json
    │   └── ...
    ├── plots/
    │   ├── <virus>_histogram.png
    │   ├── SuperExpressor_<virus>.png
    │   ├── umap_binary.html
    │   └── umap_continuous.html
    └── results/
        ├── viral_summary.tsv
        ├── per_cell_viral.tsv
        └── cell_type_enrichment.tsv
```

`cell_type_enrichment.tsv` is present only when `--cell-types` is supplied.
UMAP files are present only when `--umap` is supplied. `host_filtered/` is
present only when `--host-filter` is supplied.

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

## `cell_type_enrichment.tsv`

Tab-separated, one row per detected virus and labeled cell type. Written only
when `--cell-types cell_types.csv` is supplied.

| Column | Description |
|--------|-------------|
| `virus` | Virus name |
| `cell_type` | Cell-type label from the CSV |
| `n_infected` | Infected labeled cells in this cell type |
| `n_total` | Total labeled cells of this type |
| `pct` | `n_infected / n_total × 100` |
| `OR` | One-sided Fisher exact odds ratio |
| `pvalue` | Raw Fisher exact p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |

Input CSV requirements:

```csv
barcode,cell_type
AAACCCAAGAGT-1,T cell
AAACCCAGTGCA-1,Monocyte
```

Barcodes must match `adata.obs_names`. If no barcodes overlap, ViralScan skips
the enrichment table and logs a warning.

---

## AnnData files (`.h5ad`)

The AnnData objects can be loaded with [scanpy](https://scanpy.readthedocs.io/):

```python
import scanpy as sc

adata = sc.read_h5ad("output/sample/kb-python/counts_unfiltered/adata_multimap.h5ad")
print(adata)
# Layers: counts_original, counts_corrected
```

Key layers:

| Layer | Description |
|-------|-------------|
| `counts_original` | Raw kb count matrix (unique-mapping reads) |
| `counts_corrected` | Extra multi-mapped read share (additive correction) |

`adata.X` = `counts_original + counts_corrected` (combined count matrix).

---

## Multiple samples

When `--sample1` and `--sample2` contain comma-separated FASTQ lists,
ViralScan processes each pair separately:

```bash
viralscan \
  -t t2g.txt -i index.idx -o output/ \
  -s1 A_R1.fastq.gz,B_R1.fastq.gz \
  -s2 A_R2.fastq.gz,B_R2.fastq.gz
```

Expected directories:

```text
output/A/
output/B/
```

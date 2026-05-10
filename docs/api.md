# API Reference

This page documents the importable ViralScan API. Most users should use the
`viralscan` CLI; import these functions when you are writing tests, notebooks,
or small automation around ViralScan outputs.

Snakemake worker scripts such as `scripts/detection.py`, `scripts/multimap.py`,
and `scripts/umap.py` depend on Snakemake globals (`snakemake.input`,
`snakemake.config`, and similar). Treat those files as workflow implementation
details rather than a stable Python API.

---

## Common imports

```python
from viralscan.enrichment import cell_type_enrichment, write_cell_type_enrichment
from viralscan.scripts.ncbi_fetch import fetch_reference
from viralscan.scripts.build_reference import build_combined_reference, fetch_host_cdna
from viralscan.data_fetch import fetch_viral_data, ensure_viral_data
```

---

## Cell-Type Enrichment

Use this API when you already have a ViralScan AnnData object and want to
compute or re-run enrichment by cell type outside the full Snakemake workflow.

### `cell_type_enrichment`

```python
def cell_type_enrichment(
    adata: anndata.AnnData,
    group_by_virus: dict[str, list[str]],
    cfg: dict[str, Any],
) -> pandas.DataFrame: ...
```

Computes one-sided Fisher exact tests for each virus and cell type. A cell is
treated as infected for a virus when the sum across that virus's listed
features is greater than zero.

Minimal example:

```python
import scanpy as sc
from viralscan.enrichment import cell_type_enrichment

adata = sc.read_h5ad("output/sample/kb-python/counts_unfiltered/adata_multimap.h5ad")
group_by_virus = {
    "Influenza A virus": ["INFLUENZA_A_PB2", "INFLUENZA_A_NS1"],
}
cfg = {"cell_types": "cell_types.csv"}

df = cell_type_enrichment(adata, group_by_virus, cfg)
print(df.sort_values("padj").head())
```

`cell_types.csv` must contain `barcode` and `cell_type` columns:

```csv
barcode,cell_type
AAACCCAAGAGT-1,T cell
AAACCCAGTGCA-1,Monocyte
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `adata` | `anndata.AnnData` | AnnData object whose `obs_names` are cell barcodes. |
| `group_by_virus` | `dict[str, list[str]]` | Virus name to feature names in `adata.var_names`. Missing features are ignored. |
| `cfg` | `dict[str, Any]` | Must include `"cell_types"` with the CSV path. |

Returns an empty `DataFrame` when the CSV is missing, unreadable, empty, lacks
overlapping barcodes, or none of the requested viral features are present.

Returned columns:

| Column | Description |
|--------|-------------|
| `virus` | Virus name from `group_by_virus`. |
| `cell_type` | Cell-type label from the CSV. |
| `n_infected` | Infected labeled cells in this cell type. |
| `n_total` | Total labeled cells of this type. |
| `pct` | `n_infected / n_total * 100`, rounded to 4 decimals. |
| `OR` | One-sided Fisher exact odds ratio. |
| `pvalue` | Raw Fisher exact p-value. |
| `padj` | Benjamini-Hochberg adjusted p-value. |

### `write_cell_type_enrichment`

```python
def write_cell_type_enrichment(
    cell_type_df: pandas.DataFrame,
    outputpath: str,
) -> str | None: ...
```

Writes `cell_type_df` to `{outputpath}/results/cell_type_enrichment.tsv`.
Pass the sample run directory, not the top-level multi-sample output root.

```python
from viralscan.enrichment import write_cell_type_enrichment

path = write_cell_type_enrichment(df, "output/sample")
print(path)
```

Returns the written path as a string, or `None` when the input table is empty.

### `_bh_adjust`

```python
def _bh_adjust(pvals: list[float] | numpy.ndarray) -> numpy.ndarray: ...
```

Internal helper for Benjamini-Hochberg FDR correction. It is documented for
contributors, but downstream code should prefer `cell_type_enrichment()`.

---

## Reference Downloads

These functions support the CLI reference modes and are useful in automation
when you want to prepare reference files before running `viralscan`.

### `fetch_reference`

```python
def fetch_reference(
    accessions: Iterable[str],
    out_dir: str | os.PathLike[str],
    email: str | None = None,
    api_key: str | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
) -> tuple[pathlib.Path, pathlib.Path]: ...
```

Downloads FASTA and minimal GTF files for one or more NCBI nucleotide
accessions, merges them into `reference.fasta` and `reference.gtf`, and caches
per-accession downloads.

```python
from viralscan.scripts.ncbi_fetch import fetch_reference

fasta, gtf = fetch_reference(
    ["NC_002021.3", "NC_045512.2"],
    out_dir="refs/ncbi_virus",
    email="you@example.org",
)
print(fasta, gtf)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `accessions` | `Iterable[str]` | required | RefSeq or GenBank nucleotide accessions. |
| `out_dir` | `str \| PathLike` | required | Directory for merged `reference.fasta` and `reference.gtf`. |
| `email` | `str \| None` | `None` | NCBI contact email. Falls back to `NCBI_EMAIL`. |
| `api_key` | `str \| None` | `None` | Optional NCBI API key. Falls back to `NCBI_API_KEY`. |
| `cache_dir` | `str \| PathLike \| None` | `None` | Cache directory. Defaults to `~/.cache/viralscan/ncbi/`. |

Raises `NCBIFetchError` for missing accessions, missing email, failed
downloads, parse failures, or empty merged outputs.

### `fetch_host_cdna`

```python
def fetch_host_cdna(
    species: str,
    out_dir: str | os.PathLike[str],
    cache_dir: str | os.PathLike[str] | None = None,
) -> tuple[pathlib.Path, pathlib.Path]: ...
```

Downloads the current Ensembl cDNA FASTA and GTF for a supported host species
and copies them into `out_dir`.

```python
from viralscan.scripts.build_reference import fetch_host_cdna

host_fasta, host_gtf = fetch_host_cdna("human", "refs/human")
```

Run `viralscan build-ref --list-species` for supported names.

### `build_combined_reference`

```python
def build_combined_reference(
    host_species: str,
    virus_accessions: list[str],
    out_dir: str | os.PathLike[str],
    email: str | None = None,
    api_key: str | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
    run_kb_ref: bool = True,
) -> dict[str, pathlib.Path | None]: ...
```

Builds a host + virus reference in the same way as `viralscan build-ref`.
When `run_kb_ref=True` and `kb` is on `PATH`, the returned dictionary includes
the kallisto index and `t2g.txt` paths.

```python
from viralscan.scripts.build_reference import build_combined_reference

ref = build_combined_reference(
    host_species="human",
    virus_accessions=["NC_045512.2"],
    out_dir="ref_human_sars",
    email="you@example.org",
)

print(ref["index"], ref["t2g"])
```

Returned keys:

| Key | Description |
|-----|-------------|
| `fasta` | `combined.fa` |
| `gtf` | `combined.gtf` |
| `index` | `index.idx`, or `None` when `kb ref` was skipped or failed |
| `t2g` | `t2g.txt`, or `None` when `kb ref` was skipped or failed |

---

## Viral Panel Cache

Use these functions when application code needs to locate or populate the
external bundled GTF panel.

### `fetch_viral_data`

```python
def fetch_viral_data(
    cache_dir: str | pathlib.Path | None = None,
    archive_url: str | None = None,
    expected_sha256: str | None = None,
    force: bool = False,
) -> pathlib.Path: ...
```

Downloads the Zenodo archive, verifies checksums, extracts `.gtf` files, and
writes a manifest. Returns the data directory, usually
`~/.cache/viralscan/data/`.

```python
from viralscan.data_fetch import fetch_viral_data

data_dir = fetch_viral_data()
```

### `ensure_viral_data`

```python
def ensure_viral_data(cache_dir: str | pathlib.Path | None = None) -> pathlib.Path: ...
```

Returns the cached GTF directory if it is valid. Raises `ViralScanDataError`
with the user-facing `viralscan data fetch` instruction when the cache is
missing or incomplete.

---

## Shared Helpers

### `configure_logging`

```python
def configure_logging(verbose: bool = False, quiet: bool = False) -> None: ...
```

Configures the `viralscan` logger. `quiet=True` takes priority over
`verbose=True`.

### `setup_script_logging`

```python
def setup_script_logging() -> logging.Logger: ...
```

Initializes logging for Snakemake worker subprocesses and returns the
`viralscan` logger.

### `load_config`

```python
def load_config(path: str | pathlib.Path) -> dict[str, Any]: ...
```

Loads a YAML configuration file and returns a plain dictionary. Raises
`ValueError` when the top-level YAML object is not a mapping.

---

## Constants

### `VIRUS_NAME_MAP`

```python
VIRUS_NAME_MAP: dict[str, str]
```

Maps viral gene-ID prefixes used in GTF files to human-readable virus names.

```python
from viralscan.constants import VIRUS_NAME_MAP

print(VIRUS_NAME_MAP.get("INFLUENZA_A", "INFLUENZA_A"))
```

### `ENSEMBL_SPECIES`

```python
ENSEMBL_SPECIES: dict[str, tuple[str, str]]
```

Maps short species names such as `human` or `mouse` to Ensembl species names
and assemblies used by `build_combined_reference()`.

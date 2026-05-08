# API Reference

This page documents the importable public API of ViralScan.

The Snakemake worker scripts (`scripts/analysis.py`, `scripts/detection.py`,
`scripts/multimap.py`, `scripts/umap.py`, `scripts/createconfig.py`) run
inside Snakemake's execution environment and rely on its magic globals
(`snakemake.input`, `snakemake.config`, etc.).  They are not directly
importable and have no stable public API surface; they are therefore not
documented here.  See the source files for implementation details.

---

## `viralscan.utils`

Small shared helpers used across the Snakemake worker scripts and the CLI.

### `configure_logging`

```python
def configure_logging(verbose: bool = False, quiet: bool = False) -> None: ...
```

Configure the root `viralscan` logger.

Called once from `menu.py` after CLI argument parsing.  Each Snakemake worker
script runs in its own subprocess and calls this through `setup_script_logging()`
below.  Priority ordering: `quiet` (WARNING) › `verbose` (DEBUG) › default (INFO).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `verbose` | `bool` | `False` | Set log level to DEBUG. |
| `quiet` | `bool` | `False` | Set log level to WARNING, suppressing INFO messages. |

**Returns:** `None`

---

### `setup_script_logging`

```python
def setup_script_logging() -> logging.Logger: ...
```

Minimal logging bootstrap for Snakemake worker scripts.

Each script runs in a fresh interpreter process and must configure its own
handler.  Calls `configure_logging()` with default settings and returns the
`viralscan` logger ready to use.

**Returns:** `logging.Logger` — the `"viralscan"` logger instance.

---

### `load_config`

```python
def load_config(path: Union[str, Path]) -> dict[str, Any]: ...
```

Read a YAML config file and return it as a plain `dict`.

All ViralScan worker scripts use this central loader so that YAML parsing and
boolean-normalisation logic live in one place.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | `str \| Path` | required | Path to the YAML config file. |

**Returns:** `dict[str, Any]` — parsed YAML mapping.

**Raises:** `ValueError` if the file does not contain a YAML mapping at the top level.

---

## `viralscan.constants`

### `VIRUS_NAME_MAP`

```python
VIRUS_NAME_MAP: dict[str, str]
```

Maps the abbreviated GTF gene-ID prefix to the human-readable virus name.
Contains entries for all 195 viruses bundled with ViralScan.

Previously duplicated verbatim in `detection.py` and `umap.py`; both modules
now import this single source of truth.

**Example:**

```python
from viralscan.constants import VIRUS_NAME_MAP

print(VIRUS_NAME_MAP["INFLUENZA_A"])  # "Influenza A virus"
```

---

## `viralscan.enrichment`

Cell-type enrichment utilities extracted from `scripts/detection.py` so they
can be imported and unit-tested without triggering Snakemake magic globals.

### `cell_type_enrichment`

```python
def cell_type_enrichment(
    adata: anndata.AnnData,
    group_by_virus: dict[str, list[str]],
    cfg: dict[str, Any],
) -> pd.DataFrame: ...
```

Compute per-virus enrichment by cell type using one-sided Fisher exact tests.

Reads the barcode-to-cell-type mapping from the CSV path stored in
`cfg["cell_types"]`.  For each detected virus and each cell type, it tests
whether infected cells are over-represented in that cell type relative to all
other labeled cells.  Adjusted p-values use the Benjamini-Hochberg procedure.

Returns an empty `DataFrame` if `cfg.get("cell_types")` is falsy, if the CSV
cannot be read, or if no barcodes overlap between the AnnData and the CSV.

| Parameter | Type | Description |
|-----------|------|-------------|
| `adata` | `anndata.AnnData` | AnnData object whose `obs_names` are cell barcodes. |
| `group_by_virus` | `dict[str, list[str]]` | Mapping of virus name → list of gene/feature names in `adata.var_names`. |
| `cfg` | `dict[str, Any]` | ViralScan config dict; must contain `"cell_types"` key with path to the CSV. |

**Returns:** `pd.DataFrame` with columns:

| Column | Description |
|--------|-------------|
| `virus` | Virus name. |
| `cell_type` | Cell type label. |
| `n_infected` | Infected cells in this cell type. |
| `n_total` | Total labeled cells of this type. |
| `pct` | Percentage of infected cells in this type. |
| `OR` | Fisher exact odds ratio. |
| `pvalue` | One-sided Fisher exact p-value. |
| `padj` | Benjamini-Hochberg adjusted p-value. |

---

### `write_cell_type_enrichment`

```python
def write_cell_type_enrichment(
    cell_type_df: pd.DataFrame,
    outputpath: str,
) -> str | None: ...
```

Write the cell-type enrichment table to `{outputpath}/results/cell_type_enrichment.tsv`.

Skips writing (returns `None`) if `cell_type_df` is `None` or empty.

| Parameter | Type | Description |
|-----------|------|-------------|
| `cell_type_df` | `pd.DataFrame` | DataFrame returned by `cell_type_enrichment()`. |
| `outputpath` | `str` | Root output directory (same as `--output`). |

**Returns:** `str` path to the written TSV, or `None` if nothing was written.

---

### `_bh_adjust` *(internal)*

```python
def _bh_adjust(
    pvals: list[float] | npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]: ...
```

**Internal helper — not part of the public API.**

Benjamini-Hochberg FDR correction for a 1D array of p-values.  Used by
`cell_type_enrichment()` to compute `padj`.  Documented here for contributors.

| Parameter | Type | Description |
|-----------|------|-------------|
| `pvals` | `list[float] \| NDArray[float64]` | Raw p-values. |

**Returns:** `NDArray[float64]` — BH-adjusted p-values clipped to `[0, 1]`.

---

## `viralscan.scripts.ncbi_fetch`

NCBI E-utilities downloader that fetches a FASTA + GTF reference for one or
more nucleotide accession numbers without requiring Biopython.  Downloads are
cached under `~/.cache/viralscan/ncbi/` by default.

### `NCBIFetchError`

```python
class NCBIFetchError(RuntimeError): ...
```

Raised when an NCBI download or parse step fails (invalid accession, HTTP
error after retries, empty response, etc.).

---

### `fetch_reference`

```python
def fetch_reference(
    accessions: Iterable[str],
    out_dir: str | os.PathLike[str],
    email: str | None = None,
    api_key: str | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
) -> tuple[Path, Path]: ...
```

Download FASTA + GTF for one or more NCBI nucleotide accessions and merge
them into a single `reference.fasta` and `reference.gtf` in `out_dir`.

Multiple accessions are concatenated in order, mirroring ViralScan's
comma-separated `--fasta`/`--gtf` semantics.  Per-accession files are cached
so re-runs do not re-download.  Retries with exponential back-off on transient
HTTP errors (429, 5xx).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `accessions` | `Iterable[str]` | required | Nucleotide accession strings (RefSeq or GenBank), e.g. `"NC_002021.3"`. |
| `out_dir` | `str \| PathLike` | required | Directory where merged output files are written. Created if absent. |
| `email` | `str \| None` | `None` | Contact email required by NCBI terms of service; falls back to `$NCBI_EMAIL`. |
| `api_key` | `str \| None` | `None` | NCBI API key for higher rate limits; falls back to `$NCBI_API_KEY`. |
| `cache_dir` | `str \| PathLike \| None` | `None` | Per-accession cache directory; defaults to `~/.cache/viralscan/ncbi/`. |

**Returns:** `tuple[Path, Path]` — `(fasta_path, gtf_path)` pointing to the
merged files, ready to pass to `kb ref`.

**Raises:** `NCBIFetchError` if no accessions are provided, if no email is
available, if NCBI returns an error after retries, or if the merged output is
empty.

# Inspect results

After `viralscan validate-run --run-dir output/` reports `"ok": true`, open
`output/<sample>/report.html` for the overview and then read the tables below.
This guide covers the main results, the units they are in, and the denominator
behind every per-cell rate. Exact column definitions are in the
[Output reference](output_reference.md).

## Which file answers which question

| Question | File | First fields to read |
|---|---|---|
| What did this run find? | `results/viral_summary.tsv` | `virus_name`, `viral_molecules_total_est`, `infected_called`, `n_called_cells`, `pct_infected_called`, `detected` |
| Which cells carry a candidate? | `results/per_cell_viral.tsv` | `barcode`, `virus_name`, `viral_molecules_total_est`, `viral_fraction`, `is_called_cell` |
| What is the per-cell denominator? | `results/called_cells.tsv` | `barcode` (one row per called cell) |
| Could a negative be real? | `results/sensitivity.tsv` | `virus_name`, `observed_molecules`, `detected`, `capture`, `lod95_molecules`, `informative_negative`, `negative_blockers` |
| What does the report show? | `report.html` | Run metadata, QC summary, per-virus detection table, multimapping evidence, histograms |
| What did the run actually do? | `run_manifest.json`, `run_complete.json` | Options, input hashes, software identity, completion marker and artifact hashes |
| Which reference was used? | `results/reference_provenance.json` | Index, `t2g`, GTF, technology, multimapping settings, viral accessions |

`report.html` is self-contained and works offline. Rows with `detected = false`
are excluded from its virus overlays. `summary.txt` gives the same per-virus
figures as plain text and marks a candidate that misses the called-cell gate as
`CANDIDATE ONLY`.

## Units: molecules, not reads

ViralScan counts corrected **(cell barcode, UMI) molecules**, not reads. PCR and
sequencing duplicate reads for the same molecule are collapsed before counting;
BUS `count` is read multiplicity and is audit-only.

`viral_molecules_total_est` has two components:

- **unique support** — the molecule's compatible gene set had exactly one gene;
- **allocated ambiguous support** — the molecule was compatible with several
  genes and the selected multimapping method allocated it. The allocation is a
  model estimate and moves if the method changes.

In the count matrix, `counts_unique` and `counts_ambiguous_allocated` are
disjoint and sum to `adata.X`. The separate diagnostic layers in the AnnData
file (below) keep the raw components visible.

## `detected` is a configured threshold result

`detected = true` means the virus passed the gate configured for this run:

- at least `--detection-min-cells` called cells are infected with the virus
  (default 1), **and**
- the called cells hold at least `--detection-threshold` of its molecules in
  total (default 1).

A cell is infected for a virus when its summed allocated molecules reach
`--infected-min-molecules` (default 1.0). Weaker cells are counted in
`infected_subthreshold_cells`, never dropped.

`detected = false` means *candidate only*: the virus had signal in the
all-barcode screen but missed the called-cell gate. It is not proof of absence.
Passing the gate is candidate evidence, not a confirmed infection.

## Denominators: all barcodes versus called cells

| Field | Numerator / denominator | Use |
|---|---|---|
| `pct_infected` | infected cells over **all barcodes** | Screening view; diluted by empty droplets |
| `pct_infected_called` | infected called cells over `n_called_cells` | Within-run per-cell candidate-support rate |
| `pct_infected_comparable` | Infected cells within the comparable set over called barcodes clearing an absolute host-UMI floor (200 host molecules) | Comparing runs that used different host-filtering strategies |

`pct_infected_called` divides by the called cells of *this run*. In a two-step
host-subtraction run, cell calling happens after filtering, so the denominator
moves with the strategy: host filtering can change called cells and viral
molecules by different factors. Do
not compare `pct_infected_called` across host-filtering strategies.

`pct_infected_comparable` uses a denominator that host filtering cannot move. Its
denominator, numerator and percentage are **blank (unavailable), not zero**, when
no cell clears the host-depth floor, when host depth is unavailable, and for
virus-only matrices or external overrides without a host matrix. Cross-arm
comparisons should use the same frozen reference cell set for every arm.

None of these fields is proof of infection. The field names are retained for
schema compatibility.

## Reading one row

The values below are **illustrative**, not from a real sample.

| virus_name | viral_molecules_total_est | infected_called | n_called_cells | pct_infected_called | background_umis | detected |
|---|---:|---:|---:|---:|---:|---|
| Epstein-Barr virus | 412.0 | 37 | 9,850 | 0.38 | 19.5 | true |
| Human herpesvirus 6B | 6.0 | 0 | 9,850 | 0.00 | 6.0 | false |

Read the first row as: 412 estimated viral molecules in total, of which 37
called cells carry support, over 9,850 called cells (0.38 %), with 19.5 molecules
outside the called set. It passed the gate, so `detected` is true.

Read the second row as: 6 estimated molecules, all of them in barcodes outside
the called set (`background_umis = 6.0`), and no infected called cell. It is a
candidate that failed the called-cell gate. That pattern is consistent with
ambient RNA, barcode swapping or empty droplets; it does not establish that the
virus is absent.

For each row you intend to use:

1. Filter on `detected` to separate the detected set from candidates.
2. Read `n_called_cells` and `background_umis` before any percentage.
3. Check `role`, `reference_risk_flags`, `artifact_risk` and
   `sibling_crossmap_note`. These are labels, never filters.
4. Check `host_homology_status`; `not_measured` is not zero homology.
5. Before treating a call as biology, run the checks in
   [Assess the evidence](assess_evidence.md).

## The AnnData file

For downstream analysis, load the corrected matrix:

```python
import scanpy as sc

adata = sc.read_h5ad(
    "output/sample/kb-python/counts_unfiltered/adata_multimap.h5ad"
)
print(adata)
```

`adata.X` is the complete selected-method matrix. `obs` carries per-cell
metadata and `layers` keeps the components, including `counts_unique`,
`counts_ambiguous_allocated`, the alternative allocation layers, and the
host–virus ambiguity diagnostics. The [Output reference](output_reference.md)
lists every layer and column.

## Missing is not zero

An empty field means the measurement was unavailable, not that the value was
zero or that the result was clean. Examples: `pct_infected_comparable` when host
depth cannot be established, `capture` when no control certifies the row, and
the anellovirus alignment columns when that branch did not run. Report these as
unavailable.

# ViralScan examples

Task-oriented, runnable examples of what ViralScan is for. Each maps to a real
subcommand or feature and uses only v3 APIs and output columns. The six
CI-executed notebooks run on small synthetic matrices and quote no real-library
numbers; each includes a negative or ambiguous example next to the positive one.

These notebooks illustrate mechanics. They are **not** real-data performance
validation: none of them measures sensitivity or specificity on a real library.
The two heavy notebooks use public data and are not executed in CI.

| # | Example | Answers | Key command / API | Data | CI |
|---|---------|---------|-------------------|------|----|
| 1 | [Quickstart: from FASTQ to per-cell viral load](quickstart_fastq_to_viral_load.ipynb) | "How do I quantify virus per cell?" | `viralscan` (default), `check-whitelist` | Public data (not run in CI) | `[skip-ci]` |
| 2 | [Building a reference: fetch a panel or build your own](building_a_reference.ipynb) | "Where does the reference come from?" | `data fetch`, `build-ref` | Public data (not run in CI) | `[skip-ci]` |
| 3 | [Multimapping correction](multimapping_correction.ipynb) | "Why does it find more virus than unique counting?" | `--multimap-method em-global`, `rerun-multimap`, `em_gene_abundances`; includes an ambiguous host-virus case | Synthetic | ✅ |
| 4 | [Cell-calling & denominators](cell_calling_denominators.ipynb) | "Is my infection rate a denominator artifact?" | `--cell-calling external`, `compute_stats`; wrong-list failure example | Synthetic | ✅ |
| 5 | [Cell-type enrichment](cell_type_enrichment.ipynb) | "Which cell types are infected?" | `--cell-type-labels`, `cell_type_enrichment` | Synthetic | ✅ |
| 6 | [Specificity: true negative](specificity_true_negative.ipynb) | "Will it cry wolf when the virus is absent?" | `check-whitelist`, `whitelist_match_rate` | Synthetic (negative example) | ✅ |
| 7 | [QC & read evidence](qc_and_read_evidence.ipynb) | "Can I trust this call?" | `evidence`, `check_sibling_crossmapping` (ambiguous EM-bleed example), breadth/ambiguity | Synthetic | ✅ |
| 8 | [Host-response with depth control](host_response_depth_control.ipynb) | "Do infected cells have a host program (really)?" | `hostresponse`, depth-confound diagnostics (synthetic null) | Synthetic | ✅ |

## Suggested reading order

- **New users:** 1 → 2 → 4 → 5.
- **Evaluating the method:** 3 (multimapping allocation) and 4 (denominators) are
  the scientific differentiators; 6 shows specificity.
- **Advanced / analysis:** 7 (trust a call) and 8 (host-response, and how to
  avoid the sequencing-depth trap).

## Running them

```bash
# CI-runnable notebooks need only the Python stack (numpy, pandas, scipy,
# anndata, matplotlib) plus ViralScan on the path:
PYTHONPATH=src jupyter nbconvert --to notebook --execute \
    docs/vignettes/multimapping_correction.ipynb
```

The `[skip-ci]` notebooks additionally require `kb` (kb-python) and `snakemake`.
The quickstart notebook caches its public 1M-read subsample in
`VS_FASTQ_CACHE` (default `~/.cache/viralscan/vignettes/`) and only prints the
download commands when the cache is empty, so the FASTQ is fetched once rather
than on every run. Point `VS_FASTQ_CACHE` at a shared read-only directory to
reuse one copy across users. CI executes only the six synthetic notebooks; the
two `[skip-ci]` notebooks never download data in CI.

See [vignette design plan](https://github.com/mdmanurung/ViralScan/blob/main/docs/vignettes/VIGNETTES_PLAN.md) for the design rationale and data strategy.

# ViralScan vignettes

Task-oriented, runnable examples of what ViralScan is for. Each maps to a real
subcommand/feature and uses only v3 APIs and output columns. All executed cells run
on small synthetic data and quote no real-library numbers. Each CI notebook includes
a negative or ambiguous example next to the positive one.

Most notebooks **execute in CI** on small synthetic or committed data. Two build
a kallisto index / run `kb count` and are marked `[skip-ci]`; they carry
copy-paste public-data commands (no institutional paths).

| # | Vignette | Answers | Key command / API | CI |
|---|----------|---------|-------------------|----|
| 1 | [Quickstart: FASTQ → viral load](quickstart_fastq_to_viral_load.ipynb) | "How do I quantify virus per cell?" | `viralscan` (default), `check-whitelist` | `[skip-ci]` |
| 2 | [Building a reference](building_a_reference.ipynb) | "Where does the reference come from?" | `data fetch`, `build-ref` | `[skip-ci]` |
| 3 | [Multimapping correction](multimapping_correction.ipynb) | "Why does it find more virus than unique counting?" | `--multimap-method em-global`, `rerun-multimap`, `em_gene_abundances`; includes an ambiguous host-virus case | ✅ |
| 4 | [Cell-calling & denominators](cell_calling_denominators.ipynb) | "Is my infection rate a denominator artifact?" | `--cell-calling external`, `compute_stats`; wrong-list failure example | ✅ |
| 5 | [Cell-type enrichment](cell_type_enrichment.ipynb) | "Which cell types are infected?" | `--cell-types`, `cell_type_enrichment` | ✅ |
| 6 | [Specificity: true negative](specificity_true_negative.ipynb) | "Will it cry wolf when the virus is absent?" | `check-whitelist`, `whitelist_match_rate` | ✅ (negative example) |
| 7 | [QC & read evidence](qc_and_read_evidence.ipynb) | "Can I trust this call?" | `evidence`, `check_sibling_crossmapping` (ambiguous EM-bleed example), breadth/ambiguity | ✅ |
| 8 | [Host-response with depth control](host_response_depth_control.ipynb) | "Do infected cells have a host program (really)?" | `hostresponse`, depth-confound diagnostics (synthetic null) | ✅ |

## Suggested reading order

- **New users:** 1 → 2 → 4 → 5.
- **Evaluating the method:** 3 (multimapping allocation) and 4 (denominators) are the
  scientific differentiators; 6 shows specificity.
- **Advanced / analysis:** 7 (trust a call) and 8 (host-response, and how to avoid
  the sequencing-depth trap).

## Running them

```bash
# CI-runnable notebooks need only the Python stack (numpy, pandas, scipy,
# anndata, matplotlib) plus ViralScan on the path:
PYTHONPATH=src jupyter nbconvert --to notebook --execute \
    docs/vignettes/multimapping_correction.ipynb
```

The `[skip-ci]` notebooks additionally require `kb` (kb-python) and `snakemake`,
and download a public FASTQ sample; edit the `WORKDIR` variable to any local
directory.

See [`VIGNETTES_PLAN.md`](VIGNETTES_PLAN.md) for the design rationale and data strategy.

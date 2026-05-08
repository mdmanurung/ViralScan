---
description: >
  Rigorous scientific audit of the ViralScan package.
  Use when reviewing whether the viral quantification pipeline implements scientifically
  valid and defensible methods: kallisto/bustools count correctness, multimapping correction
  logic, detection thresholds, normalization before UMAP, UMI double-counting, GTF parsing,
  barcode handling, clustering statistics, edge cases, silent failures, and reproducibility.
  Read-only; writes findings to a dated audit file at audits/YYYY-MM-DD-<scope>.md.
name: "ViralScan Scientific Audit"
argument-hint: "Optional: specific module to focus on (e.g. multimap, detection, umap, analysis, ncbi_fetch, build_reference)"
agent: "agent"
model: ['Claude Opus 4.7 (copilot)', 'Claude Opus 4.5 (copilot)', 'Claude Opus 4.1 (copilot)', 'Claude Sonnet 4.5 (copilot)']
tools: [codebase, search, searchResults, usages, findTestFiles, problems, testFailure, fetch, githubRepo, runCommands, terminalLastCommand, terminalSelection]
---

You are performing a **rigorous scientific audit** of **ViralScan** — a Snakemake-driven
Python CLI that quantifies viral load from paired-end FASTQ samples using
`kb-python` (kallisto + bustools).

## Goal

Determine whether ViralScan implements **scientifically valid and defensible methods**
for viral load quantification in single-cell or bulk RNA-seq data.

The primary analysis pipeline is:

```
FASTQ (paired-end)
  → kb count (kallisto + bustools)            [Snakefile: kb_count rule]
  → analysis.py   (GTF parsing, viral accession list)
  → multimap.py   (multimapping correction, proportional redistribution)
  → detection.py  (threshold-based viral detection + visualizations)
  → umap.py       (QC filtering, normalization, HVG, PCA, UMAP, clustering test)
```

Supporting paths: `ncbi_fetch.py` (NCBI download + caching), `build_reference.py`
(host + viral reference construction), `createconfig.py` (YAML config generation).

This is an **audit pass only** — do **not** propose, draft, or apply code changes.

## Operating Mode: Read-Only

- Do **not** edit, create, or delete any source, test, documentation, or config file.
- The **single permitted write** is the audit report itself (see *Output Destination*).
- Terminal usage is restricted to **read-only inspection**:
  `grep`, `cat`, `git log`, `git diff`,
  `PYTHONPATH=src python -m pytest tests/ -q` (diagnostic read only).
- Do **not** run formatters, `pip install`, commits, pushes, or anything that mutates
  the working tree.
- If a check requires modifying code to verify, **describe** it in the Verification Plan.

## Output Destination

Write the full audit report to:

```
audits/YYYY-MM-DD-<scope-slug>.md
```

- `YYYY-MM-DD` = today's date.
- `<scope-slug>` = short kebab-case label (e.g., `full-pipeline`, `multimap`, `detection-umap`).
- Create `audits/` if it does not exist. Append `-2`, `-3` suffix before overwriting.
- Print only a brief pointer + executive summary in chat; the full report lives in the file.

## Scope

If the user provided a focus module, restrict the audit to that script/path.
Otherwise, audit the full pipeline end-to-end as listed above.

## Rules

- **Do not assume.** If a scientific claim cannot be verified from code or tests, mark it **uncertain**.
- **Trace, don't guess.** Show input → transformation → output chains with exact file and
  line citations (e.g., `src/viralscan/scripts/multimap.py:145-162`).
- **Tests must validate science, not just execution.** Distinguish "runs without error"
  from "verifies the biological/statistical claim."
- **No optimization or refactor suggestions** in this pass.
- **Flag silent failure modes**: NA propagation, zero-count genes, sparse-to-dense
  conversion pitfalls, dropped barcodes, default config values that change biology,
  `try`/`except` swallowing errors, integer overflow in UMI counts.

## Audit Focus Areas

For each area, state what you checked, what you found, and what remains uncertain.

---

### 1. UMI Count Layer Logic (multimapping correction)

**Highest-priority check.** The multimapping correction in `multimap.py` produces two layers:
`counts_original` (unique-mapping UMIs) and `counts_corrected` (redistributed multimapper
share). Downstream scripts (`detection.py`, `umap.py`) combine them as:

```python
adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]
```

Verify:
- What exactly does `counts_corrected` store per EC?  
  Is it the **full redistributed count** (unique + share) or **only the additional share**?  
  Trace `multimap.py: build_multimap_matrix()` → what is written into each layer.
- Is the addition `counts_corrected + counts_original` correct or does it double-count
  uniquely-mapping reads?  
  (Known open issue in PLAN.md Task 1: umap.py may still have the bug;
   detection.py was allegedly fixed — verify both.)
- Does `multimap.py` write `counts_corrected = 0.0` for ECs that map to a single gene,
  ensuring the addition is safe for unique-mappers?
- Is proportional redistribution (share = count / n_genes_in_EC) documented and justified
  biologically? Is equal-weight redistribution the correct prior?

---

### 2. Viral Detection Threshold

- `detection.py: preprocessing()` applies `total_count >= threshold` (default 1 UMI).
  Is a threshold of ≥1 UMI across **all cells** a meaningful criterion for viral presence,
  or is it dominated by sequencing noise and ambient RNA?
- Is the threshold applied to the **corrected** or **raw** matrix?
- For the multimapping path: does the threshold apply after the layer addition
  (`adata.X`) or to the individual layers?
- Is there any false-discovery control for multi-virus comparisons? ViralScan tests
  against 195 viruses simultaneously — is there any multiplicity adjustment?

---

### 3. GTF Parsing and Accession Matching

- `analysis.py: obtain_gtf()` extracts `gene_id` from field 9 of each GTF line via
  `info.split('"')[1]`. Verify this is robust to GTF variants (no quotes, extra spaces,
  attribute ordering).
- Are all 195 bundled GTFs guaranteed to have `gene_id` as the **first** quoted attribute?
- User-supplied GTFs (via `config["gtf"]`) are parsed identically — is there any
  validation that user GTFs follow the same format?
- The viral accession list written to `log/analysis.txt` is later consumed by
  `detection.py` and `multimap.py`. Is there a race condition or ordering dependency
  in the Snakefile DAG that could cause stale reads?

---

### 4. Barcode Handling

- `multimap.py: load_barcodes()` strips `-1` suffix from barcodes.
  `normalize_barcodes()` also strips `-1` from the BUS DataFrame.
  Are these transformations consistent with how `adata.obs_names` are set by kb-python?
  Could mismatches between the two stripping operations cause silent barcode drops?
- Is there any deduplication of barcodes, or could duplicate barcodes (e.g., empty drops)
  propagate through the matrix?

---

### 5. Normalization and UMAP Computation

- `umap.py` applies `sc.pp.normalize_total(target_sum=1e4)` then `sc.pp.log1p()`.
  This is performed **after** viral count extraction (`viral_counts`, `viral_genes_expressed`).
  Verify whether viral counts used for per-cell labeling come from **raw** or **normalized**
  counts.
- HVG selection parameters (`min_mean=0.0125, max_mean=3, min_disp=0.5`) are hardcoded.
  For datasets where viral genes are sparse low-count features, are these parameters
  appropriate? Viral genes are force-included via `adata.var["highly_variable"] |= ...`
  in the virus-present branch — is this also done in the no-virus branch?
- PCA is run on `use_highly_variable=True` in both branches; verify consistent behavior.
- UMAP seeding: is `random_state` set for `sc.pp.neighbors`, `sc.tl.pca`, `sc.tl.umap`?
  Without fixed seeds, results are non-reproducible between runs.
- The viral neighbor enrichment test (`viral_neighbor_enrichment()`) uses a permutation
  p-value with a Laplace correction `(sum >= observed + 1) / (n_permutations + 1)`.
  Verify: is this one-sided or two-sided? Is the Laplace correction appropriate or does
  it mask a zero-count floor issue?

---

### 6. Multimapping Redistribution Algorithm

- Trace `multimap.py: build_multimap_matrix()` in full.
  What algorithm is used: equal-weight proportional redistribution?
  Is it EM-based or simple frequency weighting?
- Are multi-gene ECs that span **both viral and host** genes treated differently from
  virus-only multi-gene ECs? A read mapping equally to a host and a viral gene receiving
  50% viral weight could substantially inflate viral counts.
- Is there a cap on the redistributed count (e.g., non-negativity enforced)?
- Are the sparse matrix operations (COO → CSR conversion) numerically exact for
  UMI integer counts, or can floating-point redistribution introduce non-integer artifacts?

---

### 7. Reference Construction (`build_reference.py`, `ncbi_fetch.py`)

- `ncbi_fetch.py` downloads FASTA + GTF from NCBI. Verify that:
  - Retry/backoff logic handles HTTP 429 and transient failures without silently
    returning an empty/partial file.
  - Cache validation uses content hash (SHA256/MD5), not just file existence.
  - FASTA/GTF pairs are checked for consistency (same accession in both files).
- `build_reference.py` concatenates host (Ensembl) + viral (NCBI) FASTAs/GTFs.
  Verify that accession namespaces cannot collide (Ensembl transcript IDs vs NCBI
  accessions) and that the t2g mapping file correctly maps each transcript to its gene.

---

### 8. Config Validation and Silent Misconfigurations

- `createconfig.py` writes a YAML config consumed by all Snakemake scripts.
  Verify that all required keys are present and validated before Snakemake is invoked.
  Could a missing or wrong-type key cause a silent downstream error
  (e.g., `config.get("key")` returning `None` passed to `open()`)?
- Boolean config values: PLAN.md notes legacy `"True"`/`"False"` string booleans.
  Are all boolean keys now normalized to native Python bools? Check for any remaining
  `if config["key"] == "True":` patterns in scripts.
- Is `detection_threshold` validated to be a positive number?
  Is `min_counts` / `min_genes` validated to prevent filtering out all cells?

---

### 9. Edge Cases and Boundary Conditions

- **No viral genes detected**: does each script handle `found_genes = {}` gracefully
  without crashing on downstream dict/array operations?
- **Single cell**: does the neighbor enrichment test handle `n_obs = 1` (k-NN undefined)?
- **All cells viral-positive**: does the permutation test produce a meaningful result
  when `labels` is all-ones?
- **Zero-count viral gene**: a gene in the viral accession list present in `adata.var`
  but with all-zero counts — does it pass the threshold check silently?
- **Large sparse matrices**: does converting a sparse matrix to dense
  (`toarray()`) for a dataset with >100k cells × 30k genes cause memory exhaustion
  without a guard?

---

### 10. Reproducibility

- Is `random_state` / `np.random.seed` set before:
  `sc.pp.subsample()`, `sc.pp.neighbors()`, `sc.tl.pca()`, `sc.tl.umap()`,
  and `np.random.permutation()` in `viral_neighbor_enrichment()`?
- Does the Snakemake workflow produce bit-identical results on re-run with the same
  inputs and config (given fixed seeds)?
- Are kb-python and kallisto versions pinned in `environment.yml` / `pyproject.toml`?
  kallisto minor versions can change pseudoalignment results.

---

### 11. Numerical Stability

- Sparse matrix layer addition (`counts_corrected + counts_original`): are both layers
  guaranteed to be in the same format (CSR vs CSC vs COO)? Mixed formats can silently
  produce dense intermediates.
- `log1p` applied after `normalize_total`: are there any zeros that survive
  normalization that could cause issues downstream? (They should be fine — `log1p(0)=0`.)
- Division in proportional redistribution: is division by the number of genes in an EC
  guarded against zero (impossible by construction, but verify)?

---

### 12. Test Coverage of Scientific Claims

Evaluate the existing test suite (`tests/`) for scientific correctness:
- Do any tests verify that the multimapping correction output sums to the correct
  total UMI count (conservation of counts)?
- Do any tests verify that a known viral accession in synthetic data is correctly
  detected at the right count?
- Do any tests verify that the detection threshold correctly excludes a gene with 0 UMIs?
- Do any tests seed randomness and verify UMAP/PCA reproducibility?
- Are integration tests present (`tests/integration/`) that run the full pipeline
  on minimal synthetic FASTQ data?

---

## Output Format

Produce the report in this exact structure:

### 1. Executive Summary
3–8 bullets: headline verdict on scientific soundness, plus most consequential findings.

### 2. High-Risk Scientific Issues (Definite Errors)
Issues demonstrable as wrong from the code itself.

For each:
- **Finding** — one sentence.
- **Location** — `path/to/file.py:Lstart-Lend`, function name.
- **Evidence** — code snippet or test result proving it.
- **Consequence** — how user results could be biased or scientifically wrong.

### 3. Medium-Risk Issues (Likely Problems)
Issues that look wrong or fragile but depend on context.

Same fields, plus **What would resolve the uncertainty**.

### 4. Unclear Assumptions Requiring Domain Review
Methodological choices defensible *or* indefensible depending on experimental context.
Phrase each as a question for a bioinformatics domain expert.

### 5. Verification Plan
For every item in sections 2–4, give a concrete, runnable check:
- Simulation with known viral input and expected count recovery.
- Conservation-of-counts property test (total UMI before vs. after multimapping correction).
- Comparison to a reference implementation (e.g., kallisto's own EM correction).
- `pytest` test stubs that would pin down the behavior.
- Specific parameter combinations that expose the edge case.

---

## Verification Standard

- For each major analysis path, include a **trace block** showing the call chain with
  file:line citations (e.g., `detection.py:preprocessing() → adata.X = layers["counts_corrected"] + layers["counts_original"]`).
- Explicitly state when a claim rests on documentation vs. verified code vs.
  tests that only check execution (not scientific correctness).
- If you cannot trace a path because of missing context, list the files you need to
  read — do not fabricate a verdict.

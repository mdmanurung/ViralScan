# Review — EBV host-response analysis (hostresponse_ebv_matched) — 2026-07-01

**Scope**: The EBV host-response classifier analysis — method (`src/viralscan/scripts/hostresponse.py`),
driver (`scripts/hostresponse_ebv_matched.py`), this session's MCC addition, and the
comprehensive report (`analysis/hostresponse_ebv_matched/reports/`).
**Files reviewed**: 6 (+ context)
**Sub-agents run**: 6 (stats-causal, data-pipeline-leakage, bioinformatics, llm-failure-modes, doc-schema-fidelity, code-quality)

## Key decisions in this analysis

- **HVG feature selection before the train/test split** — `sc.pp.highly_variable_genes`
  runs on all 1906 cells, then cells are split per seed. See **F1**.
- **EBV-positive label = ≥10 EBV UMI** — single operational threshold, no independent
  ground truth; specificity gap interpreted as low-level EBV activity. (informational)
- **Uncertainty = SD across 6 fixed seeds** of an 80/20 balanced, depth-filtered split
  — not a k-fold CV; test set ≈146 cells/seed, not 1906. See **F4**.
- **`detection_threshold` default differs** between module (1) and driver (10). See **F2**.
- **Stability selection seed** — code uses `seeds[0]` (=0); report says 42. See **F3**.
- **Baseline is chance only** (AUC 0.5 / MCC 0); GSE158275 is data source, not a numeric
  comparator. (informational; correct framing)

## Questions for the analyst

1. Is this report a **paper figure** (headline numbers quoted downstream) or an internal
   sanity check? If the former, **F1** (HVG leakage) should be fixed before the numbers
   circulate, since it likely inflates AUC/MCC.
2. Is within-sample classification the intended estimand (cells as the unit), or do you
   want a statement that generalizes across LCLs? The current design supports only the former.
3. Do you want the reported ± to represent finite-sample uncertainty (→ bootstrap/k-fold)
   or just split-shuffle stability (current)? This changes how the SD should be described.
4. Is the ≥10-UMI positive label fixed by design, or open to the sensitivity sweep already
   in the todo?

## Findings

### Data pipeline & leakage
#### Major
##### F1. Highly-variable-gene selection is fit on all cells before the train/test split
`src/viralscan/scripts/hostresponse.py:405-406` (and `_select_features` 120-127)
```python
_detect_and_normalize(host_adata)                        # all 1906 cells
X, feature_names = _select_features(host_adata, use_hvg)  # HVG on the full matrix
weights_df, metrics = _run_l2_regression(X, virus_presence, depth, seeds, feature_names)
# _run_l2_regression then calls _balanced_split — the split happens AFTER HVG selection
```
**Why it matters here**: EBV+ cells are 62% of the data and drive variance in EBV-responsive
genes, so the HVG set is biased toward the most discriminative features — selected using cells
that later appear in the test set. The headline **AUC 0.845 / MCC 0.539 are upward-biased by an
unknown amount**. The report's "evaluated on held-out cells" is true for rows but not for the
feature columns. (Confidence: high.)
**Fix**: Move HVG selection inside the per-seed loop — select on `X_train` only, apply the mask
to `X_test` — and re-run. Expect the headline metrics to drop somewhat.

#### Minor
##### F6. Stability selection runs on all cells (no held-out partition)
`src/viralscan/scripts/hostresponse.py:441-443`
```python
stab_probs = _run_stability_selection(X, virus_presence, n_stab_iter, seed=seeds[0] if seeds else 42)
```
**Why it matters here**: The 15-gene "stable" set is a full-data estimate, so "stably selected"
is slightly optimistic. Does **not** affect AUC/MCC (those use per-seed held-out cells).
**Fix**: Run on a fixed training partition, or disclaim as a full-data estimate.

### LLM coding antipatterns
#### Major
##### F2. `run_hostresponse` defaults `detection_threshold=1`, driver uses `10`
`src/viralscan/scripts/hostresponse.py:353` vs `scripts/hostresponse_ebv_matched.py:215`
```python
# module:  detection_threshold: int = 1,
# driver:  detection_threshold: int = 10,   # always passed explicitly
```
**Why it matters here**: The published 1179/727 split and all metrics use 10 (the driver passes it
explicitly, so the report is correct). But any caller invoking `run_hostresponse` as a library
function or via the module CLI silently gets a **10× looser positive call** and a different analysis.
(Confidence: high.)
**Fix**: Align the module default to 10 (or make it a required argument) and note it in decisions.md.

#### Minor
##### F8. `_raw_depth` computed before the `log1p`-in-uns early return (latent)
`src/viralscan/scripts/hostresponse.py:87-95`
```python
host_adata.obs["_raw_depth"] = raw_depth
if "log1p" in host_adata.uns:      # pre-normalized inputs return here
    return                          # _raw_depth would then hold log-sum, not UMI sum
```
**Why it matters here**: Does **not** affect this analysis (kb-python h5ads are raw counts, no
`uns["log1p"]`). Latent bug: a pre-normalized input would give a depth-ranking that's a biased
proxy, corrupting the top-50%-depth filter in `_balanced_split`.
**Fix**: Compute `_raw_depth` from a raw layer, or only in the normalization branch.

### Statistics & causal inference
#### Minor
##### F4. Seed-SD is presented as an "honest spread"; effective test N ≈146/seed, not 1906
`analysis/hostresponse_ebv_matched/reports/hostresponse_ebv_matched-report.tex` (Methods + Fig caption)
```latex
Performance was averaged over six random seeds to give an honest mean and spread ...
... on \SciVal{\NCellsMatched}{1906} matched cells.
```
**Why it matters here**: The 6 test sets overlap (same depth-filtered pool), so the SD captures
split-shuffle noise, not finite-sample uncertainty. Each seed evaluates ~146 cells (top-50% depth,
balanced 50/50, 20% held out) — the "on 1906 cells" caption implies precision the test set doesn't
have. Direction of the result stands; the precision framing is overstated.
**Fix**: State the per-seed effective test size, or move to stratified k-fold; drop "honest" or call
the ± "split-shuffle variability."

##### F5. MCC's stated rationale ("more conservative when classes differ in size") doesn't apply to the balanced test set
`...-report.tex` Methods §Definitions
```latex
It is more conservative than accuracy when the two classes differ in size, which is why we report it alongside AUC.
```
**Why it matters here**: `_balanced_split` makes the test set 50/50, so there is no imbalance to be
robust to; on a balanced test, MCC and balanced accuracy are near-redundant. The code comment
("robust to the class balance produced by `_balanced_split`") is accurate; the report's rationale
refers to the natural 62/38 split that is never tested.
**Fix**: Soften to "a single-number confusion-matrix summary reported alongside balanced accuracy,"
or compute MCC at natural prevalence from calibrated probabilities.

### Documentation & schema fidelity
#### Major
##### F3. Report says stability selection used "seed 42"; the code uses `seeds[0]` = 0
`analysis/hostresponse_ebv_matched/reports/hostresponse_ebv_matched-report.tex:119-121`
```latex
Randomized-lasso stability selection (seed 42) scored each host gene ...
```
**Why it matters here**: `_run_stability_selection(..., seed=seeds[0] if seeds else 42)` with
`DEFAULT_SEEDS[0]=0` means the actual seed is **0**; 42 is the unreachable fallback. A reader
reproducing the 15-gene set with "seed 42" would get a different list. The "42" was inferred from
the parameter default, not the call site. (Confidence: high.)
**Fix**: Change the report text to "seed 0."

#### Minor
##### F7. Provenance line numbers are wrong/untraceable
`analysis/hostresponse_ebv_matched/outputs/numbers.json` and `reports/.manifest.json`
```
numbers.json:  "computed_at": ".../<stdin>:L3"     # from an interactive heredoc
manifest.json: mcc_mean and mcc_sd both -> hostresponse.py:L218  # per-seed append, not aggregation
```
**Why it matters here**: Defeats automated freshness/drift auditing — the whole point of the
registered-value provenance. Values themselves are correct.
**Fix**: Re-register from a script path (real line numbers); point mcc mean/sd at the aggregation
(≈L239-243).

### Bioinformatics
#### (no Major/Minor beyond F8 above; positive verifications below)

### Code quality
#### (no findings — the MCC addition is internally consistent)

## What was checked but is fine
- **Double dipping**: host feature matrix verified to contain only `ENSG…` human IDs — no viral
  accession leaked into features; EBV labels and host features come from separate h5ads. Clean.
- **Numbers vs source**: all 11 report numbers (AUC 0.845, MCC 0.539±0.059, sens 0.824, spec 0.710,
  bal-acc 0.767, n=1906, 1179/727, 15 genes, threshold 10) match `hostresponse_summary.tsv` /
  `numbers.json` exactly at stated precision.
- **Within-seed split**: train/test genuinely disjoint (`replace=False`, sequential split).
- **Normalization**: `normalize_total`+`log1p` are per-cell — not leakage.
- **Barcode matching**: `-1` suffix stripped consistently on both sides; dedup handled.
- **sklearn API**: `matthews_corrcoef(y_true, y_pred)` correct; returns 0 (not ValueError) on
  degenerate splits, so the MCC line needs no try/except.
- **Pseudoreplication / single sample**: acknowledged in the report caveats.
- **AUC try/except asymmetry**: unreachable given `_balanced_split` guarantees; fine.

## Resolution (2026-07-01, same day)

- **F1 (HVG leakage) — FIXED.** HVG selection moved inside the CV fold (`_hvg_mask`
  on training cells only). Re-ran leakage-free. **Surprise:** corrected metrics went
  *up*, not down — AUC 0.845→0.866, MCC 0.539→0.570, specificity 0.710→0.744. The
  global HVG was majority-class-dominated; per-fold HVG on the balanced training set
  gives better EBV-contrast features. Verified (per-fold path ran, 15 stable genes
  identical, no test leakage). Report regenerated.
- **F2 (threshold default) — FIXED.** Module default `detection_threshold` 1→10.
- **F3 (seed 42 vs 0) — FIXED.** Report now says "seed 0."
- **F5 (MCC rationale) — FIXED.** Reworded for the balanced test set.
- **F7 (provenance) — FIXED.** `numbers.json` re-registered via a committed script
  (`analysis/hostresponse_ebv_matched/scripts/register_report_values.py`); manifest
  mcc line numbers corrected.
- **F4 (effective N) — ADDRESSED.** Methods now states ~150 cells/seed and frames the
  ± as split-shuffle variability.
- **F6 (stability on all cells) / F8 (`_raw_depth` latent bug) — documented minors,
  left as-is** (don't affect the headline; F6 is descriptive, F8 doesn't fire on
  kb-python inputs).
- 30 module tests pass after the change.

## Notes
- **F1 is the load-bearing finding.** It and F3 both touch the just-delivered report: F1 (if fixed)
  changes the headline numbers; F3 is a factual error in the report text. Recommend fixing F3 (cheap,
  unambiguous) now and deciding on F1 explicitly.
- F1 + F6 share one remediation shape (move feature/gene selection inside the CV split).
- No security, no fabricated numbers, no hallucinated APIs.

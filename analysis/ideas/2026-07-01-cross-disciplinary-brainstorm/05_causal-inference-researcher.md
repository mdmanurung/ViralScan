# Sequencing-Depth as Confounder: Bounding the Causal EBV→Host Effect via Sensitivity Analysis

## Persona
**Causal Inference Researcher** — identifying and bounding backdoor paths through measured confounders before any EBV→host claim.

## Motivation
The most credible rival to "EBV drives host genes" is **sequencing depth**: a cell with more total UMIs has more EBV UMIs AND more host counts, generating association through a common cause (cell size/library complexity). Until depth is blocked, the 15 genes could be confounding artifacts. `_raw_depth` is already stored, making this tractable.

## Connection to Existing Data
`_raw_depth`, EBV UMI, per-gene host counts in the same matrix; the 15 genes + importances; n=1906 (1179/727) supports stratified analysis.

## Approach
1. Draw the DAG: depth is a fork into EBV_UMI and every host gene (+ ambient RNA, doublets unmeasured).
2. Re-run the classifier adding log(_raw_depth) and log(total_host_UMI) as covariates; AUC/MCC before vs after.
3. E-values (Ding & VanderWeele) per gene: minimum confounder strength to explain away the OR; E<1.5 fragile, E>3 robust.
4. Depth-quintile stratified DE: associations surviving within strata weaken the depth story.

## Expected Insights
A number for how much of AUC 0.866 survives depth adjustment (turns the "association not causation" caveat quantitative); a depth-robust short list of genes; whether the sensitivity/specificity asymmetry is itself a depth artifact.

## Feasibility
- **Effort**: Low | **Data ready**: Yes | **Methods**: Standard (scikit-learn, statsmodels, scanpy)
- **Key risk**: Verify `_raw_depth` is host-only vs total (host+viral) — the wrong denominator reintroduces the confound.

---

# Reference Strategy as Instrument: Aligner Choice as an IV

## Persona
**Causal Inference Researcher** — exploiting the designed benchmark as a natural experiment.

## Motivation
Aligner choice (kallisto vs STARsolo) changes measured EBV UMI but has no direct biological effect on cell state — plausibly satisfying IV conditions (relevance, exclusion, independence). It can instrument "measured EBV load" to estimate the causal EBV→host effect, separating measurement noise from biology.

## Connection to Existing Data
`results/reference_strategy_benchmark.tsv` (per-run/strategy/aligner); the HSV-1 denominator artifact is exactly a measurement confound IV handles; P22.10 matched barcodes give the first-stage data; HHV-6A/6B ambiguity is a falsification-test natural experiment.

## Approach
1. Formalize: Z=aligner, D=measured EBV UMI, Y=host gene/classifier score, U=cell size/ambient/true intensity; IV → LATE for "complier" cells.
2. First stage: regress EBV UMI on aligner (fixed effects for strategy/run); F>10 check.
3. 2SLS per gene; compare to OLS (2SLS≫OLS → confounding down-biased naïve; collapse → up-biased).
4. Over-identification: add reference strategy as a 2nd instrument, Sargan–Hansen test; inconsistency implicates the HSV-1 artifact.

## Expected Insights
A novel aligner-as-instrument design in single-cell virology; a defensible causal estimate (or bound) robust to depth confounding; the HSV-1 artifact operationalized as a formal test; possibly a reframing from "EBV drives host" to "depth-associated host features predict EBV."

## Feasibility
- **Effort**: High | **Data ready**: Mostly (verify per-cell barcode-resolution UMI across strategies, not just aggregate) | **Methods**: Standard (statsmodels IV2SLS, linearmodels, AER::ivreg)
- **Key risk**: Exclusion restriction fails if aligners differ in host-gene quantification; test host-gene concordance across aligners in EBV− cells (instrument should be inert).

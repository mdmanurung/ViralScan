# Hill Coefficient Mapping of Single-Cell EBV Viral Burden

## Persona
**Pharmacologist** — dose-response fitting, Hill coefficients, therapeutic-window ID, with viral burden as a continuous "dose."

## Motivation
The pipeline binarizes EBV at 10 UMI. But is the host response a step at 10 UMI or a sigmoid with a definite Hill coefficient? n>1 → steep, cooperative (threshold meaningful); n≈1 → hyperbolic (binary discards gradient); n<1 → sub-linear saturation. Fitting Hill per gene gives each an EC50 and Hill coefficient — a pharmacological fingerprint of the host response.

## Connection to Existing Data
Per-cell EBV UMI (n=1906, continuous) = dose; 15 stable genes' normalized expression = effect. The 10-UMI threshold can be checked against the fitted EC50s. The P22.5 GMM infrastructure shows the project already models mixtures.

## Approach
1. Extract EBV UMI (log1p, dose) + 15-gene normalized expression (effect).
2. Per gene, fit 4-parameter Hill via scipy.curve_fit (bounded); bootstrap (500) CIs.
3. Rank by Hill coefficient: n>2 switch-like, n≈1 rheostat; overlay curves with the 10-UMI line.
4. Define a data-driven EC50-consensus threshold; compare to 10 UMI; re-evaluate classifier at the empirical EC50.

## Expected Insights
Cooperative vs graded host sensing; validation (or not) of the 10-UMI cutoff — possibly lifting MCC above 0.570; the EC10–EC90 "therapeutic window" (dynamic range for burden→transcription).

## Feasibility
- **Effort**: Low | **Data ready**: Yes | **Methods**: Standard (scipy Hill, bootstrap, matplotlib)
- **Key risk**: Sparse 1–9 UMI cells may under-fill the sub-EC50 regime, making n noisy for low-EC50 genes.

---

# Reference-Strategy Selectivity Index: Off-Target Profiling for Viral Quantification

## Persona
**Pharmacologist** — selectivity index, off-target activity, therapeutic index applied to aligner/reference "formulations."

## Motivation
Each (aligner × reference) combo is a "formulation." Target = accurate EBV UMI. Off-targets = spurious HHV-6B (contig-ambiguous) and HSV-1 (artifact) UMIs. Compute Selectivity Index SI = on-target ÷ off-target per strategy, rank, and find the "therapeutic window" in reference-strategy space.

## Connection to Existing Data
`results/reference_strategy_benchmark.tsv` (18 conditions) has per-strategy per-virus UMI + positive fractions. EBV = on-target; HHV-6B/HSV-1 = off-target (known ambiguity/artifact); the 2042-anellovirus panel = a non-specific "background binding" floor.

## Approach
1. Per condition: on-target = mean EBV UMI/EBV+ cell; off-target = Σ(HHV-6B+HSV-1 UMI)/cell; SI = ratio; average over the 3 runs.
2. Selectivity landscape scatter (potency vs off-target) with iso-SI hyperbolas at SI=1/10/100.
3. Reference complexity (human_only=1…combined=3) as an ordered "dose" vs SI — monotonic or a non-monotonic optimum?
4. Anellovirus detection fraction as the non-specific floor → the assay's limit of selectivity.

## Expected Insights
A ranked, quantitative reference recommendation replacing qualitative discussion; the cost of the combined reference's off-target load; possibly a targeted (human+EBV) reference beating the inclusive one; an intrinsic noise floor for future low-abundance-virus claims.

## Feasibility
- **Effort**: Low | **Data ready**: Mostly (one anellovirus aggregation step) | **Methods**: Standard (pandas/matplotlib/scipy)
- **Key risk**: Only 3 runs → SI variance may swamp strategy differences; use bootstrap CIs and report overlap.

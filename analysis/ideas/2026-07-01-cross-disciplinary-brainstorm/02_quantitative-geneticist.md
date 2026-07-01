# Viral Load Variance Components: Decomposing EBV Transcript Heterogeneity

## Persona
**Quantitative Geneticist** — partitioning phenotypic variance (V_P = V_A + V_D + V_I + V_E) into additive/epistatic/environmental parts.

## Motivation
In a clonal LCL, V_A ≈ 0 (shared genome), so per-cell EBV-load variance is transcriptional/stochastic. Do the 15 stable host genes act **additively** (independent dosage) or **epistatically** (some combinations far more predictive than their marginals)? The logistic classifier conflates these; decomposing them distinguishes a linear-dosage latency model from a switch-like epistatic circuit.

## Connection to Existing Data
15 stability-selected genes with per-cell expression (n=1906, 6 seeds); the AUC-0.866/MCC-0.570 gap is where interaction terms live; continuous per-cell EBV UMI is the quantitative phenotype.

## Approach
1. Additive baseline: marginal + joint L1-logistic variance-explained per gene; flag super/sub-additive pairs.
2. Pairwise epistasis: all 105 gene-pair product terms, LRT vs additive, BH-corrected; epistatic variance fraction.
3. Variance decomposition: LMM on log(EBV UMI+1) with a transcriptional "kinship" K = XXᵀ÷p → additive transcriptional heritability vs stochastic residual.
4. GxE analog: gene × reference-strategy interaction — are the 15-gene effects stable across human_only/all_virus/combined?

## Expected Insights
Whether the panel is a robust linear sensor (portable score) or a multiplicative circuit; a heritability upper bound explaining the AUC plateau; whether GxE exposes reference-alignment artifacts.

## Feasibility
- **Effort**: Medium | **Data ready**: Yes | **Methods**: Standard (scikit-learn, statsmodels, limix)
- **Key risk**: n=1906 with 105 interactions underpowered for small epistasis; "transcriptional heritability" is analogy, not literal GREML.

---

# Cross-Viral Additive Score: Does Co-infection Load Sum Additively or Saturate?

## Persona
**Quantitative Geneticist** — multi-locus (multi-virus) phenotype modeling: additivity vs dominance vs epistasis across viral "loci."

## Motivation
Each virus is a "locus" with a per-cell UMI dosage; the host transcriptome is the phenotype. Does an EBV+HHV-6B co-infected cell show the sum of individual effects (additivity), one virus dominating, or a qualitatively new response (viral epistasis)?

## Connection to Existing Data
The 3×3×2 benchmark gives per-cell EBV and HHV-6B UMI; P22.10 confirms barcode concordance; the 15 EBV genes are a ready response signature. (Exclude HSV-1 — denominator artifact.)

## Approach
1. Stratify into EBV±/HHV6B± quadrants; tabulate counts (double-positive count = power constraint).
2. Additive prediction from single-virus strata vs observed double-positive expression.
3. Two-way ANOVA per gene with EBV×HHV6B interaction, FDR-corrected; interaction variance fraction.
4. 2D dose-response surface (EBV UMI × HHV-6B UMI); additive vs multiplicative fit by AIC.

## Expected Insights
Additivity → per-virus scores are separable/summable; epistasis → genes at the intersection of two viral programs; even a null validates that per-virus UMI separates real signals.

## Feasibility
- **Effort**: Low | **Data ready**: Mostly (merge per-cell EBV+HHV-6B UMI onto the gene matrix by barcode) | **Methods**: Standard (scipy/statsmodels)
- **Key risk**: Double-positive cells may be too few in one sample; HHV-6B contig ambiguity adds noise. Run as a power check first.

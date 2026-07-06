# Mutual Information Bottleneck for EBV Status Prediction

## Persona
**Information Theorist** — compressing host transcriptomes through an information bottleneck to find the minimal sufficient statistic for EBV status.

## Motivation
"Which genes survive selection" ≠ "how many bits of EBV information host expression carries, and how efficiently we extract them." The IB (Tishby) traces the compression–relevance frontier: how many bits of the 15-gene signature actually distinguish EBV±, and is there a lower-dimensional representation with the same I(Z;Y)?

## Connection to Existing Data
n=1906 (1179/727) labeled at ≥10 UMI; 15 HVGs = X, EBV status = Y; I(host;EBV) is estimable nonparametrically on the existing matrix; AUC 0.866/MCC 0.570 is the performance ceiling; documented non-Gaussian (GMM) marginals justify nonparametric MI.

## Approach
1. Estimate I(X;Y) via KSG (mutual_info_classif proxy, then NPEET) → upper bound on predictive info.
2. Sweep IB β with a variational IB / β-VAE (1–8-bit bottleneck); record I(Z;Y), I(X;Z); plot the IB curve.
3. Find the "kink" (phase transition) = minimum sufficient statistic dimension.
4. Shapley-MI per gene vs the current 15-gene list — same genes, or does IB surface variance/correlation-missed genes?

## Expected Insights
A concrete bit-count ("EBV status = X bits; classifier captures Y bits, Z% efficient"); whether 3–5 genes suffice; whether the AUC gap is irreducible channel noise vs unused information; redundant vs synergistic encoding.

## Feasibility
- **Effort**: Medium | **Data ready**: Yes | **Methods**: Standard (scikit-learn MI, pytorch β-VAE, NPEET)
- **Key risk**: KSG degrades in high dimension (15-gene space fine; full HVG needs reduction); IB curve depends on the Z parameterization.

---

# Reference Strategy as a Noisy Channel: Capacity and Information Loss

## Persona
**Information Theorist** — each (aligner, reference) as a communication channel transmitting the true viral signal with noise/distortion.

## Motivation
The 18-condition benchmark (same cells, 3 refs × 2 aligners × 3 runs) is a channel-capacity experiment. C = max I(X;Y) quantifies how much each pipeline preserves. The HSV-1 denominator artifact (dynamic-range compression) and HHV-6A/6B ambiguity (symbol crosstalk) are concrete channel noise. "Which pipeline has higher capacity?" has a precise answer.

## Connection to Existing Data
`results/reference_strategy_benchmark.tsv` (per-cell per-virus UMI, all 18 conditions); same cells → consensus/STARsolo-combined as approximate ground truth; documented HSV-1 artifact = systematic bias; EBV ≥10-UMI labels (n=1906) as the transmitted message.

## Approach
1. True signal T (consensus) vs received R_c per condition; estimate I(T;R_c) over binned UMI (0,1–3,4–9,10–49,≥50).
2. Transition matrices P(R_c|T) as heatmaps; the HSV-1 artifact = asymmetric off-diagonal mass in human_only.
3. Decompose capacity loss: aligner (same ref), reference completeness (same aligner), interaction.
4. Model HHV-6A/6B as a binary symmetric channel with crossover p; capacity 1−H(p) bits = precise cost of the ambiguity per cell.

## Expected Insights
A ranked bits-of-fidelity table per pipeline (beats count-ratio comparisons); whether the HSV-1 artifact costs more bits than reference incompleteness; the exact bit-cost of HHV-6A/6B ambiguity; EBV-specific vs pan-viral optimal channels.

## Feasibility
- **Effort**: Low | **Data ready**: Yes | **Methods**: Standard (numpy/scipy discrete MI, matplotlib)
- **Key risk**: "True signal" is circular without external ground truth; run a sensitivity analysis using each pipeline in turn as the gold standard.

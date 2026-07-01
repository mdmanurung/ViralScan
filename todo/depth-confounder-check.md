# Depth-confounder check for the EBV host-response classifier

- **Priority**: high
- **Status**: open
- **Category**: validation
- **Date**: 2026-07-01
- **Author**: mdmanurung
- **Source**: ideation session 2026-07-01 (Causal Inference persona, idea 5a)

## Motivation
The headline result — host transcriptome predicts per-cell EBV status at AUC 0.866 /
MCC 0.570 — is explicitly labeled association, not causation. The most credible rival
explanation is **sequencing depth as a common cause**: a cell with more total UMIs has
more EBV UMIs *and* more host-gene counts, manufacturing an association through cell
size / library complexity rather than EBV biology. `_raw_depth` is already stored, so
this is cheap to test and would convert a verbal caveat into a number.

## Plan
1. Draw the DAG: depth (and unmeasured ambient RNA / doublets) as a fork into EBV UMI
   and every host-gene count.
2. Re-run the classifier adding `log(_raw_depth)` and `log(total_host_UMI)` as
   covariates; compare AUC/MCC before vs after. A large drop implicates depth.
3. Compute per-gene E-values (Ding & VanderWeele 2016): the minimum confounder strength
   needed to explain away each gene's odds ratio. Flag E<1.5 (fragile) vs E>3 (robust).
4. Depth-quintile stratified EBV± differential expression for the 15 genes — do
   associations survive within strata?

## Definition of done
- A number for how much of AUC 0.866 survives depth adjustment.
- A depth-robust subset of the 15 genes for biological follow-up.
- A note on whether the sensitivity/specificity asymmetry is itself a depth artifact.

## Key risk
Verify `_raw_depth` is host-only vs total (host+viral) depth — the wrong denominator
reintroduces the very confound it is meant to block.

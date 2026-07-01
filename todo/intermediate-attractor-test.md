# Intermediate-attractor test for the EBV specificity gap

- **Priority**: high
- **Status**: open
- **Category**: analysis
- **Date**: 2026-07-01
- **Author**: mdmanurung
- **Source**: ideation session 2026-07-01 (convergent across Evolutionary Biologist 1a,
  Stem Cell Biologist 4a/4b)

## Motivation
Four independent disciplinary lenses converged on the same hypothesis: the classifier's
**specificity gap (0.744 < 0.822 sensitivity)** may not be noise but a real, biologically
distinct **primed / intermediate cell state** between EBV-latent and EBV-lytic. The
misclassified cells (EBV+ called negative, and vice-versa) are the candidate intermediate
population. If they form a distinct host-transcriptome cluster, the binary EBV label is
too coarse and a 3-class (latent / primed / lytic) model is warranted.

## Plan
1. UMAP on the full host transcriptome (post-HVG, PCA) for all 1906 cells; overlay
   per-cell EBV UMI as a continuous gradient. HDBSCAN clusters without a preset count →
   two density peaks (bistable) vs three (intermediate).
2. Locate the misclassified cells (EBV+ classified −, EBV− classified +) in the embedding.
   Co-localized in a distinct region = intermediate-attractor evidence; scattered at the
   two-cluster boundary = classification noise.
3. If an intermediate cluster exists: DE latent vs intermediate vs lytic; check for a
   pre-lytic program (BZLF1 low, LMP1 high, BCL2/MCL1 up = survival-first).
4. Attractor depth = within-cluster ÷ between-cluster transcriptional variance → a stable
   state vs a transition.

## Definition of done
- A clear verdict: bistable (2 states) or tristable (intermediate exists).
- If tristable: intermediate-state marker genes + a recommendation to retrain with a
  3-class label.
- Either way it is publishable ("no intermediate" confirms true bistability).

## Key risk
n=1906 may be too small / heterogeneous to resolve a stable intermediate from noise; the
lytic population may be sparse. Pair with the depth-confounder check so a cluster isn't a
library-size artifact.

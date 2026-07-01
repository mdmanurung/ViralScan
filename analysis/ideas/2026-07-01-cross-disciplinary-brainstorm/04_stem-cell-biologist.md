# EBV Latency-to-Lytic Commitment Landscape: Pseudo-Viral-Load Ordering

## Persona
**Stem Cell Biologist** — per-cell EBV UMI as a latent→lytic commitment trajectory, like a differentiation pseudotime.

## Motivation
The EBV latent→lytic switch (BZLF1/Zta-driven, bistable, epigenetically remodeled) is analogous to lineage commitment. The continuous per-cell EBV UMI is a natural "viral pseudotime." Ordering cells along it reveals whether a "primed" intermediate host state precedes viral commitment.

## Connection to Existing Data
Continuous EBV UMI (n=1906); the 15 genes likely mix early-response (primed) and late (lytic) markers; the ≥10-UMI cut is arbitrary vs the distribution's true shape; SoRelle reported replication-state heterogeneity (bistability signature).

## Approach
1. GMM (1/2/3-component) on log EBV UMI + Hartigan dip test → bistable (bimodal) vs graded; valley = data-driven threshold.
2. Order cells by EBV UMI (or PC1 of the 15-gene space); Monocle3/PAGA trajectory anchored at near-zero cells.
3. GAM per gene along pseudotime → sigmoidal inflection = commitment-point genes; early graded = primed.
4. Epigenetic-memory overlay: if EBV+ LCL ATAC/H3K27ac exists in GEO, test whether primed genes sit near accessible BZLF1 enhancers; else correlate with published EBV chromatin maps.

## Expected Insights
Bistable vs stochastic-graded switch; a data-driven commitment threshold; a ranked list of "priming" genes that may gate reactivation; whether sensitivity>specificity arises because the classifier tags a primed intermediate.

## Feasibility
- **Effort**: Medium | **Data ready**: Mostly (chromatin overlay needs public data) | **Methods**: Standard (scikit-learn GMM, scanpy PAGA, pygam)
- **Key risk**: Few true-lytic cells → unreliable trajectory; fall back to 1D GMM+GAM.

---

# Host Transcriptional State-Space Geometry: A Detectable Intermediate Attractor?

## Persona
**Stem Cell Biologist** — Waddington-landscape reconstruction: two attractors (latent/lytic) or three (latent/primed/lytic)?

## Motivation
scRNA-seq revealed invisible metastable intermediate states in pluripotency/haematopoiesis — the cells that respond to perturbation. If EBV+ false-negatives occupy an intermediate attractor whose host program hasn't committed to lytic, that explains sensitivity>specificity and is testable in the host transcriptome alone.

## Connection to Existing Data
n=1906 with full host transcriptomes + EBV UMI labels; the ~211 EBV+ cells missed by the classifier and ~202 EBV− called positive are the candidate intermediate cells; the 15 genes span survival/activation/differentiation (SoRelle framing).

## Approach
1. UMAP on the full host transcriptome; overlay EBV UMI; HDBSCAN clusters (no preset count) → 2 vs 3 density peaks.
2. Locate misclassified cells in the embedding — co-localized (intermediate attractor) vs scattered at the boundary (noise)?
3. If an intermediate cluster exists, DE latent vs intermediate vs lytic; check for pre-lytic program (BZLF1 low, LMP1 high, BCL2/MCL1 up = survival-first).
4. Attractor depth = within/between-cluster variance ratio → true stable state vs transition.

## Expected Insights
Bistable vs tristable at the host level; whether to retrain with a 3-class (latent/intermediate/lytic) label; intermediate-state marker genes for mechanistic follow-up; an LCL reframing analogous to naive/formative/primed pluripotency.

## Feasibility
- **Effort**: Low | **Data ready**: Yes | **Methods**: Standard (scanpy UMAP+HDBSCAN, rank_genes_groups)
- **Key risk**: n may be too small/heterogeneous to resolve a stable intermediate; "no intermediate" is itself a valid finding (confirms bistability).

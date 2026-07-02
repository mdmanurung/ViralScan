---
topic: host-transcriptome-encodes-viral-infection-state
description: Whether a cell's own (host) gene expression carries enough signal to predict its intracellular viral infection status at single-cell resolution.
created: 2026-07-01
last_updated: 2026-07-01
status: active
---

# Host transcriptome encodes viral infection state

## F-001: A host-gene classifier predicts per-cell EBV status, but the headline is substantially depth-confounded
**Status:** contradicted (a major confounder overturns the naive interpretation; a weaker residual signal survives)
**Claim:** Within a single EBV-transformed LCL sample (SRR12682296, GSE158275), a cross-validated host-gene classifier separates EBV-positive from EBV-negative cells (EBV+ = ≥10 viral UMI) at AUC 0.866 / MCC 0.570 (n=1906; 1179+/727−; 15 stable genes). **However, a follow-up confounder analysis shows this is substantially driven by sequencing depth, not host biology:** the ≥10-raw-UMI label is strongly depth-dependent (EBV+ rate rises 27.5%→96.3% across host-depth quintiles); sequencing depth *alone* predicts EBV status at AUC 0.803 (all cells) and **0.967 within the classifier's own balanced+depth-filtered design — higher than the 0.866 host-gene model.** Because host features are depth-normalized, the classifier recovers a degraded proxy of the depth-driven label. Yet a residual genuine signal remains: **5 of the 15 stable genes retain a depth-adjusted association (E-value ≥2, p<1e-7)**, while the two top-stability genes are fully explained away by depth (p≈0.8 after adjustment).
**Implications:** The naive reading ("host transcriptome strongly predicts EBV status") is overstated — much of the apparent signal is a sequencing-depth artifact introduced by a raw-UMI positivity threshold. The pipeline's top-50%-depth + class-balancing step does NOT fix this (it matches class *counts*, not the between-class depth *distributions*; the filter actually widens the depth gap). A real but modest host-response signal exists in ~5 genes. Any downstream use (biomarker panel, mechanistic follow-up) should restrict to the depth-robust genes and re-define the label depth-independently.
**Tags:** virology, scrna-seq, host-response, ebv, classifier, single-cell, confounding, sequencing-depth

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-01 | mycelium-lifecycle-2026-07-01 | SRR12682296 (GSE158275, LCL) | ViralScan | AUC 0.866 / MCC 0.570 (leakage-corrected, per-fold HVG) | supports |
| 2026-07-01 | depth-confounder-check | SRR12682296 (GSE158275, LCL) | ViralScan | depth-alone AUC 0.803 (all)/0.967 (headline design); EBV+ rate 27%→96% by depth quintile; only 5/15 genes E≥2 | contradicts |
| 2026-07-01 | depth-matched-reanalysis | SRR12682296 (GSE158275, LCL) | ViralScan | depth-matched case-control (1006 cells, depth p=0.99, depth-alone AUC 0.48): host genes STILL predict at AUC 0.72 (15 genes) / 0.68 (5 robust genes) | refines |
| 2026-07-01 | cpm-label-crosscheck | SRR12682296 (GSE158275, LCL) | ViralScan | depth-normalized EBV-CPM label (corr with depth −0.09; depth-alone AUC 0.53): host genes predict at AUC 0.64 — 2nd independent de-confounding method agrees | refines |

**Resolved reading**: the naive headline (AUC 0.866) is depth-inflated, but a real,
depth-independent host-response signal survives — **AUC ~0.64–0.72** across two
independent de-confounding methods (depth-matched case-control; depth-normalized CPM
label), both well above their chance controls. The **5 depth-robust genes are
biologically coherent** (mapped via mygene.info): LTA (NF-κB/LMP1 target, up), ING1
(p53, up), MT-ND4L (OXPHOS, up; MT-QC caveat), EVI2B (down), MACROD2 (down) — aligning
with EBV-modulated host programs (SoRelle 2021). Not the strong readout headlined; not
an artifact.

### Open Questions
- Do the 5 depth-robust genes form a coherent biological program (pathway/GO)?
- Does the depth-matched AUC ~0.72 signal generalize across the other LCL lines / donors in GSE158275? (single-sample result)

## F-004: The host transcriptome is a continuum, not discrete EBV states, apart from a small lytic tail
**Status:** preliminary
**Claim:** Testing the (four-persona-convergent) hypothesis that the classifier's specificity gap reflects a discrete "primed/intermediate" host-cell state: it is not supported. Density-based clustering (HDBSCAN) leaves most cells unassigned and KMeans silhouette is low (~0.15–0.17 for k=2–4) both on all cells and in the depth-matched subset — the host transcriptome is largely a **continuum**, not discrete latent/primed/lytic attractors. A naive forced 3-way split of all cells stratifies mainly by **sequencing depth** (EBV+ rate 0.37→0.64→0.75 tracking median depth 12.6k→15.6k→27.7k), i.e. the apparent "intermediate" is the depth artifact. The one genuinely distinct, depth-independent population is a small **high-viral-burden (lytic-like) tail**: ~108 cells (5.7%) with EBV UMI >1000 and median viral fraction 0.20, forming a separate low-depth cluster in the depth-matched design.
**Implications:** EBV state in this LCL is continuous with a minority lytic fraction (consistent with SoRelle et al. 2021's heterogeneity framing), not a discrete tristable host-cell landscape. The "primed intermediate" reading of the specificity gap was appealing across four disciplinary lenses but is explained by depth. Any state-transition/pseudotime modeling should treat this as a continuum + rare lytic population, and single-cell "attractor" claims here need depth control first.
**Tags:** scrna-seq, ebv, cell-state, continuum, clustering, lytic, confounding

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-02 | intermediate-state-test | SRR12682296 (GSE158275, LCL) | ViralScan | HDBSCAN no clusters; silhouette ~0.15 (all + depth-matched); forced k=3 tracks depth; distinct high-burden lytic tail (5.7%, fraction 0.20) | supports |

### Open Questions
- Is the ~108-cell lytic tail transcriptionally distinct in the HOST program (BZLF1/lytic markers), or only in viral burden?
- Would a diffusion/pseudotime axis (rather than discrete clusters) reveal a graded host-response gradient after depth control?

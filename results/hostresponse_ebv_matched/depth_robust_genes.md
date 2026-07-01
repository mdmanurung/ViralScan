# The 5 depth-robust EBV host-response genes

Genes from the 15-gene stable set that retain a depth-adjusted association with EBV
status (E-value ≥ 2, p < 1e-7 in `depth_confounder_check.py`). Symbols mapped via
mygene.info. Direction = sign of the depth-adjusted odds ratio (per +1 SD expression);
OR > 1 = higher in EBV-positive cells.

| Ensembl | Symbol | Name | adj. OR | Direction (EBV+) | E-value |
|---------|--------|------|---------|------------------|---------|
| ENSG00000226979 | **LTA** | lymphotoxin alpha | 1.41 | up | 2.17 |
| ENSG00000153487 | **ING1** | inhibitor of growth family member 1 | 1.47 | up | 2.30 |
| ENSG00000212907 | **MT-ND4L** | mitochondrial NADH dehydrogenase subunit 4L | 1.35 | up | 2.04 |
| ENSG00000185862 | **EVI2B** | ecotropic viral integration site 2B | 0.62 | down | 2.61 |
| ENSG00000172264 | **MACROD2** | mono-ADP-ribosylhydrolase 2 | 0.65 | down | 2.44 |

## Biological coherence (narrative — n=5 is too small for formal GO enrichment)

The depth-robust signature is **not random noise**; it aligns with known EBV biology and
with the host programs SoRelle et al. (2021) reported (survival, activation, oxidative stress):

- **LTA (up)** — lymphotoxin-α is a TNF-superfamily cytokine and a canonical **NF-κB target**.
  EBV's LMP1 is a constitutive NF-κB activator, so LTA up-regulation in EBV-positive cells is
  the most mechanistically direct hit — an expected readout of active EBV latency-III signaling.
- **ING1 (up)** — a **p53-pathway** tumor suppressor driving growth arrest/apoptosis; consistent
  with an EBV-induced stress/survival response (EBV is known to modulate p53).
- **MT-ND4L (up)** — a mitochondrial **oxidative-phosphorylation** subunit; fits EBV-driven
  metabolic reprogramming and the oxidative-stress axis noted in the source study.
- **EVI2B (down)** — a myeloid/granulocyte-**differentiation** membrane gene; its down-regulation
  is consistent with EBV pushing proliferation over differentiation in transformed B cells.
- **MACROD2 (down)** — an ADP-ribosylhydrolase in the DNA-damage response.

## Caveats
- **MT-ND4L is a mitochondrial gene.** Mitochondrial-read fraction is also a scRNA-seq QC
  covariate (stressed/dying cells have high MT content), so this hit is the one most in need of
  a mito-fraction control before mechanistic claims. The other four are nuclear-encoded.
- n=5 genes: this is a narrative characterization, not a statistically-powered enrichment. A
  GO/pathway test would need the broader depth-robust gene set (re-selected on a depth-independent
  label), not just the top 5.
- Directions are from a single LCL sample; generalization untested.

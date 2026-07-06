# Cross-disciplinary brainstorm — ViralScan — 2026-07-01

**Personas** (7): Evolutionary Biologist, Quantitative Geneticist, Pharmacologist,
Stem Cell Biologist, Causal Inference Researcher, Ecologist, Information Theorist.
**Focus**: EBV host-response finding + reference-strategy benchmark, equally.
**14 ideas**, grounded in real project data (EBV classifier AUC 0.866 / MCC 0.570,
15 stable genes; 3×3×2 reference benchmark; 2042-accession anellovirus panel).

## Ideas by feasibility

### Low effort / data ready (start here)
| # | Persona | Idea | Data | One line |
|---|---------|------|------|----------|
| 1 | Evolutionary Biologist | Fitness landscapes of viral latency | Yes | GMM on per-cell EBV UMI → bistable attractors vs flat drift; is the specificity gap a real intermediate state? |
| 5 | Causal Inference | Depth-as-confounder E-value bounding | Yes | Re-run classifier with `log(_raw_depth)` covariate + E-values → how much of AUC 0.866 survives depth adjustment? |
| 4b | Stem Cell Biologist | Host state-space intermediate attractor | Yes | UMAP+HDBSCAN on host transcriptome; do misclassified cells form a 3rd (primed) attractor? |
| 3a | Pharmacologist | Hill-coefficient mapping of EBV burden | Yes | Fit dose-response per gene → is the 10-UMI threshold near the EC50? Cooperative vs graded response. |
| 3b | Pharmacologist | Reference-strategy Selectivity Index | Mostly | SI = on-target EBV ÷ off-target (HHV-6B+HSV-1) per strategy → ranked, quantitative reference recommendation. |
| 2b | Quantitative Geneticist | Cross-viral additive vs epistatic co-infection | Mostly | EBV±/HHV6B± quadrants + interaction test → do co-infection host effects sum or interact? |
| 6a | Ecologist | Per-cell virome alpha-diversity vs cell state | Mostly | Shannon/Simpson of the anellovirus community per cell vs EBV state → permissive vs marginal habitats. |
| 7b | Information Theorist | Reference strategy as a noisy channel | Yes | Bits of I(true; measured) per pipeline; exact bit-cost of the HHV-6A/6B ambiguity and HSV-1 artifact. |

### Medium effort
| # | Persona | Idea | Data | One line |
|---|---------|------|------|----------|
| 2a | Quantitative Geneticist | Viral-load variance components + epistasis | Yes | Decompose the 15-gene signal into additive vs epistatic; transcriptional "heritability" upper bound on AUC. |
| 4a | Stem Cell Biologist | Latent→lytic commitment pseudotime | Mostly | Order cells by EBV "pseudotime"; GAMs find commitment-point vs priming genes. |
| 6b | Ecologist | Keystone-virus co-occurrence network | Needs prep | Hub anelloviruses in a per-cell co-occurrence network → pioneers/succession, niche partitioning. |
| 7a | Information Theorist | Mutual-information bottleneck for EBV | Yes | How many bits does host expression carry about EBV; can 3–5 genes match the 15-gene signature? |

### High effort / needs new data
| # | Persona | Idea | Data | One line |
|---|---------|------|------|----------|
| 1b | Evolutionary Biologist | Anellovirus constraint mapping → reference compression | Needs prep | Per-site entropy across 2042 accessions → curated constrained-region reference vs the full panel. |
| 5b | Causal Inference | Aligner-as-instrument (IV / 2SLS) | Mostly | Use aligner choice as an instrument for measured EBV load → a defensible causal EBV→host estimate + Sargan test of the HSV-1 artifact. |

## Notable cross-cutting themes
- **The specificity gap (0.744 < 0.822) recurs** in 4 ideas (1, 4a, 4b, 5) as a possible real intermediate/primed state rather than noise — a convergent, testable hypothesis worth prioritizing.
- **The 10-UMI threshold is questioned** by 3 ideas (1, 3a, 4a) — a threshold sensitivity sweep (already a todo) would feed all three.
- **The reference benchmark reframes cleanly** as selectivity (3b), a noisy channel (7b), and an instrument (5b) — three lenses on the same 18-condition table.
- **The anellovirus panel is underused**: ecology (6a/6b) and evolution (1b) both open it up.

## Details
Per-persona idea write-ups: `01_*.md` … `07_*.md` in this directory.

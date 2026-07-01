# Fitness Landscapes of Viral Latency: Mapping EBV Load Distribution as a Selection Surface

## Persona
**Evolutionary Biologist** — fitness landscapes, constraints, and the geometry of viral persistence strategies

## Motivation
The bimodal EBV UMI distribution across 1906 LCL cells reads not as a measurement artifact but as the fitness landscape EBV navigates. The threshold split (1179/727 at ≥10 UMI) plus specificity (0.744) < sensitivity (0.822) suggests "EBV-negative" cells are a low-expression attractor rather than truly virus-free. Is low-level expression a bet-hedging strategy (robustness against immune surveillance) or neutral drift through a flat landscape region? The 15 stable host genes are the phenotypic correlates of the attractor the virus has settled into.

## Connection to Existing Data
The classifier (AUC 0.866, MCC 0.570) already names 15 discriminating genes; per-cell EBV UMI from SRR12682296 gives the viral-load axis. The false-negative stratum (EBV-UMI-low but classifier-positive) is the biologically interesting transitional population. `results/hostresponse_ebv_matched/` and the STARsolo concordance give orthogonal measurements.

## Approach
1. Fit a Gaussian mixture (2–4 components) to per-cell EBV UMI to test bimodal/trimodal vs heavy-tailed unimodal — the number of modes = number of attractors.
2. Project the 15 genes into PCA/UMAP, color by EBV load; discrete clusters = deep valleys, gradient = flat landscape.
3. Characterize the EBV-UMI-low-but-classifier-positive cells; a distinct intermediate cluster = transitional population.
4. Test bistability: per-component CV of EBV UMI vs a Poisson-at-observed-depth null.

## Expected Insights
Whether EBV heterogeneity reflects a bistable landscape (replicative vs latent attractors) or a continuous robustness landscape — predicting whether the specificity gap is noise or a real intermediate state.

## Feasibility
- **Effort**: Low
- **Data ready**: Yes
- **Methods available**: Standard tools (scikit-learn GMM, scanpy UMAP, scipy)
- **Key risk**: 727 low-UMI cells may be too few to resolve a rare intermediate mode; single-sample design can't separate cell-intrinsic bistability from copy-number heterogeneity at infection.

---

# Co-evolutionary Constraint Mapping Across the 2042-Accession Anellovirus Panel

## Persona
**Evolutionary Biologist** — neutral networks, constraint detection, and evolvability of viral detection signal

## Motivation
ViralScan bundles 2042 anellovirus accessions. Which genomic regions are purifying-selected (robust alignment anchors) vs hypervariable (unreliable, causing the HHV-6A/6B-style ambiguity already seen)? Detection sensitivity is a function of which regions the index includes — the tool's own "fitness landscape."

## Connection to Existing Data
The 2042 accessions live as GTF in `src/viralscan/data/`. The reference-strategy benchmark (3×3×2) already quantifies how reference composition affects EBV/HHV-6B/HSV-1 detection — the design needed to extend to anelloviruses. The HHV-6A/6B finding proves the reference-strategy effect is real.

## Approach
1. Build an MSA of a phylogenetically representative subset (~200 accessions by NJ on pairwise ANI); compute per-site entropy.
2. Partition genomes into constrained (<0.2 bits/site) vs hypervariable (>1.5 bits/site); measure k-mer (k=31) coverage per partition in the current index.
3. Run ViralScan on the 3 benchmark runs with (a) all-accession, (b) constrained-only, (c) one-representative-per-species indices; compare detection via the existing benchmark framework.
4. Build a constraint map: per-ORF conservation vs coverage trade-off; the Pareto front defines the evolutionarily optimal reference.

## Expected Insights
Which anellovirus regions are universal detection anchors, and whether the 2042-accession strategy is over-specified (diluting unique k-mers, worsening multi-mapping) — an evidence-based reference-compression recommendation.

## Feasibility
- **Effort**: High
- **Data ready**: Needs preprocessing (fetch FASTAs via `viralscan build-ref`/`ncbi_fetch.py`; 200-accession subset tractable)
- **Methods available**: Standard (MAFFT/minimap2, custom entropy, jellyfish); benchmark framework exists
- **Key risk**: Downloading/aligning 2042 FASTAs is heavy + NCBI-rate-limited; benchmark runs may lack anellovirus reads, needing an anellovirus-positive control.

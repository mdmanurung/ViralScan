# Per-Cell Virome Alpha-Diversity as a Cell-State Predictor

## Persona
**Ecologist** — diversity indices as ecological fingerprints of habitat quality.

## Motivation
Each cell is a habitat hosting a viral community; the anellovirus panel is a rich commensal community. Does host cell-state (EBV latency program, cycle phase) predict the viral alpha-diversity (Shannon/Simpson/richness) a cell hosts — permissive vs marginal habitats?

## Connection to Existing Data
Per-cell EBV UMI (SRR12682296) + the 2042-accession anellovirus panel = a barcode×accession OTU-table analog; the 15 host genes include latency markers (cell-state is latent in the feature matrix); P22.10 gives a joint host+viral matrix.

## Approach
1. Extract barcode×anellovirus UMI; filter ≥10 anellovirus UMI/cell; rarefy to equal depth (skbio).
2. Per-cell Shannon H, Simpson evenness, observed richness.
3. Assign state from the 15-gene classifier (EBV-high/low) + optional cell-cycle scoring.
4. Kruskal-Wallis / LM (H ~ EBV score + total UMI); violin + scatter; effect size, depth-adjusted.

## Expected Insights
Diversity ↑ with EBV load → permissive-habitat model; ↓ → competitive exclusion/ISG clearance; null → diversity is cell-extrinsic (inoculum-driven). Links host gene modules to viral community composition.

## Feasibility
- **Effort**: Low | **Data ready**: Mostly (confirm anellovirus output; classifier scores exist) | **Methods**: Standard (skbio, scanpy, scipy)
- **Key risk**: Anellovirus UMI may be too sparse per cell for rarefaction (median <5 → presence/absence richness only).

---

# Keystone Virus Detection via Co-Occurrence Network Hubs

## Persona
**Ecologist** — keystone-species identification through co-occurrence network topology.

## Motivation
A keystone is the most *connected*, not most abundant, species (sea otters, Pisaster). Are there anellovirus accessions whose presence predicts many others — pioneers facilitating a succession within a permissive EBV+ niche? Hub nodes in a per-cell virus–virus co-occurrence network are candidate keystones.

## Connection to Existing Data
Cell×anellovirus presence/absence from ViralScan (SRR12682296); the HSV-1 artifact/HHV-6B ambiguity inform which signals are trustworthy; EBV load + 15-gene scores let hub identity be correlated with host state.

## Approach
1. Presence/absence matrix (UMI>0), keep accessions in ≥1% of cells.
2. SparCC (fastspar) or phi-coefficients with 1000× label-permutation FDR.
3. Network (nodes=accessions, edges=significant co-occurrence); degree/betweenness/eigenvector centrality → top-5 hubs.
4. Validate: richness gap in hub-absent vs hub-present cells; hub taxonomy (TTV/TTMV/TTMDV); correlate with EBV load / classifier score. Flag negative (competition) edges.

## Expected Insights
Candidate keystone anelloviruses (biologically special vs phylogenetically central); succession signal (hubs enriched in EBV-high cells); modularity = niche partitioning; mutual-exclusion edges = intracellular competitive exclusion (rarely documented).

## Feasibility
- **Effort**: Medium | **Data ready**: Needs preprocessing (barcode×accession matrix; check 2042-accession sparsity) | **Methods**: Standard (fastspar, networkx)
- **Key risk**: Extreme single-cell sparsity makes co-occurrence unreliable; may need pseudo-bulk aggregation (losing the single-cell framing). Pilot: per-cell richness distribution first.

# Dedicated Viral scRNA-seq Tool Comparison (scaffold — not yet run)

**Status:** scaffold only. No comparator has been run; the manuscript does **not**
claim a dedicated-tool comparison (scope narrowed, see `docs/manuscript_draft.md`).

Scope when executed: compare ViralScan against one dedicated viral single-cell RNA-seq
detector (Venus, Luebbert et al. 2024 — already cited; or ViralTrack) on the
SRR12682296 EBV LCL matched barcodes.

Primary matched cell set: the 1,906 GSM4796271 LCL_777_B958 barcodes used in
`results/matched_barcode_comparison.tsv`.

Primary metrics:

- EBV-positive cells at >=1 viral UMI (or tool-equivalent positive call).
- EBV-high cells at >=10 viral UMI (or tool-equivalent high-confidence call).
- Lytic-marker-positive cells (BZLF1, BRLF1, BHRF1) where the tool exposes gene-level calls.
- Per-cell viral-burden rank concordance (Spearman) where both tools output continuous burden.
- Runtime and peak memory from scheduler accounting.

Acceptance rule: the manuscript may claim a dedicated-tool comparison only after a
`results/dedicated_tool_comparison.tsv` has non-empty rows for **both** ViralScan and the
comparator on the same cell set. Install the comparator in a separate conda environment
and write its raw outputs under `analysis/dedicated_tool_comparison/raw/` (do not overwrite
ViralScan benchmark artifacts).

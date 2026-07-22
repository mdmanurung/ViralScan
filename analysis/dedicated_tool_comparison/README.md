# Dedicated viral scRNA-seq comparison status

**Status: not yet run.** No dedicated v3 comparator result is eligible for a
claim. The authoritative tasks are `CMP-04` and `CMP-05` in
[`../../PLAN.md`](../../PLAN.md).

Scope when executed: compare ViralScan against **both** version-pinned Venus and
Viral-Track in isolated, unmodified environments. Run at least EBV, HHV-6B, and
HSV-1 with outcome-independent barcode anchors.

Primary matched cells will use the checksum-pinned, outcome-independent anchor
defined by the v3 validation protocol. No ViralScan-positive or legacy matched
cell set may define the comparison denominator.

Primary metrics:

- Cells with at least one unique EBV molecule on the shared feature universe.
- Cells with at least ten unique EBV molecules as a prespecified secondary burden tier.
- Lytic-marker-positive cells (BZLF1, BRLF1, BHRF1) where the tool exposes gene-level calls.
- Per-cell viral-burden rank concordance (Spearman) where both tools output continuous burden.
- Runtime and peak memory from scheduler accounting.

Acceptance rule: `results/dedicated_tool_comparison.tsv` has matched rows for
ViralScan, Venus, and Viral-Track, or an explicit reproducible failure record for
any planned row. Inputs, references, barcode anchors, feature intersections,
count layers, commands, tool versions, scheduler accounting, and unavoidable
differences must be audited. Raw outputs stay under
`analysis/dedicated_tool_comparison/raw/`; no comparator algorithm is silently
patched.

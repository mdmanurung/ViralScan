# ViralScan 3.0 molecule-counting contract

<!-- viralscan-claim:v3-molecule-conservation-synthetic status=validated_v3 -->
Input is the corrected, sorted BUS stream produced by the documented kb/bustools
stage. A molecule key is corrected `(CB, UMI)`. BUS `count` is PCR/read multiplicity
and is audit-only.

For each molecule, every observed EC is projected from transcripts to a set of
distinct genes. The intersection across its EC observations is the compatible gene
set. An empty intersection is an unresolved collision and receives no allocation.
One compatible gene is unique gene evidence, including multi-transcript ECs of the
same gene. More than one is ambiguous evidence.

Allocation conserves one unit: `equal` splits over genes; `host-conservative` sends
mixed host-virus mass only to compatible hosts and otherwise splits normally;
`unique-weighted` uses cell-local unique support plus a documented pseudocount;
`em-global` uses one sample model; `em-cell` uses cell-local unique support shrunk
toward the sample model. If a cell has no local unique support for any gene
compatible with its ambiguous molecules, `em-cell` retains equal-split ambiguity
instead of letting the positive sample prior force an assignment. Disjoint
molecules are unresolved and receive no allocation.

The opt-in `sibling-weighted` variant uses cell-local unique gene support plus
the configured pseudocount only when every compatible gene is viral and belongs
to the same non-empty `sibling_group` in the run's Virus Identity table. It
normalizes weights over those compatible genes only. Mixed host-virus sets,
cross-group sets and genes without a curated sibling group retain equal splits.
Absent identity metadata fails closed for this variant. This method does not
change the default `host-conservative` allocation or claim to identify an absent
sibling. No sample-wide prior or non-compatible gene receives molecule mass.

`counts_unique` contains one unit for each resolved singleton. The selected
method's ambiguous mass is in `counts_ambiguous_allocated` (`counts_corrected`
is its legacy alias); their sum is `X`. `counts_multimap_sibling_weighted` is
the diagnostic ambiguous-only layer for the opt-in variant. None of these
layers includes unresolved collisions, discarded barcode records or PCR
multiplicity.

The corrected-molecule host-virus boundary is a resolved CB-UMI whose final
intersection contains at least one host and one viral gene, before allocation.
The audit counts each such molecule once, independent of method and candidate
gene count. Per-gene `counts_host_viral_ambiguous` distributes one unit over
its compatible viral genes; it is a diagnostic upper mass, not the selected
allocation. This boundary cannot measure fragments discarded by STAR host
filtering: their viral compatibility is unobserved in the retained BUS stream.
SW-07's prefilter host-virus fragment boundary therefore remains open until a
joint host/viral prefilter observation and fragment identity contract exists.

Gene roles come from the same run Virus Identity table used for viral status
and sibling membership: host genes are `host`, curated viral roles are retained,
and an empty viral role becomes `unknown`. `index_kind` records `host_only`,
`combined` or `virus_only`; legacy runs infer this partition from analysis.txt
and label their metadata source explicitly. Counts per 10,000 quantified
molecules must name this indexed-gene universe. A virus-only matrix cannot
provide host depth or a host-molecule comparable-cell denominator; these require
the matched external host matrix/cell universe. The index kind describes
composition, not proof of biological absence or complete catalogue curation.

Invariants are:

`unique + ambiguous + unresolved = input molecules`

`sum(counts_ambiguous_allocated) = ambiguous molecules`

`X = counts_unique + counts_ambiguous_allocated`

All matrices must be finite and non-negative. Ordering and processing chunk size
must not change the result. The implementation streams sorted molecule records in
bounded buffers; golden/property tests enforce order and buffer-size invariance,
including the synthetic molecule-conservation claim registered for this contract.

`tests/test_multimapping_properties.py` checks all of the invariants above, for every
method, against an independent oracle that re-derives the contract from the raw
records (40 seeded random fixtures with off-list barcodes, UMI collisions and
unresolved molecules). It also checks that mass reaches only compatible genes, that
`sibling-weighted` moves mass only inside one sibling group, and that record order,
read multiplicity and streaming buffer size change nothing.

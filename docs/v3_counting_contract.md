# ViralScan 3.0 molecule-counting contract

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

Invariants are:

`unique + ambiguous + unresolved = input molecules`

`sum(counts_ambiguous_allocated) = ambiguous molecules`

`X = counts_unique + counts_ambiguous_allocated`

All matrices must be finite and non-negative. Ordering and processing chunk size
must not change the result. The implementation streams sorted molecule records in
bounded buffers; golden/property tests enforce order and buffer-size invariance,
and the retained EBV baseline verifies conservation on 103,145,071 BUS records.

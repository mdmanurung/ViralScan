---
id: F-014
date: 2026-09-27
topic: reference
tags: [anellovirus, low-complexity, k-mer, dustmasker, kallisto-dlist, CAT-17]
status: confirmed
---

# The low-complexity failure is a k-mer property, not an N-masking property

`CAT-17` was written as "mask the panel". Measuring the shipped artifacts showed
that framing is wrong in a way that matters, because the two properties are
independent and only one of them is what breaks detection.

## The measurement

Screened 2,000,000 R2 reads of the EBV LCL `SRR12682296` (EBV-immortalised LCLs,
10x v2; R1 26 bp barcode, R2 98 bp UMI+insert) against exact 31-mers from the
panel, then re-scored the captured hit reads against a low-complexity-masked
panel.

| Probe | Result |
|---|---|
| EBV `NC_007605.1`, 288,566 unique 31-mers | **29,207 reads = 1.46 % of R2** |
| TTV prototype `NC_002076.2`, 6,392 unique 31-mers | **0** |
| All 2,042 anelloviruses, no complexity filter | 7,175 / 500,000 = **1.44 %** |
| …same panel, low-complexity masked | **0 / 6,437** captured hit reads |

89.5 % of the captured hit reads carried a homopolymer run ≥ 31 bp (median
longest run 49, max 97). The k-mers that actually matched were literally
`A`*31, `AAAATAAAAAAAAAAAAAAAAAAAAAAAAAA`, `T`*31 and similar.

Method validation: exact recovery on synthetic controls at 1 %, 0.1 % and 0.02 %
abundance (1000/1000, 100/100, 20/20), zero false positives among ~100,000 random
reads.

## Why "was the panel masked?" is the wrong question

| Panel | N bases | pure-homopolymer 31-mers | all low-complexity 31-mers |
|---|---:|---:|---:|
| shipped `viral_genome.dedup.fa` (2,312 recs) | 785 / 9,925,822 = **0.0079 %** | **170** (9 recs; worst `NC_001479.1` ×85) | 7,236 / 9,853,675 = 0.0734 % |
| upstream hardmasked 2,023 reps | 284,762 / 5,927,006 = **4.80 %** | **0** | 3 / 5,337,346 = 0.0001 % |

Our panel is 99.99 % unmasked *and* contributes 170 pure-homopolymer 31-mers.
An N-content check reports it as essentially clean while the k-mer space is
still capable of manufacturing >1 % of reads. The upstream panel is 600× more
masked by N and 1,268× cleaner by k-mer. A record with a single 31-A run is
"0.003 % masked" by N and fully lethal by k-mer.

This is why the gate is an **absolute count**, not a fraction: 170 k-mers spread
over a 9.8 Mbp panel is a rounding error as a fraction and a 1.4 % false-read
floor in practice, because the reads are not rare even though the k-mers are.

## Consequence for the D-list defence

`CAT-19`/F-005 established the host genome D-list as the fix for host-homologous
reads. It does **not** help here. A D-list filters k-mers present in the host
sequences; these k-mers are present only in the viral panel, and are self-similar.
The two failure modes need two different defences, and the covid build shipped
the D-list defence without the masking one.

## What landed

`validate_reference_records` gained two opt-in gates (`max_pure_homopolymer_kmers`,
`max_low_complexity_fraction`), both default-off. `build_anellovirus_reference`
enforces them (strict when masking is requested). `scripts/build_bundled_panel_ref.py`
— the builder that made the shipped index and never called dustmasker — now
writes a viral-only `viral.fa` and gates before `kb ref`. The manifest records
`low_complexity_kmers` and `low_complexity_kmer_fraction` per sequence.

## Open

The covid `viral_summary.tsv` files still publish `Alphatorquevirus` 1,167,103
and 1,605,631 UMI, and `main` @ `4fcd748` closed F-005 as accession-level
artifact. The retraction is `CAT-30` and is blocked pending explicit user
approval.

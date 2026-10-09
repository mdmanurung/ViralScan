# F-028: The twostep HPV77 call in hpv16 SRR19537339 is human chr1 reads STAR did not align

**Date**: 2026-10-06
**Status**: result recorded (one sample; the mechanism is the leading explanation, not tested on a second sample)
**Found by**: `DSR-15` redetect, then `DSR-02` evidence (job 25720142, code `ef0c0c8`)

## Question

`twostep_v2` calls Human papillomavirus 77 in 263 of 10,270 reference cells (2.56 %, 280 molecules) in
`hpv16/SRR19537339`. `combined_off` and `combined_artefact` call only HPV29 (11 cells). Is HPV77 real?

## Result

No. `viralscan evidence --virus "Human papillomavirus 77" --blast` on the twostep_v2 run:

- All reads that aligned in the competitive step (GRCh38 plus panel) went to host; none went to HPV77.
  366 reads (313 molecules, 296 cells) sit in one window, `chr1:153607043-153607122`, at about 99.8 % identity.
- BLAST: all 353 reads with a hit matched human chr1 (mostly 100 % identity); no viral hit.
- Flags: `host_homology` flagged (1.0), `sibling_or_host_ambiguity` flagged (1.0).

## Why the arms differ

Twostep removes reads STAR aligns to the host, then maps the rest with kallisto against a virus-only index.
Reads from a human locus that STAR failed to align have no host target there, so the closest viral
sequence (HPV77) takes them. The combined arms keep host transcripts in the index and send them to host.

## Consequence

- The twostep row is a false positive; the combined arms are right for this sample.
- Any twostep call the combined arms do not make needs `viralscan evidence` before it is reported.
- Not yet checked: why STAR missed this locus (parameters `pinned`, `outFilterMismatchNmax 4`,
  `outFilterMatchNminOverLread 0.9`), and whether the other twostep calls are affected.

## Mechanism (added 2026-10-06)

- **TSO enrichment.** In the first 2 M reads of `SRR19537339` R2, 14.0 % carry the 10x TSO
  (`AAGCAGTGGTATCAACGCAGAGTACATGGG`) in the input, 49.2 % after the STAR host filter (about 3.5x).
  A 31-nt TSO prefix leaves under 90 % of the read alignable, so `outFilterMatchNminOverLread 0.9` rejects it
  and the read passes into the twostep quant.
- **Why HPV77.** All 384 target reads share one 31-mer with HPV77 (Y15175.1): a `GGGG(CAG)n` tract. All 384
  carry the CAG repeat and 199 (52 %) the TSO. The human locus is a CAG repeat.
- **Same call in two more samples** (twostep_v2 vs combined_off, cells): `hpv16/SRR19537341` 115 vs 0,
  `kshv_gse190558/GSM5725695` 38 vs 2. Evidence (jobs 25720309, 25720310, `ef0c0c8`):
  - SRR19537341: 167 reads, all HOST (chr1 167), `host_homology` flagged, no viral hit.
  - GSM5725695: 42 reads, all HOST (nine chromosomes, MAPQ 0, 100 % identity, `complex_body_fraction` 1.0),
    pure `CAGCAG…` repeat reads of about 55 nt. BLAST returned 0 reads (low-complexity masking), so
    `host_homology` and `low_complexity` both read `not_flagged`. The alignment table, not the flags, shows it is host.
- **Verdict rule consequence.** `dsr02_verdicts.py` must read `host_competitive_fraction` and
  `complex_body_fraction` from `alignment_qc.tsv`, and treat an empty BLAST table as "not assessed".
- **NCBI web BLAST** (12 reads, `blastn -remote -db nt`): pending at time of writing.

## D-list does not fix it (2026-10-07)

- All 6 shared 31-mers are in the genome D-list (`dlist.fa`), yet the read still pseudoaligns, on the cat42d
  and the genome+cDNA (`dlist2`) virus-only index alike (EC 3670,3671 = Y15175 E2/E4).
- One read x200 through `kallisto bus -x 0,0,16:0,16,28:1,0,0`: TSO alone 0/100, CAG alone 0/100, CAG+tail
  0/100, TSO+CAG (58 bp) 100/100. The junction `ACATGGGGCAGCAG…` k-mers carry the match.
- A kallisto D-list does not delete target k-mers that also occur in the decoy; a larger D-list is the wrong
  lever. Remaining levers: trim the TSO before the twostep quant, or require `viralscan evidence` on every
  twostep-only call (user decision, HANDOFF open items).

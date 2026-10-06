# F-026: The HHV-6A residual in the HHV-6B sample is HHV-6B reads, not HHV-6A

**Status:** confirmed (read-level). 2026-10-05. PLAN `DSR-05`.
**Sample:** SRR20710641 (HHV-6B, 10x 5′), combined run `ebv_latest_ref_2026-09-27/runs/combined/SRR20710641` (`-x 10xv3`, 3.0.0.dev0).
**Work dir:** `/exports/para-lipg-hpc/mdmanurung/ViralScan/hhv6a_in_hhv6b_2026-10-05/` (`run.sh`, Slurm job 25701668).

## Claim
HHV-6A (`NC_001664.4`) has 55 barcode-UMI records in 32 ECs: 10 records in HHV-6A-only ECs and
45 in ECs shared with HHV-6B (`NC_000898.1`); none involves a host transcript. Extracting the 95
read pairs (UMI exact, CB ≤ 1 mismatch) and aligning R2 (150 nt) with minimap2 against GRCh38 +
the panel, and against 6A + 6B alone:
- **Shared records (82 reads):** 52 align better to 6B, 25 tie, 1 favours 6A, 3 no hit, 1 host.
- **6A-only records (13 reads):** 5 align better to 6B, 3 tie, 4 favour 6A by one or two
  mismatches only (alignment score 290-300 vs 280-290), 1 short host hit (chr2, 44 nt).
- No read is telomeric (contrast F-025); no read carries a 6A-specific stretch with several
  mismatches against 6B.

## Interpretation
The EM share given to HHV-6A is bleed from shared regions (as `possible_em_bleed` already flags), and
the few 6A-only k-mer hits are single-base differences between this patient's HHV-6B strain and the
6B reference (Z29), not a 6A co-infection. There is no read-level evidence for HHV-6A in SRR20710641.

## Implications
- Keep the 6A residual as bleed; do not report HHV-6A co-detection.
- A 6A-only EC is not proof of 6A: 31-mer exclusivity breaks on single strain-level SNPs. A sibling
  call needs reads with several sibling-specific differences, which `DSR-02` should report per read.

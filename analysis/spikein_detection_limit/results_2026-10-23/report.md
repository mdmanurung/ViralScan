# Spike-in computational detection thresholds: results

Grid complete: 66 conditions scored from `scores_v2.tsv`. Rules and definitions: see `aggregate.py` docstring and README.

## Flank-index verdict: **RETIRE**

Rule: keep iff flank threshold <= old threshold in every background x layout AND flank adds <= 5 unsupported UMIs per null library (paired by replicate, max).

| background | verdict | max added unsupported UMIs (flank - old) | flank absolute max | notes |
|---|---|---|---|---|
| D03 | retire | 12.5714 | 12.5714 | c005: flank 10 <= old 10; c100: flank 10 <= old 10; flank adds up to 12.5714 unsupported UMIs per null library > 5 |
| R01 | retire | 9.0 | 9.0 | c005: flank <=3 <= old <=3; c100: flank <=3 <= old <=3; flank adds up to 9 unsupported UMIs per null library > 5 |

## Computational detection threshold (source, background and layout-specific; lowest level met in >= 2 of 3 seeds, and every higher level met)

A threshold that is not established (">300 (not established)") is not a negative result; "<=3" is left-censored. Last column: seeds with the call by level 3/10/30/100/300.

### D03

Call: combined

| arm | layout | threshold | non-monotone | seeds called by level |
|---|---|---|---|---|
| flank | c005 | 10 |  | 3:0;10:3;30:3;100:3;300:3 |
| flank | c100 | 10 |  | 3:0;10:3;30:3;100:3;300:3 |
| old | c005 | 10 |  | 3:0;10:3;30:3;100:3;300:3 |
| old | c100 | 10 |  | 3:0;10:3;30:3;100:3;300:3 |

Call: star

| arm | layout | threshold | non-monotone | seeds called by level |
|---|---|---|---|---|
| star | c005 | <=3 |  | 3:3;10:3;30:3;100:3;300:3 |
| star | c100 | <=3 |  | 3:3;10:3;30:3;100:3;300:3 |

Call: kb_only

| arm | layout | threshold | non-monotone | seeds called by level |
|---|---|---|---|---|
| flank | c005 | 10 |  | 3:0;10:3;30:3;100:3;300:3 |
| flank | c100 | 10 |  | 3:0;10:3;30:3;100:3;300:3 |
| old | c005 | 10 |  | 3:0;10:3;30:3;100:3;300:3 |
| old | c100 | 10 |  | 3:0;10:3;30:3;100:3;300:3 |

### R01

Call: combined

| arm | layout | threshold | non-monotone | seeds called by level |
|---|---|---|---|---|
| flank | c005 | <=3 |  | 3:3;10:3;30:3;100:3;300:3 |
| flank | c100 | <=3 |  | 3:3;10:3;30:3;100:3;300:3 |
| old | c005 | <=3 |  | 3:2;10:3;30:3;100:3;300:3 |
| old | c100 | <=3 |  | 3:2;10:3;30:3;100:3;300:3 |

Call: star

| arm | layout | threshold | non-monotone | seeds called by level |
|---|---|---|---|---|
| star | c005 | <=3 |  | 3:3;10:3;30:3;100:3;300:3 |
| star | c100 | <=3 |  | 3:3;10:3;30:3;100:3;300:3 |

Call: kb_only

| arm | layout | threshold | non-monotone | seeds called by level |
|---|---|---|---|---|
| flank | c005 | <=3 |  | 3:3;10:3;30:3;100:3;300:3 |
| flank | c100 | <=3 |  | 3:3;10:3;30:3;100:3;300:3 |
| old | c005 | <=3 |  | 3:2;10:3;30:3;100:3;300:3 |
| old | c100 | <=3 |  | 3:2;10:3;30:3;100:3;300:3 |

## Mean recovery (recovered / spiked) by level

| background | arm | layout | 3 | 10 | 30 | 100 | 300 |
|---|---|---|---|---|---|---|---|
| D03 | flank | c005 | 0.56 | 0.78 | 0.77 | 0.79 | 0.79 |
| D03 | flank | c100 | 0.56 | 0.78 | 0.77 | 0.78 | 0.78 |
| D03 | old | c005 | 0.56 | 0.57 | 0.69 | 0.65 | 0.64 |
| D03 | old | c100 | 0.56 | 0.57 | 0.69 | 0.65 | 0.64 |
| D03 | star | c005 | 1.00 | 1.00 | 1.00 | 0.99 | 1.00 |
| D03 | star | c100 | 1.00 | 1.00 | 1.00 | 0.99 | 1.00 |
| R01 | flank | c005 | 1.00 | 1.00 | 0.96 | 0.97 | 0.96 |
| R01 | flank | c100 | 1.00 | 1.00 | 0.96 | 0.97 | 0.96 |
| R01 | old | c005 | 0.89 | 0.87 | 0.78 | 0.83 | 0.84 |
| R01 | old | c100 | 0.89 | 0.87 | 0.78 | 0.83 | 0.84 |
| R01 | star | c005 | 1.00 | 1.00 | 0.98 | 1.00 | 1.00 |
| R01 | star | c100 | 1.00 | 1.00 | 0.98 | 1.00 | 1.00 |

kb arms: NET kb molecules (unique + allocated ambiguous, spike barcodes, truth accessions, minus the unspiked replicate of the same background subsample) per spiked read, not capped. star: any alignment on a truth contig.

## Background baseline of the truth virus in the spike barcodes (unspiked replicate)

| background | arm | spiked conditions | mean | max | flag (max >= 3) |
|---|---|---|---|---|---|
| D03 | flank | 30 | 0 | 0 | ok |
| D03 | old | 30 | 0 | 0 | ok |
| R01 | flank | 30 | 0 | 0 | ok |
| R01 | old | 30 | 0 | 0 | ok |

## Null floor (unspiked replicates)

| background | arm | metric | libraries | mean | max |
|---|---|---|---|---|---|
| D03 | old | kb_unsupported_molecules | 3 | 0.333 | 1 |
| D03 | old | n_false_calls | 3 | 0 | 0 |
| D03 | flank | kb_unsupported_molecules | 3 | 10.2 | 12.5714 |
| D03 | flank | n_false_calls | 3 | 0 | 0 |
| D03 | star | star_nonspike_reads_total | 3 | 6.67 | 8 |
| D03 | star | star_contigs_with_reads | 3 | 1.33 | 2 |
| R01 | old | kb_unsupported_molecules | 3 | 0.333 | 1 |
| R01 | old | n_false_calls | 3 | 0 | 0 |
| R01 | flank | kb_unsupported_molecules | 3 | 7 | 9 |
| R01 | flank | n_false_calls | 3 | 0 | 0 |
| R01 | star | star_nonspike_reads_total | 3 | 0 | 0 |
| R01 | star | star_contigs_with_reads | 3 | 0 | 0 |

## STAR gate evidence: fate of spikes by class (all levels, seeds and layouts pooled per background)

Pooled over levels, seeds and layouts: seeds are resamples of ONE source library per virus (pseudo-replicates, not independent libraries), the high levels dominate the pool, and class composition differs between backgrounds, so pooled recovery is not comparable across backgrounds. n is shown beside recovery with a 95% Wilson interval; classes with n < 30 are marked insufficient.

| background | class | spiked (n) | on truth contig | recovery | 95% CI | n check | host | other virus | other | absent |
|---|---|---|---|---|---|---|---|---|---|---|
| D03 | all | 2658 | 2646 | 0.995 | [0.99, 1.00] | ok | 0 | 25 | 0 | 12 |
| D03 | antisense_in_gene | 474 | 472 | 0.996 | [0.98, 1.00] | ok | 0 | 12 | 0 | 2 |
| D03 | flank3_unannotated | 252 | 252 | 1.000 | [0.98, 1.00] | ok | 0 | 3 | 0 | 0 |
| D03 | intergenic | 92 | 92 | 1.000 | [0.96, 1.00] | ok | 0 | 0 | 0 | 0 |
| D03 | sense_in_gene | 1840 | 1830 | 0.995 | [0.99, 1.00] | ok | 0 | 10 | 0 | 10 |
| R01 | all | 2658 | 2654 | 0.998 | [1.00, 1.00] | ok | 0 | 0 | 0 | 4 |
| R01 | antisense_in_gene | 26 | 26 | 1.000 | [0.87, 1.00] | insufficient | 0 | 0 | 0 | 0 |
| R01 | flank3_unannotated | 6 | 6 | 1.000 | [0.61, 1.00] | insufficient | 0 | 0 | 0 | 0 |
| R01 | sense_in_gene | 2626 | 2622 | 0.998 | [1.00, 1.00] | ok | 0 | 0 | 0 | 4 |

Fates are the PRIMARY alignment (host = host competition; other virus = primary on another viral contig); on-truth counts any alignment on the truth set.

## NOT measured

- One source virus per background: HHV-6B reads into D03, EBV (type 1 + type 2) reads into R01_SOC4. No other virus, no
  other source library, and no source-background chemistry crossing was tested, so thresholds do not transfer to other viruses.
- NH policy: the EBV source keeps reads with NH <= 2 whose records all lie on the two EBV genomes (`--max-nh 2`, an
  unvalidated choice); HHV-6B keeps NH = 1 only. The effect of that choice on the threshold was not measured.
- No error-free simulation arm in this grid (the `synth` probes were not run here): recovery losses are not separated
  into read-sequence effects versus index effects.
- "First/second source" caveat: spikes come from ONE real library per virus; seeds resample reads and barcodes from that
  library, not new libraries, so seed agreement reflects sampling within one source only.
- Per-class recovery is available for the STAR arm only: kb molecules are summed per barcode and gene, not per spike, so
  the kb arms have class "all" only. The kb arms get their STAR half of the combined call from the old-index run.
- STAR "absent" lumps reads unmapped for any reason (including too short after the 90 nt trim) and reads lost before the
  BAM (host-filter loss); the two are not separated here. The product host-filter loss column was not extracted.
- Thresholds are on a 5-level grid (3, 10, 30, 100, 300): a threshold is a grid level, not an interpolated value; "<=3" is left-censored and ">300 (not established)" is NOT a negative result (it means the grid did not establish a threshold). The null floor
  is 3 libraries per background; maxima of 3 values bound nothing beyond those libraries.
- Host and virus used the pinned indexes of this repository; depth-matching to real libraries (~45 M pairs) is by subsample.

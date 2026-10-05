# F-018 — Since SW-13, a Drop-seq run without `-w` is pre-filtered to kb's data-derived cell list

**Status:** confirmed (2026-09-29) · **Tags:** barcode-correction, Drop-seq, kb, cell-calling, denominators, SW-13

## Observation

The SW-13 fix stopped passing `kb count -w None` when no on-list is given. Job
25666447 reran EBV and HSV-1 on the v1 combined index, so both runs used the
same index before and after the fix. `kb_info.json` shows the command
sequence per chemistry:

- **10xv2 (EBV, SRR12682296):** `bustools correct -w 10x_version2_whitelist.txt`,
  the official on-list.
- **Drop-seq (HSV-1, SRR8315713):** `bustools allowlist -o whitelist.txt`, then
  `bustools correct -w whitelist.txt`. The on-list is built from the data by a
  knee.

| | uncorrected (`-w None`) | corrected (current code) |
|---|---:|---:|
| EBV molecules | 906,202 | 885,933 (−2.2 %) |
| EBV barcodes (all / infected) | 793,308 / 71,497 | 409,752 / 57,762 |
| HSV-1 molecules | 32,404 | 23,170 (**−28 %**) |
| HSV-1 barcodes (all / infected) | 1,809,992 / 10,227 | **5,468** / 1,127 |
| HHV-2 `possible_em_bleed` on HSV-1 | 54.8 | 34.8 |

For SRR8315713, `inspect.json` reports 5,468 barcodes on the allowlist (0.30 %).
The `summary.txt` headline now matches `viral_summary.tsv` in both runs
(885,939 and 23,204.7; SW-16).

## Interpretation

- **10x behaves as intended.** Barcodes off the official on-list are dropped or
  corrected. The −2.2 % change in EBV molecules is that correction.
- **For Drop-seq, the allowlist acts as a cell caller.** Every barcode outside
  kb's knee is discarded before ViralScan counts molecules or calls cells:
  - the "all barcodes" molecule totals lose the ambient and low-count barcodes;
  - the infected-cell denominator becomes kb's list of 5,468;
  - ViralScan's own cell calling (knee or EmptyDrops) then runs on
    pre-selected cells.
- SW-20 makes multimap read kb's corrected BUS, so the multimap layers and
  detection are pre-filtered too.

## Action

- `PLAN SW-21`, which the chemistry module decides (`MECH-D`): for technologies
  without an official on-list, either skip correction (kb's bypass value) or
  keep the allowlist and document the pre-filter.
- Until then, Drop-seq absolute numbers from current code are not comparable
  with pre-SW-13 runs, and 10x numbers move by about 2 %.

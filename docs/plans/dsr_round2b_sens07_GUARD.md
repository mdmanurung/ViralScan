# DSR round 2b (SENS-07 planted sensitivity) - comparability guard (binding on every agent)

Template. Copy to `dsr_round2b/GUARD.md` after the user has answered Q48 and Q49 and `config.tsv` is filled in;
never edit the copy afterwards. Mirrors `docs/plans/dsr_round2_GUARD.md` (PLAN `DSR-16`, `SENS-07`).
Every run is exploratory: labelled `explore`, unscored, one library per dataset (n = 1), no scored-gate (G5b) and
no VAL-01 claim. No pass/fail threshold anywhere.

Pinned code: /exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/35de249 (clean `main` commit; never edit it).
**Deviation from round 2 (user, Q49):** round 2 ran on `vs_pinned/e7c1b08`; round 2b uses `35de249` because the
`SENS-07` scorer-group and exact-sequence receipt repairs landed after `e7c1b08`. The planted run is a separate
measurement, not a re-run of round 2's cells, so the code difference does not change any round-2 cell; the ledger
note and every row's `code_sha.txt` record it.
Reference: `cat42d` (`viral_ref_cat42d/build`), unchanged from rounds 1 and 2.
Rows: dsr_round2/manifest.tsv, read-only (array index i = data row i, 0-based); chemistry:
dsr_round2/chem/<sample>.json, read-only -> `-x` and `--strand` for both runs of a row.
Per-row choices (run yes/no, expected family, n_background): `dsr_round2b/config.tsv` (frozen copy of
`config/dsr_round2b_sens07.template.tsv`). Rows with `run=no` are recorded n/a, not run and not "fixed":

| Row | Reason (n/a unless the user decides) |
|---|---|
| `SRR20710641` (hhv6b) | R1 is 150 nt, not barcode + UMI (26 nt); `plant` refuses it |
| `GSM4682311` (kshv_ebv_gse154900) | R1 is 28 nt, documented chemistry 10xv2 needs 26 |

Per row (scripts/slurm_dsr_round2b_sens07_array.sh): `plant` -> planted/ + unplanted/ twin FASTQs of the first
`n_background` pairs; both through `viralscan` with round 2's `combined_off` flags (`-i panel.idx`, `--cell-calling
emptydrops`, `--anellovirus-gene-ids`, no filter flag); `score` (planted minus twin, layer `X` of
`adata_multimap.h5ad`). Plant parameters (Q48 defaults): exact genomes only (no held-out), families
`Orthoherpesviridae Herpesviridae Papillomaviridae Polyomaviridae` plus the row's expected family, `--per-family 2`,
`--seed 20261012`.

Environment (set before every sbatch):
  ROUND_DIR=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round2b
  CODE_DIR=/exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/35de249
Pass `-o $ROUND_DIR/logs/%x_%A_%a.log` on every sbatch line (the in-script #SBATCH -o is only the default).

## Explore ledger (mandatory)
Every finished `viralscan` run (arms `combined_off_planted`, `combined_off_unplanted`) is appended to
`$ROUND_DIR/explore_ledger.tsv` by the run script (`scripts/explore_ledger.py append`, label `explore`). Round 2b has
its own ledger and chain; it never writes to `dsr_round2/explore_ledger.tsv`. Do not edit, sort, truncate or
hand-write the ledger. A run that finished but was not logged is reported, not back-filled. The script refuses a
manifest row naming a Track B accession (`config/track_b_denylist.txt`); if it does, STOP and report. After the last
array finishes: `python scripts/explore_ledger.py --ledger $ROUND_DIR/explore_ledger.tsv check` must exit 0.

## Throttle
`--array=0-16%2`. The user has other long jobs running: always use the throttle given in your assignment.

## Reporting
Per dataset: `score/family_summary.tsv` (capture, Wilson interval, `capture_status`, `lod95_true_molecules`) and
`score/recovery_by_condition.tsv`. Public development set and institutional datasets (`covid`, `sfl_tonsil`) are
reported separately, never pooled. State "n = 1 library, internal-fragment forward-strand plants, first
`n_background` read pairs" next to every number.

## You MAY change
--mem, --time, --cpus-per-task, partition, array throttle (%N, never above 2), and resubmit a failed array index.
On OOM or TIMEOUT: resubmit that index ONCE at 2x mem/time, memory capped at 480G (the plant step holds every
(CB, UMI) of the first `n_background` pairs in memory: about 10 GB per 50M pairs).

## You MAY NOT change
-x, -w, --strand, --cell-calling, --read-filter, the reference/index, the multimap method, the code path,
`PLANT_FAMILIES`, `PER_FAMILY`, `SEED`, `n_background`, the manifest, the chem JSONs, the config, the ledger, or any
file under the repo / vs_pinned. No git commands. No package-code edits.
A failure that seems to need one of these: STOP that cell, leave it, and report it (n/a: <reason> or the user).
emptyDrops failure => report `n/a: emptydrops failed`. NEVER switch to none/knee.
A planting error (too few barcodes for the conditions, geometry mismatch): report it; do not lower
`--cells-per-condition` or trim reads to make it pass.

## Write only under dsr_round2b/ (runs/, status/, logs/). Never write under dsr_round1/ or dsr_round2/ (read-only). Do not write inside the repo tree.
## Submit from the repo dir: /exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
## Wait with Monitor or a Bash until-loop (run_in_background), checking sacct no more than once per 5 minutes.

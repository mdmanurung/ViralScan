# DSR round 2c (SENS-07 planted sensitivity, per-family seeds + trimmed HHV-6B) - comparability guard (binding on every agent)

Template. Copy to `dsr_round2c/GUARD.md` after round 2b has finished, the `PIN` below is filled and `config.tsv` is
frozen; never edit the copy afterwards. Mirrors `docs/plans/dsr_round2b_sens07_GUARD.md` (PLAN `DSR-16`, `SENS-07`).
Every run is exploratory: labelled `explore`, unscored, one library per dataset (n = 1), no scored-gate (G5b) and
no VAL-01 claim. No pass/fail threshold anywhere.

Pinned code: /exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/PIN (a clean `main` commit that contains
`--selection-seed-mode` and `--trim-r1-to-geometry` in `scripts/sens07_planted_sensitivity.py`; never edit it). The
driver script is that same commit's `scripts/slurm_dsr_round2b_sens07_array.sh` (round 2c is that script with the
settings below; the script refuses a `CODE_DIR` whose planter lacks the flags it is asked to pass).
Reference: `cat42d` (`viral_ref_cat42d/build`), unchanged from rounds 1, 2 and 2b.
Rows: dsr_round2/manifest.tsv, read-only (array index i = data row i, 0-based); chemistry:
dsr_round2/chem/<sample>.json, read-only -> `-x` and `--strand` for both runs of a row.
Per-row choices (run yes/no, expected family, n_background, trim_r1): `dsr_round2c/config.tsv` (frozen copy of
`config/dsr_round2c_sens07.template.tsv`).

## Changes from round 2b (user, 2026-10-13)
* **Q64, per-family seeds.** `SELECTION_SEED_MODE=per-family` (flag `--selection-seed-mode per-family`): each family's
  genomes are drawn from a generator seeded with sha256(`SEED`, family), so a family's planted genomes are the same in
  every dataset whatever other families that dataset plants. Round 2b drew all families from one shared generator
  (`global`), so a row adding e.g. `Coronaviridae` planted different Papillomaviridae genomes than the other rows;
  round 2c removes that. Round 2b's numbers are not re-derived and not edited; compare 2b and 2c per family only with
  this difference in mind.
* **Q63, HHV-6B.** Row `SRR20710641` (hhv6b) has R1 of 150 nt (barcode + UMI = 26 nt for 10xv2) and is now planted:
  `trim_r1=yes` passes `--trim-r1-to-geometry`, which trims R1 to the first 26 nt (sequence and quality) in BOTH the
  planted and the unplanted twin (identical trim; R2 is untouched). This is a declared read change, recorded in that
  row's `plant/plant_manifest.json` (`r1_trim`: original 150, trimmed 26; `read_changes`) and reported next to every
  HHV-6B number. The trim never rescues an R1 shorter than barcode + UMI.
* **GSM4682311 stays n/a.** Its R1 is 28 nt against 26 (longer, so the same trim rule would apply), but it is a
  `FORCE_TECH` row whose chemistry `check-chemistry` could not resolve; the user's Q63 keeps it n/a.

Per row (scripts/slurm_dsr_round2b_sens07_array.sh): `plant` -> planted/ + unplanted/ twin FASTQs of the first
`n_background` pairs; both through `viralscan` with round 2's `combined_off` flags (`-i panel.idx`, `--cell-calling
emptydrops`, `--anellovirus-gene-ids`, no filter flag); `score` (planted minus twin, layer `X` of
`adata_multimap.h5ad`). Plant parameters (Q48 defaults): exact genomes only (no held-out), families
`Orthoherpesviridae Herpesviridae Papillomaviridae Polyomaviridae` plus the row's expected family, `--per-family 2`,
`--seed 20261012`, `--selection-seed-mode per-family`.

Environment (set before every sbatch):
  ROUND_DIR=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round2c
  CODE_DIR=/exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/PIN
  SELECTION_SEED_MODE=per-family
Pass `-o $ROUND_DIR/logs/%x_%A_%a.log -J dsr2c_sens07` on every sbatch line (the in-script `#SBATCH` lines are
round 2b's defaults). Check with `DRY_RUN=1` first: every row's `== row` line must show `selection=per-family`, and
row 1 (`SRR20710641`) `trim_r1=yes` with `--trim-r1-to-geometry` in its `plant` command.

## Explore ledger (mandatory)
Every finished `viralscan` run (arms `combined_off_planted`, `combined_off_unplanted`) is appended to
`$ROUND_DIR/explore_ledger.tsv` by the run script (`scripts/explore_ledger.py append`, label `explore`). Round 2c has
its own ledger and chain; it never writes to `dsr_round2/` or `dsr_round2b/` ledgers. Do not edit, sort, truncate or
hand-write the ledger. A run that finished but was not logged is reported, not back-filled. The script refuses a
manifest row naming a Track B accession (`config/track_b_denylist.txt`); if it does, STOP and report. After the last
array finishes: `python scripts/explore_ledger.py --ledger $ROUND_DIR/explore_ledger.tsv check` must exit 0.

## Throttle
`--array=0-16%2`. The user has other long jobs running: always use the throttle given in your assignment.

## Reporting
Per dataset: `score/family_summary.tsv` (capture, Wilson interval, `capture_status`, `lod95_true_molecules`) and
`score/recovery_by_condition.tsv`. Public development set and institutional datasets (`covid`, `sfl_tonsil`) are
reported separately, never pooled. State "n = 1 library, internal-fragment forward-strand plants, first
`n_background` read pairs" next to every number, and for HHV-6B also "R1 trimmed from 150 to 26 nt in planted and
unplanted".

## You MAY change
--mem, --time, --cpus-per-task, partition, array throttle (%N, never above 2), and resubmit a failed array index.
On OOM or TIMEOUT: resubmit that index ONCE at 2x mem/time, memory capped at 480G (the plant step holds every
(CB, UMI) of the first `n_background` pairs in memory: about 10 GB per 50M pairs).

## You MAY NOT change
-x, -w, --strand, --cell-calling, --read-filter, the reference/index, the multimap method, the code path,
`PLANT_FAMILIES`, `PER_FAMILY`, `SEED`, `SELECTION_SEED_MODE`, `trim_r1`, `n_background`, the manifest, the chem JSONs,
the config, the ledger, or any file under the repo / vs_pinned. No git commands. No package-code edits.
A failure that seems to need one of these: STOP that cell, leave it, and report it (n/a: <reason> or the user).
emptyDrops failure => report `n/a: emptydrops failed`. NEVER switch to none/knee.
A planting error (too few barcodes for the conditions, geometry mismatch): report it; do not lower
`--cells-per-condition`, do not add `trim_r1` to another row and do not trim reads yourself to make it pass.

## Write only under dsr_round2c/ (runs/, status/, logs/). Never write under dsr_round1/, dsr_round2/ or dsr_round2b/ (read-only). Do not write inside the repo tree.
## Submit from the repo dir: /exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
## Wait with Monitor or a Bash until-loop (run_in_background), checking sacct no more than once per 5 minutes.

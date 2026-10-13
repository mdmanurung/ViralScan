# DSR round 2 - comparability guard (binding on every agent)

Template. Copy to `dsr_round2/GUARD.md` once `<SHA>` is filled in; never edit the copy afterwards.
Adapted from round 1's `dsr_round1/GUARD.md` (grill session 3, 2026-10-12; PLAN `DSR-16`).
Every run is exploratory: labelled `explore`, unscored, no scored-gate (G5b) claim.

Pinned code: /exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/<SHA> (clean `main` commit; never edit it).
Reference: `cat42d` (`viral_ref_cat42d/build`, virus-only index `build_virus_only`), unchanged from round 1.
Manifest (frozen, read-only): dsr_round2/manifest.tsv = round 1's manifest minus the two rows below
(data row i = array index i, 0-based):

| Dropped | Reason (recorded n/a) |
|---|---|
| `ERR13027122` (tonsil_atlas) | n/a: HTO/ADT library, not a gene-expression library |
| `GSM4682312` (kshv_ebv_gse154900) | n/a: no barcode read |

Chemistry JSON (frozen, reused from round 1, never re-run): dsr_round1/chem/<sample>.json copied to
dsr_round2/chem/<sample>.json  -> -x and --strand for every arm.

Arms differ ONLY by flags (implemented in scripts/slurm_dsr_rerun_array.sh, env ARM=):
  combined_off | combined_artefact (--read-filter artefact) | twostep (--host-filter starsolo, virus-only index)
`twostep` calls cells on the STARsolo host matrix in the same run (DEF-05); there is no `twostep_v2` redetect pass.
The round-1 baseline for `twostep` is round 1's `twostep_v2`.

Environment (set before every sbatch):
  ROUND_DIR=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round2
  CODE_DIR=/exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/<SHA>
  (REF_BUILD defaults to the cat42d build; do not override it.)
Pass `-o $ROUND_DIR/logs/%x_%A_%a.log` on every sbatch line (the in-script #SBATCH -o is round 1's).

## Explore ledger (mandatory)
Every finished run is appended to `$ROUND_DIR/explore_ledger.tsv` by the run script
(`scripts/explore_ledger.py append`: timestamp, dataset, sample, arm, code sha, reference id, label `explore`).
Do not edit, sort, truncate or hand-write the ledger. A run that finished but was not logged is reported, not
back-filled. The run script refuses a manifest row naming a Track B accession (`config/track_b_denylist.txt`);
if it does, STOP and report. After the last array finishes: `python scripts/explore_ledger.py --ledger
$ROUND_DIR/explore_ledger.tsv check` must exit 0.

## Throttle
Combined arms (`combined_off`, `combined_artefact`): `--array=<rows>%4`. `twostep`: `--array=<rows>%2`.
The user has other long jobs running: always use the throttle given in your assignment.

## Reporting
The public development set is reported on its own. Institutional datasets (`covid` x213/x216 and
`sfl_tonsil` x223) are reported separately and never pooled with it. Headline cells: the DSR-16 common-cell set
(`scripts/dsr_common_cells.py`), not each arm's own `pct_infected_called`. Round 2 vs round 1:
`scripts/dsr_round_compare.py`, filtered on `detected` (AUD-21).

## You MAY change
--mem, --time, --cpus-per-task, partition, array throttle (%N, never above the limits above), and resubmit a failed
array index. On OOM or TIMEOUT: resubmit that index ONCE at 2x mem/time, memory capped at 480G.

## You MAY NOT change
-x, -w, --strand, --cell-calling, --read-filter, the reference/index, the multimap method, the code path,
the manifest, the chem JSONs, the ledger, or any file under the repo / vs_pinned. No git commands. No package-code edits.
A failure that seems to need one of these: STOP that cell, leave it, and report it (becomes n/a: <reason> or goes to the user).
emptyDrops failure ("insufficient unique points" etc.) => report `n/a: emptydrops failed`. NEVER switch to none/knee.
A chem JSON that is missing, exit!=0 with chemistry null, or says unresolved => the cell is n/a, report it.

## Write only under dsr_round2/ (runs/, status/, logs/). Never write under dsr_round1/ (read-only). Do not write inside the repo tree.
## Submit from the repo dir: /exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
## Wait with Monitor or a Bash until-loop (run_in_background), checking sacct no more than once per 5 minutes.

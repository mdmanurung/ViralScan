# Completion plan — reference-strategy 2×2 benchmark

**Goal:** finish the 2×2 (reference strategy `combined` vs `two_step` × aligner
`STARsolo` vs `ViralScan`) across the three target viruses so the
reference-strategy comparison and the Selectivity Index become computable, and
the benchmark can be honestly re-included in the manuscript (currently excluded).

**Design:** 3 viruses × 2 references × 2 aligners = 12 rows. Command manifest:
`benchmark_runs/reference_strategy_2026-06-28_fresh12/commands.jsonl`.
Summarizer: `scripts/summarize_reference_strategy.py` → `results/reference_strategy_benchmark.tsv`.

**Success criteria:** all 12 rows `status=complete`; `summarize_reference_strategy.py`
exits 0 with no `failed`/`blocked`/`incomplete`/`running_or_incomplete` rows; the
Selectivity Index (on-target ÷ off-target UMI, per aligner×reference) is computable
for EBV (the on-target virus).

---

## Current status — 4 of 12 complete

| virus | STARsolo·combined | STARsolo·two_step | ViralScan·combined | ViralScan·two_step |
|-------|-------------------|-------------------|--------------------|--------------------|
| HHV-6B | ✅ 1295 / 1876 UMI | ⏳ no final output | ✅ 1999 / 3217 UMI (11 off-tgt) | ⛔ blocked (no kallisto) |
| EBV    | ❌ failed (FASTQ) | ⏳ no final output | ⚠️ incomplete (no outputs) | ⛔ blocked (no kallisto) |
| HSV-1  | ✅ 0 / 0 | ✅ 0 / 0 | ⚠️ incomplete (no outputs) | ⛔ blocked (no kallisto) |

## Root causes (three failure classes)

1. **ViralScan `two_step` — all 4 blocked (HHV-6B, EBV, HSV-1, + the EBV/HSV-1 above).**
   `viralscan --host-filter kallisto` preflight (`_check_host_filter_tools`) requires
   standalone `kallisto` **and** `bustools` on `PATH`. The benchmark env
   (`evonk/.../test_viralscan`) ships `kb` + `bustools` but **no standalone `kallisto`**,
   and the run script's kb-python-bundled-kallisto PATH shim did not expose one.
   → **Fix: a conda env that has standalone `kallisto` (Step 1).**

2. **EBV STARsolo `combined` — failed:** `FATAL ERROR in reads input: quality string
   length is not equal to sequence length`. Geometry is correct (10xv2, CB=16/UMI=10),
   so this is a **malformed FASTQ record**, not a parameter bug. The *same* EBV sample
   (SRR12682296) processed successfully with STARsolo for the manuscript §3.3 comparison
   using a different (showcase) copy — the benchmark used a fresh ENA re-download
   (`benchmark_inputs/reference_strategy/SRR12682296/`, R1 2.27 GB / R2 8.67 GB) that
   STAR rejects. kallisto/ViralScan tolerate the record; STAR is stricter.
   → **Fix: repair or re-fetch the EBV FASTQ (Step 2).**

3. **ViralScan `combined` incomplete (EBV, HSV-1)** and **STARsolo `two_step`
   running_or_incomplete (HHV-6B, EBV):** rows started but produced no final marker /
   no `per_cell_viral.tsv`+`viral_summary.tsv`. Likely wall-time/OOM truncation on the
   deep samples, or a downstream step that never ran. → **Fix: re-run to completion and
   inspect the per-row logs if they fail again (Step 4).**

---

## Step 0 — Restore the ViralScan kallisto index — ✅ NO REBUILD NEEDED (2026-07-03)

**Discovered during preflight:** the *scratch* tree
`/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/` was deleted by scratch
cleanup — including the ViralScan combined index `-i`/`-t` that **all six ViralScan rows** use.
This — not just the missing `kallisto` binary — is why the ViralScan `combined` rows are
"incomplete".

**But a persistent ARCHIVE copy survives** (the Jun-22 build the Jun-28 benchmark used), so the
chosen "faithful" route (A) collapses to a **restore** — no `kb ref` rebuild:
`/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_showcase/viralscan_showcase/fullrun/refs/merged/`
- `index_plus_anellovirus.idx` (371 MB), `t2g_plus_anellovirus.txt` (226,005 host ENST +
  821 target-virus rows — EBV/HHV/HSV present), `transcriptome_plus_anellovirus.fa`,
  `viral_serratus_plus_anellovirus.gtf` (the `-gtf` the two_step ViralScan rows use).
Because it is the **same index** that produced the 4 complete rows, those stay valid — **only
the 8 incomplete rows need re-running**, not all 12.

**Version-compat gotcha (important):** the index was built with kb-python's bundled kallisto.
The **standalone conda `kallisto` 0.52 segfaults** reading both the combined and host indices,
but the **kb-python bundled kallisto reads both** (verified: 0.51.1 and 0.52.0 bundled). `kb count`
uses the bundled kallisto internally (fine), but `--host-filter kallisto` calls the standalone
`kallisto` on PATH — so the harness must **prepend the kb-python bundled kallisto dir to PATH**
so `which kallisto` resolves to the bundled binary (reads the host index), exactly as the original
run script did. The surviving source FASTAs/GTFs on archive remain a rebuild fallback if ever needed.

Rebuild materials (fallback only, if a *fresh* reference is ever wanted): human GRCh38-2024-A
`genome.fa`/`genes.gtf` + Serratus `fasta_split`/`split_gtf`/`Serratus_v2` + the archived
`create_final_transcriptome.slurm` recipe; or `scripts/build_bundled_panel_ref.py` (Route B).

**What survives (on archive) — the index is fully rebuildable:**
- Human GRCh38-2024-A `genome.fa` (3.0 G) + `genes.gtf` (1.6 G).
- Serratus viral source set: `.../evonk/old/intern/fasta_viruses/Serratus/{fasta_split (174),
  split_gtf (174), Serratus_v2/fasta (22), Serratus_v2/gtf (22)}` + `all_fastas.fasta`/`all_gtf.gtf`.
- The **original build recipe** `.../Serratus/create_final_transcriptome.slurm` (the exact
  `kb ref` command that produced the index).
- The **STARsolo combined index (33 G)** — intact, so the STARsolo side needs no rebuild.
- The host kallisto index (483 M) used by `--host-filter kallisto`.

**Two rebuild routes (a decision, not a default):**
- **Route A — faithful reconstruction.** Re-run the archived `create_final_transcriptome.slurm`
  `kb ref` (GRCh38 genome + Serratus fasta_split/split_gtf) and add the anellovirus panel
  (regenerable via ViralScan's bundled data) → `index_plus_anellovirus.idx` + t2g. Matches the
  Serratus sequences the surviving 33 G STARsolo index was built from, so the aligner axis stays
  matched. ~1–2 h build (96 G mem per the original recipe), then re-run all 12 rows.
  Con: the reference is a bespoke Serratus set, not something a reader can reproduce.
- **Route B — clean rebuild via `scripts/build_bundled_panel_ref.py`.** ViralScan's maintained
  builder (Ensembl human cDNA + NCBI viral accessions + bundled anellovirus) → a reproducible,
  documented reference. Con: it is a *different* viral panel than the Serratus STARsolo index,
  so to keep the aligner comparison fair the STARsolo index must be rebuilt from the same source
  too — more work, but far more defensible for the manuscript. Aligns with the repo's
  "prefer native `viralscan build-ref`" convention.

**Consequence (restore path taken):** because the archive copy **is the same index** that
produced the 4 complete rows, those stay valid — **only the 8 incomplete rows need re-running**
(array indices `1,3,4,5,6,7,10,11`). (A *rebuild* — Routes A/B above — would instead force
re-running all 12, since a fresh index is not byte-identical. Restore avoids that.)

## Step 0b — Scratch cleanup damage is broader than the index (audit, 2026-07-03)

Re-parsing the run dir with `summarize_reference_strategy.py` (which reads live output files)
plus a per-input existence check revealed the scratch filesystem
(`/exports/para-lipg-hpc`, 93 % full) is **actively deleting benchmark data**, not just the index:

- **ViralScan index**: deleted from scratch — **restorable from archive** (Step 0). ✅
- **EBV FASTQs** (`benchmark_inputs/reference_strategy/SRR12682296/`): **survive**. ✅
- **HHV-6B FASTQs** (`viralscan_showcase/data/hhv6_carT_ref/SRR20710641/`): **GONE.** ❌
- **HSV-1 FASTQs** (`viralscan_showcase/data/hsv1_fibroblast/SRR8315713/`): **GONE.** ❌
- **Row outputs**: 2 rows that were "complete" on Jun 28 have since **lost their output files**
  (`hhv6b__viralscan__combined`, `hsv1__starsolo__two_step`), while 2 rows the summarizer had
  marked failed/incomplete were actually **complete all along** (`ebv__starsolo__combined` —
  223 MB matrix, "ALL DONE"; `ebv__starsolo__two_step`) — the original TSV captured stale
  early-attempt statuses.

**Corrected current state (summarizer on live files): 4 complete, 8 to re-run.**
- Complete (survive): `hhv6b__starsolo__combined`, `ebv__starsolo__combined`,
  `ebv__starsolo__two_step`, `hsv1__starsolo__combined` (array 0,4,5,8).
- Re-run needed: array **1,2,3,6,7,9,10,11** (all 6 ViralScan rows + the 2 STARsolo two_step rows).

**Consequences / new decisions this raises:**
1. **Re-fetch HHV-6B + HSV-1 FASTQs** from ENA (public: SRR20710641, SRR8315713 — several GB each;
   `scripts/fetch_reference_strategy_fastqs.py` is the tracked helper). Without them, 6 of the 8
   re-run rows cannot run.
2. **Scratch instability**: since cleanup is ongoing, the 4 surviving "complete" rows may vanish
   too. For a stable, self-consistent result, either (a) re-run **all 12** into a fresh dir and
   summarize immediately, or (b) accept the risk of "restore + re-run 8". Writing reference
   artifacts + key outputs to **archive** (persistent) rather than scratch is advisable.
3. Given the benchmark is currently **excluded from the manuscript**, weigh whether completing it
   (re-fetch 2 datasets + full re-run, fighting active cleanup) is worth the compute now.

## Step 1 — Dedicated conda environment — ✅ DONE (2026-07-03)

Built `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench`
from the documented `environment.yml` spec + `samtools`:
Python **3.11.15**; **`kallisto` 0.52.0** (the previously-missing piece), `bustools`
0.45.1, `kb-python`, `STAR` 2.7.11b, `snakemake-minimal`, `samtools`, and the full
scientific stack; `viralscan` 2.5.0 installed editable (`pip install --no-deps -e .`).
(`seqkit` is being added for the Step 2 sanitize; not required — Step 2 Option B needs no extra package.)

Verified under a **clean activation** (what the SLURM job sees):
```bash
source /share/software/tools/miniconda/3.10/23.3.1/etc/profile.d/conda.sh
conda activate /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
```
`python kb kallisto bustools snakemake STAR samtools viralscan` all resolve to the env;
`python --version` = 3.11.15; `viralscan --version` = 2.5.0.

**Two_step fix (refined):** the conda standalone `kallisto` 0.52 *segfaults* on the older host
index, so the harness prepends the **kb-python bundled kallisto** to PATH (as the original run
script did) — `which kallisto` then resolves to the bundled binary, which reads the host index.
That satisfies `_check_host_filter_tools` **and** avoids the segfault. Verified in `viralscan_bench`.

> In an *interactive* shell the login profile may auto-activate the `codex` env and shadow
> `python`; `conda deactivate` it first, or rely on the SLURM job's clean single activation.

## Step 3 — Staged run harness — ✅ DONE (fresh12b)

Non-destructively staged at
`/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_runs/reference_strategy_2026-06-28_fresh12b/`
(preserves the original fresh12 run + provenance):
- `run_reference_strategy_array.sh` — activation → `viralscan_bench` (+ bundled-kallisto shim);
  `RUN_DIR` and `#SBATCH --output/--error` → fresh12b.
- `commands.jsonl` — 12 rows, JSON-validated; ViralScan `-i`/`-t`/`-gtf` repointed to the surviving
  **archive** `merged/` copy; all outputs → fresh12b. `--host-index` (evonk) verified present + readable.

## Step 4 — READY TO SUBMIT: EBV 2×2 (chosen scope, 2026-07-04)

Chosen scope: **EBV 2×2 only** — the on-target virus, the only one whose Selectivity Index is
computable, and completable with **no re-fetch** (EBV FASTQs survive; index restored; kallisto
fixed). Self-contained into `fresh12b` (re-runs all 4 EBV rows, so it does not depend on the
fragile surviving outputs that scratch cleanup keeps deleting). All 4 EBV rows' inputs verified
present; scratch ViralScan checkout (v2.5.0) present; outputs → fresh12b.

**Submit (stage-only — you submit):**
```bash
sbatch --array=4,5,6,7 \
  /exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_runs/reference_strategy_2026-06-28_fresh12b/run_reference_strategy_array.sh
```
Rows: 4 = ebv·starsolo·combined, 5 = ebv·starsolo·two_step, 6 = ebv·viralscan·combined,
7 = ebv·viralscan·two_step. (No EBV FASTQ sanitize needed — Step 2's premise was disproven:
STAR actually completed on this FASTQ; the "quality-string" status was a stale early-attempt artifact.)

**Summarize after the array finishes (from the archive repo checkout, viralscan_bench python):**
```bash
cd /exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench/bin/python \
  scripts/summarize_reference_strategy.py \
  --run-dir /exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_runs/reference_strategy_2026-06-28_fresh12b \
  --out /tmp/ebv_2x2.tsv --validate /tmp/ebv_2x2.tsv
```
The 4 EBV rows should read `complete`; the EBV combined-vs-two_step × STARsolo-vs-ViralScan cells
are then all populated → compute the Selectivity Index (on-target EBV UMI ÷ off-target UMI).

**Deferred (needs re-fetch; user chose not to do now):** HHV-6B (SRR20710641) + HSV-1 (SRR8315713)
rows — FASTQs deleted from scratch; re-fetch from ENA via `scripts/fetch_reference_strategy_fastqs.py`
before their rows can run. HSV-1 has no signal; HHV-6B carries count-layer/6A-cross-mapping confounds.

## Step 2 — Repair the EBV FASTQ (unblocks EBV STARsolo)

**Diagnosis (done): both EBV FASTQs pass `gzip -t`** (R1 OK, R2 OK), and record counts
match (R1==R2==127,045,580). So this is **not** a truncated download — it is a single
malformed record inside a valid gzip that STAR rejects but kallisto tolerated (ViralScan's
manuscript EBV run used this sample successfully). **Fix = sanitize, not re-fetch:**

```bash
# Option A (seqkit): locate + drop the bad record, keeping mates paired
seqkit seq --validate-seq -w0 R2.fastq.gz > /dev/null   # reports the offending record
# then write STAR-clean, still-paired copies (..._1.clean.fastq.gz / ..._2.clean.fastq.gz)
```
```awk
# Option B (no extra package): drop any record whose seq/qual lengths differ, per mate,
# then re-pair on read name. Pure zcat|awk|gzip — works with only the base env.
```
Keep R1/R2 paired (equal record counts) and point the EBV STARsolo rows at the cleaned copies.
Re-fetch is the last-resort fallback (ENA URLs in
`benchmark_inputs/reference_strategy/SRR12682296/source_urls.tsv`; no `sra-tools` needed —
these are direct `.fastq.gz` HTTPS URLs).

> Env note: `seqkit` is being added to `viralscan_bench` (slow mamba solve). Option B needs
> no extra package, so the sanitize step is not blocked on the seqkit install.

Cross-check the fix by confirming `bustools`/kallisto already accepted this sample
(ViralScan's manuscript EBV run did), i.e. STAR-specific strictness is the only barrier.

## Step 3 — Corrected run harness

Copy `run_reference_strategy_array.sh` to a `fresh12b` run dir and change the single
activation line to the new env:
```bash
conda activate /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
```
Keep the existing `commands.jsonl` (STAR paths are absolute to the `starsolo` env and
still valid; ViralScan rows use `python src/viralscan/menu.py` which now finds `kallisto`
on `PATH`). If the EBV FASTQ path changes after repair, regenerate the manifest with
`scripts/prepare_reference_strategy_benchmark.py` so the EBV rows point at the clean copy.

Pre-flight the references (all present, but confirm before a 12-h array):
`references/starsolo/{combined_GRCh38_2024A_serratus_plus_anellovirus,human_GRCh38_2024A,
all_virus_serratus_plus_anellovirus}` and the ViralScan kallisto index
`.../viralscan_showcase/fullrun/refs/merged/index_plus_anellovirus.idx` (+ the host
kallisto index used by `--host-filter kallisto`).

## Step 4 — Re-run the 8 non-complete rows

Target only the incomplete rows by SLURM array index (0-based, from `commands.jsonl`):
- `1` hhv6b·starsolo·two_step, `3` hhv6b·viralscan·two_step
- `4` ebv·starsolo·combined (after Step 2), `5` ebv·starsolo·two_step,
  `6` ebv·viralscan·combined, `7` ebv·viralscan·two_step
- `10` hsv1·viralscan·combined, `11` hsv1·viralscan·two_step

```bash
sbatch --array=1,3,4,5,6,7,10,11 run_reference_strategy_array.sh   # (fresh12b copy)
```
`--cpus-per-task=8 --mem=128G --time=12:00:00` are already set. If any ViralScan row is
truncated by wall-time/OOM on the deep EBV sample, bump `--time`/`--mem` and inspect that
row's `logs/%x_%A_%a.err` before re-submitting (do not silently retry).

## Step 5 — Verify + summarize

```bash
python scripts/audit_reference_strategy.py    # refresh reference/fastq audits
python scripts/summarize_reference_strategy.py  # must exit 0, all 12 rows complete
```
Then compute the Selectivity Index (on-target ÷ off-target UMI per aligner×reference) —
now possible because EBV rows have real signal — and refresh the tracked run packet
(`run_packets/2026-06-28_fresh12/failure_summary.tsv` → empty; add a completion note).

## Step 6 — Fair-comparison harmonization (before any manuscript claim)

The 2×2 is only publishable if the aligner axis is apples-to-apples:
- **Count-layer parity:** STARsolo reports `GeneFull.raw`; ViralScan reports
  `per_cell_viral.viral_umi`. Compare on a common layer (e.g. both unique-count, or note
  ViralScan's multimap layer explicitly) rather than cross-layer.
- **HHV-6A↔6B cross-mapping:** ViralScan's extra HHV-6B calls + small off-target may
  reflect 6A cross-mapping, not pure sensitivity — separate before claiming superiority.
- **HSV-1 has no signal** in this sample (0 under all cells), so it contributes only a
  specificity/negative check, not a sensitivity comparison.

## Step 7 — Manuscript re-inclusion decision

Only after Steps 5–6: if all 12 complete and the comparison is fair, lift the exclusion in
`analysis/reference_strategy_benchmark/REFERENCE_STRATEGY_BENCHMARK.md` and
`docs/PUBLICATION_READINESS.md`, and add the reference-strategy result to the manuscript.
Otherwise keep it as provenance-only. Update `PLAN.md`.

---

## Risks / notes

- Deep EBV sample (R2 8.67 GB) is the main wall-time/memory risk for the ViralScan rows.
- STAR paths in `commands.jsonl` point at the `starsolo` env by absolute path — fine, but
  if that env moves, regenerate the manifest.
- The env build + a full 12-h array are the long poles; everything else is minutes.
- Keep bulky outputs under gitignored `benchmark_runs/`; only lightweight provenance is tracked.

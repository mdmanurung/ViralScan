# Anellovirus remaining work — plan (2026-10-04)

Status: **superseded 2026-10-04** by `2026-10-04-package-completion-plan.md`
(user: finish package development before more experiments).
- A1 (REL-16) is now M1.2 there.
- A3 (evidence measures) is now M4.
- A2, Track B (`ANDET-09f`) and Track C are deferred experiments.

The steps below stay valid as the deferred design.
Tracker rows: PLAN.md `ANELLO-PRIOR.1–.3`, `REL-16`, `ANDET-09f`.

Two tracks run in parallel:
- **Track A** (critical path, mostly local) makes read-level evidence available
  for a *default* kallisto call and measures it on the counted molecules.
- **Track B** (Slurm) recalibrates the STAR branch's filters.
- **Track C** lists the decisions only the user can make.

## Already done (c8bc8c7)
- TSO-in-R1 audit closed. The 6 TSO "barcodes" have 0 records in the corrected
  BUS and are not in the count matrix, so they make up 0 of the 57,715 UMI.
- Correction recorded: the 19,785-read census set came from the July evidence
  job, built from the uncorrected BUS plus the raw merged FASTQ. It is not the
  counted set, so the 81/62 bound bounds that set only.
- Checked: today's `evidence_run.py:75` replays `kb_r1`/`kb_r2`, the
  host-filtered input that kb counted. No bug there.

---

## Track A — critical path

### A1. `REL-16`: resolve kallisto and bustools through kb_python, not PATH (≈ 3 h)
**Why first:** `viralscan evidence` segfaults when conda's kallisto is first on
PATH, and A2 and A3 both go through it. The bustools calls in `multimap.py` have
the same exposure: the two bustools binaries also differ (1,064,040 B vs
1,160,504 B).

Changes:
1. One resolver, `tool_binary(name)`, in an existing module (`validation.py`,
   next to `FULL_TOOLS`).
   - It returns `kb_python.config.get_kallisto_binary_path()` or
     `get_bustools_binary_path()`, which are kb's own helpers and respect
     `--kallisto`/`--bustools` overrides.
   - It falls back to PATH with a warning only if kb_python cannot be imported.
2. Route every direct call through it:
   - `evidence.py:392` (preflight `have_tools`), `:400-416` (`kallisto bus`),
     `:418-437` (`bustools capture/sort/text`);
   - `scripts/multimap.py:142-174` (`bustools correct/sort/text`).
3. `doctor_report()` shows the resolved path, not `shutil.which`.
4. Provenance: record the resolved path and md5 per binary.
   - In `evidence_manifest.json`: no fingerprint, so this is safe.
   - In `run_manifest.json`: under a new key **excluded from
     `run_fingerprint`**, as `completion_marker` is (`run_safety.py:110`).
     Otherwise every `--resume` against an older manifest breaks.
5. Tests:
   - unit test the resolver with kb_python's config patched;
   - an evidence preflight test;
   - a `--resume` fingerprint-stability test.
   - **Trap:** the integration tests (`test_exact_lineage.py`,
     `test_tiny_end_to_end.py`) build their index with PATH `kallisto index`.
     They must build it with the resolved binary, or the bundled kallisto meets
     a conda-built index and spins (the 787 % CPU case, `.living/learnings.md`).
6. Acceptance: the original repro passes — `viralscan evidence` on the cat42b
   `panel.idx`, with conda kallisto deliberately first on PATH.

### A2. Re-run the census on the counted reads (≈ half a day including queue)
Gives the bound that can be compared against the 57,715 UMI.
1. Slurm: run `viralscan evidence --virus Alphatorquevirus --viral-fasta <panel>
   --host-fasta <host>` on `covid_viralscan/results_hostfilter/LUM-SJ-x213-g`.
   - 8 CPUs, partition `all`, code from a pinned `vs_pinned/<sha>` checkout.
   - The replay reads `host_filtered/R{1,2}` (19 GB).
2. **Version trap, which gates step 3.**
   - The counted run used kallisto **0.51.1** (kb 0.29.5, `run_info.json`). The
     resolver gives 0.52.0.
   - Read-number lineage is only valid if the replay reproduces the count.
   - Pre-registered parity gate: the replay's Alphatorquevirus (CB, UMI) set
     must equal the counted set, which is 57,715.
   - If it doesn't, replay with a 0.51.1 binary (kb 0.29.5 env), passed
     explicitly; don't patch around it.
3. Run `scripts/anello_body_census.py` on the evidence BAM's anellovirus
   primary records.
   - Report k/N in reads and in **counted** CB+UMI.
   - Use the pre-registered label from the docstring as the bound; the
     post-hoc line is labelled as such.
4. Update F-019, quoting this bound, and only this one, against the 57,715.

### A3. Read measures in `viralscan evidence` (≈ 2 h)
This closes the "columns only with `--anello-align`" row for any kallisto call,
through the documented follow-up command.
1. In `_alignment_qc_from_text` (`evidence.py:601-718`), per reference:
   `complex_body_fraction` and `reagent_fraction`.
   - They come from `Alignment.read_seq()`, `has_reagent` and
     `is_complex_body(seq, query_span())`.
   - **Not** `accession_metrics`: minimap2 has no NH tag, so `weighted_reads`
     and `breadth_unique` would be wrong.
2. `r1_tso_fraction` is left out here, with a one-line note in the docs.
   Exact-lineage reads all carry corrected on-list barcodes, so the column would
   always read ~0.
3. Docs: the `alignment_qc.tsv` row in `output_reference.md`
   ("`viralscan evidence` output", line ~546), then re-pin the doc in the
   inventory.
4. Tests: pure text-parser tests with fixture SAM lines, forward and reverse.
   Reuse the real covid reads from `tests/test_anello_align.py`.
5. Automatic evidence on every default run, as a Snakefile rule, is **out of
   scope**. That would be a separate decision (Track C).

---

## Track B — `ANDET-09f` filter recalibration (Slurm, parallel to A)

### B1. Regenerate the plant (≈ 20 min job)
The committed set predates the poly-A tail fix: its `truth.tsv` has no
`body_len`/`tail_len`, and its manifest has no `polya_tail_on_3p`.
- Same arguments as job 25696086; recover them from `sacct`/the log.
- Pinned checkout at current HEAD.
- Write to a new directory (`runs_anello_plant_v2/`) so the old set stays
  comparable.

### B2. Baseline at the current filters (≈ 20 min array + 15 min report)
- Run `anello_acceptance.sbatch` (3 arms, 8 CPUs, 96 G) and then the report,
  **without** `--reuse-per-genome`.
- Expect about ¼ of 3′ reads to be lost to the 0.80 coverage filter.
- **Re-check C3:** A-rich tails may now be lost to the host filter.

### B3. Code for the sweep (≈ 3 h)
1. `align_cmd(..., overrides=None)`: an override mapping used by the report
   only. No user-facing config key until a default actually changes.
2. The sweep must move all **three** coverage filters.
   `--outFilterScoreMinOverLread 0.80` caps the other two otherwise.
3. The report reads `body_len` from `truth.tsv` and reports 3′ recovery per
   `body_len` bin (0–19, 20–39, 40–59, 60–90).
4. A sweep mode in `anello_acceptance_report.py`:
   - per setting, re-run the per-genome STAR sets (seconds each) and reuse the
     kallisto per-genome results unchanged;
   - the negative arm is STAR-only on `runs/negative/*/host_filtered/`, not a
     full viralscan run (≈ 11 min each).
5. Tests: bin assignment, override splicing, and the sweep table schema.

### B4. Pre-register the grid and decision rule in PLAN before submitting
**User confirms these.** Proposed:
- `MatchNminOverLread = ScoreMinOverLread` ∈ {0.80, 0.66, 0.50, 0.33}
- `MismatchNoverLmax` ∈ {0.08, 0.10, 0.15}
- That makes 12 settings.
- **Rule:** take the most permissive setting with **0 negative-arm molecules**
  that also passes C3. On ties, take the highest 3′ recovery in the 20–39 nt
  `body_len` bin. Then re-score C1 on the `uniform` window at that setting.

### B5. Run the sweep and report (≈ 1–2 h job)
Output: a table of setting × {C1 per genome, 3′ recovery by bin, negative
molecules, C3}.

### B6. Decide (user)
If the chosen setting passes C1 on all 8 genomes:
- revisit `--anello-align` default-off under the pre-registered `ANDET-09e`
  rule;
- set the new filter values in `ALIGN_ARGS`.

---

## Track C — decisions for the user
1. **Re-quantification of covid and SFL tonsil.** Recommendation: **defer
   until B6.**
   - At today's 0.80 filters the STAR branch gives about 0 anellovirus
     molecules on covid, so its columns would come out empty, i.e. not
     measured.
   - Read-level validation of the 57,715 comes from A2, not from a re-quant.
   - When it runs: STAR branch only, on the existing `host_filtered/` (cheap).
     Don't do a fresh full run.
2. **`ANELLO-PRIOR.1`, a real TTV-positive dataset.** The B1 plant already is
   the held-out planted control. Proposal:
   - a bounded search (≤ 1 h; GEO/SRA scRNA-seq of transplant or
     immunosuppressed cohorts);
   - stop when nothing usable turns up;
   - fall back to SRR32170409 (bulk lung-transplant donor, 2,659 anellovirus
     reads) as an orthogonal check.
   - This blocks nothing.
3. **`ANELLO-PRIOR.2`, k-mer capture on CDS vs CDS+UTR.** No UTR annotation
   exists, so it needs a method choice:
   - (a) derive UTRs heuristically (5′ = nt 1–513; polyA = first AWTAAA after
     the last CDS + 20 nt), or
   - (b) park it until real data from C2 exists.
   - Recommendation: (b).
4. **Automatic evidence on default runs** (a Snakefile rule after detection for
   anellovirus calls). Recommendation: decide after A3 shows the columns are
   useful.
5. **Housekeeping:** `scancel 25651950` (DependencyNeverSatisfied), and push
   when you choose.

## Deferred (no action unless triggered)
- rc(UMI)/rc(CB) in R2 vs R1 on read-through reads: only if read-through
  becomes common in a call.
- A host R2 fixture showing no forward-TruSeq starts.
- About 50 kept reads with a degraded TSO (≥ 3 mm): a known ceiling in
  F-019.

## Order and estimates
| Day | Track A | Track B |
|---|---|---|
| 1 | A1 (3 h) → submit A2 | B1 → B2 (jobs), B3 code (3 h) |
| 2 | A2 parity gate + census; A3 (2 h) | B4 sign-off → B5 → B6 decision |

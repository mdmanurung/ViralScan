# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-07-01] 64 GB references/ was untracked but NOT gitignored

**Category**: gotcha

**What happened**: During ingest of the reference set, `git check-ignore references` returned nothing — the 64 GB `references/` tree (STARsolo genome dirs) was untracked but not ignored. A stray `git add -A` / `git add references` would have tried to stage 64 GB.

**Why it matters**: Accidentally staging/committing multi-GB genome indices bloats the repo irreversibly (git history keeps them forever) and can hang or OOM the commit.

**Resolution**: Added `references/`, `starsolo_p22_6/`, `starsolo_p22_6b/` to `.gitignore`; registered the data via mycelium ingest with in-place pointer doc + provenance instead of committing bytes.

**Tags**: git, large-files, gitignore, ingest, bioinformatics, references

**mitigation_type**: convention

**structural_mitigation_candidate**: A pre-commit hook rejecting staged files > ~50 MB would structurally catch this class of error; not yet shipped.

### [2026-07-02] Controlling for a composite covariate that CONTAINS the tested feature is circular

**Category**: gotcha

**What happened**: A verification review flagged that the %mito control for MT-ND4L was circular — %mito was computed from a mito-gene set that INCLUDED MT-ND4L, so regressing MT-ND4L on %mito partly regresses it on itself (self-suppression), which would drop its significance even absent any real confound.

**Why it matters**: Composite covariates (%mito, %ribo, total counts, module scores, a cell-type signature) often contain the very feature being tested. Controlling for them then absorbs the feature's own variance and manufactures a null. The verdict is uninterpretable until the tested feature is excluded from the covariate.

**Resolution**: Recomputed %mito with MT-ND4L excluded (leave-one-out); MT-ND4L still dropped (FDR 4e-7→0.067), so the mito-QC confound is real and now non-circular. Fixed `go_enrichment.py`. General fix: exclude the tested feature from any composite covariate (leave-one-out) before adjustment.

**Tags**: confounding, circularity, covariate, scrna-seq, mito, review, statistics

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: When a covariate is a sum/score over features, build it leave-one-out per tested feature (or assert the tested feature ∉ the covariate's inputs).

### [2026-07-02] Mitochondrial genes need a %mito control before claiming them as biology

**Category**: gotcha

**What happened**: MT-ND4L was one of the 5 "depth-robust" EBV host-response genes (depth-adjusted FDR 4e-7). Adding a percent-mitochondrial covariate dropped it to FDR 0.07 — it was a mitochondrial-content/QC effect, not an EBV response. Across the full set, 217 of 518 depth-robust genes were lost once %mito was controlled.

**Why it matters**: In scRNA-seq, %mito tracks cell stress/quality and correlates with many conditions; a mitochondrial gene surviving depth adjustment can still be a QC artifact. Depth control is not enough — %mito is a separate confounder for mito-encoded genes (and for stress signatures generally).

**Resolution**: Added a %mito covariate alongside log-depth; reported only the depth-AND-mito-robust set (301 genes) and dropped MT-ND4L. See `go_enrichment.py`, finding F-001.

**Tags**: confounding, mitochondrial, scrna-seq, qc, gene-selection, causal-inference

**mitigation_type**: convention

**structural_mitigation_candidate**: Convention: "for scRNA-seq association/DE, control for %mito (and cell-cycle where relevant) in addition to depth; flag any mito-encoded hit for a %mito control." Fits alongside the depth-label convention.

### [2026-07-01] Raw-count positivity thresholds silently confound with sequencing depth

**Category**: gotcha

**What happened**: The EBV+ label (≥10 raw viral UMI) turned out to be strongly depth-dependent — EBV+ rate rose 27.5%→96.3% across host-depth quintiles, and depth alone predicted the label at AUC 0.80–0.97, beating the host-gene classifier. Class balancing + top-depth filtering did NOT fix it (they match class counts, not between-class depth distributions; the filter even widened the gap).

**Why it matters**: Any "positive/detected/expressed at ≥N raw counts" label makes the label a proxy for sequencing depth, so any depth-correlated feature looks predictive. This can turn a technical artifact into a headline biological result. It nearly did here.

**Resolution**: Ran a depth-confounder check (depth-alone AUC baseline, quintile rates, E-values); reported the headline as substantially confounded; recommended a depth-normalized label. See findings F-001 (contradicted) and F-003.

**Tags**: confounding, sequencing-depth, count-data, thresholding, label-definition, scrna-seq, causal-inference

**mitigation_type**: convention

**structural_mitigation_candidate**: A convention "define positive/detected labels depth-independently (CPM/fraction or depth-matched); always report depth-alone predictive baseline for any count-threshold label." Candidate for `.living/conventions.md`.

### [2026-07-01] Feature selection (HVG) before the CV split leaks into held-out metrics

**Category**: gotcha

**What happened**: Review of the EBV host-response classifier found `sc.pp.highly_variable_genes` fit on all 1906 cells before the per-seed train/test split. Because EBV+ cells drive variance in EBV-responsive genes, the HVG feature set is chosen using cells later used for evaluation — the held-out AUC/MCC are upward-biased. Documentation ("evaluated on held-out cells") was true for rows but not feature columns.

**Why it matters**: This is the single most common silent leak in scRNA-seq classifiers, and it inflated a headline metric that was about to go into a report. Any per-cell feature selection, normalization-parameter fit, or gene selection done before the split is suspect.

**Resolution**: Flagged as Major (F1) in `.living/outputs/reviews/2026-07-01-hostresponse-ebv-matched.md`; fix is to move HVG/gene selection inside the CV loop and re-run.

**Tags**: leakage, feature-selection, hvg, scrna-seq, cross-validation, classifier, review

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: A review-checklist item (and a possible convention) "feature selection / normalization params must be fit inside the CV split"; a test that asserts HVG is recomputed per fold. Candidate for promotion to `.living/conventions.md` if it recurs.

**Update (fixed 2026-07-01)**: After moving HVG inside the CV split, the corrected metrics differed from the pre-fix values by **less than one seed-SD on every metric** — i.e. unchanged within noise. HVG is *unsupervised* (variance-based, never sees labels), so this leak class is expected to be negligible, and it was. Two lessons: (1) don't assume "remove leak ⇒ lower number" — an unsupervised-selection leak can move the number either way or not at all; (2) more importantly, don't narrate a mechanism for a sub-SD shift (an earlier draft claimed per-fold HVG "tracks the EBV contrast better" — that was story-fitting noise; corrected per advisor). Fix the leak because it's methodologically correct, not because it changes the answer here.

**Category**: insight

**What happened**: Re-ran `_run_l2_regression` (fixed seeds) in a different conda env (bioenv, sklearn 1.7.1) than the original run to add MCC. Sensitivity/specificity/balanced-accuracy and the 15-gene stable set reproduced **exactly**; AUC differed in the 4th decimal (0.84456 → 0.84487).

**Why it matters**: AUC is rank-based on continuous `predict_proba` outputs, so it's sensitive to tiny BLAS/sklearn-version numerical differences; thresholded metrics (predictions at 0.5) are not. When reporting AUC across environments, quote ≤3 decimals or pin the env, or the "same" analysis looks irreproducible at high precision.

**Resolution**: Reported AUC as 0.845 (3 dp); documented the reproduction in the analysis provenance. Adopted the re-run outputs wholesale so all metrics come from one consistent run.

**Tags**: reproducibility, auc, sklearn, numerical, metrics, bioinformatics

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Pin sklearn/BLAS versions in the analysis env and record them in ENVIRONMENTS_INSTALLATIONS.md; assert AUC reproduces to 3 dp in a regression test.

### [2026-07-01] ACTIVE_CONVENTIONS.yaml is malformed after install_convention.py

**Category**: gotcha

**What happened**: After installing convention packs, `.living/conventions/ACTIVE_CONVENTIONS.yaml` contains `active_conventions: []` followed by dangling list entries. As YAML, `active_conventions` parses as an empty list and the entries below are orphaned — so a tool that reads `active_conventions` to see what's installed would see nothing.

**Why it matters**: The ingest protocol step 3 says "read ACTIVE_CONVENTIONS.yaml to see what's installed"; if parsed literally it reports no active conventions, so domain validation (bioinformatics) could be silently skipped.

**Resolution**: Worked around by reading the file's list entries directly (bioinformatics + skill-bridge are installed). Did not patch the mycelium-generated file. Consider filing a mycelium convention-gap issue.

**Tags**: mycelium, tooling, yaml, conventions, bug

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: install_convention.py should append entries under the `active_conventions:` key (or replace the `[]`), and a yaml.safe_load round-trip assertion in its test suite would catch the malformed output.

### [2026-07-02] kb ref / kallisto index have two silent FASTA-vs-GTF contracts

**Category**: gotcha

**What happened**: The covid_viralscan Stage 2 build hit two distinct failures on the
same combined host+viral reference. (1) `kb ref` **hung for 3.5 h at 0 bytes** because
`combined.gtf` carried Ensembl *chromosomal* seqnames (`1`, `2`, `X`) while the cDNA FASTA
headers are ENST transcript IDs — ngs_tools' genome-split step scanned 1.4 GB looking for
chromosome sequences that don't exist. (2) After switching to a cDNA-level GTF (seqname =
ENST) and letting the 2 h cDNA extraction finish, the final `kallisto index` **aborted**
with `Error: repeated name in FASTA file` because the source panel `viral_genome.fa`
contained accession `NC_002076.2` (Torque teno virus 1) twice (byte-identical, present in
both the anellovirus and Serratus sets), so ngs_tools extracted its 4 gene models once per
copy → 4 duplicate names.

**Why it matters**: `kb ref` exits 0 and logs "complete" even when the kallisto index step
failed with a 0-byte `index.idx` — you MUST verify by artifact (`kallisto inspect`, non-zero
size), never by exit code. Both failure modes are silent contracts: GTF seqnames must match
FASTA headers, and every target name in the FASTA must be unique.

**Resolution**: (1) generate a cDNA-level GTF via `gen_combined_cdna_gtf.py`; (2) keep-first
dedup of `cdna.fa`/`t2g.txt`/`combined.fa` (removed records byte-identical → equals a clean
build), resume `kallisto index` directly on deduped inputs, and add a dedup guard (Step 2.6)
to `slurm_build_ref.sh` so any future duplicate accession is handled automatically. Related:
the bulk `panel.idx` (P23) was separately found host-less because its June 24 build died at
the Ensembl `current_gtf` 404 (fixed 2026-07-01) — same "verify the artifact, not the log"
lesson.

**Tags**: kb-python, kallisto, kb-ref, ngs_tools, bioinformatics, reference-build, gotcha, verify-by-artifact

### [2026-07-02] bustools correct silently needs an UNCOMPRESSED whitelist

**Category**: gotcha

**What happened**: ViralScan's Snakefile passed the ngs_tools bundled 10x v3 whitelist
(`10x_version3_whitelist.txt.gz`) straight to `kb count -w`, which hands it to
`bustools correct`. bustools 0.45.1 does NOT decompress a gzipped whitelist — it reads the
compressed bytes as barcode lines and aborts with `Error: on-list file malformed;
encountered barcode length 137 on a line but barcode length 57 on another line`. kb count
swallowed this (the Snakefile captured `kb count ... 2>&1` into a shell variable that was
then discarded), exited 0, and left `counts_unfiltered/` unwritten — so the pipeline failed
one step later at `mv counts_unfiltered/` with a misleading "No such file or directory".

**Why it matters**: Two compounding silent failures — a tool that needs plain-text input but
gives a cryptic length error on gzip, and a wrapper that discards the tool's stderr. Verify
by artifact (does `counts_unfiltered/` exist?), never by exit code, and never capture a
long-running tool's output into a variable you throw away.

**Resolution**: `kb_count` Snakefile rule now decompresses any `*.gz` whitelist to a plain
file before kb count, tees kb output to `<output>kb_count.log`, and asserts
`counts_unfiltered/` exists (printing the log tail + exit 1 otherwise). Confirmed: with the
decompressed whitelist, `bustools correct` finds all 6,794,880 barcodes. Related:
[[kb-ref-kallisto-index-fasta-gtf-contracts]] (same verify-by-artifact lesson).

**Tags**: kb-python, bustools, kallisto, whitelist, gzip, snakemake, bioinformatics, verify-by-artifact

### [2026-07-02] Two panel dirs — bulk scan defaulted to the incomplete one

**Category**: gotcha

**What happened**: There are two bulk-panel reference directories:
`viralscan_bulk_gse128078/ref/` (host-less, viral-only, 2,742 entries — the June-24 build
that died at the Ensembl 404) and `viralscan_panel_ref/ref/` (the complete host+viral build,
470,533 entries, made 2026-07-02). `scripts/bulk_viral_scan.sh` hardcoded
`REFDIR=$WORKDIR/ref` = the incomplete one, so a bulk scan would have silently quantified
against a host-less index (no host reads absorbed → skewed viral specificity).

**Why it matters**: A stale/partial reference in a plausible-looking location is a silent
correctness trap — the scan runs fine and produces numbers, they're just against the wrong
index. When a rebuild lands in a new path, grep for every consumer of the old path.

**Resolution**: Repointed `bulk_viral_scan.sh` `REFDIR` (now `REFDIR=${REFDIR:-…/viralscan_panel_ref/ref}`,
overridable). Still TODO: point `bulk_viral_summarize.py --t2g` at the new `panel.t2g` before
the B4 summarize step. Related: [[kb-ref-kallisto-index-fasta-gtf-contracts]] — the same
incomplete build is why P23.op1 had to be rerun.

**Tags**: bulk, reference-panel, paths, gotcha, gse128078, silent-correctness

### [2026-07-02] Wrong 10x whitelist silently produces an all-empty-droplet matrix

**Category**: gotcha

**What happened**: The covid_viralscan quant used the 10x v3 whitelist for a library whose
barcodes are NOT in it (likely GEM-X 5′). `bustools correct` marked 96.5% of records
"uncorrected" and kept going; the resulting matrix had 163,203 barcodes but median 1 UMI and
max 5,311 — essentially all empty droplets. CellRanger called 28,922 real cells from the same
FASTQs. Only 0.4% of raw R1 barcodes match the v3 whitelist; the best bundled ngs_tools list
(v4/GEM-X) matched just 4.7%, while the R1 barcodes matched CellRanger's called cells 56.5%.

**Why it matters**: A wrong whitelist doesn't error — it silently drops most reads and yields
a plausible-looking but meaningless single-cell matrix. Any "% of cells infected" computed on
it is noise. ALWAYS sanity-check the barcode/whitelist match (correction rate, per-barcode UMI
distribution / knee, overlap with a trusted CellRanger cell set) before per-cell claims.
Diagnostic: `p_pseudoaligned` was also low (6.4%), and the "uncorrected" fraction from
`bustools correct` (visible now that kb_count.log is teed) is the smoking gun.

**Resolution**: (pending) re-run with the correct chemistry whitelist — determine the 10x
chemistry (GEM-X 5′?) CellRanger auto-detected and pass its whitelist to `viralscan -w`, or
add proper 5′/GEM-X whitelist support. Related: [[bustools-correct-needs-uncompressed-whitelist]]
(the gzip bug that had to be fixed first to even reach this stage).

**Tags**: 10x, whitelist, barcodes, bustools, gem-x, 5-prime, scrna-seq, empty-droplets, gotcha, verify-by-artifact

## [2026-07-03] Verify agent-authored plans against HEAD (and numbers against primary artifacts) before executing

**Category**: process / gotcha

**What happened**: An agent-authored publication-readiness plan contained three classes of error
caught only by reconciling against the live repo: (1) an already-completed task (package
`emptydrops.R` + its test were already at HEAD, commit `e6b1186`); (2) a hallucinated API rewrite
(`cell_type_enrichment` to a `SimpleNamespace` signature — the real API is dict-based and `api.md`
already documented it correctly); (3) a cross-design statistical pairing — it paired the manuscript's
stale all-cell AUROC 0.845 with the depth-alone 0.967 from a *different* (balanced+depth-filtered)
design. The tracked `hostresponse_summary.tsv` already carried 0.866, and `depth_confounder.txt`
explicitly pairs 0.967 with the 0.866 headline in the same n=1906 design.

**Why it matters**: Executing the plan verbatim would have re-done work, documented a nonexistent API,
and reintroduced the exact depth-confounding error the fix was meant to remove — under a "corrected"
banner. Presence of a number in an artifact is not provenance of that number; read the same-split
field, not a coincidental match.

**Resolution**: Reconcile every task step-by-step against HEAD; verify each cited number against the
primary artifact that generated it (here `depth_confounder.txt` line 1 = n=1906, lines 13–14 =
same-design 0.866 vs 0.967). Wrote `...-reconciled.md` before touching code.

**Tags**: plan-reconciliation, verify-by-artifact, host-response, depth-confound, hallucinated-api, manuscript, process

## [2026-07-03] reference-strategy 2×2 benchmark: root causes of the 8 failed rows

**Category**: bioinformatics / environment / gotcha

**What happened**: The 4/12-complete reference-strategy benchmark (combined vs two_step ×
STARsolo vs ViralScan × 3 viruses) failed in three distinct ways:
1. **All 4 ViralScan `two_step` rows blocked** — `viralscan --host-filter kallisto` preflight
   (`_check_host_filter_tools`) needs standalone `kallisto` AND `bustools` on PATH. The benchmark
   env (`evonk/.../test_viralscan`) had `kb` + `bustools` but **no standalone `kallisto`** (the
   kb-python-bundled-kallisto PATH shim did not expose one). Fix: built
   `mdmanurung/conda/envs/viralscan_bench` from `environment.yml` (has `kallisto 0.52.0`); verified
   `_check_host_filter_tools('kallisto')` now passes.
2. **EBV STARsolo `combined` failed** — `FATAL ERROR in reads input: quality string length is not
   equal to sequence length`. Geometry was correct (10xv2, CB16/UMI10); both FASTQs pass `gzip -t`
   with matching record counts → a single malformed record STAR rejects but kallisto tolerates
   (the same sample ran fine for the manuscript §3.3 STARsolo comparison via a different copy).
   Fix = sanitize the record (seqkit), not re-fetch.
3. **Incomplete ViralScan `combined` + STARsolo `two_step` rows** — started, no final outputs;
   likely wall-time/OOM on the deep EBV sample. Re-run + inspect per-row logs.

**Why it matters**: The benchmark is excluded from manuscript claims until all 12 rows complete;
these are the exact unblock steps. Plan: `analysis/reference_strategy_benchmark/COMPLETION_PLAN.md`.

**Tags**: reference-strategy, benchmark, starsolo, viralscan, host-filter, kallisto, conda-env, fastq, gotcha, verify-by-artifact

## [2026-07-03] Scratch cleanup deleted the ViralScan kallisto index — recovered from archive; kallisto version gotcha

**Category**: environment / reproducibility / gotcha

**What happened**: Completing the reference-strategy 2×2 surfaced that the ViralScan combined
kallisto index (`index_plus_anellovirus.idx` + t2g) used by all 6 ViralScan rows had been
**deleted by scratch cleanup** (`/exports/para-lipg-hpc/.../viralscan_showcase/fullrun/refs/`
is gone) — the real reason those rows are "incomplete", beyond the missing kallisto binary.
A **persistent copy survived on archive** at the (doubled) path
`/exports/archive/.../mdmanurung/viralscan_showcase/viralscan_showcase/fullrun/refs/merged/`
— the exact Jun-22 build the Jun-28 benchmark used (t2g: 226,005 host ENST + 821 target-virus
rows). So it's a **restore, not a rebuild**; being the same index, the 4 complete rows stay valid
(only the 8 incomplete rows re-run). Full rebuild materials (GRCh38 + Serratus fasta_split +
`create_final_transcriptome.slurm` recipe; or `scripts/build_bundled_panel_ref.py`) also survive.

**Kallisto version gotcha (important)**: the **standalone conda `kallisto` 0.52 SEGFAULTS**
reading these older kb-python-built indices (combined AND host), but the **kb-python bundled
kallisto reads them** (0.51.1 and 0.52.0 bundled both OK). `kb count` uses the bundled kallisto
internally (fine), but `viralscan --host-filter kallisto` calls the *standalone* `kallisto` on
PATH — so the run harness must **prepend the kb-python bundled kallisto dir to PATH** so
`which kallisto` = bundled. Verified fix.

**Why it matters**: on this HPC, scratch (`/exports/para-lipg-hpc`) is cleaned; archive
(`/exports/archive`) persists. Keep reference indices on archive. And don't assume a newer
standalone kallisto reads an index a bundled kallisto built — verify by `kallisto inspect`.

**Tags**: kallisto, kb-python, index, scratch-cleanup, archive, version-mismatch, reference-strategy, benchmark, gotcha, verify-by-artifact

## [2026-07-04] Verify-by-artifact turned the reference-strategy benchmark inside out (twice)

**Category**: verification / gotcha

**What happened**: Two "obvious" premises about the reference-strategy 2×2 were both wrong, caught
only by checking live artifacts:
1. The EBV STARsolo "failed / quality-string-length" row had **actually completed** (Log.final.out
   "ALL DONE", 127M reads mapped, 223MB matrix) — the results TSV captured a **stale early-attempt**
   status (job 25102703 vs the successful 25102837). So the planned 11GB FASTQ "sanitize" was a
   no-op (R1=0 malformed; STAR mapped everything). The advisor's insistence on *detect before rewrite*
   saved the wasted rewrite.
2. Re-parsing with `summarize_reference_strategy.py` (reads live files) showed scratch cleanup is
   **actively deleting benchmark data**: the ViralScan index (restorable from archive), the HHV-6B +
   HSV-1 FASTQs (gone → ENA re-fetch), and 2 previously-"complete" rows' outputs. Net true state:
   4 complete, 8 to re-run — a *different* 4 than the stale TSV claimed.

**Resolution**: chose the cheap high-value slice — the **EBV 2×2** (on-target; Selectivity Index
computable; EBV FASTQs survive; index restored). Staged self-contained re-run of the 4 EBV rows into
`fresh12b` (array 4,5,6,7) + a summarize command. HHV-6B/HSV-1 deferred (need re-fetch).

**Why it matters**: on a 93%-full scratch, benchmark "completeness" is not stable — re-derive status
from live artifacts, never trust a cached results TSV; and keep reference indices/inputs on archive.

**Tags**: verify-by-artifact, reference-strategy, benchmark, starsolo, scratch-cleanup, stale-status, ebv, gotcha

## [2026-07-04] `conda activate` doesn't guarantee its python is first — auto-activated env shadows it in SLURM

**Category**: environment / gotcha

**What happened**: The EBV 2×2 SLURM array (job 25144701) failed in 7s — all 4 tasks —
with `ModuleNotFoundError: No module named 'kb_python'` at the run script's kb-python shim.
Cause: the login profile auto-activates the `codex` conda env, and even after
`conda activate viralscan_bench` in the batch script, `python` still resolved to
`codex/bin/python` (Python 3.14, no kb_python). `set -e` + the failed `$(python -c import kb_python)`
substitution aborted the job before any real work. (Reproduced interactively: `which python`
after activate = codex, not viralscan_bench.)

**Fix**: immediately after `conda activate <env>`, force the env's bin to the front:
`export PATH="/exports/archive/.../conda/envs/viralscan_bench/bin:$PATH"`. Then `python` =
env python (kb_python present) and the bundled-kallisto shim works. Re-submitted as 25144705.

**Why it matters**: on this HPC, don't trust `conda activate` alone in non-interactive/SLURM
shells when a profile auto-activates another env — explicitly prepend the target env's bin, or
use absolute paths to the env's `python`.

**Tags**: conda, slurm, path-shadowing, kb-python, environment, batch, gotcha

## [2026-07-04] Benchmark "ViralScan ≫ STARsolo" headline was mostly a GTF-annotation artifact + a regex naming bug

**Category**: verification / gotcha

**What happened**: The completed reference-strategy 2×2 initially looked like ViralScan is
~18× (EBV), ~3.7× (HHV-6B), ∞ (HSV-1 STAR=0) more sensitive than STARsolo. Fair-comparison
harmonization (anchor-restricted, unique-layer) collapsed all of it except HHV-6B (~1.9×):
- **HSV-1 STAR=0 was a summarizer regex bug**: STARsolo names HSV-1 `HHV1gp…`; the HSV-1
  `target_regex` omitted the `hhv-?1` alias (EBV had `hhv-?4`, HHV-6B `hhv-?6b`), so it matched
  0 features. Fixed in `src/viralscan/reference_strategy.py` + regression test
  (`test_reference_strategy_benchmark.py`). Real STARsolo HSV-1 ≈ 19 UMI (still a GTF artifact).
- **EBV "2×" is reference-completeness**: 95% of ViralScan's EBV UMI comes from CDS-only genes
  STARsolo's combined GTF cannot count — the two tools measure near-disjoint gene sets.
- **Confirmed**: reference-strategy axis (combined vs two_step) is minor; multimap adds +51–122%.

**Why it matters**: only HHV-6B gives a defensible aligner comparison. Keep the benchmark
excluded from the manuscript, or present HHV-6B-only with the GTF-artifact caveat. Never trust
a benchmark's headline magnitudes before harmonizing count-layer + denominator + feature-naming.

**Tags**: reference-strategy, benchmark, starsolo, viralscan, gtf-artifact, regex, harmonization, verify-by-artifact, gotcha

## [2026-07-04] Multimap main-pass speedup: EC-precompute wins ~15%; CSR fancy-index gather BACKFIRES

**Category**: performance / gotcha

**What happened**: Optimizing `build_multimap_layers` (the ~42-min/method single pass over ~103M
BUS records), verified against a golden snapshot (6 seeds × 4 methods × 8 layers, exact 0.00e+00):
- **Hoisting per-EC invariants** (gene classification, conservative/selected masks) out of the
  per-record loop + iterating column arrays via `zip` instead of `itertuples`: **~15% faster**
  (50.4s vs 58.9s on 600k synthetic records), byte-identical. Committed.
- **Replacing the per-gene `_matrix_value` scalar lookups with a per-record CSR fancy-index gather**
  (`original_counts[cell, gene_list]`) **made it 1.6× SLOWER** (96.8s) — CSR is row-oriented and each
  fancy index builds a new sparse object. Reverted. The benchmark caught it; a blind "obvious"
  vectorization regressed.

**Why it matters**: (1) always benchmark a perf change against the unoptimized baseline — the
intuitive vectorization was slower. (2) The remaining hotspot (the `_matrix_value` weights lookup +
CSR triplet construction) needs the profiler's cProfile attribution (pending) to target, and any
bulk `original_counts` gather must be validated for speed, not assumed. Golden-equivalence harness:
/tmp/multimap_equiv.py (exact-match gate for any further rewrite).

**Tags**: multimap, performance, csr, sparse, benchmark, gotcha, verify-by-artifact

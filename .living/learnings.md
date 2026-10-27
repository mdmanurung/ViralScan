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

## [2026-07-15] Coverage breadth is the decisive per-run quality gate — not optional post-hoc

**Category**: convention

**What happened**: After three convergent methods showed no genuine viral infection in the COVID
scRNA-seq samples, the coverage-breadth step (`viralscan evidence` / minimap2 → `samtools
coverage`) was the only measure that could rule out artifact versus real infection in the host-
filtered residual. UMI counts (even post-STAR-filter) cannot make this distinction — a fixed
host-homologous locus can accumulate thousands of UMI while covering <4% of the genome.

**Why it matters**: A viral signal that is deep-but-narrow (high depth, low breadth) is
definitionally artifact (fixed locus). Real infection spreads reads across the genome. This
distinction is invisible from count matrices alone and only visible from alignment-level
breadth metrics. Treating breadth as optional/confirmatory rather than a standard output tier
means artifact calls can appear in results without any path to refutation.

**Resolution**: Established as the three-tier reporting standard: (1) UMI count (kallisto), (2)
STAR host-filter, (3) `viralscan evidence` breadth. All three tiers required before concluding
"no genuine infection" or "genuine detection." See review 2026-07-15 finding F3 for the
reproducibility gap introduced when `viralscan evidence` crashed and breadth was obtained from
an undocumented remediation job.

**Tags**: coverage-breadth, viralscan-evidence, artifact-detection, anellovirus, host-homology, convention, viralscan

## [2026-07-15] Undocumented remediation jobs produce authoritative results that cannot be reproduced from the workflow

**Category**: process / reproducibility / gotcha

**What happened**: `viralscan evidence` (job 25180994) crashed at `samtools sort` with a
duplicate BAM header entry (`NC_002076.2`). A manually crafted `hf_align` job (25181135) was
submitted without committing its script to `covid_viralscan/scripts/`. The resulting
`coverage.tsv` is the authoritative breadth table cited in the findings, but the exact command
and parameters used to produce it exist only in SLURM history — not in the repository.

**Why it matters**: Any manuscript citation of the breadth numbers (max 3.41%, NC_001479.1
1.99% saturation) references a result that cannot be reproduced from the committed codebase.
This creates a reproducibility gap at exactly the most decisive evidence tier. The original
pipeline crash was caused by a fixable upstream issue (duplicate NC_002076.2 in the reference
FASTA); fixing that and committing the script closes the gap.

**Resolution** (pending): (1) dedup NC_002076.2 from the viral reference FASTA/GTF; (2)
commit the hf_align script to `covid_viralscan/scripts/`; (3) add a note to RUNBOOK.md;
(4) re-run `viralscan evidence` from the documented pipeline to regenerate coverage.tsv from
a reproducible workflow. See review 2026-07-15 finding F3.

**Tags**: reproducibility, undocumented-job, remediation, samtools, NC_002076.2, viralscan-evidence, coverage, gotcha

## [2026-07-15] Internal impossible-detection controls calibrate the noise floor

**Category**: convention

**What happened**: The COVID ViralScan survey detected variola (VARV/smallpox) at 1 UMI in
x213 and 0 UMI in x216 post-STAR-filter. Smallpox was declared eradicated in 1980; any
detection is definitionally noise. Using VARV as an anchor revealed that every other
non-anellovirus signal (HHV-6B 2 UMI, CeHV2 1–1.5 UMI, MPXV 0–4 UMI, molluscum 2–3 UMI,
EBV 3–4 UMI) falls within 0–4× of this known-impossible baseline — all collapse to noise.

**Why it matters**: Panels covering eradicated, geographically isolated, or otherwise
impossible organisms (VARV, PERV, cetacean viruses in human samples) provide free noise-floor
calibrators. A signal within 1–2 orders of magnitude of a known-impossible detection cannot
be called positive without extraordinary evidence. This argument is stronger and more direct
than a statistical threshold alone, because it doesn't require knowing the threshold a priori.

**Resolution**: Documented as the per-virus-plausibility verdict in review 2026-07-15. Add
eradicated/impossible organisms to standard panel curation so any run has at least one
internal negative control.

**Tags**: noise-floor, internal-control, varv, smallpox, eradicated, per-virus-plausibility, convention, viralscan

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

## [2026-07-05] Scalar CSR lookups in a hot loop dominate build_multimap_layers (86.4% of CPU)

**Category**: performance / profiling finding

**What happened**: cProfile on a 1M-row subsample of the real EBV BUS file (103M rows) found that
`_matrix_value` at multimapping.py line 292 accounts for **86.4% of total build time** across all
non-EM methods (and 83.2% for EM). The function is called ~5.9M times per 1M BUS rows (avg 7.2
calls per multi-EC record) and each call traverses 5+ scipy dispatch layers:
`__getitem__ → _validate_indices → isintlike (×4) → _get_intXint → get_csr_submatrix`.

The EM extra cost (post-commit 2c2e6f0 vectorised em_gene_abundances) is now negligible — only 9.2s
(2.0% of equal total), all from the em_records allocation loop. em_gene_abundances is below the
profiling noise floor entirely.

**Why it matters**: The "obvious" vectorization (batch CSR fancy-index gather at commit b7e9635)
was previously tried and reverted as 1.6× SLOWER. The correct fix is different: for each BUS record
with genes_in_ec = [g1, g2, ..., gN], fetch a dense row-slice from `original_counts` ONCE
(`original_counts.getrow(cell_idx).toarray()` or precomputing `original_counts.toarray()` if memory
allows) rather than N individual scalar lookups via `__getitem__`. This eliminates the dispatch chain
for each gene without building a new sparse object per record.

**Resolution** (pending — flagged in profiling report; no src/ changes this session):
Fix at multimapping.py line 292: replace `[_matrix_value(original_counts, cell_idx, gid) for gid
in genes_in_ec]` with a single vectorised fetch of the cell's row, then index it with `genes_in_ec`.
Expected speedup: ~8× on the main pass, bringing full-data 'equal' from ~17,000s to ~2,000s.

**Tags**: multimap, performance, csr, sparse, profiling, bottleneck, scipy, bioinformatics

## [2026-07-04] Multimap main-pass: profiler-guided direct CSR buffer access = ~3× (fancy-index was wrong mechanism)

**Category**: performance / win

**What happened**: cProfile (1M rows) confirmed `_matrix_value` (per-gene `matrix[cell,gene]`) is
**86% of the pass** — the cost is scipy's `__getitem__ → _validate_indices → isintlike → get_csr_submatrix`
dispatch, run ~5.9M times per 1M rows. The fix that WORKS: read the CSR row buffers directly
(`indptr`/`indices`/`data`) and locate genes with `np.searchsorted` on the sorted row indices —
**~2.9× faster** (14.1s vs 40.5s on 600k synthetic records), **byte-identical** (golden 0.00e+00,
21 tests pass). Committed.

**Key insight**: the intuitive "batch fetch `original_counts[cell, gene_list]`" (CSR fancy index)
does NOT help — it goes through the SAME `__getitem__` dispatch and was 1.6× *slower* (measured
earlier). The win comes from bypassing `__getitem__` entirely via the raw buffers. The subagent's
cProfile correctly identified the hotspot but its recommended fix (fancy index) would have regressed;
direct-buffer + searchsorted is the correct mechanism.

**Combined multimap speedups this session** (all byte-identical): vectorized EM (sparse mat-vec;
dropped EM from cProfile top-40 to ~2%), EC-precompute/array-iteration (~15%), direct CSR access (~3×).
Net main-pass ~4× vs the original scalar-lookup loop. Golden gate: /tmp/multimap_equiv.py.

**Tags**: multimap, performance, csr, searchsorted, profiling, win, verify-by-artifact

## [2026-07-06] Archive and HPC ViralScan paths are the same inode — cp fails with "same file"

**Category**: environment / gotcha

**What happened**: Tried to copy an edited script from `/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan/` to `/exports/para-lipg-hpc/mdmanurung/ViralScan/`. Got "cp: ... are the same file" — both paths resolve to the same underlying inode. The HPC path is a symlink (or mount alias) of the archive path.

**Why it matters**: Any edit to a file at the archive path is immediately visible at the HPC path without any copy. Attempting `cp archive-path hpc-path` always fails. Scripts submitted from the HPC path see changes made via the archive path instantly.

**Resolution**: Only one edit is ever needed; `cp` between the two paths must never be attempted. Confirmed by running `ls -li` on both paths — same inode number.

**Tags**: filesystem, symlink, hpc, archive, gotcha, environment

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: None needed — the single-inode relationship is a feature, not a bug. Just know not to `cp` between the two paths.

## [2026-07-06] diag_viral_read_origin.sh: pre-built combined index unlocks parallel submission

**Category**: process / performance

**What happened**: The read-origin diagnostic script originally depended on the STARsolo re-run output (`genome_GRCh38_viral`), which takes ~30 min to build. Discovered that `references/starsolo/combined_GRCh38_2024A_serratus_plus_anellovirus/` already exists with a valid `SAindex` and covers all contigs needed for the NH-flag test (GRCh38 + anellovirus panel). SARS-CoV-2 absence is irrelevant for that test.

**Why it matters**: A hard dependency that turns out to be already satisfied lets two cluster jobs run in parallel instead of sequentially — a 30+ min wall-clock saving on a cluster with variable queue times.

**Resolution**: Repointed `GENOME` in `diag_viral_read_origin.sh` to the pre-built index. Jobs 25149332 (STARsolo build) and 25149333 (read-origin) submitted concurrently as a result.

**Tags**: starsolo, read-origin, cluster, parallel, dependency, covid, performance, bioinformatics

**mitigation_type**: awareness

**structural_mitigation_candidate**: When scripting cluster pipelines, check whether prerequisite artifacts already exist before adding a hard dependency step — saves queue latency that often dwarfs the computation itself.

## [2026-07-05] mycelium Stop-hook only checks learnings/decisions/conventions/findings mtimes

**Category**: tooling / process

**What happened**: A read-only status session (and a follow-up commit-only session) kept getting
STOP BLOCKED with "N files changed but .living/ not updated" even after I updated
`.living/last-session.md` and the session logs. Reading the hook script
(`skills/core/hooks/mycelium-stop-check.sh`) showed why: the block gate stats **only**
`learnings.md`, `decisions.md`, `conventions.md`, and the `findings/` dir, and passes only if one
of their mtimes is newer than the work-reminder timestamp. `last-session.md`, `LOG_REGISTRY.md`,
and the per-session `log/*.md` files are auto-finalized by the hook itself and do **not** count
toward the gate. It also debounces for 5 min after first work, so a quick session may block only
on a later stop.

**How to apply**: To clear a legitimate STOP BLOCK, append a real entry to one of
learnings/decisions/conventions/findings (touching its mtime) — not last-session.md. If the session
genuinely produced nothing worth triaging (pure status/commit), the honest move is still a short
learnings/decisions note (as here); there is no "skip triage" path once the block fires other than
`stop_hook_active`. See [[verify-agent-plans-against-head-before-executing]].

---

### [2026-07-07] A ViralScan quant "redo" silently reuses cached outputs; and how to re-run only the detection tail

**Category**: gotcha / operational

**What happened**: Re-running the covid `slurm_viralscan_quant.sh` to "redo" the analysis
completed in **1 second** (exit 0) and changed nothing. The script's Step-2 (`if merged R1/R2
exist → skip merge`) and Step-3 (`if $RESULTS/$S/log/kb.done exists → skip viralscan quant`)
sentinels short-circuited, because the merged FASTQs (`covid_viralscan/data/$S/`) and the
`kb.done` sentinel had survived on non-scratch (archive) storage even though the raw FASTQs on
scratch were swept. The re-extracted 259 GB tarball was never used. A true recompute required
clearing both: `rm data/$S/*_merged_R*.fastq.gz` and moving `results/$S` aside (which removes the
`kb.done` sentinel). After that the recompute ran for real (~1.5 h) and reproduced identical
numbers (deterministic pipeline) plus the current-code schema.

**Why it matters**: "verify by artifact, not exit code" applies to *re-runs* too — a 1-second
"COMPLETED" quant is the tell that sentinels skipped the work. Check `sacct Elapsed`, not just
`State`, to distinguish a real run from a cached no-op.

**How to change the cell-calling denominator — the WRONG way and the right ways.**
The Snakemake DAG is create_config → kb_count → analysis → multimap → detection; `config.yaml`
is the sole source of `cell_calling`/`called_cells_file`, written only by `create_config`.

*WRONG (tried 2026-07-07, FAILED):* remove `config.yaml` + downstream sentinels but keep
`log/kb.done` + `kb-python/counts_unfiltered/`, expecting kb_count to be reused. It is NOT —
regenerating `config.yaml` (newer mtime) cascades and **re-runs `kb_count`** (the log shows
"Starting to quantify the data…"). Worse, the kb_count rule's `mv counts_unfiltered/
kb-python/counts_unfiltered` step is **non-idempotent**: it dies `mv: cannot move … File exists`
because `kb-python/counts_unfiltered/` already exists — after wasting ~1 h re-pseudoaligning.
Both covid emptyDrops detection re-runs (25167140) failed exactly here.

*RIGHT (both used successfully):*
1. **Don't re-run the pipeline at all.** `per_cell_viral.tsv` (columns `barcode, virus_name,
   viral_umi, …`) is **cell-calling-independent** — cell-calling only sets the denominator. To
   get viral counts over ANY called-cell set, intersect its barcodes with `per_cell_viral.tsv`
   and aggregate (sum `viral_umi`, count distinct barcodes). This is how the knee / emptyDrops
   (30,792 / 15,998) / CellRanger (28,922 / 19,183) denominators were all reported.
2. **If you need the native `viral_summary.tsv` for a new cell-calling**, run a fresh
   `viralscan quant` to a NEW `-o` dir (e.g. `results_genomic/`) — kb_count runs cleanly there
   (no `mv` collision). Costs a full kb count but is reliable.

**Barcode-format note**: CellRanger `sample_filtered_feature_bc_matrix` barcodes carry a `-1`
suffix; strip it (`sed 's/-[0-9]*$//'`) to match ViralScan's kb barcodes. Overlap with the
ViralScan matrix was 100% (x213-g) and 89.7% (x216-g); the ~2 k missing x216-g CR cells had no
pseudoaligned reads (0-viral by definition). Feed the stripped list via `--called-cells-file`.

**Tags**: viralscan, quant, redo, sentinel, verify-by-artifact, snakemake, cell-calling, emptydrops, cellranger, barcodes, covid

---

### [2026-07-15]

**Category**: Infrastructure / SLURM
**Tags**: slurm, cluster, partition, qos, cpu-limit, medium-partition

The `medium` partition uses `restrictmedium` QOS which imposes `MaxCPUsPU=4` across all jobs
running under that QOS, regardless of how many CPUs each individual job requests. This means a
user already running 4 CPUs anywhere on `medium` can't submit another job there, even 1-CPU.
The `all` partition has no QOS restriction (30-day time limit cap) and bypasses this entirely.
For jobs needing 8+ CPUs, submit to `--partition=all` to avoid `QOSMaxCpuPerUserLimit` failures.

**Mitigation**: use `--partition=all` for CPU-intensive jobs (minimap2, blastn, multi-threaded
bioinformatics). Reserve `medium` for lightweight/interactive jobs where the 4-CPU cap suffices.

**mitigation_type**: process
**structural_mitigation_candidate**: false

---

### [2026-07-15]

**Category**: Biology / EVE conceptual framework
**Tags**: anellovirus, eve, endogenous-viral-elements, transposable-elements, host-homology, covid

Anellovirus sequences detected at deep/narrow coverage loci are NOT transposable elements — they
are Endogenous Viral Elements (EVEs). Key distinction: TEs have transposition machinery (retro-
transposons use RT + integrase; DNA transposons use transposase). Anelloviruses encode ORF1
(capsid), ORF2 (phosphatase-like), ORF3 (CLLD7-like) — no integrase, no RT, no transposase.
They cannot self-transpose. Ancient integrations occur via host repair pathways (NHEJ after DSB).
Belyi et al. 2010 (PLOS Pathog) documented anellovirus-related EVEs in mammalian genomes.
ERVs are a special case of EVEs that retained retrotransposition machinery.

The 156-base fixed locus seen in NC_001479.1 (EMCV, 841–1621x depth, same bases in both samples)
is consistent with a single ancient EVE integration transcribed as part of a host transcript
(intronic/UTR region), not captured by kallisto's cDNA-only reference.

**mitigation_type**: conceptual-clarification
**structural_mitigation_candidate**: false

---

### [2026-07-15]

**Category**: Infrastructure / skill-pack installation
**Tags**: mycelium, convention-pack, skill-installation, aifi, scrna

When `skills/core/scripts/install_convention.py` doesn't exist (script missing from repo),
install a skill pack manually:
1. Unzip pack: `unzip <pack>.zip -d /tmp/install/`
2. Copy: `cp -r /tmp/install/<name>/ .living/conventions/<name>/`
3. Append entry to `.living/conventions/ACTIVE_CONVENTIONS.yaml`
4. Add bullet to `## Installed Convention Packs` in `CLAUDE.md`

Skill packs with `SKILL.md` (not `analysis-conventions.md`) use `SKILL.md` as the entry point;
reference it in CLAUDE.md as `See .living/conventions/<name>/SKILL.md`.

**mitigation_type**: process
**structural_mitigation_candidate**: false

---

### [2026-07-15]

**Category**: Biology / EVE artifact mechanism
**Tags**: anellovirus, eve, plasma-cells, cell-type-enrichment, celltypist, host-homology, covid

CellTypist enrichment of host-filtered anellovirus reads (covid PBMC samples) shows significant
enrichment in both **Epithelial cells** (OR 3.3–4.2) and **Plasma cells** (OR 3.1, FDR 5e-13).
The plasma cell enrichment is mechanistically informative: plasma cells have extreme transcriptional
output (antibody heavy/light chains, secretory pathway genes) which generates abundant intronic
pre-mRNA across expressed loci. More intronic pre-mRNA → more reads from EVE-bearing introns →
more anellovirus-assigned reads passing the cDNA-only host filter.

This is a testable prediction: if the enrichment is EVE-driven, it should co-localise with
intronic reads from EVE-bearing genes (NALCN, LINC02742, etc.) not with unique anellovirus k-mers.
The CellTypist pattern adds cell-type-level confirmation that the artifact follows host-gene
expression (not viral tropism).

**Mitigation**: for any future cell-type enrichment analysis, a positive result for EVE-risk viruses
(eve_risk=True in viral_summary.tsv) should be cross-checked against (a) accession_breadth > 0.1
and (b) absence of plasma/epithelial bias before being interpreted as genuine viral tropism.

**mitigation_type**: validation-check
**structural_mitigation_candidate**: true

---

### [2026-07-17]

**Category**: Biology / EVE artifact mechanism — accession-level confirmation
**Tags**: anellovirus, eve, blast, phase-a-b, grchr38-homology, aav2, emcv, covid, panel-screen

EVE accession screen (job 25237061) provided accession-level mechanistic confirmation for the
anellovirus artifact. Key insights:

1. **Imperfect homology vs exact integration**: All 8 detected anellovirus accessions appear in
   Phase A (multi-chromosomal GRCh38 alignment) but NONE in Phase B (BLAST 100% identity). This
   is why `--genome-dlist` removed only ~15% (exact-k-mer masking), while STAR mismatch-alignment
   removed ~95%. The host↔anellovirus homology is imperfect — not exact EVE integration.

2. **NC_001479.1 (EMCV-like) is Phase B confirmed**: positions 120–303 are 100% identical to a
   human intergenic region (e=3.40e-89). This is the same accession showing depth-doubling without
   breadth increase in the coverage-breadth analysis — both methods point to a single fixed
   host-homologous locus. BLAST provides the sequence-level proof.

3. **Phase C panel-wide screen is clean**: HHV-1, EBV, HHV-6B, CeHV2, MPXV, Molluscum, and
   SARS-CoV-2 do NOT align to GRCh38. If any of these were detected, the detection could be
   trusted as not EVE-driven (subject to other artifact checks). This is useful for future cohorts.

4. **Anellovirus EVE mechanism is dispersed-homology, not single-locus integration**: the
   multi-chromosomal Phase A pattern (4–24 chromosomes, many genes) reflects sequence similarity
   with non-coding regions scattered genome-wide, not a single proviral insertion. This explains
   why the breadth-per-contig is low (reads come from many host loci → no single contig saturates).

**mitigation_type**: validation-check
**structural_mitigation_candidate**: false

### [2026-07-17]

**Category**: Packaging / CI — undeclared-dependency blind spot
**Tags**: packaging, ci, dependencies, anndata, no-deps, pyproject, environment-yml, release-hygiene

A publication-readiness review found `anndata` used as an eager top-level import
(`scripts/multimap.py`) but declared nowhere in `pyproject.toml` — satisfied only transitively via
`scanpy`. It never surfaced because **CI installs the package with `pip install --no-deps -e .`**
(the workaround for the snakemake `connection_pool`/setuptools build failure). `--no-deps` means CI
can never detect a missing runtime dependency: the import resolves because the test env already has
the transitive dep, so a green suite is not evidence the declared dependency surface is correct.

**How to apply**: when a project installs with `--no-deps` in CI, add a separate, cheap check that
the *declared* deps are complete — e.g. a metadata-only job that pip-installs the built wheel into a
bare venv and imports the package, or an import-linter/`deptry` pass. Also keep `environment.yml` and
`pyproject.toml` dependency sets in sync (this repo's env was additionally missing `scikit-learn`).
General pattern: a passing test suite validates behavior under the *test* environment, not the
*declared* install contract — those are different things and need different gates.

**mitigation_type**: process-gap
**structural_mitigation_candidate**: true

---

## Tutorials must be verified against BOTH the executed code AND the live CLI parser
<a name="vignette-cli-flags-and-runnability"></a>

**Tags**: docs, vignettes, tutorials, cli, verification, staleness, notebooks

Building the ViralScan vignette suite (2026-07-20) surfaced three staleness/verification traps that a
notebook-execution harness alone does not catch:

1. **`[skip-ci]` notebooks silently rot.** The old `cell_type_enrichment.ipynb` passed a bare `dict`
   to `cell_type_enrichment()`, but the function had since moved to attribute access (`cfg.cell_types`
   via `RunConfig`). It never failed because the notebook was `[skip-ci]` — never executed. Any
   tutorial not run in CI is presumed broken until proven otherwise.

2. **An exec harness verifies code cells, not the CLI flags in markdown bash blocks.** Flag names
   written into ```bash examples were *inferred from `RunConfig` field names* and were wrong: the
   `hostresponse` subcommand uses `--label {raw,cpm,fraction}` / `--depth-match`, NOT
   `--hostresponse-label` / `--hostresponse-depth-match` (those are the config-field spellings). The
   `evidence` example also omitted the required `--run-dir`. Wrong flags in a tutorial are the same
   defect class as the dict-vs-RunConfig staleness.

3. **"No exception" ≠ "correct value."** A green exec can still print a pedagogically broken number.
   The headline demos were value-checked (EM recovery 7.17×; B-cell enrichment OR highest, padj 3e-19),
   not just run-to-completion.

**How to apply**: for docs/tutorials, (a) make them CI-runnable wherever possible and actually execute
the code cells; (b) grep every `--flag` out of the docs and verify each against `<tool> <sub> --help`
(RunConfig field names are NOT CLI flag names — argparse renames them); (c) for headline/"utility"
examples, assert the *output value*, not just no-exception; (d) treat any `[skip-ci]` doc as
unverified. Data reproducibility: confirm example data is git-tracked (`git ls-files` / `check-ignore`)
before a notebook "reuses on-disk artifacts" — ignored large files = a local-path dependency.

**mitigation_type**: process-gap
**structural_mitigation_candidate**: true

---

## Reorganize a sprawling repo non-destructively with a gitignored symlink view
<a name="curation-symlink-view"></a>

**Tags**: repo-organization, curation, symlinks, non-breaking, tooling, gitignore

When analysis artifacts are scattered (`analysis/<theme>/`, `results/`, `scripts/`, `docs/`) but are
load-bearing for imports, scripts, tests, packaging, and CI, physically reorganizing them breaks
paths. A safe alternative: a **gitignored, regenerable directory of relative symlinks** generated
from a tracked manifest by a tracked builder. Originals stay canonical; the view is a disposable
lens. `git`/CI/packaging see nothing (`git check-ignore` + absent from `git status`). See
[[curation-symlink-view]] decision (2026-07-20) — `curation/` from `scripts/curation_manifest.yaml`.

**How to apply**: (1) relative symlinks (`os.path.relpath`) so the view survives a repo move; (2)
gitignore the view, track only manifest+builder (the reproducible recipe); (3) idempotent build +
`--clean`; (4) skip-and-warn on missing manifest targets so the view never carries dangling links as
files move; (5) guard destructive rebuild with a marker file so the builder can't clobber a
non-generated directory. Organize two ways at once — by narrative (paper results) and by artifact
type — since curation and scanning want different shapes.

**mitigation_type**: tooling
**structural_mitigation_candidate**: true

---

## Committed profiling artifacts go stale — re-profile the current code before trusting them
<a name="stale-committed-profiles"></a>

**Tags**: performance, profiling, cprofile, verify-by-artifact, multimapping, gotcha

A performance review of ViralScan (2026-07-21) nearly reported the wrong bottleneck: the committed
`analysis/multimap_profiling/outputs/cprofile_em.txt` shows `_matrix_value` / scipy sparse
`__getitem__` consuming 403s of 485s in `build_multimap_layers`. But that hotspot was **already
fixed** in commit 3c53ad7 (a vectorized CSR-buffer path taken whenever `original_counts` is sparse —
i.e. always, in production); the committed profile predates the fix. A fresh cProfile on synthetic
sparse input confirmed `_matrix_value` is no longer called; the real current cost is the pure-Python
per-record loop (`list.append`, `dict.get`, and a per-record `pd.isna`).

**How to apply**: treat a checked-in profile like a checked-in benchmark number — it reflects the code
*at capture time*, not now. Before acting on it, (1) check its date against `git log -S<hotspot>` for
later optimizations, and (2) re-run a quick profile on the current code with a representative synthetic
input. For sparse-matrix code specifically, confirm which branch runs (element-wise `matrix[i,j]`
`__getitem__` vs. direct CSR `indptr/indices/data` buffer access) — they differ by ~orders of
magnitude and a stale profile can point at a path that no longer executes.

**mitigation_type**: process-gap
**structural_mitigation_candidate**: true

---

## Collapsing linear-in-weight records is an exact, big speedup — but measure the duplication factor on REAL data first
<a name="collapse-linear-records"></a>

**Tags**: performance, multimapping, vectorization, scipy-sparse, golden-test, verify-by-artifact

`build_multimap_layers` looped over ~100M BUS records emitting per-gene shares. Every share is
`count · k(cell,ec)` (linear in count for a fixed cell+EC), and scipy sums duplicate COO `(row,col)`
entries order-independently — so pre-summing counts by `(cell_idx, ec_id)` is mathematically exact
and cut the loop ~3x (real BUS data has ~3.18x `(cell,ec)` duplication, one record per UMI). Verified
with a golden harness (all 8 output layers × 4 methods at rtol=1e-9) before trusting it.

**How to apply**: (1) the collapse key must be the RAW id (`ec_id`), not a derived key (distinct
gene-set) — different ECs can share a gene set. (2) The win is entirely a function of the real
duplication factor: a naive synthetic profiler drawing records uniformly from a huge (cell×ec) space
has ~0 duplication and will show collapse as a *loss* (groupby overhead, no reduction) — you MUST
measure `n_distinct/n_records` on an actual data file before choosing collapse vs numba vs "it's at
the floor." (3) Don't promise byte-identical for fractional layers: `k·c1+k·c2` (two COO triplets)
vs `k·(c1+c2)` (one) differ ~1 ULP; aim for rtol=1e-9 and keep integer layers exact. Cross-ref
[[stale-committed-profiles]].

**mitigation_type**: technique
**structural_mitigation_candidate**: true

---

## A "monitor to terminal state" instruction is a claim about the world — verify it before acting on it
<a name="verify-tracker-claims-before-acting"></a>

**Tags**: provenance, slurm, tracker-hygiene, verify-by-artifact, silent-correctness, viralscan

`PLAN.md` and `analysis/legacy_v2_v3/TRACKER.md` both told the next session to monitor
fresh-control arrays `25331035`/`25331037` to terminal state. Both arrays had already
terminated with all ten rows failed. The authoritative tracker was asserting something
false about the present, and a session that trusted it would have sat waiting on finished
jobs. One `sacct` call settled it in seconds.

**How to apply**: (1) when a tracker says "in flight", "running", or "monitor X", query the
scheduler/service before planning around it — the tracker records what was true when it was
written, not what is true now. (2) A failure record is not self-describing: four of these ten
rows had `workflow_exit_code: 0` with complete output trees and were still labeled `failed`,
because a wrapper's artifact check looked one directory too shallow and manufactured exit 65.
Always read the *inner* exit code, not just the aggregate status. (3) Add a `stage` field to
any status record that can fail at more than one point — `workflow` vs `artifact_validation`
is exactly the distinction that would have made this self-diagnosing. Cross-ref
[[silent-correctness]].

**mitigation_type**: technique
**structural_mitigation_candidate**: true

---

## Verify that a pinned DOI actually resolves — an unregistered one fails closed at the worst moment
<a name="verify-pinned-doi-resolves"></a>

**Tags**: provenance, packaging, zenodo, external-apis, release-gate, viralscan

`src/viralscan/data_fetch.py` pins `VIRAL_DATA_DOI = "10.5281/zenodo.20112332"` for the
195-GTF viral panel, which was moved out of the wheel. That record is not registered:
Zenodo's API returns `"The persistent identifier is not registered."` and `doi.org` returns
404, while an unrelated third-party Zenodo DOI in the same repo resolves 200 from the same
host. So `viralscan data fetch` fails for every user from a clean install — discovered only
because five cluster jobs died on it.

**How to apply**: (1) when moving bundled data to an external archive, add a cheap
reachability test to CI that resolves the DOI/URL — a placeholder identifier committed ahead
of publication is indistinguishable from a working one until someone runs it. (2) Prove the
egress theory before blaming the network: resolving a *different* DOI from the same host
separates "record missing" from "no outbound HTTPS". (3) Resist repairing this by writing the
unregistered DOI into a locally built cache manifest just because `cache_valid()` requires
DOI equality — that fabricates provenance for content that never came from the archive.
Cross-ref [[verify-tracker-claims-before-acting]].

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## A schema shipped in two places will drift — pin the copies with a test, not a convention
<a name="pin-duplicated-schema-copies"></a>

**Tags**: schemas, packaging, duplication, silent-correctness, viralscan

`validation_protocol.schema.json` exists twice: `schemas/v3/` (what
`scripts/validate_v3_protocol.py` actually reads, via `REPO_ROOT / "schemas"`)
and `src/viralscan/schemas/v3/` (what ships in the wheel). They were byte-
identical, and nothing enforced that. Editing the packaged copy produced a
validator that still rejected the document, because the validator was reading the
other file — the edit was invisible to the thing it was meant to change.

**How to apply**: (1) when a config or schema is duplicated for packaging, add a
test asserting the two files are byte-identical; a convention in someone's head
does not survive the next session. (2) Check which copy the consumer actually
resolves before editing — `grep` the default path constant, don't assume the
`src/` copy is canonical. (3) Avoid `json.load`/`json.dumps` round-trips to edit a
hand-formatted JSON file; the reformat buries the real change in an 18 KB diff.
Do a targeted text insertion, then parse once to fail closed on malformed output.
Cross-ref [[verify-tracker-claims-before-acting]].

**mitigation_type**: structural
**structural_mitigation_candidate**: true

---

## An integrity rail is only as good as the list of things it covers — enumerate the exclusions
<a name="enumerate-what-a-rail-excludes"></a>

**Tags**: provenance, integrity, review, silent-correctness, preregistration, viralscan

A digest-and-ledger rail was added so no frozen protocol section could be edited
without a recorded amendment. It covered five sections named in a dict. An
independent reviewer then changed an already-frozen random seed and an
already-frozen factor level and got **zero validation errors**, because seeds and
factors were not sections and nobody had asked what the map left out. The
protocol's promise that seeds are never chosen after seeing outcomes had no
evidence behind it, in the same commit that claimed to make amendments
tamper-evident.

**How to apply**: (1) when you build a rail keyed on an explicit list, write the
complement down — what is *not* on the list, and whether any of it is load-bearing.
The five covered sections were the ones with a convenient `status` and
`contract_sha256` field; coverage followed the existing data shape rather than the
actual risk. (2) Ship the exclusions as a field in the artifact itself, not as a
commit-message aside, so the next reader sees the boundary. (3) Have the
adversary run the experiment: this was found by someone *mutating the file and
rerunning the validator*, not by reading the code. A reviewer who only reads
resolution notes confirms intent; one who tampers confirms behaviour. Cross-ref
[[verify-tracker-claims-before-acting]].

**mitigation_type**: structural
**structural_mitigation_candidate**: true

---

## Securing a layer moves the attack one layer down — keep reviewing until a round passes clean
<a name="secure-one-layer-attack-moves-down"></a>

**Tags**: integrity, provenance, review, adversarial-testing, preregistration, viralscan

Three independent review rounds on the same preregistration, three verdicts of
does-not-pass, and **two of the blockers were in fixes written earlier in the same
session**:

- Round 1: frozen sections could be edited and re-hashed with no record. Fixed by
  a digest-and-ledger rail over five sections.
- Round 2: seeds and factors were outside that rail; a frozen split seed could be
  rewritten with zero errors. Fixed by hashing them into a `frozen_inputs` digest.
- Round 3: the ledger the rail depends on was not append-only in fact. Editing a
  record in place, rather than appending, left no trace. Fixed by chaining each
  record's hash to its predecessor's.

Each fix was correct and each moved the weakness down one level: section → the
inputs the section trusts → the record of changes to both. This is the normal
shape of integrity work, not a sign of sloppiness, but it means **one round of
review is never enough** and "the reviewer confirmed my fix" is not the same as
"the property now holds."

**How to apply**: (1) after building any integrity rail, ask what the rail itself
trusts, and whether *that* is checkable — a digest trusts the record of digests; a
record trusts its own immutability. (2) Budget for at least three adversarial
rounds and stop only when one passes clean, not when the previous round's findings
are closed. (3) Insist reviewers **mutate the artifact and rerun the check**;
every blocker here was found that way, and none by reading resolution notes. (4)
Make each fix state its own trust boundary in the artifact, so the next reviewer
starts where the last one stopped instead of rediscovering it. Cross-ref
[[enumerate-what-a-rail-excludes]].

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## Adversarial attention goes where you point it — scope a review at the thing you have not been reviewing
<a name="scope-review-at-the-neglected-layer"></a>

**Tags**: review, experimental-design, preregistration, adversarial-testing, viralscan

Five independent review rounds on the same preregistration. Rounds 1 to 4 all
landed on the integrity machinery — digests, ledgers, hash chains — because the
first round happened to find a defect there and each subsequent round inherited
the frame. Each found something real, and the layer genuinely hardened. But four
rounds of rail-chasing meant the **experimental design** had received almost no
adversarial attention.

Round 5 was scoped explicitly to the science and told not to re-litigate the
rails. It returned three blockers immediately, including two that would have
wasted the artifact about to be built: evidence tiers defined only through a
tool-specific layer with no comparator adapter, so 72 of 90 planned rows could
never produce the metric they existed to produce; and a five-factor stratified
split that computes to ~324 strata against ten available biological samples.

**How to apply**: (1) after two or three rounds converge on one layer, that is
evidence the *other* layers are unexamined, not that they are sound. Re-scope
rather than continue. (2) Say explicitly in the brief what is out of scope — "the
integrity layer has a defensible resting point, do not re-litigate it" — or the
reviewer follows the same gradient as the last one. (3) Review a specification
*before* implementing it: the stratification defect was a text change at review
time and a full rebuild once the generator had run. (4) Ask reviewers to
**compute** rather than assess — "count the actual samples available" produced a
number that no amount of reading the prose would have surfaced. Cross-ref
[[secure-one-layer-attack-moves-down]].

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## Fixing by addition leaves the contradiction in place — amend the old text, do not append a correction beside it
<a name="amend-do-not-append-corrections"></a>

**Tags**: review, documentation, specification, silent-correctness, viralscan

Round 5 of a protocol review found that five of six declared comparison axes were
uncomputable. The fix added a new field, `predecessor_incomparable_axes_rule`,
saying so. It did not touch `predecessor_comparison_rule` four lines above, which
still listed all six. The document now asserted both readings at once, and the
commit message said the finding was resolved.

Round 6 caught it by **diffing the commit** rather than reading the message, and
found the identical pattern a second time in the same batch: an estimator amended
to disclaim a stratification while the denominator beside it still promised one.

The pull toward appending is strong because appending is safe — it cannot break
anything that parses. That is exactly why it is dangerous in a specification: the
old text keeps its authority, and a later reader has no way to know which clause
is current.

**How to apply**: (1) when resolving a finding in a document, grep the whole
artifact for every statement of the thing you are correcting, not just the field
the finding named. (2) Treat "I added a field that says the opposite" as an
unfinished fix. (3) Have someone diff the change against the claim, since a commit
message asserting resolution is not evidence of it — both misses here were
invisible from the message and obvious from the diff. (4) Where a number is
asserted, recompute it: the same batch pinned a per-stratum floor of two that the
document's own apportionment rule makes degenerate, which one line of arithmetic
would have caught. Cross-ref [[scope-review-at-the-neglected-layer]].

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## One claim spread across three fields gets corrected three times — grep the claim, not the field
<a name="grep-the-claim-not-the-field"></a>

**Tags**: review, specification, documentation, silent-correctness, viralscan

A protocol said its limit-of-detection fit was stratified by chemistry and
homology. That single claim lived in three places: the estimator, a denominator's
`population`, and an endpoint's `unit`. Round 5 corrected the estimator. Round 6
found the denominator still contradicting it and called it a blocker. Round 8
found the endpoint unit still contradicting both and called it a blocker again.
Three rounds, three blockers, one claim.

The same shape produced a separate finding: a truth population that three
denominators referenced had no dataset, no asset, and no columns, because each
earlier inventory had enumerated *artifacts that existed* rather than *populations
the metrics require*.

**How to apply**: (1) when a finding names a field, search the artifact for the
*claim* that field encodes, not the field name — every place the same fact is
asserted must move together. (2) Enumerate from the consumers backwards: list what
every metric and denominator requires, then check each requirement has a declared
source, rather than listing what exists and assuming it is sufficient. (3) When you
fix one instance, write the regression test against the *invariant* ("no field
promises a stratification the estimator disclaims"), not against the one string you
changed. (4) Expect the third instance. Two corrections of the same claim is
evidence of a third, not evidence of completeness. Cross-ref
[[amend-do-not-append-corrections]].

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## Fix the reason a defect class recurs, not the instances — and expect the fix itself to need one round
<a name="fix-the-class-not-the-instance"></a>

**Tags**: review, integrity, specification, adversarial-testing, viralscan

Ten adversarial rounds on one preregistration. Rounds 6, 8, and 9 each found a
stale contradiction in the *same section*, corrected it, and moved on. Round 9
finally asked why that section specifically, and found the answer: `endpoints` and
`hypotheses` carried no digest and appeared in no coverage list, so nothing could
ever detect drift there. Three rounds of instance-fixing, one round of
class-fixing.

Two other observations from the same sequence:

**Expect the fix to need a follow-up round.** Rounds 9 and 10 each found
second-order defects in the previous round's fixes — a digest check that matched a
value anywhere in the file rather than in the section its scope named, a digest
change omitted from its own record's bookkeeping. These were caught one round
later rather than three, which is the improvement; they were not avoided.

**A reviewer can be right about the defect and wrong about the direction.** One
blocker was reported as "delete this word from the unit"; the denominators showed
the *estimator* was the wrong side. Fixing as instructed would have propagated a
wrong requirement into the artifact about to be built. Check the reasoning behind
a proposed fix against the rest of the document before applying it — and record
the disagreement, because the next reviewer will re-derive it.

**How to apply**: when the same section yields findings twice, stop fixing
instances and ask what makes that section undetectable. Budget one confirmation
round after any fix to the checking machinery itself. Treat a reviewer's proposed
remedy as a hypothesis, not an instruction. Cross-ref
[[grep-the-claim-not-the-field]].

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## A frozen constant no call site passes is decorative — the protocol side of a contract cannot show you the code side

**Date**: 2026-07-27
**Context**: mycelium six-agent review of `codex/viralscan-v3` vs `main`, after
eleven rounds of SCI-05 independent protocol review had converged (rounds 10 and
11 returning no blocker).

The clearest Major the review found: `protocol.yaml` declares
`harmonization.cell_universe.host_only_emptydrops.seed_source: seeds.cell_calling`
with the frozen value `20260727002`. The plumbing is complete on both ends —
`emptydrops.R:65` calls `set.seed(seed)` from `args[[6]]`, and `emptydrops_cells`
appends `str(seed)` as the sixth argv element. It is severed in the middle:
`call_cells:208` passes `rscript`/`fdr`/`lower`/`niters` and stops, so emptyDrops
always runs at the signature default `seed=100`.

Eleven adversarial rounds could not have caught this. Those rounds read
`protocol.yaml`, and a protocol reviewer reads the promise; only a code reviewer
reads the delivery. It is invisible from the artifact too — the run succeeds and
produces plausible cells, exactly as a working implementation would.

**The counter-lesson, from the same review.** I first reported a second finding of
apparently the same shape — `verify_frozen_fastq` emitting
`audited_storage_sha256` after checking only `st_size`, which I wrote up as a field
naming a guarantee the code never delivers. Wrong. `audit_fastq_pair.py:63` streams
a real SHA-256 at freeze time; the docstring says "verify runtime size against an
identity established by a full-stream audit"; the field `runtime_size_verified`
says exactly what it does. It is a documented cost trade-off with a narrower
residual gap (no re-verification between audit and run), not a lie. Pattern-matching
the second finding to the first is what produced the overstatement — I had a
compelling shape and stopped checking.

**Why**: contracts split across a declaration and an implementation have a seam,
and review scoped to either side alone never crosses it. The declaration side is
the one that gets reviewed, because it is where the science is written down. But
"the code doesn't do what the name says" is also the most seductive finding shape
available, so it attracts false positives at the same rate.

**How to apply**: for every entry under `seeds:` and every verification claim in
`analysis/v3_validation/protocol.yaml`, grep for the consumer and confirm the value
reaches the call. Before writing up a name-versus-behavior finding, read the
docstring and grep for the function that *establishes* the value — half the time
the guarantee is delivered somewhere else and the name is honest. Cross-ref
[[grep-the-claim-not-the-field]] — that learning was about one claim in three
places within a document; this one is about one claim on both sides of a
document/code boundary.

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## An enforcement rail that binds only the final state lets the draft drift, and the divergence surfaces at the worst moment

**Date**: 2026-07-27
**Context**: fixing the 2026-07-27 code-review findings. A one-sentence
clarification added to `partitions.sibling_pair_rationale` changed that section's
`contract_sha256`. I wrote in `PLAN.md` that no ledger record was needed "because
the section is `pending`" — and only checked after an advisor pushed.

`validate_amendment_ledger` walks the digest chain of every section whose status
is `frozen`, and skips the rest. That is the right rule for a section still being
drafted. It is the wrong rule for a section that *already has ledger records*:
`partitions` carried an eight-link chain reconciled at `0dae092a…`, and my edit
orphaned it. Nothing reported this. The draft gate stayed green, the append-only
check stayed green, and the divergence would first have surfaced at the moment
`SCI-03` flipped the section to frozen — when the chain would suddenly assert a
history that never happened, with no way left to reconstruct which edit broke it.

The fix was a record (`DEV-019`, which also re-digests the frozen `frozen_inputs`
because declaring the record changes it), plus a new
`pending_section_ledger_drift` check that blocks the training phase.

Two failures of mine worth separating:

**I asserted rather than checked.** "The section is pending, so no record is
required" is a claim about the enforcement rule, and I wrote it into the tracker
without running the two-line query that would have settled it. Exactly the shape
of the finding I had just written up an hour earlier, about fields that name a
guarantee nobody verified.

**My first check was wrong and looked convincing.** The ad-hoc script I wrote to
survey the drift read only each record's top-level `protocol_sha256_after`,
missing links recorded in `additional_digest_changes`. It reported *two* orphaned
sections and implied both predated me. Using the codebase's own
`_ledger_chain_for_section` gave the true answer: one section, orphaned by my own
edit. A quick script that reimplements a subtlety the real code already handles
will confidently produce a wrong picture — and a wrong picture that blames
history is more comfortable than one that blames this turn's edit, which is
exactly why it deserves the extra minute.

**Why**: rails scoped to a terminal state defer their cost to the transition into
that state, which is the point of highest commitment and lowest ability to
reconstruct what happened.

**How to apply**: when a check binds on `status == "frozen"` (or `published`, or
`released`), ask what accumulates in the states it skips. If those states carry
records that assert history, bind the *existence of records*, not the status.
Before writing "no record/migration/approval is required" into a tracker, run the
query. When surveying repo state, call the repo's own resolver rather than
reimplementing it. Cross-ref [[frozen-constant-no-call-site]].

**mitigation_type**: structural
**structural_mitigation_candidate**: true

---

## Two rules writing one file with mode "w" is a race decided by DAG order, and it is invisible from either call site

**Date**: 2026-07-28
**Context**: auditing `SW-04` (does `rerun-multimap` regenerate every
method-dependent artifact?) one artifact at a time.

`multimap.py` and `detection.py` both opened `{output}/summary.txt` with mode
`"w"`. detection runs later, so multimap's three totals — including "Total viral
molecules (selected method)", the headline quantity of the tool — were computed
and truncated away on every run. Never published, since the beginning.

Neither call site is wrong on its own. `open(path, "w")` in a script that owns
its summary file is ordinary. The defect exists only in the pair, plus the rule
order, and neither is visible from the file you happen to be reading.

**The reason it survived**: nothing consumed the output. A repo-wide search for
the literal strings found no test, no doc, no notebook, no parser. A missing
output cannot fail a test that was never written, and it looks identical to an
output nobody asked for. Absent consumers are not evidence that an artifact is
unneeded — here they were evidence that nobody had ever seen it.

**The counter-lesson, again, from the same audit.** I opened this audit by
asserting a *different* bug: that the rerun layer swap corrupts
`molecule_audit.allocated_ambiguous_mass`. I built a fixture, ran it, and watched
the contract break — convincing. It was wrong. `host_conservative` divides
`count / sum(cons_eligible)` across eligible genes rather than dropping mass, so
allocated mass is method-invariant, and `MoleculeAudit.validate` *requires* that,
which means my fixture was a state the pipeline cannot produce. I had reproduced
a bug in my fixture, not in the code. Second time in two days that a compelling
shape outran the checking; the tell both times was that I built the demonstration
before reading the code that would have refuted it.

**Why**: pipeline-shaped code hides defects in the *relations between* steps —
write order, sentinel drops, which rule owns which file — while review attention
lands on the steps themselves.

**How to apply**: grep every `open(..., "w")` (and `to_csv`/`write_text`) against
a shared output directory, group by path, and check for more than one writer; if
there is one, the DAG order decides the content. When a "regenerate everything"
requirement appears, enumerate the artifacts and audit each one — the answer per
artifact is usually different, and three of the six here needed no work at all.
And when a reproduction is cheap to build, read the code that would refute it
first. Cross-ref [[frozen-constant-no-call-site]],
[[pending-section-ledger-drift]].

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## A test that has never failed has not been tested — mutate the fix before flipping the checkbox

**Date**: 2026-07-28
**Context**: closing `SW-05`, the integration test asserting that no stale
artifact survives a `rerun-multimap` method change.

Eight assertions, all green on the first run. That is exactly the state in which
a test is worth least: green against code you just wrote, having never
demonstrated it can go red. Two of this week's defects survived *because* the
tests around them could not fail — `test_fresh_control_refuses_fastq_storage_size_drift`
used a fabricated `"a"*64` digest against a path that never hashed anything, and
multimap's three summary totals had no test at all because nothing consumed them.

So before flipping the row, I reintroduced each defect and confirmed the specific
assertion that should catch it does:

- `clear_stale_virus_outputs` forced to return `[]` → fails
  `test_the_demoted_virus_leaves_no_hostresponse_csv_behind`
- manifest path reverted to `sample_dir.parent / RUN_MANIFEST` (the original bug)
  → fails `test_provenance_names_the_new_method`

Both failed on the named assertion and nothing else, which is the useful signal:
the test is specific, not merely sensitive.

**Why**: a green suite measures agreement between code and test, and writing both
in one sitting guarantees agreement regardless of whether either is right. The
mutation is the only cheap evidence that the test constrains the code rather than
describing it.

**How to apply**: for any test written to close a specific defect, re-break the
defect and watch that test fail before marking the work done. Two minutes. If the
defect cannot be re-broken by a small edit, that is itself informative — it
usually means the test is asserting on a different mechanism than the one that
failed. Also state in the test's docstring what it does *not* cover: this one
does not invoke Snakemake, so it proves nothing stale survives when rules re-run,
not that the DAG re-executes. Cross-ref
[[two-rules-one-file-mode-w]], [[frozen-constant-no-call-site]].

**mitigation_type**: process
**structural_mitigation_candidate**: true

---

## I shipped a regression from an unverified agent finding — running the thing is the only check that would have caught it

**Date**: 2026-07-29
**Context**: `SW-10`, the tiny end-to-end run. It immediately falsified a fix I
had committed the day before.

A multi-agent review reported that `rerun-multimap` rewrote `run_manifest.json`
at the tree root while "the manifest lives beside each sample's `config.yaml`",
making the rewrite a silent no-op. The finding was detailed, cited file:line,
quoted real code, named a plausible trigger, and had survived an adversarial
verifier prompted to refute it. I spot-checked five of the review's findings.
This was not one of them. I implemented the fix, wrote four tests encoding the
claimed layout, and committed.

Running the pipeline end to end took nine minutes and settled it in one command:

    $ find out -name run_manifest.json
    out/run_manifest.json          # the root — beside out/<sample>/, not inside

The original code was right. My fix moved the rewrite to a path that never
exists, so the manifest stopped being updated at all — strictly worse than the
behaviour it "fixed". My tests passed because I had built the fixture from the
finding's description rather than from a real run's output.

**What made this fail closed nowhere.** The finding was internally consistent: if
the manifest *were* per-sample, every sentence in it would be true. Nothing in
the codebase contradicts it in a single file — the layout emerges from
`prepare_output_directory` writing to `--output` and `createconfig` creating a
subdirectory beneath it, two facts in two modules. Static reading can support
either conclusion; only the filesystem after a run distinguishes them.

**Why**: a fixture built from a claim tests the claim, not the system. Both my
unit tests and my integration test derived their directory layout from the review
text, so they agreed with it perfectly and proved nothing. This is the same
failure as the fabricated-`"a"*64` digest test, one level up: there the test
could not fail; here it could not fail *for the right reason*.

**How to apply**: before implementing a fix to pipeline layout, ordering, or
filesystem behaviour, produce the state from a real run and look at it —
`find`, `ls`, one command. Never build a test fixture from a bug report's
description of the layout; build it from a run's actual output, or from the code
that creates the directories. Adversarial verification is not a substitute for
execution: a verifier reading the same files reaches the same wrong conclusion,
and its agreement raises confidence without adding evidence. Cross-ref
[[mutation-test-before-flipping-the-checkbox]] — mutation testing proved my tests
constrained my code, which was true and irrelevant, because the code and the
tests were wrong together.

**mitigation_type**: process
**structural_mitigation_candidate**: true

## [2026-09-27] A case-insensitive marker match inflated the EBV latent headline 2.5x

**Tags**: gene-programmes, ebv, catalogue, silent-correctness, verify-by-artifact

**What happened**: the gene-programme generator resolved curated names with
`gene_attr.lower() == needle`, so the EBV latent marker `BARF1` (BamHI-A)
also matched `BaRF1` (BamHI-a, the lytic ribonucleotide reductase). BARF1 itself
is latent only in epithelial cancers, not B cells. On the EBV LCL run
`BARF1.2` carried 13,668 uniquely-placing molecules and 0 allocated ones, and
it became the docs' showcase for the unique-layer design. Removing both
markers moved unique-layer latent calls from 2,240 to 895 and shrank the
unique layer's claimed latent-sensitivity gain from +75 % to +4.6 %.

**Why it matters**: the number that "proved" the design was produced by the
mislabel it should have been checked against. Nothing crashed; every test
passed, including two that could not fail (`PROG-10`).

**How to apply**: gene symbols are case-sensitive identifiers (BARF1 vs
BaRF1, ORF vs orf). Curated biology tables need a primary citation per row,
and a result that dramatically favours a design should be re-derived with its
single largest contributor removed before it is quoted.

**mitigation_type**: code
**structural_mitigation_candidate**: false

## 2026-09-27 — Unmapped 5′ scRNA reads are mostly TSO and homopolymer junk; filter before any viral alignment

In the SFL tonsil 5′ v3 library, 79 % of the 118.9 M unmapped GEX reads fail a
simple prefilter:
- 62.6 M are shorter than 50 bp once the TSO or adapter is trimmed;
- 30.3 M carry a homopolymer of 20 or more;
- 0.95 M have one base above 60 %.

Unfiltered, 1 M of these reads produced 486,044 viral "hits" and 72 M
alignments with `minimap2 -N 200`, and a collect-all parser ran out of memory
at 24 GB. The hits were:
- TSO + poly(T) reads on the HCV 3′ poly(U) tract;
- poly(C) reads on the EMCV poly(C) tract;
- poly(A) + TSO-rc reads on A-rich anellovirus regions.

After filtering, the whole library ran in 8.5 minutes.

Rules:
- Trim the TSO (`AAGCAGTGGTATCAACGCAGAGTAC`), its reverse complement, and the
  Illumina adapter.
- Filter on complexity.
- Parse aligner output as a stream, grouping records by read.

Two more checks:
- A 64-bp windowed DUST does **not** catch CAG trinucleotide repeats (score
  ≈ 9). A minimum aligned length of 50 bp is what removed them.
- RefSeq `NC_001526.4` (HPV16) is linearised at E1, not at the K02718 origin,
  so literature coordinates such as p97 must be shifted (p97 → 7139).

### [2026-09-29] Snakemake: an unconsumed extra output of a rule never triggers its rerun; a rewritten input does
- Category: pipeline gotcha
- What happened: adding results/virus_identity.tsv to rule analysis did not reschedule analysis in old run dirs (no rule requests it), while rewriting config.yaml (input of kb_count) rescheduled kb_count and everything downstream.
- Why it matters: "skip step X" rerun commands silently rerun X if they touch any of its inputs; new outputs are not backfilled until consumed.
- Resolution: probe with `snakemake -n` on a `cp -r --attributes-only --preserve=timestamps` skeleton copy (empty files, real mtimes). Logged as PLAN SW-22.
- Tags: snakemake, rerun, mtime, dry-run

### [2026-09-30] Agent worktrees can start from a stale base commit
- Category: git-workflows
- What happened: two of three `isolation: worktree` subagents found their worktree on 2ef70f2 (an old main merge), not the orchestrator's HEAD 815838a, and had to `git reset --hard 815838a`.
- Why it matters: an implementer on the wrong base silently diffs against old code; merges then drag in or revert unrelated history.
- Resolution: state the exact base sha in every implementer brief and have the orchestrator check `git merge-base HEAD <sha>` per worktree before merging.
- Tags: git, worktree, subagents, orchestration

### [2026-09-30] Worktree test runs miss untracked data; compare skip counts in the main tree
- Category: testing-patterns
- What happened: in worktrees, tests needing untracked bundled GTFs failed or skipped (test_index_reconciliation, test_ncbi_fetch), while the main tree ran 0 skipped. Golden tests with absolute paths still ran.
- Why it matters: a worktree "green" is not the gate; data-dependent regressions only show in the tree that has the data.
- Resolution: record the baseline pass+skip count in the main tree before launching, and gate every merge on that count there.
- Tags: pytest, worktree, untracked-data, baseline

### [2026-09-30] Structural guards built on "every index measured" miss third-party index builds
- Category: bioinformatics
- What happened: the MECH-A guard (a GTF gene whose t2g column 5 is an index transcript is host) held on 465,769 host rows / 0 viral rows across our combined indexes, but evonk's covid index (kallisto 0.51.1 build) has a viral gene (HUM_HERP6B_DR1) whose column 5 is its own transcript ID, so it flipped to host (1 UMI).
- Why it matters: a rule validated only on indexes we built encodes our builder's row shape, not a property of kb t2g files.
- Resolution: count-parity diff against the previous code on stored third-party runs caught it; logged as PLAN MECH-A step 4a.
- Tags: kallisto, t2g, guard, parity-check, third-party-reference

### [2026-09-30] Snakemake ignores a missing intermediate input until a downstream job must run
- Category: data-pipelines
- What happened: after the consumer rules gained `results/virus_identity.tsv` as an input, a complete pre-35940ec run dir dry-ran as "Nothing to be done" even though the table was missing. With `detection.done` and `umap.done` removed (as `rerun-multimap` does), snakemake scheduled analysis → multimap → detection → umap. After backfilling the table with `analysis.txt`'s mtime, only detection and umap were scheduled.
- Why it matters: a new rule input in an old run dir only hurts once something downstream reruns. Then it cascades from the producer, and the rerun is slow and can change the counts.
- Resolution: backfill the missing intermediate in-process, set its mtime equal to a sibling output of the same rule, and check with a `cp -r --attributes-only --preserve=timestamps` skeleton dry run.
- Tags: snakemake, rerun, mtime, intermediate-files, backfill

### [2026-09-30] Self-named t2g rows (tx = gene = column 5) are viral records, not host cDNA
- Category: bioinformatics
- What happened: the 4a false positive was wider than logged — all 97 HUM_HERP6B_* genes of the covid index are self-named (VIRTUS-style: one FASTA record per gene), not just DR1. Host cDNA rows that point at an index transcript never have transcript ID = gene ID (0 of 465,769 host rows in each of 5 stored combined t2g files).
- Why it matters: exempting on "column 3 set" would have been wrong (3,037 viral rows in the final panel have it empty); tx == gene == col5 is the discriminator. Un-guarding alone would also have split one genome into 97 accession-keyed viruses.
- Resolution: self-named rows record no accession and are named via the legacy prefix map (f9cc725). Covid parity: every multimap layer identical to 815838a.
- Tags: kallisto, t2g, guard, virus-identity, third-party-reference

### [2026-10-01] A panel rebuild can rename gene IDs; the old run GTF then no longer matches the index
- Category: bioinformatics
- What happened: the CAT-42 rebuild wrote 124 RefSeq anellovirus gene IDs with an accession prefix (`TTV3_gp1` → `NC_014081.1_TTV3_gp1`), so `hpv16_gse189670/ref/viral_panel.gtf` (cut from the 09-28 build) no longer matched. The builder commit responsible is not yet traced.
- Why it matters: `-gtf` decides which genes are viral; a stale GTF silently drops renamed genes. Identity itself is unaffected (t2g column 5 is still the accession).
- Resolution: cut the viral GTF from each build's own `combined.gtf` (seqnames of `viral.fa`) → `viral_ref_cat42/viral_panel.gtf`.
- Tags: gtf, kb-ref, gene-id, reference-rebuild, anellovirus

### [2026-10-01] Anellovirus gene models stop at the ORF ends, short of the polyA site
- Category: bioinformatics
- What happened: on TTV-1 (NC_002076.2) the AATAAA is at nt 3,073 but the ORF1/ORF2 models end at 2,901/2,875 (only spliced gp1 reaches 3,077). Across 1,992 ORF-model anellovirus genomes, the gap from the last exon end to the polyA site is median ~20 nt, 90th percentile ~220 nt (AWTAAA heuristic). 48 genomes are whole-genome pseudo-transcripts.
- Why it matters: 3′ 10x reads sit just upstream of the polyA tail; reads in the uncovered 3′ UTR have no transcript to pseudoalign to, a sensitivity loss for true anellovirus.
- Tags: anellovirus, 3prime, utr, polyA, sensitivity, gtf

### [2026-10-03] A low-complexity N-mask can *raise* kallisto's k-mer count
- Category: bioinformatics
- What happened: cat42b, the cat42 index plus an N-mask of 1,586 bases, has 84,860,192 k-mers versus 84,860,149 for cat42.
- Why it matters: kallisto replaces N with pseudorandom bases rather than breaking k-mers, so masked windows become unique junk k-mers that no read matches. That is harmless, but "masked → fewer k-mers" is not a valid sanity check.
- Tags: kallisto, masking, index

### [2026-10-03] Worktree agents can start from the default branch, not the current HEAD
- Category: tooling
- What happened: one of two `isolation: worktree` agents launched together branched from `main` (2ef70f2), not `codex/viralscan-v3`. It reported "no REF-13 row, no config/". The other agent got the right base.
- Resolution: check `git merge-base <agent commit> HEAD` before merging. The REF-13 change was two new files, so a cherry-pick onto v3 was clean.
- Tags: git, worktree, subagents

### [2026-10-03] The mycelium Stop hook sha256s every untracked file in the repo
- Category: tooling
- What happened: Stop got very slow. `session_file_changes.py` runs `git status --untracked-files=all` and `_stat_fingerprint` → `_content_fingerprint` SHA-256s each listed file in full. The 5′ strand reruns wrote ~246 GB of untracked kb output into the repo, so each Stop read 380 GB (28,170 files) over the network FS.
- Resolution: listed the run-output dirs in `.git/info/exclude` (local, uncommitted). The hook now sees 61 files, 0 GB, and the helper runs in 0.63 s. Upstream fix to propose: hash only when size/mtime changed, or skip files above a size cap.
- How to apply: write large run outputs outside the repo, or exclude them before submitting the jobs.
- Tags: mycelium, hooks, git-status, performance, untracked

### [2026-10-03] Test env and the three-field inventory re-pin
- Category: tooling
- What happened: `python` is not on PATH and `.venv` has no pytest, so CLAUDE.md's `PYTHONPATH=src python -m pytest` fails. A re-pin that changed only the commit and sha256 still failed `test_artifact_inventory`.
- Resolution: run tests with `/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin/python`. A re-pin of `docs/output_reference.md` in `analysis/v3_artifact_inventory.tsv` changes three fields: the git commit, the sha256 of `git show <commit>:<path>`, and the `git show <commit>:<path>` retrieval command.
- Tags: testing, governance, conda

### [2026-10-03] The Stop hook charges a read-only session for a concurrent session's edits
- Category: tooling
- What happened: a read-only grill session (no edits) was blocked at Stop for "12 files changed". Those changes came from another session working in the same tree at the same time: commits ee9c278 and 8a4db0c (chemistry, DEF-02) plus uncommitted PLAN.md/Snakefile/chemistry.py edits.
- Why it matters: the hook diffs against a session-start baseline, not against the session's own tool calls. Concurrent sessions in one working tree block each other, and they invite `.living/` entries about work the session never did.
- Resolution: check `git log`/`git status` against your own tool calls before triaging. Record only what this session did. Run concurrent work in separate worktrees.
- Tags: mycelium, hooks, concurrency, worktree

### [2026-10-03] A published "contigs" FASTA can be ORF1 only
- Category: bioinformatics
- What happened: `spyros-lytras/anellovirus-diversity` ships `data/Modha_contigs.fas`, described as the 829 assembled genomes. All 829 of its records match the metadata's `ORF1_len` exactly, and none matches `wg_len`. It holds ORF1 nucleotide sequence, not genomes. The whole genomes are in `Modha_genomes_annotated.gbk` (829 `ORIGIN` blocks, every length equal to `wg_len`).
- Why it matters: planting reads from it would have sampled only the hypervariable ORF1, left the 5'/3' windows undefined, and biased every divergence estimate downward — ORF1 identity to the nearest panel genome is systematically lower than whole-genome identity.
- Resolution: take sequence and CDS coordinates from the GenBank file, which also keeps both in one coordinate system. Check a sequence file's lengths against the paper's own metadata before using it.
- Tags: anellovirus, reference-data, provenance, ORF1, gotcha

### [2026-10-03] kb uses its own bundled kallisto, and a foreign index can spin for ever
- Category: tooling
- What happened: a tiny end-to-end smoke run built `index.idx` with the `viralscan_bench` env's `kallisto`, then `kb count` ran `kallisto bus` from `site-packages/kb_python/bins/linux/kallisto/kallisto`. It burned 787 % CPU for 13 minutes on **one** read pair before being killed.
- Why it matters: this is the version lock the anellovirus pangenome repository warns about, and it does not fail loudly — it looks like a slow job. Any hand-built index fed to `kb count` is exposed.
- Resolution: build indexes through `kb ref` (as `build_bundled_panel_ref.py` does) so one kallisto writes and reads them, or point `kb` at the matching binary. The cat42b `panel.idx` is unaffected.
- Tags: kb-python, kallisto, index, version-lock, gotcha

### [2026-10-03] Snakemake's preamble makes a `__future__` import a SyntaxError in the rule
- Category: pipeline gotcha
- What happened: the ANDET-09 acceptance run died in the `anello_align` rule with `SyntaxError: from __future__ imports must occur at the beginning of the file`, after the host filter and kb count had already finished (job 25696097_2, ~6 min of cluster time per arm). Snakemake copies a `script:` module with its preamble (`import sys; ...; snakemake = pickle.loads(...)`) inserted at the **top**, above the docstring, so a `from __future__` import is no longer the first statement.
- Why it matters: it cannot be caught by importing the module, by `ast.parse`, or by any unit test, because the file is valid Python on its own. It only fails inside the rule, at the end of the expensive part. Three of the nine `script:` modules had one: the new `anello_align.py` plus `gene_programs.py` and `hostresponse.py`, so `--gene-programs` and the host-response rule carried the same latent crash.
- Resolution: delete the import from any `script:` module (PEP 585 generics evaluate natively on the supported versions); convert `X | Y` annotations to `Optional[...]` where the declared minimum is 3.9. `tests/test_snakemake_script_modules.py` reads the `script:` targets out of the Snakefile and fails on any future import.
- Tags: snakemake, gotcha, scripts, ci, python

### [2026-10-04] Tests that read git-ignored data pass locally and fail on every clean checkout
- What: `src/viralscan/data/*.gtf` has been git-ignored since 6bb5c64, but stale local copies (May) stayed on the developer tree. Two tests (`test_index_reconciliation` Retroviridae and `test_ncbi_fetch` EBER) globbed or read them, so they passed here and failed in every clean worktree (and in CI). One PLAN closure (CAT-27) had also cited such a file as evidence.
- Evidence: a full suite on `git worktree add --detach <scratch> HEAD` gave 2 failures; the main tree gave 0. Fixed in 227a2f6 (hermetic fixtures, tracked tables).
- When useful: before declaring a gate green, and before citing a repo file as evidence, run `git ls-files --error-unmatch <path>` or run the suite in a clean worktree.
- Scope: any repo with ignored-but-present data files.

- 2026-10-04 (MECH-B): re-pinning a row in `analysis/v3_artifact_inventory.tsv`
  means changing three fields: awk `$8` (git_sha), `$11` (sha256), and the
  `git show <sha>:<path>` retrieval command. Separately, "unchanged at threshold
  1" holds for integer count layers only, because EM layers are fractional, so a
  grouped sum can add calls there.

### [2026-10-05] Switching to a stale local `main` silently deletes the gitignored bundled GTFs
- Category: gotcha
- What happened: local `main` was 350 commits behind and still tracked the 195 `src/viralscan/data/*.gtf`, which 6bb5c64 untracked and gitignored ("local copies kept"). `git checkout main` overwrote the ignored local copies without a warning, because git overwrites ignored files freely. The fast-forward to 319fd77 then deleted them, since the new main no longer tracks them. Two tests failed: `test_every_marker_resolves_in_the_bundled_panel` and `test_bundled_panel` (933 genes against an expected >2000).
- Resolution: restore each missing file with `git show 6bb5c64^:<path>` (never `git checkout <rev> -- path`, which stages it). Old main's GTFs were byte-identical to 6bb5c64^. No other ignored path was affected.
- Prevention: update a stale branch without checking it out (`git fetch origin main:main`), or diff `git ls-tree` between the old and new tips for ignored paths before switching.
- Tags: git, gitignore, bundled-data, gtf, checkout

### [2026-10-06] em-cell stores a dense all-gene theta per barcode and OOMs at full depth
- Category: gotcha
- What happened: `rerun-multimap em-cell` on SRR12682296 (793,308 barcodes, 46,238 genes) was OOM-killed at 64 GB after ~20 min in `multimap` (job 25702291, MaxRSS 63.6 GB, TotalCPU 21:45 over 23:36 elapsed, so about 1 core used of 8). `em_cell_abundances` (`multimapping.py:326`) returns a dense length-`n_genes` array, and `cell_theta` (`multimapping.py:679`) keeps one per barcode with ambiguous records: 370 KB per barcode, up to about 293 GB. The source host-conservative run peaked at 9 GB.
- Implication: em-cell cannot run on a 2,000-cell VAL-01 sample at unfiltered-barcode scale without large memory. Keep only the compatible genes per cell (sparse), or restrict it to called cells. Grid jobs that include em-cell need `--mem` sized for this, or em-cell dropped from the grid.
- Tags: em-cell, multimap, memory, OOM, Q10d, rerun-multimap

### [2026-10-06] Stop hook flagged `prompt.MD` on a question-only session
- Category: gotcha
- What happened: session 6339e1f7 only answered a clarifying question (user pasted the "00 — Current repository baseline" spec, no request). The Stop hook still reported `prompt.MD` as changed and demanded a `.living/` update. The file is not in `git ls-files -o -m` and this session did not edit it; it was probably a user-side file that existed briefly.
- Implication: no scientific content to triage. This entry exists only to satisfy the hook. If the spec in `prompt.MD` is later adopted, record that as a decision then.
- Tags: mycelium, hooks, stop-check

### [2026-10-08] `viral_ref_final/build/viral.fa` is not in this checkout; synthetic-test pitfall
- Category: learning
- `viral_ref_final/` does not exist under either working directory, so any script that takes `--panel-fasta` has no built panel to read here. The NCBI cache (`~/.cache/viralscan/ncbi`) holds 2,346 per-accession FASTAs that cover the 2,345 `shipped` rows and is the offline fallback (unmasked).
- k-mer scoring of 26 genomes against 2,319 panel genomes: 22 s, 0.9 GB with numpy (`np.lexsort` on kmer then group); the full 973-candidate run should fit the same envelope plus the candidates.
- Test pitfall: `random.Random(seed).choice(...)` inside a generator re-seeds on every call and yields a constant sequence; create the `Random` once. It made a sharing test pass for the wrong reason before it was caught.
- Tags: panel, k-mer, testing, numpy

### [2026-10-09] `build_bundled_panel_ref.py` reconciliation can be run offline and stopped before `kb ref`
- Category: tooling
- What happened: with all NCBI/Ensembl files cached, the builder reaches step 7 (reconciliation) in about 5 minutes with no network calls; it only demands a non-empty `NCBI_EMAIL`. The background run is stopped by step 8 only if `kb` is off PATH, so I let it fail at step 7 under `--strict-reconciliation`, which exits before `kb ref`.
- Why it matters: it is the cheap acceptance test for any catalogue or GTF change. It found 2 pre-existing misses (TTMDV12 `NC_038359`, `AB303562`): the builder's anellovirus fetch takes only the 2,020 clareaulab accessions, never the `viralscan-refseq` rows in `anellovirus_accessions.tsv`.
- Resolution: open; needs a user decision (fetch those rows or record exclusions).
- Tags: builder, reconciliation, anellovirus, CAT-05

### [2026-10-11] Overlap-group transitive closure chains through short exons
- Category: gotcha
- What happened: regenerating gene_programs.tsv from REF-13-normalised GTFs put EBV BMRF1 (productive) in one group with latent EBNA genes (23 genes, g5). A 37-bp EBNA-1.2 exon at 50134-50170 overlaps BPLF1/BOLF1; build_overlap_groups unions transitively (BOLF1 -> BORF1 -> BaRF1.1 -> BMRF1).
- Why it matters: latent reads would count as productive breadth; tests/test_gene_programs.py:173 fails.
- Resolution: shipped TSV left stale; grouping-rule decision (min overlap length, or split latent/productive first) pending with user.
- Tags: prog-07, ref-13, overlap-groups

### [2026-10-11] Exon overlap is not pseudoalignment ambiguity (Biomni review)
- Category: method
- What happened: A-01 proposed min-overlap / pairwise-independence grouping (T=31/91 bp). Biomni max review: 31 bp aggregate overlap need not hold one shared 31-mer; another informative k-mer can resolve a read; strand and missing k-mers matter; greedy independent set is maximal, not maximum.
- Why it matters: any overlap-group rule must be built on k-mer/read compatibility under the real strand setting.
- Resolution: PLAN AUD-01 rewritten; findings.tsv A-01 fix revised. Several title-based citation rejections (A-05) were wrong: judge support against the claim, not the title.
- Tags: prog-07, overlap-groups, audit, biomni

### [2026-10-11] Check interrupted agent worktrees before starting a PLAN row
- Category: process
- What happened: AUD-20 was implemented by hand while `codex/aud-h` already held a finished WIP commit of the same row; the work was duplicated and thrown away. `git rev-list --count HEAD..codex/aud-*` plus `git worktree list` and per-worktree `git status` finds unmerged or dirty agent work in seconds.
- Why it matters: stopped agents leave committed-but-unmerged branches that never touched PLAN.md. Also, 3 tests (artifact_inventory, claim_registry stale `build_reference.py` pin, docs_consistency `infected_subthreshold_cells`) already failed on `codex/audit-fixes-2026-10` before any of this work; the stop-hook "mycelium-data-tracker" processes hold the session lock while agents run, so Stop stays blocked until they finish.
- Resolution: run that check first when resuming.
- Tags: agents, worktrees, process

### [2026-10-11] Stop hook re-blocks while background agents are active
- Category: process
- What happened: `mycelium-stop-check.sh` only accepts a `.living/` update newer than the last recorded work timestamp, and every tool call by a background agent bumps that timestamp, so Stop keeps failing until the agents finish.
- Resolution: finish or wait for agents first, write the `.living/` and `last-session.md` update last.
- Tags: mycelium, hooks, agents

### [2026-10-12] anndata 0.13 iterates layers as None + X; argparse BooleanOptionalAction help differs <3.11
Category: gotcha
What happened: CI py3.12 resolved anndata 0.13.4/pandas 3.0.6 where set(adata.layers) contains None (the X matrix) and sorted() raised, hidden by a broad except as attempt_contract. Separately, BooleanOptionalAction appends "(default: X)" to help before py3.11, so generated docs differed across versions. macOS APFS is case-insensitive, so r1.fastq.gz and R1.fastq.gz collide.
Why it matters: broad except hid the real error; reproduce with a venv pinned to the CI versions from the Successfully installed log line.
Resolution: skip None layer keys; strip the suffix in the generator; rename test inputs.
Tags: ci, anndata, argparse, macos

### [2026-10-18] Perf wave 1 (B24/B23/B04/B06): card candidates need a combined gate, and governance pins follow the source

**Category**: gotcha

**What happened**: Four accepted byte-identical perf cards were combined in worktree `vs_work/perf_wave1`. Combined checks found two things the isolated cards could not: mypy `[index]` on B23's blocked `unique[start:stop]` (`_marker_matrix` returns `Optional`), and the artifact-inventory/claim-registry tests failing on the B24 sha256 of `build_reference.py` (pins must be re-pinned in a later commit because the pinned `git_sha` has to contain the new blob). `tests/test_performance_exactness.py` also needs the untracked `tests/fixtures/performance_original.py`; without it the whole suite aborts at collection. A baseline from a non-git copy differs only in provenance `kind: installed` vs `checkout`, so compare against a git-checkout baseline.

**Why it matters**: component evidence is not combined evidence; load noise (load 8-61) produced 2x swings on unchanged code paths (an untouched B06 path shows ratio 0.44), so only ratios beyond baseline-vs-baseline spread count.

**Resolution**: B23 `None` pass-through, governance re-pin commits, fixture carried over; see `analysis/computational_bottlenecks_2026-10-17/followups/20261017-122240-G00/c00/C00_REPORT.md`.

**Tags**: performance, governance, mypy, timing-noise

### [2026-10-18] A tracker that exists in three checkouts is not a tracker

**Category**: gotcha

**What happened**: the perf plan doc was untracked in the shared `codex/drop-py39` checkout and was then committed on `codex/perf-wave1-b24-b23-b04-b06` (reconciled, 06:05). Three copies now exist: the branch's reconciled copy, a stale root copy still reading "implementation not started" (Oct 17 12:54), and the wave-2 branch's fork of the branch copy. The wave-2 copy will conflict on rebase, and the branch's own header, a "pending execution" bullet and an empty `C00_REPORT.md` §3b were still stale even after the card bodies were updated.

**Why it matters**: the reader who opens the stale copy is misinformed — they may redo G00 or re-pick a finished card — and half-reconciled trackers are worse than obviously stale ones because the checkboxes and the prose disagree.

**Resolution**: reconcile once on the branch that owns the work, record the reconciliation as its own tracked row (`PERF-DOC-01`), let later branches rebase onto it, and treat the root copy as stale until the PR lands.

**Tags**: performance, governance, documentation-drift, worktrees

### [2026-10-18] Pattern pkill on a shared account kills other sessions' runs

**Category**: gotcha

**What happened**: while restarting its own suite, the wave-1 C00 session ran `pkill -f` and `ps | grep pytest | xargs kill` patterns that were not restricted to its own processes; pytest runs belonging to other sessions in roughly 04:50-05:00 may have been terminated, and the session's own shell died twice (exit 144). See `C00_REPORT.md` §5.

**Why it matters**: on a shared HPC account an aborted run looks like a bug or a flake in the other session's work, and the cause is invisible unless it is written down.

**Resolution**: never use pattern-based kills on a shared account; address the exact PID. The incident is recorded here so a session whose run vanished in that window can attribute it instead of re-debugging.

**Tags**: hpc, shared-account, process-safety
### [2026-10-18] A rebase can fuse two independently-correct features into a crash, and an exactness grid only covers the index spaces its fixtures encode

**Category**: gotcha

**What happened**: On the wave-1 rebase onto main, wave-1's budget-gated column view (layers
rebound to a view over `selected`, per-gene access remapped through `gene_positions`) and main's
Q123 per-virus positive-cell count (`unique[:, idx]` with `idx` built from
`adata.var_names.get_loc`, i.e. absolute) landed in the same loop. Each was correct on its own
branch. Combined, the absolute index was applied to the narrow view. Every exactness fixture —
the 90-case method x layout x dtype grid and the frozen-oracle suite — used leading genes
(`g0..g3`), where `selected` coincides with the absolute indices, so a standalone comparison
against main's own implementation reported 90/90 byte-identical and the merge looked clean.

**Why it matters**: Two green-looking signals ("all tests pass", "N/N identical to main") were both
true while the merged code raised `IndexError` on the real layout. `build_reference.py`
concatenates the host GTF first, so viral genes are always a trailing column block (2,000 host +
40 viral columns reproduced it, and the real call site groups every viral gene). The failure has
two shapes: an index outside the view crashes loudly, but a non-identity `selected` whose width
exceeds the absolute index reads the wrong column silently. Perf work that rebases onto a moving
main is structurally exposed, because the optimization changes *which object* a later feature's
index refers to.

**Resolution**: Probe the index space, not only the values: include a fixture whose referenced
columns are a trailing (or interior) block, and compare against main's implementation on it before
trusting a rebase. Fixed by resolving the index through the same `gene_positions` mapping the
sibling access uses and accumulating the mask per integer column.

**Tags**: rebase, byte-identical, fixture-blindspot, index-space, performance, evidence, multimapping

**mitigation_type**: ambient-awareness

### [2026-10-18] CI lint (mypy) passes locally but fails in CI because of unpinned numpy stubs

**Category**: gotcha

**What happened**: PR #73's "Lint (ruff + mypy)" job failed with two errors in
`src/viralscan/scripts/umap.py`:
`Incompatible types in assignment (expression has type "list[Never]", variable has type
"ndarray[Any, Any]")` at the `counts = []` and `counts_perm = []` lines of the
`viral_neighbor_enrichment` binary-label fast path. The same command
(`mypy src/viralscan`) reported "Success: no issues found" in the local
`viralscan_test_full` environment (numpy 2.4.6). The CI environment (numpy 2.5.4) resolves
`np.mean(..., axis=1)` to `ndarray[Any, Any]`, so the same name being reassigned `[]`
later in the branch becomes a type error; numpy 2.4.6's stubs resolve it to something
compatible with `list`. CI installs lint tools unpinned (`pip install ... numpy pandas ...`),
so the drift is silent and only surfaces in CI.

**Why it matters**: A lint fix that is verified only against the local env can ship a red
CI run (PR #73 was blocked on exactly this). The fast path is also performance-sensitive,
so the tempting "make both branches produce the same type" restructure would risk the
byte-identical-output contract and the isolated speedup; the correct fix is a type
annotation only, which has zero runtime effect.

**Resolution**: Annotated the two accumulators as `np.ndarray | list[float]` at their
declaration sites in `viral_neighbor_enrichment` (no runtime change; byte-identical
outputs preserved; 125/125 `test_performance_exactness.py` and 26/26 `test_umap.py` green).
Verified in a CI-faithful environment built with `uv` (numpy 2.5.4, mypy 2.4.0) *and* in
the local environment (numpy 2.4.6); both clean. Future lint work should be validated in a
CI-faithful env (e.g. `uv venv` + the pinned versions from the "Install lint tools"
step) rather than the shared conda env, and pinning the lint tools in CI would remove the
drift at source.

**Tags**: mypy, numpy-stubs, CI-drift, lint, byte-identical, performance, annotation

**mitigation_type**: ambient-awareness

### [2026-10-18] Perf wave 2: an unused eager read is still the only validator, and `splitlines()` is not a newline iterator

**Category**: gotcha

**What happened**: Two wave-2 gates turned on validation semantics rather than arithmetic. (1) B10 — for the disabled UMAP rule the counts H5AD is read but never used except for one `uns["count_schema_version"]` lookup, so replacing `sc.read_h5ad` looked free; the read is in fact the only validator, and its failure surface spans h5py open errors, per-element anndata checks (`obs` group presence, dense X 2-D, layer shape, contents/row-count agreement) *and* anndata's encoding registry (`IORegistryError` on an unknown `encoding-version` on X/obs/uns), plus read-tolerated damage (deleted X/var/obsm groups, zeroed CSR `indptr`, duplicate `obs` index) that a replacement must keep tolerating. Driven through the unmodified `umap.run`, `uns_only` silently accepted 7 damaged inputs, `backed='r'` 3 (and reordered a two-fault report), a hand-rolled h5py validator 5, and `anndata.experimental.read_lazy` needs `xarray` and errored on valid input in the test env. (2) B25 — `read_text().splitlines()` splits on eight separators beyond newline iteration (`\v \f \x1c \x1d \x1e \x85 \u2028 \u2029`), and chunked incremental UTF-8 decoding reports `UnicodeDecodeError` positions *relative to the failing chunk*, so the same corrupt file yields chunk-size-dependent messages (file offset 12 reported as 12 for 64-byte chunks but 0 for 1- and 4-byte chunks).

**Why it matters**: "the object is unused, so skip the read" is a validation change, not an optimization; a streaming text reader changes the failure contract unless it re-fabricates the codec's error with a global offset. What decides these gates is the *class* of failure (registry-level, element-level, contents-level), not the count of cases in a finite frozen matrix — so freeze the two-fault ordering cases too, they pin precedence.

**Resolution**: B10 recorded no compliant fast path (`cards/b10/decision.md`), keeping the 37-run freeze as the disabled rule's regression net. B25 implemented only the byte-concatenation half (chunked copy with trailing CR/LF holdback; peak flat at two chunks instead of 2–3× the largest input, 63.7× lower at a 64 MB input, ~2× faster) and deferred streaming headers/GTF with measurements attached (`cards/b25/{feasibility,decision}.md`).

**Tags**: performance, anndata, validation, io

### [2026-10-18] Perf wave 2 B18: a fused reducer can be a wash on its own, and tracemalloc invents CPU regressions

**Category**: gotcha

**What happened**: B18's fused one-pass molecule reducer (seen-QNAME set + per-molecule maxima, replacing `best_scores`'s `{qname: {side: score}}`) removes 300 B/read of per-read dicts but holds the seen set at the same time as the molecules — while the frozen two-phase shape frees the per-read dicts before rows accumulate. Measured as `list(molecule_rows(...))` the prototype is 0.968x frozen at 1 read/molecule and 0.471x at 8: the fused parse alone merely swaps which structure peaks. Pairing it with a writer that streams rows from the molecules (instead of a materialized row list) is what removes the trade: 0.728x/0.641x/0.571x/0.526x at 1/2/4/8 reads per molecule, 0.468x at 1 M reads, versus 0.970x at 1 read/molecule for the fused parse *with* rows materialized. Separately, the first wall-clock numbers said 1.06-1.13x slower; without tracemalloc the same code is 1.034x (2/5 trials faster) and 0.936x (5/5 faster) — tracemalloc charges per allocation, so it systematically penalizes whichever arm allocates more per line.

**Why it matters**: an optimization that only changes *which* container peaks has not reduced anything, and a peak-based gate must be driven through the production entry point, not through a convenience API that materializes results. Timing comparisons of allocation-heavy variants must run without tracemalloc enabled, or the measurement reports the profiler's overhead as the change's cost.

**Resolution**: shipped the fused parse *with* the streaming writer in `molecule_verdict.py` (parse completes before the gzip is opened, so the parse-failure boundary and byte stream are unchanged; no temp file because gzip's `FNAME` header embeds the output file name). Evidence: `analysis/computational_bottlenecks_2026-10-17/followups/20261018-112322-wave2-g00/cards/b18/`. Mutant "writer materializes rows" is caught only by the 1-read-per-molecule memory guard, which is where the row list is a significant fraction.

**Tags**: performance, memory, profiling, io

### [2026-10-18] Perf wave 2 C00: two checkouts at different commits cannot be byte-compared without excusing the checkout sha

**Category**: gotcha

**What happened**: The wave-2 combined qualification compared a throwaway git checkout at the base commit against the branch tip and found 147 of 148 harness artifacts byte-identical, with `b24/prepare_reference_inputs.json` differing at exactly three JSON paths: `build_receipt.software_identity.git_commit` (`faf5667...` vs `6de5a65...`). Wave-1's C00 never hit this because its baseline was a distribution-installed copy (whose `kind` differed instead). The field is the one thing that *cannot* be equal between two commits, so a blanket "one file differs" would have been misleading in one direction and a blanket excuse in the other.

**Why it matters**: an exactness gate has to distinguish "provenance" from "behaviour" mechanically, not by eyeballing the diff, or the excuse becomes the hole. The comparer therefore (a) walks both JSON documents and lists the differing paths, (b) only excuses paths matching a single declared suffix, and (c) requires the two documents to be equal *after* scrubbing that field — so a clock change, a reordered key set or a shifted array would still fail even if it landed in the same file.

**Resolution**: `cards`-level receipts stay byte-exact; the combined receipt records the excused path and the scrub proof (`c00/exact_compare.json`, `C00_WAVE2_REPORT.md` section 1). Also worth knowing: the full suite at the tip passed 5,028 tests in 457.9 s, so a combined wave of two small cards is comfortably re-qualifiable in one session.

**Tags**: performance, exactness, provenance, testing

### [2026-10-26] st_dev is a per-host mount id; use inode alone to prove a rename
**Category**: testing
**What happened**: R1 tree tests failed on node exe032 ("inode changed") though the tree was the same; st_dev was 136 vs 105 at capture.
**Why it matters**: A [dev, ino] check false-fails when the test runs on another login or compute node.
**Resolution**: compare st_ino only (test_regression.py:36).
**Tags**: testing, hpc, migration

## L-78 (2026-10-26) Relative `Path.cwd().parents[n]` paths break silently when a tree is renamed; path-rewrite regexes never see them
**What**: After `viralscan_work` moved under `data/runs/work`, `slurm_finish_evidence.sh:18` (`Path.cwd().parents[1]/'vs_pinned/...'`) resolved to `data/runs/vs_pinned` and the queued job failed with FileNotFoundError. The old-path scan and rewrite tools match literal path tails, so they reported the file clean.
**Evidence**: job 25882793 log; `sibling_specific_reads.py:29` same pattern; `confirm_hits.py:411` was dangling even before the move.
**Rule**: before and after any tree move, `grep -nE "parents\[[0-9]\]|Path\.cwd\(\)|\.\./\.\."` the live scripts and fix those by hand; add them to the per-group diff list.
Tags: migration, paths, slurm, testing, gotcha

<!-- 2026-10-27 batch 4e (D64/D69): 28 experiment-only entries live only in the reproducibility repo copy (https://github.com/mdmanurung/ViralScan-reproducibility/blob/main/.living/learnings.md); this copy keeps package and shared entries -->

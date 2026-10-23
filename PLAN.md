# ViralScan 3.0 execution tracker

Status: **active**

Branch: per-unit `codex/<unit>` branches off `main`

Last reconciled: 2026-09-27

Release target: `3.0.0rc1`, then `3.0.0`

Scientific target: methods-grade manuscript generated only from frozen v3 artifacts

This is the single operational tracker. The accepted scientific and product
contracts are in
[`docs/plans/2026-07-22-viralscan-v3-correctness-and-publication.md`](docs/plans/2026-07-22-viralscan-v3-correctness-and-publication.md).
The former pre-v3 tracker is retained only in the local governance archive and
excluded from Git and source distributions because it contains institutional
paths and private-analysis notes. It must not be used to establish v3
completion.

## Next action

- **2026-10-23:** spike-in grid done, flank index retired (F-028, SPIKE-VERDICT). Next: STAR-GATE-01 (user), DISC-12 freeze (#103), rc1 release prep.
- [~] `DISC` Kraken2 cell-resolved taxonomic discovery: plan approved 2026-10-21 (all nine §16 decisions).
  `DISC-R1` (grouping perf) done. `DISC-09` (R2: provenance, fingerprint, publish, resume) done. `DISC-R5` (per-UMI rule: deepest-consistent) done. `DISC-R6` (polyA/T trim before Kraken2) done. `DISC-R7b` (Kraken2 defaults 0.05 / 4) done. `DISC-01` (R3: discovery env, doctor profile) done. `DISC-10` (R4: real-Kraken2 integration test) done. `DISC-11` (R7: `discover-collapse` rank views) done. `DISC-R8` (minimal `candidates.tsv`) done. `DISC-R9` (docs: `docs/discovery.md`) done. `DISC-R10` (EBV real-data validation, 3 arms) done. `DISC-R11` (close-out) done. Defaults 0.05 / 4 confirmed by the user after R10 (2026-10-22). Next: push and PR (user); then `DISC-12`.
- [x] `F027-01` finding F-027 (`.living/findings/hhv6b-flank-reads-missed-by-cds-only-index.md`): HHV-6B reads in the positional 3' flank class (8.9% of viral reads, GSM6416022, 5' v2) are outside a CDS-only kallisto index.
  - Note (2026-10-22): Lareau comparison done against the public reference at caleblareau/hhv6-reactivation@dd2eb77: 103 CDS-only entries, same NCBI CDS coordinates, no flank/UTR, no masking.
- **Focus: STARsolo host+virus readiness.** Phases 1-4 are committed on `codex/star-confirm`; next are SPIKE-01 (phase 5, unrun), the
  STAR-GATE-01 decision, and the phase 6 EBV acceptance run.
- [x] REVIEW-A (WP-1): called-cell summary denominator and host-only plot null.
- [x] REVIEW-B (WP-2): bounded coverage zero filling; preserve depth policy.
- [~] REVIEW-C (WP-3): coverage-map code/fixtures pass; real-directory qualification pending.
  - Note: real Track A BAMs (`VIRUS|U31784.1`) vs a bare-accession GTF gave `annotation_mismatch`;
    `resolve_reference_annotation` now also accepts `ACC` for viral-class references only and
    records `annotation_seqid` / `annotation_seqid_normalised` in the coverage-map JSON.
- [x] REVIEW-D (WP-4): opt-in candidate-complete source-labelled review.
- [x] REVIEW-E (WP-5): pool-scoped tier-2 molecule fingerprint.
- [x] REVIEW-F (WP-9): docs, acceptance ledger, parity/benchmark receipts.
- [x] STAR-HOSTVIRUS-01: STARsolo alignment of host-unmapped reads to one host + full-viral-panel index
  (`viralscan.star_hostvirus`: `write_reference` HOST|/VIRUS| FASTA+GTF, `build_star_hostvirus_index`,
  `align_host_unmapped`; sorted BAM with CB/UB for `plot_coverage_maps`). Code + synthetic tests done; wired into the CLI
  by STAR-03 (`viralscan star-index`) and STAR-04 (`viralscan run --star-confirm`). Real-library runs (EBV GSM4796271,
  HHV-6B GSM6416022) remain Track A diagnostics.
  - Note: multimap cap 100 and relaxed score filters are deliberate (keep host competition visible); STAR 2.7.11b
    still segfaults in Solo with `--soloCBwhitelist None` when zero reads map (ANDET-09-STAR139); handled since
    STAR-02 (shared `star_run.run_star`). Docs (`cli_reference`, `output_reference`, `quickstart`, CHANGELOG) done in phase 4;
    SPIKE-01 (spike-in detection limit) is **not run**, so the arm has no measured sensitivity or specificity, and phase 6
    (real-library acceptance on EBV GSM4796271) is not done.
- [x] STAR-04: `viralscan run --star-confirm --star-hostvirus-index DIR` (branch `codex/star-confirm`): rule `star_confirm`
  (`scripts/star_confirm.py`) aligns the kb input pair (host-filtered reads) with `star_hostvirus.align_host_unmapped(prune=True)`,
  writes `star_hostvirus/viral.bam(.bai)`, `results/star_support_by_virus.tsv` (per `virus_name`) and `log/star_confirm.done`;
  `detection` appends `star_reads_primary, star_reads_any, star_molecules, star_host_competing, star_other_virus` to
  `viral_summary.tsv` only when the flag is on. Labels only: no call, tier or threshold changes.
  - Note: menu dies without the index, without `--host-filter`, on a failing `validate_index_manifest`, on missing STAR/samtools,
    or on an index without the flag. `summarise_by_virus(virus_of=...)` aggregates per virus (a read on several accessions of one
    virus counts once; `star_other_virus` needs a different virus). Review fixes: mypy on `read_index_manifest`, wheel/sdist
    allowlists carry the three STAR modules, menu wiring/prune-TSV/manifest/run_star tests, docs state the HHV-6/7 telomere
    100-locus loss and the Unique-mode `star_molecules` lower bound. Viruses STAR saw nothing for read 0; alignment-only rows stay blank. Strand = `--strand` or `unstranded` (kb's
    effective default). Docs/CLI-reference/output-reference rows landed in phase 4. No real STAR run.
- [x] STAR-ACCEPT-01 (2026-10-21): real-library acceptance of `star_confirm` on EBV GSM4796271
  (`analysis/star_confirm_acceptance/`, Slurm 25860995, 10 min, exit 0): the real Snakefile rule + script
  (`--allowed-rules star_confirm`, so no host_filter/kb_count rerun) on the Track A host-filtered reads, the existing
  index adopted through a symlink view plus a generated manifest. All 2,343 VIRUS| contigs have the same alignment
  record count as the earlier manual STAR run (0 mismatches), so on-list correction and strand change only barcode
  tags, as designed. `star_support_by_virus.tsv`: Epstein-Barr virus 1,170,981 primary / 941,111 molecules;
  Human gammaherpesvirus 4 type 2 502,137 / 724,380; 1,204,497 reads place on both (`star_other_virus`), the known
  EBV type 1/type 2 shared-contig signal; `star_host_competing` is 0 because the input is host-filtered.
  Second library (Slurm 25861544, 24 min, exit 0): R01_SOC4, 3-prime 10xv3, forward, RAW reads (172.4 M pairs), same
  rule: 2,343/2,343 contigs identical to the pilot STAR run (0 mismatches). It exercises the 10xv3 on-list (6.8 M
  barcodes gunzipped), --strand forward, a 12-nt UMI and the zero-viral-reads path (header-only table, no crash), but
  R01_SOC4 has no viral reads, so per-virus counts on 3-prime v3 are untested until the SPIKE-01 EBV spikes. Not yet
  covered: multi-accession viruses other than EBV.
  Exit-139 fallback verified on REAL STAR 2.7.11b (Slurm 25864228, 1 min, exit 0): 2,000 random pairs with the Drop-seq
  chemistry (no on-list, `--soloCBwhitelist None`, zero reads mapped) crash STAR in Solo; `star_run.run_star` retried with a
  12-nt dummy whitelist, logged the fallback, left `star139_fallback.txt` and the pruned BAM, and wrote the header-only
  per-virus table (`analysis/star_confirm_acceptance/fallback139.py`).
- [ ] STAR-GATE-01: let STAR evidence gate or tier calls. Decided only after SPIKE-01 (user decision pending); STAR-04 adds columns only.
  - 2026-10-23: SPIKE-01 is done (F-028): STAR arm threshold <=3 spiked UMIs in both backgrounds, 0-8 non-spike STAR reads in nulls. The decision is now the user's. `fate_host` is not a measure of host competition (pruned BAM); see F-028.
- [x] STAR-03: `viralscan star-index` verb (branch `codex/star-confirm`): `star_hostvirus.star_index` wraps
  `write_reference` + `build_star_hostvirus_index` (reference in `<out-dir>/reference/`, index in `<out-dir>/index/`) and
  writes `star_hostvirus_index.json` (`write_index_manifest` / `read_index_manifest` / `validate_index_manifest`).
  - Note: manifest is written last (absent = unfinished build). Refuses a non-empty index dir; `--force` deletes it only
    if it holds a manifest. No real STAR build run (tests monkeypatch it). Next: STAR-04 (phase 3) validates the manifest.
- [x] STAR-02: harden `star_hostvirus` (branch `codex/star-confirm`): on-list default via
  `resolve_whitelist` (kb on-list gunzipped; literal `None` only for list-less chemistries, explicit list wins),
  `--soloStrand` from a viralscan strand (`SOLO_STRAND`, out of `ALIGN_ARGS`), exit-139 fallback lifted from
  `anello_align` into `viralscan/star_run.py` (both callers; dummy barcode = `cb_len` As), `summarise_by_virus`
  / `write_summary` / `STAR_SUMMARY_COLUMNS`, `prune_to_viral` (`prune=` option, default off).
  - Note: `align_cmd` stays pure (whitelist `None` -> literal `None`; the on-list default lives in
    `align_host_unmapped`) so the existing `align_cmd` test holds. Summary keys are the VIRUS| contig without
    its prefix. A pruned BAM no longer holds the HOST| partners of viral multimappers, so coverage-map
    competition categories need the full BAM or must be taken from the summary TSV.
  - Review fixes: `anello_align._run_star` takes `cb_len` from `cb_umi_geometry(technology)` (was hardcoded 12,
    broke the 139 retry for 16-base chemistries); `star_molecules` documented as Solo-assigned (CB, UB) only
    (not read support; use `star_reads_any`); tests for the on-list default and `cb_len` wiring through
    `align_host_unmapped`.
    still segfaults in Solo with `--soloCBwhitelist None` when zero reads map (ANDET-09-STAR139), not handled here.
- [x] SPIKE-03 (2026-10-22; branch `codex/spikein-posthoc`, from tag `spikegrid-code-26c1ac5`; code only, grid not re-run):
  post-hoc scorer fix, rescore, recovery ratio, aggregate hardening in `analysis/spikein_detection_limit/`.
  - BUG fixed: `virus_identity.tsv` also lists ~41k host ENSG rows (`viral=false`) and contaminant rows, and the extractor and
    `score_condition.kb_stats` treated every row as viral, so the floor columns (`kb_unsupported_molecules`,
    `kb_total_viral_molecules`, null D03-NULL-r1 ~33,399) were host counts. New `score_condition.viral_identity` keeps
    `viral` == true (case-insensitive) and role != contaminant, refuses a table with no `viral` column; used by the scorer
    (`kb_stats`, `evidence_tier`) and `run_condition.extract_kb`. Truth-side columns verified unchanged on D03-L100-s1-c100 and
    D03-NULL-r1 (before/after in `viralscan_work/spikein_grid_2026-10/rescore_check_2026-10-22/`).
  - `rescore.py`: writes `bundle/scores_v2.tsv` (never `scores.tsv`), refuses conditions without DONE; `aggregate.py
    --scores-name` (default `scores.tsv`) so one scorer covers all conditions.
  - `aggregate.py` hardening: threshold labels `<=3` (left-censored) / `>300 (not established)` / `undetermined`; `non_monotone`
    flag + report warning; a level with missing seeds that could still be met is undecided, so missing data is `undetermined`
    and never not-reached or a retire reason; verdict backgrounds from `grid.BACKGROUNDS`; blanks no longer become 0 in
    `recovery_rows` or the null floor; docstring contradiction fixed (old not reached, incl. both, = undetermined).
  - `recovery_ratio.py`: flank/old net kb recovery per background x layout x level (RATIO_MIN 0.80), flags
    `needs_user_decision` only; never changes the Q15 verdict.
  - Q16 wording ("computational detection threshold (source, background and layout-specific)") in README, report, TSV headers
    (`lod.tsv` keeps its file name; column `limit` is now `threshold`). Decision: an existing retire from a complete layout still
    wins over an undetermined one (overall = retire > undetermined > keep, unchanged).
- [x] SPIKE-VERDICT (2026-10-23, grid 25864985 66/66, rescored `scores_v2.tsv`): flank index v2 **retired** by the locked
  rule (Q3/Q15): same threshold as old (D03 10, R01 <=3) but up to 12.57 / 9 unsupported UMIs per null library (limit 5);
  0 false calls in either arm; flank recovers 1.12-1.38x more net kb molecules (no ratio < 0.80, so no ask-back).
  STAR arm threshold <=3 in both backgrounds. Finding F-028; receipts `analysis/spikein_detection_limit/results_2026-10-23/`.
  - PR #109 review: STAR fates came from the pruned VIRUS|-only BAM, so `fate_host` cannot show host competition (host-primary
    spikes read as `absent`); STAR detection is unaffected. `count_pairs` now uses `gzip -dc` (macOS `zcat` wants `.Z`);
    `SHA256SUMS` no longer lists itself.
- [x] SPIKE-01 (2026-10-21, exploratory under the Q21/LOCK-01 A2 override; branch `codex/spikein-detection-limit`,
  stacked on the unpushed `codex/star-hostvirus-align`): empirical detection limit of the old kb index, the flank
  index v2 and STARsolo host+virus, by spiking real HHV-6B reads (GSM6416022 BAM, stratified by flank-audit class)
  into null background libraries. Tooling written and unit-tested (`analysis/spikein_detection_limit/spikein.py`,
  `score.py`, `smoke.sbatch`; `tests/test_spikein_detection_limit.py`, 4). NOT run: the smoke job
  (1 background D03 SRR29320037, levels 10/100, seed 1, 5 cells, first 5 M pairs) is unsubmitted. Next: smoke,
  then the 60-run grid (levels 3/10/30/100/300 x 3 seeds x 2 layouts x D03/R01_SOC4 at ~45 M pairs).
  - Pre-set rules (user "approve all" 2026-10-21): call = >= 3 kb molecules AND >= 1 STAR read on the same contig;
    limit = lowest level recovered in >= 2 of 3 seeds; keep the flank index only if its limit <= the old index's and
    it adds <= 5 unsupported UMIs per null library (3 unspiked replicates per background).
  - SPIKE-02 / task I1 (2026-10-22): spike source generalized to a truth-contig SET so the EBV source (GSM4796271,
    type 1 NC_007605.1 + type 2 NC_009334.1, forward, R01_SOC4 background 3' v3 forward) works. `spikein.py spike
    --contig A,B --max-nh N` (default 1 = old behaviour: NH=1 primary records); N>1 also keeps reads with NH<=N whose
    records inside the set equal NH (reads with a HOST| / other-virus alignment are excluded; counts printed to stderr).
    Class and src_pos come from the alignment on the FIRST truth contig (type 1 for EBV) so labels do not depend on
    which contig STAR called primary. Truth table adds `contigs`, `nh_stratum`, `src_contig`; `contig` is now the
    primary contig. `score.py star --contig` is optional (read from truth table `contigs`, else meta json): n_source_contig_any
    = any alignment on any truth contig; n_other_virus = primary on a VIRUS| contig outside the set. Fixed a latent bug:
    `samtools view -s N.0` (frac 1.0, pool smaller than target) kept nothing; now `-s` is omitted. Tests: 13 (9 old + 4
    new, synthetic SAM->BAM). NOT done: per-stratum recovery reporting in score.py (truth carries `nh_stratum`; kb scoring
    via multimap not touched), 98->90 nt trim effect on NH unmeasured, nothing run on real libraries. For EBV/R01 pass
    `--source-orient forward --bg-orient forward --cb-len 16 --umi-len 12 --r1-len 28 --r2-len 90`.
  - SPIKE-02 / task I2 (2026-10-22) [~]: grid tooling in `analysis/spikein_detection_limit/`: `grid.py` (CONFIG block + 66-row
    `conditions.tsv`: 60 spiked + 3 unspiked replicates x 2 backgrounds; ids have no underscore), `run_condition.py` (one
    condition end to end with the REAL product: seeded cached background subsample, `spikein.py spike`, `viralscan run` for the
    old then the flank index with `--host-filter starsolo`, `--star-confirm` on the old run only, small bundle, scoped deletion,
    DONE last, resumable, `--dry-run`), `grid.sbatch` (array, `--row $SLURM_ARRAY_TASK_ID`), `score_condition.py` (rows per
    condition x arm old/flank/star: kb molecules = unique + ambiguous_allocated over spike barcodes and truth accessions,
    evidence tier, STAR any-alignment on the truth set, the three call booleans, unsupported-molecule floor).
    `tests/test_spikein_grid.py` (9, synthetic bundles + dry runs). NOT run for real: nothing submitted, no FASTQ/STAR/kb touched.
    Decisions: spiked bg replicate = spike seed, so null replicate s is the baseline of spike seed s (`--baseline-bundle`);
    called list frozen per bg replicate = top 100 on-list barcodes (`--cell-calling external`, no emptyDrops); EBV `--max-nh 2`
    (my pick); host filter is repeated for the flank run (product path), STAR columns come from the old run. Unverified:
    spike read names surviving the product host filter into viral.bam (runner fails loudly at level >= 100 with 0 STAR reads),
    `viralscan run` flags exactly as in the explorer report (not executed), multimap h5ad layer names read from code only.
    Not done: aggregation across conditions (limit per index, keep/drop rule), per-NH-stratum recovery, host-filter-loss column.
  - SPIKE-02 / task I3 (2026-10-22) [~]: `aggregate.py` (subcommands `check`, `report`), `tests/test_spikein_aggregate.py` (13, synthetic score rows and bundles), `README.md` (run order, budget, decisions). Limit = lowest level met in >= 2 of 3 seeds with every higher level met; null floor from the unspiked replicates; flank verdict by the pre-set rule per background x layout; STAR fates by class from the bundles. `report` refuses (exit 2, writes nothing) unless all 66 conditions are DONE with valid scores. Decisions: "adds <= 5 unsupported UMIs" = per-replicate flank minus old, max (absolute flank max reported beside it); both limits not reached = tie, flagged; overall = retire if any background retires, else undetermined if any is, else keep; kb arms have class "all" only. NOT run on real data; budget in the README is the sbatch settings, not measured (explorer report not available to I3). The aggregation item of the I2 "Not done" line is now done.
  - SPIKE-02 / review fixes (2026-10-22) [~]: 7 reviewer findings checked against the code; all confirmed except that "extract_kb index is not guaranteed to be the gene id" was right only in part (multimap.py stores gene_id as a var column, index = kb matrix names; now matched on `var['gene_id']`). (1/3) kb call is now NET of the unspiked replicate: `aggregate.apply_baseline` subtracts the null replicate (same bg subsample, bg_rep = seed) molecules in the same spike barcodes and truth accessions, floored at 0, and recomputes call_kb/call_combined; done in aggregate, not `finish()`, because nulls are the LAST array rows and may not exist when a spiked condition finishes. `check` now requires the null kb tables; the report has a baseline section flagging any background whose baseline >= KB_MIN. (2) `extract_kb` fails when no viral column matches or (level >= 100) no spike barcode is in obs_names; strips a `-N` barcode suffix. (4) `verdict`: old index not reached in any background x layout is `undetermined` (never keep), headline says so. (5) STAR class table prints n, 95% Wilson interval, "insufficient" for n < 30, and states the pooling is over resamples of one source library. (6) `count_pairs` takes a per-background flock and writes the cache atomically. (7) new tests: DONE/stage skips, zero-STAR gate at level >= 100, called-list assert, frac < 1, count_pairs cache, extract_kb gates, baseline subtraction (50 spike-in tests pass). Mutations checked: removing the DONE skip, stage skip, called-list assert, level gate, frac subsampling each fail a test. NOT run on real data; the on-disk `adata_multimap.h5ad` layout is still from code and a synthetic AnnData only.
  - SMOKE RESULT (2026-10-21, job 25860984, 20 min, exit 0; D03 first 5 M pairs, HHV-6B spikes at 10 and 100, seed 1, 5 cells;
    one library, one seed, so descriptive only). STAR host+virus: 10/10 and 99/100 spikes on the source contig, 0 on host,
    1 lost to "unmapped: too short" (a spliced source read trimmed to 90 nt); the real STAR host filter kept 10/10 and
    100/100 (no host-homology loss); unspiked base: 0 viral STAR reads. kb, gene-level matrix: old index 6/10 and
    49/100, flank index 4/10 and 25/100. Per-spike BUS fate (`kb_fate.py`): pseudoaligned old 8/10 and 95/100, flank
    10/10 and 100/100, but multi-gene equivalence classes (dropped from the unique-gene matrix) old 46/100 vs flank
    75/100. So the flank index is MORE sensitive at pseudoalignment and LESS at the unique-gene matrix; ViralScan's
    multimap step redistributes ambiguous molecules, so the kb arm of the grid must be scored at BUS level
    (unique + ambiguous), not the raw matrix, or the verdict is an artefact of the metric (the pilot's
    `viral_umis.py` Arm-1 numbers are matrix-level). Unspiked-base viral UMIs: old 0, flank 3 (three anello genes).
  - `spikein.py synth` (circular-origin and short-fragment probes, Q4 option b): error-free reads from a circular
    genome record, origin-spanning vs internal control, and 30/50-nt short reads; unit-tested, not yet run on STAR/kb.
  - Note: one spiked read = one UMI (conservative; real data carry PCR duplicates). Smoke uses the first 5 M pairs,
    not a random subsample. The host-filter arm runs the real `_starsolo_filter` on the spike FASTQ alone (no viral
    competitor), compared with the full host+virus run. The R01 reverse-complement path is unit-tested but not yet
    run on real data.
  - Phase 5 code fixes (tooling only, still no compute run): `score.py star` adds `n_source_contig_any` (any
    alignment incl. secondary on `--contig`, the headline for the ">= 1 STAR read on the same contig" rule) and
    `n_other_virus`; `spike` takes `--whitelist` (on-list barcodes only) and avoids background UMIs of the chosen
    barcodes (one extra pass over bg R1, ~45 M pairs); `bg` gains seeded random `--frac`/`--seed` (`--n-pairs` stays
    first-N, smoke only); smoke.sbatch runs the spike-only arms first with a spike-barcode whitelist (STAR 139
    guard), prints host-filter accounting, and keeps VIRUS| + SPIKE BAMs before deleting full BAMs.
  - Caveats to settle before the grid: the NH:1-only source pool makes the limit optimistic (multimapping viral
    reads are excluded); the R01 arm places 5-prime-profile reads (GSM6416022 is 5-prime 10xv2, reverse) into a
    3-prime v3 background, so the +1 kb 3-prime flank class is under-exercised (decide before the grid); one spike
    read = one UMI; the smoke background is first-N, not random.
- WP-6/WP-8 deferred; WP-7 read-only. No scientific execution authorized.


- [~] `COVMAP-01` (2026-10-21, `codex/spq-coverage-map`): the evidence step's `coverage.raw_vs_deduplicated.png` is replaced by
  scPathoQuant-style coverage maps (`plot_coverage_maps` in `src/viralscan/evidence.py`, following
  `spq.visualization.viz.generate_coverage_maps`): per viral reference and BAM layer (`raw`, `deduplicated`) a `bamCoverage` bigWig,
  the GTF `gene` rows converted by BEDOPS `gtf2bed`, `make_tracks_file`, the same ini edit (labels on, grey becomes green) and
  `pyGenomeTracks` at 200 dpi, written as `<ref>_coveragemap.<layer>.{png,pdf}`. Two extras per layer: `<ref>_coveragemap_competition.<layer>`
  (one coverage track each for `unique`, `other_virus` and `host_competing` reads, classified from each read's `SA` partners, since the
  competitive BAM keeps no secondary alignments) and `<ref>_coverage_stacked.<layer>` (matplotlib stacked histogram, sqrt y-axis, packed
  gene track). New tools: `pygenometracks`, `deeptools`, `bedops` in the optional `environment.coverage.yml` (doctor `--profile coverage`, not `FULL_TOOLS`: CI libmamba could not co-solve them with `bustools=0.45.1` (hdf5 1.14.3 vs mpich/mvapich) and `r-base=4.5.3` in `environment.yml`, so they were moved out, as discovery did); the plot is diagnostic and a
  failure only logs a warning. Cross-virus multimapping is not drawn: Track A evidence BAMs have 0 secondaries because the pinned code
  (97a64f6) predates `AUD3-04` (`--secondary=yes`); see `AUD3-G3-EVIDENCE`. Interactive multimap functions (option C) are deferred.
  STAR/STARsolo BAMs (`codex/spq-coverage-map-star`, option A): STAR writes no `SA`; when the target's records have `NH`>1 and no `SA`,
  the target's multimapper QNAMEs go to a file, one `samtools view -N` pass collects their other records, and the same precedence
  applies (any `HOST|` partner, else another `VIRUS|` contig, else unique; same-contig hits stay unique). Streams records; memory is the
  multimapper QNAME set. Tested on a synthetic NH/HI BAM.
  **Open:** layout review by the user, a pyGenomeTracks run on a real library beyond HHV-6B AF157706.1, a `environment.coverage.yml` solve in CI.
- [~] `REF-11-PREP` (2026-10-23, `codex/panel-archive-script`, user: keep the script): `scripts/make_panel_archive.py` builds the
  deterministic viral-annotation panel archive for Zenodo (zip in the `viralscan data fetch` layout + `PANEL_MANIFEST.json`,
  `SHA256SUMS`, DRAFT `zenodo_metadata.draft.json` with TODO values). Preparation only: nothing uploaded, no DOI reserved.
  Checked: two builds of the local 393-GTF panel at 3.0.0rc1 gave byte-identical zips (sha256 5c5b5fea...49d9), corpus
  check passed; `tests/test_make_panel_archive.py`. **Open (user):** Zenodo deposit and the TODO metadata.
- [~] `RC1-01` (2026-10-23, `codex/release-3.0.0rc1`, decision Q15): version `3.0.0.dev1` -> `3.0.0rc1` in every file
  `scripts/check_versions.py` compares (`__init__`, Dockerfile, Singularity.def, conda recipe, CITATION.cff x2) plus the
  cli_reference and release_policy pins; CHANGELOG section renamed. Historical receipts (release_checks, PLAN history) unchanged.
  **Open (user):** merge after #109 + green main, then push tag `v3.0.0rc1` (TestPyPI per release.yml; not PyPI).
- [x] `VIG-RUN-01` (2026-10-23, `codex/quickstart-real-run`, decision Q10): the quickstart vignette now runs `viralscan run` when
  `VS_INDEX`/`VS_T2G` (and, for an index without a build manifest, `VS_GTF`) are set; the whitelist preflight needs `VS_WHITELIST`.
  Real Slurm run 25867622 (1M-read SRR12682296 cache, cat42d index + GTF, 8 CPUs, 3 min): 975 EmptyDrops cells, EBV in 199
  (20.4 %), `validate-run` ok at the tree root. Found: `kb info` exits 1 on success (kb-python 0.30.2), so the version
  cell no longer checks kb's status; an index not made by `build-ref` needs its GTF or `viralscan data fetch`. Method names
  (`em-global`) left as is by decision (accepted silently by the CLI).
- [~] `PATHS-01` (2026-10-22, `codex/strip-institutional-paths`): `PLAN.md` site paths replaced by `$SITE_ROOT` (archive work root), `$SITE_SCRATCH` (scratch checkout root) and `$COLLAB_ROOT` (collaborator root); 10 lines. **Open:** about 190 more hits in `.living/` and about 300 in tracked analysis outputs, manifests, docs and Slurm scripts (`git grep -c -E "/exports/|lumc"`); scope of a repo-wide scrub is the user's call. **Done (2026-10-22, `codex/strip-paths-docs-scripts`):** tracked `docs/`, `scripts/` and `covid_viralscan/` now use `$SITE_ROOT`/`$SITE_SCRATCH`/`$COLLAB_ROOT` placeholders in prose and the fail-fast env vars `VIRALSCAN_SITE_ROOT`, `VIRALSCAN_SITE_SCRATCH`, `VIRALSCAN_COLLAB_ROOT` in scripts (`#SBATCH` log paths are now relative; `--chdir` lines dropped, submit from the repo root). `.living/`, `analysis/`, `results/` and other outputs deliberately untouched; no history rewrite.
- [x] `RESTREE-01` (2026-10-21, `codex/research-tree-gitignore`): research-tree hygiene before any analysis/ trees are committed. `.gitignore` gains bulk rules (`analysis/computational_bottlenecks_*/{execution,followups}/`, `trackA_report_*/outputs/`, `*.rds *.mtx *.pkl *.stdout *.stderr *.sam *.bam *.bus` under analysis/, `/tmp/`, `.living/archive/`, `/biomni_bundle/`, root `HANDOFF*.md`/`SESSION_RESUME.md`); root-anchored where the name also exists tracked (`.living/HANDOFF_*.md`). New `scripts/check_tracked_file_sizes.py` + `tests/test_tracked_file_sizes.py` fail on any tracked/staged `analysis/` or `results/` blob over 1 MB unless listed in `config/large_file_allowlist.txt` (empty; passes on main as-is). Not hooked into pre-commit: the existing `check-added-large-files --maxkb=500` is already stricter repo-wide; CI runs the pytest. `sources/viral_panel`, `migration_plan/`, `viralscan_lshwcp_spec/` deliberately not covered.
- [x] `ANDET-09-STAR139` (2026-10-20, `codex/anello-star139`): the `anello_align` rule crashed on any run with
  `--soloCBwhitelist None` (chemistries without an on-list, SW-21) when no read mapped: STAR 2.7.11b segfaults (exit
  139) in the Solo step. Found on 30 hsv1 Drop-seq libraries in Track A; reproduced on a 500k-read subset (no-Solo and
  a dummy whitelist both exit 0; Unique-only and `Gene` still crash; peak RSS 1.7 GB, so not memory).
  `scripts/anello_align._run_star` retries once with a one-barcode dummy whitelist and accepts the result only if
  `Log.final.out` reports zero mapped reads (otherwise it raises, so no barcode is dropped from a library with
  alignments); the retry is recorded in `anello_align/star139_fallback.txt`. Any other failure, and exit 139 with a
  real whitelist, still raises. Tests: `tests/test_anello_star139.py` (4). The Track A wrapper
  `trackA_anello_align.py` carries the same logic and was confirmed end to end on hsv1 GSM3511319 (smoke job
  25852483: crash, retry, 22.0M reads, 0 mapped, classify completed). NOT tested: a no-whitelist library that has
  mapped reads and still crashes (it raises); STAR upstream bug report not filed.

- [x] `INV-REPIN-01` (2026-10-20, `codex/fix-inventory-pins`): three `analysis/v3_artifact_inventory.tsv` rows pinned
  pre-rebase commits that the rebase-merges of #77 and #78 left out of `main`'s history (`build-reference-command`
  04562cf, `claim-registry` 4c352b7, `public-readme` ab47b5a), so `test_repository_artifact_inventory_is_valid`
  fails in CI for any branch cut from `main` (local runs pass because the old objects still exist locally).
  Re-pinned to the equivalent `main` commits (51854ee, 16a05ec, e94b036); artifact sha256 values are unchanged.
  Lesson: a row pinned to a PR commit goes stale on rebase-merge; re-pin to the merged commit in a follow-up.
  Check: every `git_sha` in the inventory is an ancestor of `origin/main`.

- [x] `DOCS-CAVEATS-01` (2026-10-22, `codex/reference-caveats`, Q17): docs-only reference caveats. `docs/reference_panel.md` gains a
  'Known limitations' section (HSV-1 LAT two-exon model and poly(A) caveat; `kmer_twin_of` is a heuristic; blank `risk_class` = unknown);
  `docs/output_reference.md` `risk_class` row now says blank = unknown. No code change: code leaves `eve_risk`/`artifact_risk` unflagged
  for blank rows and falls back to a genus-name rule only for uncatalogued viruses, consistent with 'unknown'.
- [x] `DOCS-PERTPY-01` (2026-10-20, `codex/docs-pertpy`, stacked on `codex/orig-api-docs`): docs site follows
  pertpy's (https://pertpy.scverse.org/en/stable/) theme and landing structure, on user request. Theme
  `scanpydoc` (needs `sphinx-book-theme`, which scanpydoc does not declare), `sphinx-design` cards,
  `sphinx-copybutton`; landing page is intro, six cards, short getting-started, with the sidebar grouped
  General / Tutorials / About. Replaces `sphinx_clarity_theme` and its `beautifulsoup4==4.13.3` pin (pixi docs
  feature). Verified: `sphinx -W` build clean in a scratch venv (scanpydoc 0.15.4, Sphinx 9.0.4, book theme
  1.4.0). Prose rewritten (wording only; code fences, links, tables and headings verified identical):
  `index`, `installation`, `quickstart`, `migration`, `vignettes/README`, `api`,
  `reference_panel`. Deliberately NOT rewritten: `faq` (claim/inventory-pinned), generated `cli_reference` and `workflow_comparison`,
  and the specification pages (`output_reference`, `v3_counting_contract`, `panel_integrity`,
  `release_policy`, `reference_maintenance`, `changelog`), so no pinned claim changes meaning. NOT done: a
  logo, a funding block, `pixi lock`/`pixi install` for the docs env, a container rebuild. Possible doc bugs
  noticed, not fixed: `reference_panel.md` repeats `viralscan data fetch` under "Check the cache location",
  and lists `--anellovirus` where the `build-ref` example uses `--reference-panel anellovirus`.

- [x] `DOCS-ORIG-01`: README and Sphinx site restructured on the original upstream
  (mahfouzlab/ViralScan) skeleton: README Introduction / Installation / Input Data / User Guide
  (original-style short-flag examples, multi-sample) / License, with the v3 notice, claim markers
  and Limitations retained. Sphinx nav now Home / Tutorials / API / CLI and outputs; added
  `docs/api_reference.md` (autodoc, NumPy docstrings) and a `docs` extra (also installed by
  `.readthedocs.yaml`). `build_combined_reference` is not autodoc'd: its docstring is not
  reST-clean and `build_reference.py` is sha-pinned by claim `v3-reference-homology`. Vignettes are
  rendered unexecuted (`nbsphinx_execute = "never"`). The draft CAT-30 notice is in
  `exclude_patterns`, so `sphinx -W` is clean. Inventory row `public-readme` re-pinned to the new
  README (sha256 + git_sha `ab47b5a`).
  Follow-up: `docs/index.md` now mirrors the upstream docs/index.md heading order (Introduction /
  Installation / Input Data / User Guide / License), keeping the toctrees, migration/SUPPORT/SECURITY links.
  Theme switched to `sphinx_clarity_theme` (PyPI 2.1.0, Sphinx>=8; MIT/community licence requires the
  footer attribution, which the theme emits by default: do not remove). `pixi.toml` (mirrors environment.yml pins,
  dev/docs environments) and `Singularity.pixi.def` (dev image, re-solves) added with installation.md Options 3-4.
  Verified 2026-10-19: `apptainer build Singularity.pixi.def` (module apptainer 1.3.6, unprivileged) -> 1.5 GB SIF;
  in the image `viralscan --version` 3.0.0.dev1, kb-python 0.30.2 (kallisto 0.52.0, bustools 0.45.1), STAR 2.7.11b,
  samtools 1.23.1, minimap2 2.31, CD-HIT 4.8.1, DropletUtils 1.30.0, `doctor --profile full` ok, and `check-chemistry`
  on a 100k-pair slice of SRR20710641 calls 10xv2 5' (88.2% on-list). Build exposed that `pixi install` solves every
  environment: the image installs `-e default`, and `[feature.docs]` pins beautifulsoup4==4.13.3 (a Clarity dependency
  pin). NOT verified: `pixi lock` (uv cache lock held by another process); no pixi.lock committed.

- [x] `DOCS-REV-01` (2026-10-20, `codex/docs-revision`): documentation revision into three tasks (Run ViralScan,
  Inspect results, Assess the evidence), ported from the uncommitted opencode edits onto `main`. Keeps `main`'s
  `scanpydoc` theme and `api`/`api_reference` toctree; adds `inspect_results`, `assess_evidence`, `reference_help`,
  `project_info` (also in `MANIFEST.in` and `config/public_ship_scope.json`); `docs.yml` builds with
  `-W --keep-going`. `cli_reference.md` was not ported: `main`'s generated copy is current
  (`gen_cli_reference.py --check` clean). Strict `sphinx -W --keep-going` build passed locally only with the `alabaster` theme and without `sphinx_design`/`sphinx_copybutton` (not installed; no env has `scanpydoc`), so the first CI `docs.yml` run is the theme check. Plan: `docs/plans/2026-10-20-documentation-revision-plan.md`.
  2026-10-21 rebase onto `main` after CLI-UNIFY-01 (#87): every command example in the revised pages, README table and
  quickstart vignette moved to the `run`/`--r1`/`--r2`/`--cache-dir`/`--cell-type-labels`/`validate-run --run-dir` grammar;
  `faq-doc`/`output-reference-doc`/`public-readme` inventory rows re-pinned in a final commit. Lesson (INV-REPIN-01): pins
  name the last content commit of this branch, which survives only a merge-commit PR; after a rebase/squash merge re-pin to the merged commit.

- [~] `DOCS-PAGES-01` (2026-10-20, `codex/docs-pages`): publish the Sphinx docs on GitHub Pages (user choice over
  Read the Docs). `.github/workflows/docs.yml` builds on push to `main` (and by hand), installs pandoc for nbsphinx
  and `.[docs]`, and deploys with `actions/upload-pages-artifact` v5.0.0 / `actions/deploy-pages` v5.0.1 (SHA-pinned,
  `check_actions_pinned.py` clean). Not built locally (no env has `scanpydoc`). Left for the user: Settings -> Pages ->
  Source "GitHub Actions", push, merge; `[x]` once the first deploy is live. `.readthedocs.yaml` is kept.

- [x] `CLI-UNIFY-01` (2026-10-20, `codex/cli-unify`): one CLI grammar (user decisions 2026-10-20: hard break,
  `viralscan run` subcommand, post-hoc outputs inside the run). PR A: `viralscan run`; `--r1/--r2` replace
  `-s1/-s2`; no single-dash long flags; `BooleanOptionalAction` replaces `_OptionalBool`; `--cache-dir`,
  `--cell-type-labels`, `--hostresponse-mito-control`, `--use-hvg/--no-use-hvg`, `--min-umi`, `--run-dir` on
  validate-run, `--json PATH` everywhere, `--cores/-c` and `--verbose/--quiet` on every subcommand. Removed spellings
  error with "use X instead"; a bare `viralscan` warns and runs `run`. `dest` names unchanged (config.yaml,
  fingerprint, resume unaffected). Decisions: `build-ref --anellovirus` stays (it adds the anellovirus table to a
  host+virus build; `--reference-panel anellovirus` builds an anellovirus-only one, so not a duplicate);
  `--cores` added only to `build-ref` (no flag that does nothing elsewhere); `run` makes `-o/--r1/--r2`
  argparse-required (the `_die` checks went); `hostresponse -o` (really an input run dir) is left to
  CLI-UNIFY-02. Not touched: `scripts/run_fresh_control.py` (drives frozen pinned v2/v3 binaries), historical
  `scripts/slurm_*.sh`, `tests/_frozen_oracles/`, `.living/`/`analysis/` notes, `docs/write-docs-prompt.md`.
  Supersedes `CLI-ORIG-01`. Tests: `tests/test_cli_unify.py`; migration table in `docs/migration.md`.
- [x] `CLI-UNIFY-02` (2026-10-21, `codex/cli-unify-02-layout`, stacked on `codex/cli-unify-02-hostresponse`):
  `hostresponse` takes `--run-dir` (input) + optional `--output`; user decision 2026-10-21: post-hoc outputs
  default INSIDE the run in a reserved `<run-dir>/posthoc/` (one folder per verb), `--output` overrides,
  `layout_version` 2. Design: `outputs.posthoc_dir/is_posthoc/SUPPORTED_LAYOUT_VERSIONS`. `evidence` ->
  `posthoc/evidence/<virus>` (`review_candidates` for `--review-candidates`, whose TSV moved out of `results/`);
  `rerun-multimap` -> `posthoc/multimap_<method>` (copy skips `posthoc/`; non-empty output still refused; an
  explicit `-o` inside the run but outside `posthoc/` is refused); `rerun-programs` -> `posthoc/programs`
  (`gene_programs.main(results_root=)`; no `log/gene_programs.done`, no marker re-stamp); `hostresponse` ->
  `posthoc/hostresponse` (no marker re-stamp); `evidence-plot` unchanged (writes beside its evidence dir).
  `validate-run` accepts layout 1, 2 and absent, errors on others, and its rglobs skip `posthoc/`. `run
  --hostresponse` / `--auto-evidence` are part of the run and keep their paths. Tests:
  `tests/test_posthoc_layout.py` (hashes the run before/after). Docs: migration, output_reference, faq, CHANGELOG,
  cli_reference. Not re-pinned: `analysis/v3_artifact_inventory.tsv`, claim registry.
  - Note (2026-10-22, PR #99 review): `run` now refuses a derived sample ID `posthoc` (it would be read as
    post-hoc output by `validate-run` and the `rerun-multimap` copy); `menu._check_sample_ids` holds this and the
    duplicate check. Hand-written `cli_reference` sections for evidence / rerun-multimap / hostresponse /
    rerun-programs updated to the layout 2 defaults.

- [x] `OUT-LAYOUT-01` (2026-10-20, `codex/output-layout`, stacked on `codex/orig-api-docs`): output
  layout kept close to evonk's 2.2.0 (user decisions, grill 2026-10-20). `count_audit.tsv` moves into
  `results/` (path owned by `viralscan.outputs.V3Outputs`); `run_manifest.json` gains `layout_version: 1`
  (outside the fingerprint); `combined_summary.tsv` at the run root; `summary.txt` opens with the original's
  block, v3 lines after. Not done by decision: no `v3/` wrapper, manifest stays at the root (pinned by
  `test_run_manifest_is_at_the_tree_root`), banner untouched (already identical; `prog` and help text still
  differ from upstream). Decision: candidate-only viruses get no "has a viral load of" line in the original
  block (would overstate an ungated virus). Old runs: `count_audit.tsv` stays at the sample root there.

- [x] `CLI-ORIG-01` (superseded 2026-10-20 by `CLI-UNIFY-01`, user chose a hard break; the forms now fail
  with the replacement named and `tests/test_original_cli_compat.py` is removed): accept the original upstream (mahfouzlab/ViralScan) invocation forms:
  `-reference` alias and optional True/False values on `--reference`, `--visual`, `--multimapping`
  (strict `_OptionalBool` action, never `type=bool`). Defaults and v3 contracts unchanged. Tests:
  `tests/test_original_cli_compat.py`. `-umap` stays a bare flag. Known limit: a bare flag directly before a subcommand name takes it
  as the value (`-v data` is a usage error); documented, pinned by test.

- [~] `SCAUTOQC-01` (2026-10-19, branch `codex/scautoqc`, ported from `codex/scautoqc-wip` 92dd88e): optional 2026-10-22: rebased onto origin/main via tmp/wave2-base-2026-10-22; tests adapted to the `run` subcommand.
  `--scautoqc` SCTK cell/cluster QC flags added in place to the count H5ADs and UMAP outputs; cells are flagged,
  never removed. Opt-in rule/flag/config key are absent when off. Conservative defaults on the two open design
  questions (user to confirm): (a) the legacy n_counts/n_genes filter in `umap.py` is UNCHANGED with the flag
  (the WIP skipped it); (b) the rule runs after multimap and before cell calling, so it scores every barcode in
  the kb H5AD, and that cell set is recorded in `adata.uns['scautoqc']` and `log/scautoqc.done` (JSON).
  SCTK ships as the optional extra `viralscan[scautoqc]`, not in `full`. Tests: flag-off, all-cells-fail guard,
  flags-only, missing-sctk error, DAG dry-run. **Open gate:** real FASTQ-to-UMAP validation on SRR8315713 (needs
  kb/kallisto/bustools env, reference and user-approved compute) is NOT done; clusterwise multi-resolution QC
  (absent in SCTK 0.2.0) and ribo/globin marker sensitivity are untested. Flip to `[x]` only after that gate.
  Rebased onto origin/main + shifted-log + tchem (wave-2 prep stack, not upstream): `viralscan.scautoqc` is under mypy
  strict and `viralscan.scripts.scautoqc` under the relaxed Snakemake-script override; legacy filter kept.
  Review follow-up: scoring all barcodes (empty droplets included) is a KNOWN LIMITATION, not mitigated; flags on a
  raw 10x matrix are not meaningful and clustering may be slow until the rule is moved after cell calling (gated
  on the real run). Added BioMart failure handling + low-mapping warning, case-insensitive `MT-`, extra tests;
  virus-only/no-host index still hard-fails (documented). Clustering path needs `scanpy[leiden]` (the extra).
  - Note (2026-10-22, PR #100 review): globin markers matched human symbols only, so every mouse `--scautoqc` run
    aborted on the mandatory hb check; `GLOBIN_PATTERN` now also covers mouse (Hba-a1/a2, Hba-x, Hbb-bs/bt/b1/b2,
    Hbb-bh1, Hbb-y, Hbq1a/b), with a regex test that needs no SCTK.
  - Note (2026-10-22, grilling Q8): CI now installs `scanpy[leiden]`, `sctk==0.2.0` and `pybiomart` in the py3.11 Ubuntu
    pytest job only, so `tests/test_scautoqc.py` runs in CI instead of skipping. Not py3.12: sctk 0.2.0's numpy-groupies
    pin fails to build there (`configparser.SafeConfigParser` removed), so the `[scautoqc]` extra is py<=3.11 in practice.

- [x] `PERF-PLAN-02`: draft and review actionable cards for all 27 remaining
  performance candidates; preserve completed B05/B07/B20. Plan:
  `docs/plans/2026-10-17-remaining-performance-candidates.md`. (Row text carried over
  from the main checkout's uncommitted PLAN.md; not re-reviewed here.)

- [x] `PERF-EXEC-01`: first-wave B07/B20/B05 complete; E01–E05 passed and E06
  assessed/deferred. 1,055 tests passed, one R skip; 125 independent exactness tests
  reverified. Combined Python workload 4.21x faster; no whole-run speedup claim.
  Plan: `docs/plans/2026-10-17-performance-triage-byte-identical-outputs.md`. (Carried
  over from the main checkout's uncommitted PLAN.md, authored in an earlier session;
  not re-reviewed here.)

- [~] `PERF-EXEC-02` (2026-10-18, branch `codex/perf-wave1-b24-b23-b04-b06`): wave-1 B24
  (reference low-complexity counters once), B23 (blocked programme calls), B04 (local
  CSR column view) and B06 (row-blocked dense reduction) implemented byte-identical and
  qualified together (C00, component level). Combined vs baseline `file_manifest`
  identical (148 artifacts), ruff/mypy clean. Benefit in one session on a loaded host: B24 gate ~0.36x time;
  B06 0.09-0.42x isolated RSS on wide inputs, time neutral; B23 sparse-negative 0.20x traced
  peak, time inconclusive; B04 time inconclusive, +4 MB view cost. **Rebased onto `main`**
  (`6910c2d`; the original base was 188 files behind) and landed there: one real source conflict in
  `build_reference.py` (AUD3-05 x B24), the tracker/governance conflicts, and both governance
  digests recomputed at the merged bytes (`6c8cb96`, `f852214`). That rebase surfaced a
  **merge-induced defect**: `main`'s Q123 positive-cell count indexed the layer that B05/B07/B20's
  budget-gated view had rebound to a narrow span, which raised `IndexError` on any real run (the
  reference puts host genes first, so viral genes always trail) — fixed in `5d7de93`, with the
  frozen evidence oracles re-pinned onto `main` in `b6bec78` (user-approved; nothing loosened).
  Rebased-tip suite **5,512 passed, 41 skipped, 139 deselected, 6 failed**; the six are the
  pre-existing, environmental `test_release_static.py` failures that also fail on a pure
  `origin/main` checkout. Not done: real
  dustmasker/`kb ref`/R integration, a second timing session, whole-run measurement,
  remaining cards. Tracker and evidence:
  `docs/plans/2026-10-17-remaining-performance-candidates.md`,
  `analysis/computational_bottlenecks_2026-10-17/followups/20261017-122240-G00/c00/C00_REPORT.md`.
  **Next up:** review of the PR; a quiet-host second timing session; land wave 1. B08 is
  deferred on its gate and B10 plus wave 2 (B25, B27, B18) are in flight in the
  `20261018-112322-wave2-g00` packet, on `codex/perf-wave2-b25-b27-b18`; that session must re-check
  B18's oracle against `main`'s post-merge `molecule_verdict.py` before trusting it.

- [x] `PERF-DOC-01` (2026-10-18, same branch): reconcile the tracker and wave-1 plan against the
  receipts — header status, the stale "G00/G01/G02/C00 not started" bullet, the certified final
  tip-suite result, and the empty `C00_REPORT.md` §3b; record B08's deferral and the wave-2 in-flight
  scope without rewriting cards the wave-2 branch owns. No code, golden, comparator or gate verdict
  changed; no card newly claimed as benefit-qualified. Tracker:
  `docs/plans/2026-10-17-remaining-performance-candidates.md` (Decision and evidence log).

**2026-10-15 (audit-fix integration, `codex/aud3-integration`):** the independently verified audit findings are fixed in per-unit `codex/aud3-*` branches merged here (`AUD3-01a/01b/02/03/04/05/06a/06b/07/09a/09b`), with one governance re-pin (`AUD3-GOV`). Verification report and per-finding matrix: `analysis/audit_2026-10-verification/`. Open: `AUD3-02-OPEN` (viral-free UMAP null still above nominal), `AUD3-FLAKY`. `AUD3-08` is done (see its `[x]` row below). PR #71 merged. AUD3-04 and the `harmonization`-adjacent items land before the G3 signature.

- [x] `AUD3-01a` — Virus Identity table vs matrix reconcile (audit P9-02, verified). One shared helper
  `virus_grouping.reconcile_identity_to_matrix` (raises `IdentityMatrixMismatch`) now runs in `detection.preprocessing`
  and `multimap` and replaces the inline CELLS-05 check in `cellcalling.viral_partition` (same message, wrapped in
  `CellCallingError`). A matrix column the table does not list stops the run (a table viral gene renamed in the matrix,
  e.g. `vB_g1` vs `vB_g1.2`, previously dropped the detection silently on `--cell-calling none`/external/two-step); viral
  table genes absent from the matrix log one warning, no error (no counts can exist). Detection thresholds, layers and
  denominators are unchanged. Docs: `CLAUDE.md`/`PLAN.md` branch line now says per-unit `codex/<unit>` branches (the
  `codex/viralscan-v3` branch never existed); `docs/output_reference.md` sentence. Governance pins to re-pin (not edited
  here): inventory row 17 `output-reference-doc` (docs/output_reference.md; `tests/test_artifact_inventory.py` fails on the stale sha256 until the final re-pin); CLAUDE.md/PLAN.md/CHANGELOG.md have no row.
  Tests: `tests/test_identity_matrix_reconcile.py`.
- [x] `AUD3-01b` — detection gate float tolerance (audit P9-03). `detect_genes` and the called-cell gate in `compute_stats` compared summed viral mass with a bare `>=`; both now use `threshold * (1 - INFECTED_RELATIVE_TOLERANCE)` like `infected_masks`
  (no new constant; locked decision 2026-10-15). A molecule split 6 or 7 ways sums to just under 1 in float64 and was
  silently dropped everywhere (no row, no CANDIDATE ONLY). Thresholds, layers and denominators are unchanged: only float
  under-shoot is removed. Other gates reviewed: `infected_masks`, `tier_measures`, `evidence_tiers` already carry the
  tolerance; `sensitivity.py` default `detected` inference and `hostresponse._virus_presence_label` (per-cell label, not
  the detection gate) are left bare on purpose. Docs: `docs/output_reference.md` (`detected` row). Governance pins to
  re-pin (not edited here): inventory row `output-reference-doc` (docs/output_reference.md).
  Tests: `tests/test_detection_gate_tolerance.py`.
- [x] `AUD3-03` — `report.html` run summary and main table headline `pct_infected_called` (comparable-cell fallback when there are no called cells), all-barcode rate secondary, badge on the primary rate; in two-step (`--host-filter starsolo`) runs `pct_infected` and `viral_molecules_per_10k_est` are NA (blank in the TSV) with a new `denominator_note` column (audit P3-03, P4-02). Changes: `detection.py` (`compute_stats`, `primary_rate`, `_virus_summary_line`, `write_tsv_outputs`, `generate_html_report`), `report.html.j2`, `docs/output_reference.md`, `docs/faq.md`, CHANGELOG. Tests: `tests/test_headline_denominators.py`. `pct_infected_called` / `pct_infected_comparable` / called-cell anchor and the combined-mode values are unchanged; no viral_summary schema exists. Pinned docs touched (re-pin in the inventory unit): inventory rows `faq-doc` and `output-reference-doc`. Not changed: the summary.txt headline count `x/y all barcodes` still shows the host-filtered barcode count in two-step runs (a count, not a rate).

**2026-10-15 (`SE-01`):** super-expressor summary script added (Lareau 2023 definition, >=10 viral UMIs); it runs on every Track A library after arm N. Track A full-depth round 1 is set up under `viralscan_work/trackA_round1`.

**2026-10-14 (host-only cell calling, `codex/host-only-cell-calling`):** emptyDrops calls cells on the host
transcriptome only, one call per library reused by every analysis (Q153); combined mode now calls on the host
gene columns of the kb matrix. Done as `CELLS-02`; `CELLS-01` amended. G3 review of #64 found the seed overflow
(Q154); fixed as `CELLS-05` (`codex/cell-calling-review-fixes`). Next: merge, then the G3 signing procedure.

**2026-10-14 (independent UMI/Q137 review):** verify the resumed UMI branch
without editing Claude's shared worktree. Q137 trimmed-query coverage must count
only aligned M/=/X positions inside the trimmed body; adaptor and tail positions
must not inflate its numerator. Thresholds and the flag-only policy are unchanged.

- [x] `AUD3-05` — reference near-duplicate screen at build time (audit P2-01, P2-04, P2-10): strand-aware k-mer twin screen between reference records (reuses `KMER_TWIN_BELOW` / `canonical_kmers` / `sharing`), warn-only, `--strict-reference` fails on non-sibling pairs, visible host-homology-not-measured warning.
  Reuse: `canonical_kmers`, `sharing`, `_CODE` and `KMER_TWIN_BELOW` moved to `src/viralscan/kmer_sharing.py`
  (re-exported by `scripts/panel_kmer_sharing.py` / `panel_candidates.py`, their tests unchanged); no new cut.
  `validate_reference_records(near_duplicates="warn"|"error"|"off")` (default warn); `build-ref` and the manifest
  writer screen once (viral-involving pairs only), `--strict-reference` on `build-ref` and
  `build_bundled_panel_ref.py` raises `NearDuplicateError`; manifest/bundled build log "host homology NOT measured".
  Decision: the 10%-diverged copy is NOT flagged (about 9% of 31-mers survive vs the 95% cut); deferred to `REF-08`
  (test pins "not flagged"). Only the top partner is judged; declared `sibling_group` pairs are exempt.
  Tests: `tests/test_reference_near_duplicate.py`. Docs: `docs/panel_integrity.md`, `cli_reference.md` (regenerated),
  `quickstart.md`, `output_reference.md` L246 (pinned by inventory row 1.1.0 `output-reference-doc`; re-pin needed).
- [x] `EVID-Q137-COVERAGE` — correct trimmed-query coverage; regressions cover
  leading TSO, trailing poly-A, reverse strand and an intact aligned body.
  Validation: 187 focused unit/integration tests passed; lint and typing passed.
  The 14 saved torque-teno calls remain at zero passing molecules; saved HPV16
  remains 2,733 reads / 687 molecules / 679 passing provisional gates. Source
  BAM/FASTA/manifest/lineage hashes were unchanged. Full-suite qualification is
  recorded separately; no calibration, G3 signature or detection threshold changed.
- [x] `AUD3-04` — audit P3-02: the competitive `minimap2 -ax sr` had `--secondary=no` (sr preset), so exact host/virus ties
  could never be `tie` (200 perfect ties: 105 host_best / 95 virus_best / 0 tie before, 200 tie after). `evidence.MINIMAP2_ALIGN_OPTIONS`
  = `-ax sr --secondary=yes -p 0.8 -N 5` now used by `align_reads_to_viral`; `molecule_verdict.py` docstring and
  `docs/output_reference.md` say what the code does. minimap2 2.31 `--help`: `-p` = min secondary-to-primary score ratio,
  `-N` = retain at most N secondaries, so the docs' "limit 5" meant `-N 5` (the `sr` preset itself is `-p 0.5 -N 20`).
  The evidence manifest gets a top-level `minimap2_options` and the minimap2 version string now also for host-verdict runs
  (it was `null` there; only competitor mode filled it). `verdict_of`, tier cutoffs, thresholds and the flag-only policy untouched.
  Tests: integration `test_minimap2_secondary_ties.py` (fails before), `test_align_reads_to_viral_minimap2_argv_is_pinned` (fails
  before), `test_manifest_records_minimap2_version_and_options_without_competition`; 236 existing evidence/verdict/read-support/tier
  tests pass. Replay of the 15 saved `dsr_round2` calls (old flags reproduce the saved BAMs record-for-record): no verdict or
  passing-molecule change (HPV16 2,733 reads / 687 molecules / 679 passing, 680 virus_best, 0 tie; 14 torque-teno calls 0 passing);
  BAM records 2,736 -> 2,805 for HPV16. Evidence made before this row has no ties and its manifest lacks `minimap2_options`:
  re-run before feeding `molecule_verdicts.tsv.gz` to tier measure S/G3. Pinned inventory row affected:
  `analysis/v3_artifact_inventory.tsv` `output-reference-doc` (`docs/output_reference.md`); re-pin left to the final unit.
- [x] `AUD3-GOV` — governance re-pin after the AUD3 units: inventory rows `docs/faq.md`, `docs/output_reference.md`, `README.md`, `docs/v3_counting_contract.md` and `src/viralscan/scripts/build_reference.py` re-pinned (`artifact_sha256`, `git_sha` and the `git show` retrieval command only); `scripts/validate_artifact_inventory.py` passes. `claims/registry.json` unchanged (no registered claim changed).
- [x] `AUD3-08` — resume after a purged `kb-python/`: verified read-only (2026-10-15) that `viralscan --resume` with `log/kb.done` present and `kb-python/` deleted exits 0, runs no `kb count`, and rewrites `run_complete.json` without the counts entry (`Snakefile` `kb_count` declares only `log/kb.done`; `run_safety.prepare_output_directory` checks the manifest fingerprint only). Smallest fix: on resume require `kb-python/counts_unfiltered/adata.h5ad` for every sample with `log/kb.done`, else stop ("delete log/kb.done or use --overwrite"); or declare the file as a `kb_count` output. Decided 2026-10-15 (user): option (a), the resume-time guard, as its own small unit after the audit-fix PR. **Implemented 2026-10-15** (`codex/aud3-08-resume-guard`): `run_safety.prepare_output_directory` refuses `--resume` when `<sample>/log/kb.done` exists without `kb-python/counts_unfiltered/adata.h5ad`, before the completion marker is cleared; `tests/test_run_safety.py::test_resume_refuses_purged_counts`. `viralscan validate-run` already reports it.
- [ ] `AUD3-FLAKY` — `tests/test_cellcalling.py::test_real_r_frozen_protocol_seed_runs_and_is_reduced_modulo_int_max` failed once when R tests ran concurrently with other pytest processes and passed on rerun; investigate (R/DropletUtils temp state or timeout) before relying on parallel suite runs.
- [ ] `AUD3-G3-EVIDENCE` — evidence made before `AUD3-04` (the 35 `trackA_round1` evidence directories and the `dsr_round2` calls) used `minimap2 -ax sr` without secondary alignments, so it has no `tie` verdicts and its manifests lack `minimap2_options`. Replaying the 15 saved `dsr_round2` calls showed no change in verdict counts or passing molecules (none holds an exact tie), but re-run the evidence before `molecule_verdicts.tsv.gz` feeds tier measure S at G3 preparation.
- [ ] `AUD3-REF08-DIVERGED` — audit P2-01 residual: the build-time near-duplicate screen (`AUD3-05`) reuses the curation twin cut (>= 95% shared canonical 31-mers), so a host homologue diverged by about 8% (about 90% of its k-mers unshared at k=31) is not screened; `host_homology_status=not_measured` stays the honest label. A homology threshold or alignment-identity cut belongs to `REF-08` (CAT-14 deliberately has none); do not invent one here.
- [x] `DOC-API-01` — `docs/api.md` cell-type enrichment example is wrong in two ways (found 2026-10-15 during the analysis-layout grill): (1) it passes a dict as `cfg` (`docs/api.md:54-56`), but `cell_type_enrichment()` reads `cfg.cell_types` and `getattr(cfg, "infected_min_molecules", ...)` (`src/viralscan/enrichment.py:49`, `:114-115`), so the example raises `AttributeError`; pass `RunConfig` (as the CI vignettes already do); (2) the text says a cell counts as infected when the virus sum is greater than zero (`docs/api.md:40-41`), but the code applies the `infected_min_molecules` threshold shared with detection. Docs-only fix; no behaviour change. **Fixed in `docs/api.md` (2026-10-19, folded into `codex/orig-api-docs`; the standalone `codex/plan-doc-api-enrichment` branch is retired).**

**2026-10-14 (UMI deduplication and Q137 read-support flags, `codex/umi-dedup-q137`):** every count that feeds a verdict, flag, threshold, tier measure, summary number or plot is a UMI-deduplicated molecule count (corrected CB+UMI, deduplicated AFTER the competitive alignment so the representative read keeps its host and virus records); raw read counts stay only where labelled reads (alignment_qc `raw` layer, `read_verdicts`, host_capture, anello `*_reads`). Audit table in the PR body. Q137 (user-approved Biomni gates) are flags in `molecule_support.tsv` / `read_support_summary.tsv`, never removal. Next: PR review; `protocol.yaml` and `v3_score.py` untouched (their counts are already molecule counts); thresholds stay provisional until the calibration gate.

- [x] `EVID-UMI-01` — alignment_qc decisions, `interpretation_flags`, BLAST sample, DSR verdicts/cards/adjudicator and the figure use the deduplicated layer; regression `tests/test_read_support.py` (5 reads of one molecule + 1 of another = 2 molecules everywhere).
- [x] `EVID-Q137` — `read_support.py`: aligned length / trimmed-query coverage, DUST (own symmetric DUST, equals sdust), A/T warning, competitive human, hotspot (data TSV `read_support_loci.tsv`) and DNAJB6 decoy, per-call `known_artefact`; `evidence-plot` panel D (aligned-length histogram, 60-base gate), `--layer umi|raw`, lanes only up to 60 molecules.

**2026-10-14 (arm P prep, `codex/translated-search-arm-p`):** arm P is prepared for every DSR library (`AA-03`) and piloted on one 5' library (`AA-04`); the shared index's host mask was found to be the plain one (`AA-04`), so **next up: decide `AA-05` (rebuild with the encoded mask) before any further arm-P runs**; `AA-02` stays `[~]`.

**2026-10-14/15 (combined PalmDB reference, `codex/translated-search-palmdb`, Q155):** `AA-05` is superseded; the arm-P reference is PalmDB RdRP + the annotated proteins in one kallisto 0.50.0 index (`AA-06`, built and smoke-tested on EBV). **Next up:** the full-depth arm-P runs on all libraries against this index (submission pending; see `AA-06`), then read-level confirmation of suspicious hits.

**2026-10-14 (user-requested translated protein search):** prepare an all-protein
reference for the viruses in the local GTF corpus (RNA and DNA viruses), and
port the kallisto translated-search tutorial with explicit attribution.
- [x] `AA-01` — extract accession-matched, catalogue-hash-verified GenBank CDS
  translations; preserve missing annotations in a coverage table; prepare
  Snakemake host-masked indexing/counting and verify a small real-tool example.
  This is a separate exploratory reference, with no change to v3 detection claims.
  **Preparation verified 2026-10-14:** 4,124 proteins from 389/393 GTF
  accessions cover 296/299 catalogue labels (111 RNA, 185 DNA); four annotation
  gaps are retained in `coverage.tsv`. Four focused tests, lint/format, and data
  governance pass. The Snakemake reference/index/count workflow passes a real
  kallisto 0.52.0/kb-python 0.30.2 smoke test: 20/20 synthetic UMIs each for
  Aichi virus and HSV-2, using a synthetic host mask. The real human-mask build
  passes dry-run only; real samples, specificity/sensitivity, and reactivation
  qualification remain unrun. Attribution and run instructions:
  `analysis/translated_search/README.md`.
  **Snakefile provenance:** `validation_2026-10-14.json` and the earliest ledger rows describe the
  pre-kallisto-0.50.0 `Snakefile` (sha256 `acafbd0a…`); the committed `Snakefile` (sha256 `3f9ba67c…`,
  edited for the kallisto 0.50.0 retry) is the version that built the host-masked index and ran
  AA-02 (SLURM 25772166, both counts COMPLETED). The JSON is intentionally unchanged.
- [!] `AA-02` — **Deferred 2026-10-20 (user parked translated-search work behind Track A on 2026-10-19; a full-depth run is not authorized).** user-authorized bounded exploratory sanity check: build the
  human-genome/cDNA-masked protein index, then search the first 1,000,000 intact
  paired reads from public EBV-positive SRR12682296 (10xv2) and HSV-1-positive
  SRR8315713 (Drop-seq). Keep immutable subsets, input/reference hashes,
  exploration ledger, expected-virus and off-target summaries. No full-depth
  analysis, held-out negatives or calibration is authorized here. (The Codex-era "no
  commit/push" line is superseded for the code only; see the takeover note in the
  "VR pilot" section.) **Note 2026-10-14 (`AA-04`):** the index behind these counts used the plain, untranslated
  host mask, so its host masking was largely inactive; the AA-02 numbers are unmasked.
- [x] `AA-03` — arm P (translated search) prepared for every DSR library (`codex/translated-search-arm-p`,
  stacked on the VR-pilot PR; nothing merged). `analysis/translated_search/make_configs.py` writes one config per
  technology/strand/whitelist group from the DSR manifest and chem JSONs (chem JSON is authoritative for
  technology, so GSM4682311 is 10xv2; each R1/R2 pair is opened and R1 must have the chem JSON's `r1_length`);
  optional Snakefile `whitelist` key adds `kb count -w` (default behaviour unchanged; needed by the three
  user-whitelist libraries); `summarize_translated.py` (generic AA-02 summary, `--cells`, expected taxa from
  `expected_taxa.tsv`) reproduces the EBV 1M numbers (2,312 UMIs, 798 barcodes); `check_barcodes.py` (overlap and
  format differences vs a frozen cell list, abort below 90%). Eight tests in
  `tests/test_translated_search_tools.py`. Decision: the manifest already lists ERR13027027's `_2` file (28 bp) as R1, the
  barcode read, so no swap is applied.
- [!] `AA-04` — **Deferred 2026-10-20 (user parked translated-search work behind Track A on 2026-10-19; a full-depth run is not authorized); full-depth arm-P runs need explicit approval.** first 5' pilot of arm P (first 1,000,000 intact pairs of one 10xv3-reverse-strand institutional
  library; local only, results not committed): `kb count --aa --parity single -x 10xv3 -w <user whitelist>` runs
  to completion on 2 CPU / 8 GB in minutes with R1 = 28 bp barcode+UMI, so `--parity single` works for 5' reads.
  **Index provenance finding:** the shared kallisto 0.50.0 index used the plain nucleotide host mask, not the
  encoded one (the retry note is wrong); 0.50.0 does not translate the D-list, so host masking is largely inactive
  (details: `analysis/translated_search/README.md`, "Index provenance"). Per-virus output of this index is
  therefore unmasked and exploratory; it is not a specificity result.
- [!] `AA-05` — rebuild the protein index with the CFC-encoded host mask under kallisto 0.50.0. **Superseded
  2026-10-14 (Q155):** both original workflows (the kallisto tutorial and `LSCHWCP_2023` notebook
  `7_virus_host_captured_dlist_cdna_dna`) use a plain nucleotide cDNA+DNA FASTA as the D-list, so the plain mask is
  the intended one; the encoded mask is not built.
- [x] `AA-06` — combined translated-search reference: PalmDB RdRP (the tutorial's three files, sha256 pinned;
  BSD-2-Clause repo, CC0 PalmDB source) + the 4,124 annotated proteins in one index. `scripts/combine_palmdb_reference.py`
  (+ `tests/test_combine_palmdb_reference.py`) writes the FASTA, t2g (`taxon_<taxid>` / `palmdb_<rep_ID>`), a taxonomy table
  with catalogue taxon ids and the overlap report; `count_library.py` runs `kb count --aa --num` for one library.
  Index: kallisto 0.50.0, `-k 31`, plain host mask (sha256 `8c788453...601b`), 606 flanking k-mers, SLURM 25778661
  (35 min, 4.6 GB). EBV 1M smoke: 1.4 % pseudoaligned, `taxon_10376` 2,311 UMIs / 797 barcodes (v0500: 2,312 / 798);
  PalmDB adds 344 UMIs, none EBV (PalmDB has no DNA viruses). Not a validation. Numbers: `analysis/translated_search/palmdb_provenance.json`
  and the README. **Open:** full-depth runs on all 17 libraries (jobs prepared, not submitted: the permission
  system blocked the submission) and read-level confirmation of suspicious hits.
- [x] `AA-07` — read-level confirmation and artefact checks for arm-P hits (`codex/translated-search-confirm`, stacked on
  `codex/translated-search-arm-p`; pilot data only, 1M pairs). `analysis/translated_search/confirm_hits.py`:
  `kallisto bus --aa --num` keeps the read number in the BUS flag (works in 0.50.0; the corrected BUS file has the same
  (barcode, UMI, EC) set and read total as the original), then per group up to 2,000 reads are aligned (minimap2) to the
  host genome, host transcriptome and viral genomes, blastx'd against viral and human proteins (Ensembl 116 peptide),
  checked for DUST/poly-A/G/TSO, given a per-read and per-group verdict, and summarised by barcode pattern; host
  capture is a measurement (share of reads that also align the host), not a filter. `recurrence` lists proteins hit in
  several libraries. Findings: the EBV expected-virus controls (BHLF1, BNLF2b, UL43) are viral_supported; HCMV pp71 and
  HPV83 in the HSV-1 pilot are host-derived; the dominant 5' institutional hit (HPV83 E1) is TSO/poly(T) adapter
  artefact that capture would not remove. Decisions: DIAMOND and bwa are not installed, so blastx and minimap2 are used;
  minimap2 runs with `-I 300M` (1G exceeded 8 GB); `adapter_artefact` and `protein_only_viral` verdicts were added to the
  requested set. Negative-control subsets (HPV16 parental, HNSCC, tonsil atlas: FASTQs present) are prepared, not run.

**2026-10-14 (tier-cutoff measures, `codex/tier-measures`, Q140-Q148):** `TIER-02` builds the three frozen `probable` measures (m, F, S; Q141) in `src/viralscan/tier_measures.py`; F is wired into `viral_summary.tsv`. Training blocker `tier_cutoff_measures` is NOT lifted (no evidence-run driver per sample x virus yet; see `TIER-02`). Next: that driver, then the emptyDrops anchor (`real_background_anchor`, Q142), after the round-15 PR merges.

**2026-10-14 (SCI-05 round 15 fixes, `codex/sci05-r15-fixes`, `DEV-032`, Q139):** round 15 (independent agent: PASS with conditions, 0 Blocker, 0 Major, N1-N6; Biomni: 1 residual Major) is folded into one PR: E8 reports a model-free bracket beside every not-estimable interval (Q139 a), the LOD range check is made in log space (Biomni Major), exact-dataset libraries without truth rows are refused and a `run_provenance_binding` blocker names V3 and row authenticity (N2), and the G3 signing procedure is written below under the SCI-05 row, "Round 15". Next: user review and merge, then the user follows that procedure to sign G3 (nothing is frozen or signed here; status `draft`, execution flags false).

**2026-10-14 (SCI-05 round 14 fixes, `codex/sci05-r14-fixes`, `DEV-031`):** round 14 (independent agent: 6 Majors, 13 Minors; Biomni: 5 Majors) returned `does-not-pass`; this branch fixes every Major and every real Minor in one PR (D6 not applicable, recall over all planted truth, `score_production` bound to the verified protocol and an EmptyDrops run config, `--production` guards and the `real_background_anchor` blocker, `frozen_inputs` downgrade forbidden, stratified `f1_interval`, BCa/LOD failure policy, probit quasi-separation, S population). Declared limitations are listed under the SCI-05 row, "Round 14". Next: user review of the PR, then the round-15 re-review; G3 stays unsigned, top-level status `draft`, execution flags false.

**2026-10-14 (LOCK-01 prerequisites, `codex/lock01-prereqs`, Q123-Q133):** `LOCK01-PRE` lands the package side of U1: the adjudicator now computes the D1 locus count and the D2 positive-cell count as LOCK-01 defines them (new `viral_loci.py`, two appended `multimap_evidence.tsv` columns, `results/viral_loci.tsv`), and the study repo marks consumable VAL libraries (`consumable` column, U6). Next: user review and merge; the merge SHA is then pinned in LOCK-01 section 11 (calibration may not start before that). The study-repo rows and Q123-Q133 are recorded below.

**2026-10-13 (SCI-05 round 13, protocol side, `codex/sci05-r13-protocol`, Q118-Q122):** `DEV-030` answers the round-13 protocol findings (agent M6, M7, M8, protocol-text Minors; Biomni BM1-BM4 text and Minors) with the user's decisions: G3 freezes the `probable` tier only (strong, AUPRC and the M8 confidence out of the freeze; the three tier-cutoff measures are defined by Q135/Q136, but their code is pending, so `tier_cutoff_measures` still blocks training), real-data validation is descriptive only (`E11_real_data_descriptive`), unknown-origin calls on real backgrounds are never false positives, and the generator crosses abundance by homology per donor library. The code findings (M1-M5, BM5, BM6) are on `codex/sci05-r13-code`. Next: round 14 once both branches merge; build the code for the Q135 measures (per-cell tier input, features per sample, called-cell support fraction via an evidence run per sample and virus) before any training run.
**2026-10-13 (SCI-05 round 13 code side, `codex/sci05-r13-code`):** round 13 (independent agent plus Biomni) returned `does-not-pass`; this branch fixes the code-side Majors M1-M5 and Biomni BM1, BM3-BM6 plus the code Minors listed under the `SCI-05` row "Round 13". M6, M7, M8 and the protocol wording (BM2, BM1/BM3 YAML) are on `codex/sci05-r13-protocol`. **Next up:** merge both branches (renumber the ledger record `DEV-029` if the protocol branch also appends one), re-run round 14, then the user signs G3.

**2026-10-14 (planter strand fix, `codex/planter-strand`, Q117):** planted reads are now sense to the indexed transcript (and antisense for `--strand reverse` libraries); round 2b `SENS-07` capture is understated. Next: user review of the PR, then round 2c on a fresh pin after round 2b finishes (not submitted).

**2026-10-13 (probable/strong tier code, `codex/evidence-tiers`, Q96/Q104):** `TIER-01` builds the tier code to the frozen calibration-grid contract and adds the one level-to-tier mapping column to `adjudicate_candidates.py`. Next: user review of the contract gaps listed under `TIER-01` (a `strong` rule, the support-fraction definition) before `SCI-03` freezes calibration; then `REF-09`.

**2026-10-13 (Q83 choices pasted, `codex/protocol-fill-uc`):** the user accepted all 25 USER CHOICE recommendations; the pending-phase blocks are pasted, re-based onto the post-#46 protocol, as ledger record `DEV-028` (W2 now covers all 13 datasets, 98 rows; DSR-13 retry rule; comparator version policy; independent-review pass bar and sidecar path; and the real-data rule: Q103/Q106 real-data arms report agreement with presumed labels descriptively, accuracy claims only from synthetic truth). **No G3 flip was done.** **Next up:** SCI-05 round 13 (independent agent plus Biomni), then the user signs G3; before the freeze the tier rule must expose the three grid knobs (UC-C1) and the multimap rerun cost should be measured (UC-D2). See the `SCI-05-Q83` row.

**2026-10-13 (SCI-05 round 12, protocol side, `codex/sci05-r12-protocol`):** `DEV-027` lands M1, M5, M6 (partial), the V03 context row (`T2-AMEND-01`), m1-m4 and the structured B15 fields; see the `SCI-05` row for per-finding status. User decisions 2026-10-13: Q86 (exposure on context-only data is not an outcome; no flag flipped), Q85 as released (the grill-2 HNSCC/tonsil hold-out is withdrawn; their runs are inventoried as `PFR-13`), Q87 (scored-T2 class deferred; V03 stays a context row). Still open: CAL-11 gating in the study repo. The code-side M2/M3/M4/m5 are PR #43.

**2026-10-13 (multimap CLI, branch `codex/multimap-3-options`):** `Q71-Q76` recorded under `CLI-17`; `EMC-01` re-measure passed, so `em` = `em-cell`. Next: user review of the PR; after G5 the default follows the frozen protocol rule (Q75).

**2026-10-13 (`SCI-05` round 12, branch `codex/sci05-r12-code-honouring`):** round 12 returned `does-not-pass`; the code-honouring findings M2/M3/M4/m5 are fixed here (see the `SCI-05` row). Next: the protocol-side findings (B1, M1, M5, M6, M7) on a branch that may edit `protocol.yaml`/`deviations.yaml`/schemas, including the structured `empty_barcodes_per_sample` field M4 still needs.

**2026-10-13 (Claude continuation):** validate the combined package branches in an isolated checkout. `SENS-07` scorer-group and exact-sequence receipt corrections landed via `codex/land-resume-fixes` (with the host-response docs/pin refresh, `REL-23` lock-audit workflow and `REL-24` maintainer checklist from the resume-qualification branch); real per-dataset measurement stays pending. Round-2 tooling and discovery registries resume from their saved work. The approved publication steps await the current sandbox approval; scientific gates remain unchanged.

**2026-10-13 (grill session 4, user-confirmed Q53-Q70; branch `codex/trackb-gates-q53-q70`):** *Process:* agents run downloads and installs under a GUARD; the user runs the `PANEL-01` WP5 build, owns the `SCI-05` sign-off and branch protection; an agent places the `DEF-07` notice (Q53, Q70). Before G3 only data prep (downloads, `PRE-01`, pilots, T2T index); no quantification until G3 and WP5 (Q54). *Data:* calibration sets now (V01 local, V02 `GSE226438` 3+3, V03 if raw GEX exists, V04 5 CHB + 5 FC, `CAL-08` backgrounds, N01) plus a 1-2 run pilot per discovery/replication set; bulk streamed after `LOCK-02` (Q55); layout `viralscan_work/discovery/{fastq,runs,checksums}`, 90 % disk cap (Q67); T2T-CHM13v2.0 `GCF_009914755.1` with a minimap2 `sr` index, HPRC for L1+ candidates only (Q68). *Pre-registration:* locus, positive-cell and sibling-bleed definitions, the rare-super-expressor route, `CAL-13` strata and decoy ratio go into `LOCK-01`, set from calibration data only (Q56); `CAL-13` decoys run once in calibration (V positives + 5-10 real negatives; Q57); the lock cites the sha256 of the protocol file bytes (Q59); the lock schema is checkout-only (Q60); lane = `flowcell:lane` (`PRE-01`, study repo, Codex; Q61). *Tools:* `ADJ-01` L3 reads the nine herpesviruses through a taxid map (Q58, done in this branch); `LONG-01` defaults accepted: `min_detected_level` L1, pooled libraries excluded from trajectories but counted for contamination, depth table from registry + run outputs (Q62); round 2c plants HHV-6B `SRR20710641` with R1 trimmed to barcode + UMI, `GSM4682311` stays n/a (Q63); per-family seeds from round 2c on (Q64). *Reference and review:* WP5 builds now on the current shipped tier, WP3 gates as acceptance, remaining WP2 promotions go to protocol v2 (Q65); `SCI-05` round 12 = independent agent review + Biomni external review, user signs G3 (Q66); `VAL` truth panels are generated by an agent after WP5 against the frozen reference (Q69). **Next up:** the user's WP5 build (command under `PANEL-01` WP5), then `LOCK-01` from calibration data.

**2026-10-12 (grill session 3, user-confirmed Q19-Q40):** two public-dataset tracks. *Track A, exploratory round 2 (`DSR-16`, new):* full rerun of every arm on a new pinned `main` commit, `cat42d` reference, the 22 frozen round-1 chemistry JSONs, arms `combined_off`, `combined_artefact`, `twostep` (new DEF-05 cell calling; compare with round-1 `twostep_v2`); public development set plus COVID x213/x216 and SFL x223 reported separately; `ERR13027122` (HTO/ADT) and `GSM4682312` (no barcode read) dropped as n/a; agents submit under `dsr_round2/GUARD.md` (resources only), throttle `%4` combined / `%2` two-step; `SENS-07` per-dataset planting as round 2b on the same commit. Every run logged `explore`. *Track B, broad virus discovery study (`DISC-01`, new):* lives in a local private study repo next to the package (seeded from the untracked Biomni plan bundle, which is never committed here); generic tools land in ViralScan as package features. Shipped-tier reference after the user's `PANEL-01` WP5 build; `DEC-02`..`DEC-04` confirmed (hash-pinned panel, combined host+virus, fixed `unstranded`); calibration, then `LOCK-02`, then blind discovery, then replication; no discovery run before the lock. Dataset facts (NCBI/ENA metadata, 2026-10-12): all series open access except `GSE317012` (raw data never deposited: replacement positive to be searched, fallback a recorded two-positive deviation); about 3.8 TB GEX FASTQ; `GSE151671` is 10x 3' v2 with R1/R2 as separate runs (pair by sample; no Drop-seq gate); `GSE190338` is BAM-only (`bamtofastq`); most discovery sets are 10x 5'. Influenza `GSE336913` also joins the package's scored T2 rows via an amendment in the SCI-05 review. Storage: calibration FASTQs kept until the lock, discovery/replication streamed (download, run, keep outputs + checksums, delete FASTQ). Host-homology screen: T2T-CHM13 now, HPRC only for L1+ candidates (user downloads). Healthy 3' and 5' PBMC backgrounds outside the package holdout, picked by an agent and approved by the user.

- [x] `PERF-B13` (2026-10-19) — bounded invocation-local raw-sequence classifier
  cache accepted: 256 entries, maximum 512 key characters, long reads bypass it.
  Raw case and barcode/UMI span retained; validation/serial output order unchanged.
  47 focused tests pass, including exact gzip bytes and late failures. Two seven-pair
  complete-filter sessions: duplicate wall 0.0527/0.0525x, reagent 0.1356/0.1331x;
  unique 1.0232/1.0173x, long-read fallback 0.9944/0.9944x. Synthetic stage benefits;
  no production whole-run claim. Independent baseline d056986; combined local/native qualification passes.

- [x] `PERF-B28` (2026-10-19) — parse captured immutable GenBank text once for
  CDS/noncoding groups, retaining CDS-first emission, ordinals and all independent
  integrity reads. 96 focused tests pass (one network test deselected); full offline
  stale-GTF preparation passes two seven-pair benefit sessions (0.8964/0.8863x wall).
  Warm-cache control 1.0243/0.9697x. Raw output/cache tree signatures match; no network
  requests. Mutable-path digest reuse rejected by a sidecar-mutation counterexample.
  Combined local/native qualification passes; no remote-service or production whole-run claim.

- [x] `PERF-ASSESS-04` (2026-10-19) — remaining cards have explicit dispositions:
  B21/B13/B28/B16 accepted subject to final combined C00. B01 weighted-buffer
  prototype wall1.0209x/peak0.9241x and B02 index reuse wall0.9308x/peak1.0012x
  fail benefit gates. Real R4_51/DropletUtils A/A passes; B12 adapter not material
  on the bounded fixture. B22 repeated preparation0.696% of profiled real fits.
  B03/B09/B14/B15/B17/B26/B29 retain observable rereads/failures/validation; B19 has
  no repeated BLAST DB workload in current orchestration. B28 optimizes immutable
  parsing only. G01 retrospective mutation receipt passes with historical goldens
  unchanged. B11/B30: total-budget policy approved and implemented 2026-10-19 on `codex/cpu-budget` (`--cores` bounds kb/minimap2/CD-HIT/rule threads; `build-ref --cores` new, unset = legacy; numerical libs pinned at 1; manifest worker fields still open). 2026-10-20 review-fix requalification (local, no SLURM): plumbing tests added for every CLI-to-native call site (`tests/test_legacy_worker_requests.py`: menu `main`, `build_ref_main` anellovirus and combined paths, kb/minimap2/cd-hit-est argv all carry `--cores`); DAG test now asserts rule expressions with a scheduler above the budget (`--cores 64`) so a literal `threads: 8` fails; multimap BLAS pin checked at runtime; deleting the kb_count OMP/OPENBLAS export fails 12 of 16 `test_count_legacy_request_and_library_environment` cases. Mutation checks: every site and the literal-threads mutation fail a test. The overnight tiny-fixture native matrix matched all scientific bytes, but an earlier retained one-core finish-reference index differs by 55 bytes from its raw-identical baseline A/A. Forty subsequent direct native builds (20 each at one/eight workers) matched the baseline; they do not explain or supersede the earlier mismatch. Four downstream count runs using the retained differing indexes match scientific files on the two-read fixture. Logs/manifests and the harness's own `config.yaml` also differ as recorded. No comparator relaxed; raw-index compatibility remains unresolved and merge stays blocked pending that decision. See docs/plans/2026-10-20-cpu-budget-qualification-status.md. PR77 CI exposed an unrelated Rscript dependency in the CLI plumbing fixture; select knee calling explicitly and assert the expected successful exit. The Rscript-absent reproduction fails before this test-only correction; the affected CLI/worker suites pass 184 tests (16 integration cases deselected) with Rscript hidden after the fix. Ruff and formatting pass; production code and defaults are unchanged. Full grids stopped
  after early rejection, not claimed complete. See the qualification/disposition
  report and resource proposal under docs/plans/2026-10-19-performance-*.

- [x] `PERF-RESOURCE-01` (2026-10-19, user explicitly approved local qualification)
  — qualified explicit legacy worker requests only: kb count/ref8 and homology
  minimap2 3; unchanged effective workers, Snakemake thread declarations, CD-HIT
  -T0/-M8000, search parameters and numerical-library environment. Independent
  baseline9f0c397; native A/A before candidate edits, then A/B raw-artifact and
  error/partial-output checks plus 1/6/8/32-core DAG controls. No scheduler jobs,
  new data, scientific changes or total-admission policy implementation. Resource
  provenance differences must remain truthful and explicit; no performance claim.
  Native six-stage A/A/A-B and supplemental count artifacts pass; 237+48 default
  regressions,34 scheduler/failure tests and32 paired DAG cases (64runs) pass.
  PR76 initial CI found stale build-reference evidence hashes; refresh the
  claim/inventory pins against immutable source commits without changing runtime
  or claim calibration status. The registry inventory pin follows the immutable
  registry commit1b731be; preserve its non-release provenance tag through rebase.
  Final local native count/strand/lineage/end-to-end54passed with no skips;
  governance/worker27passed. PR76 final CI remains the merge gate.
  Full B11/B30 total-admission policy remains open.

- [x] `NORM-SHIFT-01` (user note, 2026-10-19) — use shifted-log normalization
  from [Single-cell best practices](https://www.sc-best-practices.org/preprocessing-visualization/normalization/#shifted-logarithm):
  `normalize_total(target_sum=None)` (median library depth), followed by `log1p`;
  do not use a fixed TP10K target. Four current sites: `umap.py:host_only_pca`,
  both main UMAP preprocessing paths, and `hostresponse.py:_detect_and_normalize`.
  Keep raw counts and normalization provenance. This changes scientific features
  and requires explicit downstream validation; no normalization edits are in PR76.
  Implementation started after user "next": isolated branch from6bd1652, four
  fixed targets replaced with Scanpy median-depth scaling and log1p. Preserve
  raw input files/layers, raw-depth covariates, preprocessed-input handling,
  detection/counting and thresholds. Validate dense/sparse/zero-count behavior,
  both UMAP branches, host-only PCA and fitted host-response outputs; this is
  an intentional scientific preprocessing change, not a byte-identical speedup.
  Local234regressions pass (including13new normalization cases); historical
  weight golden retained on explicitly preprocessed legacy input, and no-HVG
  fixtures made constant-expression without changing thresholds/assertions.
  Paired original A/A and candidate A/A fit artifacts are raw deterministic;
  expected A/B differences include synthetic AUC0.5→0.375 and stable genes8→5.
  This is not evidence of predictive improvement. Scanpy1.11.5 has different
  dense/CSR default medians with empty rows; exact target_sum=None semantics
  are retained and version/procedure provenance recorded. Native CLI and final
  full PR CI remain merge gates. No production inputs or previous outputs altered.
  Documented native CLI21passed/no skips, including reference build, visual
  run, evidence artifacts and validate-run. Governance21passed.
  Public output documentation is re-pinned to implementation commitdc7855f;
  preserve non-release tag provenance/shifted-log-2026-10-19 through rebase.
  No claim status, historical golden values or model/QC thresholds changed.

- [~] `PERF-EXEC-04` (2026-10-19) — authorized performance code delivery:
  B21/B13/B28/B16 accepted and combined on frozen baseline d056986 (tree6e70eda).
  Each has independent exactness and two seven-pair benefit sessions. Combined:
  154 new oracle tests,53 real offline native integration tests,12/12 retained
  reducer raw files against pre-performance main6910c2d; bounded fitted cell/group/
  strata outputs27/30raw equal, three precise source-provenance manifests. Source
  freeze326files; original oracle function bodies match baseline raw source and AST.
  Ruff/format275files, mypy60 and data governance pass. G01 historical mutation
  receipt recovered without golden edits; all17 legacy component timing rows now
  repeated in two7pair sessions.
  PR #75 merged as9f0c397 after all12CI checks passed; code delivery complete.
  Total13accepted/15deferred cards; rejected edit grids stay unchecked. G02 production
  attribution/whole-run measurement and B11/B30 total-admission policy retain
  separate authorization gates; local legacy qualification is approved above, with concrete proposals in docs/plans/2026-10-19-
  performance-{whole-run-proposal,resource-policy}.md. No scheduler jobs, downloads,
  scientific/worker changes or retained-data deletion. Qualification report and raw
  receipts retain all failed runner/prototype attempts. This row remains partial
  until the gated plan work is resolved; code delivery is a separate completion claim.
  **Deferred 2026-10-20:** G02 production attribution/whole-run measurement and the B11/B30 CPU-budget follow-ups (manifest worker fields, total-admission policy) are closed as deferred: the user parked performance work behind Track A on 2026-10-19, and the Slurm measurement needs approval before any job is submitted. Row stays `[~]` for the delivered code.

- [~] `PERF-EXEC-03` (2026-10-18, original branch `codex/perf-wave2-b25-b27-b18`,
  historical baseline `faf5667`; finalized on `codex/perf-wave2-final` over merged main
  `5e05a27` after PR #73): wave-2 B08/B10 gate decisions (B08 deferred, B10 no compliant fast path — no
  source edits), B27 negative (whitelist I/O dominates — no source edit), B25 (bounded
  reference-input concatenation) and B18 (fused molecule reducer + streaming writer)
  implemented byte-identical, all qualified together as C00-wave2. Combined vs a `faf5667`
  git checkout: 147/148 harness artifacts byte-identical, the 148th differing only in its
  own `build_receipt.software_identity.git_commit` (equal after scrubbing that field), tip
  A/A identical; B25/B18 card harnesses re-run identical at the tip; 5,028 tests pass (19
  skipped, 138 deselected), governance pins consistent on both arms (no re-pin needed),
  ruff/mypy clean. Component benefits: B25 peak flat at 2.01 MB vs 7.96x/31.84x/63.68x
  higher at 8/32/64 MiB inputs and 0.71-0.83x wall (0.49x at 64 MiB); B18 traced peak
  0.728x/0.641x/0.571x/0.526x at 1/2/4/8 reads per molecule (0.468x at 1M reads), wall 0.925x
  median over 7 pairs, 7/7 faster. Both timed in one loud window (load 11.5-19.1); no
  whole-run claim. Not done: a quiet-host second timing session, whole-run measurement,
  wave-3 B29/B17 and the remaining cards. Landing authorized by the user's merge/finalize
  request. Reconciled source preserves main's newer deduplication helpers and both
  unchanged frozen oracles; focused tests 97 passed, one sdust skip; lint/format/mypy
  pass. Final full-suite and CI results are recorded on the finalization PR. Historical evidence:
  `analysis/computational_bottlenecks_2026-10-17/followups/20261018-112322-wave2-g00/c00/C00_WAVE2_REPORT.md`
  and the branch's `docs/plans/2026-10-17-remaining-performance-candidates.md` card rows.

**2026-10-12 (grill session 2, user-confirmed Q1-Q18; branch `codex/grill2-decisions-2026-10-12`):** *Datasets:* exploratory runs (labelled `explore`, logged in an explore ledger) on the development set HSV-1, EBV, HHV-6B, COVID x213/x216, KSHV `GSE190558`, HPV16 `GSE189670` and BC-1 `GSE154900` (context-only rows are exempt from the partition rules); HNSCC `GSE164690` and the tonsil-atlas negative subset `E-MTAB-13687` are held out untouched until G3. Only the HPV16 raft/parental pair and iSLK.219 latent/induced become scored T2 rows; EBV, HHV-6B, HSV-1, BC-1 stay descriptive. *Build:* max EC size 9,921 is the `CAT-40` baseline; no `kb count --mm`; `pbmc_10x_healthy_v3` is the D1 background and T4 negative (user downloads and checksum-pins it; label 'presumed negative'). *Code:* `B-19` zero-fill host-called cells missing from the viral matrix, drop and log viral barcodes not host-called, error if no host cell set; `C-12` exclude `cores` from the run fingerprint (record it as information) and name differing options in the error; cell-type labels only from author annotations, never from ViralScan output. *Content:* DEF-07 notice approved as drafted (beside, never over); `CAT-10` phased, v1 = 9 herpesviruses + HIV-1/2 + HTLV-1/2 + HBV with `verify-references` and user review before commit. *Claims:* no wet controls required for the first release (affected claims stay 'screening support'); no numeric specificity pass/fail: report the upper bound with n; >= 3 donor-distinct libraries for any T2/T3 sensitivity claim; >= 10 independent T4 libraries before 'specificity assessed'; LOD95 per genus reported, not pass/fail. SCI-05 review round 12 starts when the software queue has merged; G3 freezes right after the user's sign-off. *Process:* branch protection on `main` (user, GitHub settings); weekly non-blocking lock-audit workflow; prune merged worktrees/branches (dirty old worktrees saved as patches under `.living/archive/` first); external labs install from a pinned commit, no `v*` tags before the RC gates pass; an agent drafts `MAINTAINER_CHECKLIST.md` (user fills it; nothing published or reserved).

**2026-10-12 (grill session, user-confirmed; branch `codex/grill-decisions-2026-10-12`):** goal = finish all software so exploratory runs on public datasets can start. Exploratory (unscored, labelled `explore`) runs come first; scored G5b runs wait for the G3 freeze. Decisions are recorded under `AUD-12` (M-03/M-04), `CAT-21`, `CAT-29`, `DEF-04`, `DEF-07`, `ANDET-04`, `ANDET-07` and in `.living/decisions.md`. Software queue, in order: batch 1 `M03-IMPL`, `M04-GEN`, `CAT-29`+`DEF-07` (parallel worktree agents, 3 max); batch 2 `M04-SCORE`, `HR-01..05`, DSR/verdict tooling (`DSR-02`, `DOSSIER-01..03`, `VERDICT-01`, `CARD-01`, `SENS-07`, `DLIST-02`), `REL-04/05`, AUD-13 minors, ANDET-01/02 follow-ups, `CAT-40`. **Deferred, not dropped:** protocol decisions (which context-only positives become scored rows; acceptance thresholds), DEF-00/SCI-05 sign-off and the G3 freeze, wet controls (same-pool negative, reagent control, orthogonally confirmed positive). First exploratory batch once software is done: HSV-1 `SRR8315713`, EBV `SRR12682296`, HHV-6B `SRR20710641`, COVID x213/x216 (standard arm set each), then KSHV `GSE190558` and HPV16 `GSE189670`. PRs: one `codex/*` branch per task, green CI, `PLAN.md` in the same commit; the user merges (`! gh pr merge <n> --merge`).

**2026-10-12 (branch `codex/ci-green-2026-10-12`, follow-up to PR #13):** the four red-CI groups are fixed and awaiting a CI run on the PR. (1) `DEV-006` (user decision: archive ref, ledger untouched): tag `archive/dev-006` on `291ff14` pushed, so the commit resolves under a full-history checkout; it still sits on no branch, and `check_ledger_append_only.py` only needs it to resolve. (2) `compare_legacy_v2_v3_attempts.py`: anndata >= 0.13 yields `None` (the X matrix) when iterating `adata.layers`, which broke `sorted()` and was swallowed as `attempt_contract`; `None` is now skipped (reproduced on CI's anndata 0.13.4 / pandas 3.0.6, py3.12 only). (3) Python 3.9: `zip(strict=...)` removed from four scripts, five call sites, one multi-line and found only by an AST scan (the length invariant in `prepare_fresh_controls.py` is now an assert); `gen_cli_reference.py` strips the `(default: X)` suffix that `BooleanOptionalAction` adds to help before 3.11, so the doc is identical on 3.9 to 3.12. (4) macOS: `test_mate_mismatch_fails` wrote its inputs to `r1/r2.fastq.gz`, which are `R1/R2.fastq.gz` on a case-insensitive filesystem, so the filter truncated its own input; inputs renamed. `pyproject.toml` build requirement raised to `setuptools>=77` (SPDX `license` string) and the pytest job installs `setuptools>=77`, fixing the missing `setuptools.config.pyprojecttoml` on the macOS 3.9 runner. Local: 3,388 pass, ruff/mypy/governance/ledger clean.

**2026-10-11 (Codex parallel backlog resume, branch `codex/audit-fixes-2026-10`, baseline `8d03f17`):** recover the existing uncommitted `codex/aud-15b`, `codex/aud-21b` and `codex/aud-13` work in separate worktrees; finish the reference, detection and minor-finding software before integrating. The `aud-m2` merge is already complete. `AUD-17/18/22` are merged; ANDET-01 breadth and ANDET-02 per-locus coordinates are implemented. ANDET-04 already has an opt-in rule; default-on remains conditional on calibration and a user decision. Continue agent-doable packaging, documentation, release and validation-tooling backlog in subsequent parallel batches. Refresh governance pins after the final source/doc edits and verify one combined checkout. Scientific outcomes remain behind the existing package/preregistration gates; external metadata, DOI registration, publishing and infrastructure configuration cannot be inferred from software checks.

**2026-10-10 (wave 4, branch `codex/wave4-package-a`, plan in `~/.claude/plans/thoroughly-plan-to-finish-mellow-starlight.md`):**
the agent-doable list-A rows are done or advanced. **Closed:** DSR-08, SW-07 (boundary count = `mixed_host_virus_manifest.tsv` read_ids),
CAT-06 (viral IDs from the build manifest, catalogue sha256 pinned), PROG-12, PROG-13, CAT-28 (report only). **Advanced `[~]`:** ANDET-01,
ANDET-02 (orientation defect fixed, demotion waits on REF-08), ANDET-04 (opt-in `--auto-evidence` rule), REF-07/03, REF-04 (compare tool),
DEF-08/REL-03 (explicit lock `environment.explicit.txt`), REF-13 (build-time normaliser `gtf_normalise.py`; integrity 2,626 raw -> 7 normalised).
**Left for the user, in order:** (1) `gene_programs.tsv` is stale against the normalised index (63 of 75 rows' overlap groups differ; EBV
independent groups 5 -> 2). The generator's HHV-6A name is fixed (2026-10-11), so a full regeneration now changes only `overlap_group`
and `non_overlapping`; but it puts EBV BMRF1 (productive) in one group with the latent EBNA genes (a 37-bp EBNA-1.2 exon at 50134-50170
overlaps BPLF1/BOLF1 and transitive closure chains to BMRF1), which fails `test_cited_biology_corrections_2026_10_02`. Shipped TSV left
unchanged pending `AUD-01` (k-mer-compatibility grouping; Biomni review rejected the min-overlap option); (2) `PANEL-01` WP5 cluster build and freeze, then
regenerate `gtf_corpus_manifest.tsv`; (3) re-run `kb count` for runs made before the PROG-13 t2g merge; (4) re-measure `em-cell` memory (EMC-01);
(5) replay `environment.explicit.txt` with `conda create --file` and the CI runner result of `snakemake9-dryrun`; (6) GitHub Action SHA pins,
docker, `REF-11`, `REL-15`; (7) design calls: `DEF-04`, `reactivating` gate beyond HCMV, `BBLF2` alias, the unwired `build-ref` t2g merge and normaliser.

**2026-10-10 (wave 3, branch `codex/wave3-package-a`, plan in `~/.claude/plans/thoroughly-plan-for-all-swift-elephant.md`):**
**done 2026-10-10 except what is marked `[~]`/`[!]`:** full suite 2,878 pass, 70 integration tests pass (real kb/snakemake tiny run),
mypy (48 files), `ruff check .`, `ruff format --check .` and the data-governance check are clean. **Left for the user, in order:**
(1) ~~decide the `NC_002076.2` GTF~~ done 2026-10-10: `Alphatorquevirus homin1.gtf` kept (REF-13); (2) `PANEL-01` WP5 cluster build and freeze; (3) re-measure
`em-cell` peak memory on the 793 k-barcode run before lowering `--mem` (EMC-01); (4) GitHub Action SHA pins, docker build, `REF-11`, `REL-15`;
(5) design call on `DEF-04` gene-level `whole_genome` exclusion; (6) VAL generator rows stay behind G3. Wave 3 rows, one commit each:
- [x] W3-01 tracker reconciliation (DEF-05, CHEM-01, CAT-18 closed; AB303562.1 rows, GTF-manifest check noted done).
- [x] W3-02 `--n-seeds` / `hostresponse_n_seeds` outside 1..6 is rejected (was silently truncated; 0 became 6).
- [x] W3-03 `pct_infected_comparable` blank, not 0.0, when no cell is comparable (kallisto and alignment-only rows).
- [x] W3-04 DEF-05 hardening (shared host-matrix validator before R, exact `matrix.mtx` in `emptydrops.R`, one `COMPARABLE_CELL_MIN_UMI`; Rscript+DropletUtils exist in the conda env, so the R path is tested for real).
- [x] W3-05 DEF-06 property tests (`tests/test_multimapping_properties.py`, 388 cases; claim/inventory pins added in the close-out governance commit).
- [x] W3-06 EMC-01: em-cell keeps each cell's fitted abundance only at its compatible genes (identical results, dense-reference test; tracemalloc test fails on the old code).
- [x] W3-07 CAT-38/39 done; REF-13 offline: astrovirus feature type fixed; `[!]` duplicate `NC_002076.2` (two GTFs, different gene-ID conventions: which one the builder consumes is a user call), exon-less CDS normaliser, 6 genes with neither CDS nor exon, `curated` list; DEF-04 deferred (see row).
- [x] W3-08 VAL-06/07 core scorer + shared partition (code only, no outcome runs; both rows `[~]`, see them).
- [x] W3-09 SW-03 (optional `results/molecule_assignments.tsv.gz`, `--multimap-molecule-assignments`; weights reconstruct `corrected` per cell/gene on 20 random seeds x 6 methods; `rerun-multimap` forces a full pass when the file exists).
- [x] W3-10 MECH-C (`menu._write_run_config` is the one writer; Snakemake reads `--configfile`; the `create_config` rule, `createconfig.py`, `_build_config_args` and `to_snakemake_config_args` are deleted; 396 focused + 70 integration tests pass incl. the tiny end-to-end run).
- [~] W3-11 SW-09: the 13 parser builders moved verbatim to `src/viralscan/cli/parsers.py` (menu.py 2,530 -> 1,196 lines; `gen_cli_reference.py --check` byte-identical; `menu` re-exports every name). Handlers stay in `menu.py`: they call names the tests patch on `menu` (`_check_required_tools`, `_swap_multimap_layer`, ...), so moving them would break those patches for no gain.

**2026-10-10 (review of wave 2, pushed):** the counting contract, two-step cell calling and host-response work
were reviewed; three defects fixed with regression tests: `run_manifest` schema lacked `sibling-weighted`
(`validate-run` would have rejected such runs); `cell_calling_input.json` recorded no input digests
(now SHA-256 of the GeneFull matrix, called-cells file and `emptydrops_cells.tsv`); host-response reported a
fold-subset depth baseline as all-fold and masked per-stratum failures in `viruses_status`. Uncertain, not
changed: `pct_infected_comparable` is 0.0 (not unavailable) when no cell is comparable; HTML "% Infected" is an
all-barcode figure in two-step runs; host-response differential p-values are cell-level with no donor
adjustment and `--n-seeds` above 6 is silently truncated. `uv.lock` is now tracked and consistent (`uv lock
--check` passes; snakemake>=9 is gated on python>=3.11). Docker build is not possible on this host (no docker).

**2026-10-10 (wave 2 finished and committed):** full suite 2,265 pass; ruff check and
mypy (47 files) clean; the only failures are the two governance pins, re-pinned last. Closed: `SENS-CORR-01`,
`EVID-COMP-01..05`, run-provenance (`software_identity` stays outside the resume fingerprint), the panel
count test (4 stale `AB303562.1` rows removed: 2,511 rows / 1,994 accessions), and three silent-pass paths in
`panel_integrity.py` / `prepare_reference_inputs`. Counting contract, two-step cell calling and the recovered
host-response work are covered by their new tests (passing) but were not re-reviewed this session.
Left to the user: action SHA pins, `uv.lock`, docker build,
NCBI fetches, EC-size gate on a real index, `REF-11`, `PANEL-01` WP5.

**2026-10-09 (resumed after Claude session limit):** wave 1 is integrated through `654bca0`:
`PKG-01`, typing, `SENS-CORR-03` and `EVID-CORR-01`.
Installs, push, CI, `REL-15`, `REF-11` and the `PANEL-01` WP5 build stay the user's.

**2026-10-09 (package continuation):** `PANEL-01` WP2 role plumbing and `SW-26` scratch isolation
are implemented and verified. Catalogue labels reach the TSV/HTML report without changing counts;
native reference builds preserve the checkout's existing `tmp/`. Source commits: `e665003`, `dc989a0`.
Verification: 1,929 unit passes; the sole stale registry-pin failure was repaired, then all 21
inventory/claim checks passed. CLI/reference checks: 260 passes; affected native integration checks:
25 passes. Changed-file lint/format and data governance pass; two pre-existing lint findings remain
in `test_chemistry_check.py` and `test_snakefile_dag.py`. No scientific claim status changed.
Next package work: remaining `PANEL-01` promotions (SV40/AAV1/7/8 need NCBI records), then WP3 gates;
WP4 GTF packaging remains deferred and WP5 build/freeze remains user-owned.

**2026-10-07 (latest, dossier):** inventory first (`DOSSIER-01..03`): run
`python scripts/dsr05_dossier.py <round>` and review `<round>/dossier/round_index.md` with the user before any
regeneration, barcode-rank plot or planted-read job. Evidence regeneration on `vs_pinned/14ca266`
(`evidence_v2/calls.tsv`, 140 calls) is enumerated but not submitted.

**2026-10-07:** TSO handling landed (`TSO-01`, uncommitted). DSR round 1: redetect twostep rows 5-7,
14-16 running (25725712); x223 row 18 needs the sample-name fix in `slurm_dsr_twostep_redetect.sh`
(manifest `x223` vs directory `LUM-SJ-x223`). Then rebuild `common_cells_summary.tsv` with `--reference-cells`
and `--verdicts`, re-enumerate calls, second evidence batch, commit.

**2026-10-06 (latest): DSR round 1 is running** (`DSR-13`). Every downloaded
dataset is rerun on one pinned commit (`bbf1821`) with one arm set, driven by
a frozen manifest and `check-chemistry` JSONs; subagents submit and monitor the
arrays. The evidence OOMs were a `samtools depth -a` plot bug, not BLAST (`DSR-14`, fixed).
Next: pin the fix as `vs_pinned/<fixsha>`, rerun the F-025 gate at 64 G (the alignment layer
already reproduces F-025: 14/14 reads telomere repeats), then enumerate calls
(`scripts/dsr02_enumerate_calls.py`), run the evidence array, and record F-027… and the arm matrix.

**2026-10-05 (latest): science track, Phase B.** The VAL-01 design
(`docs/plans/2026-10-05-val01-generator-design.md`) is fully decided: the user
closed D1–D7. Next:
1. ~~evaluate REF-06~~: done 2026-10-05; all gates pass, so cat42d is the
   latest reference;
2. ~~REF-07 host-homology table~~: done 2026-10-05 (F-023);
3. ~~`DEF-00`~~: written into `protocol.yaml` v0.3.0 on 2026-10-05 (`DEV-020`);
   §F decided as `none`/`repeat_homology`. Grid membership was decided on
   2026-10-05; `COST-01` is done (F-024); all 5 grids written (DEV-021,
   DEV-022, DEV-023; `DEF-01` done). D18 reconciled (DEV-025, 2026-10-06).
   Next: check the grid cost against Q10d, then the independent review
   (SCI-05) and the user's G3 sign-off;
4. only then code VAL-01.

**2026-10-04 (latest): package completion before more experiments** (user
directive).
- Plan: `docs/plans/2026-10-04-package-completion-plan.md`. M0 tracker
  reconciliation → M1 G1 green + `3.0.0.devN` → M2 G2 distribution → M3 G6a
  non-numeric docs → M4 schema-changing features. All experiments are
  deferred, including `ANDET-09f` and the census on the counted reads.
- Baseline (HEAD c8bc8c7):
  - unit suite green (1,681);
  - ruff check/format **red** (18 / 26, mostly Markdown and vendored code);
  - integration gate **hangs** (`test_tiny_end_to_end`: REL-16 kallisto
    mismatch), 1 fail (`test_anellovirus_chain` without `dustmasker`), 2 tool
    skips.
- `REL-16` is therefore a G1 blocker.
- **User, 2026-10-04:** focus on the package; the user does network installs,
  push and CI afterwards. D4 was not answered, so the plan's recommendation
  stands (only schema-changing features before rc).
- Progress:
  - [x] M0 tracker reconciliation (69e2975).
  - [x] M1.1 lint. ruff excludes now match the "package + tests" intent
    (`extras`, notes, `*.md`); `ruff check` and `ruff format --check` are
    green.
  - [x] M1.2 `REL-16` (d352927).
  - [~] M1.3 test tool skips: no change needed. The triage showed
    `test_anellovirus_chain` fails on a stale fixture, not a missing
    dustmasker (ANELLO-14). The 22 tool skips wait on the user's full env.
  - [x] M1.4 triage: 7 confirmed defects, **all fixed** (2026-10-04):
    SW-22, ANDET-05, CAT-41, PROG-14, CAT-37, REF-13 residue, ANELLO-14,
    MECH-F off-list count. The two tests that passed only with git-ignored
    `data/*.gtf` are now hermetic (227a2f6). Caveat: the EBER test checks a
    committed fixture, not live `ncbi_fetch` output.
  - [x] M1.5 G1 gate run. **Green 2026-10-05 in the full env**
    `viralscan_test_full` (micromamba; every gate tool present):
    - `pytest -m "integration and not network"`: **73 passed, 0 skipped**;
    - `pytest tests/`: 1,715 passed; `ruff check` and `ruff format --check`: clean;
    - the full env exposed two non-hermetic test files
      (`test_multimap.py`, `test_bus_correction_audit.py`). With `kb` on PATH,
      `tool_path` ran `kb info` through their fake `subprocess.run`. They now
      stub `tool_path`.
    - Earlier partial run (2026-10-04, HEAD d627924):
    - `pytest tests/`: 1,704 passed;
    - on a clean checkout: 1,697 passed, 7 skipped (need git-ignored data);
    - `ruff check` and `ruff format --check`: clean;
    - `pytest -m "integration and not network"`: 50 passed, 23 skipped, no
      hang.
    - `DEF-09` followed (local tag `v3.0.0.dev1`).
  - Also landed in parallel:
    - M4 schema items ANDET-03 (`claim_scope`) and ANDET-01 (index-gene
      breadth);
    - M2.1 (REL-01/02, DOC-04);
    - Wave 2: the evidence-route measures; CAT-25 pins, CAT-13 size gate,
      the DEF-08 CI dry-run and REL-03 prep; DOC-09 coverage. HEAD is green:
      1,710 passed.
    - M3 doc fixes (ANDET-06, OPS-01/02 drafts; DOC-07 stale note).
    - M4 `MECH-B` (2026-10-04): detection thresholds each virus's summed
      count; `sensitivity.tsv` gets zero rows for undetected indexed viruses.
  - **Pre-rc package work is complete** (M4 recommended list done; `CAT-18`
    waits until after rc). G1 is green and `DEF-09` is done (local tag
    `v3.0.0.dev1`, 2026-10-05). Next is user-side: push (the tag triggers the
    PyPI/ghcr release), lock resolution, Sphinx and CI.

**2026-10-04 (night): plan for the remaining anellovirus work** (superseded by the
package-completion plan above; its experiment steps are deferred).
- The TSO-barcode audit closed at 0 of 57,715.
- **Caveat:** the 81/62 census bound covers the July evidence set (raw FASTQ,
  uncorrected BUS), not the counted molecules, so do not quote it against
  57,715 UMI.
- Critical path:
  1. `REL-16`;
  2. census on the counted reads;
  3. the measures in `viralscan evidence`.
- In parallel on Slurm: `ANDET-09f` (regenerate the plant, then the sweep).

**2026-10-04 (evening): reagent leak fixed.** `reagent_fraction` (junction-context
TSO + forward TruSeq) and `r1_tso_fraction` replace `tso_fraction`. On the
covid census, kept reads fall 6,171 → 790 with sensitivity unchanged. Next: the
TSO-barcode audit under `ANELLO-PRIOR.3`, then `REL-16` and `ANDET-09f`.

**2026-10-04 (later): body census done, and the read-side measures leak.**
Over all 19,785 covid anellovirus-aligned reads, the genuine component is
≤ 81 reads / 62 CB+UMI (CP95 0.51 % / 2.84 %), and on inspection every one is a
low-complexity tract. A post-hoc refinement leaves 13, none with a body
> 26 nt. But `has_tso` misses edge-truncated TSO,
and `read_body` scores a 5′ TSO/TruSeq fragment as a complex body. Next: the
reagent-leak fix under `ANELLO-PRIOR.3`, which needs one user decision on
column semantics. Then `REL-16` and `ANDET-09f`.

**2026-10-04 (latest): `ANDET-09` closed; the covid TTV signal is an artefact.**
The STARsolo anellovirus alignment branch is implemented and tested, and ships
**off** (`ANDET-09e` criterion 1 failed 1/8 genomes: `SRR2037085_NODE_7436`,
kallisto 14, alignment 0; criteria 2 and 3 passed). The commissioned adversarial
review returned and *upheld the artefact conclusion while correcting its stated
reasoning* — 0/30 read bodies align to any anellovirus genome, and the
pileup/poly-A/`NM = 0` triad is non-diagnostic
(`.living/findings/covid-anellovirus-signal-is-polyg-reads.md`, update
2026-10-04). Next: `ANDET-09f` (filter sweep on the existing plant set),
`REL-16` (the kallisto segfault), and the open user decision on scope for the
covid and SFL-tonsil re-quantification.

**2026-10-03 (night): `ANDET-09` started.** A STARsolo anellovirus
alignment branch that does not depend on kallisto, from a grill of an external
bulk/rustar plan.

**2026-10-03 (latest, evening): DEF-02 chemistry module landed.** `-x` now
defaults to read-based detection and a contradicted `-x` stops the run
(`--force-technology` overrides); SW-21 closed (Drop-seq gets `-w None`); the
3′ strand pilot picks forward. Waiting on HSV-1 job 25695959 for the
end-to-end check. Per-chemistry default *values* stay with DEF-00.

**2026-10-03 (latest): CAT-42 closed; cat42b is the current panel.**
Every swap gate passed (see `CAT-42`). The covid 5′ strand reruns failed in
emptyDrops, because the `test_viralscan` env's `Rscript` has no Matrix. They
were resubmitted with `--resume` as 25695488. Still open in this pass:
GOV-06 closed 2026-10-03. The 5′ reruns finished and `--strand auto` landed on
2026-10-03 (see DEF-02 and F-020). TTMDV reverse-strand read check: done, an artefact (`EXPL-TTMDV`, F-020) (reconciled 2026-10-04).
PROG-07 was re-measured but reopened (`PROG-17`: layer 2 ignored the
called-cell set; both closed later the same day). The REF-13 manifest landed 2026-10-03.
PROG-08 closed 2026-10-03: KSHV is `complete` (LANA cluster, K12, LANA2);
real-data validation is `PROG-18`. PROG-17 and PROG-07 closed 2026-10-03:
layer 2 scores called cells only; EBV LCL = 932 cells, 526 latent, 0 inversions.

**2026-10-02 (latest): finish-pass over every started (`[~]`) row.**
The reviewed plan is at `~/.claude/plans/read-last-handoff-document-moonlit-noodle.md`.
- W0 (housekeeping and blocker notes) is done.
- W1 runs three implementers: `v3/sw24-mask` (SW-24 plus the CAT-42 composition mask), `v3/cli` (SW-06, DOC-02, later `--strand auto`) and `v3/gov06` (LVC-11/12, protocol 1.2.0).
- User decisions: CAT-42 replaces the shipped panel if every check passes; GOV-06 v3 gets an explicit `-gtf` (protocol 1.2.0); no PyPI installs; no push.

**2026-10-01 (latest): CAT-42 panel rebuild running (job 25691642,
homopolymer mask ≥20).** Checks 2 and 3 under `CAT-42` are queued behind it
(25691838 GSE189670, 25691839 EBV/HSV-1). When the build ends, read
`kallisto inspect` (check 1), then the two runs, then return to the
5′ `--strand` reruns and `DEF-00`. `MECH-A2` landed while it ran. Handoff:
[`.living/HANDOFF_2026-10-01.md`](.living/HANDOFF_2026-10-01.md).

**2026-09-29 (latest): every run without `-w` skipped barcode correction
(`SW-13`, fixed).** Snakemake turned an empty whitelist into `kb count -w None`,
which means "bypass". Rerun any result that needs absolute numbers.

The mechanism review (`WP1D`) ranks the next work:

1. ~~`SW-14`…`SW-20` defect fixes~~, all landed 2026-09-29.
0. **Grill decisions confirmed 2026-09-29** — see `WP1E`. They bind all work
   below. On 2026-09-30, MECH-A step 4, `DEF-03` and `--strand` landed.
   Steps 4a–4c (the guard false positive, the NC_000898.1 catalogue row,
   old-run backfill) also landed, which closes MECH-A. Next: the
   5′ reruns with `--strand`,
   then `DEF-00`, the protocol amendment.
2. `MECH-A`, the Virus Identity table, is done (2026-09-30). Step 1 (the catalogue
   merge, 4,128 accessions with taxid and `panel` scope) and step 2
   (`virus_identity.py`, golden-tested on 3 real indexes) landed on 2026-09-29.
   Step 3 (the `analysis` rule writes `results/virus_identity.tsv`) also
   landed. Step 4, re-pointing the consumers, landed 2026-09-30, and so did
   steps 4a–4c. Both guard rules are decided (grill Q8 = B, adding an
   index manifest; Q8b = A2, genus plus "Anelloviridae (genus unassigned)").

The 4,127-genome max panel (`viral_panel_max_2026-09-28/`) built and passed
every gate, but must not become the default before `MECH-A`/`MECH-B`. EBV
type 2 took 16 % of EBV molecules through the equal split of shared k-mers
(F-017).

**2026-09-28: the final reference is BUILT — 2,343 genomes.** `CAT-31`
through `CAT-35` shipped together; the panel is 323 GTF-backed genomes plus
2,020 anelloviruses, and the CAT-31 guard closed 32 of 34 catalogued gaps inside
that build. Full measurements in
[`.living/findings/final-reference-panel-2343-genomes.md`](.living/findings/final-reference-panel-2343-genomes.md)
(F-016). Four things now matter, none of them "add more viruses":

1. **`CAT-40` — gate the index on `kallisto inspect`** (max EC size, discarded
   ECs) before trusting it. The HPV k-mer metrics are computed outside kallisto
   and cannot see the pseudo-inverse threshold effect, which is `CAT-19`'s real
   failure mode. This is the one open item that could still invalidate the build.
2. **`CAT-36` — the catalogue is behind the panel by 96 rows**, so 96 indexed
   genomes have no provenance, family or tier. `CAT-31`'s reverse check will keep
   reporting them.
3. **`CAT-38`/`CAT-39` — three whole-genome pseudo-transcripts** entered the
   panel (`NC_001489.1` HAV, `NC_001488.1` HTLV-2, `M12737` HPV-8). This is the
   exact shape that collapsed 99.8 % of anellovirus UMI into one bucket, so they
   need a deliberate decision, not a silent default.
4. **`CAT-30` is still the only row blocked on you** — retract the published
   Alphatorquevirus number, or leave the published results frozen? The science
   is closed; the retraction never happened.

**2026-09-27 (later): F-015 — the catalogue is not the index, and it outranks
every reference-import question above.** Measured against the built
`viral_genome.dedup.fa`: **13 of 16 catalogued HPV are absent** (including
HPV18 and HPV31), and **all 5 Retroviridae — both HIV-1 and HIV-2 — are
absent**, plus 7/30 influenza, 4/7 coronaviruses, 3/8 polyomaviruses. A
catalogued virus that is not indexed is *undetectable*, so this silently caps
sensitivity today and would do the same to any import.

1. **`CAT-31` shipped 2026-09-28 (guard only — the 34 misses are still open).**
   The build now reconciles the catalogue against the panel it actually emits
   and writes `catalogued_not_indexed.tsv`, so the gap is reported on every
   build instead of being invisible. The guard is deliberately **red**: it
   reports all 34 catalogued-but-unindexed accessions as `unexplained` and fails
   under `--strict-reconciliation`, because none of them is a documented
   exclusion — they are F-015's measured sensitivity loss. Each one now needs a
   genome or a named decision, not a silent default. Next: work the report.
2. **`CAT-32` is the cheapest real EBV win.** EBV type 2 (`NC_009334.1`) is
   missing; 31-mer overlap with type 1 is 79.9% shared, leaving ~20%
   type-discriminating k-mers. One accession turns EBV detection into EBV
   *typing*. The `CAT-31` report will confirm it lands.
3. **`CAT-33` reframes the HPV question.** VIRTUS2 offers 92 HPV *types* vs our
   16 → 76 genuinely new, but as weak-provenance `lcl|` records; source the same
   types from curated RefSeq `NC_` accessions. Start with the **13 already
   catalogued** types in the `CAT-31` report — that is 13 real genomes, not 76
   speculative ones.
4. **`CAT-34` closes a worry.** No `chrEBV`/viral contig exists in our actual
   host build (GRCh38-2024-A, 194 contigs; GRCh38.116 GTF), so the VirDetect
   silent-invisibility hazard does not apply here.
5. **VIRTUS2's other 400 records are not worth importing** (`CAT-35`, done):
   legacy influenza lab strains and 27 anelloviruses we already cover better.
   Take the *method* (host subtraction, strand-aware counting), not the genomes.

**2026-09-27: `CAT-17` shipped and the panel expansion was adversarially
reviewed — five new rows (`CAT-25`…`CAT-30`), four retractions, 41 external
files vendored.** Three things now need a decision rather than more work:

1. **`CAT-30` is the most urgent and is blocked on you.** `main` @ `4fcd748`
   already closed F-005 as an accession-level artifact, but
   `covid_viralscan/results/*/results/viral_summary.tsv` **still publishes**
   `Alphatorquevirus` 1,167,103 and 1,605,631 UMI at 17.6–19.0 % of cells. The
   published number and the closed finding contradict each other today.
2. **`CAT-05` has a live blocker.** The new gate rejected *both* the upstream
   2,023-rep panel and the shipped panel on one duplicate —
   `NC_038359.1` / `AB303562.1`, byte-identical. It is an upstream
   dereplication failure we inherited, and a duplicate sequence has broken
   `kallisto index` before.
3. **`CAT-20` stays gated** on the F-013/F-011 reconciliation, per the
   previous session's own instruction. Adopting the 2,023 reps is the same
   unmeasured question as restoring the 1,522, approached from the other side.

An advisory note on sequencing: the panel question was framed as "cover more
viruses", but for anelloviruses — 88 % of the panel — more genomes is provably
*not* the answer, because the ceiling is set by upstream dereplication. The
half of the request with a clear answer is breadth outside Anelloviridae
(`CAT-21` SARS-CoV-2 lineages, `CAT-22` influenza clades, `CAT-27` EBER,
`CAT-29` HSV-1 latency), and that half is cheap.

**Earlier, superseded:** WP4H gained a design document and 8 rows
(`CAT-17`…`CAT-24`) — read
[`docs/plans/2026-09-27-viral-reference-panel-expansion.md`](docs/plans/2026-09-27-viral-reference-panel-expansion.md)
before touching `CAT-09`…`CAT-16`.** It corrects three rows as
under-specified and adds six failure modes that were not tracked:

1. **`CAT-17` is the next thing to build.** A naive whole-genome anellovirus
   panel manufactures **1.44 % of R2 reads as false hits** in the EBV LCL
   `SRR12682296` — 6,437/6,437 captured hit reads had *zero* genuine
   anellovirus k-mer once low-complexity k-mers are masked. It is cheap, has
   no dependencies, and every panel expansion silently reintroduces it.
2. **`CAT-12` is not executable as written.** RefSeq holds **exactly one**
   SARS-CoV-2 genome and **zero** lineage-labelled RefSeq records
   database-wide, so "one per WHO variant lineage" cannot be built from
   RefSeq. → `CAT-21`, blocked on a user decision (document §8 Q1).
3. **"One influenza strain per subtype" is the wrong axis** — the diversity is
   in the clade, and the 8 current IAV RefSeq records are all segments of *one
   1934 lab strain*. → `CAT-22`.
4. **Panel expansion is not monotone-good** (`CAT-19`): `kb count` *discards*
   multimapping UMIs by default in the scRNA-seq path, and the one benchmark
   that measured this found F1 falls ρ = −0.73 with reference-set size for
   >99 %-identical genomes.
5. **The anellovirus panel is structurally incomplete** (`CAT-20`): genus-name
   queries reach 21,283 of 42,755 family records; SENV sits outside every
   genus; ICTV is now 37 genera / 243 species, not the 8 the shipped table
   encodes.
6. **The Serratus anchor becomes checkable** (`CAT-23`): the 18 latent-DNA-virus
   accessions of Lareau et al. *Nature* 2023 are already in the panel, so they
   are an acceptance criterion rather than new breadth.

`CAT-18` (two-index architecture) is the structural fix for the 99.8 %
single-bucket failure that `CAT-01` did not address: whole-genome and real-CDS
transcripts for the same virus still share one equivalence-class space.

**New, 2026-09-27: WP4H gains a breadth-and-diversity expansion
(`CAT-09`…`CAT-16`).** A census of the current reference found 2,313 records but
only ~106 non-anellovirus species, nearly all single-genome, with SARS-CoV-2,
HIV-1/2, HTLV-1/2, 3 of 4 seasonal coronaviruses, hMPV, bocavirus, influenza D
and every high-risk HPV except 16 **absent entirely**, and influenza A present
only as one 1934 lab strain. `CAT-01` (build-ref discards the real GTF) and
`CAT-03` (nothing groups segments) block all of it and come first. Diversity
targets are set by measured 31-mer capture, not by quota.

**New, 2026-09-27: `WP4J` — the SFL tonsil pool shows no TTV and no HPV in any
of the 24 donors (`TONSIL-01` done).** The screen used the host-subtracted
cellranger BAM: 118.9 M unmapped GEX reads, and 23 reads survived to the end,
none of them anellovirus or HPV. A positive-control plant recovered 99 %. A
held-out plant recovered only 8–40 % of reads from TTV strains ≤ 84 % identical
to any reference, so divergent low-level TTV is not ruled out. Only
`x223` is gene expression (`x225` is ADT), and the library is 5′ v3 R2-only.
None of these zeros is an informative negative (`SENS-06`). A native `viralscan`
run waits for `TONSIL-02` (strand and whitelist support). The WP4I → WP4G → WP4H
order below is unchanged.

**New, 2026-09-27: three work packages opened, in this order — `WP4I`
(latent/lytic), `WP4G` (anellovirus detection), `WP4H` (human-virus catalogue).**
`PROG-10` is done; next is `PROG-11` (catalogue
biology), because the `complete` labels and the `PROG-07` numbers are unreliable
until they land. `CAT-01` (build-ref discards the real GTF) blocks every natively
built index, so it precedes any `CAT-08` build and `ANELLO-13`. Do not flip the
anellovirus default (`ANDET-07`) before `ANDET-01`–`ANDET-04`. Housekeeping that
preceded this: the WP4B2–WP4F work was squashed into one commit so every commit
passes its own suite, and the schema 1.1.0 governance migration was finished
(see the evidence log).

**New, 2026-09-26: WP4F (`ANELLO-01`–`ANELLO-13`) landed** — the Anelloviridae
panel now carries real NCBI gene structure instead of one placeholder gene per
genome. The root cause was a *discarded* GenBank GTF, not a missing one:
`_fetch_one` already wrote a real annotation for 1,995 of 2,042 accessions and
both panel builders threw it away and rebuilt a placeholder from the FASTA. A
second, independent defect would have survived that fix — gene IDs were the bare
`/gene=` value, so `ORF1` (150 genomes) would have collapsed 1,995 genomes onto
2,316 columns. `extras/build_anellovirus_genes.py` generates
`src/viralscan/data/anellovirus_genes.tsv`; measured CDS coverage is 205/206
(99.5 %) across all eight genera, and the low CDS *count* is biology, not a gap:
75 % of the panel is Betatorquevirus TT-mini genomes that genuinely carry one
ORF. NCBI carries **no** genogroup for Anelloviridae, so that column ships empty
rather than inferred. **WP4F also fixes the two things `HPV-09` was blocked on**
(genome-scoped `_genbank_to_gtf` IDs, public `fetch_genbank()`); that row is left
for its owner to close against their own tests. `ANELLO-12`/`ANELLO-13` record
the limit: the 99.8 %-in-one-bucket covid artifact is now explainable and
testable, but the test needs `kb` and a re-run, and per-genotype anellovirus
quantification is still not defensible for the one-ORF majority of the panel.

**New, 2026-09-26: WP4E (`HPV-01`–`HPV-10`) landed** — HPV ORFs are now named.
All 14 high-risk genotypes plus HPV1/HPV2, 124 ORFs in
`src/viralscan/data/hpv_genes.tsv`, names taken from each record's own `/gene` or
`/product` qualifier rather than from a coordinate table, because papillomavirus
genomes are linearised circles cut at the submitter's choice of point and a
coordinate table is wrong for 15 of the 16 genotypes here. Two follow-ups are
recorded rather than done: `HPV-09` needs a `ncbi_fetch.py` fix (out of scope
for that change) and `HPV-11` needs the real index rebuilt to *measure* the L1
cross-mapping this row currently only predicts. **Do not publish a per-type HPV
number before `HPV-11` closes.**

**New, 2026-09-26: WP4D (`PROG-01`–`PROG-07`) landed** — layer 2 gene-programme
inference, opt-in via `--gene-programs`, with the EBV LCL regression test
pinning the design. `PROG-08`/`PROG-09` are the open scope questions. WP4E does
not change `PROG-09`: HPV still has no latency/lytic dichotomy, so it stays out
of the gene-programme catalogue. Also in this change: WP4C (`SENS-01`–`SENS-05`)
and WP4B2 (`HOST-01`–`HOST-04`) landed alongside the `REF-01`/`REF-13`
reference-visibility fixes.** An audit of detection sensitivity found that
`limit_of_detection` existed only as an unrun endpoint in the protocol schema,
that the bundled 20-genome TTV panel captured
**1.36 %** of the real anellovirus 31-mer space (median per-genome coverage
**0.00 %**, 85.8 % of genomes sharing zero 31-mers), and that 49.2 % of gene IDs
in the covid panel resolved to no virus name — so 9 of 17 published
`viral_summary.tsv` rows were bare gene IDs rather than viruses.

The immediate consequence: **`informative_negative` is `false` on essentially
every run, by design**, because certifying a negative needs a measured k-mer
capture term and no shipped workflow plants a control. `SENS-06` is the row that
fixes that, and it should be built together with `VAL-01` — a generator that
plants a target at known abundance *is* the positive control `SENS-04` consumes.
Until then, no negative result from this package may be reported as an absence.

**Governance follow-up required at commit time (not done here, deliberately).**
`analysis/v3_artifact_inventory.tsv` has no rows for the artifacts this work
introduces, and adding them before the commit would break
`test_artifact_inventory`: `governance_utils` resolves each row's `path` from
the **git tree at the recorded `git_sha`**, so a row pointing at an untracked
file can never validate. At commit, add rows for
`src/viralscan/sensitivity.py`, `src/viralscan/gene_programs.py`,
`src/viralscan/data/gene_programs.tsv`, `src/viralscan/scripts/gene_programs.py`,
`results/sensitivity.tsv`, `results/positive_control.json`,
`results/gene_program_summary.tsv` and `results/gene_program_cells.tsv` (and
re-derive the five already-failing rows:
`artifact-inventory-schema`, `claim-registry`, `claim-registry-schema`,
`output-reference-doc`, `v3-counting-contract-doc`, which fail today only
because those files are modified-but-uncommitted in the working tree).

Pre-existing and still first in the queue: **`G0` archive-build evidence, then
`G1`.** `GOV-03`,
`GOV-04`, and `GOV-05` now have fail-closed local evidence: the sanitized
artifact inventory, claim graph and public ship-scope validators pass, and a
wheel plus sdist built directly through the installed setuptools backend match
the positive allowlist. `G0` remains partial because the locked Python lacks the
PyPA `build` frontend required by the frozen validation command. No dependency
was downloaded or installed without approval.

### Do now

**Superseded 2026-10-04.** Follow
`docs/plans/2026-10-04-package-completion-plan.md` (M0 → M4). The user supplies
network installs (the `build` frontend, Sphinx, test tools, locks), the push and
CI. The original list is kept below for its G0 commands.

1. Supply a compatible PyPA `build` frontend to
   `benchmark_runs/legacy_v2_v3/env_full` through an explicitly approved network
   install or a user-provided offline artifact.
2. Run the prescribed `python -m build --no-isolation` command and the exact
   wheel/sdist/context membership check; rerun the full local gate and frozen
   dirty-path hashes before promoting `G0`.
3. Treat the separate integration-marker failures as `G1` evidence: the real CLI
   currently requires a `conda` executable for `--use-conda`, and this runtime
   does not provide one; the installed Snakemake also emits no rule listing with
   the tests' `--quiet` dry run.

## Software backlog reconciliation (2026-10-09)

Merged from the user's 59-item backlog (checkout `b1a494e`) and the open-row survey of this
tracker. Software only: tiny fixtures, install checks and regression tests count; dataset reruns,
threshold calibration, benchmarks and holdout evaluation do not. Rows already tracked below keep
their ID and are only pointed at; new rows are `QA-*`/`PKG-*`. Order is priority: tiers 1-3 first.
Verified 2026-10-09: ruff 2 findings, `ruff format --check` 16 files (green on 2026-10-05, so a
regression), mypy 166 errors in 21 files, 1,930 unit passes. `DEF-05` reading `Gene/raw`
was checked and closed 2026-10-10 (`cellcalling.solo_raw_dir` rejects it; `tests/test_cellcalling.py:465`). Unverified: `EVID-CORR-01` (`evidence.py:1025` still has
`dedup="umi"`/`strand_aware` parameters).

**Tier 1 — gates and correctness**
- [x] `QA-01` — fix the 2 ruff findings (`test_chemistry_check.py` computed default,
  `test_snakefile_dag.py` membership expression). **Done 2026-10-09:** `ruff check .` clean; the B008 default
  is now a module-level `_SHARED_TAIL` (same one-tail-per-import behaviour); E713 auto-fixed.
- [~] `QA-02` — `ruff format` the 16 files and fix the 166 mypy errors without blanket ignores. **Format done
  2026-10-09** (15 files; `tests/test_dsr02_enumerate_calls.py` is the user's dirty file and is left
  unformatted; 1,930 unit passes after). **Left:** 166 mypy errors (`evidence.py` 81, `anello_align.py` 30,
  `chemistry.py` 25, `multimap.py` 13, ...), assigned per file to the agent that owns that file.
- [x] `PKG-01` — integrated `a192ef6` (2026-10-09): wheel/sdist and Docker-context
  membership reconciled, tracked-tree regression added. Claude's offline artifact verification
  refuted the wheel-import premise; the Docker context and allowlists were the actual defects.
  Untracked user-owned `reads.py` is outside this tracked-tree completion claim.
- [x] `SENS-CORR-03` — integrated `69deaef`/`11ca0da` (2026-10-09): canonical
  `SensitivityRecord` drives TSV/JSON/text/HTML; capture 0.5 distinguishes 2.996 expected
  observed molecules from 5.991 required true molecules. Fresh combined verification pending.
- [x] `SENS-CORR-01` (residue) — integrated 2026-10-10: `exact_sequence` certifies only with a
  checksum-bound receipt (accession, sequence, index and manifest digests, evidence); targets resolve through
  `VirusIdentityTable`; immutable `control_receipts/<sha256>.json`; report/receipt schemas cover verified/
  failed/unavailable and reject forged "certifying" reports. The 8 tests that failed built a measured control
  with no receipt, so they now use the shared `tests/_exact_control.py` helper; new tests prove no receipt or a
  changed index digest certifies nothing. 222 tests pass. The stale "0.32/0.06" text is already gone
  (`output_reference.md:304` reads 0.063/0.013). Left: document the receipt fields in `output_reference.md`.
- [x] `EVID-CORR-01` — integrated `b45d1b6` (2026-10-09): one deterministic
  representative per corrected CB/UMI/reference, shared by coverage, starts and QC.
  Fresh combined verification pending.

**Tier 2 — chemistry, counting, provenance**
- [~] `CHEM-01` (closed 2026-10-10: the bulk path stays advice-only), `DEF-05` (done: reader's matrix verified, missing host matrix rejected),
  `DEF-06`, `DSR-08`, `SW-07`, `MECH-F` (index-aware denominators, gene-role contract).

**Tier 3 — packaging and reference availability**
- [ ] `REL-04`/`REL-05` repeat on the finished source (clean installs, `doctor`, `validate-run`);
  `REL-03`/`DEF-08`/`CAT-25` (review the untracked `uv.lock`); `CAT-06` (index carries its GTF and
  catalogue); `REF-11` and `PANEL-01` WP4 bundling (deferred; the DOI stays external).

**Tier 4 — reference construction and panel**
- [~] `REF-03` (2026-10-10: `prepare_reference_inputs` now rejects a reference with no viral records when masking
  is requested, with a test; the gtf-manifest check in the NCBI-fetching builder run is network-pending),
  `REF-02`/
  `REF-05`/`REF-07` (real provenance, retrieval vs build time), `REF-04` (deterministic rebuild).
- [~] `PANEL-01` WP2 **done 2026-10-10** (SV40 `NC_001669.1`, AAV1/7/8 `NC_002077.1`/`NC_006260.1`/`NC_006261.1` fetched from NCBI,
  catalogued as `shipped` + `contaminant`, GTFs generated locally and listed in `gtf_corpus_manifest.tsv`; 2,272 tests pass) and **WP3 gates**
  (self-consistency matrix, EC-size gate, `CAT-40` `kallisto inspect`); iciHHV-6 flag and CAR-T vector
  filter under `CAT-04`. B19 genotype 2 stays deferred (no complete reference).
- [~] `REF-13`, `CAT-36`, `CAT-38`/`39` (done 2026-10-10), `DEF-04` (deferred, see its row) (annotation integrity); stale `AB303562.1` rows in
  `anellovirus_genes.tsv` (removed in `f4320d2`, 2026-10-10); GTF-manifest check in the build lifecycle
  (done: `build_bundled_panel_ref.py:_check_gtf_manifest` runs `write_gtf_manifest.py --check` before any fetch;
  only the network-backed builder run is left).

**Tier 5 — competitive evidence and diagnostics**
- [x] `EVID-COMP-01..05` — done 2026-10-10: manifest validator, class-encoded FASTA, tie-aware
  BLAST with no-hit rows, sampled-read vs molecule accounting, CLI/schema/docs. Competitor mode omits the
  two-sided HOST/VIRUS molecule verdict, recorded as `competition.molecule_verdict`. 119 unit and 2 live
  integration tests pass. `molecule_verdict`, the verdict files and the competitor-mode outputs are now in
  `docs/output_reference.md` (2026-10-10); governance re-pin of that file is the follow-up commit.
  Remaining rows in the old line: (manifest validator → reference classes → complete BLAST ties/no-hits →
  lineage accounting → CLI/schema/docs); `ANDET-01`/`ANDET-02` (genome breadth, host-homologous
  concentration); `ANDET-04` (auto-run rule; trigger is a user decision); `CARD-01`; `CAT-04`,
  `PROG-09`/`PROG-12`/`PROG-13`, `CAT-28`.
- [ ] Added from the survey, not in the user's list: DSR tooling that is code (`DSR-02`, `DOSSIER-01..03`,
  `VERDICT-01`, `SENS-07`, `DLIST-02`), `ANDET-07`, `TONSIL-02`, `CAT-19`, `DEF-07`,
  `CAT-14`/`15`/`16`.

**Tier 6 — host response** (`HR-01..05`: group-disjoint folds, same-fold baselines, structured null,
cell-type restriction, provenance and entry-point parity). Spec: `viralscan_lshwcp_spec/04_HOSTRESPONSE_ROBUSTNESS.md`
and `SOURCE_AUDIT_IMPLEMENTATION_PLAN.md` (in the `$SITE_SCRATCH/ViralScan` checkout). The core was recovered
into `hostresponse.py` by `c4e7d2e` (2026-10-10) without rows here; reviewed against the spec on 2026-10-12 (branch `codex/hr-host-response`):
- [x] `HR-01` — group-disjoint folds (`--cv group --groups COL --cv-folds N`), canonical row-order-invariant folds, fold tables, in-fold HVG,
  frozen aggregation (`AUC_STATISTIC`), no fallback to cell CV. **Verified 2026-10-12:** acceptance items map to existing passing tests in
  `tests/test_hostresponse_groups.py` (no group leakage `test_no_group_in_train_and_test`/`test_no_leakage_in_saved_membership`; reorder invariance
  `test_row_order_*`; invalid train/test classes `test_invalid_fold_classes_after_filtering_fail_explicitly`; metadata failures
  `TestMetadataValidation`/`test_missing_group_column_fails_without_fallback_or_cleanup`; fold-local features `test_hvg_selection_never_sees_held_out_cells`;
  legacy cell mode `test_cell_mode_is_default_and_unchanged`). The 9 host-response test files pass (247). Parity of the HR-01 flags across the
  standalone parser, subcommand, main CLI and Snakefile is `HR-05`.
- [x] `HR-02` — structured null (`--permutations`, `--permutation-unit`, `--permutation-block`; full refit per replicate, failed replicates retained, pinned
  empirical-p formula). **Done 2026-10-12:** the null-spec rules (declared unit, block consistency, group CV) were only checked inside the `hostresponse` rule,
  i.e. after `kb count`; they now live in `src/viralscan/hostresponse_contract.py` (`validate_evaluation`, import-light) and run in `RunConfig.from_snakemake_config`,
  `menu.errorhandler` (clean exit before the output directory exists) and `run_hostresponse`. Tests: `tests/test_hr_entrypoint_contract.py::test_config_rejects_an_invalid_null_specification`,
  `::test_errorhandler_rejects_an_invalid_null_specification`, `::test_a_declared_group_null_passes_the_config` (the first two fail on the previous config, which accepted all four
  bad combinations); null/permutation behaviour itself is covered by `tests/test_hostresponse_null.py` (empirical-p pin, block preservation, determinism, full refit, invalid-null fixtures).
  Naming note: the summary column is `observed_cv_auc` (spec text says `observed_group_auc`); not renamed.
- [x] `HR-03` — cell-type strata (`--cell-type-column`, `--cell-types`; separate cohorts, per-stratum support checks, collision-checked `cell_type_*` directories, status rows;
  existing tests in `tests/test_hostresponse_strata.py`). **Done 2026-10-12:** `--cell-types` was comma-split (no strip: `"T, B"` gave `" B"`, which then failed as an unknown type) on
  the standalone parser but `nargs="+"` on the subcommand and main CLI, while the FAQ shows the comma form (`viralscan hostresponse --cell-types T_cells,B_cells` arrived as one unknown
  type). All three now use `hostresponse_contract.CellTypesAction`: space- and/or comma-separated, stripped. **Ambiguity resolved conservatively:** a name containing a comma cannot be
  given on the command line (use the config); an unknown name still fails loudly listing the available types. `cell_types` without `cell_type_column` is now rejected in the config and
  `menu.errorhandler` (before kb count), not only inside the rule. The standalone-parser test that pinned `" B"` now pins `"B"`. Tests: `tests/test_hr_entrypoint_contract.py::test_cell_types_syntax_is_the_same_on_every_entry_point`,
  `::test_cell_types_without_a_column_fail_before_kb_count`. `docs/cli_reference.md` generated region regenerated (help text).
- [x] `HR-04` — same-fold baselines (`depth_only`, `cell_type_only`, `depth_plus_cell_type`, optional in-fold `panel_depth` via `--panel-in-fold`; training-only scaler/category
  encoder/panel, unseen test levels map to zero indicators; headline columns `expression_auc`/`depth_auc`/`cell_type_auc`/`depth_cell_type_auc`/`expression_minus_best_baseline_auc`).
  **Verified 2026-10-12** against the spec by the existing `tests/test_hostresponse_baselines.py` (identical folds/cohorts `test_every_model_uses_identical_folds_and_cohorts`,
  `test_same_fold_holds_with_depth_filtering`; training-only objects `TestBaselineDesign`, `TestInFoldPanel`; unseen cell type `test_unseen_test_cell_type_is_handled`; cell mode on the
  expression model's splits `TestCellModeBaselines`). **Done 2026-10-12:** `panel_in_fold` without group CV was only caught inside the rule (after kb count); now also rejected in the config
  and `menu.errorhandler` via `hostresponse_contract.validate_evaluation`. Test: `tests/test_hr_entrypoint_contract.py::test_panel_in_fold_without_group_cv_fails_before_kb_count`.
- [x] `HR-05` — provenance and entry-point parity. **Verified/done 2026-10-12:** `hostresponse_manifest.json` records every spec field (`hostresponse_schema_version`, `sklearn_version`,
  `cv_mode`, `group_column`, `n_folds`, `seeds`, `permutation_count`, `permutation_unit`, `feature_mode`, `depth_matching`, `cell_type_filter`) plus input and output SHA-256, fold digests and
  per-virus/stratum/skip/error statuses (`hostresponse_status.tsv`); owned-output cleanup covers the new TSVs and `cell_type_*` directories (`tests/test_hostresponse_strata.py`,
  `tests/test_hostresponse_groups.py`); the settings fingerprint, rerun inheritance and config round trip are in `tests/test_hr_entrypoint_contract.py`. Parity gaps closed in `HR-02`/`HR-03`/`HR-04`
  (one shared validator, one `--cell-types` syntax). Also: parser help and `docs/cli_reference.md` cited "Luebbert et al. 2026" (now 2025; the module header already said 2025) and the standalone
  `--differential` help said "ALL features" although the table covers the HVG set unless `--no-use-hvg` (B-18); module docstring now lists the manifest, status, fold, null and stratum outputs.
  New module `hostresponse_contract.py` added to the wheel and sdist allowlists in `config/public_ship_scope.json` (`test_ship_scope_tracked_tree`).
  Test: `tests/test_hr_entrypoint_contract.py::test_help_cites_luebbert_2025_and_states_the_differential_family`. No hostresponse completion-marker coverage was added (artifact-local manifest only, as specified).
  **Documentation complete 2026-10-13:** the host-response section now covers the manifest, status, CV, permutation and stratum outputs;
  barcode alignment describes the approved host-cell zero-fill policy, and the changelog cites 2025. The output-reference pin is refreshed after the documentation commit.
  A host-response manifest schema remains optional and is outside this acceptance contract.
- **B-18/B-19 (AUD-13):** B-18 done (`fdr_family`, `n_genes_tested`, `inference_unit`, `donor_adjusted` columns in `<virus>_differential.csv`);
  B-19 is implemented and documented: host cells absent from the viral matrix are zero-filled, viral-only barcodes are dropped, and the manifest records both counts and fractions.

**Tier 7 — containers, release, docs**
- [ ] `REL-06..14`, `RC-01..03` (one OCI artifact, SIF from its digest, parity, product-environment audit,
  SHA-pinned Actions, exact-revision gating); `DOC-01` (Sphinx `-W`), `DOC-05`, `DOC-07`, `OPS-01`/
  `OPS-02`, `REL-15`, `STB-04`. Maintainer credentials, `REL-15` and `REF-11` stay user-owned (`[!]`).

**Tier 8 — deferred architecture and scale**
- [~] `SW-03` (done 2026-10-10), `MECH-C` (done), `SW-09` (parsers done, handlers deferred), `MECH-E` (open), `CAT-18` (closed, not adopted);
  `EMC-01` (done 2026-10-10; the cluster re-measure is the user's); optional chunked BUS ingestion (`QA-03`, **closed 2026-10-12 without a chunked mode**: ingestion is already one streaming pass over `bustools text` rows (SW-03), so no chunk boundary can split a CB-UMI molecule; the acceptance proof is `tests/test_qa03_bus_ingestion.py` (46 seeded property tests: 7 read-buffer sizes incl. mid-row cuts x every multimap method give identical layers and audit, unsorted text is refused, in-memory row order is irrelevant) on top of the DEF-06 properties in `tests/test_multimapping_properties.py`; a real chunked reader is only worth building if a deep-sample RSS measurement shows the stream is not enough).

**Tier 9 — validation tooling (software only, deferred until protocol readiness)**
- [~] `VAL-01`, `VAL-02`, `VAL-03`, `VAL-06`/`VAL-07`, `VAL-08`/`VAL-09` (audited 2026-10-12: `VAL-07` closed, the rest `[~]` with what needs data, the cluster or a user decision, see the WP5 rows); conditional `VAL-RA-01..04` and
  `VAL-RA-CAL` (keep `virus_key` certification rejected until matching evidence exists).
  - **`VAL-RA-01..04` / `VAL-RA-CAL` audit 2026-10-12 (`codex/val-software-audit`):** these ids have no row text beyond this line and
    `VAL-RA-DECISION-01`; they are conditional on measured transfer evidence (`SENS-CAL-01`, cluster, and the user's protocol
    amendment), so nothing is implemented and they stay `[ ]`. The guard is verified and one bypass is closed: `virus_key` was refused
    by `RunConfig.from_snakemake_config` but accepted by `RunConfig.from_yaml`, which detection and the `rerun-*` commands use, so a
    hand-edited `config.yaml` with `positive_control_scope: virus_key` would have let `control_claim`/`scoped_capture` certify the
    target. `from_yaml` now raises the same `VIRUS_KEY_REJECTED` error. Tests `tests/test_val_ra_guard.py` (the YAML test fails on the
    old code) and `tests/test_positive_control.py::TestScopeConfig::test_virus_key_scope_rejected_without_transfer_calibration`.
    `sensitivity.scoped_capture` stays a pure function (its `VIRUS_KEY` branch is tested for the future calibrated case); the gate is
    config-level. Lift the guard only together with an approved calibration and a reviewed amendment.

**Tracker hygiene (same commit as the first row above):** the LSCHWCP `IMPLEMENTATION_BACKLOG.tsv`
still lists `SENS-CORR-01/02` as TODO (both done 2026-10-06); `SESSION_RESUME.md:230` and `HANDOFF.md`
(says branch `main`; current is `codex/panel-reference-roles`) are stale. Not missing, do not restart:
role reporting, scratch isolation, the exact substitution model, chemistry geometry, read extraction,
molecule verdicts, TSO trimming, the three modelled no-CDS GTFs.

## How to use this tracker

1. Take the first unchecked item whose dependencies are all `[x]`.
2. Add or update its test before changing production behavior.
3. Run the stated acceptance command and save the named artifact.
4. Mark `[x]` only with dated evidence; use `[~]` for partial and `[!]` for an
   external blocker.
5. Update **Next action**, the work-package row, and the claim registry in the
   same commit.

Status legend: `[x]` verified; `[~]` partial; `[ ]` ready/not started; `[!]`
blocked on a person, credential, private datum, external service, or unavailable
compute.

## Critical path

| Order | Work package | Status | Depends on | Exit gate |
|---:|---|:---:|---|---|
| 0 | Governance and legacy freeze | `[~]` | none | `G0` |
| 1 | Core software contracts | `[~]` | `G0` | `G1` |
| 2 | Install, lock, and artifact parity | `[~]` | `G1` | `G2` |
| 3 | Preregister scientific validation | `[~]` | `G0` | `G3` |
| 4 | Freeze production references | `[~]` | `G3` | `G4` |
| 5 | Build and validate the truth panel | `[ ]` | `G3`, `G4` | `G5a` |
| 6 | Run public positives and comparators | `[ ]` | `G2`, `G4`, `G5a` | `G5b` |
| 7 | Score, audit, and freeze results | `[ ]` | `G5b` | `G5` |
| 8 | Regenerate docs and claims | `[~]` | `G1`; numbers require `G5` | `G6a` |
| 9 | Regenerate and review manuscript | `[ ]` | `G5`, `G6a` | `G6` |
| 10 | Publish and test `3.0.0rc1` | `[ ]` | `G2`, `G5`, `G6a` | `G7` |
| 11 | Publish stable `3.0.0` and submit | `[ ]` | `G6`, `G7` | `G8` |
| 12 | Observe and maintain | `[ ]` | `G8` | ongoing |

Work packages 1-4 may proceed in parallel. Do not start truth-panel outcome
analysis before `G3`, regenerate quantitative documentation before `G5`, or tag
a release before its gate.

## Verified baseline

- [x] `BASE-01` — v3 branch and locked correctness roadmap exist.
- [x] `BASE-02` — molecule-safe streaming allocation, distinct-gene EC
  projection, five explicit methods, contracted H5AD layers, and exact mass
  audits pass synthetic/property tests.
- [x] `BASE-03` — retained EBV input completed with 103,145,071 BUS records,
  57,957,364 input molecules, and exact partition/matrix conservation. Evidence:
  [`analysis/v3_ebv_baseline/README.md`](analysis/v3_ebv_baseline/README.md).
- [x] `BASE-04` — exact-target evidence completes extraction, competitive
  alignment/BLAST, deduplication, QC, coverage, BAM indexing, and IGV on the tiny
  real-tool fixture.
- [x] `BASE-05` — the 2026-07-22 local gate passed 720 tests with 21 deselected;
  Ruff check/format, data governance, and protocol validation also passed. Rerun
  the full gate after each remaining software slice rather than treating this
  development-branch result as release evidence.

## WP0 — Governance and legacy freeze

Objective: prevent unsafe pre-v3 results, private data, or unsupported claims
from entering the release. Estimated remaining effort: 1-2 days.

- [x] `GOV-01` — mark pre-v3 corrected/combined counts scientifically
  incompatible and ineligible for v3 claims.
- [x] `GOV-02` — add a ship-scope data-governance check and exclude private
  clinical outputs from package/container/source-distribution scope.
- [x] `GOV-03` — inventory every claim-bearing input, reference, intermediate,
  result, and scheduler record with software version and SHA-256. Write
  `analysis/v3_artifact_inventory.tsv`.
- [x] `GOV-04` — define the public ship-scope allowlist, then remove or quarantine
  institutional defaults from every included script and manifest.
- [x] `GOV-05` — expand `claims/registry.json` into a validated claim graph with
  source location, artifact digest, Git SHA, input/reference hashes, schema,
  layer, denominator, generation command, scope, and status.
- [x] `GOV-06` — **closed 2026-10-03 (diagnostic-only, outcome-ineligible).**
  LVC-11 to LVC-14 are `[x]` in `analysis/legacy_v2_v3/TRACKER.md`. The
  interpretation is in `FRESH_CONTROLS.md` (merged from `v3/gov06-close`,
  17f59f3).
  - **v2:** all five fresh rows reproduce the archive on every record within
    1e-6, including the attempt-3 highmem row SRR6825024.
  - **v3:** all five fresh rows completed with `validate-run` 0. They used
    the frozen wheel and an explicit `-gtf` (protocol 1.2.0, packet
    `66918c50…`). Exact-read evidence succeeded for every expected target.
  - **Fresh v3 vs archived v3:** 836 of 862 records match. The 26
    mismatches:
    - 20 HIV reference differences (the arms differ on HIV);
    - 3 Cercopithecine herpesvirus 1→0 molecules (unresolved);
    - 3 EBV: 2 float noise, and SRR12682298 at −2 molecules (unresolved).
  - **HIV:** v2's HIV result is structurally negative, not a sensitivity
    comparison.
  - **EBV:** v3 is 0.46–0.59× v2. These are different quantities and are
    not rescaled. EBV recovery is qualitative positive-control evidence
    only.
  - **Caveats:** the frozen v3 wheel predates SW-13 and carries SW-24. The
    GTF/t2g parity residual was accepted fail-closed. No truth, comparator,
    calibration, release or publication gate closes on this.
  - **Retained failure:** v3 evidence attempt 1 failed because `kallisto`
    and `bustools` were not on PATH; it is kept.
  - **Open:** for the HIV controls, the evidence helper's "largest
    non-target" is another HIV gene, so no non-HIV candidate was traced.
  - Previous text follows as history. Execute the outcome-ineligible
    ViralScan 2.2.0 versus v3
  diagnostic in `analysis/legacy_v2_v3/`. The identical-BUS arm is complete
  (44/44 valid rows) and the five-control input gate is closed. Fresh
  matched-FASTQ attempt-2 arrays `25331035` and `25331037` are **terminal with
  all ten tasks failed** across three independent causes: a nested-output
  validator defect that manufactured exit 65 on four genuinely successful v2
  rows, an unpopulated Zenodo viral-data cache that killed all five v3 rows
  before quantification, and one genuine v2 out-of-memory on `SRR6825024`.
  Attempt 3 is fully wired and tested but **not frozen**, because the cache pin
  it requires is blocked on `REF-11`; `LVC-13`–`LVC-14` remain unstarted and have
  no comparison tooling yet. A same-input head-to-head with evonk's released
  2.2.0 is tracked separately as `CMP-06` (WP6B).
  - 2026-10-02 (`v3/gov06`, e3a08b9).
    - **Tooling:**
      - highmem tier wired in, via `--highmem-task`;
      - explicit `-gtf` path, so the cache is optional;
      - GTF/t2g parity check;
      - LVC-12 exact-read evidence step;
      - `compare_legacy_v2_v3.py fresh-vs-archive`. The four attempt-2 v2 rows match the archive on all 300 records, within 1e-6.
    - **LVC-11:** the attempt-3 v2 highmem packet is submitted as job 25694919.
    - **LVC-12:** blocked on a user decision.
      - The parity check of the v2-arm GTF (195 files, 2,692 genes) against the frozen index t2g fails:
        - 89 GTF-only genes;
        - 200 non-Ensembl t2g genes not in the GTF, including all of HIV-1.
      - The v2 panel has no HIV GTF, so the archive could never call HIV.
      - The protocol 1.2.0 text is drafted, not applied, until that decision.
    - **LVC-12 unblocked (2026-10-02, `v3/gov06-hiv`, 325cc58).** The user decided to add an HIV GTF.
      - Source: evonk `Serratus_v2/gtf/Human_immunodeficiency_virus_NC_001802.gtf`. Its gene IDs and coordinates equal the index's IMMUNO_HIV1gp1–10.
      - Parity after the addition: 2,613 of 2,702 genes in the t2g. The residual (89 GTF-only, 190 t2g-only) is accepted fail-closed (`--accept-parity-residual`).
      - Protocol **1.2.0** applied, `outcome_triggered: true`; the arms now differ on HIV.
      - Packet `fresh_control_packet_attempt3b_v3/` (5 v3 rows) frozen and submitted.

`G0` passes when the governance scan is green, every public quantitative claim
is registered, all pre-v3 quantitative claims are rejected or historical, and
the inventory identifies enough retained BUS/reference/FASTQ material to rebuild
each eligible result.

Evidence to record:

```text
python3 scripts/check_data_governance.py
python3 scripts/validate_claim_registry.py --coverage
sha256sum analysis/v3_artifact_inventory.tsv claims/registry.json
```

## WP1 — Core software contracts

- [x] `SW-26` — (2026-10-09, found during PANEL-01 verification) native `--reference` must not use
  kb-python's current-directory `tmp`. The checkout's unrelated `tmp/` makes the documented tiny
  reference-build sequence fail before counting. Allocate a fresh per-build temporary workspace;
  preserve the existing directory. Implemented in shared `build_reference._run_kb_ref` for all three
  native kb-ref call sites. Regression covers successful/failed shell-outs, scratch cleanup, fresh
  retries and preservation of an existing `tmp/`; CLI/reference unit checks: 260 pass. Acceptance met:
  `TestDocumentedSequence` plus the anellovirus integration chain: 25 pass in this checkout.

Objective: close remaining correctness and workflow-consistency gaps before
scientific-scale execution. Estimated remaining effort: 4-7 engineering days.

### WP1A — Schemas, assignments, and reruns

- [x] `SW-01` — move/package all v3 schemas inside the installed `viralscan`
  distribution, load them with `importlib.resources`, and fail closed when a
  required schema is missing. Add wheel and sdist tests.
- [x] `SW-02` — enforce every public JSON/TSV/H5AD v3 schema at write and
  `validate-run` boundaries; remove generic silent skips. Three of the six
  shipped schemas had no reader at any boundary: `count_audit`,
  `reference_manifest`, and `evidence_manifest` shipped without ever validating
  a document. `h5ad_contract.json` is not a JSON Schema — it declares no
  keywords, so handing it to a validator would have accepted everything while
  looking like enforcement; `validate_json_schema` now refuses it by code
  (`not_a_json_schema`) and `_matrix_issues` reads `required_layers` and
  `required_uns` from it instead of from literals, closing a two-sources-of-truth
  gap that had left `quantification_unit` and `multimap_method` unchecked.

  Write boundaries raise (`SchemaContractError`) rather than returning issues;
  `validate-run` reports. The split is deliberate and is about authorship: at
  `validate-run` the artifact is input and a violation is a finding, while at a
  write boundary ViralScan is the author and a violation is a defect in this
  code, so publishing the file anyway would ship it under a schema it does not
  meet. Do not re-litigate this into a uniform policy.

  Note the h5ad count invariants were already fail-closed at *construction*
  (`multimapping.py:698` raises); `_matrix_issues` re-derives them from bytes on
  disk, which is a different guarantee — it catches a truncated write, a
  hand-edited file, or an artifact from another version rather than a compute
  bug. Two audit fields (`resolved_molecules`, `ignored_read_multiplicity`) and
  the contract's third invariant (unique mass equals audited unique molecules)
  had no reconstruction-side check at all; they do now.

  Silent skips removed: the `if manifest:` guard that let an empty run manifest
  pass schema validation and fingerprint checks, and the `ImportError` branch in
  `validate_json_schema` that turned a broken install of a hard dependency into
  a soft finding.
- [x] `SW-03` — **done 2026-10-10.** Add optional compressed molecule-assignment evidence containing
  CB, UMI, ECs, distinct genes, ambiguity class, method, weights, and exclusion
  reason without changing default matrix mass.
  - 2026-10-10 (wave 4): doc note: `docs/quickstart.md` now lists
    `--multimap-molecule-assignments` in the multimap paragraph and the "Useful flags" table.
- [x] `SW-04` — make `rerun-multimap` regenerate every method-dependent artifact
  in a new result tree: matrix/layers, count audit, summaries, evidence tiers,
  UMAPs, and host-response inputs.

  Audited 2026-07-28, one artifact at a time, because "regenerate everything" is
  not checkable without knowing which artifacts are actually method-dependent:

  - **matrix/layers** — regenerated. `_swap_multimap_layer` rewrites `X`, both
    compositional layers, and `uns["multimap_method"]`.
  - **count audit** — *not* stale, and does not need regenerating. Every field in
    `molecule_audit` is method-invariant: the molecule counts describe resolution
    rather than allocation, and `allocated_ambiguous_mass` is invariant because
    `host_conservative` divides `count / sum(cons_eligible)` across eligible genes
    rather than dropping mass — total allocated mass equals the ambiguous molecule
    count under every non-EM method, which `MoleculeAudit.validate` requires
    anyway. An earlier reading of this row assumed the swap corrupted the audit;
    it does not, and the fixture used to "reproduce" it was a state the pipeline
    cannot produce.
  - **summaries** — was broken, and not only on rerun. `multimap.py` opened
    `summary.txt` with mode `"w"` and wrote three totals; `detection.py` runs
    later in the DAG and opened the same path the same way, so those totals were
    destroyed on **every** run and never published. Nothing consumed them — no
    test, doc, or notebook. detection is now the sole writer and recomputes them
    from the H5AD, so they are both published and correct after a layer swap.
  - **evidence tiers, UMAPs, host-response inputs** — audited 2026-07-28 by a
    multi-agent review. Three further defects, one shared root cause: the
    command's invalidation list was incomplete, so `shutil.copytree` left old
    artifacts in a tree labelled with the new method.
    - `run_manifest.json`: **the reported defect was not real, and the fix for it
      was a regression.** The review claimed the manifest lives beside each
      sample's `config.yaml`, so the rewrite at `output_dir/` (the tree root) was
      a no-op. Running the pipeline end to end for `SW-10` showed the opposite —
      a completed run of `viralscan -o out` produces `out/run_manifest.json`
      alongside `out/<sample>/config.yaml`, exactly one manifest, at the root.
      The original code was correct. The per-sample rewrite shipped in `e6315cb`
      moved it to a path that never exists, so the manifest stopped being updated
      at all. Reverted 2026-07-29, with the layout now pinned by
      `test_run_manifest_is_at_the_tree_root`. The finding was accepted without
      being run; five of the review's findings were spot-checked and this was not
      one of them.
    - `plots/` was never cleared. Which viruses clear `detection_threshold` is
      method-dependent, and `generate_html_report` globs the directory, so a
      demoted virus's figure was re-embedded into a report whose own table no
      longer listed it. Only detection-owned patterns are cleared; `umap.py`
      writes into the same directory and its output is left alone.
    - `hostresponse/` was never cleared and `log/hostresponse.done` was never
      dropped, so it re-ran only if snakemake happened to judge it stale by
      mtime, and a virus falling below `MIN_VIRUS_CELLS` kept its old CSVs.

  Closed by `SW-05`.
- [x] `SW-05` — add an integration test proving no stale artifact survives a
  method change and the source result remains untouched.
  `tests/integration/test_rerun_no_stale_artifacts.py`, eight cases over the
  scenario that actually orphans files: a virus that clears the detection
  threshold under the source method and falls below it under the new one, so it
  is skipped rather than rewritten. Asserts the demoted virus leaves no plot or
  CSV, the surviving virus is regenerated rather than merely kept, provenance
  names the new method, artifacts owned by other rules survive the cleanup, every
  method-dependent sentinel is dropped, and the source tree is byte-unchanged.

  Both fixes were mutation-tested before the row was flipped: making
  `clear_stale_virus_outputs` a no-op fails
  `test_the_demoted_virus_leaves_no_hostresponse_csv_behind`, and restoring the
  old root-level manifest path fails `test_provenance_names_the_new_method`.

  Scope limit, stated in the module docstring: Snakemake is not invoked, so this
  proves nothing stale survives *when the rules re-run*, not that the DAG
  re-executes. That needs `SW-10`.

### WP1B — Workflow safety and architecture

- [x] `SW-06` — `doctor`, `validate-run`, fingerprints, resume/overwrite safety,
  and atomic manifests exist; finish whole-workflow staging and an atomic
  completion marker.
  - Done 2026-10-02 (`v3/cli`, deaf5c9). The user accepted that the
    completion marker serves as the whole-workflow staging.
    - `run_complete.json` holds the fingerprint, version and samples, plus the
      sha256 of `results/{viral_summary,virus_identity,multimap_evidence}.tsv`
      and `adata_multimap.h5ad` when present. It is written atomically after
      the sample loop and after rerun-multimap.
    - `prepare_output_directory` clears it on resume, overwrite and the
      rerun-multimap copy. rerun-programs and hostresponse re-stamp it.
    - `completion_marker: true` is added to the manifest after the
      fingerprint hash, so old `--resume` still matches.
      `validate_run` errors only for manifests that declare the marker; for
      old runs a missing marker is only a warning.
    - Schema `run_complete.schema.json` is in `REQUIRED_V3_SCHEMAS`.
    - Tests: 7 unit tests, plus an e2e check that deleting the marker fails
      validate-run.
- [x] `SW-07` — exact-fragment STAR filtering and mate synchronization exist;
  emit a reason for every retained/removed fragment plus lost-truth and
  host-virus-ambiguous boundary counts.
  - 2026-10-02 (`v3/tests-rel`, cbd5274):
    - `fragment_lineage.tsv.gz` now has one row per input fragment:
      `retained/host_unmapped` or `removed/host_mapped`. Memory is
      O(retained), by streaming input IDs against the retained set.
    - `lost_truth_counts()` is a Python hook (no CLI flag) that gives
      D15/D16 counts. Tests include a real STAR integration test.
    - Documented in `output_reference.md`.
    - Still open: the host–virus-ambiguous boundary count (deferred until
      defined), and the real lost-truth numbers, which need VAL-01.
  - 2026-10-10 (wave 4): **boundary count defined and done.** The truth contract has no `mixed` value in the `label`
    column (only `viral`/`host`); mixed host-virus fragments are the rows of `mixed_host_virus_manifest.tsv`. So
    `lost_truth_counts(..., mixed_manifest_tsv=)` adds `mixed_fragments` and `removed_host_virus_mixed` (the manifest's
    `read_id`s that the boundary removed; `None` without the manifest; reads the lineage never saw are not counted).
    No new vocabulary. Real numbers still wait on VAL-01, which generates the manifest. Needs the two keys added to the
    `output_reference.md` host-filter paragraph (wave-4 doc commit).
- [x] `SW-12` — range-check the EM parameters. `multimap_pseudocount` was
  guarded; `multimap_em_max_iter` and `multimap_em_tol` were not. A budget of
  zero makes `range(1, max_iter + 1)` empty, so `em_gene_abundances` and
  `em_cell_abundances` return their pre-loop seed weights while the H5AD records
  `multimap_method` as an EM method. Mass conservation still holds, so
  `MoleculeAudit.validate` cannot see it, and the only trace is a
  `converged: false` diagnostic that nothing reads. Both now raise in
  `RunConfig.from_snakemake_config`.
- [x] `SW-08` — remove unsafe kallisto CB-UMI-wide host filtering from the stable
  CLI because exact fragment identifiers are unavailable.
- [~] `SW-09` — **parsers done 2026-10-10 (`src/viralscan/cli/parsers.py`); handlers deliberately left in `menu.py`, see W3-11.** Split the oversized CLI into thin parsers plus importable service
  functions; convert Snakemake scripts to minimal wrappers without changing
  outputs.
- [x] `SW-10` — run one tiny paired-end fixture through documented CLI commands:
  preflight, reference, quantification, molecule allocation, cell calling,
  summaries, evidence, BAM/BLAST/plots/IGV, and `validate-run`.
  Executed for real on 2026-07-29 against `tests/data/evidence_tiny` (two read
  pairs, two reference sequences) using the tool binaries in
  `benchmark_runs/legacy_v2_v3/env_full/bin` (kallisto 0.50.1, bustools 0.43.2,
  kb_python 0.28.2, snakemake 8.20.5). `create_config` → `kb_count` → `analysis`
  → `multimap` → `detection` completes with exit 0, publishes all thirteen
  expected artifacts, and `validate-run` returns `ok: true` with zero issues —
  which exercises every `SW-02` schema check against a real artifact rather than
  a fixture. Codified as `tests/integration/test_tiny_end_to_end.py`, 17 cases,
  skipped when the binaries are absent.

  Two things the run established that no unit test could:
  - **`REF-11` is a hard blocker on the default path, demonstrated rather than
    inferred.** Without `-gtf` the run dies in `analysis`: the bundled panel is
    fetched from the unregistered Zenodo DOI. The fixture now ships its own
    minimal `viral.gtf`. When `REF-11` resolves, add a variant that drops `-gtf`.
  - **`run_manifest.json` lives at the tree root**, beside the per-sample
    directories rather than inside one. See `SW-04` below.

  Not covered: the evidence/BAM/BLAST/IGV leg, which needs `blastn`,
  `makeblastdb`, and `minimap2` — none present in this environment. Those legs
  are exercised by `test_exact_lineage.py` and `test_evidence_chain.py`. The row
  stays `[~]` until they run in one sequence.
  - Done 2026-10-02 (`v3/tests-rel`, 8149423). `TestDocumentedSequence` runs
    the tiny fixture end to end:
    1. `doctor --profile full`.
    2. The documented `-ref -fasta -gtf` build (`kb ref`). `build-ref` needs
       a network download, so it isn't used. The fixture needs a distinct
       viral record name and `transcript` rows.
    3. A run with visuals (plots and HTML).
    4. `--cell-calling external`.
    5. `viralscan evidence --blast` on the real run (all 15 artifacts).
    6. `validate-run` before and after evidence.

    Correction to the text above: env_full now has blastn, makeblastdb,
    minimap2, samtools, Rscript, STAR and cd-hit-est. The no-`-gtf` variant
    still waits on REF-11.
  - 2026-10-10 (wave 4, `CAT-06`). `TestManifestRunWithoutGtf` runs the tiny
    fixture with no `-gtf` against an index that has a build manifest; the
    integration run passed (19 cases). A manifest-less run without `-gtf` (the
    true default path) still needs the bundled panel, so `REF-11` still blocks that.
- [x] `SW-11` — make production cell calling fail closed: caller exceptions,
  zero-match external lists, invalid barcode geometry, and canonical collisions
  must never silently turn every barcode into a cell; `none` remains explicit.
  `detection.py` caught every exception and continued with `called_mask=None`,
  which `compute_stats` expands to all-ones, so a caller failure silently
  replaced the called-cell denominator with every barcode and still labelled the
  result a called-cell rate. All four paths now raise `CellCallingError`;
  `--cell-calling none` remains the explicit way to report over all barcodes.
  Completed 2026-07-27 by the review follow-ups: `emptydrops_seed` and
  `emptydrops_niters` are declared configuration with CLI flags rather than
  function-signature defaults, so the protocol's frozen `seeds.cell_calling`
  actually reaches `set.seed()` in `emptydrops.R`; and `Rscript` is preflighted
  whenever the resolved caller is `emptydrops`, because failing closed at the end
  of a multi-hour run for a knowable reason is the wrong place to fail.

- [x] `SW-13` — resolve `snakemake --config` wire values in Python, not in shell
  templates (2026-09-29).
  - Snakemake parses an empty `whitelist=` as `None`, and `{config[whitelist]:q}`
    rendered it as the literal `None`. **Every run without `-w` therefore called
    `kb count … -w None`, which kb reads as "bypass barcode error correction".**
    Seen in the SRR12682296 `kb_info.json`: bus → sort → inspect → count, no
    `bustools correct`, 793,308 barcodes. All reference-strategy and max-panel
    runs since then are uncorrected and need rerunning for absolute numbers.
  - Booleans: Snakemake converts only `True`/`False`, so `gene_programs=false`
    arrived as the non-empty string `"false"` and the rule was always planned.
  - Now a rule `params` lambda resolves the whitelist, and the gene_programs gate
    compares `str(...).lower() == "true"`.
  - Regression tests in `tests/test_snakefile_dag.py`: a template check, plus
    `snakemake -n -p` renders `WL=` when unset and `WL=/fake/wl.txt` when set.
    Both fail on the old Snakefile.
  - Root-cause removal is `MECH-C` (single Run Config writer).
  - Follow-up: `tests/data/evidence_tiny/R1.fastq` carried off-list barcodes
    (`AAAA…`/`GGGG…`), and the tiny e2e run had passed only because correction
    was bypassed. With correction on, kb drops them and cell calling fails
    closed on 0 barcodes. The fixture now uses two v3 on-list barcodes.
- [x] `SW-20` — (2026-09-29) the other half of `SW-13`. ViralScan's own matrix still skipped
  barcode correction when no `-w` is given.
  - `multimap.prepare_resolved_bus` sorts the **raw** `output.bus` and runs
    `bustools correct` only for a user whitelist. The multimap-derived X, which
    feeds detection, is therefore still uncorrected, while kb's
    `counts_unfiltered` is now corrected.
  - Use kb's corrected and sorted `output.unfiltered.bus`, or its copied on-list
    (e.g. `10x_version3_whitelist.txt`), when `-w` is absent.
  - Count the dropped off-list records in the audit.
  - Test with a fixture barcode one mismatch from the on-list.
  - Fixed. `multimap.select_bus_input` picks kb's `output.unfiltered.bus`
    (`KbCountOutputs.kb_corrected_bus`) when no user on-list is given, and warns
    when neither kb nor the user corrected. Checked end to end on the tiny
    fixture plus one read with barcode `GAACCCAAGAAACACT`, which is one
    mismatch from on-list `AAACCCAAGAAACACT` and has a unique neighbour:
    - fixed code: the read is corrected, 3 input molecules, 2 viral molecules;
    - previous code: the raw barcode is dropped **before** the audit counts it,
      2 input molecules, 1 viral molecule, and the loss is invisible.
  - Still open, under `MECH-F`: count off-list drops in the audit rather than
    before it.
- [x] `SW-21` — **closed 2026-10-03** (user decision: bypass). See DEF-02 (c).
  Original text: (found 2026-09-29, F-018) since `SW-13`, a Drop-seq run
  without `-w` is pre-filtered to a data-derived cell list.
  - Cause: Drop-seq has no official on-list. `kb count` therefore runs
    `bustools allowlist`, a knee on the data, and corrects against it. Every
    other barcode is discarded, both in kb's matrix and in the corrected BUS
    that `SW-20` feeds to multimap.
  - Measured on HSV-1 SRR8315713 (job 25666447, v1 index):
    - kept barcodes: 1,809,992 → 5,468;
    - HSV-1 molecules: 32,404 → 23,170 (−28 %);
    - infected barcodes: 10,227 → 1,127.
  - The 5,468 barcodes become the "all barcodes" denominator, and ViralScan's
    own cell calling runs after this implicit one.
  - 10x runs behave as intended. EBV SRR12682296 (10xv2, official on-list)
    keeps 885,933 of 906,202 molecules (−2.2 %), and the `summary.txt` headline
    now matches, at 885,939.
  - Decision needed, under `MECH-D`, the chemistry module: a technology with no
    official on-list either keeps no correction (pass kb its bypass value) or
    keeps kb's allowlist and documents the pre-filter. This is not changed in
    the Snakefile until then.
- [x] `SW-22` — (found 2026-09-29) `rerun-multimap` does not skip `kb_count`.
  - **Done 2026-10-04.** `_run_rerun_multimap` restores config.yaml's mtime after
    the rewrite, so the sentinels alone choose what reruns. Pinned by
    `test_rewritten_config_keeps_its_mtime_so_kb_count_is_not_rerun`. The
    MECH-C refactor (delete the k=v round trip) stays open as architecture,
    not a G1 defect.
  **Closes with `MECH-C` (reconciled 2026-10-04).** It is a G1 defect candidate (package plan
  M1.4).
  - Cause: `_run_rerun_multimap` rewrites the copy's `config.yaml`, which is an
    input of `kb_count`. Under snakemake's default mtime trigger, `kb_count`
    (kallisto on the FASTQs), `analysis` and `multimap` all rerun, which also
    undoes the fast layer swap.
  - Shown by dry-run on a skeleton copy of `combined_corrected/SRR8315713`.
    Without the config rewrite, only detection, umap and all are scheduled.
  - This predates MECH-A. A fix (for example keeping the config's mtime, or
    moving the multimap parameters out of `kb_count`'s inputs) belongs with
    `MECH-C`, the single Run Config writer.
- [x] `SW-14` — (2026-09-29) the Snakemake invocation in `menu.main` and `_run_rerun_multimap`, now one helper, `menu._snakemake_run_command`, covered by `tests/test_cli.py::TestSnakemakeRunCommand`:
  - put the `all` target before `--quiet`, because snakemake 9 lets
    `--quiet [...]` consume the target;
  - drop the unconditional `--use-conda`, since conda is not required at run
    time.
- [x] `SW-15` — (2026-09-29) STAR host-filter geometry. Done in `host_filter.starsolo_barcode_args` (pure) and `_plain_whitelist`, with tests in `TestStarsoloBarcodeArgs`. Merging with the benchmark builder is left to `MECH-D`:
  - drop the invalid `--outSAMflag None`;
  - add `--soloBarcodeReadLength 0`, because 10x 5′ with a 150 bp R1 aborted
    with "barcode length 150 ≠ 28" on SRR20710647's library;
  - decompress a gzipped whitelist;
  - share `reference_strategy.starsolo_geometry_args`.
- [x] `SW-16` — (2026-09-29) the `summary.txt` headline reported "Viral molecules … 0; Cells
  with viral reads 0/793308" beside 906,202 EBV molecules in
  `viral_summary.tsv`. `_headline_totals` filters virus *names* against
  `var_names`. Fixed: the call site passes the detected gene IDs, and
  `_headline_totals` raises if handed IDs that match no `var_name`. A new tiny-e2e
  assertion (`test_summary_headline_matches_viral_summary_total`) compares the
  headline with the `viral_summary.tsv` total; it fails on the old code.
- [x] `SW-17` — (2026-09-29) the `evidence --virus hhv6a/hhv6b/hhv8/kshv` selectors raised. `VIRUS_ALIASES` now lists every label each virus resolves to (legacy, catalogue, RefSeq) and matches case-insensitively, until `MECH-A` replaces it with taxid lookup. The cause was:
  `VIRUS_ALIASES` overwrites `hhv8`, and there is a `6B`/`6b` case mismatch.
- [x] `SW-18` — (2026-09-29) `anellovirus.gtf_text_for` returned `''` for uncatalogued
  accessions when `fasta_texts` is absent, so genomes silently lose their GTF.
  Its docstring promises a 1 bp exon.
  - Fixed: it now raises `ValueError` naming the unsized accessions.
  - It also iterated `accessions` twice, so a generator argument silently lost
    every placeholder. The argument is now materialised once.
  - Tests are in `tests/test_anellovirus_reference.py`.
- [x] `SW-19` — (2026-09-29, reproduced then fixed) `rerun-multimap` rewrote the *source* run's h5ad, because
  `from_yaml` runs before the `output` rewrite (menu.py ~595 vs 616). Write an
  e2e test through `_run_rerun_multimap` first, and fix only if it reproduces.
  - **Reproduced.** `TestRerunLeavesSourceUntouched` drives the real function
    with snakemake stubbed, and the source h5ad's sha256 changed.
  - Fix, part 1: the fast path resolves the h5ad under the copy's sample
    directory.
  - Fix, part 2: `cfg["output"]` keeps its trailing separator. Without it the
    Snakefile's `f"{config['output']}log/…"` paths would have been
    `…SAMPLElog/`.
  - Fix, part 3: the `--config` args come from
    `RunConfig.to_snakemake_config_args()` instead of a private `_arg_val`.

`G1` passes when all count invariants, schema checks, safety scenarios, rerun
consistency, and the full tiny workflow are green. No known correctness or data-
loss defect may remain.

Required gate:

```text
NUMBA_CACHE_DIR=/tmp/viralscan-numba-cache PYTHONPATH=src python3 -m pytest tests/ -q
python3 -m ruff check .
python3 -m ruff format --check .
PYTHONPATH=src python3 -m pytest -m "integration and not network" -q
```

### WP1D — Mechanism consolidation (architecture review, new 2026-09-29)

A read-only review found three walks of the code whose scientific decisions are
re-derived at every call site from string formats:

- virus identity and detection;
- the counting path and Run Config;
- reference construction.

The helpers pass their tests while the call sites produce wrong outputs.
`SW-13`…`SW-19` are the defects it observed. The user chose `MECH-A` first,
and settled three decisions in review:

- **A "virus" row** is one NCBI organism with rollup columns. The key is the
  taxid; segmented viruses are keyed by (species, strain or isolate).
- **An uncatalogued indexed gene** counts as viral if it is in `--gtf`, with a
  warning.
- **The catalogue** is merged into the packaged `virus_catalog.tsv`.

The review report (candidates A–F, the defect table, the builder × gate
matrix and MECH-A progress) is published as a private artifact:
https://claude.ai/artifact/QXWXwBk3BiJNSBnioYiUKH (2026-09-29).

- [x] `MECH-A` — (2026-09-30) per-Run **Virus Identity table**, `src/viralscan/virus_identity.py`,
  built once by the `analysis` rule. It maps gene_id → genome accession (t2g
  column 5) → catalogue row → viral status, virus key, name, family, sibling
  group and risk class.
  - [x] Step 1 (2026-09-29): catalogue merge. `virus_catalog.tsv` now holds
    4,128 accessions (504 species, every row with a taxid), up from 2,249.
    `build_virus_catalog.py` adds `taxid`, `organism`, `strain` and
    `serotype`, and carries forward the overlay columns `common_name`,
    `sibling_group`, `role` and `panel`. `extras/seed_catalogue_overlays.py`
    seeded them:
    - `common_name` on 303 rows, from the legacy prefix map;
    - `sibling_group` by taxid: HSV, HHV-6, and HHV-4 (EBV-1/2, F-017);
    - `risk_class=eve` on 3,448 rows;
    - `role=decoy` on 9 rows;
    - `panel`: `shipped` on 2,345 rows, `max` on 1,783.
  - Decision: CAT-31 now reconciles only `panel=shipped` rows.
    `catalogue_detection_targets(panel="shipped")` is the default, and
    `panel=None` checks every row. The max-panel rows are there so the identity
    table can name them; they are not claims of the shipped panel.
  - Fix: the builder keeps the newest version of each accession. The pre-merge
    catalogue had `NC_006312.1`, but the bundled GTF uses `.2`.
  - Finding for step 2: strain text is not a usable key. Two of the eight PR8
    segments spell the strain "A/Puerto Rico/8/1934(H1N1)". The isolate-level
    taxid (211044) is the same on all eight, but 24 Influenza A rows carry only
    the species taxid 11320. The segmented key therefore needs a normalised
    strain, with the `(HxNy)` suffix stripped.
  - [x] Step 2 (2026-09-29): `src/viralscan/virus_identity.py` and
    `tests/test_virus_identity.py` (32 tests).
    - `build_identity_table(t2g, gtf_gene_ids)` gives one row per indexed gene,
      with status catalogued / uncatalogued / host / legacy_prefix.
    - It fails when no gene is viral, and warns with a count for uncatalogued
      genes.
    - Golden tests pass on `t2g_v2`, the final `panel.t2g` and `t2g_max`, and
      skip where those files are absent:
      - every viral gene is catalogued, and every host gene is ENSG (41,145);
      - EBV-1 and EBV-2 are two viruses in sibling group HHV-4;
      - HPV45's two records form one virus;
      - the 18 SARS-CoV-2 genomes form one virus;
      - PR8 is one virus with 8 segments;
      - no virus holds a segment twice;
      - no two viruses share a display name.
    - The virus counts are 113 (v2), 236 (final) and 421 (max).
  - Decision, revising the "segmented = (species, strain)" key: segments are
    keyed by **taxid**. A taxid is split by normalised strain only when it holds
    the same segment twice.
    - Reason: strain text is unreliable inside one RefSeq set. It is empty on
      some segments ("Hantavirus Z10" M, "Pichinde" S) and spelled "…/1/96" or
      "…/1/1996" on others (goose/Guangdong, taxid 93838).
    - Only the species-level taxid 11320 ("Influenza A virus", 4 isolates)
      splits.
  - Decision (user, grill Q8b = A2, 2026-09-29): anelloviruses are keyed by
    **genus** (`ANDET-05`) and not by taxid. Genomes with no assigned genus are
    named "Anelloviridae (genus unassigned)", because a bare family row beside
    genus rows reads as their total.
    - 4 of 120 anellovirus taxids span several genera. They are the catch-all
      bins that hold most rows:
      - 2055263 "Anelloviridae sp." (Alpha, Beta, Gamma);
      - 68887 "Torque teno virus" (Alpha, Beta, Gamma, Samek);
      - 93678 "TTV-like mini virus" (Beta, Gamma, Het);
      - 432261 "Torque teno midi virus" (Gamma, Mem, Samek).
    - Catalogued anellovirus rows carry no `common_name`. (An earlier count of
      "15" treated the table's genus label and the same genus from the
      catalogue as different labels.)
    - The genus comes from the anellovirus accession table. Otherwise it falls
      back to the catalogue genus, when that is an ICTV genus name (one word
      ending in "virus"), then to "Anelloviridae". This keeps the non-genus
      lineage token "Small anellovirus" out of the keys.
  - Step 3 caveat for step 4: a run directory made before 35940ec has no
    `results/virus_identity.tsv`. Snakemake does not rebuild it, because no rule
    consumes it yet (checked by dry-run on a skeleton copy of
    `combined_corrected/SRR8315713`: only detection, umap and all are
    scheduled). Once step 4 makes it an input, `rerun-multimap` and resume must
    build it in process for old runs, or the `analysis` rule reruns.
  - Guard (user, grill Q8 = B): keep the structural guard below and add an
    index build manifest that records the host and viral gene sets. The
    manifest overrides `--gtf`, and a run is refused when `--gtf` would mark a
    manifest-host gene as viral. The manifest, fixing `reference_strategy.py`'s
    `-gtf`, and detection reading the table are all step 4 work.
    Structural guard: a `--gtf` gene whose t2g
    column 5 is a transcript of the index is host. This is the host-cDNA row
    shape of every combined index measured: 465,769 host rows, 0 viral rows.
    - Reason: `reference_strategy.py:704` passes the combined GRCh38+viral GTF
      as `-gtf`. Today's `analysis.py` makes every host gene in that path viral.
      Our SLURM runs passed viral-only GTFs and are unaffected.
  - [x] Step 3 (2026-09-29): the `analysis` rule writes
    `results/virus_identity.tsv` as a declared output, next to
    `log/analysis.txt`, which is unchanged.
    - `analysis.write_identity_table` builds the table from `config.transcripts`.
    - A run whose index has no viral gene now stops at this step.
    - The column reference is in `docs/output_reference.md`.
    - The tiny e2e (20/20) checks that the fixture's 3-column t2g takes the
      legacy fallback.
    - No consumer reads the table yet.
  - [x] Step 4 (2026-09-30, branch `v3/mech-a-consumers`: 63c3de1 phase 1,
    160a4db phase 2): the consumers read the table.
    - `multimap.py` takes its viral/host partition from the `viral` column.
    - `detection.py` groups by `virus_key`, names by `virus_name`, and takes
      the sibling note from `sibling_group`. The weaker members of a group
      are flagged against its dominant member, which generalises the old
      pairs. `eve_risk` comes from `risk_class`.
    - `umap`, `hostresponse`, `evidence` (`--virus` resolves by key, taxid,
      handle, name, organism or family) and `gene_programs` (by `virus_key`)
      also read the table.
    - The consumer rules declare `virus_identity.tsv` as an input.
    - `SIBLING_VIRUS_PAIRS`, `EVE_RISK_GENERA` and `merged_name_map` remain
      only as the fallback for runs with no table. CONTEXT.md has the new
      terms.
    - Tests: `test_mech_a_consumers.py` and `test_mech_a_presentation.py`.
      They cover every retired sibling pair and EVE genus, per-gene →
      per-group conservation including `genus:` keys, the partition,
      selectors and programme lookup. Main tree: 1,519 passed, 0 skipped.
    - Decision (EVE, user-confirmed 2026-09-30): a virus is flagged when any
      of its genes has `risk_class=eve`. An empty `risk_class` on a catalogued virus means no
      EVE risk (680 rows are empty; flagging them would flag EBV).
      Uncatalogued and legacy viruses fall back to the genus-name test, with
      a warning.
    - Count parity was checked against 815838a on copies of the EBV, HSV-1,
      HHV-6B and covid x213 runs (`viralscan_work/parity_mechA/`).
      `count_audit`, `found_genes`, `analysis.txt` and `positive_control` are
      identical in all four; every multimap layer is identical for EBV,
      HSV-1 and HHV-6B.
    - Explained naming diffs:
      - `Anelloviridae` → "Anelloviridae (genus unassigned)", 274 genes;
      - "Torque teno virus" and raw `D1P6x_gpN` IDs → `Alphatorquevirus`;
      - covid `VARVgp184`: VZV → Variola virus, which is correct (1 UMI);
      - `n_viral_accessions_in_reference` 9193 → 5093, because only indexed
        genes are counted.
    - Name checks pass: EBV "Epstein-Barr virus"; HHV-2 `possible_em_bleed`
      on HSV-1 at 665:1; the EBV-1/EBV-2 split is unchanged.
  - [x] Step 4a (2026-09-30) — **structural-guard false positive** (covid
    x213): gene `HUM_HERP6B_DR1` has t2g column 5 equal to its own
    transcript ID. The guard therefore marked a viral gene as host: 1 UMI in
    1 cell, `counts_unique_viral` −1.
    - t2g shape, measured: the covid index has 97 **self-named** rows
      (transcript = gene = column 5), all `HUM_HERP6B_*`. They are VIRTUS-style
      records, one per gene, so DR1 was only the one with a UMI. Host cDNA
      rows that point at a transcript never have transcript = gene (`ENST` ≠
      `ENSG`): 0 of 465,769 host rows in each of the five stored combined
      t2g files (covid, v1 full panel ×2, final 2,343, max).
    - Fix (`virus_identity.read_t2g`): a self-named row names no genome. Its
      accession is `""`, and it is not structural-host. `_resolve_genes`
      names a GTF gene with no accession by the legacy prefix map (status
      `legacy_prefix`, key `name:<virus>`), so the 97 genes are one virus,
      "Human herpesvirus 6b", not 97 accession-keyed ones. The exemption
      does not use column 3: 3,037 viral rows of the final panel have it empty.
    - Covid parity (`viralscan_work/parity_mechA/covid/final4a/`, job
      25684236): the 97 genes form one viral row, and HHV-6B gets its 1 UMI
      back. Every multimap layer, including `counts_unique_viral`, is now
      identical to 815838a. The HHV-6 note does not fire: HHV-6 has 38 UMI against HHV-6B's
      1, which is under the 50:1 threshold, and legacy-prefix rows carry no
      `sibling_group`. Detection warns that `eve_risk` comes from the genus
      fallback, as designed for uncatalogued viruses.
    - Tests: two in `tests/test_virus_identity.py` (a self-named row, and
      one prefix-named virus next to an `ENST`/`ENSG` self row that stays
      host). Full suite: 1,527 passed.
  - [x] Step 4b (2026-09-30) — **catalogue row for NC_000898.1** (HHV-6B,
    taxid 32604, `common_name` "Human herpesvirus 6b", `sibling_group`
    HHV-6). The catalogue held HHV-6B only as AF157706.1. On the stored
    SRR20710641 index HHV-6B therefore read `NC_000898.1`, and the HHV-6A/6B
    note did not fire.
    - `panel=legacy`, a new value, not `shipped` as first written here.
      `NC_000898.1` is in neither built panel. The shipped `viral.fa` has
      HHV-6B only as AF157706.1, and max-panel dedup dropped it as an
      `exact_hash` duplicate of AF157706.1 (`dedup.tsv`), the CAT-05 pattern.
      Only pre-v3 stored indexes (VIRTUS2-sourced) carry it. A `shipped` or
      `max` value would claim a panel that does not index it, and a
      reconciliation of that panel would report it as catalogued but not
      indexed. Reconciliation skips any other `panel` value. The identity
      table does not filter by `panel`, so naming is unaffected. It is the
      only catalogue row missing from its panel's FASTA.
    - Checked on a copy of the stored run with the working tree
      (`viralscan_work/parity_mechA/hhv6b/final/`, job 25683797). It now
      reads "Human herpesvirus 6b" (5,596 UMI), and the note fires on
      "Human herpesvirus 6" at 323:1.
  - [x] Step 4c (2026-09-30): **in-process backfill for old run dirs**.
    - `virus_identity.backfill_identity_table(config)` builds the table from
      the run's own `log/analysis.txt`, not a fresh GTF glob, so catalogue
      or GTF changes since the run cannot move its partition. The file takes
      `analysis.txt`'s mtime.
    - `rerun-multimap` calls it after rewriting the config, and so does
      `--resume` before snakemake starts. `write_identity_table` moved from
      `scripts/analysis.py` into `virus_identity.py` (re-exported), so
      `menu.py` can call it without the script's logging setup.
    - Tests in `tests/test_identity_backfill.py`, including a byte-identical
      match with the table the `analysis` rule writes.
    - Dry run on a timestamp skeleton of the stored SRR8315713 run, with
      `detection.done` and `umap.done` removed (what `rerun-multimap`
      does): without backfill, analysis + multimap + detection + umap are
      scheduled; with it, only detection + umap.
    - A complete old run with every target present schedules nothing either
      way, because snakemake does not rebuild a missing intermediate for
      up-to-date targets.
  - MECH-A stays `[~]` until 4a and 4b land and the HHV-6B name check passes
    on the stored index. All three were met on 2026-09-30, and MECH-A is
    closed.
  - It replaces 7 prefix-matching call sites, the GTF-only viral/host
    partition, the name-keyed `SIBLING_VIRUS_PAIRS`, and the substring EVE test.
  - Observed failures it fixes:
    - the HHV-6 sibling note is empty at 297:1 (SRR20710641);
    - one genome appears under two names (`TTVgp1` vs `NC_002076.2_gene1`);
    - the 1,912 max-panel accessions would report as raw IDs;
    - EBV-2 bleed shows as 80 unflagged per-gene rows (F-017).
- [x] `EXPL-R2.10`: exploratory reruns on the current reference (grill
  R2.10 = D, 2026-09-29). These do not tune anything (R2.2, R2.4). Submitted
  on 2026-09-29:
  - 25672273 `scripts/slurm_quant_max_corrected.sh`: F-017 redone on the max
    panel with barcode correction on, for EBV and HSV-1. Output goes to
    `viral_panel_max_2026-09-28/runs/combined_max_corrected/`.
  - 25672274 `scripts/slurm_quant_hhv6b_5p.sh`: HHV-6B SRR20710641 as 5′
    `-x 10xv2` on the v1 index. Its R1 has the TSO at base 27, so the barcode
    is 16 bp and the UMI 10 bp. Output goes to
    `ebv_latest_ref_2026-09-27/runs/combined_corrected_5p/`.
  - 25672275 `scripts/slurm_quant_covid_exploratory.sh`: covid x213 and x216
    with current code, the published index, and the CellRanger whitelist. R1 is
    28 bp with no TSO, so `-x 10xv3` geometry is correct. Barcode correction
    was already on in the published run, and there is no read filter yet, so
    the F-019 poly-G signal is expected to come back. Output goes to
    `covid_viralscan/results_v3_exploratory/`.
  - 25672276 `scripts/slurm_strand_test.sh`: 4M read pairs of each 5′ library
    (HHV-6B, x213, x216), run through `kb count --strand`
    forward/reverse/unstranded. Output goes to
    `viralscan_work/strand_test/*/strand_summary.tsv`.
  - Results (2026-09-30):
    - **F-017 is confirmed with barcode correction on:** EBV-2 takes 16.1 %.
    - **HSV-1 max panel:** +0.05 %.
    - **HHV-6B as 10xv2:** 5,596 molecules in 3,556 cells. It still reports as
      `NC_000898.1_gene1`, the MECH-A naming defect.
    - **Covid on current code:** 1.08M / 1.50M Alphatorquevirus, which
      reproduces F-019.
    - **F-020, the strand default:** kallisto's forward default keeps only
      6.5–8.9 % of reads on every 5′ library. `--strand` (DEF-02) is now a
      priority.
    - **Tonsil read origin:** reproduces the F-019 artefact, with 3.02M reads
      and 0 clean full-length viral reads.
- [x] `EXPL-TTMDV` — **Result 2026-10-03: artefact** (F-020 "EXPL-TTMDV
  result"). The reproduction gate passed exactly (x213 55,046 / 24,658;
  x216 7,391 / 1,829).
  - 0 clean, full-length, viral-best reads in either sample. 100 % carry
    homopolymers and 96 % carry the 10x TSO: they are TSO–oligo-dT–poly-G
    concatemers.
  - They hit a 26-nt poly-A tail after `AATAAA` in the covid index's
    AB303552.1/AB303557.1, which is the F-021 sink class. cat42b masks it.
  - gget blast: 0 of 14 reads hit an anellovirus.
  - It says nothing about donor TTV carriage, which is commensal and
    expected. My earlier "no homopolymer over 12 nt" was measured on the
    masked cat42b FASTA, and is corrected.
  - Original row: (2026-10-03, user-requested) read-level check of the
  reverse-strand Gammatorquevirus signal in covid x213/x216 (F-020 update):
  55,046 + 24,658 unique molecules on AB303552.1/AB303557.1 under
  `--strand reverse`. Exploratory; it tunes nothing. Anelloviruses are
  commensal (F-022), so the check is neutral on the prior.
  - Method, reused from F-019 (`viralscan_work/f005_readorigin/`), outputs in
    `viralscan_work/ttmdv_readcheck/<sample>/`:
    1. `kallisto bus -n --rf-stranded` (0.51.1, same index, `-x 10xv3`).
    2. Capture the anellovirus ECs, then correct.
    3. Extract the exact R2 reads.
    4. Competitive `minimap2 -ax sr` against GRCh38 + the covid panel.
    5. `classify`, then a spliced alignment to the TTMDV genomes.
    6. gget blast (login node) of the per-sample consensus and about 20
       position-stratified clean reads.
  - **Reproduction gate (before classifying):** the captured molecules on
    AB303552.1/AB303557.1 must match the reverse run (x213: 55,046 /
    24,658). Stop if they do not.
  - **Pre-registered reading:**
    - **Real TTMDV:** most molecules clean (no homopolymer ≥15, no TSO, not
      low entropy), full-length and viral-best; reads spread across the
      transcribed region (not one hotspot); blast top hits are anelloviruses.
      Variants recurring across independent CB-UMIs point to a donor strain.
    - **Artefact:** most molecules low-complexity, TSO or host-best; a single
      hotspot, especially in the GC-rich non-coding region; blast hits human
      or vector.
    - A consensus identical to the reference at every covered site is a
      contamination/reference signature.
    - Not evidence: strand (guaranteed sense by `--rf-stranded`), and
      identity ≥ 0.9 to AB303552 (donor strains diverge).
  - **Gate run 1 failed (job 25695673, 2026-10-03).** Re-counted vs run,
    in unique molecules:
    - x213: 57,889 vs 55,046 (+5.2 %) and 25,109 vs 24,658 (+1.8 %);
    - x216: 7,947 vs 7,391 (+7.5 %) and 1,956 vs 1,829 (+6.9 %).

    Always high, which fits a gate defect: the capture kept only
    target-EC records, so a molecule with other-EC reads looked
    target-unique. The tolerance is unchanged (2 %). Run 2 (job 25695740)
    gates on every record of each target-touching UMI (`bustools capture
    -u`). If it still fails, stop.
    - Descriptive only: identity distribution; called vs empty-droplet
      barcodes (ambient); shared CB-UMIs between x213 and x216. Donor
      identity of x213/x216 is unknown (batch 1/2), so a cross-sample
      consensus comparison is descriptive.
- [x] `EXPL-HPV16` — (2026-10-01) HPV16 positive control on GSE189670 (Bedard et al.,
  Nat Commun 2023, PMID 37031202), user-requested 2026-09-30. Isogenic NIKS
  keratinocyte rafts, 10x 3′ v3: SRR19537341 (GSM5705760, "HPV16 infected
  keratinocytes") vs SRR19537339 (GSM5705759, "Normal keratinocytes"), one
  run each. Shipped 2,343-genome panel (`viral_ref_final/build/panel.idx`),
  `-x 10xv3`. Exploratory; it tunes nothing.
  - Truth is qualitative only: the authors see HPV16 early and late genes in
    the HPV16 rafts and none in the parental rafts, with no per-cell % given.
    Pass = HPV16 present in the HPV16 sample, and ~0 in the normal one.
  - ENA serves R2 only for these runs, so reads come from SRA with
    `fasterq-dump --include-technical --split-files`.
  - Gene identity: the bundled GTF's `gene` attribute maps gp1=E6, gp2=E7,
    gp3=E1, gp4=E2, gp5=E1^E4, gp6=E5, gp7=L2, gp8=L1. The genes overlap,
    and 3′ reads pile up at the early and late polyA sites, so per-gene
    counts reflect polyA position, not ORF identity.
  - Scoping of the other user-supplied datasets (GSE164690, GSE208653,
    CELLxGENE HPV/CMV, E-MTAB-13687 tonsil):
    `viralscan_work/dataset_scoping_2026-09-30/`.
  - Result (2026-10-01, job 25684772, HEAD f40f72a): **HPV16 PASS on
    presence.** SRR19537341 (HPV16) has 7,725 molecules in 2,442 called
    cells. SRR19537339 (normal) has 8 molecules in 6 cells, about 970:1.
    - Per gene, HPV16 sample, as unique + allocated:
      E5 2,374; E1^E4 1,501 (1,478 allocated); E7 1,440; E1 1,100 (all
      allocated); E2 696; E6 451; L1 112; L2 51. Early genes dominate, which
      matches the authors' "mostly early genes". E5 sits next to the early
      polyA site, as expected for 3′ reads.
    - Both array tasks were marked FAILED only because the script's final
      `test -s $OUT/results/...` was wrong. viralscan nests results under
      `$OUT/<sample>/`. Snakemake completed every rule.
    - The cell-calling knee is not usable: total ≥10 gives 60,475 /
      1,569,733 (HPV16) and 86,135 / 1,292,165 (normal) barcodes. This is
      `SW-23`.
    - Indicative per-cell result, cut at the barcode-rank steepest descent
      (total ≥952 / ≥739, computed outside viralscan): 1,869 / 8,023 HPV16
      cells are HPV16+ (23.3 %), holding 7,034 of 7,725 molecules. The
      normal raft has 0 / 11,045, and all 8 of its molecules sit in empty
      droplets.
    - **emptyDrops rerun** (2026-10-01, job 25689583, detection only on
      copies in `runs_emptydrops/`; FDR 0.01, lower 100, seed 100):
      | | HPV16 raft | normal raft |
      |---|---:|---:|
      | called cells | 16,079 (knee 9,895, inflection 1,052) | 11,884 |
      | HPV16+ called | **2,012 (12.5 %)** | **0** |
      | HPV16+ comparable (≥200 host UMI) | 1,927 / 9,086 (21.2 %) | 0 / 9,793 |
      | Gammatorquevirus+ called | 2,552 (15.9 %) | 3,314 (27.9 %) |
      These are the reportable numbers. The rerun reused the run's stored
      `virus_identity.tsv`, so it still shows the old label "16,18"; new runs
      get the fixed name. The anellovirus row is the F-021 artefact, now
      27.9 % of the normal raft's called cells.
    - **Specificity failure, anellovirus.** Gammatorquevirus has 5,688
      molecules (HPV16) and 10,037 (normal) in an anellovirus-free cell
      line. In the normal sample all 10,037 sit on one placeholder gene,
      `KP343824.1_gene1` (an "UNVERIFIED" TTV isolate S57 record whose
      first 29 nt are poly-T; the builder masks only runs ≥31). This is the
      same single-bucket sink as F-005/F-019. Read-level check: of the
      first 20M R2 reads, 27,098 carry a 31-mer from that poly-T head. They
      are host mRNA ends running into **poly-A tails**, not the poly-G reads
      of F-019. See F-021 and `CAT-42`. KP343824.1 is also in the max
      panel, yet the 3′ EBV (10x v2) and HSV-1 (Drop-seq) runs put 0 on it,
      so the capture depends on the library, not on 3′ chemistry alone.
    - The label was wrong: HPV16 reported as "Human papillomavirus 16,18".
      Fixed in this commit (catalogue `common_name` and both legacy prefix
      maps). `HUM_PAP_1618_*` genes are HPV16 only; HPV18 is NC_001357
      under `HPV18_*`.
    - Low-level calls present in both samples at similar levels (HPV118 15–18,
      HPV29 17, HHV-6 11–12) are background, not HPV16-specific.
  - **Read-level validation (2026-10-02, job 25694035; scripts in
    `hpv16_gse189670/readcheck/`): HPV16 is real.**
    - The scan found 13,119 R2 reads with a non-homopolymer HPV16 31-mer
      in the HPV16 raft, against 10 in the normal raft.
    - Competitive minimap2 against GRCh38 + the viral panel:
      - 11,349 (86.5 %) are clean, full-length and HPV16-best;
      - 38 (0.3 %) are host-best;
      - 33 are low-complexity.
    - 99.2 % read in the mRNA-sense direction (11,254 vs 95), so this is
      RNA, not DNA.
    - The pile-up sits just upstream of the early polyA site: 4,173 reads
      start at nt 3,000–3,250 of the indexed layout, which numbers from the
      E1 ATG, so E6 = 7,125.
    - **Canonical HPV16 splice junctions** (K02718 numbering, on a genome
      rotated to start at nt 7,500):
      - 880^3358 (E1^E4): 81 reads;
      - 880^2709: 29;
      - 226^409 (E6*I): 25;
      - 880^3391/3361: 22;
      - 226^526: 1.
    - gget BLAST of 12 reads sampled along the genome: 11 are full-length
      HPV16 at e = 7e-38, and 1 is a 71-nt HPV16 chimera.
    - The normal raft's molecules are index hopping: 4 of its 9 clean
      HPV16 CB-UMI pairs also occur in the HPV16 raft.
    - Caveat: ~5,000 reads sit in E6/E7/E1, far from the early polyA. They
      are still sense-strand and clean, so they are most likely internal
      priming on viral RNA. Not investigated.
- [x] `SW-24` — `MoleculeAudit.validate` (`multimapping.py:72`) compares
  allocated ambiguous mass with an absolute `atol=1e-9`. On EBV SRR12682296
  with the CAT-42 index, 10,250,998 ambiguous molecules sum to
  10,250,998.000000002 (diff 1.9e-9, relative 1.8e-16, float rounding), and
  the run fails (job 25691839, probe 25694067). Fix: a relative tolerance,
  plus a test with ~10M fractional shares. Blocks CAT-42 check 3.
  - Done 2026-10-02 (`v3/sw24-mask`, 08ee19a): `rtol=1e-9, atol=1e-9`.
    `tests/test_multimapping_tolerance.py` covers 3.3M thirds in a csr matrix,
    a +0.5 mismatch that still raises, and the observed 10,250,998 + 2e-9
    case, which `rtol=0` would reject.
- [x] `SW-23` — `--cell-calling knee` puts the knee in the empty-droplet
  tail (found 2026-10-01). `cellcalling.knee_cells` takes `argmin` of the
  signed distance (the point furthest *below* the chord). On a barcode-rank
  curve the cell plateau lies above the chord and the empty tail below, so
  the knee lands at about `knee_min_umi` (10).
  - Every knee run this week logged "knee at total>=10": HHV-6B 15,800,
    covid 118,061, HPV16 60,475, normal 86,135 cells.
  - On GSE189670 the steepest descent sits at rank ~8,000 / ~11,000, and
    `argmax` (above the chord) at rank 3,061 / 4,690.
  - Every `infected_called` / `pct_infected_called` from a knee run is
    inflated in its denominator.
  - **Decision (user, 2026-10-01): emptyDrops for every reported number.**
    `knee` stays in the code only because the frozen protocol lists it
    under `sensitivity_only_callers`. Its estimator is not fixed.
    (2026-10-05: DEF-00 A7 removed `knee` from `sensitivity_only_callers`, so
    the protocol no longer needs it. The CLI option stays, warned as
    sensitivity-only. Removing or fixing it is now a separate decision.) Every knee
    run now logs a WARNING (sensitivity-only, SW-23), and the help text and
    docs say so (`cli_reference.md`, `faq.md`, `output_reference.md`).
  - The three tracked SLURM scripts (`slurm_quant_covid_exploratory.sh`,
    `slurm_quant_hhv6b_5p.sh`, `slurm_quant_max_corrected.sh`) now use
    `--cell-calling emptydrops`. They put
    `conda/envs/R4_51/bin` (R 4.5.1, DropletUtils 1.30.0) at the end of
    `PATH`, because the bench env has no R. The quickstart vignette is also
    switched.
  - `TestKneeCells` passes with the bug because its data is a trivial
    bimodal split. It was left unchanged, and the new test only asserts the
    warning.
  - Check before quoting: the covid `pct_infected_called` figures under
    `HOST-03` (143,243 → 28,921 called cells) may come from knee runs.
- [ ] `ANELLO-PRIOR` — anelloviruses are a commensal virome, detectable in
  most people without disease (user, 2026-10-01). So an anellovirus call in
  human tissue or blood is biologically expected and must not be filtered as
  noise by family. Any gate has to be read-level: are the assigned reads
  low-complexity (homopolymer/poly-A/poly-G), and do they land on viral
  sequence along the genome, near the viral polyA, or only on a homopolymer
  tract? This applies to F-005, F-019, F-021, CAT-42 and DEF-01 wording.
  - First application (2026-10-01): on GSE189670 normal raft, 0 of 27,098
    KP343824.1-hitting reads (20M scanned) carry any non-homopolymer viral
    31-mer, so that call is the poly-A sink (F-021). The scan script
    (whole-genome k-mer position profile) is the template for the gate.
  - Detection audit (2026-10-01, F-022), against Kane et al. 2026 (plasma
    DNA prevalence 79 % young / 100 % older). No human anellovirus call
    has been read-validated yet. The default call is 1 UMI on one gene,
    with no read-level gate. The CDS-only models leave out the 5′ region
    (median 513 nt) and, for ~10 % of genomes, ≥200 nt before polyA. The
    tonsil whole-genome null (F-010) still says cellular mRNA is rare.
    Sub-steps, in order:
    - [ ] `ANELLO-PRIOR.1` plant positive control: held-out genomes,
      3′/5′-end reads into a real host library. Includes a TTV-high
      public dataset if one exists.
    - [ ] `ANELLO-PRIOR.2` rerun `measure_kmer_capture.py` on CDS exons
      vs exons + UTR. Its F-011/F-012 numbers were measured on whole
      genomes. Extend models to the mRNA extent only (never whole genome),
      judged on the negative controls too.
    - [~] `ANELLO-PRIOR.3` default read-level **diagnostic** (not a gate).
      - **Re-specified 2026-10-04.** The original spec — "homopolymer/entropy
        fraction and position along the genome" — names two measures the
        adversarial review ruled non-diagnostic. 10x 3' chemistry sees only
        transcript ends, so a genuine anellovirus read is poly-A rich by
        definition; and an artefactual poly-A read lands on the genome's
        longest templated A-tract, which in an anellovirus *is* the polyA
        site. Position and homopolymer fraction cannot separate the
        hypotheses. See F-019, update 2026-10-04.
      - **User decision, 2026-10-04: tune for not missing a real infection.**
        Flag, never filter; the confirmatory checks come afterwards. This
        follows `ANELLO-PRIOR` — anelloviruses are commensal, a call is
        expected, and every homopolymer-based filter removes precisely the
        genuine 3'-end reads needed to prove one. Such evidence *bounds* a real
        infection; it never excludes one.
      - **Done 2026-10-04.** Three measured columns, all labels:
        - `alignment_complex_body_fraction` — the decisive measure. The body
          5' of the first >=15-nt run must be >=20 nt at >=2.0 bits
          dinucleotide entropy. Thresholds measured, not assumed: over 5,000
          random ACGT draws the genuine minimum is 2.21 bits at 20 nt, while
          the artefact classes top out at 1.70 (poly-A 0.00, AC 1.00, CAG
          1.58, A-rich 1.70). 2.0 sits in that gap. The review's 3.5 was
          measured on 90-nt reads and would discard half of all genuine 25-nt
          bodies. `MIN_BODY_LEN` is 20, not the review's 25, because a 3'UTR
          ending in A's merges them into the tail (<=6 nt lost, measured).
        - `alignment_tso_fraction` — TSO in R2 in either orientation at <=2
          mismatches. It cannot occur in a genuine molecule: 10x 3' chemistry
          sequences only fragments carrying the bead-oligo end. Entropy does
          not catch the TSO (H = 3.56), which is why both columns exist.
        - `alignment_median_query_coverage` — identity without the length it
          was measured over is not evidence. **Needed before `ANDET-09f`,**
          which sweeps the very coverage filter that currently bounds it.
      - Also fixed: `scripts/plant_anello_10x.py` gave its `3p` plants no
        untemplated poly-A tail, so any poly-A-sensitive measure scored on that
        set would have looked free. `truth.tsv` now records `body_len` and
        `tail_len` per read, so recovery can be scored against body length.
        Pinned by `tests/test_plant_anello_tails.py`.
      - Caught in review, fixed before the columns were trusted: SAM stores
        SEQ reverse-complemented on a reverse-strand record, which **inverts**
        both read-side measures — a genuine `[body][poly-A]` arrives as
        `[poly-T][rc body]` and scores as artefact, while a `[poly-A][TSO]`
        chimera arrives TSO-first and scores as complex. With `--soloStrand
        Unstranded` that is about half of all records. `Alignment.read_seq()`
        now restores sequencing orientation; pinned by a test that fails
        against the unfixed code.
      - [x] Remaining, and **the four arms as specified will not give a
        reading**: covid x213 and the synthetic negative yield ~0 aligned reads
        through the branch (its 0.80 coverage filter already rejects the
        chimeras), so their columns come out *empty*, which is "not measured",
        not "artefact"; and HPV16 is not an anellovirus, so the branch never
        touches it. Validate instead by running `is_complex_body` / `has_tso`
        directly over the 19,785 covid anellovirus-aligned reads already
        extracted — which is also the full-set body census that turns the
        review's 0/30 into a bound. Un-revcomp any read taken from BAM SEQ.
        - **Done 2026-10-04** (`scripts/anello_body_census.py`; F-019 update
          "full body census"). Pre-registered bound: ≤ 81/19,785 reads
          (CP95 0.51 %), ≤ 62/2,791 CB+UMI (CP95 2.84 %). Every one is a
          low-complexity G/C/A tract. Post-hoc refinement: 13, none with a
          body > 26 nt. The refinement leans on `is_complex_body`, so it is
          not the bound.
          **The measures failed validation**: `has_tso` misses TSO truncated
          at the read edge, and `read_body` takes a 5′ TSO/TruSeq fragment as
          the body (H ≈ 3.5). 6,171 reads kept, 98.9 % of them with a body that
          places nowhere. Fix tracked in the row below.
      - [x] Fix the reagent leak found by the census (2026-10-04). `has_tso`
        should match a TSO overlapping the read edge (flag-only, fits the
        locked decision). Pin it with the census's real reads as fixtures, and
        re-run `scripts/anello_body_census.py` to show the leak closed.
        **User decision:** do TruSeq R1 and other reagent hits go into
        `tso_fraction`, or into a renamed `reagent_fraction`? And should
        `read_body` strip a leading reagent before measuring the body?
        - **Done 2026-10-04.** The user approved the five changes from the
          Biomni review `tsk_015vCaV6Sg0UDc0WwVBVOeVb` ("accept with
          changes"):
          - `reagent_fraction` replaces `tso_fraction`. It flags forward
            TruSeq R1, or a TSO in junction context (no complex body after
            it). TSO followed by a complex body is a genuine full-length
            molecule, per 10x.
          - New `r1_tso_fraction`, from raw `CR`/`UR`.
          - `read_body` strips a leading reagent.
          - Alignment overrides entropy.
          - The 1.70-bit claim is corrected.
          - Census re-run: kept 6,171 → 790; the bound and sensitivity
            (69/81) are unchanged; the 4,762 edge-TSO + 459 TruSeq leak is
            closed. Known ceiling: ~50 kept reads carry a degraded TSO
            (≥ 3 mm). F-019 "evening" update.
      - [x] Audit the TSO-in-R1 "barcodes": 815 reads, 4 CB+UMI. Were they
        whitelist-corrected and counted in the kallisto/bustools 57,715-UMI
        call? If so, the call is attributed to real cell barcodes (Biomni
        review, point b).
        - **Done 2026-10-04: no, they contribute 0.** None is within Hamming
          distance 1 of the whitelist. There are 0 records in the corrected
          `output.unfiltered.bus` and none in `cells_x_genes.barcodes.txt`.
        - Found on the way, a correction: the census evidence set was built
          from the uncorrected BUS plus the raw merged FASTQ (July job
          `hf_evidence_25180994`), not from the counted reads. The 81/62 bound
          therefore bounds the evidence set and must not be quoted against
          the 57,715 UMI (F-019, "Correction").
      - [ ] Re-run the body census on the **counted** reads. Steps:
        `viralscan evidence` with exact lineage (`kb_r1`/`kb_r2`) for
        Alphatorquevirus on covid x213, then `scripts/anello_body_census.py`.
        That gives the bound comparable to the 57,715. Blocked by `REL-16`.
      - [ ] Deferred from the review: check that rc(UMI)/rc(CB) in R2 match R1
        on read-through reads; add a fixture of host R2 reads from the same
        library showing no forward-TruSeq starts. Only needed if
        read-through reads become common in a call.
      - [x] The columns are reachable only with `--anello-align`, which ships
        - **Done 2026-10-04.** `viralscan evidence` writes `complex_body_fraction`
          and `reagent_fraction` per reference in `alignment_qc.tsv`. They
          come from the same pure functions, on reads restored to
          sequencing orientation, and they are labels only.
          - `r1_tso_fraction` is omitted by design: exact-lineage reads carry
            corrected on-list barcodes only.
          - Reverse/forward parity is pinned by a test.
          - Any kallisto call can now get the columns from one documented
            follow-up command. An automatic run is a separate decision.
        off, so a default run's anellovirus call still carries nothing but
        `artifact_risk`. Apply the same pure functions in `viralscan evidence`,
        which `output_reference.md` already names as the route for checking a
        kallisto call's reads. **This puts `REL-16` (the segfault that breaks
        `viralscan evidence`) ahead of `ANDET-09f`.**
      - [ ] Expect the regenerated planted arm to lose roughly a quarter of its
        3' reads: a read starting >228 nt into the 300-nt window carries a tail
        of >=19 nt, dropping coverage below 0.80. That is the branch filter, not
        the new columns, and it is exactly the sensitivity cost the permissive
        decision rejects — so score 3' recovery by `body_len` bin from
        `truth.tsv` and make that the axis of `ANDET-09f`. Re-check criterion 3
        as well: A-rich tails may now be lost to the host filter.
    - [x] `ANELLO-PRIOR.4` (user decision) relabel `eve_risk` for
      Anelloviridae; its basis is the F-005 mechanism F-019 revised.
      - **Done 2026-10-03** (user chose "clear + interim note" in a grill).
        Anellovirus `eve_risk` is now False.
        `viral_summary.tsv` gains `artifact_risk`, set to `low_complexity`
        (F-019/F-021). It is a diagnostic label, never a filter, until
        `ANELLO-PRIOR.3`'s measured metric replaces it.
      - Implementation:
        - All 3,448 catalogue `risk_class=eve` rows became
          `low_complexity`. They are Anelloviridae only, and no other column
          changed.
        - `EVE_RISK_GENERA` is now empty. The mechanism is kept for a real
          EVE family; inherited ciHHV-6 is parked.
        - The new `LOW_COMPLEXITY_RISK_GENERA` seeds the label, and the
          genus-name fallback covers uncatalogued names.
      - Basis: a quick Europe PMC search found no germline human anellovirus
        EVE; the one integration is somatic, in the SKNO-1 cell line (PMID
        42671192).
      - cat42b probe: 8 anellovirus groups are labelled, the other 222
        viruses are not, and `eve_risk` is False everywhere.
    - MECH-B (group → sum → threshold) must close before any anellovirus
      sensitivity claim.
- [x] `CAT-42` — homopolymer mask misses runs under 31 nt (F-021,
  2026-10-01). KP343824.1 starts with 29 T and captured poly-A tail reads
  in GSE189670 (10,037 molecules in an anellovirus-free cell line). The 3′
  EBV and HSV-1 max-panel runs put 0 on it, so the capture is
  library-dependent.
  - Panel scan (2026-10-01): 38 runs of ≥20 identical bases in 38 of 2,343
    genomes (32 poly-A, 3 T, 2 G, 1 C; 9 are 25–30 nt). Every
    Gamma/Betatorquevirus molecule in both GSE189670 rafts (5,688 + 941 and
    10,037 + 727) landed on one of those genomes. The ~250 "genus
    unassigned" molecules did not.
  - `build_bundled_panel_ref.py --homopolymer-run-length` default is now
    20 (was 31). Rebuild job 25691632 → `viral_ref_cat42/build/`. After it
    finishes:
    1. check `kallisto inspect` against CAT-40 (max EC size, ECs
       discarded);
    2. rerun GSE189670 (both rafts) for anellovirus ~0 and HPV16
       unchanged;
    3. rerun the EBV/HSV-1 regression (molecules within ±5 %).
  - Queued 2026-10-01 with `afterok:25691642`:
    `viral_ref_cat42/quant_gse189670.sbatch` (25691838, emptyDrops, output
    in `hpv16_gse189670/runs_cat42/`) and `viral_ref_cat42/regression.sbatch`
    (25691839, output in `viral_ref_cat42/runs/`). Both use
    `viral_ref_cat42/viral_panel.gtf`, cut from the new `combined.gtf`. In
    it, anellovirus gene IDs carry their accession prefix, so 124 gene IDs
    differ from the old panel GTF; HPV16 is unchanged.
  - The rename hits 124 RefSeq anellovirus genes (`TTV3_gp1` →
    `NC_014081.1_TTV3_gp1`, `D1P64_gp1` → `NC_038336.1_…`). Identity is not
    affected, because t2g column 5 is still the genome accession, so those
    genes stay catalogued by accession. The output gene IDs do change.
    Before this index ships, find the commit that changed the builder's
    gene-ID naming (not yet traced).
  - Pass criterion for check 2, judged by accession, not by genus label:
    molecules on KP343824.1 and on the other 37 genomes with a ≥20-nt run
    are ≈0, and HPV16 stays 7,725 molecules / 2,012 cells. *Corrected
    2026-10-02:* the ~250 "genus unassigned" molecules are **not** a
    positive signal. They appear equally in the normal raft, an
    anellovirus-free cell line, so they are artefacts too.
  - **Composition mask (2026-10-02, `v3/sw24-mask`, 665b274):**
    - `build_bundled_panel_ref.py --lowcomplexity-kmer-mask`, off by default.
      It N-masks each 31-nt window with ≤2 distinct bases and a base ≥20, or
      any base ≥28, in **Anelloviridae records only**. A review simulation of
      the panel-wide rule masked EBNA-2 and HSV-1 s-gene sequence. It also
      writes `lowcomplexity_mask.tsv`.
    - Dry run on `viral_ref_cat42/build/viral.fa`: 1,586 bases masked across
      46 records. HPV16, EBV-1/2, HSV-1 and HSV-2 lose 0 bases. KP343822.1
      loses 34 bases, all inside nt 2405–2460. KP343842.1 loses 35 of the 38
      nt at 3306–3343.
    - The rebuild goes into `viral_ref_cat42b/`. It is checked against cat42,
      so the mask is the only variable.
  - **Results (2026-10-02):**
    - Build 25691642: `BUILD_OK`, 1,346 bases masked. `kallisto inspect`:
      max EC 9921, 0 ECs discarded, 471,944 targets. **Check 1 passes.**
    - GSE189670 (25691838): Gamma/Betatorquevirus 10,037 → **0**. HPV16
      is still 7,725 molecules (2,007 / 15,502 called cells). "Genus
      unassigned" stays at 246 / 273 in both rafts (KP343822.1,
      KP343842.1), so the artefact persists below the 20-nt mask.
      **Check 2 passes for KP343824.1, with one residual artefact.**
    - Regression (25691839): HSV-1 completed. EBV failed in multimap, on
      the conservation check, with `SW-24`. **Check 3 is blocked.**
  - **cat42b results (2026-10-03; build 25694852, checks 25694853/25694854,
    code d01043c). Every gate passes, so cat42b is now the current panel**
    (user decision 2026-10-02: replace if every check passes). Compared with
    cat42, so the composition mask is the only variable:
    - `kallisto inspect`: 0 ECs discarded, 471,944 targets, 84,860,192 k-mers
      (cat42: 84,860,149; kallisto fills N with pseudorandom bases).
    - GSE189670: "genus unassigned" 246 → **0** (HPV16 raft) and 273 → **0**
      (normal raft). KP343822.1 220 → 0, KP343842.1 26 → 0. HPV16 stays at
      **7,725** molecules (2,006 / 15,615 called cells; cat42 2,007 / 15,502.
      emptyDrops runs on the changed matrix, so the called set moved). The only anellovirus left is 1 Samektorquevirus molecule in the
      normal raft, as on cat42.
    - EBV SRR12682296: EBV-1 747,531 → 747,532, EBV-2 142,960 → 142,960.
      Gene-level changes: KP343822.1 65 → 0, KP343821.1 1 → 0, BGLF3/3.5/4
      +0.33 each. EBNA-2 unchanged. Against v1 (885,933): EBV-1 + EBV-2 =
      890,492.
    - HSV-1 SRR8315713: 23,187 → 23,187 (v1 23,170). KP343822.1 19 → 0. HHV-2
      is still `possible_em_bleed` (666:1). HSV-1 s-genes unchanged.
  - **Swap:** the golden identity tests now include `cat42b`
    (`tests/test_virus_identity.py`, `final` kept). CMP-06's "latest
    reference" is now `viral_ref_cat42b/build/` (see the note there).
  - **Gene-ID rename, traced:** the `{accession}_{token}` gene-ID rule
    (`ncbi_fetch.py:350`) came in with 3379b7c (2026-09-27). The 241
    gitignored bundled GTFs were regenerated with it at 2026-09-28 22:46,
    during the a97cc00 final-panel work, after `viral_ref_final` was built.
    So cat42 and cat42b carry prefixed IDs for 124 RefSeq anellovirus genes,
    and `final` does not. Identity is unchanged, because t2g keeps the
    genome accession.
- [~] `CAT-41` — catalogue display names that contradict their genome
  - 2026-10-04: the three real errors are fixed.
    - Influenza D: 7 rows.
    - NC_001664.4 → "Human herpesvirus 6A". The gene_programs.tsv HHV-6A
      markers and the 6A/6B sibling map were renamed with it, so no marker
      silently unresolved.
    - KC138720.1 → HPV type XS2.

    Pinned by `test_cat41_display_names_match_their_genome`. Smaller exemplar
    names (HPV61/2/14/68/6, Enterovirus D) were not reviewed.
  (found 2026-10-01 from the HPV16 "16,18" label, which is fixed). Still
  open:
  - all 7 Influenza D segments (NC_036615–21) display as "Influenza D virus
    segment 7";
  - HHV-6A (NC_001664.4) displays as "Human herpesvirus 6" (changing it
    changes sibling-note text);
  - XS2 (KC138720.1) displays as "Human papillomavirus 78";
  - smaller species/type mismatches: HPV61 for Alphapapillomavirus 3, HPV2
    for Alphapapillomavirus 4, HPV14 for 14D, HPV68 for 68a, HPV6 for 6b,
    "Human enterovirus 68, 70" for Enterovirus D.
  Re-run the scan in `CAT-41` after any catalogue merge.
- [x] `MECH-A2` — prefix-named rows lack catalogue decisions (residual of
  MECH-A step 4a, 2026-09-30). On a VIRTUS-style index (covid x213) HHV-6B
  resolves as `legacy_prefix`, key `name:Human herpesvirus 6b`, with no
  `sibling_group` or `risk_class`. The same virus is `taxid:32604` on the
  stored SRR20710641 index. So the HHV-6A/6B note can never fire there, and
  `eve_risk` takes the genus fallback. Candidate fix: adopt the catalogue
  row when the prefix name matches exactly one taxid's `common_name`.
  - Done 2026-10-01: a prefix name equal to a catalogue display name (these
    are unique per key) takes that key's row, keeping status
    `legacy_prefix`. The key's eve row is used when it has one, so
    `eve_risk` stays fail-closed. `virus_facts` uses the genus fallback
    only for `name:` keys, because `risk_class` is blank on 681 non-EVE
    rows. Covid x213 parity (job 25691885, against final4a): all 12
    multimap layers, X and `viral_summary.tsv` are identical. The 97
    `HUM_HERP6B_*` genes are now `taxid:32604` / `HHV-6`, and the
    genus-fallback warnings went from 1 to 0.
- [x] `MECH-B` — virus-level Detection: group, then sum, then threshold.
  `accession_breadth` becomes coverage over reference genes (today it is always
  1.0), and `sensitivity.tsv` gets zero rows.
  - Done 2026-10-04. `detect_genes(..., groups=)` thresholds each virus's
    summed count over its index genes and reports every nonzero gene of a
    passing virus. Viruses spread thinly over genes are no longer missed, and
    `viral_molecules_total_est` no longer drops their sub-threshold genes.
    At the default threshold of 1, calls are unchanged on integer layers
    only. On a fractional (EM) layer, sub-1 genes can now sum to a call, and
    trace bleed genes of a passing virus enter `found_genes`, which raises
    infected-cell and headline totals. That is the intended "do not miss"
    direction (corrected 2026-10-04; ce7d201's message says "unchanged").
    Without an identity table, unmatched genes group under their own ID, so
    their zero rows are named by gene ID (left as is). Breadth landed
    earlier as `ANDET-01`. `build_sensitivity_table(index_viruses=)` adds an
    `observed_molecules` = 0 row per undetected indexed virus. Expect many
    rows on large panels; nothing reads the file programmatically.
- [x] `MECH-C` — **done 2026-10-10.** A single Run Config writer; delete the Namespace → `k=v` →
  YAML round trip and `createconfig`. This is the root cause of `SW-13`.
- [~] `MECH-D` — **tracked under `DEF-02` (reconciled 2026-10-04): one module, one row;
  this closes with DEF-02.** (module landed 2026-10-03) Chemistry module: one geometry for kb, STARsolo and the
  preflight, plus an R1-length / polyT check. The 10x v2 run as `-x 10xv3`
  gave 1.78M "cells" with no error.
- [~] `MECH-E` — the reference pipeline: source adapters → one gate module →
  masking → GTF emitter → index with a recorded D-list. The prototype is
  `viral_panel_max_2026-09-28/01–04`.
  Software side done under `AUD-15` (2026-10-11/12): source adapters, `_prepare_viral_fasta` masking gate, `_emit_viral_gtf` emitter, `_finish_reference_build` index boundary with `genome_dlist_applied` in the manifest. Left: governance re-pin and a real (non-mock) index build.
- [~] `MECH-F` — the corrected-BUS boundary, index-aware denominators, and one
  - 2026-10-04 (1203250): **off-list drops are now counted.** The count audit gains:
    - `barcode_correction`;
    - `bus_records_raw` / `bus_records_after_correction` /
      `offlist_dropped_records`, and the matching reads fields;
    - `bus_totals_source`.

    Computed from `bustools inspect`, for both the ViralScan and the kb
    correction paths. The fields are optional in schemas/v3 (both copies
    identical). Index-aware denominators and the single gene-role catalogue
    stay open as architecture. Not yet run against real bustools; the tests
    fake the inspect JSON.
  gene-role catalogue for every family.

### WP1E — Grill decisions for 3.0 (user-confirmed 2026-09-29)

The user answered a structured grill (41 questions, three rounds) and confirmed
the summary on 2026-09-29. These decisions are binding for the rows below and
for WP3–WP11. A `DEF-*` row is the implementation work a decision creates.

| # | Decision | Answer |
|---|---|---|
| Q1 | Release before measured defaults | Internal pre-release only; outputs stamped `defaults_status=provisional` |
| Q2 | Coordination | One orchestrator writes PLAN.md, `.living/`, inventory, claim registry, protocol |
| Q3 | MECH-A owner | This main session (user: "Continue MECH-A here") |
| Q4 | Default panel | Current default stays until MECH-A, MECH-B and the sibling experiment land |
| Q5 | 10x 5′ | In scope, as `-x 10xv2`/`10xv3` plus strand handling |
| Q6 | Chemistry declaration | Auto-detect from reads; fail on ambiguity or a mismatch with `-x` (override flag) |
| Q7 | Host strategies | Combined and STAR two-step both candidates; D2 decides |
| Q8 | Host genes in `--gtf` | Index build manifest of host/viral genes plus the column-5 structural guard; manifest wins; contradictions refused |
| Q8b | Anellovirus unit | Genus; unassigned = "Anelloviridae (genus unassigned)" (done, a566dad) |
| Q9 | Canonical environment | `env_full`, locked (REL-02/03); CI adds a snakemake 9 dry-run; retire py3.8 env for repo code |
| Q10a | Subagent commits | Own `v3/*` branches; only the orchestrator merges; every merge gated |
| Q10b | Concurrent implementers | 3; shared files (`menu.py`, `Snakefile`, `multimap.py`) serialised |
| Q10c | Downloads | Approved: chemistry subsamples, PBMC backgrounds, comparator envs; checksum-pinned manifests |
| Q10d | Tuning compute | 20k core-hours |
| Q11 | Covid anellovirus numbers (CAT-30) | Files frozen; retraction notice "poly-G no-signal reads ≈90 %, host-homologous ≈10 % (F-005, F-019); low-level divergent anellovirus not excluded"; claim retracted in registry |
| Q12 | Whole-genome features (CAT-38/39) | Virus-level detection only; tagged `whole_genome`; excluded from gene-level outputs and programmes; low-complexity filter required |
| R2.0 | Read-artefact filter | Default pre-count filter (homopolymer ≥15 nt, TSO/adapter, low complexity) with audit table, plus reference masking |
| R2.1 | Preregistration | `defaults_selection` in `protocol.yaml` (grid, objective, hard constraint, 1-SE conservative tie-break, sample bootstrap); reviewed, signed off, frozen before tuning |
| R2.2 | Already-seen data | Excluded from tuning and holdout; context only |
| R2.3 | Specificity constraint | Zero probable/strong calls on host-only and planted-homology negatives binds every default |
| R2.4 | Tuning reference | No tuning run until G4 passes; reference chosen after G4 |
| R2.5 | Sibling allocation (F-017) | Within-sibling-group variant (unique support or group EM) in D3; counting contract and conservation claim amended first |
| R2.6 | Strand | `--strand` option; otherwise inferred from a 1M-read pilot; recorded in manifest |
| R2.7 | Per-chemistry defaults | Allowed; live in the Chemistry module; preregistered per chemistry |
| R2.8 | Cell calling under two-step | `--called-cells` if given; else STARsolo host `GeneFull` counts into EmptyDrops (SCI-02) |
| R2.9 | Generator backgrounds | Real PBMC primary plus synthetic GRCh38-transcript control |
| R2.10 | Reruns now | All exploratory reruns now (see `EXPL-R2.10`) |
| R3.1 | Compute split | Shared defaults once on 10x v3 (~55 %); chemistry-specific per chemistry (~25 %); per-chemistry confirmation (~10 %); reserve (~10 %) |
| R3.2 | Decision rules | Per-default objective: recall (host strategy, STAR, read filter); weighted allocation L1 per Index Kind (multimap); cell Jaccard vs truth labels (cell calling, never viral outcomes); SCI-03 F1 (thresholds, tiers). No passing option → most conservative, reported "screening only" |
| R3.3 | Comparators | Full SCI-04 matrix, both index arms; ~10k core-hours, separate budget |
| R3.4 | Holdout blinding | Seeded partition; truth outside repo, owner-only, hash in protocol; only the H1 brief names the path; unblinding needs user's written sign-off after the defaults freeze |
| R3.5 | Releases | `3.0.0.devN` git tags only; `3.0.0rc1` after G2/G5/G6a; `3.0.0` after G6/G7. User sign-off at: protocol freeze (G3), defaults freeze, holdout unblinding, rc1 publication, 3.0.0 tag, any retraction |

Implementation rows:

- [~] `DEF-00` — amend `analysis/v3_validation/protocol.yaml`:
  - add the `defaults_selection` section (R2.1, R3.1, R3.2);
  - add 5′ to the chemistry scope (Q5);
  - add the within-sibling D3 variant (R2.5);
  - add the per-chemistry defaults (R2.7);
  - add the both-background generator (R2.9);
  - add the blinding layout (R3.4);
  - decide the knee sensitivity arm: `sensitivity_only_callers` requires
    `knee`, whose estimator is defective (SW-23, 2026-10-01). Fix it or
    replace it in the amendment.

  Then an independent review, user sign-off, and the freeze (G3).
  - 2026-10-05: **written**, protocol v0.3.0, ledger `DEV-020` (all 25 rows
    plus §F; row-to-location map in
    `docs/plans/2026-10-05-def00-protocol-amendment.md` §G). Every factor and
    section stays `pending`.
    - §F: `host_virus_homology` = `none`/`repeat_homology`. Decided by Claude at
      the user's instruction, because REF-07 found the homology is
      low-complexity (F-023).
    - Knee: removed from `sensitivity_only_callers`, not fixed (SW-23).
    - B2: real-background datasets split by donor library within chemistry,
      with the full 7 × 2 grid nested in each library. Synthetic-background
      datasets keep the generated sample as the unit.
    - Schema, validator and packaged schema now know D25, `five_prime` and
      `defaults_selection`. `defaults_selection` is digested under
      `frozen_inputs`.
    - Gates: draft valid; `--phase training` shows only freeze blockers;
      ledger append-only; 150/150 protocol tests pass.
    - Two governance tests still fail until commit: they pin the sha256 of
      `protocol.yaml` and `deviations.yaml`, so these need a re-pin commit.
    - Remaining for G3:
      - [x] fill the empty `defaults_selection` grids. Three are written
        (DEV-021, 2026-10-05): host_strategy, multimap_allocation, and
        cell_calling (one point). `star_alignment` is `[pinned, star-default]`
        (DEV-022, 2026-10-05), with star-default named the conservative point
        except on multimapping (user confirmed 2026-10-05).
        `read_artefact_filter` is `[off, artefact]` (DEV-023, 2026-10-05; user).
        artefact is the strictest point;
        - [x] STAR-parameter knob: `--host-filter-star-params {pinned,star-default}`
          (`host_filter.py` `STAR_FILTER_PARAM_SETS`, 2026-10-05). `star-default`
          passes STAR 2.7.11b's defaults explicitly. Unset leaves the run manifest
          unchanged (`omit_when_unset`), and the audit records `star_param_set`;
      - [x] reconcile D18's "exact truth-labelled" population with twin-difference
        specificity. **Done 2026-10-06 (DEV-025, user-confirmed).** On a real
        background, D18 and D19 drop every anchor cell the unplanted twin calls
        for the same virus, at any tier, and report those cells separately as
        background signal. This is the same per-cell exclusion as D25. D19 is
        included because E9 scores both endpoints on one population.
        Synthetic GRCh38 backgrounds have no twin and are unchanged;
        - [x] DEV-026 (2026-10-06): the `datasets` planting_rule no longer says
          "planted-minus-twin difference" (it read as rate subtraction); the
          DEF-00 amendment doc marks D18 resolved;
      - run SCI-05.
- [x] `COST-01` — measure the real per-sample run cost for kb count and the STAR
  two-step, in core-h and peak RAM, before sizing the `defaults_selection` grids
  (VAL-01 D5, user 2026-10-05).
  - **Done 2026-10-05 (F-024).** Measured on 7 existing full-depth jobs on real
    libraries rather than the pilot panels; the user approved this because
    the panels need VAL-01, which waits on G3. Per 50 M pairs: combined about
    1.3 core-h, two-step about 5.6. The 2 × 4 grid costs 7.4–18.2 core-h per
    training sample. Script: `scripts/cost01_per_sample_cost.py`.
- [x] `SW-25` — `reference_strategy.py:723` builds the two_step row with
  `--host-filter kallisto`, but the CLI accepts only `starsolo`
  (`menu.py:1478`), so the harness two_step row cannot run.
  - **Done 2026-10-07.** The row now passes `--host-filter starsolo --host-index
    <references.starsolo.human_only.genome_dir>`. `tests/test_reference_strategy_benchmark.py::
    test_viralscan_two_step_row_parses_with_the_real_cli` parses the built argv with
    `menu.build_parser()` (fails without the fix). The manifest's
    `references.viralscan.human_only.kallisto_index` is now unused by this row; left in place.
- [x] `DEF-01` — read-artefact filter before `kb count`, with an audit table
  (R2.0, F-019). Reference homopolymer/low-complexity masking stays under
  CAT-17.
  - **Done 2026-10-05.** `--read-filter {off,artefact}`, opt-in and off by
    default (user, 2026-10-05; keeps the 2026-10-04 flag-not-filter decision).
    `scripts/read_filter.py` runs after any host filter. It drops, in this
    order: `r1_tso` (TSO in the CB+UMI span), `r2_reagent` (`has_reagent`), and
    `r2_no_complex_body` (`is_complex_body` on R2 *or* its reverse complement).
    All three reuse the `anello_align` classifiers. The either-orientation rule
    was added the same day (DEV-024): 5′ R2 is antisense, so a genuine
    `[poly-T][body]` read was being dropped. It writes `read_filtered/read_filter_audit.tsv`, with pinned
    constants, and `fragment_lineage.tsv.gz` in the host filter's schema, so
    `lost_truth_counts` works on it.
  - R2.0's literal "drop any read with a homopolymer of 15 nt or more" is not
    used: it removes genuine `[body][poly-A]` reads. Unset leaves the run
    manifest unchanged. `anello_align` now takes its inputs from `kb_r*`, which
    closes a race once the filter is on.
  - Speed: exact pigeonhole shortcuts in `_find` and `_edge_len` (equivalence
    tests). The filter takes about 80 µs per pair on one core, which is about
    1.15 core-h per 50 M pairs.
- [~] `DEF-02` — **priority after F-020**, because every 5′ run so far used
  about 7–9 % of its reads. The Chemistry module (`MECH-D`), covering:
  - auto-detection and fail-closed checks (Q6);
  - `--strand` and the pilot inference (R2.6, strand test job 25672276);
  - per-chemistry defaults, including SW-21 barcode correction (R2.7).
  - [x] `--strand {forward,reverse,unstranded}` plumbing: 2026-09-30, branch
    `v3/strand` (d8bf730). The option reaches `kb count` through
    `params.strand`; unset omits the flag, so kb keeps its per-technology
    default. The default stays unset until DEF-01, because unstranded
    re-admits the forward-only poly-G artefact (F-019).
    - `run_manifest.json` records `options.strand` only when it is set. So an
      old manifest still resumes when `--strand` is not given. An old manifest
      with `--strand X`, or a different explicit strand, refuses with a
      reason. Tiny e2e asserts that `kb_info.json` records `--strand
      unstranded`.
    - `reverse` fails on the 2-read fixture (kallisto exit 1), which is why
      the e2e uses `unstranded`.
  - [ ] Strand pilot inference and the per-chemistry default (the rest of
    DEF-02, MECH-D). Next: rerun the 5′ positives with `--strand reverse` and
    `unstranded`. These reruns are exploratory; they do not set the default.
    - **Reruns done (2026-10-03; full runs, emptyDrops, code d01043c;
      results in F-020's update).**
      - Covid x213/x216, reverse vs unstranded: 47.8/53.6 % and 46.2/54.4 %
        pseudoaligned, against 6.4/8.9 % with no strand.
      - HHV-6B: 2,405 molecules reverse vs 5,596 unstranded.
      - **Correction:** the HHV-6B baseline already ran unstranded in effect
        (60.2 %; kallisto 0.52.0, `-x 10xv2`, no strand flag). The covid
        baseline behaved as forward (kallisto 0.51.1, `-x 10xv3`, no flag).
        Whether version or chemistry sets the default is not separated. So
        the pilot must pass all three strands explicitly.
      - **Resolved by `EXPL-TTMDV` [x] (artefact) (reconciled 2026-10-04).** Original note:
        under reverse, covid x213 puts 55k + 25k
        molecules on two TTMDV placeholders (AB303552.1, AB303557.1), in
        64.9 % of called cells. These genomes have no homopolymer over 12 nt.
        Needs a read-level check before any claim. Logged, not started (no
        new tasks in this pass).
      - `--strand auto` check: reverse/unstranded = 0.85–0.89 on all three
        5′ libraries, forward/unstranded = 0.11–0.16. So τ = 0.8 picks
        reverse. For HHV-6B, reverse is the choice that loses 57 % of viral
        molecules.
    - **`--strand auto` landed (2026-10-03, 8ff6bec, opt-in; the default
      stays unset).**
      - `src/viralscan/strand.py`: the pilot runs on the first 1M pairs, with
        three **explicit** strands, and reads `p_pseudoaligned`.
      - Rule `infer_strand`, τ = 0.8: the larger qualifying ratio wins; a tie
        or unstranded ≤ 0 gives unstranded.
      - The manifest keeps `options.strand="auto"` and adds a per-sample
        `strand_inference` block (choice, rates, pilot_reads, tau) outside the
        fingerprint. Resume reuses the recorded choice; auto vs an explicit
        strand refuses.
      - Tests: `tests/test_strand_auto.py`. The three measured F-020 rate sets
        pick reverse; the forward and unstranded branches are synthetic only.
      - Still open: a measured 3′ ratio (done 2026-10-03, see (d) below) and
        a real-`kb` end-to-end pilot run.
  - **Chemistry module, 2026-10-03 (user decisions the same day).**
    - Decisions: SW-21 → bypass kb's barcode correction for a technology with
      no official on-list, so ViralScan's cell calling is the only caller.
      Drop-seq is inferred when no 10x on-list matches and R1 is 20 bp.
      Per-chemistry default *values* stay with DEF-00 (R2.7 preregistration);
      today's behaviour (strand unset) ships as provisional.
    - [x] (a) `src/viralscan/chemistry.py`: one table (`CHEMISTRIES`, kb name,
      CB/UMI length, ngs_tools on-list); `evidence.cb_umi_geometry` re-exports
      it (the unused `10xv3_5p` entry is gone). `detect()` samples 100k R1
      reads, streams each on-list once, and reads the UMI length from the
      TSO/poly-T start (26/28) or a trimmed R1 length. `resolve()` fails on
      ambiguity, on samples that disagree, or on a mismatch with `-x`;
      `--force-technology` keeps an explicit `-x`. Tests:
      `tests/test_chemistry.py`.
      Real libraries (19 s for all four): EBV SRR12682296 → 10xv2 (97.3 % v2
      list); HHV-6B SRR20710641 → 10xv2, 5′ (88.2 %, TSO at 26); covid x213
      → refused without `-w` (best bundled list v4 4.4 %), 10xv3 geometry
      with Cell Ranger's list (67.2 %); HSV-1 SRR8315713 → dropseq. EBV run
      as `-x 10xv3` is refused (the MECH-D failure). No bundled list covers
      GEM-X 5′. No local 10xv3 3′ library exists to check v3/v4 overlap;
      two lists above 0.5 fail closed. Bounded from the lists instead: v3 ∩ v4
      = 68,254 barcodes (1.0 % of v3), v2 ∩ v3 = 77,142 (10.5 % of v2), so a
      10xv3 library cannot clear 0.5 on two lists.
      Fix (same day): Drop-seq with its own `-w` list resolves to dropseq
      (a 20 bp R1 gives no UMI length, so it used to be refused).
    - [x] (b) run preflight, 2026-10-03: `-x` defaults to detection,
      `--force-technology` keeps an explicit `-x` (excluded from the
      fingerprint). Per-sample `chemistry_detection` blocks go into
      `run_manifest.json` outside the fingerprint (`record_manifest_block`,
      shared with `strand_inference`). The advisory `_whitelist_preflight`
      is gone; `check-whitelist` stays. Resume: the resolved `-x` is what is
      fingerprinted, so an old run made with the old `10xv3` default resumes
      when its reads are 10xv3, and refuses (with the reason) when they are
      not. Docs: `docs/faq.md`, `docs/cli_reference.md`. Tests:
      `TestResolveChemistry`, `test_force_technology_does_not_change_the_fingerprint`.
      Note: slurm scripts that run HHV-6B as `-x 10xv3` now stop; it is 10xv2.
    - [x] (c) SW-21 bypass, 2026-10-03: `chemistry.kb_whitelist_arg` gives
      `kb count -w None` for a chemistry with no on-list and no `-w`. kb 0.29.5
      then skips `bustools correct` and writes no `output.unfiltered.bus`, so
      multimap starts from the raw BUS (`correction="none"`, already logged).
      Tests: `TestKbWhitelistArg`, a Snakemake dry-run render. Not yet rerun
      on HSV-1 SRR8315713; expect ~1.8M barcodes and ~32k HSV-1 molecules.
    - [x] (d) measured 3′ ratio, 2026-10-03 (job 25695946,
      `scripts/slurm_strand_pilot_3p.sh`): EBV SRR12682296, 10xv2 3′, cat42b
      host+viral index, first 1M pairs, real `kb`. Pseudoaligned: forward
      72.6 %, reverse 2.8 %, unstranded 73.8 %. forward/unstranded = 0.98,
      so τ = 0.8 picks **forward**, as a 3′ library should. The forward
      branch is no longer synthetic-only.
    - [~] End-to-end check of (b)+(c): HSV-1 SRR8315713 run with **no `-x`**
      (job 25695959, `scripts/slurm_quant_sw21_hsv1.sh`, output
      `runs/combined_sw21/`). Pass: detection → dropseq, kb log shows
      `-w None`, and barcodes/HSV-1 molecules return to ~1.8M / ~32k (vs
      5,468 / 23,170 with the knee allowlist, job 25666447).
- [x] `DEF-03` — index build manifest of host and viral gene sets, overriding
  `--gtf` (Q8). Fix `reference_strategy.py:704` so it passes a viral-only
  GTF. Folds into MECH-A step 4.
  - 2026-09-30, branch `v3/def-03-build-manifest`. `build-ref` and
    `--reference` write `<index>.build_manifest.json` (sorted host and viral
    gene IDs, plus provenance sha256) after a successful `kb ref`. The
    `analysis` rule passes it to `build_identity_table(build_manifest=...)`.
  - Contradiction rule (user-confirmed 2026-09-30): compare de-versioned IDs,
    restricted to the index t2g.
    The viral set the no-manifest resolution would produce must equal the
    manifest's viral set, and every index gene must sit in one manifest set.
    Otherwise raise `BuildManifestContradiction` (a `ValueError`) in the
    `analysis` rule. `--gtf` is a union of panels, which is why the rule
    compares the resolved set and not the raw GTF set.
  - No manifest (every existing index: v2, final, max) → a logged warning,
    then the old path, byte-identical. `reference_strategy.py:704` now
    passes the all-virus GTF for both strategies.
- [~] `DEF-04` — *(2026-10-10: no shipped GTF carries a `whole_genome` gene any more, `panel_integrity.py` reports 0; the gene-level exclusion code is deferred: it only affects FASTA-only user references where the whole-genome gene is the virus's only unit, so it needs a design call, not plumbing.)* `whole_genome` tag in the catalogue and its exclusions (Q12).
  - 2026-10-12 (grill, user-confirmed): deferred past the exploratory public-dataset runs; no shipped reference is affected.
- [x] `DEF-05` — **done (f0dd5b3, c4e7d2e, 1841e5c; tracker reconciled 2026-10-10).** Two-step cell calling: STARsolo `GeneFull` into EmptyDrops,
  `--called-cells` override, and fail closed when neither exists (R2.8,
  under MECH-F).
- [x] `DEF-06` — **done 2026-10-10** (amendment in `docs/v3_counting_contract.md`; property tests in `tests/test_multimapping_properties.py`). Counting-contract amendment and property tests for the
  within-sibling allocation (R2.5). Must land before the D3 grid is frozen.
- [x] `DEF-07` — CAT-30 retraction notice beside the covid results, and the
  - 2026-10-13 **done** (user approved placement, Q70/Q78): `RETRACTION_NOTICE_CAT-30.md` (section 2 text of the draft, dated 2026-10-13,
    sha256 `90052893879a002f0d2a0e281081c86abc638a1ff6b5234ac97d7ceb36f55a7c`) placed beside the four `viral_summary.tsv` and next to
    `SURVEY_SUMMARY*.md` in the untracked `covid_viralscan/results/`; no published file edited. `claims/registry.json` `ttv-biological-signal`
    `reason` rewritten as drafted; `analysis/v3_artifact_inventory.tsv` row `claim-registry` re-pinned in the follow-up commit.
  - 2026-10-12 (grill, user-confirmed): user will sign off the notice wording (R3.5) before it lands; add the notice beside the covid results, never edit the published TSVs.
  claim marked retracted in `claims/registry.json` (Q11). Needs the user's
  sign-off on the wording (R3.5).
  - 2026-10-12: **draft only, nothing published or applied.** Draft notice and proposed registry diff:
    `docs/notices/CAT-30_retraction_DRAFT.md`. Findings that change the brief: the registry **already** holds the claim as
    `rejected` (`ttv-biological-signal`, `claims/registry.json:94-110`), so the diff is a one-field `reason` rewrite (its
    mechanism, "cDNA-only host-homology", is the one F-019 revised; re-pin `analysis/v3_artifact_inventory.tsv` when applied);
    the 27,098-read figure is the HPV16 raft (F-021), not covid, so the covid evidence cited is F-019 (1,269,176 reads, 90.1 %
    viral-best homopolymer, 9 clean all host-best); 1,167,103.0 is x213 only (both x216 dirs publish 1,605,631.0); 17.63/19.05 %
    is over all barcodes (56.64/64.77 % of called cells); the notice retracts the published magnitude, not the possible presence
    of anellovirus (F-019 `:31`, `ANELLO-PRIOR`). UNVERIFIED items are at the foot of the draft. Next: the user approves or
    edits the wording; then place `RETRACTION_NOTICE_CAT-30.md` beside the four `viral_summary.tsv` files and `SURVEY_SUMMARY*.md`,
    apply the registry diff, re-pin, and flip to `[x]`.
  - 2026-10-13 (Q70, user-confirmed): an agent places the notice (the wording was approved as drafted on 2026-10-12). **Not placed yet, deliberately:** the target,
    `covid_viralscan/results/` (4 sample dirs + `SURVEY_SUMMARY*.md`), is untracked and outside git (only `covid_viralscan/{README,RUNBOOK}.md` and `scripts/` are
    tracked), so a PR cannot carry it. Needed: (1) an explicit go to write `RETRACTION_NOTICE_CAT-30.md` into those 5 untracked locations of the shared checkout
    (no published file edited); (2) the date for the `2026-10-DD` placeholder; (3) a follow-up commit, after this branch is pushed, that applies the one-field
    registry `reason` rewrite and re-pins `analysis/v3_artifact_inventory.tsv` with the repo's own pin script (commit-bound pins must reference pushed commits).
- [~] `DEF-08` — lock `env_full`, and add a snakemake 9 dry-run to CI (Q9).
  - 2026-10-04: CI gains a `snakemake9-dryrun` job (`snakemake>=9,<10` +
    `test_snakefile_dag.py`). It has not run on a runner yet; watch for the
    `connection_pool` build issue. Locking env_full is still open.
  - 2026-10-10 (wave 4): the lock target is `viralscan_test_full`, not
    `benchmark_runs/legacy_v2_v3/env_full` (a stale kb 0.28.2 env). The explicit
    lock is `environment.explicit.txt` (see REL-03). `tests/test_snakefile_dag.py`
    passes 13/13 (incl. the 9 integration-marked dry-runs) with that env's
    snakemake 9.23.1. The CI-runner result of `snakemake9-dryrun` stays the user's.
- [ ] `CMP-06` (WP6B) — same-input comparison with evonk's original
  ViralScan 2.2.0 on the **same latest reference** as v3, so only the
  implementation differs (user request 2026-09-30).
- [x] `DEF-09` — release tagging: `3.0.0.devN` after G1, with outputs stamped
  `defaults_status=provisional` in `run_manifest.json` (Q1, R3.5).
  - Done 2026-10-05. `run_safety.DEFAULTS_STATUS = "provisional"` is
    stamped after hashing, as `completion_marker` is, so `--resume` still
    matches old manifests. Both schema copies list the enum
    `provisional`/`final`. The version is now `3.0.0.dev1` and the local tag
    is `v3.0.0.dev1`.
  - **Not pushed.** Pushing a `v*` tag runs `release.yml`, which publishes to
    PyPI and overwrites ghcr `:latest` with this dev build. That is the
    user's call.

### WP1C — Simplification pass (ponytail audit, new 2026-09-27)

A read-only audit on 2026-09-27 ("ponytail") covered all 58 k tracked Python
lines. The user chose groups 1–3.

- **Kept out of scope, as provenance of registered numbers:** the
  `analysis/multimap_profiling` scripts (`fast_profile*.py`,
  `interpret_cprofile.py`).
- **Not taken:** `extras/build_anello_table.py`. It is the only regenerator of
  a frozen TSV.
- **Rejected:** replacing stdlib `urllib` with `requests`. That runs the wrong
  way on the stdlib-first rule.

- [x] `SIMP-01` — delete verified-dead code:
  - `multimap.normalize_barcodes`: no caller, and it references an undefined
    `output`.
  - `host_filter._kallisto_filter`: unreachable.
  - `sensitivity.capture_reference`: an alias with no caller.
  - `menu._config_value`: only its test calls it.
  - `KbCountOutputs.bus_txt`: the legacy path, only its test reads it.

  Kept: `evidence.blast_identity`. Its only caller is the BLAST integration
  test that `BASE-04` cites, and that test cannot run here.
- [x] `SIMP-02` — collapse duplicated helpers in `src/` to one copy each:
  - `sha256` of a file: `src/` keeps `run_safety.sha256_file`.
  - the packaged-TSV path resolver: 4 copies.
  - the gzip-aware text opener: 3 copies.
  - the TSV writer in `evidence_run.py`: 6 copies.
  - the `--verbose`/`--quiet` argparse flags in `menu.py`: 7 copies.
  - `detection._sum_axis0/1`.

  **The `scripts/` copies of `sha256`, the TSV readers/writers and the atomic
  JSON writers stay duplicated on purpose.** The legacy v2/v3 and fresh-control
  tools copy `benchmark_v3_multimap.py`, `compare_legacy_v2_v3.py` and
  `run_fresh_control.py` into hash-frozen packets and run those copies alone
  (`tests/test_benchmark_v3_multimap.py`, `tests/test_freeze_fresh_control_packet.py`).
  An import of `governance_utils` would break a frozen copy. The audit missed
  this.
- [x] `SIMP-03` — dependencies:
  - Drop `pyfiglet`. `menu.py` already falls back to plain text.
  - Drop `seaborn`. Its 2 plots move to matplotlib.
  - Swap `enrichment._bh_adjust` for `scipy.stats.false_discovery_control`,
    with `scipy>=1.11` pinned and an equivalence check against the old loop.
  - Apply the change in `pyproject.toml`, `environment.yml` (inventoried, so it
    needs a re-pin), `conda-recipe/meta.yaml` and CI.

## WP2 — Install, lock, and artifact parity

Objective: make the pip tier honest and make conda, OCI, and Apptainer execute
the same tested build. Estimated remaining effort: 5-8 engineering days plus
runner time.

### WP2A — Pip and locked environment

- [x] `REL-01` — move full-workflow-only dependencies such as Snakemake out of
  - **Done 2026-10-04 (7d4bae0).** New `full` extra: snakemake, scanpy,
    scikit-learn, plotly, seaborn, pyfiglet. The core gains jsonschema, which
    was missing.
  - `doctor --profile pip` no longer needs snakemake, and a missing snakemake
    or kb names `pip install "viralscan[full]"`.
  - Clean-install check pending (REL-04/05, user-run network step).
  mandatory pip runtime dependencies; make `doctor --profile pip` check only
  Python/API/reporting/validation capabilities.
- [x] `REL-02` — define Python 3.11 on `linux-64` as the canonical full-workflow
  - **Done 2026-10-04 (7d4bae0)**, as a declaration only: Linux classifier,
    environment.yml marked canonical python=3.11 linux-64, docs. No lock was
    resolved; that is REL-03.
  toolchain; advertise other platforms only after the same workflow passes.
- [~] `REL-03` — generate and commit reproducible runtime/development lockfiles
  - 2026-10-10: `uv.lock` committed (`c5c4f9c`); `uv lock --check` passes. Runtime conda lock still open.
  - 2026-10-10 (wave 4): runtime conda lock `environment.explicit.txt` committed: 487
    conda packages (URL + md5, sorted), from `viralscan_test_full` (matches `environment.yml`:
    kb-python 0.30.2, snakemake 9.23.1, STAR 2.7.11b), written by
    `scripts/write_explicit_lock.py` from `conda-meta` because no conda binary is on PATH.
    Equivalent to `conda list --explicit --md5`; pip-only packages are not listed (as there).
    The file name follows the command already documented in `environment.tools.lock.txt`.
    `environment.tools.lock.txt` has no generator in the repo (hand-written sha256 rows) and was
    left as is. Not done: no `conda create --file` replay test (needs the network); the
    `environment-file` CI job still builds from `environment.yml` (proposed: switch it to
    `environment-file: environment.explicit.txt` once the user has replayed the lock). Status
    stays `[~]` until that replay and REL-06.
  - 2026-10-04 prep: `environment.tools.lock.txt` lists the verified binaries
    (version, build, sha256) and the user's lock commands (conda-lock or
    pixi, `uv lock`). Resolving the locks needs the network (user).
  and an exact external-tool version manifest.
- [~] `REL-04` — build wheel and sdist once, run `twine check`, inspect packaged
  assets, install each in a clean environment, and write `SHA256SUMS`.
  - 2026-10-02: built once from a `git archive` export with the codex
    env (`build --no-isolation --skip-dependency-check`, no network).
    - Outputs: wheel 587,519 B and sdist 663,361 B. `twine check` passes
      both.
    - Both install offline into a `--system-site-packages` venv and report
      `viralscan 3.0.0.dev0`.
    - `SHA256SUMS`: whl `9061229c…`, sdist `d718456a…`; `sha256sum -c` is OK.
    - The first `check_ship_scope` run failed: 9 members were missing from
      `config/public_ship_scope.json` (`virus_identity.py`, `sensitivity.py`,
      `index_exclusions.tsv`, `run_complete.schema.json`). They were added in
      this merge; the re-check is recorded below.
    - Still open: clean installs from an index (needs the network; REL-05)
      and the `release.yml` steps (REL-13).
  - Re-check after the allowlist fix (2026-10-02, HEAD a968a63 export):
    the build succeeds, `twine check` passes on both artifacts,
    "ship-scope validation passed", and `sha256sum -c` is OK.
  - 2026-10-12 repeat (HEAD `4cd916f`, `git archive` export in a clean temp dir): isolated `python -m build`
    (build 1.5.1; the network was available, so no `--no-isolation`) -> `viralscan-3.0.0.dev1` wheel 770,724 B and
    sdist 895,952 B. `twine check` PASSED both; `check_dist_size.py --max-mb 60` exit 0; `check_ship_scope.py`
    "ship-scope validation passed"; no `*.gtf`/`*.fa`/`*.fastq` member; `sha256sum -c` OK (wheel `f27f17c2...`,
    sdist `c582af61...`). Both install into separate Python 3.11 venvs (no system site packages, deps from PyPI;
    `pip check` clean, `viralscan --version` 3.0.0.dev1). Evidence: `analysis/release_checks/2026-10-12/`
    (`SHA256SUMS`, member lists, `rel04_build_and_venvs.txt`). Still `[~]`: to be repeated on the finished source
    (these hashes belong to `4cd916f`, not the release), and the `release.yml` steps are REL-13.
- [~] `REL-05` — run `doctor --profile pip` and `validate-run` against a packaged
  v3 fixture from both clean installations.
  - 2026-10-12: `scripts/release_package_smoke.py` (packaged-asset check, synthetic conserved v3 fixture, then
    `viralscan.menu doctor --profile pip` and `validate-run`) exits 0 from both clean installs (wheel venv and sdist
    venv, `viralscan` imported from site-packages): doctor `ok: true`, 11/11 schemas load, `validate-run`
    `ok: true`, `issues: []`. Log: `analysis/release_checks/2026-10-12/rel05_doctor_validate_run.log`. Source was
    `4cd916f`; repeat on the finished source. Not covered: the bundled-panel run from a clean install is still
    blocked by `REF-11`, and `doctor` (pip profile) checks no external tools by design.

### WP2B — Containers and parity

- [!] `REL-06` — Miniforge is digest-pinned; rebuild Docker from the committed
  lock and the exact tested wheel rather than from a mutable source install.
  - Blocked (2026-10-02): needs the committed lock (`REL-03`, `[ ]`) and the tested wheel (`REL-04`). No docker daemon here, and verification needs a CI build job.
- [ ] `REL-07` — publish a versioned OCI candidate, record its digest and tool
  manifest, and never assign `latest` to an RC.
- [ ] `REL-08` — build Apptainer/SIF from that OCI digest rather than performing
  a separate dependency solve; record the SIF SHA-256.
- [ ] `REL-09` — run the identical tiny workflow in locked conda, Docker, and
  Apptainer and compare validated molecule counts and deterministic output hashes.
- [ ] `REL-10` — save `analysis/release_parity/parity_report.json` with commands,
  versions, digests, hashes, and explained nondeterministic files.

### WP2C — Supply chain and release workflow

- [x] `PY310-01` — drop Python 3.9 (user decision 2026-10-12, taken from the `REL-11` lock audit: seven pins on the `python<3.10` branch had no fix supporting 3.9). Done 2026-10-12 (`codex/drop-py39`): `requires-python >=3.10`, 3.9 classifier and CI jobs removed, ruff `target-version py310` (`B905` ignored: `zip(strict=)` is allowed, not required; two `UP035` `Callable` imports moved to `collections.abc`), `uv.lock` regenerated (no `python<3.10` branches; `scripts/audit_locks.py`: 164 registry packages, 0 advisories, saved as `analysis/release_checks/2026-10-12/rel11_lock_audit_py310.txt`), README badge and docs updated. Option 2 (relax 3.9-only coding rules): `X | None` and `zip(strict=True)` are now fine everywhere. **Kept:** Snakemake `script:` modules must still not import `__future__`: that ban is Snakemake's preamble above the docstring, not a Python-version rule (`CLAUDE.md` Style says so). Left: the pytest matrix stays on `ubuntu-24.04` until 3.10-3.12 on Ubuntu 26.04 are verified (ubuntu-latest moves 2026-10-19); no CI gate was added for the lock audit (advisories are reported, accepting one is a maintainer decision).
- [~] `REL-11` — security CI exists; audit the locked product environment rather
  than the scanner job, add OCI scanning, a full SBOM, dependency/data licence
  report, and reviewed vulnerability exceptions.
  - 2026-10-12 (software part): `scripts/audit_locks.py` audits the locked environments against OSV
    (network, no credentials); output in `analysis/release_checks/2026-10-12/rel11_lock_audit.txt`.
    `uv.lock`: 201 PyPI-registry packages, all hashed, no git/url/path sources, `uv lock --check` passes.
    Found 10 vulnerable pins; bumped the three fixable ones (`urllib3` 2.7.0 -> 2.8.0, `gitpython` 3.1.59 -> 3.2.0,
    `virtualenv` 21.7.4 -> 21.14.5, dev/transitive only; `uv lock --upgrade-package`, no other pin moved).
    **Left for a maintainer decision (reviewed exceptions):** seven pins on the `python < 3.10` branch have no fixed
    release that supports 3.9 (pillow 11.3.0, requests 2.32.5, urllib3 2.6.3, setuptools 82.0.1, pytest 8.4.2,
    filelock 3.19.1, soupsieve 2.8.4; pillow/requests/urllib3 are runtime). Either drop Python 3.9 from
    `requires-python`/CI (3.9 is end-of-life) or record each as a reviewed exception; not decided here.
    `environment.explicit.txt` (487 conda packages, 429 conda-forge + 58 bioconda, md5 on every row, no duplicate
    names) has no name-matched OSV advisory (approximate: OSV has no conda ecosystem; urllib3 2.8.0, requests 2.34.2,
    pillow 12.3.0 are in it). Still open: OCI scan/full SBOM/licence report run only on a tagged release
    (needs docker), no exception file exists, and no scan was wired into CI as a gate because the 3.9 branch
    would turn it permanently red until the decision above is made.
- [x] `REL-23` — implement the grill-session-2 weekly non-blocking lock audit. Reuse `audit_locks.py`,
  pin Actions, retain per-lock logs, exit codes and skipped/failed step outcomes; keep it independent of
  CI/release gates. Software workflow preparation only, no dependency installs or lock updates.
  - 2026-10-13: `.github/workflows/lock-audit.yml` schedules Monday 07:23 UTC and supports manual dispatch;
    uv/conda audit jobs retain logs, lock hashes, exit codes and actual step outcomes as artifacts and summaries.
    Advisory or audit failure is tolerated, never a release/CI prerequisite; skipped/failed is not a clean result.
    House Action pins reused. Offline static/shell checks simulate exit 0/1/2 and skipped setup; 30 focused
    release tests and 2 tracked-tree ship-scope tests passed; ruff check/format, Action pins, data governance
    and diff whitespace passed. No live OSV or GitHub run.
  - 2026-10-13 fix: job-level `env` used `${{ runner.temp }}`, which GitHub rejects there ("workflow file issue"), so every push showed a failed lock-audit run. The receipt dir is now set by a first step via `$GITHUB_ENV`; `tests/test_release_static.py` 30 pass.
- [x] `REL-24` — draft `MAINTAINER_CHECKLIST.md` for `REL-15`, `RC-02`, `OPS-01` and `OPS-02`.
  Maintainers fill factual metadata, publishing settings and policy decisions; drafting grants no publication
  authority, reserves no identifier and makes no new support promise (grill-session-2 Q1-Q18).
  - 2026-10-13: root checklist draft links identity/contact, repository/publisher settings, policy approvals,
    audit review and candidate/archive ownership to their existing tracker items. All factual fields remain
    maintainer fill-in; `REL-15`, `RC-02`, `OPS-01` and `OPS-02` stay open. Both new files are outside the
    positive public ship-scope allowlists; no package/reference bytes or lock pins changed.
- [~] `REL-12` — pin GitHub Actions by commit SHA and generate provenance/
  attestations for wheel, sdist, OCI, SIF, locks, and checksums.
  - 2026-10-10: `check_actions_pinned.py` and the prerelease-keeps-`latest` rule are verified by mutation
    (`tests/test_release_static.py`). Pending (user): resolve action SHAs (network), review `uv.lock`, docker
    build, wire `check_versions`/`check_actions_pinned` into CI and `release_sif_definition.py` into `release.yml`.
  - 2026-10-12 (software part): every `uses:` in `ci.yml` and `release.yml` (40) is now a 40-hex commit SHA with a
    `# vX.Y.Z` comment, resolved with `gh api repos/<action>/commits/<tag>` (`scripts/check_actions_pinned.py`
    passes). `checkout` v7.0.1, `setup-python` v7.0.0, `upload-artifact` v7.0.2 and `download-artifact` v8.0.2 are
    first-party Node 24 majors (`action.yml` `using: node24` verified); the third-party actions keep their current
    major at its latest release (codecov v4.6.0, setup-micromamba v1.11.0, docker/setup-buildx v3.12.0,
    build-push v6.19.2, login v3.7.0, anchore sbom v0.24.3 / scan v7.4.2, pypi-publish v1.14.2; the `release/v1`
    branch is now an immutable SHA). The docker, codecov and micromamba actions at those majors are still Node 20:
    moving them to their Node 24 majors (codecov 7, micromamba 3, buildx 4, build-push 7, login 4) changes inputs and
    defaults, cannot be exercised without a CI run or a tag, and is left to the maintainer. `check_actions_pinned`
    and `check_versions` now run in the `ci.yml` lint job. Provenance: `pypa/gh-action-pypi-publish` v1.14.2 has
    PEP 740 `attestations` defaulting to true under trusted publishing, so wheel and sdist are attested at
    publication. GitHub artifact attestations for OCI/SIF/locks/checksums were NOT added (need `id-token` and
    `attestations: write` plus a tag run to verify). Still `[~]`.
  - 2026-10-12 CI hygiene: `ubuntu-latest` becomes Ubuntu 26.04 on 2026-10-19. `actions/python-versions` has no
    `linux26.04` build of Python 3.9 (3.10-3.12 have one), so the `pytest` matrix now runs on `ubuntu-24.04`
    (codecov condition updated). No apt packages, hard-coded tool paths or Python versions elsewhere in
    `ci.yml`/`release.yml`; the other jobs stay on `ubuntu-latest` (3.11/3.12 and micromamba/docker on 26.04 cannot
    be verified before the migration; watch the first run after 2026-10-19).
- [ ] `REL-13` — make release publication depend on green CI for the exact tagged
  SHA; build distributions once and make all downstream jobs consume them.
- [ ] `REL-14` — add protected TestPyPI/PyPI environments, stable approval, GitHub
  release assets, and separate RC/stable container tag behavior.
- [!] `REL-15` — maintainers must provide final author, licence-holder,
  maintainer/contact, ORCID, affiliation, support-window, and trusted-publisher
  metadata.

- [x] `REL-16` — **a standalone `kallisto` on `PATH` segfaults on a `kb ref`-built
  - **Done 2026-10-04 (d352927).** `validation.kb_tools()` parses `kb info` for the
    binaries the `kb` on PATH runs. `tool_path()` routes evidence (kallisto
    bus, bustools capture/sort/text), multimap (bustools correct/sort/text),
    doctor and the integration tests through them. Path, version and sha256
    are recorded in `run_manifest.json` (`tool_binaries`, outside the
    fingerprint) and in `evidence_manifest.json`.
  - Trap found: the test env's importable kb_python (0.29.5) differs from the
    `kb` CLI on PATH (0.30.2), so resolution follows the CLI, not the import.
  - `test_tiny_end_to_end` + `test_exact_lineage`: 24 passed (previously hung
    at 760 % CPU).
  index, and `viralscan evidence` picks it up.** Found 2026-10-03 while running
  `viralscan evidence` in `conda/envs/viralscan_bench`: exit `-11` (SIGSEGV) on
  the cat42b `panel.idx`. The two binaries that both report **version 0.52.0**:

  | Binary | Size | md5 | On a `kb ref` index |
  |---|---|---|---|
  | `$ENV/bin/kallisto` (conda) | 2,795,304 B | `c9cb3b19f7dbf09e85b70b45e1d50c97` | SIGSEGV |
  | `site-packages/kb_python/bins/linux/kallisto/kallisto` | 8,518,336 B | `695d3d418a11263adbc451dd63ea1233` | works |

  The version string cannot be used to tell them apart, so a preflight
  `shutil.which('kallisto')` check passes and the failure surfaces as a crash
  mid-run. Related and already recorded: the mirror case where a *hand-built*
  index fed to `kb count` makes the bundled kallisto spin at 787 % CPU for
  minutes on one read pair (`.living/learnings.md`, 2026-10-03) — same
  version-lock, opposite direction, and neither fails loudly.
  - Fix: resolve kallisto through `kb_python`'s bundled path (or an explicit
    `--kallisto`), not `PATH`, wherever ViralScan shells out to it directly;
    record the resolved binary's md5 in the run manifest so a parity run can
    detect the swap.
  - Until then: run `viralscan evidence` with kb's `bins/linux/kallisto` first on
    `PATH`.

`G2` passes when a clean lock-created environment plus Docker and Apptainer run
the same tiny workflow with identical validated counts; installable schemas and
assets are present; security/licence reports contain no unresolved actionable
finding; and the exact-SHA release workflow is dry-run verified without
publishing.

## WP3 — Preregister scientific validation

Objective: freeze outcome-independent decisions before producing new v3
scientific results. Estimated effort: 1-2 days plus reviewer sign-off.

- [x] `SCI-01` — create schema-valid `analysis/v3_validation/protocol.yaml` and
  `README.md` with hypotheses, supported scope, datasets, factors, seeds,
  chemistry, input/reference hashes, and exclusions.
- [x] `SCI-02` — freeze the barcode universe/cell-calling rule, feature
  intersection, primary unique-only cross-tool layer, secondary ambiguity-aware
  layers, denominators, and failure-reporting policy.
  Frozen subsection SHA-256:
  `699854b71169222e74d26c2119f1c9d7742d5962ac4ed8c02dfa9f288e3339dc`;
  dependent hypotheses/endpoints/exclusions digest:
  `fd818661a8d06f28c5c23e78bfb1c5b48b151023f5f749453260d023adbff333`.
- [~] `SCI-03` — freeze training/holdout partitions, evidence-tier calibration
  rules, metrics, limit-of-detection method, and biological-sample bootstrap unit.
  All 14 round-1 specification findings against these sections are resolved:
  largest-remainder apportionment replaces ceiling rounding, the sort key is
  pinned, an empty eligible grid has a declared failure outcome, the tie-breaker
  names its standard error, the bootstrap refuses degenerate small samples, and
  two denominators no longer let failed rows flatter a result. Both sections
  remain **`pending`**: `F1` and `F2` need `VAL-01` and `REF-08` factor levels,
  and dormant validator rules now refuse a freeze without them.
  - Blocked (2026-10-02): F1/F2 factor levels need `VAL-01`/`REF-08`, and the freeze follows `DEF-00`. A review round before `DEF-00` would be wasted.

- [~] `SCI-04` — enumerate every ViralScan, STARsolo, traditional alignment,
  Venus, Viral-Track, and VIRTUS row with exact environment/reference
  requirements. Twelve workflows over 78 rows. Both fairness defects are fixed:
  the only exact-truth dataset now reaches all three dedicated comparators, and
  each runs a native-published and a matched-accession-index arm so a tool
  difference is separable from a reference difference. Breadth and tuning
  asymmetries are disclosed rather than removed. Sections remain **`pending`**
  until `SCI-05` round 2 passes; `tool_environments` still blocks on `REL-03`.
  - Blocked (2026-10-02): `tool_environments` needs `REL-03`. Freeze after a clean `SCI-05` round.

- [~] `SCI-05` — obtain an independent protocol review, resolve findings, then
  record the protocol SHA-256 and Git SHA before outcome-generating runs. Three
  rounds are complete and recorded, each verifying the previous round's fixes by
  tamper experiment rather than by reading resolution notes. Round 1: 24 findings.
  Round 2: found two round-1 fixes cosmetic and reproduced a hole in the amendment
  rail. Round 3: confirmed those fixes real, then found the ledger itself was not
  append-only in fact. Rounds 4-9 progressively hardened the digest scope, the
  hash chain, and the git anchoring, each round finding second-order defects in
  the previous round's fix. Rounds 10 and 11 returned no blocker, and two
  reviewers independently judged the specification buildable. All rounds are
  dispositioned and every reproduced tamper now fails closed. `verdict` stays
  `does-not-pass`: no round has yet passed clean with zero open Majors, which is
  the bar for flipping it. `R11-F2` — the last open item, the ledger checker's
  fail-open branches — is closed as of 2026-07-27.

  **Round 12 (2026-10-13, independent review at `bb7a34f`): 1 Blocker, 7 Major, 5 Minor; does not pass G3.**
  Protocol side, branch `codex/sci05-r12-protocol`, ledger record `DEV-027` (one record, user Q79/Q80):
  - [x] `M1` — `frozen_inputs` now hashes whole dataset objects (asset `pending_reason` excepted),
    `protocol_version`, `seeds.derivation`, `partitions.applies_to_dataset_roles` and
    `generated_sample_structure`, and `freeze`. The reviewer's tamper (hhv6b `analysis_role`, version bump)
    now fails twice: in memory and, re-digested, against the ledger. Side effect: the G3 freeze and the
    holdout hand-over each need their own ledger record.
  - [x] `M5` — `synthetic_mixed_host_virus` is in W1 (`expected_row_count` 90 to 91; W2 left alone, E6 is
    ViralScan-primary). The validator requires every required dataset to have a row and a zero-row dataset to
    carry `no_workflow_rows_reason` (a stale or required-dataset marker fails). Endpoint-level coverage
    (each endpoint's datasets in a primary row) is not checked: left open.
  - [~] `M6` — outcome defined and a 13-entry hashed pre-freeze inventory added (`pre_freeze_runs`,
    `PFR-01..13`; repo evidence re-hashed, work-tree evidence recorded only). Rule enforced: an
    outcome-generating run may touch only a context-only dataset, and a pending dataset stays unseen. No
    recorded run touches `pbmc_10x_healthy_v3` or a synthetic dataset. DEV-020..026 justifications restated by
    `field_corrections`; **no `outcome_triggered` value was flipped.** User decisions
    (2026-10-13): (Q86) exposure on context-only data is not an outcome, so `outcome_triggered_rule` is now stated as
    decided and DEV-020..026 stay false; (Q85, then released) HNSCC `GSE164690` and tonsil-atlas `E-MTAB-13687` ran in DSR
    rounds 1 (2026-10-06) and 2 although grill 2 held them out until G3; their exact libraries and arms are itemised in
    `PFR-13` (round 2b ran none at inventory time). The user first had the rest cancelled and quarantined unread, then
    withdrew the hold-out ("just analyse all datasets and describe what we find; we can iterate always") and the tasks
    were resubmitted: they are ordinary development data, recorded in `DEV-027` as a decision, not a breach.
    **Consequence for real-data validation:** labels on real data are presumed, never known (user, 2026-10-13); truth
    exists only for synthetic plants, so any real-data validation claim needs other unseen datasets, or none.
    Still open: CAL-11 in the study repo must be gated on the V03 freeze commit (Codex session).
  - [x] `B1`/`M7`/Q80 — `influenza_gse336913` (GSE336913, V03) registered as a CONTEXT row with pending assets
    (`required: false`, `freeze_required: false`, no workflow row), not a scored T2 row. Scoring needs PRE-01, a
    5-prime read model, an outcome-independent read-structure-only GEX-experiment rule, influenza A in the frozen
    panel, a comparator strand statement, a whole-tier promotion rule and counts-only reporting. B1's tier
    decision (scored T2 role, E11/D26, `design-known`) is **deferred** (Q87, user 2026-10-13): V03 stays a context row.
  - [x] `m1` — CI append-only check compares a push to main with `github.event.before`.
  - [x] `m2` — schema hashed (`frozen_inputs.schema_sha256`); the two-copy identity test already existed
    (`test_packaged_schema_matches_the_canonical_schema`) and now also checks the pin. Residual: the schema pins
    harmonization digests as `const`, so every amendment of them is still a protocol + ledger + two-schema edit.
  - [x] `m3`/`m4` — stale blocker text and the DEF-00 doc header and section G.
  - [x] code-side request (PR #43, M4): `partitions.generated_sample_structure` gains `empty_barcodes_per_sample` 3000 and
    `empty_median_umi` 50 (DEF-00 B15, previously prose only), inside the `frozen_inputs` digest.
  - `M2`, `M3`, `M4` code paths and `m5` are the code-side PR (#43). If it also appends a deviation record, DEV
    numbering, the chain and `recorded_deviations` collide with `DEV-027`: renumber the later one.
    (PR #43 appended no record; the next number was `DEV-028`.)
  - 2026-10-13: `DEV-028` pastes the user-accepted Q83 choices before round 13 (see `SCI-05-Q83`) and records the
    real-data rule (Q103, Q106). Round 13 (independent agent plus Biomni) is next; G3 is not flipped.

  **Round 13 (2026-10-13, independent review of `origin/main` @ `de39995` plus Biomni `tsk_014nVmVBNhsmW0pqnbWIvB0b`; both
  `does-not-pass`: 0 Blocker, 8 + 6 Major).** Code side, branch `codex/sci05-r13-code`, ledger record `DEV-029` (M2 only):
  - [x] `M1` — `require_execution_allowed` (and through it `score_production` and `require_protocol_run_config`) refuses any
    phase except `training`/`holdout`; the draft phase validates cleanly, so it must not open a run gate.
  - [x] `M2` — `frozen_inputs_sha256` hashes whole `references` objects (asset `pending_reason` excepted), the `status`/`freeze_required` of every seed
    (plus `seeds.derivation`) and `planned_freeze_sections`, so `host_release`, `viral_panel`, `purpose`, `assets[].freeze_required`, the planned
    sections' `status`/`freeze_required` and `seeds.*.freeze_required` leave a trace. `frozen_inputs.limitation` reworded;
    `DEV-029` re-baselines the digest (before `9544411...`, after the live value). No schema edit (`covers` enum untouched).
  - [x] `M3` — `v3_partition.apportion` implements the frozen algorithm: quota `size * f`, floors, then largest remainders up to
    `ceil(N * f)` seats (exact `Fraction` arithmetic). The production shape (42 strata of 4) now gives 51 holdout seats, not 50.
    Frozen text unchanged. A one-library chemistry is training in `split_libraries`.
  - [~] `M4` — `--production` is refused with `--background-r1` (`production_design_errors`): the frozen host-only EmptyDrops
    anchor and the B15 empty rule (host-only UMI total in [10, `lower`)) need inputs the generator does not have (anchor
    barcodes, host-only UMI totals). **Left:** make them required inputs and drop the refusal; read-count ranking stays a
    golden-scale stand-in only.
  - [x] `M5` — `production_design_errors` requires every frozen abundance and chemistry level, every `training-and-holdout`
    dataset and `cycle` covariates with replicates >= levels. `--donor-libraries` calls `split_libraries` for a real background,
    so a donor library can be holdout (it was always training). The 10-libraries-per-chemistry floor is not checked: real
    background production is refused (M4), so there is no production path to check it on.
  - [x] Biomni `BM1` (Q122 a) — `precision_recall_f1(real_background=True, known_false=...)` reports calls outside planted truth as
    `n_unknown_origin`, not FP; only sibling/homolog confusions are FP; recall unchanged; `twin_excluded`/`cell_of` apply the D18/D19
    twin exclusion to E3/E4. `BM3` — `f1_from_counts` = `2TP/(2TP+FP+FN)`, precision `None` kept separate, `f1_interval` bootstraps
    the same statistic. `BM4` — `v3_lod.bca_interval` withholds constant, non-finite or bias-undefined intervals with a
    `reason_code` (LOD and F1 both). `BM5` — `score_by_truth_status` reads `truth_status` from the frozen registry: only `exact`
    datasets get false-call reports, all others a descriptive agreement report (`score_production` routes through it).
    `BM6` — `fit_probit` returns `not_converged` when BFGS fails and the gradient is not stationary; `lod_interval` keeps failed-fit
    reason codes (`failed_fit_statuses`). Quasi-complete separation is still undetected (documented ceiling).
  - [x] Minors: `empty_median_umi` read from the protocol (was 50), `D25` in `expected_denominator_fields` (m2), undeclared appended
    ledger records refused (m5), freeze sidecar opened and its protocol/training/threshold hashes re-checked (`freeze_sidecar_errors`;
    layout documented in its docstring, the file is still written at G3). Not done: m1 registry hygiene, m9, m10, m11, and the
    protocol-text Minors.
  - Tests: `tests/test_sci05_round13_code.py` (one per Major and per Biomni item). `test_v3_scorer.py` seat-total assertion and
    `test_protocol_code_honouring.py::test_m4_scorer_entry_scores_when_the_gate_opens` updated to the new rules.

  **Round 14 (2026-10-14, independent review of `origin/main` @ `1d7b842` plus Biomni; both `does-not-pass`: 0 Blocker, 6 + 5 Major).**
  One branch, `codex/sci05-r14-fixes`, ledger record `DEV-031` (the protocol text; code and tests in the same PR):
  - [x] `R14-M1` D6 — harmonization D6 is `not_applicable` at G3 (`count_layer` deferred, `missing_output` not-applicable-and-row-not-failed), consistent
    with `calibration.freeze_scope`; it stays in the E4 map so the map is complete, and no E4 row fails on it.
  - [x] `R14-M3`/Biomni 3 recall — `precision_recall_f1` computes recall over ALL planted truth (Q122); the twin exclusion and unknown-origin rule
    touch only precision-type counts, and the F1 pool takes the twin exclusion on calls and truth alike (stated in the protocol).
    The tests that preserved recall 1/3 were rewritten.
  - [x] `R14-M4`/Biomni 2 — `score_production(phase, libraries, viruses, *, run_config)`: no `protocol=` override (registry is the on-disk file the gate
    verified), `require_protocol_run_config` plus an EmptyDrops-only check (D2), a library's dataset read from its truth rows (a contradicting
    `dataset` flag is an error), `real_background_from_manifest`/`score_planted_calls` read the background from `twin_sample_id`, the
    negative-only domain is explicit (`domain: negative_libraries_only`; planted positives are refused), and the ungated report calls
    T2-T4 calls "unexpected", not "false". Libraries without manifest rows can still only declare their dataset (declared limitation).
  - [x] `R14-M5`/M4 — `--production` refuses `synthetic_factorial` without a real background; new blocker `real_background_anchor` (EmptyDrops
    anchor, B15 empties, a frozen donor-library list) separate from `simulator`. Production is therefore blocked end to end until the anchor lands.
  - [x] Biomni 1 digest — `frozen_inputs.status` is `const frozen` in both schema copies; the digest is recomputed and the declared-record binding applies
    whatever the status; a downgrade with a ledger baseline is an error.
  - [x] Biomni 4 `f1_interval` — `strata` is a required argument; units are resampled within stratum (donor libraries indivisible). Protocol
    `uncertainty.stratification_unit_rule`: chemistry and abundance for a generated-sample unit, chemistry alone for a donor-library unit (every library
    carries the full abundance grid). Point estimate, bootstrap and jackknife share one statistic.
  - [x] Biomni 5 / `V1` / `V2` BCa and LOD — no infinite replicates (a replicate that does not reach the LOD makes the interval not-estimable,
    `lod-not-reached`), no silent zero acceleration (`acceleration-undefined`), finite/ordered/positive-width bounds (`invalid-interval`); reason codes are the
    protocol's hyphenated vocabulary, now listed in `degenerate_interval_policy`. **Consequence:** with the frozen no-tolerance policy the LOD
    interval is often not estimable when the LOD sits near the top tested level (probe: 4 of 6 seeds); a tolerance would be a ledgered user choice.
  - [x] `R14-M6`/`P1` probit — `fit_probit` flags quasi-complete separation (the groups meet at one abundance; exact in one dimension) as
    `quasi_separation`; the agent's probe (0/4, 0/4, 2/4, 4/4, 4/4 at 1, 3, 10, 30, 100) is a test. Protocol `limit_of_detection.failed_fit_policy`.
    The old test that called a tie "a finite MLE" was wrong and is replaced.
  - [x] `R14-M2`/`Q135` S — the S population is the counted molecules of the called cells with at least one alignment to an index gene of v.
  - [x] Minors fixed: `m4` (`frozen_at` accepts a date or a date-time in every section), `V1`, `V2`, `V4` (stale "round 13" text; the append-only ledger keeps
    DEV-030's date and its dangling `review_reference`, and DEV-031 says so), `V5` (the 140-sample bound is stated as a unit count, replicates multiply samples).
  - Declared limitations (not defects of this branch): `m1` (sigma 0.5 hard-coded in the generator; apportion string not in `FROZEN_KEY_CONSUMERS`),
    `m6` (refresh `pre_freeze_runs` and the inventory at the G3 commit: it moves with every explore row until then), `m7` (the freeze sidecar does not exist
    until G3), `m10` (minimum three units is the user's frozen choice), `m11` (the explore ledger is cooperative), `m12` (accepted in round 13),
    `V3` (the adjudicator does not record the cell-calling method, the GTF hash is not checked against a frozen panel, gene ids are not virus-keyed;
    a single-row origin-spanning gene raises on malformed GTF; LOCK-01 A10 stays procedural, while `score_production` now requires an EmptyDrops
    run config), authenticity of truth-manifest rows for the dataset binding (the executor's pin).
  - Tests: `tests/test_sci05_round14_code.py`; `test_sci05_round13_code.py`, `test_v3_lod.py`, `test_protocol_code_honouring.py` updated to the new rules.

  **Round 15 (2026-10-14, independent agent review of `origin/main` @ `21572b1`: PASS with conditions, 0 Blocker, 0 Major, 6 new Minors N1-N6; Biomni: 1 residual Major).**
  One branch, `codex/sci05-r15-fixes`, ledger record `DEV-032` (protocol text, code and tests in the same PR); decisions in `.living/decisions.md` (Q139):
  - [x] `N1`/Q139 a — E8 keeps zero tolerance; a not-estimable interval is reported with its reason code AND a model-free bracket (`v3_lod.lod_bracket`): per virus and
    chemistry, from the highest tested level whose pooled detection is below 95 percent to the lowest level above it at or above 95 percent, no fit; lower bound only
    when no level reaches 95 percent. Attached by `lod_interval` whenever the interval is not `ok`. Protocol: `limit_of_detection.failed_fit_policy` and E8 `reporting`.
    Tests use the review's probe design (10 libraries x 7 levels x 2 samples per level, probit truth).
  - [x] Biomni Major (LOD overflow) — `lod95` classifies against the tested range in log space before `10 ** x`; a NaN log LOD is status `nonfinite_lod`. Tests: tiny positive
    slope on the point and the bootstrap path.
  - [x] `N2` — new training blocker `run_provenance_binding` (SCI-04 executor: bind each library's dataset to the run manifest, record and check the cell-calling method per
    library (V3), hash the GTF or panel behind D1; truth-row authenticity). Code: `score_by_truth_status` refuses an exact-dataset library without truth rows
    (`exact-dataset-without-truth-rows`), which also covers Biomni's rowless-library point. `score_production`'s EmptyDrops check remains an attestation by the caller's run
    config, not a check of the libraries' runs (stated in the blocker).
  - [x] `N3` — `failed-fit` (now naming `below_range` and `nonfinite_lod` replicates) and `point-estimate-failed` are in `failed_fit_policy` and `degenerate_interval_policy`.
  - [x] `N4` — `lod_interval(strata=...)` is required (a missing label is an error), as for `f1_interval`.
  - [~] `N5` — fixed: `precision_recall_f1` returns `n_truth` and `n_recalled` beside the twin-filtered `tp`/`fn`; an exact dataset requires tier T1. Declared, not changed: the public
    `real_background` keyword of `precision_recall_f1` is a pure metric parameter used by the tests; protocol runs derive it from the manifest through `score_planted_calls`
    (renaming it would not make it unforgeable).
  - [x] `N6` — the G3 signing procedure below, and `scripts/check_g3_diff.py` (test `tests/test_sci05_round15_code.py`).
  - Declared limitations unchanged: `m1` (register sigma 0.5 and the apportion string at VAL-09), `m10`, `m11`, `m12`. `m6` and `m7` are discharged in the signing commit (steps 3 and 4).
  - [x] `G3-PREP` (2026-10-14, user-approved, `codex/g3-prep-tests`) — the G3 content (draft PR #61) failed step 5: 10 tests asserted that the live protocol is still draft, and
    the signing commit may not touch tests or docs. Eight `tests/test_validation_protocol.py` tests now check invariants that hold before and after the freeze (a section is
    frozen only with its planned entry, `frozen_at` and the `partitions_metrics` removal), and the rule probes build their own pending input (`_pending`). The workflow-comparison
    generator no longer reads `workflow_matrix.status` (frozen is not executed), so the page reads the same after G3. Verified on main's draft protocol and with PR #61's protocol,
    ledger and sidecar overlaid: 155/155 both ways. Then main goes into #61, pins are refreshed, and step 5 is re-run.

  **G3 signing procedure (user signs; agents never sign or freeze).** The reviewed draft and the signed commit differ by construction, so the signing commit must differ from the
  reviewed `origin/main` commit `R` (the merge commit of the DEV-032 PR, parent of the signing commit) in exactly the fields below and nothing else:
  1. Start from `R` on a branch off `origin/main`; `R` must have a clean round (agent PASS, Biomni no open Major) recorded here. The sidecar records both `R` and the round-15
     reviewed SHA `21572b1`, and states that `R` versus `21572b1` is DEV-032 and nothing else (user confirms).
  2. Edit only these `protocol.yaml` paths (the allowlist of `scripts/check_g3_diff.py`): `factors[*].levels|status|pending_reason` (the eight pending factors);
     `seeds.generation.status|value|pending_reason`; `{partitions,calibration,defaults_selection,workflow_matrix,failure_and_deviation_reporting}.status|frozen_at|contract_sha256`
     (`defaults_selection` has no digest); `planned_freeze_sections.{partitions,calibration_metrics_lod,workflow_matrix,failure_and_deviation_reporting,independent_review}.status`;
     `freeze.freeze_manifest|frozen_at`; `pre_freeze_runs.inventory|inventory_sha256` (m6: refresh from the explore ledgers, recording `tonsil_atlas` and `sfl_tonsil` as seen if they
     become candidate backgrounds); `frozen_inputs.contract_sha256` and an append-only `frozen_inputs.recorded_deviations`; and the removal of exactly the training blockers
     `partitions_metrics`, `independent_review_round_2`, `independent_review`. Not changed: top-level `status` (stays `draft`), `freeze.training_execution_allowed|holdout_execution_allowed`
     (false), `training_allowed|holdout_allowed` (false), every other blocker (`real_background_anchor`, `simulator`, `tier_cutoff_measures`, `execution_gate_unimplemented`,
     `run_provenance_binding`, `data_hashes`, `tool_environments`, ...), `harmonization`, the schema (both copies), `thresholds_sha256`/`training_results_sha256` (null), and all code.
  3. Digests recomputed with the repo's own functions in dependency order; one ledger record per digest scope that moved (append-only); then refresh the artifact inventory and claim
     registry pins (content commit, push, pin, push, re-pin the claim registry, push).
  4. Write `analysis/v3_validation/freeze_sidecar.yaml` (m7): `protocol_file_sha256` over the committed protocol bytes, reviewer identities and findings, `R` and `21572b1`, the
     explore-ledger sha256, and the user's sign-off and date (the user's words). The sidecar is bound by `freeze_sidecar_errors`; the sign-off is a digest-bound part of the signing commit.
  5. Mechanical checks, all required and recorded in the commit body: `python scripts/check_g3_diff.py R HEAD` (fails on any path or file outside the lists above; allowed files: the protocol,
     the ledger, the sidecar, the inventory, the claim registry, PLAN.md, `.living/decisions.md`); `python scripts/check_ledger_append_only.py --revision R`; `cmp` of the two schema copies;
     `validate_v3_protocol.py` (valid at the draft phase, and `--phase training` must still fail closed with the remaining blockers), `validate_artifact_inventory.py`,
     `validate_claim_registry.py`; ruff, mypy and the full test suite. Any failure means the commit is not signed.

  Standing limitation, established by the 2026-07-27 code review: every one of
  these rounds read `protocol.yaml`. None could have caught a frozen value that
  no call site consumes, because that is invisible from the protocol side and
  from the output side alike. Protocol review does not substitute for verifying
  that the code honours the contract.
  - 2026-10-13 (Q66, Q53, user-confirmed): round 12 = an independent agent review plus a Biomni external review; **the user signs G3** and owns this
    sign-off. The `GSE336913` influenza T2 amendment (row `T2-AMEND-01`) is reviewed in the same round.

  **Round 12 (2026-10-13, independent review of `origin/main` @ `bb7a34f`; verdict
  `does-not-pass`: 1 Blocker scoped to the V03 amendment, 7 Major, 5 Minor).** The
  review was the first per-key audit of the standing limitation above. Dispositions
  (user approved the fixes 2026-10-13, Q79; branch `codex/sci05-r12-code-honouring`
  takes the code-honouring ones, the protocol-side ones need `protocol.yaml`,
  `deviations.yaml` and both schema copies and belong to another branch):
  - `[x]` **M2** (frozen `seeds.cell_calling` vs the shipped EmptyDrops seed 100):
    `validate_v3_protocol.run_config_errors` / `require_protocol_run_config` read
    the seed named by `host_only_emptydrops.seed_source`, plus `fdr`, `lower` and
    `niters`, and refuse a mismatch after `require_execution_allowed`. The non-protocol
    CLI/`defaults.py` default stays 100 (wrong only under a protocol); a test pins the
    other three defaults to the protocol. The guard's real call site is the future
    SCI-04 executor's two `required_call_sites`; nothing calls it today.
  - `[x]` **M3** (frozen values with no consumer): `FROZEN_KEY_CONSUMERS` registry
    (`consumed` / `pinned` / `not_yet_consumed: <stage>`) plus
    `frozen_key_coverage_errors`, which fails when a numeric frozen leaf under `seeds`,
    `harmonization`, `partitions` or `calibration` has no entry, an entry is stale, or
    a named consumer does not mention the key. `v3_lod.REPLICATES` / `MINIMUM_UNITS`
    are pinned to `calibration.uncertainty`. Honest status: `seeds.root` (VAL-01;
    `seeds.generation` is pending), `seeds.evidence_sampling` and `seeds.bootstrap`
    (SCI-04 executor), `barcode_normalization`, `holdout_evaluations_allowed` and the
    calibration `grid` are `not_yet_consumed`; `seeds.root` was not wired into
    `Panel.rng` because that changes every generated byte.
  - `[~]` **M4** (gate had no call site): `--production` now runs
    `require_execution_allowed("training")` first (re-raised as `ValueError`), then
    enforces `cells_per_sample`, `reads_per_cell` and `seeds.split` from the protocol;
    `v3_score.score_production` is the gated scorer entry; a test lists every script
    that touches the protocol/scorer/generator and requires a gate call or a named
    exemption. **Partial:** the B15 empty-barcode count (3,000, median 50) is prose in
    `harmonization.cell_universe.synthetic_rule`; the guard reads a structured
    `partitions.generated_sample_structure.empty_barcodes_per_sample` and refuses a
    synthetic background while it is absent. **Needs a protocol edit** (add that field,
    digest scope included). Still no workflow-row executor, so the gate is wired on the
    generator and scorer only.
  - `[x]` **m5** (`knee` selectable): protocol runs accept only `emptydrops` or
    `external` (`auto` falls back silently, `none` has no anchor); the CLI keeps `knee`
    as the warned sensitivity-only option (SW-23, DEF-00 A7).
  - Not in this branch: B1, M1, M5, M6, M7 and m1, m2, m4 (protocol/schema/CI/doc
    side). m3 (stale `execution_readiness` blocker text) is digested protocol text, so
    it needs a ledgered DEV record, not code.

`G3` passes when the protocol validates against its schema, has no unresolved
review finding, is hashed, and there is no recorded pre-freeze outcome on scored
or protected pending-row datasets (attestation wording, user, 2026-10-13).

- [~] `SCI-05-R13P` — SCI-05 round 13, protocol side (`DEV-030`, 2026-10-13, base `b43fedf`; Q118-Q122)
  - Round 13 (agent: 8 Majors, 12 Minors; Biomni: 6 Majors) returned does-not-pass. This row covers the protocol-side
    findings; M1-M5, BM5 and BM6 and the code Minors are on `codex/sci05-r13-code`.
  - `[x]` **M6/BM2 (Q119):** `calibration.freeze_scope`. G3 freezes the probable tier (cumulative tiers) only. Strong,
    AUPRC and the cell score (D6 not applicable; `cell_score_requirement` deferred) and M8 (removed) are outside the
    freeze. The measures behind the three grid axes are defined by **Q135/Q136 (user, 2026-10-14)** in
    `calibration.freeze_scope`: a cell is `probable` when m >= M and F >= f and S >= s*; m = summed host-conservative X
    per cell (fractions kept); F = index genes of the virus with summed called-cell X >= 1 (not > 0), per sample x virus;
    S = called-cell virus_best / (virus_best + host_best + tie) counted molecules, zero denominator not-estimable; F and
    S are broadcast from sample x virus to its cells. Their implementation is pending, so training blocker
    `tier_cutoff_measures` stays until the code exists (`evidence_tiers.py` needs a per-cell loop and float molecules;
    S needs an evidence run per sample and virus).
  - `[x]` **M7 (Q120):** `E11_real_data_descriptive` (D10, D12) and `H8`; E10, the anellovirus rule and the T4 gate
    reworded (restates and supersedes DEV-028 UC-C5); `docs/plans/real-data-sensitivity-specificity.md` carries a
    supersession note.
  - `[x]` **M8 (Q121):** generator crosses abundance by homology per donor library on real backgrounds, with a test;
    protocol text states it, and that low_complexity is not planted into a real background.
  - `[x]` **BM1 (Q122), BM3, BM4:** `derived_metric_rules` (unknown origin, count-form F1), calibration objective
    (selection unit, strata, micro pooling), `degenerate_interval_policy`.
  - `[x]` Minors: m2 (D25 expected fields), m3 (stale prose), m9 (SE of the best eligible point), E5/D25 "not assessed",
    G3 attestation wording; Q118 recorded as a known limitation in `principal_risks`.
  - Not done here (code or later): m1, m4-m8, m10-m12, BM5, BM6 and the Biomni validator/sidecar minor.
  - Merge note: rebased after #51 (code side, `DEV-029`): this record is `DEV-030`, digests recomputed on top of the
    post-#51 chain and the pin chain redone.

- [x] `SCI-05-Q83` — Q83 pending-section drafts: user choices accepted and pasted (`DEV-028`, 2026-10-13)
  (`docs/plans/2026-10-13-protocol-pending-sections-draft.md`, 2026-10-13, base `88f02a1`).
  - Covers `partitions`, `calibration`, `defaults_selection`, `workflow_matrix`,
    `failure_and_deviation_reporting`, `independent_review`. Five already had full
    bodies, so each section gets a delta, a source table (86 rows, classed user / agent /
    none) and a USER CHOICE list; `independent_review` has no schema key, so it gets the
    planned-section entry and a sidecar skeleton.
  - Blocks verified on scratch copies only (`protocol.yaml` and the schema are untouched):
    the pending-phase edits and the full G3 flip validate at the draft phase.
  - **Answered 2026-10-13: all 25 USER CHOICE items accepted as recommended** (UC-P4: `pbmc_10x_healthy_v3`
    stays external-only; UC-I5: the G3 date and sign-off wording stay a placeholder for the user).
    Pasted on `codex/protocol-fill-uc` as ledger record **`DEV-028`**, re-based onto the post-#46 protocol
    (the draft's blocks predate `DEV-027`), every digest recomputed with the repo's own functions. The G3
    freeze flip is **not** done: no section status, factor, `freeze.*` or top-level status moved.
    - Applied to protocol text: UC-D1 (6-candidate joint grid) and UC-D4 (no chemistry-specific default at G3)
      in `defaults_selection`; UC-W1/W2 (version policy, pins by one REL-03 amendment) and UC-W3 (W2 extended
      from 6 to all 13 W1 datasets, `expected_row_count` 91 to 98, `docs/workflow_comparison.md`
      regenerated) in `workflow_matrix`; UC-F1 (re-based: points at `pre_freeze_runs.outcome_definition`, explore
      ledger sha256 goes in the sidecar) and UC-F2 (DSR-13 retry rule for protocol rows) in
      `failure_and_deviation_reporting`; UC-I1..I4 in `planned_freeze_sections.independent_review.description`.
    - Decision only, no protocol text (`DEV-028` `user_decisions_recorded`): UC-P1..P4, UC-C1..C6, UC-D2, UC-D3,
      UC-I5, UC-X1. `partitions` and `calibration` are untouched.
    - Not pasted: the draft's outcome definition (`DEV-027` already has one) and its "no outcome-free deviation if a
      ledger lists a run" clause (contradicts Q86); the `g3-freeze` stanza; the sidecar file
      (`analysis/v3_validation/freeze_sidecar.yaml`, written at G3; `freeze.freeze_manifest` stays null).
    - **DEV-028, real-data rule (Q103, Q106):** real-data validation arms report agreement with PRESUMED labels
      descriptively only, with unit counts; no accuracy, sensitivity, specificity, precision, recall or PPV claim from
      any dataset whose `truth_status` is not `exact`; accuracy claims come only from synthetic truth (planted
      panels). Recorded as a new item in `failure_and_deviation_reporting.prohibited`. No endpoint, hypothesis, dataset
      or schema edit, so the harmonization digests and both schema copies do not move; H2/E2 on synthetic exact
      negatives keep the zero-call rule (UC-C4).
    - **Follow-ups (not blockers in the protocol):** (1) UC-C1: the tier rule that emits probable/strong does not exist
      in shipped code; before G3 it must expose exactly the three grid knobs (molecules, features, support) so the frozen
      45-point grid is consumable. (2) UC-D2: measure the multimap rerun cost on the pilot panel and re-check the grid
      cost against Q10d (about 3.9k to 10.2k of the 11k shared core-h). (3) UC-C5: the claim registry has no tier or
      gate field (its `policy` is a string), so the T2/T3 (>= 3 donor-distinct libraries) and T4 (>= 10 independent
      libraries) minimums are recorded in `DEV-028` and the Track B lock, not in `claims/registry.json`.
    - **Next:** SCI-05 round 13 (independent agent plus Biomni, Q82), then the user signs G3 and the freeze commit
      flips the six sections, eight factors and `freeze.freeze_manifest`, with one ledger record per frozen digest scope.

## WP4 — Freeze production references

Objective: make reference contents reproducible and calibrate host-homology
safeguards without holdout leakage. Estimated effort: 3-5 engineering days plus
about 8 cluster hours per full GRCh38 build.

### WP4A — Profiles and provenance
  - Blocked (2026-10-02): run the next review round after the `DEF-00` amendment lands.

- [x] `REF-01` — profile names and expanded anellovirus behavior exist; the
  expanded panel is now the **default** rather than opt-in *(correction
  (reconciled 2026-10-04): the CLI default is off, see `ANDET-06`/`ANDET-07`; only the
  profile behaviour landed)*, and its gene IDs
  reach detection. Accession lists are frozen for `anellovirus-representative`
  and `anellovirus-expanded` via `src/viralscan/data/anellovirus_accessions.tsv`
  (2,042 accessions, 2,042 unique, 0 duplicates). `curated` and
  `broad-discovery` still carry **no** frozen accession list and remain
  label-only — see the new `REF-13`.
  - **2026-09-26 evidence for the default flip.** Measured 31-mer coverage of
    the bundled 20-genome RefSeq TTV panel against the 2,042 real human
    anellovirus genomes (`anellovirus.fa`): median **0.00 %**, and **85.8 %** of
    genomes share *zero* 31-mers with the panel. Positive control: EBV
    NC_007605.1 against itself = 100.0 % over 144,283 31-mers, so the method is
    sound. Per genus, zero-coverage genomes: Betatorquevirus 98.4 % (n=1,542),
    Gammatorquevirus/Samektorquevirus/Hetorquevirus/Gyrovirus/Memtorquevirus
    100 %. Leave-one-out capture for a *novel* strain: bundled panel 0.04 %
    (P(90 bp fragment captured)=0.024) vs expanded panel 20.5 % (P=1.0000).
    Cost of expanding: 5,994,773 bp total, 4,888,291 distinct 31-mers, only
    **1.21x** k-mer space inflation, 91.1 % of 31-mers unique to one genome,
    median genome 20.5 % redundant. CD-HIT-style downsizing does not help:
    20 genomes retain 1.1 % of the k-mer space, 204 retain 11.3 %, 2,042 retain
    100 %. Conclusion: keep all 2,042.
- [~] `REF-02` — fail-closed fetches and manifests exist; complete accession
  version, taxonomy, snapshot, retrieval date, SHA-256, length, licence, cluster,
  representative status, rationale, and missing-accession fields.
  - Blocked (2026-10-02): the licence fields need `REF-05`. Wiring the manifest into the production builder is `REF-03` (`[ ]`).
- [x] `MASK-01` — (found 2026-10-07 while starting `REF-03`) `build-ref --anellovirus` masking was a
  no-op or a hard failure. `_run_dustmasker` ran one window and wrote a **soft** (lowercase) mask; the
  k-mer gate and kallisto upper-case, so with `mask=True` (gate limit 0) any panel with a poly-A
  stretch failed its own gate (reproduced on a synthetic record), and a panel that happened to pass was
  not masked. The real panel was masked by `viral_panel_max_2026-09-28/03_mask_gate.py`.
  - **Done 2026-10-07.** `build_reference.mask_low_complexity` is that recipe: dustmasker windows 64 and
    30 (level 30), union of lowercase positions to `N`, then every N-free k-mer window the gate's
    classifier flags is masked in full. Input is upper-cased first, so an already soft-masked FASTA is
    not turned into `N`. The dedicated builder now calls it; `_run_dustmasker` keeps its signature plus
    `window`. Parity: 60 of 60 sampled prototype genomes (of 1,912) come out byte-identical to the
    prototype's masked records. Regression test `test_masking_a_poly_a_stretch_passes_the_gate_...`.
  - The targeted pass is a per-window Python scan (`ponytail:`), fine for curated panels, slow on the full
    2,343-genome panel. `REF-03` then wires the combined path to it.
- [~] `REF-03` — apply identical masking, duplicate-ID/sequence validation, and
  manifest generation to dedicated and combined build paths.
  - **Masking and gate done 2026-10-07.** `build_combined_reference` takes `mask=True` and runs
    `mask_low_complexity` (`MASK-01`) on the viral FASTA, then the same `_enforce_low_complexity_gate`
    the dedicated builder uses (extracted; limit 0 masked, 2 / 5 % with `--no-mask`). Host cDNA is not
    masked or gated (it is full of poly-A); N keeps lengths, so GTF coordinates are unchanged.
    `build-ref` preflights `dustmasker` before any download unless `--no-mask`; `--no-mask` now applies
    to every build path. **Default `build-ref` therefore needs BLAST+ and produces an N-masked viral
    panel.** Duplicate-ID/sequence validation and `write_reference_manifest` were already shared.
  - 2026-10-10 (wave 4): the old "Still open" note is stale. Since `f4320d2`, `menu._build_kb_ref` calls
    `prepare_reference_inputs(mask=True)` (dustmasker, gate, per-sequence `reference_manifest.json`). Remaining
    limitation: that manifest records `host_species="unspecified"` for `--reference`, because a user-supplied
    FASTA/GTF carries no species.
- [~] `REF-04` — make frozen inputs rebuild byte-identical panel FASTA/GTF/t2g
  contents and save a reproducibility audit.
  - 2026-10-10 (wave 4): comparison tooling done. `scripts/compare_reference_reproducibility.py A.json B.json`
    (function `compare_reproducibility`) exits 0 when two `reference_reproducibility.json` files agree on
    `content_sha256` and every per-file digest, else 1 with one `field: a != b` line per difference. The
    binary `.idx` is not claimed. Test: `tests/test_compare_reference_reproducibility.py`. **Left (cluster
    task):** rebuild the real frozen-input panel twice and save the comparison as the audit.
- [ ] `REF-05` — replace vague source-data licence text with reviewed terms for
  every redistributed or fetched reference source.
- [!] `REF-11` — publish the viral annotation panel archive and register its
  Zenodo DOI. `src/viralscan/data_fetch.py` pins
  `VIRAL_DATA_DOI = "10.5281/zenodo.20112332"`, but that identifier is **not
  registered**: `https://zenodo.org/api/records/20112332` returns
  `{"status": 404, "message": "The persistent identifier is not registered."}`
  and `https://doi.org/10.5281/zenodo.20112332` also returns 404, while an
  unrelated third-party DOI referenced elsewhere in the repo resolves normally.
  The 195 GTFs remain in `src/viralscan/data/` in the source tree but are absent
  from the installed package, so `viralscan data fetch` — and therefore every
  bundled-panel run from a clean install — fails for all users. Blocks `REL-05`,
  `DOC-05`, and `G2`; currently blocking the `GOV-06` attempt-3 cache pin.
  Owner action: publish the archive and register the DOI, or correct the pinned
  record identifier.
  - 2026-10-20: `scripts/make_panel_archive.py` prepared (not run, not shipped; not in
    `MANIFEST.in` or ship scope). Builds a deterministic zip + `PANEL_MANIFEST.json` +
    `SHA256SUMS` + draft Zenodo metadata after checking the GTFs against
    `analysis/panel_expansion/gtf_corpus_manifest.tsv`. Publishing stays the user's.

### WP4B — Full-genome competition and anellovirus

- [~] `REF-06` — the public `--genome-dlist` path and raw annotations exist;
  build the production curated human-plus-virus index with a checksum-pinned
  GRCh38 D-list.
  - **2026-10-03: grill approved, build and measure.**
    - Finding: the shipped cat42b has no genome D-list. kb's default D-list
      is the positional `combined.fa`, which is 465,769 host cDNA records
      plus the 2,343 viral genomes; that gives 632,261 D-list k-mers.
    - Literature: Sullivan et al., NAR 2025 (PMID 39657125) for the
      distinguishing flanking k-mer mechanism; Luebbert et al., Nat
      Biotechnol 2026 (PMID 40263451) use a host genome + transcriptome
      D-list, but for translated search.
    - Biomni review stalled with no output, so an internal review was used.
  - **Design (single variable):** `viral_ref_cat42d/` is cat42b's inputs
    and pinned code d01043c, with D-list = cat42b `combined.fa` + GRCh38-2024-A
    `genome.fa`, sha256 pinned in `dlist.sha256`.
    - Gates are one-sided: HPV16 may lose at most 1 %, EBV/HSV-1 at most
      5 %, HHV-2 must still be flagged as bleed, and per-gene losses are
      checked. Any increase is investigated.
    - Swap rule (default): cat42d replaces cat42b if every gate passes.
  - **2026-10-05: gates assessed. All pass, so the swap rule applies and
    cat42d replaces cat42b as the latest reference.**
    - Script: `scripts/ref06_compare_cat42b_d.py`, run over
      `multimap_evidence.tsv`.
    - The index has 3,094,720 D-list k-mers, against 632,261 before.
      `combined.fa` is byte-identical to cat42b's.

    | run | cat42b unique / est | cat42d unique / est | gate |
    |---|---|---|---|
    | HPV16 SRR19537341 | 4,769 / 7,725 | 4,771 / 7,726 | ≤ 1 % loss: pass, per-gene order unchanged |
    | EBV SRR12682296 | 184,485 / 747,532 | 184,825 / 749,314 (+0.2 %) | ≤ 5 % loss: pass |
    | HSV-1 SRR8315713 | 14,275 / 23,187 | 14,413 / 23,419 (+1.0 %) | ≤ 5 % loss: pass; HHV-2 still `possible_em_bleed` |
    | covid x213, `--strand reverse` | 74 / 110 across 12 viruses | 4 / 15 | specificity, not a gate |

    - **Per-gene, genes with ≥ 20 unique molecules:**
      - EBV, 46 genes: median 0.0 %. The worst gene is `BaRF1.1`, 1,001 →
        950 (−5.1 %). It is already dropped as a marker, and RRM2 homology is
        a plausible cause, though not checked.
      - HSV-1, 20 genes: worst −1.4 %.
    - **Host-derived background removed.** Spurious low-level unique calls
      went to 0: HHV-6 `p23` (11 → 0 in HPV16 rafts; 18 → 0 in covid), plus
      Monkeypox, MOCV, HPV9, HHV-1 and EBV `EBNA-2` in non-EBV data.
    - **The increases were investigated.** Unresolved molecules fall on
      every run (EBV 1.55 M → 1.39 M, HSV-1 1.10 M → 0.96 M). This is
      consistent with UMIs whose reads conflicted between host-genome and
      viral ECs now resolving to viral, once the host-genome reads are
      D-listed out. Consistent, not proven per molecule.
    - **Side effect: the called cells move.** emptyDrops, with the same
      seed, calls fewer cells:
      - HPV16 rafts: 15,615 → 12,186 (−22 %) and 12,251 → 10,165;
      - EBV and HSV-1: about −1 %.

      The lost barcodes are low-count: median 130 UMIs, against 2,608 in
      the kept barcodes. Kept barcodes lose about 8.5 % of their host UMIs,
      because intronic and intergenic reads are now discarded. HPV16
      infected-called cells fall only 2,006 → 1,964, so `pct_infected` rises
      from 12.8 % to 16.1 %. A protocol consequence: a reference swap moves a
      kb-derived cell anchor. The anchor must therefore be frozen *after* the
      reference (G4), or come from a reference-independent source.
    - Still open for `[x]`: the curated reference manifest (`REF-07`) and
      the G4 freeze. cat42d is the latest reference, not yet the frozen
      production one.
  - **Jobs:**
    - build 25695872;
    - regressions 25695873 (GSE189670) and 25695874 (EBV/HSV-1), both
      afterok;
    - covid x213 `--strand reverse` on cat42b (25695875) vs cat42d
      (25695876).
  - **Side check, tiny HPV16 index (kallisto 0.50.1):**
    - 90-nt reads crossing a CDS end into unmodelled sequence: 0/21
      aligned with the whole-viral-genome D-list (11 D-list k-mers), 1/21
      without it. Fully-inside controls: 18/21 both ways.
    - So these reads are lost by pseudoalignment itself, not by the D-list.
      Removing viral sequence from the D-list is not needed.
    - This supports F-022: gene models must reach the mRNA ends.
- [~] `REF-07` — save the reference manifest, host-homology/low-complexity table,
  index/t2g/GTF, commands, versions, checksums, and build resource accounting.
  - **2026-10-05: the host-homology table is in progress** (only that part; the
    manifest and accounting are untouched). It exists to freeze VAL-01's
    `host_virus_homology` levels, which VAL-01 D7 put ahead of generator code.
    - Script: `scripts/ref07_host_homology_table.py`, output
      `viral_ref_cat42d/ref07/host_homology_loci.tsv`.
    - **Per locus, not per genome.** `build_reference.measure_host_homology`
      keeps the single longest hit per viral record, which answers "is this
      genome host-homologous". VAL-03 needs every host locus, because each is a
      *source* of challenge reads.
    - **Dead end: genome as the minimap2 reference.** `-x asm10/asm20` against
      GRCh38 returns **zero** alignments for all 2,343 genomes — the asm presets
      want long colinear blocks, while viral/host homology is short diverged
      patches — and costs 78 s plus 13 GB just to index. Indexing the 11 MB
      viral panel and streaming the genome as the query uses 0.5 GB, runs ~70x
      faster, and finds the real signal (chr21 alone: 743 hits, HHV-6A/6B
      telomeric and HSV-2).
    - Levels are chosen from **training** clusters only, with low-complexity
      hits excluded, because `low_complexity` is its own VAL-01 factor. Hits are
      clustered by their *viral* interval and whole clusters are assigned, so a
      multi-copy family cannot straddle the split (D4).
    - Sensitivity control: the genomes whose spurious calls the REF-06 D-list
      removed must each appear in the table
      (`scripts/ref07_level_proposal.py`).
  - **2026-10-05: the table is built, and the result changes `VAL-01`
    (finding F-023).** 5,971 alignments over 51 of 2,343 genomes, in 89
    clusters (union-find over viral and host intervals; 70 training, 19
    holdout). An earlier draft said 112, which was the count before the union-find
    fix.
    - The sensitivity control passes 8/8, so the result is not an alignment
      failure.
    - **98.8 % of the viral bases that align to GRCh38 are low-complexity**,
      while dustmasker covers only **4.50 %** of the panel's bases. Exactly
      **1** of 89 clusters carries a 31-mer outside low-complexity, and **0**
      of those are also exonic. Only **12** hits touch a transcript 3' window,
      which is what a 3' chemistry samples.
    - **The sensitivity control is synthetic, not the genome-level one.**
      Dust-free human exons mutated to a known identity and run through the
      same command recover 20/20 at 0.95 and 0.90, 16/20 at 0.85, 12/20 at
      0.80 (`scripts/ref07_recovery_control.py`). Genuine non-repeat homology
      at 90-95 % would have been found without exception. The earlier 8/8
      genome-level check is too weak to carry this: HHV-6 passes on 1,439
      telomeric hits whatever happens at `p23`.
    - Clusters join on a shared viral segment **or** a shared host locus. The
      first version merged viral intervals per record only, which put
      `chr9:109765432-109765517` (MPXV and Cowpox) and the HHV-6A/6B/7
      telomeres into both partitions. After the fix: 89 clusters and 0 of
      1,721 1-kb host windows in both partitions; realized holdout share
      0.213 against a nominal 0.30, recorded rather than corrected because
      these are templates, not samples.
    - The top contributors are repeat-driven: HHV-6B (1,751 hits, 100 %
      low-complexity), HHV-6A (1,439), Cowpox (690), HPV19 (388), MPXV (364).
      HHV-6's are its telomeric `TTAGGG` arrays against human telomeres.
    - **Consequence:** `host_virus_homology` and `low_complexity` are not
      separable factors, and the identity/aligned-length bands the factor
      describes have no material. `DEF-00` section F puts three options to the
      user; the recommendation is two levels, `none` and `repeat_homology`,
      with each challenge read carrying its own low-complexity fraction in the
      truth manifest. That makes the strata 7 x 2 x 3 = 42 and the floor 168
      samples plus twins.
    - `REF-09` consequently tests mostly low-complexity reads, which `DEF-01`'s
      filter already targets. The protocol should say so rather than let the
      row read as independent evidence.
    - Still open for `[x]`: the reference manifest, the index/t2g/GTF copies,
      and the build resource accounting.
  - 2026-10-10 (wave 4): the dead end above is now fixed in the package too: `build_reference.measure_host_homology`
    uses the viral-panel-as-reference orientation (see `ANDET-02`).
  - 2026-10-10 (wave 4): resource accounting and tool versions are now in `build_receipt`. `_run_kb_ref` records
    `peak_rss_kib` (`getrusage(RUSAGE_CHILDREN).ru_maxrss`, KiB on Linux) with `peak_rss_scope`: it is the
    high-water mark over **every** child process the build has reaped (dustmasker, minimap2 ... as well as kb),
    not kb alone. `record_reference_build` writes `build_receipt.tools` (path, version, SHA-256 for kb, kallisto,
    bustools, dustmasker, cd-hit-est; `null` when absent) via `validation.tool_provenance` plus version probes for
    the three tools `kb info` does not report. `menu._build_kb_ref` now passes `resources=` and `index=` (it
    previously discarded both). Left for `[x]`: the curated-panel manifest itself, the index/t2g/GTF copies, and
    the real build.
- [ ] `REF-08` — calibrate homology/complexity exclusion thresholds using only
  preregistered training controls and freeze `thresholds.json`.
- [ ] `REF-09` — prove planted human-homology reads cannot reach probable/strong
  evidence on holdout; retain raw measurements and all excluded calls.
  - 2026-10-13 (Q96, `codex/evidence-tiers`): **depends on `TIER-01`** (the probable/strong code, now built) and on
    `SCI-03` freezing (a) a `strong` rule and (b) the exact-fragment support-fraction definition. Until a supplied
    calibration result exists, `viralscan.evidence_tiers.assign_tier` returns `uncalibrated`/`no_tier`, so the test has
    nothing to measure yet.
- [~] `TIER-02` — **new, 2026-10-14 (Q135, Q141; `codex/tier-measures`).** The three measures behind the frozen `probable`
  grid, per `calibration.freeze_scope`. Software only; no data read, no calibration run, `protocol.yaml` untouched (so no
  digest, no DEV record, no collision with `DEV-032` on `codex/sci05-r15-fixes`).
  - `[x]` `src/viralscan/tier_measures.py` (in both lists of `config/public_ship_scope.json`): `molecules_per_cell` (m,
    summed `adata.X` over the index genes, fractions kept), `features_supported` (F over called cells, >= 1 molecule, with
    `n_index_genes`; `cutoff_unattainable`), `support_fraction` (S from an evidence directory), `measure_virus` and
    `probable_cells` (m >= M and F >= f and S >= s*, inclusive, F and S broadcast). `evidence_tiers.TierEvidence` takes a
    float molecule count and `_meets` uses the same relative slack as `utils.infected_masks`.
  - `[x]` F wired into `detection.compute_stats` and `viral_summary.tsv` (`n_features_supported`, `n_index_genes`, called
    cells only, `adata.X`); documented in `docs/output_reference.md`. m is already `viral_molecules_total_est` /
    `per_cell_df` (`is_called_cell`).
  - `[x]` S input is explicit: `support_fraction(evidence_dir, virus, called_barcodes)` reads `molecule_verdicts.tsv.gz`
    from a `viralscan evidence` directory. Population: counted molecules in called cells with a viral alignment (host-only
    molecules excluded and counted); a repeated `(cb, ub)` collapses; zero denominator is `no_counted_molecules`; missing
    directory, missing verdict file and competitor mode are separate reason codes. A manifest naming another target, a hash
    mismatch, or no counted barcode in the matrix raise `EvidenceMismatch` (never a guess).
  - Dependencies, read-only: `molecule_verdict.py` file format as on `main`; the `(cb, ub)` key is the corrected CB-UMI
    that `codex/umi-dedup-q137` audits. The verdicts carry no per-gene coordinates, so "aligned to an index gene" is the
    evidence run's extraction plus a present viral alignment score.
  - Decision: the F >= 1 and m >= M comparisons carry `INFECTED_RELATIVE_TOLERANCE` (1/3 + 1/3 + 1/3 in float32 counts as
    one molecule); sums are float64.
  - `[ ]` Blocker `tier_cutoff_measures` stays: nothing runs `viralscan evidence` for every sample x scored virus
    including negatives (S would be not-estimable on every negative, so the zero-false-call constraint would pass
    vacuously), and competitor mode writes no verdicts. Lifting needs that driver plus a ledgered amendment (next free
    number after `DEV-032`: `DEV-033`), digest refresh and pin chain (r15 has merged; rebase on its digests). Blocker text is partly stale
    (the measure code now exists; the driver does not).

- [~] `TIER-01` — **new, 2026-10-13 (Q96, Q104).** `src/viralscan/evidence_tiers.py`: probable/strong tier code built to
  the frozen `calibration.threshold_search` contract (UC-C1: the 45-point grid is a contract on this code; UC-C6: no
  eligible grid point = `no_tier`). Software only; no calibration was run and no data read.
  - Done (`codex/evidence-tiers`): `grid_from_protocol` parses the grid prose (strict; exactly the three sets in
    order; `protocol.yaml` is untouched, so no digest moves); `select_cutoffs` takes each grid point's training F1, its
    resampling-unit bootstrap SE and the host-only / planted-homology false-call counts and applies the contract
    (eligible = zero calls in both classes; within one SE of the best F1 pick the lexicographic maximum of
    `(molecules, features, support)`; none eligible = `failed`); `calibration_from_dict` refuses an off-grid point and a
    `strong` that is not >= `probable` on every cutoff; `assign_tier` returns `uncalibrated` (no result), `no_tier`
    (failed), else `not_detected`/`candidate`/`probable`/`strong`. `validate_v3_protocol.FROZEN_KEY_CONSUMERS` now lists
    the grid as `consumed: viralscan.evidence_tiers.grid_from_protocol` (was `not_yet_consumed`); the module is in both
    wheel and sdist lists of `config/public_ship_scope.json`. Pins refreshed in the same PR (`output-reference-doc`,
    `validation-protocol-command`, the registry artifact and claim, then the `claim-registry` row).
    Q104 mapping table (`STUDY_LEVEL_TO_PACKAGE_TIER`, documented in `docs/output_reference.md#evidence-tiers`): `L0`
    candidate, `L1` probable, `L2` strong, `L3` strong; `scripts/adjudicate_candidates.py` appends one column,
    `package_tier_analogue` (an analogue by position, not an equivalence and not a package-assigned tier). The two
    scales stay separate: the study adds read-level evidence on top.
  - **Left / needs the user (contract gaps, nothing invented):** (1) the grid and objective define only
    *probable-or-stronger*; `strong` is an optional caller-supplied grid point nested above `probable`, never emitted
    without one; (2) "exact-fragment support fraction" and "distinct supporting molecules per cell / features per
    sample" are inputs (`TierEvidence`) that nothing in the package measures yet; comparisons are inclusive (`>=`);
    (3) the tie-breaker's "one standard error" is read as the SE of the best eligible point; (4) the grid is a prose
    string, so a structured form is a freeze-time item (it would change `calibration.contract_sha256`); (5) not wired
    into `detection`, the Snakefile or any output column (`evidence_tier` keeps its `candidate_*` vocabulary); (6) the
    F1/SE computation and the bootstrap run belong to the `SCI-03`/`SCI-04` executor.
- [!] `REF-10` — an orthogonally confirmed anellovirus-positive sample is needed
  for real sensitivity claims; without it, ship screening support only and label
  every result accordingly.
- [~] `REF-13` — **new, 2026-09-26.** Detection-side reference visibility and
  - 2026-10-10 (wave 4 audit): a fresh `scripts/panel_integrity.py` run reports **2,626** errors (2,469 transcripts with no exon,
    103 HHV-6B CDS rows with no `transcript_id`, 54 reused `unassigned_transcript_N`) plus 47 catalogue gaps; the drop from 2,629 is the
    removed duplicate `NC_002076.2` GTF. **5** genes (not 6) have neither CDS nor exon (`TTV3_gp6`, `D1P65_gp2`, `EPSTEIN_HHV4_BARF0`,
    `HUM_ADENO_HAdVC_gp12`, `HUM_PARVO_B19V_gp4`). kb-python 0.30.2 (ngs_tools 1.8.6) ignores CDS rows, so each CDS-only gene is
    indexed as one gene-span transcript. The normaliser design (build-time, intron removal only for single-protein spliced CDS,
    flatfile joins for the 13 collapsed HHV-6B joins) is in `~/.claude/plans/thoroughly-plan-to-finish-mellow-starlight.md`.
  - 2026-10-10 (wave 4): **normaliser implemented, branch `codex/wave4-ref13`** (user decision: build-time, shipped GTFs untouched).
    `src/viralscan/gtf_normalise.py::normalise_viral_gtf` (idempotent) is called from `scripts/build_bundled_panel_ref.py` (the
    `combined.gtf` assembly), `scripts/panel_integrity.py` (default view; `--raw` gates on the shipped files) and
    `extras/build_gene_programs.py::parse_gtf`. Rules: gene span = one transcript (ID = `gene_id`) so cDNA is byte-identical; introns
    removed only for one-protein CDS and only for gaps >= 11 bp (the influenza A PA-X 1 nt gap is a frameshift); the 13 collapsed HHV-6B
    joins come from the cached `.gb`; adenovirus-style multi-protein genes keep the span; HHV-6B terminal-repeat copies get one
    transcript each (`G`, `G-c2`); HBV P/S origin-wrapping genes keep ngs_tools' last-gene-row behaviour (`wrap_last_row`, cDNA
    unchanged, because `ngs_tools` concatenates exons in genomic order); UL111A also keeps its span transcript (`-span`).
    **Deviation from the brief:** B19V `HUM_PARVO_unassigned_gene_1` is not given a new gene row; its 11 kDa CDS rows span exactly
    the gene-only gene `HUM_PARVO_B19V_gp4` (4890-5174), so the locus is already indexed and a second gene row would duplicate it.
    The CDS rows are re-keyed onto `gp4`. Tests: `tests/test_gtf_normalise.py`. QC artefact and counts: next commit.
  - 2026-10-10 (wave 4, QC): `scripts/ref13_normaliser_qc.py` -> `analysis/panel_expansion/ref13_normaliser_qc.tsv` (4,147 genes: 4,040
    cDNA-identical, 107 changed: 101 spliced incl. UL111A, 6 HHV-6B repeat-copy genes). `panel_integrity` on the normalised view:
    **7** errors (was 2,626 raw; 7 are CDS transcript IDs with no exon row in genes that already have exons), `docs/panel_integrity.md` updated.
    viral-only kallisto index: targets 4,206 -> 4,213, k-mers 4,653,849 -> 4,638,071, max EC 17 -> 17, 0 discarded.
    UL111A intron retention not verified (Jenkins 2004 abstract says only "different splicing pattern"): span transcript kept.
    **Open:** `gene_programs.tsv` not regenerated (overlap groups change: EBV 5 -> 2 independent groups); HBV P/S first part still
    unindexed; the unpinned `build-reference` path was not wired.
  - 2026-10-10 (wave 3): astrovirus `NC_001943.1` column-3 `Non structural gene` fixed to `CDS`
    (`scripts/fix_thin_gtfs.py`; integrity count 2,628 -> 2,629 because those two CDS rows are now visible and
    exon-less like the other legacy records). **Duplicate `NC_002076.2` resolved (user decision, 2026-10-10):** the builder consumes
    `Alphatorquevirus homin1.gtf` (gene IDs `TTVgp1..`); the local gitignored `Torque_teno_virus_NC_002076.gtf` (`TTV_TTVgp1..`) was removed
    and its row dropped from `gtf_corpus_manifest.tsv` (393 GTFs). `anellovirus_genes.tsv` keeps `NC_002076.2_TTVgp1..` (accession-prefixed;
    not checked against the builder's GTF gene-ID handling, so confirm at the WP5 build). **Not done, needs the user:** the exon-less CDS normaliser (changes counts); the 6 genes with neither CDS nor exon; the `curated` list.
  - 2026-10-04 (1166f42): duplicate FASTA IDs now fail loudly in
    `package_starsolo_viral_references.validate_pair` and in
    `evidence.write_competitive_fasta`, before samtools/STAR; this was the
    July `NC_002076.2` crash.
    - The packaged `references/starsolo/.../viral_genome.fa` on disk still
      carries the duplicate and must be rebuilt.
    - GTF-side dedup in `copy_many` is not implemented.
  name resolution. Two defects made most of the panel inert while the runs
  still looked clean:
  - `SENS-01` (**fixed**) `scripts/analysis.py` globbed only the packaged panel
    directory, but the expanded anellovirus GTFs are materialized into the
    *built index* by `build-reference`, so 2,022 of 2,042 genomes (91 %) were
    countable and never reportable. `anellovirus.candidate_gene_ids()` now
    derives the `{accession}_geneN` IDs the builder emits, and
    `--anellovirus-gene-ids` (default on) controls it. Over-inclusion is the
    safe direction: `detect_genes` only reports IDs that are real columns.
  - `SENS-02` (**fixed**) 49.2 % of gene IDs in the panel the covid runs
    actually used (`references/starsolo/.../viral_genome.gtf`, 4,650 genes)
    resolved to no virus name, so each became its own row in
    `viral_summary.tsv` — the covid run published 9 of 17 rows as bare gene IDs
    (`HHV1gp00p39`, `CeHV2gUL24`, `MPXV_gp132`), which silently broke
    `accession_breadth` (1.0 by construction), sibling cross-mapping and
    `eve_risk` for exactly the herpesvirus calls. Fixed by adding the
    underscore-delimited tokens to `VIRUS_NAME_MAP` and a new
    `VIRUS_GENE_ID_ALIASES` tier for concatenated schemes. 49.2 % -> 13.5 %,
    **0 regressions** across 7,443 real gene IDs. The boundary rule was *not*
    weakened (the `AICHIX`/`BORF1`/`BUNYAMW` guards are now regression tests).
  - Still open under this row: `curated`/`broad-discovery` accession lists;
    per-GTF SHA-256 + retrieval dates for the 195 packaged GTFs (no manifest
    row exists for any of them; *2026-10-03:* `scripts/write_gtf_manifest.py
    [-o OUT] [--check OUT]` now writes and checks a per-GTF sha256/size
    manifest. It hashes 324 GTFs on this machine. `source_date` is always
    `unknown`, because retrieval dates are not recoverable. No manifest is
    committed, since the GTFs are gitignored, and the check is not wired into
    CI); the 6 genes with neither CDS nor exon; the
    duplicated `NC_002076.2`; the malformed astrovirus feature column; and the
    2,520/2,692 gene IDs (93.6 %) with no `exon` record, which STARsolo
    comparators cannot count at all.

`G4` passes when panel contents reproduce byte-for-byte, manifests validate,
duplicates are absent, GRCh38 competition is operational, and holdout host-
homology negatives cannot become probable/strong calls.

### WP4B2 — Host-filtering design and reporting (new 2026-09-26)

- [x] `HOST-01` — pin every STAR parameter the host filter depends on
  (`STAR_FILTER_ARGS` in `scripts/host_filter.py`) and record each in
  `host_filter_audit.tsv`. The command previously set **no** alignment or filter
  options, inheriting whatever the installed STAR defaulted to. Two defaults
  were wrong for viral subtraction: `outFilterMismatchNmax 0` rejects any read
  with one host mismatch, so paralogues and allele variants escaped as
  "unmapped" and reached the viral index; `outFilterMultimapNmax 1` reports a
  multi-mapping read as unmapped, which is the dominant false-positive route
  for host repeats and EVEs. Now `4` / `20`, with `outFilterMatchNminOverLread`
  lowered `0.66` -> `0.9`. Closes the parameter half of G8 step 4; the
  per-fragment removal-reason half remains open (see `HOST-02`).
- [x] `HOST-02` — record `pct_retained` in the audit so a re-run that retains a
  different fraction is attributable to a parameter change. Per-fragment removal
  reasons are still **not** recoverable: `--outSAMtype None` discards the SAM and
  `fragment_lineage.tsv.gz` logs retained reads only.
- [x] `HOST-03` — add `pct_infected_comparable`, a strategy-independent
  denominator (absolute 200-molecule host-UMI floor intersected with the called
  set). `pct_infected_called` is **not** comparable across host-filtering
  strategies: on one covid PBMC sample, called cells fell 143,243 -> 28,921 while
  Alphatorquevirus UMI fell 1,167,103 -> 57,715, so `pct_infected_called` *rose*
  from 56.64 % to 62.89 % and inverted the comparison. `pct_infected_called`
  stays the within-run primary.
- [x] `HOST-04` — document the measured three-way host-control comparison in
  `docs/cli_reference.md` and **retract** the "~4x more sensitive" two-step claim
  in `BENCHMARK_COMPARISON.md`. That figure came from a 1M-read subsample of an
  implementation that no longer exists (kallisto CB-UMI-wide subtraction, since
  removed); at full depth on the same sample the arms differ by **3 %**
  (1,479,894 vs 1,434,619 UMI), and the original two-step arm was `blocked`, not
  completed. `--genome-dlist` is now documented as the *weakest* option: it cut
  `p_unique` 2.1 % -> 0.6 % (3.5x) to remove only 14 % of the anellovirus
  artifact, because a k-mer D-list cannot see diverged host sequence.
- [ ] `HOST-05` — **open.** The sensitivity cost of the v3 STARsolo path is still
  unmeasured. It blocks `MS-02`/`CMP-01`–`CMP-03` and the preregistered
  `D15`/`D16` endpoints; no number should be quoted for it until then.

## WP4C — Detection sensitivity and negative-result claims (new 2026-09-26)

Objective: make a negative result self-describing. Added after an audit that
found `limit_of_detection` present only as an unrun endpoint in
`schemas/v3/validation_protocol.schema.json` and in no code path, while
`min_counts`/`min_genes` gated only the UMAP and never detection. A run that
detected nothing carried no depth caveat, so "nothing there" and "did not look
hard enough" were indistinguishable from the output.

- [x] `SENS-01` — `src/viralscan/sensitivity.py`: Poisson detection
  probability, depth-only LOD95, and the exact-match k-mer capture curve
  `P = 1 - (1 - (1-d)^31)^(L-30)`. Depth is the sum of the count matrix, not
  raw reads, because only quantified molecules can be detected.
- [x] `SENS-02` — every run writes `results/sensitivity.tsv` (LOD95, band,
  `depth_sufficient`, `capture_measured`, `informative_negative`) and states the
  limit in `summary.txt` and `report.html`. `viral_summary.tsv` gains
  `pct_infected_comparable`.
- [x] `SENS-03` — `informative_negative` requires depth **and** a *measured*
  capture term, so it is `false` on almost every run by construction. This is the
  load-bearing design choice: depth was ample in every real run here (LOD95
  0.0003–0.0056 per 10k host UMI, all `informative`), while the covid samples
  called SARS-CoV-2 = 0 at 21.6 M quantified molecules. Reference capture — not
  depth — is the binding limit, and capture is unmeasurable without a control.
- [x] `SENS-04` — positive control: `--positive-control-gene` +
  `--positive-control-molecules` (required together), `--require-positive-control`
  to fail closed, `results/positive_control.json`, and bisection inversion of the
  capture curve to an implied divergence. `failed` (control invisible),
  `over-recovered` (not spike-in-specific) and `gene-not-in-reference` are all
  distinct, reported states.
- [x] `SENS-05` — LOD semantics calibrated against molecule-level downsampling
  of the bundled EBV LCL run (103,145,071 molecules, 1,636,934 EBV): P(detect)
  stayed 1.0000 down to 1,270 downsampled reads and first reached 0 at 127,
  i.e. the observed floor is the Poisson floor. `fragment_capture` cross-checked
  at 90 bp and 150 bp for 5–30 % divergence.
- [~] `SENS-06` — **open.** `SENS-01`–`SENS-05` make the limit *reportable*;
  they do not make a negative *certifiable* in practice, because no shipped
  workflow plants a control. Closing this needs either a spike-in recipe in
  `docs/vignettes/` or integration with `VAL-01`'s generator, after which
  `E8`/`D17` LOD95 can be estimated by the preregistered probit fit rather than
  reported as an analytic floor.
- [~] `SENS-CORR-02` — exact substitution-only model. **Done (2026-10-06):**
  `sensitivity.fragment_capture_exact` (O(Lk) DP, absorbing-state mass, validated
  inputs) now *is* the capture model: `fragment_capture` is a thin deprecated alias
  that delegates to it, so every caller moved to the corrected curve (90 bp, k=31:
  5 % -> 0.708, 10 % -> 0.254, 15 % -> 0.063, 20 % -> 0.013; the retired
  overlapping-window formula gave 1.0 / 0.903 / 0.323 / 0.058). It is labelled a
  "substitution-only heuristic" in docstrings, notes, log text and
  `positive_control.json` (`implied_divergence_note`) and never feeds
  `capture_measured`. `detection._implied_divergence` is now
  `substitution_model_implied_divergence` (no alias kept: no other consumer; the
  JSON key `implied_divergence` is unchanged). The "0.32 at 15 %, 0.06 at 20 %"
  prose in `sensitivity.py`/`detection.py` is now computed from the exact model
  (`sensitivity.capture_cliff_text`), and `faq.md`/`CHANGELOG.md` are updated.
  Pinned tests were re-pinned with a comment saying why (0.9031 -> 0.2537 etc.;
  the old formula is kept as a test-only reference to document the overstatement).
  **Left:** `docs/output_reference.md` still says "0.32 at 15 % divergence and 0.06
  at 20 %" (replace with the exact values + heuristic label, and document the new
  `positive_control.json` fields). It is pinned by `analysis/v3_artifact_inventory.tsv`
  row 16 (git blob at `2ceea6b`), so editing it needs the governance re-pin after a
  commit; I did not edit it, to keep `test_artifact_inventory` green. Also stale:
  `SESSION_RESUME.md:230` (not mine). `substitution_only_kmer_survival` alias name
  from the audit was not added (no consumer).
- [~] `SENS-CORR-01` — scoped capture. **Done (2026-10-06):** `CaptureScope`,
  `PositiveControl`, `scoped_capture` are wired. `build_sensitivity_table` applies
  the control's capture only to rows in scope (`detection.row_capture`); every other
  row, including undetected indexed viruses, gets no capture term, so
  `capture_measured` is False and the LOD is the labelled depth-only floor (the
  `capture` column is empty, user decision 2026-10-06; the internal 1.0 is only the depth-only floor).
  A control is built from the raw ratio and only when `status == "measured"`, so the
  clamped over-recovered 1.0 certifies nothing. `positive_control.json` gains
  `scope`, `target`, `certified_targets`; `certifies_negatives` is true only when
  `certified_targets` is non-empty, and `summary.txt`/the negative statement now
  name only those targets. The report is written before the `--require-positive-control`
  raise, so a withheld run keeps its diagnostics. CLI `--positive-control-scope
  {exact_sequence,virus_key,panel_mechanics}` and `--positive-control-virus-key`
  flow through `menu.py` -> `RunConfig` (so `createconfig`/`config.yaml`) and are in
  `run_safety.omit_when_unset`, so the fingerprint of runs that do not set them is
  unchanged. A legacy control with no scope maps to `panel_mechanics` and warns once;
  it no longer certifies any virus (tests that asserted the old behaviour were
  changed with a comment). **Decision:** `virus_key` is accepted by the parser but
  rejected by `RunConfig` until VAL-RA-CAL provides an approved transfer calibration
  (audit acceptance: "virus_key without valid calibration is rejected").
  **Left:** `exact_sequence` identity is declared by flag, not verified: no versioned
  accession/sequence digest or pinned evidence tying it to the control and the
  current reference; rows are matched by the group key in the sensitivity table, not
  resolved through `VirusIdentityTable` (legacy prefix-named rows are matched by
  display name); the immutable capture-measurement record (digests, count layer,
  evidence digest); a schema-valid control-failure receipt (the JSON is now written
  before the raise but has no schema); schema/packaging decision for
  `positive_control.json`; SENS-CORR-03 rendering (HTML always prints the unmeasured
  statement; the empty `capture` column for uncertified rows is done, 2026-10-06).
- [ ] `SENS-CAL-01` — **TODO** empirical capture calibration against ViralScan's own
  index. Take real fragments (88-91 bp R2 geometry) from held-out viruses, mutate
  them on a divergence grid (0-15 %, i.i.d. substitutions first, then clustered),
  run `kb count` against the real panel + host index, compare recovery with
  `fragment_capture_exact`, and fit one effective-length (or effective-divergence)
  parameter. Training-only: it feeds VAL-RA / VAL-RA-CAL and is never scored on the
  untouched holdout. Why: Luebbert et al. 2025 Fig 1c (EBOV, i.i.d. substitutions,
  figure reading +/-3 pts, read length unverified) puts real recall near 55 % at 4.4 %
  divergence where the exact DP at L=88 gives ~76 %; the DP also cannot see
  panel/host competition. Depends: SENS-CORR-02, a frozen panel+host index.
- [x] `VAL-RA-DECISION-01` — **decision (2026-10-06, user):** the held-out genome
  partition policy is whole-taxon exclusion (not single accessions), near-duplicate
  rule ANI >= 95 % (CAMI II), and at least two sister taxa kept in the reference
  (Kraken 2 design); sparse strata are reported "not estimable", and no per-distance-
  stratum probit LoD is claimed unless pre-registered. Frozen before any recovery is
  looked at. Precedents: `viralscan_lshwcp_spec/07_LITERATURE_REVIEW_CAPTURE.md` Q3.
  Still needs the reviewed protocol amendment under `CLAIM-DECISION-01` before G3.

## WP4D — Gene-programme inference, layer 2 (new 2026-09-26)

Objective: for viruses layer 1 detected, distinguish latent from productive
expression per cell. Added after measuring that a per-gene comparison cannot do
this at all on this data.
  - Blocked (2026-10-02): needs a route decision (spike-in vignette vs `VAL-01`/`ANDET-08`). The probit fit needs VAL data.

- [x] `PROG-01` — catalogue generator `extras/build_gene_programs.py`. Hand-curated
  biology joined programmatically to bundled-panel attributes, with
  **overlap groups computed by exonic interval intersection** rather than
  hand-assigned. Two bugs were caught by generating rather than assuming: a
  monotonic sweep chained all 96 EBV genes into one group, and bounding boxes
  put LMP-2A (whose exons sit at both genome ends because LMP-2 is spliced
  across the origin and the genome carries terminal repeats) in a group with
  everything. EBV now resolves to 14 groups, largest 6.
- [x] `PROG-02` — ship `src/viralscan/data/gene_programs.tsv` (79 rows, 9
  viruses) and register it in `pyproject.toml` package-data,
  `config/public_ship_scope.json` (wheel + sdist) and `MANIFEST.in`.
  `include-package-data = false` means MANIFEST alone would not ship it.
- [x] `PROG-03` — `src/viralscan/gene_programs.py`. Evidence is the
  `counts_unique_viral` layer; breadth counts distinct non-overlapping overlap
  groups; the multimap-allocated breadth is reported alongside, never merged.
  `latent` and `mixed` are unreachable when `latency_observable_in_rna=false`.
- [x] `PROG-04` — integration: optional `gene_programs` Snakemake rule gated on
  `config["gene_programs"]`, depending on `log/detection.done` **and**
  `results/viral_summary.tsv` so it cannot run before layer 1; sentinel added to
  `rule all` only when enabled; `--gene-programs` / `--programme-min-breadth`
  CLI; `viralscan rerun-programs` operating in place (layer 2 changes no counts,
  so unlike `rerun-multimap` there is no reason to copy the run).
- [x] `PROG-05` — outputs `results/gene_program_summary.tsv` and
  `results/gene_program_cells.tsv`, plus a report section that surfaces
  `panel_completeness` and `latency_observable_in_rna` so a partial row is not
  over-read.
- [x] `PROG-06` — the EBV LCL regression test. Asserts the unique layer calls
  `latent` on a matrix built to reproduce the cross-mapping scenario, so the
  protection fails loudly if either the evidence layer or the overlap-group
  logic is removed.
- [x] `PROG-07` — **closed 2026-10-03 on the called-cell set (after `PROG-17`).**
  Same copy, `rerun-programs` with PROG-17: EBV scores **932** cells (exactly
  the called cells with marker evidence; cross-checked against
  `per_cell_viral.tsv`).
  - Unique layer: **526 latent**, 67 productive, 184 mixed, 155
    indeterminate.
  - Allocated layer: 339 latent, 12 productive.
  - Inversions (latent on unique → productive on allocated): **0**. Of the
    526 latent cells, 319 stay latent and 207 turn `mixed`.
  - The block below is the superseded first pass over 1,679 barcodes.
  Superseded first pass, 2026-10-03, on layer 2's own cell set:
  `rerun-programs` ran on a **copy** of the cat42b EBV run
  (`viral_ref_cat42b/runs_prog07/SRR12682296`, post-SW-13/SW-24, after
  PROG-11, code d3905e8). Layer 2 scored 1,679 cells. Only 932 of them are among
  the 2,763 emptyDrops-called cells, so these are **not** called-cell numbers:
  - Unique layer: **695 latent**, 78 productive, 185 mixed, 721
    indeterminate.
  - Allocated layer: 526 latent, 73 productive.
  - Inversions (latent on unique → productive on allocated): **0**. Of the
    695 latent cells, 478 stay latent on the allocated layer and 217 turn
    `mixed`.
  - Stays `[~]` until `PROG-17` lands. Then rerun on the called set; the
    command is cheap (minutes).
  - Every earlier number below came from `fresh12b`: another index, no
    barcode correction (SW-13) and no cell calling. They are not comparable
    and are kept as history only.
  - `docs/faq.md` and `docs/output_reference.md` carry a dated correction.
  - **Former status:** "under re-verification (2026-09-27); do not cite
    these numbers."
  An audit of the catalogue found lytic genes filed as latent markers for EBV
  (`BaRF1.1`, the ribonucleotide-reductase subunit, matched onto latent `BARF1` by a
  case-insensitive lookup; `BHRF1`, `BNLF2a/b`), so the latent counts below may be
  inflated by lytic reads. Re-measure with `viralscan rerun-programs` once `PROG-11`
  lands. Original entry, retained as history:
  - **Re-measured 2026-09-27 with `PROG-11`'s cited changes** (`rerun-programs`
    on `benchmark_runs/reference_strategy_2026-06-28_fresh12b/runs/ebv__viralscan__combined/SRR12682296`,
    code at `50253a6`). Unique layer: **895 latent** (was 2,240), 236 productive
    (was 102), 311 mixed (was 445), 3,094 indeterminate (was 2,968); cells with
    any marker evidence 4,536 (was 5,755). Allocated layer: **856 latent** (was
    1,277), 42 productive (was 5). So the unique layer's apparent latent gain
    shrinks from +75 % to +4.6 %: most of it was `BARF1.2`, which carried 13,668
    unique and 0 allocated molecules and is not a B-cell latency marker. About
    20 % of evidence-bearing cells reach a latent call; latent sensitivity, not
    direction, is now the open question, which is what `PROG-15`'s 3′/UTR
    measurement should explain. Still provisional: `BHRF1` and `BNLF2a/b` stay
    in the latent set until verified (`BNLF2a/b` carry 0 unique molecules, so
    they cannot inflate unique-layer calls; `BHRF1` could). The pre-change
    outputs are kept beside the new ones as `*.pre-PROG-11.tsv`.
  **measured on the real run, and the reason the design exists.**
  EBV LCL `SRR12682296`: aggregate LATENT 236,342 vs LYTIC 247,633 (ratio 1.15)
  in a cell line latently infected by construction, with `EBNA-1.1` at 920 UMI
  ~155x below `BHLF1` at 142,954 — so per-gene aggregate totals are uninformative.
  Per-marker, the uniquely-placing layer has **0** molecules on `BZLF1` and
  13,668 on `BARF1.2`, where the allocated layer has 9,308 and 0. Computed
  overlap groups: g5 = {EBNA-1, EBNA-2, EBNA-LP}, g50 = {BNLF2a, BNLF2b, LMP-1},
  g39 = {BTRF1, BcLF1}. End-to-end `rerun-programs` on that run: **2,240 cells
  latent, 102 productive, 445 mixed** on the unique layer versus 1,277 latent and
  5 productive on the allocated layer, with **0 cells inverted** between the two.
  The failure mode the unique layer fixes is lost sensitivity, not a wrong
  direction — recorded here so the claim is not overstated later.
- [x] `PROG-08` — **closed 2026-10-03** (user grill: "upgrade KSHV, close
  the rest").
  - **KSHV is now `complete`.** Its latent set:
    - ORF73/ORF72/ORF71 (LANA/v-cyclin/vFLIP): one mRNA family (PMID
      9733875), counted as one breadth unit via a new `CO_TRANSCRIBED`
      override in `extras/build_gene_programs.py`;
    - K12 kaposin, also lytic-induced (PMID 17913828);
    - vIRF-3/LANA2 (K10.5), B-cell latent (PMID 11119611).
  - K1 was removed: it is tied to lytic replication (PMID 27307571).
  - Only KSHV rows of `gene_programs.tsv` changed. Tests:
    `TestKshvLatency`, plus a calling test (cluster gives breadth 1, + K12
    gives 2).
  - The other 7 are documented as biology-limited in `docs/faq.md`. No real
    KSHV dataset exists in the repo yet; that validation is `PROG-18`.
  - Original text: **open.** Four of nine viruses have a genuine latency and
  reactivation split (EBV, CMV, HHV-6A, HHV-7). The other five are `partial`:
  HSV-1/2 latency is a single transcript (`LAT`), VZV's is inferred (ORF4), and
  HHV-6B's GTF carries no attributes at all. For those, `latent` is unreachable
  by construction — so the HSV-1 benchmark cannot demonstrate a latent call.
  Extending them needs either a fuller annotation source (the HHV-6B panel has
  no `product` text at all) or a decision to accept a weaker anchor set.
- [~] `PROG-09` — **open.** The catalogue is 9 viruses. Polyomaviruses (JC/BK/KI/WU,
  MCPyV), HPV, HBV, HDV, GBV-C, HIV-1 and HTLV-1 are deliberately excluded: they
  have no latency/lytic dichotomy representable from the panel's protein-coding
  genes, and for HIV/HPV/HBV the interesting state is DNA-level latency, which
  scRNA-seq cannot observe at all. Adding them would mean emitting
  `not_applicable` rows, which is what the code does for a detected virus with no
  model — a decision to make explicitly rather than by omission.

## WP4E — Named HPV ORFs for the oncogene-versus-capsid contrast (new 2026-09-26)

Objective: give HPV real, named genes so that transcriptional activity of the
E6/E7 oncoproteins can be told apart from passive L1/L2 capsid transcription in
oropharyngeal and tonsillar tissue. This is a **different axis from WP4D** and
does not change `PROG-09`: HPV still has no latency/lytic dichotomy representable
from protein-coding genes, so it remains out of the gene-programme catalogue and
`PROG-09`'s reasoning about DNA-level latency stands. What is new here is that
the ORFs are *named at all*, which WP4D did not require.

Trigger: the built index at
`covid_viralscan/viralscan_ref/` represents HPV with 4 accessions
(`NC_001526.4`, `NC_001356.1`, `NC_001352.1`, `NC_003461.1`) and 22 t2g rows whose
gene IDs are RefSeq `locus_tag` values — `HpV16gp1`…`HpV16gp8`, `HpV1agp1`…,
`HpV2agp1`…, `Hpv1gp01`…. Nothing in the index says which is E6 and which is L1.

- [x] `HPV-01` — **investigation first: the records already carry semantic
  names.** Every HPV complete-genome record examined annotates each CDS with
  `/gene="E6"` (RefSeq and most INSDC submissions) or
  `/product="transforming protein E6"` (the records that omit `/gene`).
  **No coordinate table is used or needed**, and the per-row `name_source` column
  records which qualifier each name came from. The bundled RefSeq GTFs already
  carried the answer in a `gene` attribute; the packager kept `locus_tag` as the
  ID and dropped `gene`, which is why the index cannot answer the question.
- [x] `HPV-02` — **a coordinate table was rejected on evidence, not taste.**
  Papillomavirus genomes are submitted as linearised circles cut at the
  submitter's chosen point, so the same E6 ORF sits at 7125-7601 in
  `NC_001526.4` (HPV16) and 105-581 in `NC_001357.1` (HPV18) — opposite ends of
  their records. `NC_001526.4` is cut inside E1 and therefore reports its ORFs as
  `E1, E2, E5, L2, L1, E6, E7` while `NC_001357.1` reports `E6, E7, E1, …`; both
  are the canonical order, differing only by where the circle was opened. The
  build's ORF-order check is rotation-tolerant for exactly this reason, and the
  rotation is asserted by a test.
- [x] `HPV-03` — **names independently confirmed against protein sequence.** The
  build aborts unless every E7 translation carries the LXCXE retinoblastoma-
  binding motif, every E6 translation carries its C-X2-C zinc fingers, E6/E7/L1/L2
  fall in their known length ranges (E6 ~150 aa, E7 ~100 aa, L1 504-569 aa,
  L2 474-525 aa), and L1 is at least 3x either oncogene. HPV16 and HPV18 E7 both
  read `…LXCYEQL…`; both E6 read `…IICVYCKQQL…`. Two apparent invariants were
  falsified by writing them first and watching all 16 records fail: **L1 is not
  the longest ORF** (E1 is, at ~650 aa — it is the replication helicase), and
  **L1 is not always longer than L2** (HPV-2 annotates L2 at 525 aa against L1 at
  511 aa; HPV-1 has them within one residue). Both are recorded in the code.
- [x] `HPV-04` — genotype coverage: **all 14 high-risk types** (16, 18, 31, 33,
  35, 39, 45, 51, 52, 56, 58, 59, 66, 68) plus HPV1 and HPV2, which are carried
  only so a rebuild does not *lose* the two types the current index already has.
  RefSeq has complete genomes for only 4 of the 14 (`NC_001526.4` REVIEWED,
  `NC_001357.1` VALIDATED, `NC_075191.1` and `NC_075233.1` PROVISIONAL); the
  other 10 are INSDC, chosen as the **oldest** complete genome carrying annotated
  CDS, since the earliest submission of a type is the prototype that genotyping
  assays and published amplicons target. Attempted and documented: 208 RefSeq
  papillomavirus complete genomes were enumerated (65 human, 57 types), and every
  one of the 10 remaining high-risk types was confirmed present in INSDC (37-572
  isolate records each) before one was selected.
- [x] `HPV-05` — `extras/build_hpv_reference.py` regenerates the TSV from NCBI,
  cache-first through `ncbi_fetch`'s own cache directory with its SHA-256 sidecar
  convention, so a record already fetched by any other ViralScan entry point is
  reused rather than re-downloaded. Refuses to write a partial catalogue: it exits
  non-zero on any unresolvable accession, unrecognised ORF symbol, failed
  invariant, or duplicate `(accession, gene)` pair.
- [x] `HPV-06` — `src/viralscan/data/hpv_genes.tsv` (124 ORFs, 16 genotypes) and
  `src/viralscan/hpv_genes.py`. Registered in `pyproject.toml` package-data and
  `MANIFEST.in`; `include-package-data = false` means MANIFEST alone would not
  ship it. Gene IDs are namespaced (`NC_001526.4_E6`), because the bare symbol
  `E6` occurs once in each of 16 genomes. Classes: 32 `oncogene` (E6+E7 per
  genotype), 1 `oncogene_locus` (HPV16 E6*), 32 `late_capsid`, 59 `early`.
- [x] `HPV-07` — **`E5` is classed `early`, not `oncogene`, and `E6*`/`E7*` are
  `oncogene_locus`, not `oncogene`.** E5 is a transforming protein several
  reviews call an oncogene, but it is not part of the E6/E7 axis and putting it
  in the oncogene class would make a positive call mean something it does not.
  E6* lacks the PDZ-binding motif and E7* lacks LXCXE, so neither is
  transforming, and their reads are indistinguishable from E6/E7 at the sequence
  level — folding them in would let E6* alone produce a confident positive
  oncoprotein call. Resolved by preferring an isoform symbol in `/product` over
  its own parent in `/gene`, which is what `NC_001526.4` requires.
- [x] `HPV-08` — `--emit-reference DIR` writes a merged FASTA (16 genomes,
  128 kB) + GTF (125 exon lines, 124 unique gene IDs) for `kb ref`, built from the
  parsed CDS rather than from `ncbi_fetch._genbank_to_gtf`. Sequences are
  generated on demand and **not committed**: reference size is a packaging
  decision governed by PR 8 and should be made once, deliberately. Structurally
  validated (FASTA/GTF seqnames agree, blocks within bounds); **not yet
  index-built**, as `kb`/`kallisto`/`bustools` are not on PATH in this
  environment.
- [x] `HPV-09` — **blocker, in `ncbi_fetch.py`.** `_genbank_to_gtf` set
  `gene_id` from `/gene=`, so in a merged multi-genome reference every
  genome's E6 collapses onto one row. Demonstrated: merging `NC_001526.4` and
  `NC_001357.1` through `fetch_reference()` yielded 16 ORFs on **9** gene IDs.
  A public `fetch_genbank()` accessor was also needed so this generator
  stopped importing `_efetch`, `_cache_valid` and `_write_cached`.
  **Both asks were implemented in WP4F (2026-09-26)** and **verified here
  (2026-09-26):** `_genbank_to_gtf` now emits genome-scoped `<accession>_<token>`
  gene IDs, and `fetch_genbank()` is public and caches the raw flatfile.
  Verification: the two cached `NC_001526.4`/`NC_001357.1` `.gtf` files did
  predate the fix (see `ANELLO-11`) — deleted and confirmed they regenerate
  from the retained `.gb` flatfile with no network call. Re-running the exact
  `fetch_reference(["NC_001526.4", "NC_001357.1"], ...)` merge that
  demonstrated the bug now yields **17 distinct gene IDs across 19 ORF lines**
  (the one repeat is a real intra-genome duplicate `locus_tag` on HPV16,
  correctly disambiguated `_dup2` by `_panel_gene_ids`, not a cross-genome
  collision). `extras/build_hpv_reference.py`'s own `fetch_genbank()` wrapper
  now delegates to the public `ncbi_fetch.fetch_genbank()` instead of
  reaching into its privates; `tests/test_hpv_reference.py` and
  `tests/test_ncbi_fetch.py` pass unchanged.
- [x] `HPV-10` — **the scientific limit, recorded so it is not over-read later.**
  What the catalogue supports: an **oncogene-versus-capsid contrast per
  genotype**. E6/E7 are early-region oncoproteins transcribed in
  carcinogen-driven HPV-positive oropharyngeal tumours, while L1/L2 are
  late-region and transcribed only in productive infection, so E6/E7 reports viral
  gene expression where L1 reports virion production. What it does **not**
  support: confident **per-genotype attribution of L1 signal**. L1 is the most
  conserved coding region in the genus — it is what pan-HPV PCR primers target —
  so L1 reads cross-map freely between all 16 genotypes and kallisto's
  multimapping will distribute them; a genotype label on an L1 count is not
  independent evidence of which type is present. E6/E7 are far more
  type-divergent and better behaved, but a cross-mapped count is still not a
  transcript count. Presence-of-HPV-transcripts and the oncogene-versus-capsid
  distinction are supportable; per-type L1 attribution is not. This is also a
  *transcriptomic* catalogue: a transcriptionally silent integrated genome — the
  common state in tonsillar crypt epithelium and the state that drives
  HPV-positive oropharyngeal carcinoma — produces no reads and is invisible
  here, not negative.
  - Closed 2026-10-02 (`v3/docs`): the limit is now in user docs
    (`reference_panel.md` and the FAQ). Presence and E6/E7 vs L1/L2 are
    supportable; per-genotype L1 attribution and low off-type calls are not.
    The measurement stays in `HPV-11` / `CAT-08`.
- [ ] `HPV-11` — next: rebuild the real index with the named HPV panel and
  measure the cross-mapping directly, rather than reasoning about it. Requires
  `HPV-09` or the `--emit-reference` path, plus `kb` on PATH. Until then the
  L1 claim in `HPV-10` is a literature-based expectation, not a measurement from
  this panel. Do not publish a per-type HPV number before this row is closed.
  - **2026-09-26, not evidence:** a local run tiled 1,837 error-free 90 bp reads
    from the panel's own E6/E7/L1/L2 CDS, pseudoaligned them against a
    panel-only `kb ref` index, and found **0** cross-genotype assignments (19
    within-genome only, E6/E8 and L1/L2 junctions). That result is expected by
    construction: reads drawn from indexed sequences match themselves, and HPV
    types are defined by ≥10 % L1 divergence. It says nothing about HPV-10's
    actual risk, which is reads from strains *not* in the panel. The valid test
    is `CAT-08`'s leave-one-out design: held-out isolates (HPV16 lineage variants,
    types outside the panel) scored for where their reads land.

## WP4F — Real gene structure for the Anelloviridae panel (new 2026-09-26)

Objective: replace the Anelloviridae panel's one-placeholder-gene-per-genome
representation with the CDS features NCBI actually annotates. This is a
**confirmed bug fix, not a research question**, and it is a different axis from
WP4D/WP4E: anelloviruses have no latency/lytic dichotomy and no oncogene/capsid
contrast, so they stay out of both catalogues. What they do need is for a hit to
be attributable to a *locus* rather than to a conservation rank.

### Root cause (confirmed in code, not inferred)

The panel never took the GenBank path, and the reason is a **discarded** GTF, not
a failed fetch:

| Site | What it does |
|---|---|
| `src/viralscan/scripts/ncbi_fetch.py` `_fetch_one` | fetches GenBank, runs `_genbank_to_gtf`, writes a real GTF, **discards the GenBank text** |
| `scripts/build_bundled_panel_ref.py:228` (was) | `anello_gtf_texts.append(_genome_as_transcript_gtf(text, acc))` — threw the real GTF away, rebuilt a placeholder from the FASTA |
| `src/viralscan/scripts/build_reference.py:1072` (was) | `_gtf_from_merged_fasta(final_fasta, final_gtf)` — same discard, in `build_anellovirus_reference` |
| `src/viralscan/scripts/build_reference.py:437` | `_genome_as_transcript_gtf` emits `gene_id = f"{accession}_gene{seq_idx}"` — the placeholder itself |

So the briefed hypothesis (that anelloviruses took
`_whole_genome_gtf_from_fasta` in `ncbi_fetch.py`) is **wrong for 97.7 % of the
panel**. Measured on the pre-existing NCBI cache: 1,995 of 2,042 cached `.gtf`
files already carried real CDS-derived exons and only 47 were placeholders. The
dominant emitter is `_genome_as_transcript_gtf` in `build_reference.py`, reached
from two callers. `ncbi_fetch._whole_genome_gtf_from_fasta` is real and does fire
for genuinely CDS-less records, but it is the minority path.

There was a **second, independent** defect that would have survived the first
fix: `_genbank_to_gtf` set `gene_id` from the bare `/gene=` or `/product=`, which
are not unique. Across the cached panel only 2,316 distinct gene IDs existed
across 1,995 annotated genomes, with `ORF1` shared by 150 genomes and `orf1` by
63. A merged index would have collapsed 1,995 genomes onto 2,316 columns with
cross-genome identity. `HPV-09` independently reports the same defect for HPV.

### Measured CDS coverage (full panel, 2,042 accessions)

Every accession was retrieved and audited; the numbers below are the generator's
own output, not a projection. **1,995 of 2,042 (97.7 %) carry real gene
structure — 2,515 genes — and 47 have no CDS feature in NCBI at all.**

| genus | panel | annotated | genes | coverage |
|---|---:|---:|---:|---:|
| Betatorquevirus | 1,542 | 1,517 | 1,699 | 98.4 % |
| Alphatorquevirus | 211 | 204 | 361 | 96.7 % |
| Anelloviridae (unclassified) | 185 | 175 | 264 | 94.6 % |
| Gammatorquevirus | 78 | 74 | 150 | 94.9 % |
| Hetorquevirus | 8 | 8 | 9 | 100 % |
| Samektorquevirus | 8 | 7 | 10 | 87.5 % |
| Gyrovirus | 6 | 6 | 17 | 100 % |
| Memtorquevirus | 4 | 4 | 5 | 100 % |
| **TOTAL** | **2,042** | **1,995** | **2,515** | **97.7 %** |

A 206-accession stratified pre-implementation survey across all eight genera gave
205/206 (99.5 %), consistent with the full run. One record was identified as
CDS-less in the survey (`KP343852.1`, a bare `source` feature); 47 are CDS-less in
full.

**Coverage is not the problem. The low CDS *count* is — and it is biology, not a
gap.** 1,740 of the 1,995 annotated genomes (87 %) carry exactly **one** CDS, and
1,713 of the 2,515 genes have the product `ORF1`. 75 % of the panel is
Betatorquevirus, whose ~2.8–3.0 kb TT-mini genomes genuinely have a single ORF
spanning the genome. There is no missing annotation to recover. The gain from
this fix is that the single gene is *named, product-labelled, and
genome-scoped* instead of anonymous, and that the 255+ multi-ORF genomes finally
get their real structure. **A further honest caveat on annotation quality:** 576
of 2,515 genes (23 %) carry the generic product `hypothetical protein`, and only
258 genes carry a `/gene` symbol at all — so "named gene" means
product-labelled for most of the panel, not functionally annotated.

- [x] `ANELLO-01` — the audit above, run **before** any code was written, so the
  decision to keep a placeholder fallback is evidence-based rather than assumed.
- [x] `ANELLO-02` — `ncbi_fetch.py` now caches the **raw GenBank flatfile** as
  `<accession>.gb` with the existing `.sha256` sidecar convention, and exposes it
  through a public cache-first `fetch_genbank()`. Re-deriving the annotation
  after a code change therefore costs zero NCBI requests. This also satisfies
  the public accessor `HPV-09` asked for; that row is left for its owner to
  close.
- [x] `ANELLO-03` — **gene naming: `<accession>_<token>`, token from
  `/locus_tag` → `/gene` → `/protein_id` → `cds<N>`.** `/locus_tag` is the
  submitter's stable locus name and wins when present — it is what makes the
  reference TTV record emit `NC_002076.2_TTVgp1/2/3`, i.e. the same
  `TTV_TTVgp1` names the bundled RefSeq GTF already uses. `/protein_id` is the
  fallback because it was observed to be reused by **zero** of the 206 sampled
  accessions, making it the only globally unique identifier NCBI offers for the
  238/275 sampled CDS features that carry no `/gene` and no `/locus_tag` at all.
  `orf2/5` is sanitised to `orf2_5`. The bare symbol is preserved in `gene_name`
  (GTF) and `gene_symbol` (TSV). **Missing and duplicated products:** a CDS with
  no identifier gets `cds<N>`; a token repeated within one genome gets `_dup2`,
  `_dup3` (not a bare ordinal, which is indistinguishable from a real `orf12`).
- [x] `ANELLO-04` — **the naming is genome-scoped on purpose.** The existing
  bundled convention is `{virusToken}_{locusTag}` (`TTV_TTVgp1`); for a
  2,042-genome panel the virus token must be the accession, because `TTV` would
  collapse the panel into one label. `{accession}_{token}` keeps the shape and
  works unchanged with `anello_name_map()`'s boundary-aware prefix rule.- [x] `ANELLO-05` — **circular topology, measured rather than assumed — and the
  panel does contain a wrap.** Anelloviridae are circular ssDNA, but NCBI
  annotates in a *linear* representation and only 906 of 2,515 genes (36 %) come
  from records that even declared `circular` on the LOCUS line, so the
  declaration is recorded (`topology`) and never used to interpret coordinates.
  Across the full panel there is **exactly one origin-spanning gene**:
  `KU243129.1` (2,824 bp, `ss-DNA`, `circular`) annotates
  `join(2677..2824,1..80)` — exon 1 at the end of the linear representation, exon
  2 back at the origin. NCBI writes the intervals in transcript order, so a
  coordinate sort would emit `1..80` first and silently transpose the gene's two
  exons into a scrambled transcript. The parser therefore **never re-sorts**; the
  only normalisation is reversing the interval list for minus-strand features,
  and `_origin_spans` flags a wrap structurally (in transcript order the first
  interval starts after the last interval ends) with `origin_spanning="true"`.
  A test asserts the real `KU243129.1` gene's exon *order*, not just its
  existence, so a regression here is a visible transposition rather than a
  silent one. No interval anywhere in the panel falls outside `[1, length]`.
- [x] `ANELLO-06` — **genogroup is not derivable and was not invented; the
  column ships empty except two verbatim NCBI values.** Corrected 2026-09-26,
  and re-verified directly against the cache rather than re-asserted: an
  earlier draft of this row claimed **zero** of 2,058 flatfiles carry a
  `/genotype` qualifier and cited two `/note`-only hits instead
  (`NC_002076.2`, `JN980171.1`). That "zero" claim was simply wrong — grepping
  the retained flatfiles for a literal `/genotype=` qualifier line
  (`grep -l '/genotype=' ~/.cache/viralscan/ncbi/*/*.gb`) finds exactly **two**
  records, `NC_014081.1` (`"6"`, `/organism` "Torque teno virus 3") and
  `NC_014094.1` (`"28"`, `/organism` "Torque teno virus 6") — both of which
  contradict their own organism species number, which is the concrete
  evidence that a genogroup must never be inferred from `/organism`. The
  `NC_002076.2`/`JN980171.1` pair is real too, and distinct: those two carry
  `genotype` only as free text inside a `/note` (confirmed separately with
  `grep -l '/note=.*genotype'`), never as a structured qualifier, so they stay
  out of `source_genotype` for the same reason the column is not back-filled
  from prose generally. Four different records, two different mechanisms —
  not a contradiction to reconcile, just two separate, now-verified facts.
  The shipped TSV column is named `source_genotype`, not `genogroup`,
  precisely because it is NCBI's own `/genotype` qualifier copied verbatim
  rather than a derived or inferred genogroup, and it is empty for all but
  the two `/genotype=` accessions. `test_source_genotype_is_never_invented`
  pins the exact pair so a future hand-fill or a broader NCBI regression
  cannot pass silently. Retained source fields, recomputed from the shipped
  TSV rather than copied from the generator's own (stale) docstring:
  `source/isolate` (2,338 genes across 1,899 accessions, laboratory sample
  codes such as `MDJHem2` or `SAfiA-468-6`) and `source/strain` (102 genes
  across 31 accessions), carried verbatim so a classifier can be fitted later
  without re-fetching.
- [x] `ANELLO-07` — `extras/build_anellovirus_genes.py` → `anellovirus_genes.tsv`
  → `viralscan.anellovirus.gtf_text_for()`, mirroring the `build_gene_programs.py`
  → `gene_programs.tsv` → `gene_programs.py` precedent. Cache-first and
  resumable (flatfiles cached, nothing written until every accession is
  attempted), `--accessions` / `--limit` / `--per-genus` for subset runs,
  `--min-coverage 0.95` so a silent NCBI regression cannot ship a
  mostly-placeholder catalogue, and a coordinate/uniqueness audit that refuses to
  write on any violation.
- [x] `ANELLO-08` — all three discard sites now consume the catalogue:
  `build_bundled_panel_ref.py` Step 4b, `build_anellovirus_reference` Step 4, and
  `build_combined_reference` Step 3 (**scoped to panel accessions only**, so the
  curated 195-genome bundled panel keeps its byte-identical whole-genome GTF and
  no existing index changes shape). Uncovered accessions still get a placeholder,
  because `kb ref` silently drops a sequence with no GTF row and the genome would
  then be neither quantified nor detectable.
- [x] `ANELLO-09` — `tests/test_anellovirus_reference.py`: 42 offline tests
  (panel integrity, no surviving `_gene1`, genome-scoped and unique gene IDs,
  spliced genes have >1 exon, coordinates within `[1, genome_length]`, strand
  `±`, ORF1/Rep present in all eight genera, genogroup stays empty, committed
  GenBank fixtures for the converter including the circular-wrap and
  minus-strand cases, generator run offline from a seeded cache) plus three
  `@pytest.mark.network` tests. The default selection needs no network and runs
  in ~4 s.
- [x] `ANELLO-10` — registered in `pyproject.toml` package-data, `MANIFEST.in`
  and `config/public_ship_scope.json` wheel **and** sdist allowlists, mirroring
  `gene_programs.tsv` / `hpv_genes.tsv`. Not added to `.dockerignore` /
  `docker_context`, matching how those two shipped.
- [x] `ANELLO-11` — **behaviour change to a private helper, recorded because it
  is observable.** `_genbank_to_gtf` gene IDs are now genome-scoped, so a
  cached `<acc>.gtf` written by an older build is stale. Delete the affected
  `~/.cache/viralscan/ncbi/<acc>/<acc>.gtf` (and its `.sha256`) to re-derive; the
  retained `.gb` means that costs no network. One assertion in
  `tests/test_ncbi_fetch.py` was updated to the new contract and the
  genome-scoping property given its own test.
- [~] `ANELLO-12` — **the honest limit, recorded so it is not over-read.** The
  bug is fixed and the 99.8 %-in-one-bucket artifact is now *explainable and
  testable*: with one anonymous gene per genome there was nothing else it could
  have been, and the decisive test is to rebuild the index and check that
  anellovirus UMI spread across many genome-scoped gene IDs instead of
  concentrating in the most conserved one. **That test is not run here** — it
  needs `kb ref` and a re-run of the COVID sample, neither available in this
  environment (`ANELLO-13`). **What remains true regardless, and is the
  important caveat: 1,740 of 1,995 annotated genomes (87 %) carry exactly one
  CDS spanning the whole genome.** For those genomes a "real gene" is still a
  whole-genome transcript, so the count is still a whole-genome count and still
  cross-maps against the rest of the panel in proportion to conservation. The fix
  makes that cross-mapping *attributable and visible* — a count can now be
  traced to a genome and a product, and a conservation-driven skew is
  distinguishable from a single-genome infection — but it does **not** make
  per-genotype anellovirus quantification defensible, and it cannot: a
  unique-sequence argument is needed, not a gene name. That is a different piece
  of work, closer to `SENS-06` (measured k-mer capture) than to this fix. The
  honest summary is: **genus-level anellovirus load becomes interpretable;
  genotype-level anellovirus load does not.**
  - Blocked (2026-10-02): the decisive test is `ANELLO-13` (`[ ]`, not started). Kept `[~]` as a recorded limit.
- [ ] `ANELLO-13` — next: rebuild the reference with the real-gene panel and
  **measure** the covid artifact rather than reasoning about it. Requires `kb` on
  PATH and a re-run of the COVID scRNA-seq sample; then compare the anellovirus
  UMI distribution across genome-scoped gene IDs against the 99.8 % single-bucket
  baseline. Until this row closes, the covid number in this section is a
  measurement of the *old* build and says nothing about the new one. Do not
  report a per-genotype anellovirus number before this row is closed.
  - **2026-09-27:** the bespoke covid build cannot run this test. It reads the
    static `references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf`,
    which still carries 6,126 `_gene` placeholder IDs, so a rebuild through
    `slurm_build_ref_v2.sh` would reproduce the old artifact. Decided path: native
    `viralscan build-ref`, which first needs `CAT-01` (build-ref currently discards
    the real GTF). The retraction question is already answered by F-005
    (`.living/findings/`); this row only validates the WP4F code fix.
- [x] `ANELLO-14` — the integration test
  - **Done 2026-10-04 (672ac08).** Seeded random fixture that passes the CAT-17
    gate. Expects `NC_002076.2_TTVgp1-3` and `AB026929.1_BAA86944-6.1`.
    4 passed.
  `tests/integration/test_anellovirus_chain.py::TestAnellovirusLabelingChain::test_build_anellovirus_reference_produces_labelable_gtf`
  fails, and the failure predates WP1C (it fails at `af4d5b3`). The test still
  asserts the pre-WP4F placeholder `NC_002076.2_gene1`, but the builder now
  emits real gene IDs (`NC_002076.2_TTVgp1`…). Update the test's expected IDs
  to the WP4F gene structure. Do not revert the builder.

## WP4G — Anellovirus detection you can trust (new 2026-09-27)

Objective: make an F-005-type host-homology artifact impossible to publish by
default, and make anellovirus sensitivity a measured quantity. An audit
(2026-09-27) found no default step that would block the covid Alphatorquevirus
call today. "Reliable" here means a genus-level LOD95 from planted reads in 10x
geometry plus host-only negatives that never reach a reported call; without an
orthogonally confirmed positive sample (`REF-10`) every anellovirus result stays
`screening_only`, and no code change lifts that ceiling.

- [~] `ANELLO-15` — (2026-10-07, user request) use `clareaulab/human_anellovirus_pangenome` correctly.
  - **Provenance corrected 2026-10-07.** The repo is the renamed `anellovirus_reference` (GitHub 301,
    repo id 1072468567), HEAD `3ed77e1`, FASTA blob `173f981` unchanged since 2025-10-09. All 2,020
    `clareaulab` rows are version-identical members of its 2,023 representatives; the other three are
    `AB303562.1` (we dropped it, CAT-05), `NC_038337.1` and `NC_038359.1` (our RefSeq rows). It has **no
    licence, paper or DOI**, so docs now cite URL + commit, claim "accession list only" and drop the
    unverifiable "Lareau 2023". Its clustering is vclust (95 % ANI, 85 % qcov), not `cd-hit-est`.
    We keep fetching sequence and CDS from NCBI: its FASTA is pre-masked (dustmasker level 20, no way
    back to the raw bases, so it would bypass `MASK-01`), its GTF is one whole-genome feature per genome
    (its `anello_t2g.txt` uses `X_transcript` IDs that do not match its FASTA/GTF), and it has no host
    screen.
  - Open: (1) ask the lab for a LICENSE, a citable release and permission to ship `orf1_genus`
    (resolves 171 of our 185 generic "Anelloviridae" rows; it also disagrees on 9 Beta and 3 Gamma rows,
    so reconcile before use); (2) optional integrity check of NCBI-fetched sequence against its FASTA
    on non-N bases; (3) its `cdhit_clusters.tsv` membership as the redundancy map for `CAT-07`.
  - **Count drift fixed 2026-10-07.** Current counts: 2,041 table accessions (2,020 `clareaulab` + 21
    RefSeq); 2,021 of them have no bundled GTF (was "2,022 of 2,042"; the "91 % of the panel" beside it
    was wrong arithmetic, it is 99 %); `anellovirus_genes.tsv` has 2,515 data rows over 1,995
    accessions, of which `AB303562.1` (4 rows, the dropped duplicate) is not in the table, so 2,511 rows
    over 1,994 accessions are live. Docs, help text and docstrings now quote these; dated measurements
    on the 2,042-genome panel (k-mer capture, CDS census) are left as measured. A test pins the numbers.
    The 4 stale `AB303562.1` gene rows were removed in `f4320d2` (2026-10-10; `tests/test_anellovirus_reference.py` pins 2,511 rows / 1,994 accessions). `reference_strategy` template `anellovirus_expected_count` 2022 -> 2021.
- [~] `ANDET-01` — `accession_breadth` is always 1.0: it is computed over
  - 2026-10-04: gene-level `accession_breadth` now uses the virus's **index** genes
    (`compute_stats(index_genes_by_virus=...)`). The per-accession
    genome-coverage breadth (F-005's ≤3.41 % gate) needs read positions; it
    belongs to the evidence route and stays open.
  - 2026-10-10 (wave 4): the alignment branch's per-accession `breadth` /
    `breadth_unique` (already in `anello_alignment_by_accession.tsv`) are now
    aggregated, read-weighted, into virus-level `alignment_breadth` and
    `alignment_breadth_unique` in `anello_alignment_by_virus.tsv` and
    `viral_summary.tsv` (via `SUMMARY_COLUMNS`). A label only, never a gate:
    5′ capture concentrates reads (see the TTV decision rules). The branch is off
    by default. Still open: the default kallisto route has no read positions, so
    no breadth for non-alignment rows; `output_reference.md` rows pending
    (hash-pinned, owner adds and re-pins; `test_emitted_columns_are_documented`
    fails until then).
  - 2026-10-12 (`codex/aud13-andet-followups`): reconciled. The per-accession breadth exists wherever read
    positions do: `breadth`/`breadth_unique` in `anello_alignment_by_accession.tsv`, the per-virus
    `alignment_breadth*` (incl. `_max`, 2026-10-11, documented in `output_reference.md`), and `breadth_1x/3x/10x`
    per reference in the evidence route's `alignment_qc.tsv`. Blocked by design, not built: the default
    kallisto route has no read positions, and the F-005 gate (<= 3.41 %) is a threshold (REF-08/user).
  `found_genes`, which are already detected (`detection.py:166`, `:501-506`).
  Compute it over every index gene of the virus and add per-accession
  genome-coverage breadth, F-005's deciding gate (≤3.41 %).
- [~] `ANDET-02` — read `host_homology_annotations.tsv` (written at
  `build_reference.py:768`, read by nothing) in detection; demote calls
  concentrated in host-homologous regions; surface `eve_risk` in the report.
  - **Measurement surfaced 2026-10-07; demotion NOT done (decision: advisor + user default narrowed).**
    `viral_summary.tsv` gains `host_homology_status` (`measured`/`partial`/`not_measured`) and
    `host_homology_max_identity`/`_max_query_coverage`/`_max_aligned_bases` (maximum over the virus's
    genomes), read by `virus_grouping.host_homology_by_virus` from `reference_manifest.json` beside the
    index (the per-sequence copy of the TSV; `not_measured` is never reported as 0). No threshold, no
    change to any call or `eve_risk`: a cut chosen now would alter round-1 calls before the calibration
    gate, and the cut belongs to `REF-08`.
  - Still open: (1) demotion, once `REF-08` freezes a threshold; (2) "concentrated in" needs per-gene or
    per-interval homology, which `_parse_host_homology_paf` discards (it keeps one maximum per genome), so
    the writer must keep coordinates (candidate source: the `REF-07` per-locus table, now untracked
    `scripts/ref07_host_homology_table.py`); (3) no shipped reference has the measurement yet
    (`viral_ref_final` has no `host_homology_annotations.tsv`), so every current run reports
    `not_measured` until a `--genome-dlist` rebuild.
  - 2026-10-10 (wave 4): **measurement defect fixed.** `measure_host_homology` ran
    `minimap2 -x asm10 <host_genome> <viral_fasta>`, which `REF-07` showed returns zero alignments for all
    2,343 genomes, so every `--genome-dlist` rebuild would have recorded `max_aligned_bases = 0` ("measured",
    no homology) for the whole panel. It now indexes the viral panel and streams the host genome as the query,
    with the flags of `scripts/ref07_host_homology_table.py` (`-c -k 15 -w 10 -N 200 -p 0.01 -s 40 -m 20`).
    `_parse_host_homology_paf` reads the swapped columns: identity, coverage and aligned bases stay per
    **viral** genome (coverage is of the viral genome; `host_homology_best_target` is the host contig).
    Tests: swapped-column PAF unit test and a real-minimap2 test with a synthetic host carrying a diverged
    viral fragment (`tests/test_build_reference.py`). Output columns are unchanged, so
    `docs/output_reference.md` needs no edit. The demotion stays blocked on `REF-08` (threshold) and the
    per-locus coordinates (2); the row stays `[~]`.
  - 2026-10-12 (`codex/aud13-andet-followups`): item (2) was closed on 2026-10-11 (see `AUD-16`: loci kept in the
    reference manifest, `results/host_homology_loci.tsv` with `eve_risk`); what was missing was the output
    documentation, now `output_reference.md` section `host_homology_loci.tsv` with a test that every column is
    documented. Still open and **not implemented on purpose**: demotion (blocked on `REF-08`'s threshold).
- [x] `ANDET-03` — `claim_scope` column (`screening_only` for Anelloviridae) in
  - **Done 2026-10-04 (9361d56).** viral_summary `claim_scope`: `screening_only` for
    Anelloviridae (catalogue family, genus-name fallback), empty otherwise;
    alignment-only rows too. Documented in output_reference. The HTML report
    renders the summary table, so the column appears there.
  `viral_summary.tsv` and the report. `REF-10`'s label exists only in prose today.
- [~] `ANDET-04` — evidence replay reads the raw FASTQs (`evidence_run.py:154-155`,
  `:176`) instead of the host-filtered `kb_r1`/`kb_r2`; `--virus ttv` resolves to
  the 185 unclassified genomes only (`evidence.py:64`). Fix both; auto-run
  read-level host confirmation for detected anellovirus genera.
  Both bugs fixed 2026-09-27. `replay_fastqs()` returns the pair `kb count`
  quantified (and reconstructs `host_filtered/` for configs that predate
  `kb_r1`), both the replay and the extraction use it, and a missing input now
  fails closed instead of silently replaying different reads. `ttv` (or
  `Anelloviridae`) now selects every genus plus the unclassified group and the
  bundled "Torque teno virus" label; reproduced before the fix as
  `('Anelloviridae', ['AB303555.1_ORF1'])` on a three-genus fixture. Still open:
  the auto-run Snakemake rule.
  - Blocked (2026-10-02): the auto-run gate depends on `ANELLO-PRIOR.3` (`[ ]`). Default on vs opt-in is a user decision.
  - 2026-10-10 (wave 4): both bug fixes ticked (`[x]`: `replay_fastqs()` in
    `scripts/evidence_run.py`; family resolution in `evidence.resolve_viral_target`). Auto-run rule done,
    `[~]`: opt-in `auto_evidence` rule (`--auto-evidence`, `RunConfig.auto_evidence`/`host_fasta`/`viral_fasta`,
    `scripts/auto_evidence.py`). OFF by default and stays opt-in until `ANELLO-PRIOR.3` closes (user decision).
    - Trigger: Anelloviridae genera (identity table) in `viral_summary.tsv`, minus `alignment_only` rows
      (the replay only sees kallisto ECs). Each runs `viralscan.reads.extract_virus_reads`, which gained
      optional `viral_fasta`/`host_fasta` and calls `run_evidence`; output `results/evidence/<genus>/`,
      `log/auto_evidence.done` lists the genera or `none`.
    - Fail closed, in `errorhandler` before any work: `--auto-evidence` needs `--host-fasta`, a viral-only
      FASTA, multimapping (the resolved BUS text), minimap2 and samtools. `RunConfig` re-checks it.
    - Viral-only FASTA: `--viral-fasta`, else `<index dir>/viral.fa`, which only
      `scripts/build_bundled_panel_ref.py` writes. `build-ref` writes none, so the flag works only with
      panel-built references or an explicit file; `combined.fa`/`reference.prepared.fa` hold host sequence and are never
      used. The new args default to None so default runs keep their `--resume` fingerprint.
    - Cost: one `evidence` run per genus (host FASTA copied, FASTQs re-read); `ponytail:` comment in the script.
    - Tests: `tests/test_auto_evidence.py`, `test_cli.py::TestBuildRunConfig`, `test_snakefile_dag.py::TestAutoEvidenceDag`.
      Proposed `output_reference.md` text (hash-pinned, not edited) is in the wave 4 report.
- [x] `ANDET-05` — one genus name per genome: bundled `TTVgp1` IDs resolve to
  - **Done 2026-10-04 (cherry-picked bc63055).**
    - `VARV` → Variola virus. It was VZV, which folded Variola counts into
      VZV.
    - `UUKU` → Uukuniemi virus.
    - The TTV prefix → Alphatorquevirus, so `TTVgp1` and `NC_002076.2_TTVgp1`
      agree.
    - Side effect: legacy `TTV_*` ids also resolve to Alphatorquevirus; their
      genus is unverified.
  "Torque teno virus" while genome-scoped `NC_002076.2_TTVgp1` resolves to
  "Alphatorquevirus". Also fix the `UUKU` and `VARV` aliases.
- [x] `ANDET-06` — correct three docs: `REF-01`'s "now the default" (the CLI
  - **Done 2026-10-04 (489ede2).**
    - docs/api.md: `include_anellovirus` defaults to False.
    - The host_filter STAR_FILTER_ARGS comment is corrected against
      `STAR --help` (mismatch and multimap defaults are 10, not 0/1).
    - cli_reference now states only the measured starsolo vs genome-dlist
      ordering.
  default is off), the STAR-defaults claim at `host_filter.py:176-192`, and the
  genome-reference ranking in `docs/cli_reference.md:200`.
- [ ] `ANDET-07` — make the expanded panel the default quantification reference,
  - 2026-10-12 (grill, user-confirmed): expanded panel stays opt-in until `ANDET-08` measures the genus-level LOD95.
  **only after `ANDET-01`–`ANDET-04` land**. `protocol.yaml`'s
  `anellovirus_expanded` entry says "never the default" and sits under the frozen
  `frozen_inputs` digest, so this needs a `DEV-0xx` deviation record.
- [ ] `ANDET-08` — genus-level LOD95: plant held-out genomes from all 8 genera in
  10x v3 geometry into a checksum-pinned healthy-PBMC background (`VAL-04`),
  probit fit per the `SCI-03` method; the unplanted background must yield no
  reported call. A minimal `VAL-01` slice; closes `SENS-06` for this family.
- [ ] `ANDET-09` — STARsolo anellovirus alignment branch that does not depend on
  kallisto (grill 2026-10-03). The plan this replaces assumed bulk and rustar.
  The user chose scRNA inside v3, STAR/STARsolo, and default-on.
  - The evidence replay (`ANDET-04`) only sees reads in kallisto ECs, and
    31-mer capture is capped near 0.95^31 (F-013). Alignment is not.
  - Input is `host_filtered/R2`+`R1`. Without `--host-filter starsolo` the
    branch records `skipped_no_host_filter` (user).
  - [x] `ANDET-09a` index: `anello_star/` next to the kb index. Anelloviridae
    records of the panel's own `viral.fa`, one gene per contig, STAR
    SAindexNbases 10 / ChrBinNbits 11 (`scripts/build_bundled_panel_ref.py`).
    Built for cat42b 2026-10-03 (job 25696001): **2,040 contigs**, STAR
    2.7.11b, genome version 2.7.4a, `manifest.json` records the FASTA sha256.
  - [x] `ANDET-09b` Snakemake rule `anello_align`. STARsolo, Unstranded
    (F-020), GeneFull, `--soloMultiMappers Unique EM`, plan §11 filters,
    MultimapNmax = SAMmultNmax = 100.
  - [x] `ANDET-09c` evidence metrics merged into `viral_summary.tsv`, plus
    `detection_source`. Any genus with ≥1 unique alignment molecule gets a
    row (user). Labels, never filters.
  - [x] `ANDET-09d` `--anello-align/--no-anello-align`, **default False**.
    **Falls back to off if `ANDET-09e` fails** (user); the DEV record is
    written only on pass.
  - [x] `ANDET-09e` acceptance (`ANELLO-PRIOR.1` slice).
    - Plant: ~8 Modha BK genomes (absent from the panel), 10/100/1000
      molecules × 5′/3′/uniform windows, into 5 M x213 `host_filtered` pairs.
    - Negative: a synthetic set of host cDNA, GRCh38 and poly-G/poly-A/CAG/TSO
      artefact reads (user).
    - Pass: STARsolo planted−unplanted ≥ kallisto for every genome and >
      kallisto below 95 % identity; 0 alignment molecules in the negative; no
      planted read lost to the host filter.
    - **Result 2026-10-03** (jobs 25696086 plant, 25696105 arms, 25696114
      report). Criterion 2 **pass**: 0 alignment molecules in the synthetic
      negative, and 0 for kallisto too. Criterion 3 **pass**: 26,640/26,640
      planted reads survived the host filter. Criterion 1 **fails 1 of 8**, so
      the flag ships off (pre-registered user decision).
      Per genome, `uniform` window (molecules of 1,110 planted reads):

      | genome | genus / identity | kallisto | alignment |
      |---|---|---|---|
      | DRR140164_NODE_2 | alpha 0.969 | 184 | 1067 |
      | SRR7167047_NODE_1 | alpha 0.876 | 425 | 1071 |
      | SRR8862005_NODE_7 | beta 0.937 | 633 | 1059 |
      | SRR2037085_NODE_5080 | beta 0.452 | 499 | 861 |
      | SRR2037083_NODE_2038 | he 0.740 | 720 | 949 |
      | SRR6316308_NODE_9 | samek 0.182 | 50 | 121 |
      | SRR2037083_NODE_1170 | mem 0.108 | 29 | 107 |
      | **SRR2037085_NODE_7436** | **gamma, no alignment** | **3** | **0** |

      The one failure is the genome minimap2 could not place anywhere in the
      panel: the branch's ≥80 %-matched / ≤8 %-mismatch filters exclude it,
      while kallisto's EM still assigns a trickle. Whether those 3 molecules are
      real or EM mis-assignment is unresolved and decides whether the criterion
      was fair. **Next: relax the two filters on a sweep and re-measure
      (`ANDET-09f`), then revisit the default.**
    - Window effect, reported apart from criterion 1: the 5′ and 3′ windows lie
      outside the panel's CDS-only anellovirus models, so kallisto recovers ~0
      there regardless of identity (alpha 0.969: kallisto 0/2, alignment
      1110/1110). That is the UTR gap, not divergence, and it is why criterion 1
      uses the `uniform` window only.
- [ ] `ANDET-09f` — recalibrate the alignment filters against the held-out
  plant: sweep `--outFilterMatchNminOverLread` and `--outFilterMismatchNoverLmax`
  and report recovery vs the negative's false-positive count per setting. The
  acceptance set is built and reusable
  (`viral_ref_cat42b/runs_anello_plant/`, `--reuse-per-genome`).
  - 2026-10-03, code landed (09a builder step, 09b, 09c, 09d flag):
    - 09a: the cat42b index **is built** (`viral_ref_cat42b/build/anello_star/`,
      73 MB, 2026-10-03) (reconciled 2026-10-04).
    - 09d's default stays provisional until 09e.
    - Found by the real-STAR test: STARsolo drops homopolymer UMIs and then
      writes `CB:Z:-`/`UB:Z:-`, so those reads never count as molecules.
      Barcodes are corrected against kb's own on-list (`chemistry.onlist_path`)
      so molecules line up with the kallisto matrix. With no list, STARsolo
      keeps the raw CB.
  - **Held-out set (2026-10-03).** Modha et al. 2025, 829 genomes, TPA
    BK068993-BK069821. **0 of 829 are in the panel** (checked by accession),
    and whole-genome identity to the nearest panel genome puts **788 below
    85 %** — where an exact 31-mer cannot survive (0.85**31 ~ 0.007). Two are
    100 % identical to a panel genome under another accession.
    - Trap: the repository's `Modha_contigs.fas` is **ORF1 only** (all 829
      lengths equal the metadata's `ORF1_len`). Planting from it would sample
      just the hypervariable ORF1 and leave the 5'/3' windows undefined.
      Sequence and CDS coordinates both come from
      `Modha_genomes_annotated.gbk` instead.
    - Selection spans 6 genera x 5 identity bands, preferring full-length
      genomes: alpha >=95 (0.969) and 85-90, beta 90-95 and <85, gamma
      no-alignment, plus samek/mem/he <85 (the under-sampled genera).
  - Acceptance is scored by `scripts/anello_acceptance_report.py`: a
    per-genome assay (planted reads pulled back out by read ID, each
    (genome, window) set run through `kb count` and the branch's own STAR
    command, so recovery is a plain fraction) plus the three whole-arm runs
    (`viral_ref_cat42b/anello_acceptance.sbatch`).

## WP4H — Comprehensive human-virus catalogue (new 2026-09-27)

Objective: a generated, frozen human-host viral catalogue with real gene
structure and deterministic names. No comprehensive index exists: the widest
build (`_misc/viralscan_panel_ref_genomic`) holds ~2,216 genomes, ~88 % of them
anelloviruses, so ~99 other species.

- [x] `CAT-01` — `build-ref` discards the real GTF (`build_reference.py:609`) and
  turns every non-anellovirus accession into one `{acc}_gene1` gene (`:704`,
  `:930-947`), which would silently disable gene programmes on any natively built
  index. Use the `_genbank_to_gtf` output; emit `exon` rows for CDS-only
  features; placeholders only for CDS-less records.
  Done 2026-09-27. `index_gtf_by_seqname()` splits the merged GTF that
  `fetch_reference` already returns, and `viral_gtf_block()` picks the best
  annotation per accession: the packaged anellovirus catalogue, then the real
  NCBI CDS structure, and the whole-genome placeholder **only** for records
  carrying neither. The build logs the three counts and warns, naming
  accessions, whenever a placeholder is used, so a placeholder-heavy index is
  visible instead of silent. Measured on the real NCBI cache: HPV16
  `NC_001526.4` now yields **9** genes and HPV18 `NC_001357.1` **8**, where
  both previously collapsed to one `{acc}_gene1` bucket; `NC_002076.2` still
  takes the catalogue path with its 3 `TTVgp` genes. The existing
  `test_combines_mocked_host_and_viral_reference_without_kb_ref` asserted the
  placeholder `NC_045512.2_gene1`; it now asserts the real gene survives, which
  is the behaviour change.
- [x] `CAT-02` — **closed 2026-10-03.** The frozen `broad-discovery` list is
  now `src/viralscan/data/broad_discovery_accessions.tsv`: 4,397 rows
  (shipped 2,345, max 1,783, broad 268, legacy 1). It is a catalogue/list
  only. Building from it (`build-ref --profile broad-discovery`) is
  `REF-03`. Original row: `extras/build_virus_catalog.py` → `src/viralscan/data/virus_catalog.tsv`,
  cache-first via `ncbi_fetch.fetch_genbank()`: NCBI Virus RefSeq complete
  genomes with human host ∪ bundled panel ∪ 2,042 anelloviruses ∪ 16 HPV
  genotypes ∪ SARS-CoV-2. Becomes the frozen `broad-discovery` list (`REF-13`).
  Generator and seed catalogue done 2026-09-27; the Tier 1 union is `CAT-09`.
  Cache-first via `fetch_genbank()`, reusing `_locus_fields` and
  `_source_qualifiers`; taxonomy comes from the ORGANISM lineage. The four
  curation columns (`tier`, `persistence_class`, `risk_class`,
  `inclusion_rationale`) are carried forward on a re-run so a regeneration
  never drops a human decision, and `--check` reports drift without writing.
  **Seed catalogue: 2,215 accessions, 204 species, 30 families, 0 failures**
  (170 flatfiles newly fetched, the rest from cache). It detects 16 segmented
  species, including influenza A (8), influenza B (8), influenza C (7),
  Rotavirus C (11) and Rotavirus A (10) — the data `CAT-03` needs. Shipped in
  the wheel, sdist and Docker context (702 KB).
- [x] `CAT-03` — names from the catalogue (genome-scoped prefix → species/genus);
  segmented viruses grouped; target 0 % unnamed gene IDs (13.4 % today on the
  covid index).
  Done 2026-09-27. `src/viralscan/virus_catalog.py` loads the generated
  catalogue and exposes `catalog_name_map()` keyed by accession, versioned and
  bare, so the boundary-aware prefix rule already in `virus_grouping` resolves
  `{accession}_{gene}` with no new matching logic. `merged_name_map()` layers
  it under the anellovirus genus map and `VIRUS_NAME_MAP`, which keep
  precedence so no existing output is renamed; `detection.py`, `umap.py` and
  `evidence.py` now use it. A missing catalogue returns `[]` rather than
  raising, so an older checkout keeps working.
  **Measured on the 2,313-record reference:** unnamed gene IDs 173 → **0
  (100 % named, the target)**, and distinct groups 182 → 107 as segments and
  strains collapse into their species. Influenza A's 8 segments now group into
  one "Influenza A virus" instead of 8 separate viruses; before the change the
  same 8 gene IDs produced 8 rows.
- [~] `CAT-04` — risk classes: exclude human endogenous retroviruses; flag
  integrated ciHHV-6, `EVE_RISK_GENERA`, and vector/reagent contaminants.
  - 2026-10-10 (wave 4 reconcile): the **flags are done** in `f4320d2`: catalogue `role` (`contaminant` AAV1-8 and SV40,
    `decoy` MLV/XMRV, `endogenous` HERV-K113) and `risk_class` feed `reference_risk_flags`
    (`virus_grouping.py:261-270`: `endogenous_reference`, `vector_reagent_reference`, `endogenous_overlap_possible`,
    `iciHHV6_possible`) into `viral_summary.tsv` and the report. **Open, needs a decision/data:** HERV exclusion (HERV-K113 is
    `broad`, not shipped, so nothing to exclude yet), a CAR-T/lentiviral vector filter (no vector sequences in the repo), and a
    real ciHHV-6 integration test (today the flag is name-based).
- [x] `CAT-05` (duplicate row) — **merged into the `CAT-05` row under "WP4H
  external-asset adoption" (reconciled 2026-10-04).** Its residue (resolving the
  NC_038359.1/AB303562.1 pair; `build_bundled_panel_ref.py` reaching the guard)
  is tracked there and in package plan M1.4. Original text: duplicate guard
  (`validate_reference_records`) on every build path, including `scripts/build_bundled_panel_ref.py` (the `NC_002076.2`
  duplicate broke `kallisto index` once already).
  Second case found 2026-09-27 and still shipping: `NC_038359.1` and
  `AB303562.1` are the RefSeq and GenBank copies of one Gammatorquevirus
  genome. Both have canonical leave-one-out capture of exactly 1.000000, which
  is the detector — a genome whose every 31-mer is also contributed by some
  other panel member is by definition redundant. The guard should key on
  sequence identity, not accession, because these two differ by accession.
- [x] `CAT-06` — the index manifest carries its own viral GTF and catalogue and
  `analysis.py` reads them, so `-gtf` no longer silently drops the panel and
  nothing depends on the unregistered Zenodo DOI (`REF-11`).
  - 2026-10-10 (wave 4). User decision: "read IDs, pin catalogue". `analysis.obtain_gtf`
    now reads the raw `viral_gene_ids` of `<index>.build_manifest.json` (versions kept, so
    they match what the GTF glob wrote to `log/analysis.txt`; `load_build_manifest` strips
    them and is not used), unions any `-gtf` gene IDs, and skips `ensure_viral_data`. An
    index without a manifest keeps the old panel path, warning and error. The manifest
    gains an OPTIONAL top-level `catalogue_sha256` (sha256 of the packaged
    `virus_catalog.tsv`; no `MANIFEST_SCHEMA_VERSION` bump, old manifests read as before),
    and `build_identity_table` warns when the packaged catalogue no longer matches it.
    The GTF itself is not embedded in the manifest: the IDs suffice for `analysis.txt`.
    Remaining dependency on the Zenodo DOI: manifest-less indexes (`-i/-t` built by hand).
- [ ] `CAT-07` — diversity-aware representatives for high-diversity families,
  chosen by leave-one-out k-mer capture (the `REF-01` method), not one exemplar.
- [ ] `CAT-08` — build with native `viralscan build-ref` and measure leave-one-out
  confusability per family: where reads from non-indexed isolates land. Also the
  valid `HPV-11` test.

### WP4H expansion — breadth and measured diversity (new 2026-09-27)

Asked for: cover as many viruses as matter for infection and reactivation, and
give SARS-CoV-2, influenza and the torque teno viruses real strain diversity
instead of one exemplar.

**Census of the current reference**
(`references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.fa`,
measured 2026-09-27):

| | |
|---|---|
| Records | 2,313 |
| Anelloviridae + Gyrovirus | 2,041 (88 %) |
| Everything else | 272 records → 107 distinct species |
| …of which HHV-6B "pseudocontig" gene fragments | 97 records, not a genome |

Non-anellovirus content is therefore ~175 records over ~106 species, and nearly
every species has exactly one genome.

**Absent entirely** (checked by accession and by name):
- **SARS-CoV-2** — no `NC_045512`. Only SARS-CoV-1, MERS and HCoV-229E.
- **HIV-1, HIV-2, HTLV-1, HTLV-2** — the canonical latent retroviruses.
- **HCoV-OC43, NL63, HKU1** — 3 of the 4 seasonal coronaviruses.
- **hMPV, bocavirus, influenza D, TSPyV, HPyV6/7.**
- **HPV18** and every high-risk type except HPV16. The 16-genotype WP4E
  catalogue exists in `hpv_genes.tsv` but is in no index.

**Present but single-strain:** influenza A is 8 segments of one 1934 lab strain
(A/Puerto Rico/8 H1N1); influenza B one 1940 strain; influenza C one 1950
strain. No H3N2 and no circulating isolate.

**Decisions (2026-09-27, user):**
1. **Two-tier catalogue.** Tier 1 = every human-host RefSeq virus species, one
   genome each, the breadth floor. Tier 2 = a literature-curated
   persistence/reactivation set that gets deep strain diversity.
2. **Diversity is measured, not quota'd** — isolate counts set by leave-one-out
   31-mer capture (the `REF-01` method).
3. **Bundle sequences in the package** (self-contained, offline). This runs
   against PR 8's "move data off the package" direction; `CAT-13` is the gate.

- [x] `CAT-09` — **closed 2026-10-03** (user grill: eukaryotic, one per
  species, catalogue only). Virus-Host DB human-host RefSeqs (1,496 taxa)
  were filtered as follows:
  - phages out: 1,434 taxa left;
  - not yet catalogued: 1,026 taxa / 1,296 accessions, all fetched, 0
    errors;
  - host re-check: 1,191 kept and 105 sent to
    `analysis/cat09_sweep/review_host.tsv`;
  - collapse on NCBI Taxonomy species-rank taxids: 232 species, of which 68
    were already catalogued;
  - result: **164 species / 268 accessions added as `panel=broad`** (not
    indexed).

  Old rows are byte-identical. The script is
  `extras/cat09_human_host_sweep.py` (taxid→species cache in
  `analysis/cat09_sweep/species_taxids.tsv`); the merge is 357cb9e.
  - Defaults:
    - segments of a human-host taxon are all-or-none;
    - an organism or common name that is already catalogued is skipped, to
      avoid a display-name collision ("Norwalk virus");
    - phages are filtered on the GenBank lineage too.
  - Original row: Tier 1 breadth floor: NCBI Virus RefSeq complete genomes with
  human host, one representative per species. Closes SARS-CoV-2, HIV-1/2,
  HTLV-1/2, OC43/NL63/HKU1, hMPV, bocavirus, influenza D, TSPyV and HPyV6/7 in
  one step.
  **First batch done 2026-09-27**, closing every gap the census named. Each
  accession was resolved by a live NCBI `esearch`/`esummary` lookup, never from
  memory, per the warning in `docs/reference_panel_research_2026-09-27.md`:
  SARS-CoV-2 `NC_045512.2`; HIV-1 `NC_001802.1`, HIV-2 `NC_001722.1`; HTLV-1
  `NC_001436.1`, HTLV-2 `NC_001488.1`; HCoV-OC43 `NC_006213.1`, NL63
  `NC_005831.2`, HKU1 `NC_006577.2`; hMPV `NC_039199.1`; bocavirus
  `NC_007455.1`; TSPyV `NC_014361.1`, HPyV6 `NC_014406.1`, HPyV7 `NC_014407.1`;
  simian foamy `NC_001364.1`; influenza D `NC_036615.1`–`NC_036621.1` (7
  segments); and the 16 WP4E HPV genotypes including HPV18 `NC_001357.1`.
  Catalogue 2,215 → **2,249 accessions, 232 species, 31 families**; influenza D
  groups 7 segments into 1 virus. Remaining for this row: the systematic
  human-host RefSeq sweep beyond the ViralZone-derived list.
  **Cross-check against the Serratus screen**
  (`sources/viral_panel/Serratus_hits_all_viruses.tsv`, 129 viruses): all 17 of
  its `reactivation_candidate` viruses were already present before this batch —
  EBV leads at 6,597 high-confidence hits, then adenovirus 2,227, HSV-1 1,224,
  HCMV 1,119. That table marks HIV and HTLV as non-candidates and omits
  SARS-CoV-2 entirely, which is why the census, not the table, drove this list.
- [ ] `CAT-10` — Tier 2 persistence/reactivation set, curated from primary
  literature with a citation per entry and checked with `verify-references`
  before commit. Each row carries a `persistence_class` (`latent-episomal`,
  `latent-integrated`, `chronic-productive`, `persistent-commensal`,
  `recurrent-lytic`). Set: the 9 human herpesviruses; polyomaviruses (BK, JC,
  MC, WU, KI, TSPyV, HPyV6/7); HIV-1/2 and HTLV-1/2; HBV and HDV; HPV; the
  anelloviruses; parvovirus B19; adenovirus; measles (SSPE); enterovirus; HCV;
  pegivirus.
- [x] `CAT-11` — promote the 31-mer capture measurement to a real tool,
  `scripts/measure_kmer_capture.py`: per-genome coverage, zero-coverage
  fraction, leave-one-out capture, P(90 bp fragment captured), using
  `sensitivity.DEFAULT_K` and `fragment_capture`. The evidence behind `REF-01`
  currently lives only in the untracked `src/viralscan/scripts/kmer3.py`, which
  hard-codes institutional paths and must not be committed. Positive control
  already validated: EBV `NC_007605.1` against itself = 100.0 % over 144,283
  31-mers.
  Done 2026-09-27 as `scripts/measure_kmer_capture.py`, with 12 unit tests and a
  `--self-check`. It reproduces every `REF-01` coverage number exactly: EBV
  144,283 31-mers; bundled-panel median 0.00 % and 85.80 % zero-coverage;
  Betatorquevirus 98.4 %; the five minor genera 100 %; expanded-panel
  leave-one-out 20.55 %. **One number does not reproduce and is corrected:**
  `REF-01`'s "P(90 bp fragment captured)=1.0000" for the expanded panel is the
  analytic independence model; walking real windows measures **0.5097**,
  because shared k-mers cluster in conserved blocks instead of scattering.
  `.living/findings/` F-011 has the per-genus table. `REF-01`'s decision is
  unaffected — the expanded panel still takes zero-coverage from 85.8 % to
  0.15 % — but the panel misses about half the reads of a strain it does not
  contain, so 2,042 genomes is not the finish line, and `CAT-12` gates on the
  measured value.
  **Corrected 2026-09-27 (same day): the metric measured the forward strand
  only.** `kmers()` built single-strand k-mer sets, but kallisto indexes
  *canonical* k-mers (a k-mer folded with its reverse complement), so a panel
  genome deposited in the opposite orientation to a target scored as a miss
  when the real index would match it. `--strand canonical` is now the default
  and `--strand forward` reproduces the superseded numbers. Effect on the
  2,042-genome anellovirus panel: median leave-one-out 0.2055 → **0.2162**,
  median p_fragment 0.5097 → **0.5349**, panel 31-mers 4,888,291 → 4,843,779.
  The aggregate shift is small but the per-genome distribution is not — 350
  genomes improve, **0** degrade, 84 by more than 0.10 absolute p_fragment and
  55 by more than 0.25. The worst case, `MH649023.1`, went from 0.0329 to
  0.9514: the old metric called it essentially undetectable. The affected
  accessions cluster in the `MH648xxx`/`MH649xxx` submission block, which is
  deposited antisense to the rest of the panel. Every per-genus number below
  and in F-011 is therefore a *lower bound* until regenerated.
- [ ] `CAT-12` — *Scope (reconciled 2026-10-04): the SARS-CoV-2 and influenza bullets are
  superseded by `CAT-21`/`CAT-22`; an experiment, deferred under the package
  plan.* Set isolate counts greedily against a stated bar: **measured**
  leave-one-out P(90 bp fragment captured) ≥ 0.95 **and** zero-coverage genome
  fraction ≤ 5 %. Use the measured `p_fragment`, never the analytic
  `1-(1-c)^60` — see `CAT-11` and F-011. Measured baseline for the current
  anellovirus panel (2026-09-27): overall P = 0.5097, and per genus
  Alpha 0.788, Het 0.838, Beta 0.498, Gamma 0.332, Samek 0.255, Mem 0.247,
  Gyro 0.064 — so every genus is below the bar today and Gyrovirus is far below.
  - SARS-CoV-2: RefSeq plus one per WHO variant lineage.
  - Influenza A: one strain per relevant subtype, **all 8 segments each**, modern
    isolates beside the 1934 reference; influenza B both lineages; C and D one each.
  - Anellovirus: the 2,042 already clear P = 1.0, so the target here is the
    zero-coverage fraction and genus resolution. Add genomes for the
    under-represented genera — the 2026-09-27 held-out test put Gammatorquevirus
    at 8–40 % read recovery with no near neighbour, against 73–88 % for a
    Betatorquevirus with one at ~94 % identity.
  - HIV-1 per group/subtype; HPV, HBV, HCV, enterovirus, adenovirus by the same rule.
- [x] `CAT-13` — **package size gate.** `src/viralscan/data/` is 8.1 MB today and
  - **Done 2026-10-04:** `scripts/check_dist_size.py` (60 MB, the limit this row
    cites) runs after `python -m build` in release.yml; 3 unit tests. Its
    first real measurement happens in CI (user).
  195 of its 201 files are GTFs already gitignored (`.gitignore:99`) for size.
  Bundling sequences pushes the wheel toward PyPI's 60 MB project limit, and
  nothing in-repo would catch it: `check_ship_scope.py` is a path allowlist and
  `release.yml` has no size check. Bundle gzipped, add a wheel/sdist size
  assertion so CI fails loudly instead of the upload, and measure before `WP10`.
  If exceeded: keep Tier 2 bundled and move Tier 1 sequences to the fetch path,
  shipping the full TSV either way.
- [~] `PANEL-01` — (2026-10-08, user request) promote as many human-relevant viruses as possible into the shipped
  panel; plan approved 2026-10-08. Decisions: relevance = catalogue host field plus a curated list; new GTFs are
  **bundled** (tracked, gzipped, size-gated; note 0 of the 324 current GTFs are tracked and the wheel ships
  only TSVs, so bundling needs a `package-data` entry and a `.gitignore` exception); vector/reagent viruses
  (AAV, SV40, MMLV) ship labelled as contaminants; an NCBI `Homo sapiens[Host]` census is a gap check.
  - **WP1 done 2026-10-08 (candidate table).** `scripts/panel_candidates.py` writes
    `analysis/panel_expansion/candidates.tsv`: 2,349 candidates (2,052 catalogue `max`/`broad`/`legacy` rows plus
    297 accessions of `evonk_candidates.tsv`); 1,491 excluded (1,406 anellovirus `max` rows, 76 HPV RefSeq
    aliases, 9 sequence twins) and 858 left: H1 380, H2_species 51, H2_name 47, proposed 37, unreviewed 343.
    Nothing is promoted; `unreviewed` is never added automatically. Deterministic output.
  - **WP1 census done 2026-10-08.** `scripts/panel_census.py` (cached; offline rerun): NCBI has 401 RefSeq viral
    records with a Homo sapiens `/host`: 28 shipped, 86 `max`, 163 `broad`, 1 `legacy` and **123 not in the
    catalogue** (8 also in the colleague list; 115 census-only: human parvovirus 4, polyomavirus 9, Astrovirus
    VA1, Zika `NC_035889`, HCV genotype 7, Marburg, enterovirus A114, Coxsackievirus B3, Norovirus GIV and more,
    plus phage/bacterial records from wrong host qualifiers). The table is now 2,464 candidates, 973 not
    excluded: H1 380, H2_species 80, H2_name 50, census_host 83, proposed 37, unreviewed 343. 4 tests.
  - **External review done 2026-10-08** (Biomni; dispositions in `analysis/panel_expansion/biomni_review_2026-10-08.md`).
    Changes the order: a separability check on the candidate superset (WP1b) comes before curation; adds a
    per-family redundancy cap, a host-virus mixed-EC budget, an explicit reagent set (phiX174, AAV, SV40, MMLV,
    Ad5...), machine-actionable `role` values, a top-down pool (ICTV VMR + Virus-Host DB + the census), and a
    frozen manifest with a before/after count diff. Unverified in code: origin-spanning k-mers of circular genomes.
  - **WP1b done 2026-10-08.** `scripts/panel_kmer_sharing.py` (+ `tests/test_panel_kmer_sharing.py`, 2 pass): canonical k=31
    sharing of each not-excluded candidate against the panel, writes `analysis/panel_expansion/kmer_sharing.tsv`
    (`frac_unique`, `frac_not_in_panel`, `frac_sibling`, `frac_other_group`; "group" = `sibling_group`, else genus,
    else species, with `group_basis` shown). 973 candidates vs the built panel
    (`$SITE_ROOT/viral_ref_final/build/viral.fa`, 2,343 records): median
    `frac_not_in_panel` 1.0 in every source; **37 have <5 % outside the panel** (25 census-new, 11 catalogue-max,
    1 evonk-new; mostly RefSeq NC_ copies of panel GenBank records the twin check could not see without FASTAs)
    and 62 have <50 %. The unmasked NCBI-cache panel gives the same numbers (max difference 0.0009). 349 candidates
    have no genus, so their sibling share is understated. Candidate FASTAs are cached under
    `~/.cache/viralscan/panel_candidates/`. Top-3 partners per low-novelty candidate go to `kmer_partners.tsv`: every
    HPV row's partner is the same-type panel genome, OC43 -> panel OC43 `NC_006213`, Betacoronavirus England 1 -> panel
    MERS-related `NC_019843`, HHV-7/8 -> their RefSeq; the "other group" flags are label differences, not cross-type
    sharing. **The 37 are now excluded as `kmer_twin_of:<panel id>`** in `candidates.tsv` (973 -> 936 not excluded);
    reruns of `panel_kmer_sharing.py` -> `panel_candidates.py` are byte-stable and the table still reproduces without
    the new input. 3 new tests. Grey zone left for curation: 8 more SARS-CoV-2 genomes (5-12 % outside the panel),
    HCMV `NC_001347` (0.31), monkeypox `NC_063383` (0.12, partner `NC_003310` shares 0.84), Micro-TTV, TT virus
    sle1957, ectromelia (mouse pathogen). Max-tier Biomni review (`tsk_0127DgImaFOLeS4lyvuYNoBL`) says: keep k=31 as a
    triage screen only; twin < 0.02, grey zone 0.05-0.5, separable > 0.5; the check lacks the human host reference and
    positional unique-k-mer density (both open); verify the shipped index used k=31; decide whether strain records can
    share one gene ID before fixing thresholds.
  - **WP1b follow-ups from the second Biomni review, 2026-10-08.** `panel_kmer_sharing.py` now adds `length_ratio` and
    `frac_of_partner` to `kmer_partners.tsv` (twin vs fragment) and, with `--baseline`, scores each panel genome against
    the rest of the panel (`kmer_panel_baseline.tsv`; panel `frac_unique` 5/50/95th percentile 0.46/0.82/1.0,
    `frac_other_group` 0/0.07/0.34). All 37 `kmer_twin_of` rows have length ratio 0.9-1.1 and share >= 0.9 of their
    partner, so no fragment reason code was needed; `FR751515` (Micro-TTV, length ratio 0.34 vs `MZ286008`) is the one
    fragment in the grey zone. 105 of 973 candidates exceed the panel's 95th percentile of `frac_other_group`, but
    that number is not trustworthy until groups come from ICTV (open). `candidates.tsv` unchanged by this step.
  - **ICTV-first grouping done 2026-10-08.** `panel_kmer_sharing.py --vmr <VMR_MSL41.v1.20260729.xlsx>` (from
    ictv.global/vmr; stdlib reader, not committed) groups by catalogue `sibling_group`, else ICTV genus (by accession,
    then by cleaned virus/species name; a name that maps to two genera is dropped), else catalogue genus, else species.
    Accession lookup alone resolves 0 of the genus-less candidates (the VMR lists exemplars only), so the name match does
    the work: of 973 candidates 628 get `ictv_genus`, 180 `genus`, 165 stay `species` (HPV types, astrovirus MLB/VA, Jorvi,
    PoSCV...; not guessed). Dry run into the scratchpad: candidates above the panel's 95th-percentile `frac_other_group`
    fall 105 -> 73 (37 of those 73 are still species-basis). `kmer_*.tsv` regenerated with `--vmr` (37 `frac_not_in_panel` < 0.05, unchanged; grouping does not affect it)
    1 new test.
  - **WP2 pool started 2026-10-08.** User: pool top-down from Virus-Host DB (+ ICTV VMR), no redundancy cap. `scripts/panel_pool.py`
    (+ `tests/test_panel_pool.py`, 3 pass) reads the Virus-Host DB table (`viral_panel_max_2026-09-28/work/virushostdb.tsv`, not committed) and
    writes `analysis/panel_expansion/pool_vhdb.tsv`: 1,496 viruses with host taxid 9606 -> 219 in the panel, 339 candidates, 38 excluded
    candidates, **900 gap** (by `tier` hint: 18 disease-named, 480 RefSeq/UniProt-evidence, 402 literature-only). The VMR has no host
    field, so it adds no relevance signal here; family/genus come from the database lineage. Gap is by accession, so GenBank-vs-RefSeq
    aliases (HCV `AJ132997`, poliovirus 1-3, rhinovirus B/C) show up as false gaps, and 236 Caliciviridae + 31 Microviridae + 19 Inoviridae
    gaps are mostly gut-phage/metagenome host noise. Not curated yet: nothing here is promoted.
    Biomni max review started 2026-10-08 (`tsk_0117fBd4IZfPHBa7JCcqtRIa`, project `prj_011fdGq9tp9AeHdJuQxbSZfD`): verdict per gap group plus the
    transplant/immunosuppression reactivation set the user added (herpes-, polyoma-, adeno-, hepatitis, anello-, respiratory and enteric viruses; which are
    missing from both tables). Inputs were two collapsed tables; the 'covered' one was condensed by hand (parenthetical labels are mine).
    **Result 2026-10-08:** of the 170 gap groups 58 alias-of-covered, 53 not human-relevant (phage/metagenome, animal viruses), 24 redundant,
    24 promote-if-distinguishable, 11 promote (La Crosse, Bayou, Toscana, PIV4b, Rocio, KFDV/Alkhumra, B virus, Ravn, Guama/Caraparu); segmented ones
    need all segments. Its "missing" and "inverted" calls were against my two tables only and were wrong against the real panel: HHV-6B is shipped
    (`AF157706`; `NC_000898` is the excluded twin), LCMV is shipped (`NC_004291/4`), West Nile is shipped (`NC_001563`) with `NC_009942` a `max` candidate,
    Puumala and Andes are `proposed` candidates. Still valid: HERV-K113 candidate belongs in an endogenous tier; iciHHV-6 flag; CAR-T lentiviral-vector
    k-mers give false HIV-1 calls; AAV2 should be relabelled contaminant; anellovirus example `NC_043415` is a patent sequence (audit anellovirus
    provenance); KIPyV is the `Betapolyomavirus tertihominis` row; resolve aliases by taxon ID, not accession. Files (`gap_verdicts.tsv`,
    `transplant_relevance.tsv`, `summary.md`) are in the Biomni task, not yet saved in the repo.
    **Curation round 1 drafted:** `analysis/panel_expansion/curation_round1.tsv` (34 rows: 29 ACCEPT, 3 DEFER, 2 REJECT, each with a reason and a
    `before_promoting` step); the user flags disagreements in `your_flag`, then the accepted rows go into the catalogue edit (nothing edited yet).
    **Round 1 accepted without changes (user, 2026-10-08); accessions resolved and k-mer screened 2026-10-08** (`round1_screen.tsv`, 34 new genomes vs the
    built panel with ICTV grouping). The screen reversed 3 rows: La Crosse `NC_077808-10` (84-91 % already in shipped `NC_004108-10`), CVA16 `U05876`
    (identical k-mers to shipped `NC_001612`) and Ravn `DQ447649` (identical to candidate `NC_024781`) are aliases, now REJECT; B19 genotype 2 DEFER (only
    partial references). Sheet is now 26 ACCEPT / 4 REJECT / 4 DEFER. Everything accepted is 84-100 % new vs the panel except Toscana (57-77 %) and the
    KFDV/Alkhumra pair (Alkhumra shares 85 % with the complex, so they need one sibling_group). Resolved: HAdV-52 `NC_006879`, HEV gt4 `LC646471` / gt7
    `KJ496143` / rat `NC_038504`, HDV gt3 `L22063`, Puumala L `NC_005225`, HBV A/B/C/E/F/G/H, Toscana 1812V `PV891900-2`. HCV genotypes 2-7 are already
    candidates.
  - **WP2 catalogue edit done 2026-10-08 (round 1).** 36 accessions added through `extras/build_virus_catalog.py` and 30 more flipped
    `max`/`broad` -> `shipped` (66 in all, 4,397 -> 4,433 rows) by `scripts/panel_promote.py` (+ `tests/test_panel_promote.py`, 2 pass) from
    `analysis/panel_expansion/promotions_round1.tsv` (accession, curation id, who decided); `broad_discovery_accessions.tsv` re-synced (the
    CAT-09 test pins it to the catalogue). Existing rows unchanged. Full suite 1,920 pass. KFDV `NC_039218` and Alkhumra `JN860200` share
    `sibling_group` KFDV in the catalogue. **Correction to the sheet:** the HBV/dengue/PIV4 sibling groups I proposed do not fit
    `SIBLING_VIRUS_PAIRS`, which needs near-identity (HBV genotypes share 12-30 % of k-mers, PIV4a/4b 3.5 %); not added.
    Not done: AAV/SV40/HERV-K113 (`role` values and report plumbing, WP2 item), HHV-6 iciHHV-6 flag and CAR-T vector filter (WP3), GTFs for the
    new accessions (WP4); the shipped index does not contain the 66 until the WP5 build, which is yours.
  - [~] `PANEL-01` WP4 — **GTFs for the 66 promoted accessions** (2026-10-08). Before this, none of the 66 had a GTF in `src/viralscan/data/` (369
    non-anellovirus shipped rows: 302 had one by name). Acceptance: the builder's `--strict-reconciliation` dry run reaches zero unexplained misses; corpus
    manifest regenerated. Blocks the WP5 build.
    - **Step 1 done:** all 66 GTFs generated from the cached flatfiles with `ncbi_fetch._fetch_one` into `~/.cache/viralscan/ncbi/<accession.version>/`:
      503 genome-scoped gene IDs (`HE974370.1_Pol`), none shared across the 66 or with the existing corpus, circular origin-spanning features flagged on HBV.
      **3 are whole-genome placeholders because their GenBank records carry no CDS at all:** Alkhumra `JN860200`, Puumala L `NC_005225`, CVA24 `D90457`
      (the single-`*_gene1` shape of the F-021/F-022 artefacts, cf. `CAT-38/39`). Recommendation, needs the user: exclude those three with a named decider
      until a CDS model exists (longest ORF checked against an annotated sibling). KFDV, Puumala S/M, EV-A71 and EV-D68 already cover the same viruses.
    - **Steps 2-3 done:** the other 63 are in `src/viralscan/data/` as plain `<Species>_<ACCESSION>.gtf` (the existing convention; the builder globs `*.gtf`),
      **untracked like the other 324** (user: keep the GTFs local for now). Tracked/gzipped bundling is deferred: it needs a `package-data` entry, a `.gitignore`
      exception, the `check_dist_size.py` gate and a builder change. The three no-CDS records are not copied. `scripts/write_gtf_manifest.py` re-run: 387 GTFs,
      `--check` passes, manifest sha256 `d0d25efffc2242c7391d1b52670112db439e5913c0a11d943e8d7908be4751f9` (pinned in the study bundle). Offline coverage check:
      of 369 non-anellovirus shipped rows only those three lack a GTF seqname (`AF157706.1` is covered by its legacy pseudo-contig GTF).
      **Side effect fixed:** `tests/test_virus_grouping.py::TestRealPanelResolution` failed with the new GTFs present (521 of 4,125 gene IDs, 12.6 %, unresolved)
      because it measured the legacy token map (`anellovirus.merged_name_map`), which does not know genome-scoped IDs like `HE974370.1_Pol`. It now uses the
      catalogue-aware no-table fallback (`virus_catalog.merged_name_map`, the map the pipeline uses; MECH-A): 21 unresolved (0.5 %, all legacy `D1P*` IDs), under
      the unchanged 2 % floor. Full suite 1,920 pass.
    - **Held three resolved (2026-10-09, user: "model them").** `scripts/model_nocds_gtfs.py` writes one single-exon modelled CDS GTF each:
      Puumala L `NC_005225.1` 37..6507 (longest ORF, 2,156 aa, 69 % identity to Hantaan L `NC_005222.1`), CVA24 `D90457.1` 751..7395 (2,214 aa, 81 % to
      poliovirus `NC_002058.3`), Alkhumra `JN860200.1` 18..10274 (tblastn span of KFDV `NC_039218.1`, 96.7 %). **Alkhumra is not a clean ORF:** the record has
      at least 5 frame changes vs KFDV (longest ATG ORF is only 1,469 aa), so its gene is the aligned span, marked `note "modelled CDS ..."`; fine for read
      quantification, not for protein-level claims. All three carry a `modelled CDS` note. Manifest re-pinned: 390 GTFs, sha256 `4c5277b6…7dd11e9`
      (stdout of `scripts/write_gtf_manifest.py`; the earlier `d0d25eff...` pin in the study bundle is now stale).
    - **Step 4 done, acceptance met (2026-10-09).** Builder run with `--strict-reconciliation` on cached NCBI + Ensembl files: 2,410 viral records, **all 66
      promoted accessions and the 3 modelled ones reach `viral.fa`**, zero unexplained misses. The first run had two, both pre-existing (TTMDV12
      `NC_038359` / `AB303562`): the builder's anellovirus fetch only took the 2,020 clareaulab rows, so the CAT-05 canonical `NC_038359.1` (source
      `viralscan-refseq; CAT-05 canonical`) was never indexed. **Fix (user: "fetch and cache"):** `scripts/build_bundled_panel_ref.py` step 4b now also
      fetches any anellovirus table row no bundled GTF covers (2,021 fetched; `NC_038359.1` cached). `AB303562` is the byte-identical duplicate dropped on
      purpose, so it is the first row in `index_exclusions.tsv` (decided_by `mdmanurung`); `tests/test_index_reconciliation.py` no longer pins the allowlist
      as empty. The dry run was not stopped before `kb ref` and entered it (kb was on PATH); I killed it, nothing was built. Next full build is WP5, yours.
    - **Step 5 done:** every segment is in the panel for Toscana (S/M/L), Bayou, Puumala (S/M/L), Andes.
  - [x] **WP2 role plumbing (2026-10-09):** catalogue roles reach virus facts, the per-virus TSV and
    HTML report, with blank meaning uncurated/unknown. `roles_round1.tsv` labels 5 existing AAV references
    `contaminant` and HERV-K113 `endogenous`; existing MLV/XMRV decoys retain `decoy`. Conflicting
    nonempty roles in one virus are rejected. Focused tests: 95 pass; changed-file lint/format clean.
    Before/after comparison: 4,433 catalogue rows and 593 virus groups retain every non-role field;
    re-rendering the 17/18 saved COVID x213/x216 summary rows retains every old column. This checks
    output invariance, not a new quantification or biological validation. Final unit run: 1,929 passes;
    its sole stale registry-pin failure is repaired and 21 inventory/claim checks pass on the finished
    metadata. Output docs, builder source and claim registry are pinned to exact committed snapshots.
    Missing SV40/AAV1/7/8 records and reference promotion remain open.
  - Open: WP2 remaining curated promotions; WP3 gates (self-consistency matrix,
    EC-size gate); WP4 packaging of the new GTFs, tests, docs; WP5 cluster build and freeze (user-owned).
  - **WP5 build command (2026-10-13, Q53/Q65; prepared, NOT run).** Q65: build now on the current shipped tier (the 66 promotions + SV40/AAV1/7/8 + the
    3 modelled GTFs); the WP3 gates are the acceptance (`--strict-reconciliation`, `--ec-baseline` and the CAT-17 low-complexity gate all exit non-zero); remaining
    WP2 promotions go to protocol v2. Verified against `build_bundled_panel_ref.py --help`: `--out` (required), `--strict-reconciliation`, `--ec-baseline`
    (an `index_inspect.json`; the baseline is the committed `analysis/panel_expansion/index_inspect_viral_ref_final.json`, max EC 9,921, 0 discarded; the
    `viral_ref_final` build directory itself holds no inspect JSON), `--scale-baseline-gtf`, `--ensembl-release`, `--genome-dlist`, `--homopolymer-run-length`,
    `--lowcomplexity-kmer-mask`. Facts the command depends on:
    (a) the 393 GTFs are untracked, so a clean pinned worktree has none: copy them from the main checkout (the builder's `_check_gtf_manifest` then
    verifies them against `gtf_corpus_manifest.tsv`, 393 rows); (b) Ensembl release pinned to 116 as in the cat42d build log; (c) flags copied from the
    cat42d recipe (`--homopolymer-run-length 20 --lowcomplexity-kmer-mask --genome-dlist <GRCh38-2024-A genome.fa>`).
    **Memory and time, reconciled from PLAN evidence:** the builder docstring's `--mem=16G` is stale (it predates the D-list; no measurement backs it). `CAT-13` and
    `build_genome_panel_ref.sh` say ~64 GB / ~8 h / 8 CPUs for the host genome D-list build; the only completed full build with a genome D-list (cat42d, job
    25695872: 8 CPUs, `--mem=180G`) peaked at **49.5 GB** (MaxRSS) in **5 h 08 min** (sacct, 2026-10-13). Request `--mem=96G --time=12:00:00 --cpus-per-task=8`
    (measured peak plus ~2x margin; 64 G would leave under 30 %).
    **D-list (Q77, user 2026-10-13): cat42d parity.** Pass cat42d's exact D-list file `$M/viral_ref_cat42d/dlist.fa` (cat42b `combined.fa` + GRCh38-2024-A
    genome, sha256 `e368c0e57c749a643d51c2ad7d1b7d404b787439d4e4d40a3a82c2a41b092c65`, `viral_ref_cat42d/dlist.sha256`), not the genome alone; the
    D-list is part of what G4 freezes (`REF-08`).
    Run from a clean pinned worktree of merged `main` (the `vs_pinned/<sha>` convention), env `viralscan_test_full`:
    ```bash
    M=$SITE_ROOT
    E=$M/conda/envs/viralscan_test_full
    GENOME=$SITE_SCRATCH/malaria_bcells/data/refgenome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
    git -C $M/ViralScan fetch origin && SHA=$(git -C $M/ViralScan rev-parse --short origin/main)
    CODE=$M/vs_pinned/$SHA && git -C $M/ViralScan worktree add --detach $CODE origin/main
    cp $M/ViralScan/src/viralscan/data/*.gtf $CODE/src/viralscan/data/        # untracked local GTF corpus (393)
    OUT=$M/viral_ref_wp5/build_$SHA && mkdir -p $OUT/logs
    export NCBI_EMAIL=your.name@example.org                                    # PLACEHOLDER: use a real address
    sbatch -J vs_wp5_ref --cpus-per-task=8 --mem=96G --time=12:00:00 --export=ALL,NCBI_EMAIL \
      -o $OUT/logs/%x_%j.log --wrap "set -euo pipefail; export PATH=$E/bin:\$PATH; cd $CODE; \
      PYTHONPATH=src python scripts/build_bundled_panel_ref.py --out $OUT --strict-reconciliation \
      --ec-baseline $CODE/analysis/panel_expansion/index_inspect_viral_ref_final.json \
      --scale-baseline-gtf $M/viral_ref_cat42d/build/viral_panel.gtf \
      --ensembl-release 116 --homopolymer-run-length 20 --lowcomplexity-kmer-mask \
      --genome-dlist $M/viral_ref_cat42d/dlist.fa --threads 8; sha256sum $OUT/panel.idx $OUT/panel.t2g"
    ```
    Afterwards: read `$OUT/index_inspect.json` (9,921 is the baseline; the 66 additions can legitimately exceed it, which is the gate working),
    `$OUT/catalogued_not_indexed.tsv` and `scale_check.json`; then freeze the index hashes and regenerate `gtf_corpus_manifest.tsv`.
- [~] `CAT-14` — host cross-talk gate. Every added genome is a fresh chance to
  call human reads viral. Measure per-accession host-homologous fraction for the
  whole catalogue into `host_homology_annotations.tsv`
  (`build_reference.py:768` writes it; nothing reads it — `ANDET-02`). A
  host-only negative must produce no reported call.
  - **Software done 2026-10-12** (`reference_gates host-crosstalk` / `host-only`): the report reads
    the table and writes a per-accession TSV (union of host-locus intervals over genome length,
    `max_alignment_coverage`; an accession with no row or status other than `measured` fails and
    is never reported as 0); `--viral-fasta --host-genome` runs the existing
    `measure_host_homology` (minimap2) first. `host_only_calls(viral_summary.tsv)` returns every
    virus with `detected` true or a called infected cell; any hit fails. No fraction threshold
    was invented (the demotion cut belongs to `REF-08`).
  - Left (cluster/data): run the measurement over the 2,343-genome `viral.fa` against GRCh38
    (no shipped reference has the table yet), and run a real host-only library through the
    pipeline and `reference_gates host-only` on its `viral_summary.tsv`.
- [~] `CAT-15` — gene programmes must survive the new index: the herpesvirus
  markers in `gene_programs.tsv` resolve to real index targets, and the EBV LCL
  `SRR12682296` re-run reproduces the `PROG-07` numbers (unique layer 895
  latent / 236 productive / 311 mixed / 3,094 indeterminate). *Stale
  (reconciled 2026-10-04): PROG-07 now reports 932 called cells (526 latent); re-read
  PROG-07 before using these numbers.*
  - **Gate done 2026-10-12** (`reference_gates programmes --t2g|--gtf`; also run after every
    bundled-panel build, writing `programme_markers.json`): every catalogue virus's markers must
    resolve through `gene_programs.resolve_markers`; any unresolved marker, or a virus with no
    rows, fails. Checked on the existing `viral_ref_cat42d` and `viral_ref_final` t2g files:
    all 9 viruses / 74 markers resolve (0 unresolved).
  - Left: the EBV LCL `SRR12682296` re-run (cluster) and its comparison with PROG-07, whose
    quoted numbers are stale (not reproduced here, so none were invented).
- [~] `CAT-16` — scale check before the full build: gene count grows ~10–50×, so
  profile `analysis.py`/`detection.py` on the new GTF, and compute
  `accession_breadth` over *all* index genes of a virus rather than only
  detected ones (`ANDET-01`; `detection.py:166`, `:501-506`).
  - **Software done 2026-10-12** (`reference_gates scale`; the bundled builder now writes
    `scale_check.json` from the viral GTF before `kb ref`, new `--scale-baseline-gtf`): gene /
    transcript / contig / exon-bp counts, growth factor against a baseline GTF (warning, not a
    failure, above the row's ~50x), `analysis.py`'s `extract_gene_ids` timed on the GTF, and the
    dense-matrix memory `detection.py` needs: `n_cells x genes in the largest virus x 8 B` per
    dense copy (it densifies one virus's index genes at a time). On `viral_ref_cat42d`'s
    `viral_panel.gtf`: 6,115 genes, 230 viruses, largest `Betatorquevirus` 1,724 genes, so 3 M
    cells need 38.5 GiB per dense copy; `analysis` parse 0.05 s. Fails only on caller-supplied
    `--max-genes` / `--max-dense-gib`.
  - `accession_breadth` over all index genes is `ANDET-01` (done 2026-10-04, not re-done).
  - Left: profile a real `detection.py` run on the full-panel matrix (cluster); how many
    dense copies are live at once is not assumed.

Blocking order: `CAT-01` and `CAT-03` first — without them a natively built
index has no gene structure outside the anelloviruses
(`build_reference.py:609` discards the fetched GTF; `:697`
`_genome_as_transcript_gtf` rebuilds one `{acc}_gene1` per accession), and
nothing groups segments (`virus_grouping.py:22-59` is a prefix match), so N
influenza strains would surface as 8 × N separate viruses.

### WP4H expansion correction and additions (new 2026-09-27)

Design document: [`docs/plans/2026-09-27-viral-reference-panel-expansion.md`](docs/plans/2026-09-27-viral-reference-panel-expansion.md).
It supplies the evidence for the rows above and corrects three of them.
Six failure modes it surfaces are not in any row above and are added here.

**Corrections to existing rows**

- `CAT-12` "SARS-CoV-2: RefSeq plus one per WHO variant lineage" is **not
  executable**. RefSeq holds **exactly one** SARS-CoV-2 genome (`NC_045512.2`,
  verified 3 ways via E-utilities) and **zero RefSeq records database-wide
  carry a `pango_lineage` attribute** — it is a BioSample attribute, not a
  nuccore one, and `NC_045512.2`'s flatfile contains no `pango`/`gisaid`/
  `lineage` string. The other 9.2 M sequences are unlabelled INSDC. Also
  Pango is depth-capped at 4 nodes (measured: 0 lineages deeper than 3 dots)
  while WHO counts >3,700 JN.1 descendants, so "one genome per lineage" is
  not well-posed at any depth. Superseded by `CAT-21`.
- `CAT-12` "Influenza A: one strain per relevant subtype" under-specifies. The
  axis carrying the diversity is the **clade within a subtype**, and the
  current 8 IAV RefSeq records are all 8 segments of *one 1934 lab strain*, so
  "one per subtype" still leaves H3N2 with no circulating isolate. Superseded
  by `CAT-22`. Clade names were renamed Feb 2023 (`3C.2a1b.2a.2a.3a.1` →
  `2a.3a.1`; current emergent `2a.3a.1 (J.2.4.1)` = "K") — **do not
  hard-code `3C` strings.**
- `CAT-13` is mis-scoped as a package-size problem. The full target panel is
  ≈30 Mbp against a 3.15 Gb host genome, <1 % by mass. The real build cost is
  the **host genome D-list** (~64 GB RAM, ~8 h,
  `scripts/build_genome_panel_ref.sh`), and the real risk is semantic, not
  byte count. See the document §7.3.

**New rows**

- [x] `CAT-17` — **low-complexity k-mer masking at index-build time.** A naive
  whole-genome anellovirus panel produced **1.44 % of R2 reads as false hits
  in the EBV LCL `SRR12682296`, 100 % attributable to homopolymer/tandem-repeat
  k-mers** — measured, 6,437/6,437 captured hit reads had *zero* genuine
  anellovirus k-mer after masking. 89.5 % of hit reads carried a homopolymer run
  ≥31 bp (median longest run 49, max 97). The k-mers a poly-A tail matched were
  literally `A`*31 and its near neighbours.
  **Done 2026-09-27** as a *k-mer-space* gate, which is deliberately not an
  N-masking property: our panel is 99.99 % unmasked (785 N in 9,925,822 bp)
  yet still carried **170 pure-homopolymer 31-mers across 9 records** and 7,236
  low-complexity k-mers total. The upstream hardmasked 2,023-rep panel carries
  **0** pure-homopolymer and 3 low-complexity k-mers — a 1,268× difference that
  the N-content check alone would have missed.
  Added `low_complexity_kmer_counts` (breaks k-mers down by
  `pure_homopolymer` / `long_run` / `few_bases` / `tandem`),
  `low_complexity_kmer_fraction`, `low_complexity_report`, and two gates on
  `validate_reference_records`: `max_pure_homopolymer_kmers` (absolute count,
  the shape the failure actually has — 170 k-mers over 0.06 % of the panel is
  >1 % of reads, because poly-A reads are not rare even though the k-mers are)
  and `max_low_complexity_fraction` (per-record). Both default to *off* so
  existing callers are unaffected; `build_anellovirus_reference` sets them to
  `0` when masking is requested and `2` / `0.05` under `--no-mask`.
  `write_reference_manifest` now emits `low_complexity_kmers` and
  `low_complexity_kmer_fraction` per sequence. `scripts/build_bundled_panel_ref.py`
  — the builder that produced the **shipped** covid index and which never called
  dustmasker — now writes a viral-only `viral.fa` and gates on
  `--max-pure-homopolymer-kmers` (default 0) *before* spending ~64 GB and ~8 h
  on `kb ref`. 15 new tests; 1,325 pass.
  The error text states that **a kallisto D-list cannot fix this**: a D-list
  filters host-homologous k-mers, not self-similarity inside a viral contig.
  **Open consequence:** the covid `viral_summary.tsv` files still publish
  `Alphatorquevirus` 1,167,103 / 1,605,631 UMI, and `main` @ `4fcd748` already
  closed F-005 as accession-level artifact. That retraction is a separate row
  (`CAT-30`) and is blocked pending explicit user approval.
  Same session's positive controls: EBV `NC_007605.1` = 29,207/2,000,000 R2
  reads (1.46 %), method validated by exact synthetic recovery at 1 %, 0.1 %
  and 0.02 % abundance.
- [x] `CAT-18` — **closed 2026-10-10: not adopted.** `docs/plans/2026-09-27-viral-reference-panel-expansion.md:32-41`
  resolved §8 Q3 as one index (placeholders are 2,043 of 470,468 t2g rows; the 99.8 % failure is within anelloviruses).
  The original proposal follows for the record. **Two-index architecture.** Whole-genome pseudo-transcripts and
  real CDS transcripts for the same virus currently share one equivalence-class
  space, which is the mechanism behind the 99.8 % single-bucket failure
  (`MW455439.1_gene1`, 1,167,103/1,169,272 anellovirus UMI reported as
  "Alphatorquevirus", `ANELLO-12`). `CAT-01` fixed annotation, not the
  collision. Build `viral_gene` (real CDS only) and `viral_genome` (one
  pseudo-transcript per accession) and run two `kb count` passes on the same
  reads; write a `tier` column (`gene`/`genome`) into the v3 output schema so
  a reader can tell which space produced a number. Per-genome abundance and
  gene programmes must never read the `genome` tier.
- [~] `CAT-19` — **`kb count` discards multimapping UMIs by default in the
  scRNA-seq path** (kallisto maintainers, pachterlab/kallisto#339), so growing a
  near-identical panel converts viral reads into *dropped* UMI rather than split
  counts. Panel expansion is therefore **not monotone-good**, consistent with the
  only benchmark that measured it (van Bemmelen et al. *BMC Genomics* 2026,
  doi:10.1186/s12864-026-12874-w: ρ = −0.73 F1 and −0.57 abundance vs set size
  for >99 %-identical genomes; the largest gain came from **geographic
  restriction**, +109 % abundance / +240 % F1, not from finer clustering).
  Concretely: use `kb count --multimapping` (not `--em`); record
  `max EC size` and `number of ECs discarded` from `kallisto inspect` on every
  build and fail on regression; and run **one experiment with kallisto
  `--distinguish`** (custom workflow, zero-indexed numeric target names) before
  committing to any per-virus genome cap. k is hard-capped at 31, so genomes
  differing by ~30–70 nt share k-mers in conserved regions.
  - **Software done 2026-10-12:** the build-time record and regression gate are `CAT-40`'s
    (`reference_gates inspect`, `--ec-baseline`). **Ambiguity, most conservative reading:**
    `kb count` has no `--multimapping` (or `--em`) flag; its flag is `--mm`. ViralScan already
    recovers multi-gene UMIs with its own multimapping step (`--multimapping`, default on, from the
    BUS file), and `kb count --mm` would make kb's `adata` a mixed layer instead of the
    unique layer. So nothing was added to the `kb count` command; `tests/test_cat19_tonsil02_contracts.py`
    pins that the Snakefile carries none of `--mm`/`--tcc`/`--em`. If `--mm` is wanted, that is a
    user decision (it changes the unique layer).
  - Left: the `--distinguish` experiment (cluster run, out of scope) and the per-virus cap decision.
- [ ] `CAT-20` — *Scope (reconciled 2026-10-04): the ICTV 37-genus refresh was retracted
  (see "WP4H external-asset adoption"); what remains is the taxid re-pull plus
  SEN/TTMV.* **Anellovirus panel is probably structurally incomplete, and has
  no usable k-mer space.** `Anelloviridae[Organism]` = **42,755** nuccore
  records but genus-name queries reach only **21,283**; the ~21,472 genus-less
  legacy records (`Torque teno virus` 12,382, `Torque teno mini virus` 6,718,
  `Torque teno midi virus` 5,735, `SEN virus` 303) are structurally unreachable
  that way, and Entrez `[Organism]` matches the lineage, so a genus-derived
  collection missed them. **RefSeqViral for this family is 0**; only 178 RefSeq
  records exist and they are exactly the ICTV exemplars — everything else is
  INSDC bulk-submitted, so panel composition is currently dictated by whoever
  submitted last, not by ICTV or prevalence. Re-pull by **family taxid 687329**,
  filtering human host at the flatfile `/host` qualifier (`"Homo sapiens"
  [Organism]` returns 0 for anelloviruses and will mislead); seed with the **243
  ICTV MSL41 exemplars**; cross-check against SCANellome V2 (3,864
  representatives, Zenodo 10.5281/zenodo.7937276); add TTMV, TTMDV, **SEN virus
  (taxid 136966, outside every genus)**, and the 4 human circovirus genomes
  explicitly. Do **not** cluster (`REF-01`/`CAT-11` already showed dedup gives
  back k-mer space). Assign species by ORF1 identity to the exemplars — ICTV
  abolished genogroups, 0/2,042 records carry `/genogroup`, and the 2 that carry
  `/genotype` contradict their own organism. **Genus-level resolution is
  impossible from a nucleotide panel**: measured, **no 31-mer is shared by
  ≥1,000 of the 2,042 genomes** (max 54 k-mers reach 400–999); genera share at
  most ~44 % ORF1 *aa* identity. The only genus-sensitive route is translated
  search on ORF1 protein (kallisto `--aa`, PalmDB precedent; `--parity single`
  mandatory).
- [ ] `CAT-21` — **SARS-CoV-2 lineage policy.** RefSeq gives exactly one genome
  - 2026-10-12 (grill, user-confirmed): ship `NC_045512.2` only (option a); multi-lineage set (option b) is not adopted. No longer blocked.
  (§ `CAT-12` correction above). Decide: (a) ship `NC_045512.2` only, or
  (b) 4–6 genomes from INSDC/ENA with lineage assignment done locally
  (pangolin/UShER), which adds a GISAID DUA dependency and a CI
  reproducibility problem. Recommended set under (b): `NC_045512.2` (ancestral,
  carries UTR coverage later lineages lack) + one JN.1 descendant + one
  **BA.3.2** (only non-JN.1 branch in circulation; saltation with ~40 spike
  substitutions) + one XFG + one pre-Omicron VOC for re-analysis of 2021–2025
  data. Current WHO designations: VOI JN.1; VUMs XFG, NB.1.8.1, PQ.16.1.1,
  BA.3.2. **Blocked on a user decision (document §8 Q1).**
- [ ] `CAT-22` — **influenza clade-level budget**, replacing the "one per
  subtype" wording. IAV H1N1pdm09 6 (pre-2009 vs 2009; 3C.3a, 6B, 6B.1, 6B.2);
  H1N1 non-pdm09 2; H3N2 8 (3C.2a, 3C.2a1, 3C.2a1b.1a, 2a, 2a.1, 2a.3,
  2a.3a.1, 2a.3a.1 (J.2.4.1)); H5N1 6 (2.3.4, 2.3.4.4, **2.3.4.4b**, human-case
  genotypes B3.13 and D1.1); H7N9 3; H9N2 4; H1N2 2; H2N2 1; IB Victoria 4;
  **IB Yamagata 1 flagged `retired_from_surveillance`** (no confirmed detection
  since March 2020, dropped from NH 2026-27 CVVs, still in the Sept 2026
  Southern Hemisphere list and an LAIV component); IC and ID 1 each.
  **8 segments per strain, always.** ⚠ **Type by sequence, not organism string** —
  `OP212288` (A/Texas/61/2022, first dairy-cattle human H5N1, clade 2.3.4.4b
  B3.13) is deposited with `(H3N2)` in its organism field. Verified B3.13 =
  `PQ468757`–`PQ468764`; D1.1 = `PQ573551`–`PQ573557`. Clustering:
  **within-subtype only, 95 % nt / 85 % coverage** (`cd-hit-est -c 0.95 -aS
  0.85`) for RNA viruses, 98 % for adenovirus. Note CD-HIT was *excluded* from
  the van Bemmelen benchmark for runtime, so no benchmark-derived CD-HIT
  threshold exists; 95/85 is the closest documented viral-derep default
  (`votuderep` ANI ≥95 % + coverage ≥85 %, the CheckV method). The
  directionally relevant caution is PalmDB *Nat Biotechnol* 2025: 99 % aa
  clustering made **67.4 % of true taxa undetectable** vs 3.3 % unclustered, so
  they grouped by taxonomy instead — which is what `CAT-03` now does.
- [ ] `CAT-23` — **Serratus reactivation panel as an explicit acceptance
  criterion.** The "Serratus activation screen" is **Lareau CA, …, Satpathy AT.
  "Latent human herpesvirus 6 is reactivated in CAR T cells." *Nature*
  623(7987):608–615 (2023), doi:10.1038/s41586-023-06704-2, PMID 37938768,
  code `github.com/caleblareau/serratus-reactivation-screen` — *read*-based
  (Serratus petabase), not host-signature-based. Input panel 129 curated human
  viruses (ViralZone / Hulo 2011); reactivation criterion = DNA-genome-only,
  known latent cycles that reactivate *in vivo* (Traylen et al. *Future Virology*
  6:451, 2011), narrowed to 17 DNA viruses from Herpesviridae,
  Polyomaviridae, Adenoviridae, Parvoviridae; positive call = **≥100 reads AND
  ≥50 % mean mapped identity** per sample–virus. **Corrected 2026-09-27: the
  anchor is 17, not 18.** The in-repo source table settles it —
  `sources/viral_panel/Serratus_hits_all_viruses.tsv` has 129 rows and exactly
  **17** with `reactivation_candidate == TRUE`; **HHV-6B is not among them**
  (only `NC_001664`, HHV-6A). The upstream "17" was right and an earlier note in
  this row that said 18 was an enumeration slip on my part. **All 17 are
  already in the panel** — T2's latent-DNA
  core needs no new accessions, only real gene structure and strain diversity,
  so it becomes the *acceptance criterion* rather than new breadth. Make it a
  checkable row: `A1` in document §7.1. **T2 must be a superset, never a copy:**
  the Serratus panel contains **no SARS-CoV-2**, only 3 HPV genotypes, filed
  HHV-6 as HHV-6A only (needing a separate HHV-6B query), and TTV appears with
  206 incidental hits never flagged reactivation-relevant.
- [ ] `CAT-24` — **expected-negative arithmetic for anelloviruses and other
  low-prevalence latent agents.** TTV in an immunocompetent adult is ~10²–10³
  copies/mL plasma (*Viruses* 17(2):140, 2025, PMID 40143262) with >90 %
  prevalence (PMID 12721794), but scRNA-seq is 3′-biased and polyA-selected,
  TTV is ssDNA, and the panel has no conserved k-mer space. Measured on the EBV
  LCL: **EBV 1.46 % of R2 reads, zero genuine anellovirus reads.** So report
  anellovirus as presence-only, binary, `screening_only` (`REF-10`); attach the
  expected-value calculation to every negative instead of publishing a bare
  zero (the discipline already required for tonsil EBV, `SENS-06`); add a
  `detection_bound` column; and **never publish per-genotype anellovirus
  abundance** — 87 % of annotated genomes (1,740/1,995) still carry a single CDS
  spanning the whole genome and are competition buckets (`ANELLO-12`).

### WP4H external-asset adoption and adversarial corrections (new 2026-09-27)

**Vendored external references.** 41 files, 30.2 MiB, 0 failures, every
repo pinned by commit, under `extras/vendor_sources/` with
`VENDOR_MANIFEST.tsv` (per-file bytes + sha256) and `VENDOR_REPORT.md`.

| Source | Pin | What it gives us |
|---|---|---|
| `clareaulab/human_anellovirus_pangenome` | `3ed77e19` | **2,023 hardmasked reps** (5,927,006 bp, 284,762 N = 4.80 %, 0 pure-homopolymer k-mers), `anello_t2g.txt`, `anello_for_kallisto.gtf`, `simple_anello_metadata_V2.csv` (**3,545** rows), `vclust/{genus,species,clusters}.tsv`, ORF1 `.contree`/`.treefile`/`.faa`, and `example/expected_output/full_SRR32170409/` |
| `caleblareau/pan-viral-reactivation` | `74136de5` | `pan_virus_annotation_plain.tsv` (**724** rows, `EC/Gene/Nuccore/Virus`, 9 viruses) with standard nomenclature (`NC_006273.2`→`RL1`/`UL1`…/`US27`/`TRS1`; `NC_007605.1`→`BNRF1`/`EBNA-1`/`LMP-1`/`BARF0`/`BGLF1`…), `pan_virus_cds.fasta`, an HIV reference, and **6 real `*.kb.txt` outputs** as regression fixtures |
| `clareaulab/ad-hsv-mapping` | `53801369` | `HSV1-LATonly.fasta`, `HSV1-coding.fasta`, `VZV-coding.fasta`, `jg_NC_001806.2.gtf` |
| `yyoshiaki/VIRTUS3` | `b7873791` | `NC_007605.1_CDS_EBER12.fa` — 96 sequences = **94 CDS + EBER1 + EBER2** (named `rna-HHV4_EBER-*`) |
| `huangyh09/ViralScan` (Apache-2.0) | `d8279c37` | `Viral_GTF_maker.py` + `viral_reference/viruses_833.fasta` — **the provenance of our 195 GTFs** |

**Accession verification** against live NCBI/ENA/Datasets:
`extras/vendor_sources/ACCESSION_VERIFICATION.md` — 11 SARS-CoV-2 verified
(5 of 7 requested lineages + Wuhan-Hu-1), ~95 influenza segment accessions over
12 complete isolates, **9/9** endemic coronaviruses, **11/11** ssDNA viruses.

**Adversarial corrections to this work package (2026-09-27).** A review of my own
design found four things wrong, and two of them had already been found
independently upstream in this repo:

1. **Retracted — "genus queries reach only 21,283 of 42,755 records".** True of
   Entrez's `[Organism]` index but irrelevant here: the panel is already
   human-host filtered, and all 16 ICTV genera absent from the panel have
   **zero** complete genomes with `"Homo sapiens"[Host]` (Lambda 408/0, Eta
   256/0, Iota 139/0, Kappa 24/0 — swine, feline, canine, tupaia, pinniped).
   The ICTV 37-genus refresh is **deleted** from `CAT-20`; it buys no
   sensitivity and adds false-positive surface.
2. **Retracted — "our genus labels are unreliable".** Backwards. The 582 panel
   genomes with no GenBank genus are the ~584 the upstream ORF1 phylogeny
   **resolved**; our labels are better than GenBank's. **Do not reconcile them
   against GenBank lineages.**
3. **Retracted — my own collision measurement.** I grouped `t2g` by accession
   prefix and reported one accession carrying both real genes and a placeholder.
   Wrong key. F-013 / session `2026-09-27-001` found a second duplicate that
   also fixes *how* `CAT-05` must work: `AB303562.1` and `NC_038359.1` are the
   GenBank and RefSeq copies of one Gammatorquevirus genome, **both shipping**,
   an upstream dereplication failure. **The guard must key on sequence, not
   accession** — and it already does; see the note on `CAT-05` below.
4. **`OP212288` is not a mislabelled H5N1.** I reported it as H3N1 deposited
   `(H3N2)`. It is genuinely H3N2 (87.7 % identity to H3N2 vs 42.4 % to H5N1,
   567 aa CDS, monobasic HA0 cleavage site); the serotype lives in the
   DEFINITION isolate string. The real hazard is a **name collision**: both it
   and the A/Texas/61/2022 dairy-cattle H5N1 case are "A/Texas/61/2022". Type by
   sequence, not by organism *or* isolate string.

**Measured this session, `CAT-05` evidence.** The new gate rejected **both** the
upstream 2,023-rep panel and the shipped 2,312-record panel on the same
duplicate: `NC_038359.1` / `AB303562.1`, byte-identical sequences, present
upstream. This is an **upstream** dereplication failure we inherited, and it is
a live `kallisto index` hazard — `CAT-05` records a previous duplicate
(`NC_002076.2`) breaking the build.

- [x] `CAT-05` — **duplicate guard keyed on sequence, on every build path.**
  `validate_reference_records` already hashed sequences, and the new
  low-complexity gate made it fire in practice: run against
  `extras/vendor_sources/anello/ref/hardmasked_*.fa` it rejected the panel on
  `NC_038359.1` / `AB303562.1`. Partially landed with `CAT-17` — the guard is
  proven to catch a real duplicate. **Remaining:** decide the resolution
  (drop the RefSeq copy, or keep one and record why) and confirm
  `scripts/build_bundled_panel_ref.py` reaches the guard on every path.
- [~] `CAT-25` — **pinned environment.** `environment.yml` declares
  - 2026-10-04: environment.yml pins match the tools the G1 gate actually ran
    (kb-python 0.30.2 with bundled kallisto 0.52.0 / bustools 0.45.1,
    STAR 2.7.11b, samtools 1.23.1, snakemake 9.23.1, python 3.11).
    minimap2, blast, cd-hit, r-base and DropletUtils carry
    `# TODO(user): verify on lock`; they are not in viralscan_bench.
  - 2026-10-05: those five were resolved in `viralscan_test_full` and pinned
    to what passed the gate: minimap2 2.31, blast 2.17.0, cd-hit 4.8.1,
    r-base 4.5.3, DropletUtils 1.30.0. The old r-base 4.4.1 pin did not solve
    against the gate pins (icu conflict). `environment.tools.lock.txt` now
    lists them with no UNVERIFIED rows. The real conda-lock/pixi lockfile is
    still REL-03.
  kallisto 0.50.1 / bustools 0.43.2 / kb-python 0.28.2, but the running env is
  **0.51.1 / 0.45.1 / 0.29.5** — so every `CAT-11` number was measured off-pin.
  `blast=2.16.0` is declared (it provides `dustmasker`, **absent from PATH**) but
  `conda-recipe/meta.yaml` declares neither blast nor cd-hit. Decision: pin down
  to the declared versions and re-measure the capture bars and EC sizes.
- [ ] `CAT-26` — **freeze the current reference before the panel moves.**
  `docs/review-clear-execute-plan.md` reruns the EBV and HIV controls *"from
  identical FASTQs and the original combined reference"* and freezes reference
  hashes. A new panel must not perturb that arm; freeze a copy and treat the new
  panel as a third arm.
- [x] `CAT-27` — **Done (reconciled 2026-10-04):** the CAT-32 outcome records
  EBER-1/EBER-2 emitted after the `ncbi_fetch.py` `key == "CDS"` fix (96 EBV
  genes). *Evidence correction, same day: `data/*.gtf` is git-ignored (6bb5c64),
  so the local `Epstein_Barr_virus_NC_007605.gtf` first cited here is an
  unshipped May copy, not evidence.* They were emitted after the `ncbi_fetch.py` `key == "CDS"` fix (CAT-32
  outcome), giving 96 genes and matching VIRTUS3. Residue for `CAT-41`: their
  gene IDs read `EPSTEIN_unassigned_gene_1`. **EBER1/EBER2 audit.** VIRTUS3 ships an EBV reference of 96
  sequences = 94 CDS + **EBER1 (167 nt) + EBER2 (173 nt)**, named
  `rna-HHV4_EBER-*`. Ours may omit the two non-coding RNAs, which are the
  highest-abundance latent EBV transcripts and arguably the best latent marker.
  Check, and add from RefSeq if absent. Cross-check already passes: VIRTUS3's
  94 CDS, `HSV1-coding.fasta` 77 and `VZV-coding.fasta` 73 exactly match the
  panviral HHV4/HHV1/HHV3 gene counts — three independent repos built against the
  same RefSeq release.
- [x] `CAT-28` — **gene-nomenclature mapping.** Map
  `pan_viral_annotation_plain.tsv` (724 rows) into `gene_programs.tsv` and
  validate the markers. **Data hazards:** the file is EC-keyed so genes repeat
  (`BWRF1` ×6, `LMP-1` ×3, `US33A` ×3, `AAV2gp06` ×2, `K14` ×2), and **line 532
  is corrupted** — `530  Jvgp6  [NC_001699.1  Jcpolyomavirus` with a stray `[`.
  A join on accession drops JCV `Jvgp6` **silently**. Also records accession
  disagreements: HHV7 `U43400.1` (we use `NC_001716`), HHV8 `MK733606.1` (we
  use `NC_009333`), and HHV6B `AF157706.1` — the pseudocontig behind our 97
  "records" that are not genomes.
  - 2026-10-10 (wave 4): **done as a validation report only** (user decision: no
    mapping into `gene_programs.tsv`). `extras/cat28_panviral_marker_check.py`
    writes `analysis/panel_expansion/cat28_panviral_marker_check.tsv` (marker,
    virus, status, detail). The input's real name is `pan_virus_...`, not
    `pan_viral_...`, and it has **724 data rows** (725 lines with the header).
    The script repairs the bracket (1 accession, counted), collapses EC repeats
    per (virus, accession, gene), and matches each marker's `refseq_gene` to the
    panviral `Gene` of the same virus by **exact, case-sensitive** equality
    (near misses go in `detail`, never to `matched`).
  - **Result on the 75-row post-`PROG-13` catalogue: 41 matched, 16 unmatched,
    18 accession_mismatch.** Matched: EBV 16 (EBNA-1's two catalogue records count
    twice), CMV 9, HSV-1 4, HHV-6B 7, VZV 5. Unmatched:
    HHV-2 6 and HHV-6A 8 (neither virus is in the panviral table, which holds
    AAV2, HHV1/3/4/5/6B/7/8 and JCV only); HSV-1 `LAT` (panviral has no LAT
    gene); EBV `BBLF2`, which panviral spells as the fused `BBLF2/BBLF3`.
    `accession_mismatch`: HHV-7 8 (`U43400.1` vs our `NC_001716.2`) and KSHV 10
    (`MK733606.1` vs `NC_009333.1`); every gene name matches, only the accession
    differs, so a name-based join works and an accession-based one fails.
  - **Not a mismatch:** HHV-6B. Our panel's GTF seqname is `AF157706.1`, the same
    accession as panviral, so the 7 markers are `matched`; the detail says the
    accession is the pseudocontig, not a genome. The known HHV-6B accession
    problem is the pseudocontig itself (the 97 non-genome records), not a
    disagreement between the two tables.
  - **Decision for the user:** accept `BBLF2` -> `BBLF2/BBLF3` (and `U5/7`, `BGRF1/BDRF1`,
    `BSLF2/BMLF1`, `EBNA-3B/EBNA-3C`) as a fused-name alias, or keep the exact rule.
    The marker count it would change is 1.
- [~] `CAT-29` — **HSV-1 latency transcripts → Tier 2** (user decision).
  - 2026-10-12 (grill, user-confirmed): re-cut to the two exons (660 nt + 5,731 nt), drop the 1,956 nt intron record, then re-evaluate HSV-1 `panel_completeness`.
  - 2026-10-12 **re-cut and re-evaluation done; the "Tier 2" label waits on `CAT-10`** (evidence and tables:
    `.living/findings/hsv1-lat-exon-recut.md`):
    - The shipped index **already** holds LAT as exons only: bundled GTF `HHV1gp00s02` has two exon rows, and the built
      panel's transcript is 6,391 nt = 660 + 5,731 (`viral_ref_final/build/cdna.fa`). The 3 vendored records equal the
      NC_001806.2 slices; nothing builds from the file, so it stays byte-identical (manifest hash) and no derived FASTA
      was written. VENDOR_REPORT section 6 records the decision.
    - k-mer check (canonical k=31, 79 HSV-1 genes): exons 5,765 distinct 31-mers, 0 shared with UL9/UL41/UL44/UL54, 32.6 %
      with ICP0 and 12.1 % with ICP34.5 (the known antisense overlap; neither is a marker), 0 of 5,765 in 465,769 human cDNA
      transcripts. The intron adds 1,825 31-mers found in no exon (39.2 % shared with ICP0), which is the false-positive
      surface the user removed. Host check is cDNA only; no genome-wide host screen was run.
    - **HSV-1 `panel_completeness`: `partial` before and after** (`latency_observable_in_rna=false`, `latent` unreachable).
      The exons are still the sole latency anchor, so `gene_programs.tsv`, its generator and the docs are unchanged; new
      pin `test_hsv1_lat_is_exons_only_and_the_sole_latency_anchor`.
    - `virus_catalog.tsv` `tier`/`persistence_class` are empty on all 4,437 rows, so "Tier 2" cannot be set for one virus;
      it follows `CAT-10`. The repeat-merge row below ends "`CAT-29` stays `[!]`", which this row supersedes.
  `HSV1-LATonly.fasta` is the thing `gene_programs.tsv` calls `partial` because
  "HSV-1's latent state is unreachable by construction". **But it is not
  spliced and not contiguous**: 3 separate records — exon 660, **intron 1,956**,
  exon 5,731 — all ICP0/LAT at 118805-127151. The intron is independently
  k-mer-countable, so despite the name this does *not* restrict to LAT mRNA and
  will add false-positive surface. Re-cut the exons before adopting, and
  re-evaluate HSV-1 `panel_completeness` afterwards.
- [ ] `CAT-30` — **retract the covid Alphatorquevirus claim.** `main` @ `4fcd748`
  closed F-005 — *"All 8 detected anellovirus accessions appear in Phase A
  (multi-chromosomal GRCh38 alignment) — accession-level artifact confirmation"*,
  *"no genuine viral infection in these COVID PBMC samples"* — but
  `covid_viralscan/results/*/results/viral_summary.tsv` **still publishes**
  `Alphatorquevirus` 1,167,103.0 UMI (x213-g, x216-g.jul2), 1,605,631.0
  (x216-g), 17.63–19.05 % of cells. The science is settled; the retraction never
  happened. **Blocked pending explicit user approval** — do not edit published
  results unilaterally. `CAT-17` makes this more urgent, not less: the shipped
  panel is the one carrying 170 pure-homopolymer k-mers.

**Sequencing note (superseded 2026-09-27):** `CAT-17` shipped. Next is
`CAT-30` (user approval), then the `CAT-05` duplicate decision, then
`CAT-20`'s adoption of the 2,023 reps — which stays **gated** on the F-013/F-011
reconciliation the previous session opened, per its own instruction *"do not act
on restore-the-1,522 until this is settled."* The full dependency-ordered
sequence is in document §9.

## WP4I — Complete latent/lytic state calling (new 2026-09-27)

Objective: make `gene_programs` biologically correct and measurable. Continues
`PROG-08`/`PROG-09`.

- [x] `PROG-10` — fix vacuous tests first: `test_measured_ratio` passes for any
  counts, and one EBV-regression assertion checks the opposite of its docstring.
  Done 2026-09-27. The directional guard was `not (state == productive and
  selected_state == latent)`, unreachable once `state == latent` is asserted two
  lines earlier; it now forbids the direction its comment names (allocated
  `productive`, unique `latent`), and
  `test_directional_guard_fires_when_allocation_is_productive_only` proves that
  combination is reachable, so the guard can actually fail. The ratio test only
  checked that the catalogue has at least as many productive groups as latent
  ones, whatever the masses; it is replaced by `test_calls_depend_on_breadth_not_mass`,
  which holds the support fixed and moves mass from 1.0/1.05 to 1000/0.01 and
  requires an identical call. Neither test reproduces the real-run numbers; that
  remains `PROG-07`'s job after `PROG-11`.
- [x] `PROG-11` — catalogue biology, each change verified against primary
  literature before editing: a kinetic-class column (`latent` /
  `immediate_early` / `early` / `late`); CMV UL122/123 and HHV-6A/7 U90/U86 are
  immediate-early, not latent; EBV `BaRF1.1`, `BHRF1`, `BNLF2a/b` are early lytic;
  KSHV's latent set lists ORF16 (vBcl-2, lytic) and `partial` is hand-set despite
  three independent latent groups; exact (not case-insensitive) name matching in
  `extras/build_gene_programs.py:888`, `:896`; `panel_completeness` derived by rule.
  Done 2026-09-27, each with a primary citation: EBV `BARF1` removed (latent only
  in epithelial cancers — PMID 32708965, 39329759), which also removes `BaRF1.1`
  (the lytic ribonucleotide reductase) that the case-insensitive match had pulled
  in; matching in `_resolve` is now case-sensitive, and regenerating dropped no
  other marker (79 → 77 rows, the diff is exactly the five intended changes).
  HCMV `UL122`/`UL123` are immediate-early, not latent, and HCMV is now
  `partial` with `latency_observable_in_rna=false`: single-cell HCMV latency has
  no restricted latency programme but mirrors a late-lytic one at much lower
  levels (Shnayder et al. 2018, PMID 29535194). KSHV `ORF16` is vBcl-2, lytic,
  not vGPCR (ORF74; PMID 20860481). The docs' BARF1.2 showcase numbers carry a
  dated correction. **Deviation:** "`panel_completeness` derived by rule" is
  dropped — HCMV has two independent latent anchors (`UL138`, `UL111A`), so the
  rule would call it complete while the biology says latency is unobservable;
  the facts stay hand-set, each with its reason. **Still unverified, so not
  edited** (literature search budget ran out): EBV `BHRF1`, `BNLF2a/b` as early
  lytic; the `BcLF1` and KSHV `ORF17` notes; KSHV additions `ORF72`, `ORF71`,
  `K12`, `K15`, `K8`, `K8.1`, `ORF57`, `ORF26`, PAN and `K1`'s class; HHV-6A/7
  `U90`/`U86`; HHV-6B `U95`; HHV-7's latency set.
  - Done 2026-10-02 (`v3/prog11`, c004f47; user confirmed both decisions).
    Sources were retrieved through Europe PMC.
    - **New columns:** `kinetic_class` (latent, immediate_early, early,
      leaky_late, late, unclassified) and `kinetic_pmid`. `unclassified`
      means no source was retrieved.
    - **EBV:** BHRF1, BNLF2a and BNLF2b are removed as markers; they are early
      lytic (29864140), and BNLF2a/b overlap LMP-1. This follows BARF1.
    - **HHV-6A, HHV-6B, HHV-7:** HHV-6A U90/U86, HHV-6B U95 and HHV-7 U90
      move to productive/IE. HHV-6A and HHV-7 are now `partial` (latency not
      observable in RNA).
    - **Note fixes:** KSHV ORF17 (protease, not MTA), VZV ORF4 (IE4, not
      IE62), EBV BcLF1 (major capsid protein).
    - **Result:** the catalogue goes from 77 to 74 rows. A new test forbids
      mixed-programme overlap groups.
    - **Flagged, out of scope:** EBV BBLF4/BBLF1/BGLF4/BALF5 notes contradict
      the CAGE table; KSHV ORF71/72 etc. were never in the catalogue.
- [x] `PROG-12` — states `latent` / `reactivating` (immediate-early only) /
  `productive` / `mixed` / `indeterminate`; symmetric breadth thresholds (latent
  needs 1 group today, productive 2); a per-marker UMI floor. (Not "honour
  `non_overlapping`": `_breadth` ignores it deliberately — the overlap group is
  the unit.) HCMV needs its own handling here: because latency mirrors
  low-level late-lytic expression, a presence-based `productive` call is not
  specific either, so it needs a per-cell quantity threshold or
  `not_applicable`.
  - 2026-10-10 (wave 4): **done.** `Marker.kinetic_class` added; `reactivating`
    is in `STATES`. **Rule:** a cell is `reactivating` when it has productive
    evidence, every productive marker carrying a molecule is `immediate_early`,
    and the cell is not already `productive`/`mixed`/`latent`. It only ever
    replaces `indeterminate`, so **no PROG-07 count other than `indeterminate`
    moves** (the old `indeterminate` splits into `indeterminate` +
    `reactivating`). Tests prove it: a frozen copy of the old `_decide`, a full
    decision-table equivalence, and random matrices over every virus's real
    markers where a changed call is always `indeterminate` -> `reactivating`.
    The next `rerun-programs` on the EBV run will therefore show the 155
    `indeterminate` as `indeterminate` + `reactivating`; 526 / 67 / 184 hold by
    construction but were **not re-run on real data here** (no H5AD in this
    worktree). `docs/faq.md`'s HSV-1 paragraph is updated; its dated PROG-07
    numbers are history and stay.
  - **Parameters** (defaults reproduce today): `programme_latent_min_breadth`
    (1) and `programme_min_umi` (0). The floor is `value > 0 and value >= X`, so
    the default still counts the *fractional* values of the allocated layer
    (`adata.X`); a floor of `>= 1` would have dropped them and moved the 339
    allocated-layer latent cells. Wired like `programme_min_breadth`: defaults,
    `RunConfig` (validated), `menu`, both parsers, `scripts/gene_programs.py`,
    `docs/cli_reference.md` regenerated and `--check` clean. Summary gains
    `n_cells_reactivating`, `latent_min_breadth`, `marker_min_umi`.
  - **HCMV** is `not_applicable` for reactivating (`REACTIVATION_NOT_APPLICABLE`,
    keyed by catalogue name): never emitted per cell (IE-only evidence stays
    `indeterminate`) and `n_cells_reactivating` reads `not_applicable`. This
    names HCMV only, per the user's decision; VZV, HHV-6A/6B and HHV-7 can
    reach `reactivating` although `latency_observable_in_rna=false`, because the
    state is a positive IE observation and needs no absence claim. **Decision for
    the user:** gate every `latency_observable=false` virus instead?
  - **Catalogue facts checked:** no overlap group mixes IE and non-IE markers.
    HHV-7 has two IE groups (U90, U42), so at `min_breadth=2` it is `productive`,
    not `reactivating`. VZV's two ORF62 rows are two IE groups only because of the
    repeat copies; `PROG-13` collapses them. A latent marker's `kinetic_class`
    (VZV ORF4, HHV-6B U94 are `immediate_early`) is ignored: IE-only is defined
    over productive markers.
  - **Bug fixed on the way:** `rerun-programs --programme-min-breadth` was
    validated and then ignored (the payload came from the run's `config.yaml`).
    Its three programme flags now default to `None` and override the run's own
    config only when given.
- [x] `PROG-13` — merge 31-mer-identical repeat copies (HSV LAT/ICP0/ICP4 in
  TRL/IRL, VZV ORF62/ORF63 in TRS/IRS) to one gene_id in t2g, so they reach the
  unique layer.
  - 2026-10-10 (wave 4): **done for new builds and a CLI for existing indices**
    (user decision: post-process the t2g). `data/repeat_gene_merges.tsv` +
    `viralscan.repeat_merges.apply_repeat_merges` rewrite column 2 of the t2g
    only (every other byte, host rows and line endings included, is preserved;
    idempotent; no survivor is also a copy). Hooked into
    `scripts/build_bundled_panel_ref.py` right after `kb ref`; existing indices:
    `python scripts/apply_repeat_gene_merges.py --index-dir DIR` (`t2g.txt` or
    `panel.t2g`). Tried on a copy of the covid index's t2g: 5 rows rewritten,
    a second run 0.
  - **Verified from the local GTFs and NCBI FASTA** (canonical 31-mers of the
    strand-aware sequence; `tests/test_repeat_merges.py`, skipped when the data
    is absent). VZV is exactly identical: ORF62 `gp63`/`gp72` 3,963 of 3,963,
    ORF63 `gp64`/`gp71` 821 of 821. **HSV-1 is nested, not identical:** the
    smaller copy's k-mers are all in the larger one, which adds a tail (LAT
    `s01` 5,388 in `s02`'s 5,795, +407; ICP0 `p76` 2,655 in `p17`'s 2,682, +27;
    ICP4 `p01` 4,228 in `p15`'s 4,229, +1). The survivor is therefore the
    longer copy (`s02`, `p17`, `p15`) so no unique sequence is lost; for VZV it
    is the first listed. The test asserts containment, not equality.
  - **IDs come in two spellings** (bundled GTFs `HUM_HERP1_HHV1gp00s01` and
    `VARICELLA_HHV3_gp63`; the VIRTUS-sourced covid index, bare `HHV1gp00s01` and
    `HHV3_gp63`), so the table lists both (10 rows); an index with neither is
    left alone.
  - **Runtime readers checked.** `kb count -g` and `multimap.load_transcripts` /
    `read_ec` take the rewritten t2g as is (an EC whose two transcripts now share
    a gene yields a duplicate index, which `multimapping` already collapses with
    `set()`); `virus_identity.read_t2g` keys on gene so the merged gene is one
    gene; the build manifest tolerates extra viral IDs (`from_manifest &
    indexed`), so a manifest written before the rewrite still validates. **The
    t2g is fingerprinted, though:** `run_safety.build_run_manifest` records
    `sha256(transcripts)` per run, so `validate-run` reports
    `input_fingerprint_mismatch` for a run made with the old t2g once an
    index is rewritten in place; and `build_reference` records
    `content_files.t2g` in the manifest of a `viralscan build-ref` index, which
    then no longer matches. The bundled builder writes no such manifest. A run
    made **before** the rewrite must re-run `kb count` with the rewritten t2g
    (a new input fingerprint, so a fresh run); `rerun-multimap` alone leaves
    the copy's column in the matrix at zero. The CLI docstring says this.
  - **`gene_programs.tsv` regenerated offline** (the generator runs against the
    local GTFs): 77 -> 75 rows. HSV-1 LAT `s01` and `s02` are one row (`s02`,
    group `g1`: the copies' exons join the survivor's before overlap groups are
    computed, so an overlap group follows the merged gene), and VZV ORF62 loses
    `gp72` (which counted as a second, independent IE group). No other row
    changed. Without this, the dropped columns would fail to resolve and raise
    the PROG-14 warning on every run. **Generator caveat (fixed 2026-10-11):** it
    named HHV-6A "Human herpesvirus 6"; `extras/build_gene_programs.py` now
    uses "Human herpesvirus 6A" throughout, so a full regeneration no longer
    reverts CAT-41.
  - **Not hooked:** `viralscan build-ref` (`scripts/build_reference.py`) writes
    `t2g.txt` itself and was not asked for; use the CLI on its output. Other
    TRL/IRL and TRS/IRS pairs exist (HSV-1 RL1/ICP34.5 `p77`/`p16`, VZV
    ORF64 `gp65`/`gp70`, same length; and HSV-2 LAT `s01`/`s02`, still two catalogue
    rows) but were not in the decision; add rows
    after the same k-mer check. `CAT-29` stays `[!]`.
- [x] `PROG-14` — unresolved markers fail loudly instead of a log line
  - **Done 2026-10-04.** The gene_program summary gains `n_markers_resolved` /
    `n_markers_unresolved`.
    - Partial resolution warns, with a caveat that panel_completeness
      describes the catalogue, not the index.
    - Zero resolved now says the reference and catalogue disagree, instead of
      the wrong "no uniquely-placing molecules".
    - It warns rather than aborts, so one bad virus does not drop the run.
  (`scripts/gene_programs.py:180-187`); overlap groups from the active index's GTF.
- [ ] `PROG-15` — measure before adding antisense latency transcripts (VLT, LAT
  intron, LUNA): strandedness (`kallisto bus` runs with no strand flag), the share
  of EBV LCL reads outside annotated exons, and which markers are real index
  targets (the KSHV GTF has 26 exon rows for 96 genes).
- [x] `PROG-19` — **closed 2026-10-03.** `summarise_programs` skips the
  placeholder row when counting, and a no-model virus keeps only its "no
  programme model" caveat. Checked on the cat42b EBV LCL copy: no-model
  viruses report 0 and EBV still reports 932. Original text: (found 2026-10-03 with PROG-17) every `not_applicable`
  virus shows `n_cells_total=1` in `gene_program_summary.tsv`: the count is
  the placeholder row (`barcode=""`) that `run_one` emits. Next to
  `n_called_cells` it reads as one real cell. It should be 0. Not started.
- [ ] `PROG-18` — (found 2026-10-03 with PROG-08) validate the KSHV latent
  call on a real dataset: a public PEL or KS scRNA-seq set (BCBL-1/BC-3 PEL
  lines are latently infected). Expect mostly latent, with a lytic minority.
  Not started.
- [x] `PROG-17` — **closed 2026-10-03.** Detection writes
  `results/called_cells.tsv`; layer 2 scores only those barcodes and adds
  `n_called_cells` to the summary. Pre-PROG-17 run dirs fall back to the run's
  own `emptydrops_cells.tsv` when its method was emptyDrops, else re-call
  cells with the run config. A list that does not match the H5AD fails closed.
  The file is not a Snakemake `output:`, so resuming an old run does not
  re-run detection. Tests: `TestCalledCellsFile`.
  Original text: (found 2026-10-03 during PROG-07) layer 2
  (`rerun-programs` / `--gene-programs`) scores every cell in the multimap
  h5ad, not the called-cell set. On the cat42b EBV LCL run it scored 1,679
  cells, and only 932 of them are among the 2,763 emptyDrops-called cells
  (`per_cell_viral.tsv` `is_called_cell`). Restrict layer 2 to called cells,
  or report both sets, then rerun PROG-07. Not started (no new tasks in the
  2026-10-02 finish-pass).
- [ ] `PROG-16` — lytic acceptance test: KSHV `GSE190558` (`RUN-04`), induced vs
  uninduced; HSV-1 `SRR8315713` expected productive with no latent calls.

## WP4J — SFL tonsil TTV and HPV screen (new 2026-09-27)

Objective: answer whether torquetenovirus (TTV) and HPV are detectable in the
SFL tonsil CITE-seq pool (24 hashtagged donors, `s1`–`s24`), without publishing
an F-005-type host-homology call.

Checked before starting:
- **Only `x223` is gene expression.** `x225` is the antibody (ADT) library,
  according to cellranger `config.csv`; its R2 reads carry tag structure, not
  cDNA.
- **`x223` is 10x 5′ v3 R2-only, not 3′.** An earlier `.living` note said 3′.
  In this chemistry R2 is antisense to the transcript, and reads cluster near
  the transcription start, not the polyA site.
- The cellranger reference was GRCh38-2024-A.
- The existing index (`panel_ref_genomic`) carries only three HPV types: HPV16,
  HPV1 and HPV2. HPV18/31/33/45 are not in it. The opening note said four; the
  fourth accession, `NC_003461`, is not a papillomavirus.
- A native run is not safe yet:
  - kb's bundled 10xv3 whitelist is the 3′ list.
  - kb passes no strand flag unless given one, and `menu.py` exposes none.

- [x] `TONSIL-01` — host-subtracted screen from the existing cellranger BAM:
  1. Take the unmapped GEX reads (ADT reads excluded) and keep CB/UB.
  2. Align them with minimap2 to the panel's viral genomes plus the 16 WP4E HPV
     types.
  3. Re-check every candidate read against GRCh38.
  4. Remove UMI duplicates. Attribute cells to donors through the cellhashr
     singlets.
  5. Add a positive-control plant in 5′ geometry.

  Reads that aligned to the host never reach the unmapped set, so F-005 is
  excluded by construction.

  Decision rules, fixed before looking at results:
  - **TTV screening-positive:** ≥ 3 UMIs in ≥ 2 singlet cells of one donor. The
    call is family level unless the genus is unambiguous.
  - **HPV positive:** ≥ 3 UMIs on one type, with ≥ 1 read in URR/E6/E7.
  - **Anything else:** "not detected at this depth", never "absent".
  - The F-005 breadth gate (≤ 3.41 %) is not used to reject. 5′ capture
    concentrates true reads near the start site by design.
  - Every anellovirus result stays `screening_only` (`REF-10`).
  - Committed outputs name donors only as `s1`–`s24`.

  **Result (2026-09-27): no TTV and no HPV in any of the 24 donors — "not
  detected at this depth".** Full write-up: `.living/findings/` F-010.

  | Stage | Reads |
  |---|---|
  | Unmapped GEX reads with a valid barcode | 118,856,604 |
  | Pass the trim/complexity prefilter | 24,969,278 |
  | Any viral alignment | 34,208 |
  | Pass the viral filters | 1,511 |
  | Survive the host re-check | 23 |

  - **Anellovirus:** 1,852 raw hits, all low-complexity partial alignments
    (20–49 bp). 1,485 of them are 28-bp matches to a CAG trinucleotide repeat;
    the rest are poly(A) plus TSO-rc. None survive.
  - **HPV:** 2 raw hits, both shorter than 50 bp.
  - **The 23 survivors are not infections:**
    - HCV reads all fall on the 3′-UTR poly(U/UC) tract, at 9435–9505.
    - The macaque *Cercopithecine herpesvirus 2* reads fall on two GC-rich
      positions shared across donors.
    - HSV-1: 2 reads at one position, in one hashtag-negative barcode.
  - **HHV-6B:** 5 host-free reads at 5 genome positions and 5 barcodes. None of
    those barcodes is a singlet, so the reads cannot be attributed to a donor.
    This is below any call and is recorded, not claimed.
  - **Positive-control plant:** 1,000 reads each for HPV16 and TTV, into 1 M
    real reads. Recovery was 99.2 % for HPV16 and 99.4 % for TTV, all
    attributed to `s1`. No planted read was lost to the host re-check.
  - **Held-out anellovirus plant (sensitivity to strains not in the
    database):**
    - One genome per genus was removed from both references, and 1,000 reads
      were planted from each.
    - Recovery was 88 % (5′ window) and 73 % (uniform) for Betatorquevirus.
      Its nearest remaining genome is 94 % identical.
    - Recovery was 16 % and 35 % for Alphatorquevirus, nearest 84 %.
    - Recovery was 40 % and 8 % for Gammatorquevirus, which has no asm20
      alignment to any remaining genome.
    - **So the TTV zero rules out strains close to the 2,042 references. It
      does not rule out low-level divergent TTV.** Reads from those strains
      mostly fail to seed or fall below 85 % identity.
  - **Circular-record check:** 1 of the 1,853 dropped hits sits at a record end.
    It is G/C-run sequence, not a read spanning the origin (42 of 90 bases match
    across the junction).

  Non-obvious points:
  - **Prefilter.** Without it, 1 M unmapped reads gave 486,044 viral "hits"
    and 72 M alignments. These were TSO/poly(T) reads on the HCV poly(U) tract,
    poly(C) reads on EMCV, and poly(A) reads on A-rich anellovirus regions. The
    job ran out of memory.
  - **Missing reads.** The per-sample BAM holds only reads with a valid
    barcode. That is 1.35 B of the 1.61 B GEX reads, so the other ~16 % were
    never screened.
  - **HPV16 coordinates.** `NC_001526.4` is linearised at E1, so p97 is at
    position 7139, not 97.
  - **Limits of the result.**
    - The plant used genomes that are in the database, so this is not an LOD.
      Sensitivity to divergent anelloviruses is `ANDET-08`.
    - The EBV zero was predicted: E[EBV+ cells] ≈ 0.02–0.15.
    - The HPV prior in benign tonsil suspensions is low.
    - So none of these zeros is an informative negative (`SENS-06`).
    - BLAST spot-checks were skipped (no BLAST on the cluster). The competitive
      minimap2 re-check against GRCh38 plus the viral set stands in for them.
- [~] `TONSIL-02` — native 5′ support so `viralscan` itself can run this library
  (in progress 2026-10-06 as `DSR-12`: `--strand` exists, pilot job 25702321):
  - a `--strand` option passed to `kb count` (overlaps `PROG-15`);
  - the cellranger cell barcodes as `-w`;
  - `-x 0,0,16:0,16,28:1,0,0`.

  Measure on a 1 M-read subsample first: host mapping rate under forward,
  reverse and unstranded, then `viralscan check-whitelist`. The native counts
  must agree with `TONSIL-01`.
  - **Software audit 2026-10-12: all three items already exist, nothing new to build.**
    `--strand` (and `--strand auto`) reaches `kb count`; `-w` takes the cellranger list;
    `-x 0,0,16:0,16,28:1,0,0` is the `10xv3` geometry (16 bp barcode + 12 bp UMI).
    `tests/test_cat19_tonsil02_contracts.py` pins the geometry equality, the chemistry call from the
    pilot's list-match rate and the strand choice from the pilot rates.
  - **Pilot measurements** (jobs 25702321/25702343, 2026-10-06, first 4 M pairs; logs under
    `viralscan_work/sfl_tonsil/`): `check-whitelist` 79.9 % of 1 M R1 barcodes on the cellranger raw
    list (10xv3); `check-chemistry` 73.1 % on 100 k, calls 5' from the strand pilot, p_pseudoaligned
    forward 8.2 / reverse 51.4 / unstranded 59.2 % so `--strand reverse`. Native runs
    (`--cell-calling none`): wall 3:00 / 4:02 / 4:05, peak RSS 4.2 GB each, **0 viral molecules in all
    three arms** (1.8-1.9 M molecules in 77-87 k barcodes for reverse/unstranded), consistent with
    `TONSIL-01`'s "not detected" but not informative at 4 M of 1.6 B pairs (LOD95 0.016 per 10 k).
  - Left (cluster, `DSR-12`): the full-depth native runs and the agreement check with `TONSIL-01`.

## WP4K — Catalogue↔index reconciliation and targeted reference adds (new 2026-09-27)

Objective: stop treating "catalogued" as "detectable". F-015 measured that the
built reference omits most of the catalogue, so the first task is not importing
genomes — it is making the build honest about what it contains.

- [x] `CAT-31` — **catalogue↔index reconciliation guard (highest value).**
  F-015: 13 of 16 catalogued HPV (incl. **HPV18, HPV31**) are absent from
  `viral_genome.dedup.fa`; **all 5 Retroviridae** (both HIV-1 and HIV-2) are
  absent; plus 7/30 influenza, 4/7 coronaviruses, 3/8 polyomaviruses. A
  catalogued-but-unindexed virus is *undetectable*, so this silently caps
  sensitivity. Add a build-time assertion in
  `scripts/build_bundled_panel_ref.py` that every catalogued accession intended
  for detection appears in the emitted FASTA, and emit an explicit
  `catalogued_not_indexed.tsv` rather than failing opaquely. Decide per family
  whether a miss is intentional (e.g. segment-only, partial CDS) and record the
  reason. **Do this before any import** — otherwise `CAT-32`/`CAT-33` reproduce
  the same invisibility at larger scale.
  - **2026-09-29 (`MECH-A`):** the guard's scope is now the catalogue's `panel`
    column. The catalogue grew to 4,128 rows so the Virus Identity table can
    name max-panel genomes. Only the 2,345 `panel=shipped` rows are detection
    targets of the shipped panel, so the 1,783 `max` rows do not turn into
    reported gaps.
  - **2026-09-28: the guard shipped; the 34 decisions did not, and that is the
    point.** `reconcile_reference_panel()` in
    `src/viralscan/scripts/build_reference.py` compares the assembled panel FASTA
    against `virus_catalog.tsv` on the version-stripped, underscore-normalised
    base accession and writes `catalogued_not_indexed.tsv`
    (`accession/family/species/status/reason`) as new Step 7/8, before `kb ref`
    so a low-sensitivity panel fails in seconds rather than after ~64 GB and
    ~8 h. Misses are `intentional` only if listed in the new
    `src/viralscan/data/index_exclusions.tsv` (`accession/reason/decided_by`),
    else `unexplained`. 45 tests in `tests/test_index_reconciliation.py`.
  - **`index_exclusions.tsv` is shipped deliberately empty.** All 34 misses are
    reported `unexplained` and the guard is red. Pre-allowlisting them would
    turn a red build green without adding a genome — the exact laundering this
    row exists to stop — and the 10 GenBank-only HPV in that set are the same
    high-risk types as the 3 RefSeq ones, so "we only index RefSeq" is not a
    real policy (2,027 non-RefSeq anelloviruses are indexed). The allowlist is
    where a reviewer records "we accept this loss"; the report is the work-list.
  - **Strict is opt-in** (`--strict-reconciliation`), off by default: the
    catalogue is 11× the bundled GTF set, so a default-on gate makes the build
    unusable rather than honest, and the report is written either way. CI turns
    it on. Unexpected records (indexed but uncatalogued) and stale allowlist
    entries are reported too, so the allowlist cannot rot into a blind spot.
  - **Gotcha for the follow-up:** `build_reference.py` is a governance-pinned
    artifact, so editing it makes `test_artifact_inventory` /
    `test_claim_registry` fail until the repo's usual 3-commit re-pin chain
    (`fd8cee1` → `be1c196` → `c31742d`) is replayed against the new commit.
  - **Left open:** all 34 rows of `catalogued_not_indexed.tsv` still need a
    genome or a named decision. `CAT-32`/`CAT-33` are the first consumers.
- [x] `CAT-32` — **add EBV type 2 (`NC_009334.1`).** F-015: the index carries
  only `NC_007605.1` (type 1 / B95-8). Measured 31-mer overlap with type 2
  (`NC_009334.1` / AG876) is **79.9% shared**, leaving 20.1% type-1-unique and
  20.4% type-2-unique — a type-2 infection is currently uncallable. One
  accession, enables EBV *typing* rather than only detection. Verify against
  `CAT-27`'s EBER audit so the two EBV entries do not double-count. *(Note: an
  earlier subagent report put the type-1-unique fraction at ~77%; that is the
  shared fraction. The 20% figure is the measured one — do not re-cite 77%.)*
- [x] `CAT-33` — **HPV expansion, sourced from RefSeq `NC_` not VIRTUS2 `lcl/`.**
  F-015: VIRTUS2 carries 92 HPV *types* vs our 16 → **76 genuinely new types**.
  The apparent "94 new accessions" is partly an artifact: VIRTUS2 stores HPV as
  `gi|…|lcl|HPV##REF.1`, which will never accession-match our `NC_` records even
  for the same type — so compare by **type**, not accession. `lcl|` records are
  RefSeq *local* submissions with weaker curation; source the same 76 types from
  curated `NC_` accessions so we get GTFs. Note this is largely redundant with
  `CAT-31`: fix the 13 already-catalogued HPV first, then decide how far to go.
- [x] `CAT-34` — **record the closed negative: no `chrEBV` in our host build.**
  VirDetect warns that some hg38 builds ship a `chrEBV` contig, which would make
  EBV silently invisible to host subtraction. F-015 checked both real D-list
  inputs — `refdata-gex-GRCh38-2024-A/fasta/genome.fa` (194 contigs) and
  `Homo_sapiens.GRCh38.116.gtf.gz` — and found **zero** viral contigs. Closed:
  the hazard does not apply. Re-run the check whenever the host reference is
  swapped rather than assuming it.
**OUTCOME 2026-09-28 — the reference is built. `CAT-31`…`CAT-35` shipped; see F-016.**

Final panel: **2,343 genomes** (323 GTF-backed + 2,020 anellovirus) from
`scripts/build_bundled_panel_ref.py` on the pinned toolchain. The reconciliation
guard ran inside that build and closed **32 of the 34** catalogued gaps; the only
two survivors are the CAT-05 duplicate pair, which is the correct outcome.

What the additions actually bought, and what it cost:

- **CAT-32** — EBV type 2 is indexed. EBER1/EBER2 are now emitted too: they were
  being dropped by a `key == "CDS"` filter in `ncbi_fetch.py`, so the highest-
  abundance latent EBV transcripts were unreachable. The regenerated EBV GTF has
  96 genes (94 CDS + 2 EBER), matching VIRTUS3's independent 96-record reference.
- **CAT-33** — HPV went from 16 catalogued types to **109 human-pathogen types**
  (all 15 IARC group-1 covered). The multimapping fear was **measured and
  disproven**: 98.68 % of k-mer space is type-discriminating, HPV16/HPV18 share
  *zero* 31-mers, and L1 is *more* type-specific than non-L1. But only 68/182
  types are curated RefSeq; **114 are INSDC-only, including 8 of the 15 group-1
  types**, so the catalogue cannot treat HPV rows as curated.
- **CAT-31** — the guard works and is wired before `kb ref`, so a lossy panel
  fails in seconds instead of after ~64 GB. Its allowlist ships *empty* on
  purpose: pre-allowlisting the 34 gaps would have turned a red build green
  without a genome being added.
- **CAT-34** — closed negative recorded; the D-list genome has 194 contigs and
  zero viral, so the `chrEBV` hazard does not apply to this build.

Three latent bugs found and fixed on the way, all in F-016: **20 duplicate
`transcript_id`s across the bundled GTFs** (a live `kb ref` crash, 301
occurrences namespaced), the EBER `misc_RNA` drop, and the CAT-17 gate itself —
**dustmasker masks only 0.01 % of this panel and cannot remove the
pure-homopolymer k-mers the gate fails on**, so the builder now masks runs of ≥31
identical bases directly (panel-wide pure-homopolymer 31-mers: 85 → 0 in the
worst record, 0 panel-wide) and the unachievable `0.0` fraction default became
`0.05`, with the absolute homopolymer gate as the real control.

**New rows, from what the build exposed:**

- [~] `CAT-36` — *(2026-10-10: reverse gap is zero, every GTF seqname is catalogued with `source`/`refseq`; the
  `oncogenic_class` ask is dropped because `hpv_genes.HIGH_RISK_GENOTYPES` already carries it; the 96-row
  count is measured on a built index and is not re-checked here.)* **the catalogue is now behind the panel by 96 rows.** The build
  reports 96 genomes indexed but absent from `virus_catalog.tsv` (the INSDC HPV
  set, plus `AF157706.1` HHV-6B). They have no provenance, family or tier, and
  `CAT-31`'s reverse check will keep reporting them. Extend the catalogue, and
  propagate `source` / `refseq` / `oncogenic_class` for the HPV rows.
- [x] `CAT-37` — **inherit a format-version stamp for the GTF cache.**
  - **Done 2026-10-04 (c3f930c).** `ncbi_fetch.GTF_FORMAT_VERSION` plus a `.fmt`
    stamp. A stale or missing stamp regenerates from the retained `.gb`,
    otherwise it refetches. Existing caches regenerate once. Bump the constant
    whenever GTF emission changes.
  `ncbi_fetch._cache_valid` only checks a file against its own sidecar, so the
  2,249 already-generated GTFs are silently reused and the EBER fix does not
  reach them until they are deleted. A cache-key change is the durable fix.
- [x] `CAT-38` — **done 2026-10-10 (single-index route; `scripts/fix_thin_gtfs.py`).** The genome-spanning `misc_RNA` of HAV and `prim_transcript`/`misc_RNA` of HTLV-2 are removed and HTLV-2's gene row is re-spanned to its real CDS (2239-5187); local GTFs, manifest re-pinned. Original: **two whole-genome pseudo-transcripts entered the panel with
  the EBER fix.** HAV `NC_001489.1` (misc_RNA spans 100 % of the genome) and
  HTLV-2 `NC_001488.1` (94 % of 8,952 nt) are exactly the shape that collapsed
  99.8 % of anellovirus UMI into one bucket. Belongs to the `CAT-18` two-tier
  split (`CAT-18` was closed as not adopted); do not ship these two in a single-index panel without deciding.
- [x] `CAT-39` — **done 2026-10-10 (modelled; `scripts/model_hpv_cds_gtfs.py`, QC `analysis/panel_expansion/hpv_modelled_cds_qc.tsv`).** HPV-8 gets E6/E7/E1/E2/L2/L1 from HPV-5 (identity 70-87 %, lengths 97-105 % of the sibling); HPV-71 keeps its annotated E1 and gets E6/E7/E2/L2/L1 from HPV-16 (identity 36-59 %, lengths 91-105 %, **family-level homology, so treat E6 and E7 as weakly supported**). Spliced E1^E4 and variable E5 are not modelled. Original: **`M12737` (HPV-8) and `NC_039089` (HPV-71) are thin.** Neither
  has usable CDS upstream, so they fall back to whole-genome pseudo-transcripts
  (`CAT-38`). Re-cut or exclude deliberately rather than by accident.
- [~] `CAT-40` — **run `kallisto inspect` on the new index** (max EC size,
  discarded EC count) before trusting it. The k-mer-overlap metrics in F-016 are
  computed outside kallisto and cannot see the pseudo-inverse threshold effect,
  which is `CAT-19`'s actual failure mode.
  - **Software done 2026-10-12** (`src/viralscan/reference_gates.py`, `python -m
    viralscan.reference_gates inspect`): `parse_inspect` (recorded text or live index),
    `index_inspect.json` written on every bundled-panel build
    (`build_bundled_panel_ref.py`, new `--ec-baseline`; exit non-zero on regression) and, as a
    one-line non-fatal hook in `build_reference.py:_finish_reference_build`, on every
    `viralscan build-ref` index build (record only; separate commit, `build_reference.py` is
    pinned and was re-pinned). Regression =
    either metric above the baseline record; no absolute threshold was invented (`--max-ec-size`
    / `--max-discarded-ecs` are caller-supplied budgets). Tests: recorded text plus a tiny
    kallisto index in `tmp_path` (integration, skipped without kallisto).
  - **Measured 2026-10-12 on the existing indexes** (`analysis/panel_expansion/index_inspect_*.json`):
    `viral_ref_final` (2,343 genomes, 471,948 targets, 84,861,133 k-mers) and `viral_ref_cat42d`
    both give **max EC size 9,921, ECs discarded 0**. Note: the conda `kallisto` on PATH aborts
    (`bad_array_new_length`) on a `kb ref` index, so inspect runs with kb's own binary
    (`validation.tool_path`, REL-16).
  - Left: the user's next full build exercises the hook; whether 9,921 is acceptable (EC size
    is driven by the near-identical anellovirus genomes) is a user decision, so no pass/fail
    claim is made. The baseline file is the `viral_ref_final` record.

- [x] `CAT-35` — **VIRTUS2/VirDetect/VIRTUS3 method extraction (research, done
  2026-09-27).** Conclusion recorded in F-015: take the **method**, not the
  genomes — host subtraction before viral quantification, per-strand counting,
  explicit multimap handling. Reject VIRTUS2's 37 influenza (legacy lab strains:
  PR8, H9N2 HK/97, H5N1 goose 1996, H3N2 NY/2004, B/Lee/1940) as no better than
  what we hold, and its 27 anelloviruses as already covered with better
  provenance by the upstream 2,023 reps. VIRTUS2 parsing hazard recorded: the
  list uses space-separated accessions (`NC 000883.2`) while its FASTA uses
  underscores, which silently drops 96 of 762 records under a naive regex.

- [x] `TAXO-01` — **NCBI lineage per virus and an aggregate-up count (user, 2026-10-16; branch `codex/virus-taxonomy-rollup`).** The
  catalogue `species` label is not a lineage: EBV (taxid 10376) reads "Human gammaherpesvirus 4",
  type 2 (12509) reads "Human herpesvirus 4 type 2", and families mix Herpesviridae and
  Orthoherpesviridae, so related viruses cannot be summed. Ship `data/virus_taxonomy.tsv` (NCBI
  Taxonomy efetch per catalogue taxid: species taxid and name, genus, subfamily, family), carry it into
  `virus_identity.tsv` as additive `*_ncbi` columns, and write `viral_summary_by_species.tsv` and
  `viral_summary_by_genus.tsv` from the multimap molecule pass: molecules unique to the group,
  every molecule compatible with it, the package total, and their ratio. Flag a virus whose own unique
  fraction is < 0.05 in `viral_summary.tsv` and point to its species row. No count changes.
  Done: 710 catalogue taxids resolved (0 missing; `extras/build_virus_taxonomy.py --check` reports drift).
  Decisions: the lineage is a frozen packaged table, not a build-time lookup, so every run is offline and
  reproducible (uncatalogued genomes have no taxid and stay their own group); NCBI "… sp." bins
  (e.g. Gammapapillomavirus sp.) carry no species; anelloviruses keep their ANDET-05 genus key at the species
  level. Group molecules are method-independent and stored in `adata_multimap.h5ad` obsm, so
  `rerun-multimap` stays valid; existing run directories get the new tables only after a multimap rerun.
  Checked on EBV GSE158275 GSM4796273 (kept combined_off files, replay): EBV 15,011 unique / 18,417
  compatible (0.815), EBV type 2 177 / 3,583 (0.049, flagged), species *L. humangamma4* 18,594 / 18,594
  (1.0); package totals equal the run's (17,193.816667 and 1,400.183333). `viral_molecules_upper_bound`
  doc corrected: it splits virus-virus molecules equally, so it is not a sibling bound.
  Tests: `tests/test_virus_rollup.py`. Governance pin to re-pin (not edited here): inventory row 17
  `output-reference-doc` (docs/output_reference.md); `tests/test_artifact_inventory.py` fails on the stale
  sha256 until the re-pin.


## WP4L — Audit 2026-10 findings (new 2026-10-11)

Three read-only agent audits (reference/GTF, detection/statistics, CLI/fetch) plus a Biomni `max` science review.
Findings: `analysis/audit_2026-10/findings.tsv` (76 rows: 0 blockers, 13 majors + 5 Biomni-raised; each row has a
`biomni_verdict`). Repro scripts under `analysis/audit_2026-10/repro_{A,B,C}/`. All findings reproduced on synthetic or bundled data,
none on a real run; Biomni did not run code. `B-22` is absent from the numbering (agent B reported 35, wrote 34).
Biomni's top five, in order: AUD-C reference coherence, AUD-A attribution, AUD-B detection/denominators, AUD-D labels, AUD-E statistics.

- [x] `AUD-01` — `A-01` overlap groups: replace transitive closure by k-mer/read-compatibility grouping (31-mer sharing from the index,
  strand-aware; sensitivity reported). Supersedes the "min overlap" option under Next action (1): a 91-bp overlap cut is not an ambiguity rule.
  Done 2026-10-11 (`codex/aud-p2`): `build_overlap_groups` now groups by canonical 31-mer sharing (threshold 0.5 of the smaller gene; marker groups are flat across 0.25/0.5/0.9, see `analysis/audit_2026-10/aud01_sensitivity.tsv`) as greedy maximal cliques (maximal, not minimum), keyed by seqname; the generator fails on a cross-programme pair sharing >= 50% in different groups. **User-approved 2026-10-11:** EBNA-1/2/LP are accepted as separate groups (read ambiguity only; EBV latent groups 3 -> 4) and HSV-2 LAT s01 is merged into s02 (catalogue 75 -> 74 rows, as HSV-1). Strand note `analysis/audit_2026-10/M01_strand_note.md`: the pinned kallisto default is unstranded (unverified against kallisto source), so grouping ignores strand. Known over-count: breadth can over-count linked genes in different cliques (EBNA-1 ~ EBNA-LP 57%); default `latent_min_breadth=1` unaffected. Docs to update: output_reference.md :436-441 and :402, faq.md :314-317, panel_integrity.md :106/:120-122, `analysis/panel_expansion/ref13_normaliser_qc.tsv` rows 527-528 (`wrap_native`).
- [x] `AUD-02` — `A-03` `gene_programs.run_one` drops cells with unique breadth 0 but selected breadth > 0 (undercounts `n_cells_*_selected_layer`).
  Fixed 2026-10-11 (`codex/aud-p1`): `run_one` skips a cell only when all four breadths (unique and selected) are 0; `n_cells_*_selected_layer` was understated. A-10 (resolve_markers fallback) and A-14 (assert -> ValueError; the never-filled `n_cells` column of `gene_program_cells.tsv` dropped) in the same bundle. Docs: `docs/output_reference.md` ~480 should say `gene_program_cells.tsv` now includes cells with unique breadth 0 and selected breadth > 0.
- [x] `AUD-03` — `A-04` curated marker notes contradict RefSeq product text (8 EBV/HHV-7 genes); regenerate from product or cite. Keep programme/kinetics separate from function.
  Done 2026-10-11 (`codex/aud-p2`): `parse_gtf` reads `product` from CDS rows; EBV/HHV-7/HSV/VZV notes corrected from RefSeq product text (programme and kinetic class unchanged); U90 cites PMID 12388818; LMP-1 and EBNA-1 wording fixed; PMIDs 20860481 and 36733461 marked unverified.
- [x] `AUD-04` — `B-01` host-response treats each viral gene ID as a virus; pool by virus.
  Fixed 2026-10-11 (`codex/aud-s1`): `hostresponse._pool_viral_genes` sums gene columns per virus (identity table when present, else the name rule); per-virus outputs are now named by virus display name. B-12/13/14/17 fixed in the same bundle. **Follow-up:** `viralscan hostresponse` (menu) still passes `log/analysis.txt`, so it uses the legacy name rule; prefer `identity_path(output)` there.
- [x] `AUD-05` — `B-02` `umap` crashes without `plots/` (negative sample or `--no-visual`).
  Fixed 2026-10-11 (`codex/aud-s1`): `umap()` creates `plots/`; B-24 all-cells-filtered now fails loudly.
- [x] `AUD-06` — `B-03`/`B-23` "infected" = count > 0 on fractional EM weights; one declared evidence rule across detection/enrichment/umap/hostresponse.
  **Decision (user, 2026-10-11):** a cell is infected for a virus when its summed allocated molecules across that virus's genes are >= `--infected-min-molecules` (default 1.0, recorded in the manifest); fractional-only cells are reported as sub-threshold, not dropped; one rule reused by detection, enrichment, umap and hostresponse.
  Done 2026-10-11 (`codex/aud-d1`): `utils.virus_molecules_per_cell`/`infected_masks`; infected = summed weight >= `--infected-min-molecules` (default 1.0, relative tolerance 1e-6 for float32 EM sums); sub-threshold cells are counted, not dropped; used by detection, enrichment and umap (B-23). B-10 (`-1` suffix) fixed. `_headline_totals` 'Cells with viral reads' still counts weight > 0. **Resume:** the new options are in `omit_when_unset`; a plain run's fingerprint is unchanged, but re-running a detection step of an old run applies the new defaults.
- [x] `AUD-07` — `B-04`/`B-16` detection gate on all barcodes; per-cell label defaults to 1 UMI via CLI. Needs background-aware design (`M-03`), not a bare N-cells cutoff.
  **Decision (user, 2026-10-11):** detected = at least N distinct called cells (default 1) and at least `detection_threshold` UMIs in called cells; empty-droplet UMIs are reported separately as background; `summary.txt` says "candidate" and flags `infected_called == 0`. B-16: new `--hostresponse-min-umi` (default 10) for the per-cell label; the sample gate stays `detection_threshold`; fingerprint-safe via `omit_when_unset`.
  Done 2026-10-11 (`codex/aud-d1`): `detected` = `infected_called >= --detection-min-cells` (default 1) and called-cell molecules >= `detection_threshold`; `background_umis` = total minus called-cell molecules; every all-barcode candidate stays in `viral_summary.tsv` with `detected` False/True; `summary.txt` says candidate and flags no-called-cell viruses (with cell calling none the gate is unchanged for integer counts). B-16: `--hostresponse-min-umi` (default 10) for the per-cell label; the `hostresponse` subcommand's `--detection-threshold` is a deprecated alias. **Docs to add** (docs/output_reference.md, pinned): `viral_summary.tsv` columns `infected_subthreshold_cells`, `infected_subthreshold_called`, `background_umis`, `detected`; infected/pct rows now use the molecule rule; `per_cell_viral.tsv` infected cells only; summary.txt wording. A passed gate is candidate evidence only; ambient RNA is not excluded.
- [x] `AUD-08` — `C-01` `--gene-programs` fails at DAG build (`results/viral_summary.tsv` has no producer); add a dry-run DAG test.
  Fixed 2026-10-11 (`codex/aud-p1`): phantom `summary` input removed from `rule gene_programs`; `TestGeneProgramsDag` (integration) plans the DAG with `gene_programs: true`.
- [x] `AUD-09` — `C-02`/`C-03` NCBI API key and email stored in `run_manifest.json` and leaked in `_efetch` error text. Fixed 2026-10-11 on
  `codex/aud-09-secrets`: both options are excluded from the manifest and fingerprint (a null is kept so runs made without them still
  resume), and `_efetch` reports the exception type, not its URL-bearing text. Two regression tests fail on the old code. **Caveat:** a run
  made *with* `--ncbi-email`/`--ncbi-api-key` before this change will not match on `--resume`; any `run_manifest.json` already written
  with a key holds it in plain text: delete or rotate.
- [x] `AUD-10` — `C-04` `validate-run` false `missing_input` from another cwd (relative paths).
  Fixed 2026-10-11 (`codex/aud-r1`): the manifest gains a resolved `input_paths` map outside the fingerprint; `validate-run` works from any cwd. Same bundle: C-07 (`--yes`/EOFError on the `--reference` prompt), C-08 (clean snakemake failure), C-09 (pre-run path checks), C-10 (numeric range checks, also in `RunConfig`), C-12 safe half (resume error names the differing options; `cores` deliberately stays in the fingerprint), C-13 (host_index directory, host_h5ad, cell_types, host_fasta, viral_fasta and fetched -acc sequences are fingerprinted). **Resume caveat:** runs using `--host-filter`/`--host-index`, `--host-h5ad`, `--cell-types`, `--host-fasta` or `--viral-fasta` started on older code will refuse `--resume` (use `--overwrite`); plain runs keep an identical fingerprint (test pinned to the aca8998 value). Docs to update: cli_reference.md :135-140, :644 (`--yes`, `--resume`); output_reference.md :145-147 (pinned) for `input_paths`/`fetched_reference`.
- [x] `AUD-11` — `C-05`/`C-06`/`C-13`/`C-14`/`C-17` reference coherence: `kb ref` without `--overwrite` in `build-ref`, Ensembl cache not release-keyed, fingerprint gaps, truncated FASTA cached, host genome treated as viral.
  Fixed 2026-10-11 (`codex/aud-r1` C-13, `codex/aud-r2` the rest): C-05 `_run_kb_ref` always passes `--overwrite`; C-06 host cDNA cached under `<key>/release-N/`, new `build-ref --ensembl-release N`, no silent fallback to release 116 (an unreachable release lookup errors; existing flat caches are re-downloaded once); C-13 (see AUD-10); C-14 FASTA length/version checked against GenBank before caching; C-17 `--reference` on a host-like GTF warns (does not refuse). Same bundle: C-15 version-aware de-dup, newline per chunk, NCBI email only on a cache miss, OSError wrapped; C-16 unwritable cache dir gives `ViralScanDataError`; A-17 one short t2g row no longer flips the file to legacy; A-13 pathlib. **Pins:** `build_reference.py` changed, so `test_artifact_inventory`/`test_claim_registry` fail on its sha256 until the end-of-wave re-pin. Docs: README.md:152 and docs/quickstart.md:82-83 (email only needed on a cache miss; README is pinned).
- [~] `AUD-12` — Biomni-raised, not audited (`M-01` strand defaults across count/extract, `M-02` per-marker polyA/capturability, `M-03` ambient/swapping controls, `M-04` real-data sensitivity/specificity, `M-05` competitive host reference).
  M-03 design note written 2026-10-12 (`codex/aud-12b`, design only, no code): `docs/plans/ambient-and-barcode-swapping-controls.md` (soup bound from `background_umis`-style ambient counts, cross-run barcode/key swap tests, controls, flags, calibration, 7 user decisions). Generator gap found: `generate_truth_panel.py` `hop()` does not preserve the barcode.
  M-04 design note written 2026-10-12 (`codex/aud-12b`, design only, no code): `docs/plans/real-data-sensitivity-specificity.md` (truth tiers T1-T5, arms P/D/N/X on `generate_truth_panel.py`/`plant_anello_10x.py` and `VAL-01..10`, `SENS-04/06/07`, `DSR-06/07`; metrics reuse protocol M2-M9/E8; claim limits; library-level bounds; 8 user decisions). Row stays `[~]`: nothing implemented.
  **M-01 done 2026-10-12** (`codex/aud-12a`, write-up `analysis/audit_2026-10/M01_strand_audit.md`): the resolved strand is written to the per-sample `config.yaml` and the evidence replay re-reads it, so count and replay agree on strand (no code change needed there); ViralScan never calls `kb extract`. One real mismatch fixed: replay passed the raw `10xv4` to `kallisto bus` (kb count rewrites it to `10XV3`; kallisto 0.52.0 rejects `10XV4`). Left as design/remaining: unrecorded replay strand in `evidence_manifest.json`, `-x ...%STRAND` suffix, the pilot reading raw FASTQs and a raw `-w`. M-02..M-05 still open, so the row stays `[~]`.
  M-01 remaining item 4 done 2026-10-12 (`codex/aud13-andet-followups`): `chemistry_check.sanity_gate` no longer blames a "forward-strand default" for a low pseudoalignment rate; it takes the run's resolved strand (`sanity_report`/`sanity_gate_from_dir`, wired from `analysis.py`), names `--strand forward` as the F-020 cause only when it was used, and otherwise says no flag is unstranded under kallisto 0.52.0. Tests in `tests/test_chemistry_check.py`.
  M-01 remaining item 1 done 2026-10-12 (same branch): `evidence_manifest.json` gains an optional `replay` block (`strand`, `kallisto_strand_flag`, `technology`, `kallisto_technology`) on both write paths (empty extraction and final); schema `schemas/v3` and its packaged copy updated together (`schema_version` stays 3.0.0, the property is optional and old manifests stay valid); `output_reference.md` one-line hunk on the `evidence_manifest.json` row (pinned: re-pin pending). Still open: `-x ...%STRAND` suffix, strand-auto pilot on raw FASTQs/`-w` (design choices).
  **User decisions 2026-10-12 (grill, all taken as recommended):** *M-03:* flags are descriptive only (`detected` unchanged; gating needs a protocol amendment); first-release flag logic is the calibration-free `rho = 1` bound with a Poisson upper confidence bound, `rho` estimate deferred; the cross-run swap step is a `scripts/` tool reading a user-written pool manifest (no new subcommand); two-step runs carry per-barcode host totals by extending `host_cells`; doublets/uptake out of scope with a claim limit. *M-04:* add a CB-preserving hop option to `generate_truth_panel.py` (default output byte-identical, then re-pin) and the `D1` real-PBMC mode (planted plus unplanted twin), as separate commits; hopped/ambient molecules in a negative are reported three ways with the unconditional count as headline; truth tiers T1-T5 reported separately, never pooled (dataset-to-tier assignment waits for the G3 round); specificity headline is the any-virus sample-level false call with per-virus rates alongside; programme LOD is direction only for the first release. Deferred: tier assignment, acceptance thresholds, which context-only positives become scored rows, wet controls.
  - [x] `M03-IMPL` — descriptive ambient/swap columns and flag in the detection outputs (`rho = 1` bound + Poisson UCB), `host_cells` extended with per-barcode host totals for two-step runs, `scripts/` cross-run swap tool with pool manifest; schema/docs/tests; no change to `detected`.
    Done 2026-10-12 (`worktree-agent-acf799a50ff8fc352`, from `codex/grill-decisions-2026-10-12`). Nine `viral_summary.tsv` columns (`ambient_viral_fraction`, `_ub95`, `n_ambient_barcodes_with_virus`, `soup_expected_cells_max`, `_ge2_max`, `soup_tail_p_cells`, `viral_depth_spearman`, `ambient_assessment`, `ambient_basis`; `detection.soup_screen`, `compute_stats`); `cellcalling.HostCellSets.host_totals` (the `read_host_matrix_totals` result, zero totals omitted) so two-step runs need no `GeneFull/raw` re-read; `scripts/swap_check.py` (tier 1 barcode overlap, hypergeometric, pool manifest `library_id/run_dir/pool_id/whitelist`, schema + example in `docs/output_reference.md`); the `contamination` evidence row carries the screen. `detected`, tiers and denominators unchanged (test). Evidence: full suite 2026-10-12 3,401 passed, 12 skipped, 1 failed (`test_artifact_inventory`, the `docs/output_reference.md` pin, re-pinned in the next commit); `ruff check`, `ruff format --check`, `mypy src/viralscan`, `check_data_governance.py`, `validate_v3_protocol.py` clean.
    Decisions and limits: (1) `ambient_assessment` is `soup_rejected_rho1` when `soup_tail_p_cells` < 0.05, else `soup_sufficient` (also when no called cell is infected); only `E1` drives it, `mu` and top-cell share are not computed. (2) Ambient barcode = not called with total (host + viral molecules) <= `--emptydrops-lower` (default 100, inclusive); two-step totals add STARsolo host molecules to the viral-matrix row total. (3) Blank (not `not_applicable`) where undefined, per `pct_infected_comparable`. (4) `not_assessed` for `--cell-calling none`/`knee`, virus-only index without host totals (a two-step virus-only run IS assessed), no ambient barcode, no called cell, two-step with no `GeneFull/raw`, alignment-only rows. (5) The `contamination` flag status stays `not_assessed` (swapping is never measured in one run), the screen goes in `metric`/`basis`; nothing downstream branches on it. (6) Tier 2 (shared CB-UMI key over BUS) not built, per design section 10. Planted-data validation of the screen waits on `M04-GEN` (CB-preserving hop); calibration of any gate waits on the G3 freeze. Doublets/uptake: claim-limit text only (`docs/output_reference.md`).
  - [x] `M04-GEN` — CB-preserving hop option in `generate_truth_panel.py` (default output byte-identical, pinned-file re-pin) and `D1` real-PBMC mode in `plant_anello_10x.py`/generator, as two commits.
    Commit 1 done 2026-10-12 (`worktree-agent-a495d8c6d6720fda3`): `--hop-preserve-cb` (BooleanOptionalAction, default off); hopped molecules keep the donor `(CB, UMI)`; `truth_manifest.tsv` gains `hop_donor_sample_id`, `hop_donor_molecule_id`, `hop_cb_preserved`, `hop_landing` (cell/empty/other) only when the option is on, and `run_manifest.json` gains `hop_mode`. Default output unchanged (AUD-18 hash test passes). Tests `tests/test_generate_truth_panel_m04.py`. Design: generator design doc section 8b. Decision: no knob steers hops into cells (landing follows the recipient's barcode space; `--onlist` makes coincidence real). Re-pin: `generate_truth_panel.py` is in neither `analysis/v3_artifact_inventory.tsv` nor `claims/registry.json` on this base, so none was needed.
    Commit 2 done 2026-10-12: `D1` mode (`--background-r1/--background-r2`, `--background-chemistry`, `--n-background`). A real 10x FASTQ pair is the background; each planted `synthetic_factorial` sample is the library's reads verbatim (qualities kept, names `p<i>`) plus planted molecules, with an unplanted twin made of the same reads (`twin_sample_id`, shared `background_library_id`). Background reads carry truth `unknown` (a new origin value in generator output); cells are the top barcodes by read count (stand-in for `cell_universe`); `plant_anello_10x.iter_background` is the shared streamer (`copy_background` output unchanged). Design: generator design doc section 8c. Resolved readings (conservative): `low_complexity` forced to 0 and no artefact reads (twin differs by planted molecules only); hops off (one library), so `realised.hopped` is 0 while `requested_hopped` is above 0 in a D1 manifest; always training (one library cannot fill a holdout seat); R1 must be exactly barcode+UMI long; one chemistry, `synthetic_factorial` only; synthetic host templates still written. Not done: streaming for a production-size library, multi-library holdout allocation, the scorer (`M04-SCORE`). Evidence 2026-10-12 on the rebased tree: full suite 3393 passed, 12 skipped, 133 deselected; `tests/test_generate_truth_panel_m04.py` 17 tests; AUD-18 default hash test passes; `ruff check`, `ruff format --check`, `mypy src/viralscan`, `check_data_governance.py`, `validate_v3_protocol.py` clean.
  - [x] `M04-SCORE` — three-number scoring of hopped/ambient molecules in negatives, T1-T5 tier vocabulary in `v3_score.py` and reports (after `M04-GEN`).
    Done 2026-10-12 (`codex/m04-score`, code only, no outcome run): `v3_score.score_negative_libraries` (unconditional any-virus false call = HEADLINE with per-virus rates; libraries excluded with hop/ambient provenance; contamination-attributable share beside the residual; `format_negative_report` prints all three), `TIERS`/`check_tier` (T1-T5 validated, mixed tiers refused, `score_negative_libraries_by_tier` has no pooled entry), `realised_shortfalls` (D1 exempt), `truth_positive_molecules` and `read_truth_manifest` accept `origin=unknown` and refuse it as planted, `hop_*` columns validated. Tests `tests/test_v3_scorer_m04.py` (hop negative, ambient-only, tier pooling, bad tier, D1 exemption, unknown origin, real generator hop/D1 manifests, pre-M04 regression: reader digest on the default panel equals the pre-change reader's). Docs: last section of `docs/plans/real-data-sensitivity-specificity.md`, CHANGELOG. No acceptance threshold; no tier assigned to any dataset; no programme LOD. Resolved readings (most conservative): (ii) excludes a library only if it holds documented hopped/ambient viral molecules and its call, if any, is fully explained by them (a call with one other or no supporting molecule stays a false call); T5 is refused for scoring; the caller supplies each call's supporting `(cb, umi)` molecules because the scorer does not re-run the gate; a negative holding cell-origin planted molecules is refused. Not done: wiring into a run-level report or the VAL-08 adapter (no real negative run exists), reading calls out of ViralScan outputs. Evidence 2026-10-12: full suite exit 0 (3428 collected = 3416 passed + 12 skipped, 133 deselected; baseline 3393 + 12 skipped, plus 23 new in `tests/test_v3_scorer_m04.py`); `ruff check .`, `ruff format --check .`, `mypy src/viralscan`, `check_data_governance.py`, `validate_v3_protocol.py` clean. `scripts/v3_score.py` is in neither `analysis/v3_artifact_inventory.tsv` nor `claims/registry.json`, so no re-pin.
- [~] `AUD-13` — all 94 original and audit-4 findings triaged on 2026-10-11; software fixes and explicit deferred scientific/policy dispositions are being integrated. The parallel close-out also covers A-12 annotation prose, A-15 partition staleness, A-18 catalogue-drift provenance, B-18/B-19 host-response family/population accounting, and the stale legacy execution-freeze regression (GOV-06 tooling only; frozen diagnostic receipts are unchanged).
  **Agent close-out 2026-10-12 (`codex/aud13-andet-followups`).** `findings.tsv` has no disposition column, so "open software fix" was taken as: class not `science`, severity `minor`, not in the 2026-10-11 progress list, and the repro/code still shows the defect. A-12 done: `gtf_normalise` docstring now states the audit's measured figures (2,522 CDS-only genes without an exon row; 101,346 bp over 105 spliced transcripts of 103 genes), labelled as local-corpus numbers (the GTFs are gitignored); `tests/test_gtf_normalise.py` pins them and a `research`-marked test recomputes them from the 393-GTF corpus. The CLAUDE.md figure (393) was already corrected.
  A-18 done: `virus_identity.write_identity_table` records `catalogue_drift` (`match`/`drift`/`not_pinned`/`unavailable`, both sha256) in the run manifest's per-sample block (outside the fingerprint, like `chemistry_sanity`); no strict flag, drift stays advisory (policy); `output_reference.md` gets one separate paragraph (pinned: re-pin pending).
  Left open on purpose (2026-10-12), all remaining non-science minors: B-10, B-11 (`detection.py`/`cellcalling.py`, M03 owns them; B-10 already has a suffix test in `test_cellcalling.py`), B-18, B-20 (host-response code, out of this bundle), C-19 (CI/release, `REL-12` user-owned), C-12 second half (`cores` stays in the fingerprint, user decision recorded under `AUD-10`), B-35 and C-18 (test-gap rows: most named sub-gaps now have a test (cpm ties, NaN BH, gini, umap `plots/`, `input_paths` from another cwd, gene-programme DAG, kb rerun, credential leak, release-keyed cache); a grep found none for `_efetch` 429/5xx retry, and the rest were not re-audited one by one; test-only rows, so no failing-before test is possible). Already landed, not re-done: A-15 (`tests/test_gene_programs_staleness.py`), B-23 (`AUD-06`), B-30 (hash suffix in `reads.py`), C-11 (staged `kb-python` move, `test_kb_count_rerun_preserves_metadata_and_previous_counts`). Science/policy rows (A-02, A-04, A-05, A-16, B-15, B-16, B-19, B-25, B-26, B-29, B-31..B-33, C-13/C-14/C-17 science halves, M-02..M-05) and D-* (`findings_D_remaining.tsv`, closed under `AUD-19`/`AUD-20`) untouched.
  2026-10-12 (`codex/b19-c12`): C-12 fixed (user decision): `cores` is no longer hashed into the run fingerprint (still recorded in `options`); a pre-C-12 manifest whose stored hash reproduces with `cores` included still resumes under any `--cores`; the mismatch error lists differing option names (never `--cores`). The plain-run fingerprint pin in `test_run_safety.py` changed accordingly.
  2026-10-12 (`codex/b19-c12`): B-19 fixed (user decision): `run_hostresponse` aligns on the host h5ad cells; host cells missing from the viral matrix are zero-filled negatives, viral-only barcodes are dropped, and `barcode_alignment` in `hostresponse_manifest.json` records counts and fractions of both (`policy` is now `host_cell_universe_zero_fill`; `n_host_excluded`/`n_viral_excluded` are gone). An empty host h5ad or zero barcode overlap now raises `ValueError` (the old zero-overlap early return would have zero-filled everything). Supersedes the 'user-owned, open' B-19 note under Tier 6.
- [x] `AUD-18` — generator realism (user-approved 2026-10-11): `scripts/generate_truth_panel.py` `--viral-fasta` + gene-span TSV hook so the panel maps to the cat42b reference, and `--onlist` for real 10x barcodes; unblocks the `kb count`/ViralScan arm of `VAL-08`.
  Done 2026-10-11 (`codex/aud-v2`, code only): new `--viral-fasta` + `--gene-spans` (+ optional `--sibling-pairs`) and `--onlist` (`.gz` ok); input hashes under `inputs` in `run_manifest.json` (none when unused). Strand-aware reads (a `-` gene is read from its reverse complement; `truth_class` ambiguous iff the footprint fits another same-strand gene); viruses outside any sibling pair are absent controls; real genomes are not mutated (`pair_identity` empty). Default output pinned byte-identical to 7b0c614 (hash test). Tests run offline on tiny inputs written to `tmp_path` (`tests/_truth_panel_real.py`). Not done: running it against the cat42b FASTA/GTF (needs the reference files, `VAL-08`). Docs: generator design doc section 8a, CHANGELOG. `gen_cli_reference.py` is unaffected (script is not a `viralscan` subcommand).
- [x] `AUD-21` — follow-ups from `AUD-07`: `gene_programs.py` and the `scripts/dsr02*` / anello-acceptance reports still treat every `viral_summary.tsv` row as detected; they should filter on `detected` (only `auto_evidence` does today). `_headline_totals` still counts weight > 0.
  Programme side done 2026-10-11 (see AUD-19). Reports done 2026-10-11: `scripts/dsr02_enumerate_calls.py` and `scripts/anello_acceptance_report.py` skip rows with `detected` = false (tests added). `detection._headline_totals` done (ab9d9da, found already landed 2026-10-12): it takes `{virus: gene IDs}` of the `detected` viruses only and counts cells with the one infected rule (`infected_masks`, `--infected-min-molecules`), union over viruses, plus `cells_with_virus_called`; tests `test_candidate_only_virus_is_not_headlined`, `test_cells_follow_the_infected_rule_not_any_weight`, `test_cells_are_split_by_called_mask_and_counted_once_across_viruses`, `test_subthreshold_weights_are_not_pooled_across_viruses` in `tests/test_detection.py`. `AUD-21` closed.
- [x] `AUD-14` — finish the audit (user-approved 2026-10-11): `evidence_run.py`, `reference_strategy.py`, `constants.py`, the PROG-12/14 gating in `gene_programs.py`, the anellovirus builder. Findings only, fixes for clear bugs.
  **Audit 4 done 2026-10-11** (`analysis/audit_2026-10/findings_D_remaining.tsv`, D-01..D-18: 5 major, 13 minor, no blockers). Majors: D-01 `reactivating` reachable for VZV/HHV-6A/6B/HHV-7 (all `latency_not_observable`; only HCMV is allow-listed), D-02 the PROG-14 reduced marker set can silently produce `reactivating`, D-06 `reference_strategy.parse_viralscan_metrics` needs the legacy `viral_umi` column (every ViralScan benchmark row fails on v3 outputs), D-07 `parsed_benchmark_rows` silently shrinks the barcode universe, D-17 missing tests for these. Fix bundles: `AUD-19` (programmes: D-01/02/03/18), `AUD-20` (reference_strategy: D-06/07/08/09/10/11/12 and evidence_run D-04/05/13). Not audited: `evidence.py` internals, `anello_align.py` pileup heuristics, STAR multimap behaviour (D-15).
- [x] `AUD-19` — programmes fixes from audit 4 (D-01, D-02, D-03, D-18): gate `reactivating` by the catalogue's own `latency_observable` flag instead of an HCMV allow-list; a reduced marker set must not silently yield `reactivating`; `n_cells_reactivating` reports `not_applicable` where unreachable; `run_one`/`summarise_programs` read one catalogue row per virus and `validate_catalogue` checks they agree.
  Done 2026-10-11 (`codex/aud-g`): D-01 `reactivating` now also requires `latency_observable_in_rna` (after this only EBV can reach it; VZV, HHV-6A/6B, HHV-7 and HCMV give `indeterminate`; **resolves the PROG-12 'gate beyond HCMV' open question**); D-02 strict rule: any unresolved productive catalogue marker degrades both absence-based states (`reactivating`, `latent`) to `indeterminate`, recorded in `marker_resolution.absence_assessable` and the caveat; D-03 `n_cells_reactivating` = `not_applicable` where unreachable; D-18 one `catalogue_virus_facts()` reader and `validate_catalogue` fails on disagreeing rows. AUD-21 programme side: viruses with `detected` False get no programme row. Docs (pinned `output_reference.md`): :463 reactivating row, :478 `n_cells_reactivating`, :472 and ~:462 latent/indeterminate notes, Layer-1 table note on candidate-only viruses. Not done: a column naming which classes were assessable.
- [x] `AUD-20` — reference_strategy / evidence_run / constants fixes from audit 4 (D-04..D-13).
  **User-approved for autonomous work 2026-10-11**, together with the EEEV output-name collision fix (`hostresponse._safe_name`: 'Eastern equine encephalitis virus' vs 'Eastern_equine_encephalitis_virus' sanitise to one name). Push and PR still require asking; merged agent worktrees and branches may be removed at the end (never the locked wave4 worktree).
  Done 2026-10-11 (`codex/aud-h`): D-04 invalid flag combinations die before the replay; D-05 manifest hashing streams via `run_safety.sha256_file`; D-06 `parse_viralscan_metrics` reads `viral_molecules_total_est` (`viral_umi` kept as legacy alias) and anchors on the called cells; D-07 a failed combined row marks the dataset's rows `incomplete` and the combined-vs-two-step delta is `nan` unless both sides completed; D-08 HSV-1/EBV regexes human-anchored; D-09 presence scan streams the whole file and provenance counts must add up; D-10 slurm array `chdir` resolved at generation time (`VS_BENCH_REPO` override); D-11 identity figure corrected to ~90 % with citation; D-12 one display name per virus (HIV-1, HTLV-1, `Human papillomavirus N`); EEEV output-name collision resolved by `hostresponse._disambiguate_names`. D-13 stays with `AUD-15`.
- [x] `AUD-22` — docs and governance pins for `AUD-06/07/11/10/19/20` and `AUD-16` (anello breadth max columns); `test_emitted_columns_are_documented` was red on `infected_subthreshold_cells` and `alignment_breadth_max`.
  Done 2026-10-11 (`codex/aud-doc`): `docs/output_reference.md` (detection gate, infected rule, programme `not_applicable`, host-response threshold, `alignment_breadth_max`/`_unique_max`), `README.md`, `docs/{cli_reference,installation,quickstart,reference_panel}.md` (panel needed only for an index without a build manifest; NCBI email only on a cache miss; resume fingerprint inputs; `--ensembl-release`; `--yes`). Governance re-pin by hand, no generator: rows `build-reference-command`, `output-reference-doc`, `public-readme` and claim `v3-reference-homology` take `git_sha` = last commit touching the file and `artifact_sha256` = `sha256sum` of it; the `claim-registry` row is pinned last, to the commit that holds `registry.json` (three commits).
- [~] `AUD-15` — `MECH-E` reference pipeline refactor (user-approved for autonomous work 2026-10-11): source adapters, one gate, masking, GTF emitter; touches pinned `build_reference.py`; mock-verified only.
  **Decision (user, 2026-10-11):** refactor plus manifest schema additions allowed (D-list record, source adapters; bump `reference_manifest.schema.json`, governance re-pin, docs); a golden test must prove byte-identical FASTA/GTF/t2g on fixture inputs for everything not intentionally changed; mock-verified only. Also fold in D-13 (the anellovirus origin-wrap gene is transposed in the `kb ref` cDNA).
  Partial 2026-10-11: one `_prepare_viral_fasta` masking gate, `_finish_reference_build` index/receipt boundary, `source_adapters` and `genome_dlist_applied` in the manifest (`reference_manifest.schema.json` now accepts 3.0.0 and 3.0.1; builder writes 3.0.1), golden test `tests/test_reference_pipeline_aud15.py` (13 cases, byte-identical FASTA/GTF/t2g on fixtures). Then ab9d9da: the GTF emitter is one boundary, `build_reference._emit_viral_gtf` (catalogue and combined adapters; `_anellovirus_gtf` and `_gtf_from_merged_fasta` are thin wrappers). **D-13 closed as a documented limitation, 2026-10-12 (`codex/aud-15c`):** `ngs_tools` (the splitter `kb ref` runs) merges exons into a sorted `SegmentCollection`, so `join(2677..2824,1..80)` of KU243129.1 is indexed as `[1:80]+[2677:2824]` whatever order the GTF gives; this cannot be fixed in ViralScan's GTF without rotating the genome or splitting the gene, both of which change coordinates and the transcript model, so it needs an `ngs_tools`/kb-python change. The rows carry `viralscan_norm "wrap_native"` (as EBV LMP-2A/2B do) and `test_origin_wrap_gene_is_tagged_and_kb_splitter_transposes_it` runs the same splitter on a fixture, keeps the GTF in transcript order (ANELLO-05, never re-sorted) and fails if kb-python ever fixes the order. Impact: one gene, junction-spanning reads only. Still open: governance re-pin (user, once at the end).
- [~] `AUD-16` — `ANDET-01/02/04` follow-ups (user-approved 2026-10-11): per-accession genome-coverage breadth, host-homology coordinates kept for later demotion; auto-evidence stays opt-in.
  **Decision (user, 2026-10-11):** ANDET-02 = keep per-locus host-homology coordinates and add an `eve_risk` column, NO demotion until REF-08 calibrates thresholds; ANDET-01 per-accession genome-coverage breadth and ANDET-04 as in the row (auto-evidence stays opt-in).
  Partial 2026-10-11 (`codex/aud-a1`): ANDET-01 per-virus best-genome breadth done (`alignment_breadth_max`, `alignment_breadth_unique_max`, circular-safe; reporting only). Still open: ANDET-02 per-locus host-homology coordinates and `eve_risk` column; ANDET-04 as in the row.
  ANDET-02 coordinates done 2026-10-11: the builder keeps every host-homology locus (`host_homology_loci` in the reference manifest, PAF validated), detection writes `results/host_homology_loci.tsv` with an `eve_risk` column; no demotion (waits for REF-08). Open: ANDET-04 as in the row.
  2026-10-12 (`codex/aud13-andet-followups`): `results/host_homology_loci.tsv` is now documented (pinned `output_reference.md`, re-pin pending) and a test pins that every column is. The row stays `[~]` (ANDET-04 opt-in, demotion on REF-08).
- [x] `AUD-17` — `M-02`/`M-05` (user-approved 2026-10-11): per-marker polyA/capturability fields in the catalogue and a design note for a competitive host+virus reference option (no M-05 code authorized).
  M-02 done; M-05 design only 2026-10-11 (`codex/aud-m2`): M-02 `gene_programs.tsv` gains `polya_capture` (yes/no/unknown) and `capture_basis` (generator `extras/build_gene_programs.py`; `yes` is RefSeq-annotation-implied, not measured; only HSV-1 `LAT` has a verified override, `unknown`, PMID 38445887; HHV-6B ID-only rows and HSV-2 `LAT` stay `unknown`; no marker is `no` because EBER1/2 are not catalogue markers). Reader/validator in `gene_programs.py` (column optional, absent = all `yes`; non-`unknown` needs a basis). Strict-rule consistency (AUD-19): a *productive* marker with capture not `yes` makes `reactivating`/`latent` read `indeterminate`; any such marker adds a `caveat`. No current call changes (EBV/HCMV markers all `yes`; viruses with `unknown` markers are already latency-not-observable). Docs: `output_reference.md` defence 4 (prose, no new columns emitted), CHANGELOG. M-05 design note only: `docs/plans/competitive-host-virus-reference.md`, no code. Regenerating the TSV also rewrote 150 lines of it (new columns). Not done: measured per-transcript polyA evidence (needs 3' read-end data). Resume review: summary assessability also honors productive-marker capture uncertainty when index-resolution metadata is omitted; six regressions cover capture values and resolution metadata (140 focused programme tests pass).
- [x] `AUD3-02` — `P5-01`/`P5-07`/`P8-01`: the UMAP neighbour-enrichment permutation test must not use an embedding that contains the viral genes that define the labels. Locked option (d): keep the viral-inclusive UMAP figure (labelled so), run the test on a host-only kNN graph (host size factors, HVGs without viral genes, PCA), reword glossary/FAQ, skip cleanly with zero HVGs.
  Done (`codex/aud3-02-umap-test-graph`): `scripts/umap.py` `host_only_pca` (viral columns dropped first, host size factors, own HVGs, PCA) feeds `viral_neighbor_enrichment`; the UMAP/figures are unchanged but titled and annotated viral-inclusive; `conclusion_text` states the tested quantity and caveat; `_n_hvg` guards zero/one-HVG PCA (warning, clean skip; with a virus but no host graph the figures are kept and the annotation says the test was not run). Glossary in `report.html.j2`, `docs/faq.md` and `docs/output_reference.md` reworded. Tests: `tests/test_umap.py::TestViralFreeTestGraph`. Pinned inventory rows affected (re-pin pending, not edited here): `faq-doc` (`docs/faq.md`), `output-reference-doc` (`docs/output_reference.md`).
- [ ] `AUD3-02-OPEN` — open item: the viral-free test is still not calibrated. On pure-noise host cells with random virus-positive labels the audit saw p < 0.05 in 4/15 runs without viral genes in the embedding (nominal 5%), cause unexplained; shown by `biomni_bundle/plans3/verification/repro/P5-01_viral_in_embedding.py` (the `viral_in_embedding=False` rows). The test is descriptive only until this is understood; do not read p < 0.05 as evidence of anything.
- [x] `AUD3-06a` — host-response warnings and relabelling (external audit 3, verified findings P5-02/P5-03/P5-05/P5-06). Warnings, docs and labels only: no metric, default, file or column name changes.
  Done 2026-10-15 (`codex/aud3-06a-hostresponse-warnings`): `hostresponse._warn_cv_leakage` warns when `--groups` is given without `--cv group` (it was ignored) and once when a donor/sample-like obs column exists under cell CV (default CV unchanged; a later unit); `_run_l2_regression` logs one aggregated warning per virus when any cell-CV test fold has < 5 positives (`MIN_TEST_POSITIVES`; ~10 positives + top-50%-depth filter gives 1 pos x 1 neg per fold); `_detect_and_normalize` docstring now says max > 10 (the code) and the skip is a WARNING; stability selection is relabelled "uncalibrated selection frequency" (help, docstrings, logs, README, FAQ, output/CLI reference; `*_stability.csv`, `stab_prob`, `stable` unchanged because `hostresponse_ebv_matched.py` reads them) with a warning above `STAB_SELECTION_MAX_POSITIVES = 300` (measured noise-gene rates at p=1500: 0 / 0 / 1.5 % / 25 % at 30 / 100 / 400 / 1000 positives, so 300 sits between the last clean and first contaminated point; not calibrated). README and `docs/PUBLICATION_READINESS.md` host-response statements marked cell-level CV. Tests `tests/test_hostresponse_warnings.py` (a golden test pins the unchanged outputs). `docs/cli_reference.md` regenerated. Pinned in `analysis/v3_artifact_inventory.tsv` and needing the later re-pin: `README.md` (public-readme), `docs/faq.md`, `docs/output_reference.md`. `claims/registry.json` untouched.
- [x] `AUD3-06b` — host-response `--cv auto` default (external audit 3, verified finding P5-02; LOCKED user decision 2026-10-15 Q3): group CV when a grouping column exists, else cell CV with a loud warning; `--cv cell` stays the explicit legacy mode, `--cv group` stays fail-closed.
  Done 2026-10-15 (`codex/aud3-06b-cv-auto-default`): new `--cv auto` is the default of `viralscan hostresponse`, `run_hostresponse(cv_mode="auto")`, the pipeline `--hostresponse-cv` / `RunConfig.hostresponse_cv` and `menu` rerun; `CV_MODES = (auto, cell, group)`. Resolution (`hostresponse._resolve_auto_cv`, before any fit): `--groups COL` -> exactly explicit group (fail-closed); no `--groups` + exactly ONE donor/sample-like obs column (`_donor_like_columns`, the detector `_warn_cv_leakage` already used, now shared) -> group on it with an INFO log; several or none -> cell + loud warning (none-found warning only for auto, an explicit `--cv cell` stays as quiet as in 06a). Fallback to cell, with a WARNING and the reason recorded, applies ONLY to the auto-selected column: cheap preconditions up front (missing/blank values, raw depth, groups < `min_groups`), and a run-level replay as cell when no virus could run group CV (`group_cv_failed`/`skipped_min_groups` and no metrics rows; one private `_auto_fallback` parameter re-enters `run_hostresponse`, so outputs are the pure cell run). Explicit `--cv group` / `--cv auto --groups` never fall back. `validate_evaluation`: `auto` + `permutations`/`panel_in_fold` requires `groups_column` (fail fast, never inferred). Output rows keep the concrete `cv_mode`/group fields; only the manifest gains `cv_mode_requested` and `cv_auto_fallback` (no schema or consumer reads these; no test compared provenance key sets). `scripts/hostresponse_ebv_matched.py` calls `run_hostresponse` without `cv_mode`: single sample, no donor column -> cell + warning, outputs byte-identical (test). Tests: `tests/test_hostresponse_auto_cv.py` (13); `test_cell_mode_is_default_and_unchanged` in `test_hostresponse_groups.py` was deliberately rewritten to `test_cell_mode_is_explicit_legacy_and_unchanged` (default is now auto -> group on `sample_id`; explicit `cell` output unchanged); `test_hostresponse_warnings.py` pins `cv_mode="cell"` (golden test unchanged); strata parser default test now expects `auto`. Docs: CLI help + regenerated `docs/cli_reference.md`, `docs/faq.md`, `docs/output_reference.md`, README (two cell-level-CV statements updated; `docs/PUBLICATION_READINESS.md` left as is, it describes a past cell-CV number). Pinned in `analysis/v3_artifact_inventory.tsv` and needing the later re-pin: `README.md` (public-readme), `docs/faq.md`, `docs/output_reference.md` (`docs/cli_reference.md` is not pinned). `claims/registry.json` untouched.
- [x] `AUD3-09b` — audit-3 `FN3` + `P4-01` + `P3-01` + D7: (1) the post-count chemistry gate no longer errors on every normal run of a virus-only index without `--host-filter` (`host_in_index` now also needs the Virus Identity table's `index_kind != virus_only`); (2) a warn-level `strand` finding for an explicit `--strand forward` where the rate cannot be judged (two-step or virus-only); (3) docs only: DO-NOT-CITE banner on the ~18x/~3.7x reference-strategy figures, the P3-01 representative-read bias in `docs/output_reference.md`.
  - Code: `analysis.run` now writes the identity table first, then `_chemistry_sanity(ctx, index_kind)`; `chemistry_check.sanity_gate` gains the `strand` warning (`elif not host_in_index and strand == "forward"`; no new thresholds, no error path, `--require-chemistry-sanity` never stops on it). Tests: `tests/test_software_identity.py` (virus-only no error; combined still errors; two-step forward warns), `tests/test_chemistry_check.py::TestSanityGate`; failed before the change, pass after.
  - Non-obvious: unknown index kind (no t2g, no table) keeps the old conservative reading (host-containing). Because the table is built first, a "no viral gene in the index" `ValueError` now precedes the sanity receipt for that run. The strand warning cannot tell 5' from 3' (no library-end evidence in a forced-forward run); it says to ignore it for 3'. The warning also fires on a virus-only index with forward strand. D7: the benchmark file's "12 of 12" (`fresh12b`, 2026-07-04) vs `protocol.yaml` PFR-01 "4 of 12" (2026-06-28 `fresh12`) is flagged in the banner, not resolved; `protocol.yaml` untouched.
  - Inventory (`analysis/v3_artifact_inventory.tsv`, not edited): `docs/output_reference.md` is pinned by row `output-reference-doc` (id 1.1.0, sha256 prefix e793dae3...), which needs re-pinning; `analysis/reference_strategy_benchmark/*` is not pinned by any row. No frozen definition or threshold changed.

### WP4M — Longitudinal compiler (Track B, discovery study M7; new 2026-10-13)

Spec lives in the discovery repo (`plan/IMPLEMENTATION_PLAN.md` M7, `CLAIM_LEVELS.md`); only the package-side tool is tracked here. Consumes the ADJ-01 adjudicated candidate table; numeric cutoffs stay out of the code (LOCK-01).

- [x] `LONG-01` — `scripts/longitudinal_compile.py`: per (subject, virus) trajectories ordered by an explicit visit order, depth-normalised burden primary (viral unique molecules / host UMIs, else host reads; positive-cell fraction and locus breadth beside it, raw counts secondary), event per consecutive visit pair: `new_detection`, `burden_increase`, `burden_decrease`, `persistent`, `cleared`, `not_testable`. Writes `longitudinal_trajectories.tsv`, `longitudinal_events.tsv`, `longitudinal_manifest.json`. Tests `tests/test_longitudinal_compile.py`.
  - Decisions: (1) `burden_fold_change`, `min_host_umis`, `min_detected_level` etc. are config keys with no baked-in pass/fail number; an unset fold change makes a both-detected pair `not_testable`, not `persistent`. (2) Libraries of one (subject, visit) are pooled; a subject with several visits and no `timepoint_order` / `biopsy_order` / `time_from_transplant_days` / numeric `timepoint` is skipped and listed in the manifest, never ordered by guess. (3) Absence is only called (`cleared`, `new_detection`) at a visit whose library is testable (host depth present, optional `min_host_umis`, no `not_testable_flags`). (4) A pair whose visits carry different `thresholds_sha256` or depth bases is `not_testable`. (5) `l4_eligible` is stricter than "no germline/screening flag and level >= L2": it is also false unless the event is `new_detection`/`burden_increase` and the contamination check came back `clear` (an event never makes L4 by itself; persistent/chronic RNA is not reactivation, CLAIM_LEVELS). Pooled/unknown-subject libraries are excluded from trajectories but still count as "other subject" neighbours.
  - 2026-10-13 (Q62, user-confirmed): defaults accepted: `min_detected_level` L1; pooled libraries are excluded from trajectories but still counted for contamination
    (`LONG-02`); the depth table comes from the registry plus run outputs. No code change.
- [x] `LONG-02` — contamination demotion in the same script: a `new_detection` or `burden_increase` becomes `demoted_contamination` (`original_event`, evidence column) when a library of another subject in the same `run_id`, or the same `lane` (unless the run ids are both known and different), holds the same `molecule_key` or `locus_id` of that virus as a library of the later visit (`contamination_match`, `contamination_min_shared` in config). `contamination_check` is `clear`, `demoted` or `not_testable:<why>` (no key table, no keys for the library, no run or lane); `not_testable` never demotes but blocks `l4_eligible`. Deterministic, input-order independent (tests).
  - Open for LOCK-01 / the user: whether sharing a `locus_id` (not a molecule key) across subjects should demote (a locus-level overlap is weak for a ubiquitously detected virus; `contamination_match` can be set to `[molecule_key]`); whether the lane id must be flowcell-qualified in the registry (`SAMPLES_SCHEMA_v2` `lane` is a bare FASTQ-header lane today); the exact ADJ-01 flag vocabulary for `not_testable_flags`.
  - Q46 (user, 2026-10-13): default `contamination_match` is `[molecule_key]`; a shared locus alone does not demote (expected for common viruses such as EBV). `[molecule_key, locus_id]` stays available in config.

## WP5 — Build and validate the truth panel

Objective: create deterministic read/molecule truth across supported chemistry
and ambiguity regimes. Estimated effort: 1-2 engineering weeks plus compute.

### WP5A — Generator

  Progress 2026-10-11: fixed A-06/07/08/10/13/14/17, A-09/A-11 tests, B-05..B-09, B-12/13/14/17, B-21, B-24, B-27/28, B-34, C-07..C-10, C-12 (safe half), C-15, C-16. Open: see `findings.tsv` (rerun the repros at close-out).
  **Decision (user, 2026-10-11):** M-02 = per-marker polyA/transcript-end/capturable fields in the catalogue (default `unknown`, populated only where verified) plus docs; M-05 = design note in `docs/` only, no code (it changes the index build).
- [~] `VAL-01` — implement a seeded generator spanning viral abundance, infected-
  cell fraction, host homology, sibling viruses, low complexity, PCR duplication,
  CB/UMI collisions, ambient/index hopping, and 10x v2/v3/Drop-seq geometry.
  - **2026-09-26:** `VAL-01` also unblocks `SENS-06` (WP4C). A generator that
    plants a target at a known abundance is exactly the positive control
    `SENS-04` consumes, so the two should be built together: it turns
    `informative_negative` from always-false into a measured quantity, and
    supplies the `E8`/`D17` probit input at the same time.
  - **2026-10-05:** design draft for review at
    `docs/plans/2026-10-05-val01-generator-design.md`. It covers proposed
    levels for all nine factors and the generation order (opaque IDs → split →
    covariates → partitioned barcode pools → seeds). It also lists seven
    protocol contradictions (D1–D7) as `DEF-00` inputs:
    - the R2.9 real-PBMC background against sample independence;
    - sibling "absent";
    - per-cell vs per-sample abundance;
    - the template-leakage rule;
    - scale: 50 M pairs/sample, with a per-read truth manifest of ~2.5 TB;
    - the location of holdout truth;
    - the REF-08 dependency.

    No code until the user decides on D1–D7.
  - **2026-10-05, D1 decided (user):** real PBMC is the primary background,
    so R2.9 stands. This brings in planted/unplanted twins, partitioning by
    donor library, and the exclusion of `pbmc_10x_healthy_v3` as a
    background. D1.1 is decided too: ≥ 10 donor libraries per chemistry
    (about 30), with a donor-level cluster bootstrap. DEF-00 amends
    `partitions.unit` and `uncertainty.resampling_unit` to match.
  - **2026-10-05, D3 and D5 decided (user):**
    - D3: abundance is total molecules per virus per sample, at levels 1, 3,
      10, 30, 100, 300 and 1000.
    - D5: keep 25 k reads per cell, with truth recorded per read for non-host
      reads and per molecule for host molecules. The compute risk is
      accepted; measure the real per-sample cost before sizing the grid.

  - **2026-10-05, D2/D4/D6/D7 decided (user); all seven are now closed.**
    - D2 adds `one_member_only`; DEF-00 adds the absent-sibling metric.
    - D4 scopes the leakage rule to host and challenge loci.
    - D6: secret holdout salt, opaque read names, shuffled order.
    - D7: REF-06 then REF-07 come first, and VAL-01 coding pauses until the
      homology levels are frozen.

    The decisions feed the DEF-00 amendment.
  - **2026-10-12 audit (`codex/val-software-audit`, code only, no outcome run):** all nine factor families are implemented at tiny
    scale (abundance, infected-cell fraction, host homology via `host_homology_manifest`, sibling regimes, low complexity, PCR
    duplication, CB/UMI collisions, ambient/index hopping, 10xv2/10xv3/Drop-seq geometry; `tests/test_generate_truth_panel.py`, 30
    tests). New `tests/test_val_generator_audit.py::test_abundance_grid_plants_every_target_at_every_level_in_every_chemistry`
    pins the protocol's identifiability requirement: levels 1/3/10/30/100/300/1000 x 6 targets x 3 chemistries, and with
    collisions and ambient/hopping off a sample at level N holds exactly N distinct planted `(cb, umi)` per virus. **Left:** the
    factor levels stay `pending` until G3 and `seeds.generation` is pending (both the user's); the production-scale streaming
    rewrite is sized after the pilot (D5); real cat42b genomes need the cluster (`VAL-08`). **Observation for the reviewer, not
    changed:** with the default contamination, at abundance 1 the single molecule per virus is relocated to ambient
    (`max(2, ceil(f * N))`), and hopped-in copies from high-abundance donors of the same chemistry and partition give a 10xv3
    abundance-1 sample 6 ambient + 49 hopped truth molecules against 0 cell-origin ones, so low-end LOD scoring must use the origin split.
- [~] `VAL-02` — emit paired FASTQs, `truth_manifest.tsv`, read/molecule truth
  tables, barcode/chemistry metadata, input hashes, and a run manifest.
  - **2026-10-12 audit (`codex/val-software-audit`, code only):** at tiny scale the generator emits everything the row names:
    paired FASTQs (gzip `mtime=0`), `truth_manifest.tsv` (one row per non-host read, one per host molecule), the host-only,
    host-homology and mixed manifests, per-sample `barcode_pool.tsv` with chemistry in every truth row, and `run_manifest.json`
    (generator/ViralScan version, seeds, protocol sha256, FASTQ/template/truth hashes, `inputs` hashes for real genomes, on-list or
    background). New tests in `tests/test_val_generator_audit.py` check that every recorded hash equals the bytes on disk, that the
    manifest lists every template and truth file and no stray directory exists, and that every sample has chemistry and
    barcode-pool metadata. **Left:** nothing software-only at tiny scale; the manifest of a production panel (holdout hashes, 50 M-pair
    samples) waits on G3 and the D5 pilot, so the row stays `[~]`.
- [~] `VAL-03` — add synthetic host-only and adversarial GRCh38-homology
  conditions plus reagent/empty-droplet controls when available.
  - **2026-10-12 audit (`codex/val-software-audit`, code only):** synthetic host-only (`synthetic_host_only`: cells plus empty droplets
    carrying ambient host molecules, no viral truth) and adversarial homology (`synthetic_host_homology`: challenge reads from
    homology clusters at identity 0.8-1, `host_homology_manifest`, partitioned whole by cluster) are generated and tested
    (`tests/test_generate_truth_panel.py::test_host_homology_and_mixed_manifests`,
    `::test_artefacts_challenge_and_host_only_samples_yield_no_viral_calls_beyond_homology`, and new
    `tests/test_val_generator_audit.py::test_host_only_dataset_has_empty_droplets_with_host_ambient_and_no_viral_truth`).
    **Left (data):** reagent / no-template controls need a selected public or wet control (protocol datasets are `pending`, user);
    homology levels use the REF-07 table at G3 (user). Not a software item.
- [ ] `VAL-04` — add a checksum-pinned 10x healthy-donor PBMC v3 presumed-negative
  observational control; never label it absolute ground truth.
- [~] `VAL-05` — add planted target and sibling-virus positives for EBV, HHV-6,
  HSV-1/2, KSHV, and a separately scored anellovirus panel.

  - **2026-10-12 audit (`codex/val-software-audit`):** planted positives for EBV, HHV-6A/B, HSV-1/2 and KSHV with the three protocol
    sibling pairs (co-infected, disjoint, one-member-only regimes) are generated and checked
    (`tests/test_generate_truth_panel.py::test_origins_and_factor_constructions`, `::test_sibling_one_member_only_plants_one_member_per_pair`,
    `::test_exact_planted_molecule_recovery_on_the_golden_tiny_panel`; the real-genome path `--viral-fasta`/`--gene-spans` is AUD-18,
    and `scripts/val08_golden.py prepare` now builds those two inputs from a reference FASTA + GTF). **Left:** (a) the genomes planted
    on the cluster are the cat42b ones (`VAL-08`); (b) the separately scored anellovirus panel needs a protocol dataset id and
    anellovirus gene models, a DEF-00 amendment that is the user's (the generator's dataset list is fixed by the protocol;
    `scripts/plant_anello_10x.py` remains the separate anello planter, ANDET). No code change; stays `[~]`.

### WP5B — Scorer and tiny gate

- [~] `VAL-06` — **core done 2026-10-10, pre-G3 at the user's request, code only, no outcome runs** (`scripts/v3_score.py`, `scripts/v3_partition.py`):
  PRF (M3/M4), AUPRC (D6), negative false-call rate (M2), sibling confusion (M5), host-virus allocation (M6), the shared
  largest-remainder split. **Left:** host-homology false positives, burden concordance, calibration, LOD scoring.
  Molecule/cell precision, recall, F1, AUPRC, sibling
  confusion, host-homology false positives, burden concordance, calibration, and
  limit-of-detection scoring.
  - **2026-10-12 audit (`codex/val-software-audit`, code only):** the `Left` list above is partly stale. Already on `main`: host-homology
    false positives (`negative_cell_false_call_rate`, D7), tier specificity/PPV with one twin-exclusion set (D18/D19), cross-tool
    relative difference (M9/D20), the negative-library three-number scoring and T1-T5 tiers (`M04-SCORE`), and the native count-layer
    adapter. New `scripts/v3_lod.py` (kept out of `v3_score.py`, which `M04-SCORE` owns) adds **LOD scoring**: probit fit of detection
    on log10 abundance, LOD95, the no-extrapolation rule (`not_reached` / `below_range`, never a number), `separation` and
    `non_positive_slope` reported rather than regularised, and the BCa unit-level (cluster) bootstrap with the protocol's 2000
    replicates / 95 % / minimum three units, an interval withheld with a reason when any replicate or jackknife fit is not estimable.
    `seed` is a required argument (`seeds.bootstrap` is frozen at G3). Tests `tests/test_v3_lod.py` (14: exactly-probit data, two-level
    hand-computable fit, status table, whole-unit resampling, reproducibility). Count conservation (M1) is checked per run by
    `scripts/val08_golden.py`. **Left:** (a) tier calibration error (M8) needs the nominal confidences and the
    `calibration.threshold_search` result (protocol `pending`, user); (b) D21 per-10k burden needs a frozen common host-plus-viral
    denominator from a real run; (c) M7 two-step loss needs two-step output on identical input (cluster); (d) wiring the LOD
    estimator to a run's per-sample detection calls waits on the tier cutoffs (G3). No acceptance threshold anywhere. Row stays `[~]`.
- [x] `VAL-07` — **closed 2026-10-12 for the scorer surface on `main` (see the audit note below the row); the 2026-10-10 text follows.** **hand-computed fixtures done 2026-10-10** (`tests/test_v3_scorer.py`, 24 cases, incl. truth-manifest rejection: missing
  contract column, planted row without virus/gene, artefact row labelled as truth, bad partition/boolean). **Left:** reject
  denominator, feature, barcode, or count-layer mismatch (needs the layer map and feature universe of a real run). Validate the scorer against hand-computed fixtures and reject
  denominator, feature, barcode, or count-layer mismatch.
  - **2026-10-12 closure evidence (`codex/val-software-audit`):** acceptance = hand-computed fixtures plus rejection of denominator, feature,
    barcode and count-layer mismatch. Fixtures: `tests/test_v3_scorer.py` (PRF, AUPRC, M2/M5/M6, partition, manifest, D7, tier
    specificity/PPV, M9), `tests/test_v3_scorer_m04.py`, `tests/test_v3_lod.py` (a two-level fit computed by hand); 86 tests pass.
    Rejections: denominator `test_native_output_adapter_rejects_frozen_declaration_mismatches[denominator]`, barcode
    `[barcodes]` (missing, duplicated) and `test_tier_scoring_rejects_barcode_denominator_mismatch[truth|calls|twins]`, feature
    `[features]` (missing, empty), count layer `[layer]` plus `test_native_output_adapter_checks_persisted_count_contract` (X mismatch,
    fractional unique layer, method, unit, missing layer, duplicate barcode), digest/fingerprint/receipt pins, and on real
    kallisto/bustools/ViralScan output `_native_host_anchor` in `tests/integration/test_truth_panel_exact.py` (2 passed, 2026-10-12).
    Metrics added to the scorer later (`M04-SCORE` extensions, M8, D21, M7) carry their own fixtures under their own rows; this row
    covers the surface that exists now. Original `Left` text above is superseded.
- [~] `VAL-08` — run the golden tiny panel end to end and prove exact planted-
  molecule recovery, count conservation, determinism, row-order invariance, and
  chunk-size invariance.
  - **2026-10-12 audit (`codex/val-software-audit`, code only, no outcome run):** the local `kb count`/ViralScan arm already exists
    (`tests/integration/test_truth_panel_exact.py::test_full_viralscan_truth_and_conservation`, the Bundle V note below is stale on
    this point): exact `(cb, umi, virus)` recovery on a collision-free golden, collision-challenge losses kept in recall's
    denominator, audit conservation, and reference/FASTQ order plus gzip-member chunk invariance. New `scripts/val08_golden.py`
    makes the cluster arm reproducible: `prepare` (reference FASTA+GTF+contig map -> `--viral-fasta`/`--gene-spans`, validated by
    the generator's loader; spliced, short and unmapped genes skipped and listed), `whitelist`, `permute` (record-reversed,
    multi-member gzip), `plan` (the 11 ordered command lines: generator, ViralScan x3, check x3, compare x2; they parse with the
    real generator and `viralscan` parsers), `check` (exact recovery, no host-to-viral, count conservation from the molecule audit,
    `--reference real` for a reference without the generator's host templates) and `compare` (determinism, row-order and
    chunk-size invariance = identical assignments, counts and audit on different inputs). Tests `tests/test_val08_golden.py` (23,
    defects injected into a stand-in run) and the integration test now also calls `check_run`/`compare_runs` on real kallisto/
    bustools/ViralScan output (2 passed, 2026-10-12, 372 s). No acceptance threshold added: exactness is the row's definition.
    **Left (cluster, user):** running `plan` against the cat42b FASTA/GTF, the on-list and the cluster index, and reading the
    `check`/`compare` output; reference-order invariance on the real index needs a second index build (not planned here).
**Bundle V (2026-10-11, `codex/aud-v`, code only, no outcome run):** `scripts/generate_truth_panel.py` (seeded tiny-panel generator with a synthetic background; required `--generation-seed`; `--production` refuses until G3), `v3_partition.unit_id`/`split_libraries`, `tests/test_generate_truth_panel.py` (30) and `tests/integration/test_truth_panel_exact.py` (real kallisto/bustools: planted (cb, umi, virus) set recovered exactly). VAL-01/02 done at tiny scale; VAL-03 lacks reagent controls; VAL-05 uses random stand-in genomes (no reference gene models, no anellovirus panel); VAL-08 still needs the `kb count`/ViralScan arm. **Decided by the user 2026-10-11:** host molecules stay in `truth_manifest` as molecule rows (`host_only_manifest` = synthetic host-only dataset only); singleton strata go to training inside the generator with `apportion` unchanged; add `--viral-fasta` + gene-span TSV (sha256-recorded) and `--onlist` options (follow-up bundle V2). Original note: (iv) host molecules live in `truth_manifest` (`record_granularity=molecule`), `host_only_manifest` holds only the `synthetic_host_only` dataset; (i) singleton strata go to training inside the generator and `apportion` is unchanged (seat total = round(0.3 * N_splittable)); optional `--viral-fasta` + gene-span hook to plant real genomes not added; sample output depends on the other samples in the run (matters for a per-sample cluster rewrite).
  - 2026-10-13 (Q69, user-confirmed): the VAL truth panels are generated by an agent after the WP5 build, against the frozen reference (not before); the cluster arm
    of this row therefore waits for WP5.
- [~] `VAL-09` — freeze generator/scorer version, seeds, manifests, and expected
  tiny outputs before cluster-scale execution.
  - **2026-10-12 audit (`codex/val-software-audit`):** a freeze is a user act, so this row stays `[~]`. What exists: the generator version
    string (`viralscan-truth-panel/1`) and the tiny-panel expected outputs are pinned by hash in
    `tests/test_generate_truth_panel_real.py::GOLDEN_DEFAULT` (default and noisy panels, truth manifests, templates, per-sample FASTQ
    digests), the run manifest records generator/ViralScan version, seeds, protocol sha256 and every output hash
    (`tests/test_val_generator_audit.py`), and `scripts/val08_golden.py` fixes the check list. **Blocked on the user:** `seeds.generation`
    and the factor levels are `pending` in the protocol until G3, so nothing can be frozen; the scorer has no version string
    (`v3_score.py` is `M04-SCORE`'s); the cluster golden receipt (`VAL-08`) supplies the real expected outputs. Do not add a freeze file
    before G3: its seed and hashes would be provisional.
- [ ] `VAL-10` — run the full training and untouched holdout panels, retaining
  failures rather than silently dropping conditions.

`G5a` passes when the golden tiny panel is exact and deterministic, the full
panel is manifest-complete, and every condition has a result or reproducible
failure record.

## WP6 — Run public positives and comparators

Objective: regenerate all eligible evidence with v3 and make cross-tool
comparisons denominator-, annotation-, barcode-, and layer-matched. Estimated
elapsed time: 2-4 weeks, dominated by downloads, queues, and comparator setup.

### WP6A — Public v3 reruns

- [ ] `RUN-01` — HHV-6B `SRR20710641` under measured 10x v3 geometry.
- [ ] `RUN-02` — EBV `SRR12682296` under measured 10x v2 geometry.
- [ ] `RUN-03` — HSV-1 `SRR8315713` under measured Drop-seq geometry.
- [ ] `RUN-04` — KSHV latent/lytic series `GSE190558`.
- [ ] `RUN-05` — KSHV plus EBV replicates `GSE154900`.

Each sample must produce a validated output tree plus input/reference hashes,
command manifest, tool versions, scheduler accounting, and failure log. Expected
virus recovery is contextual evidence, not ground truth.

### WP6B — Harmonized workflow matrix

- [ ] `CMP-00` — replace historical output-derived barcode unions/intersections
  with the frozen pre-outcome anchor/feature manifests and audited count-layer
  adapters; use `counts_unique`, retain structural zeros only from declared
  complete barcode domains, and record every missing/extra mapping.
- [ ] `CMP-01` — ViralScan combined and exact-fragment STAR two-step; keep
  kallisto two-step excluded unless exact fragment lineage becomes available.
- [ ] `CMP-02` — STARsolo combined and STAR host-filter/two-step.
- [ ] `CMP-03` — traditional host-genome subtraction followed by viral alignment,
  with and without CB/UMI retention.
- [ ] `CMP-04` — Venus in an isolated version/digest-pinned environment with no
  silent algorithm patch.
- [ ] `CMP-05` — Viral-Track in an isolated version/digest-pinned environment with
  no silent algorithm patch.
- [ ] `CMP-06` — (added 2026-09-30, user request) head-to-head with the
  **original ViralScan 2.2.0**, as released and run by evonk, on the **same
  reference** as v3.
  - **Scope.** The comparison measures *implementation* differences only. It
    runs 2.2.0 on the same inputs as v3. GOV-06, by contrast, compares
    *archived* 2.2.0 outputs.
  - **Reference (user decision 2026-09-30): one arm only.** Both versions use
    the **latest reference**: today the final 2,343-genome panel,
    `$SITE_ROOT/viral_ref_final/build/{panel.idx,panel.t2g}`.
    *2026-10-03:* the latest reference is now
    `$SITE_ROOT/viral_ref_cat42b/build/{panel.idx,panel.t2g}`
    (CAT-42: same 2,343 genomes, homopolymer and anellovirus low-complexity
    masks). Its anellovirus gene IDs carry the accession prefix.
    *2026-10-05:* the latest reference is now
    `$SITE_ROOT/viral_ref_cat42d/build/{panel.idx,panel.t2g}`,
    which is cat42b plus a GRCh38 genome D-list (REF-06 gates passed).
    If the frozen G4 reference supersedes it, both versions switch to that.
    - Give both versions the same viral-only GTF, extracted from
      `combined.gtf`.
    - Never pass the combined host+viral GTF as `-gtf`: 2.2.0 would then count
      every host gene as viral (the Q8 issue).
    - There is no native-Serratus-index arm. Reference differences are out of
      scope.
  - **What 2.2.0 is.** The install is in env
    `$COLLAB_ROOT/conda/envs/test_viralscan`
    (`viralscan-2.2.0.dist-info`); evonk's run pattern is
    `$COLLAB_ROOT/viralscan/run_viralscan.sh`.
    Pin the version and the environment spec before any run. First check that
    2.2.0 accepts the 8-column kb-ref `panel.t2g` and its gene IDs; record
    whatever adapter is needed.
  - **Hold everything else equal:**
    - identical FASTQ hashes;
    - identical `-x` and whitelist;
    - the same strand. 2.2.0 has no `--strand`, so compare at kallisto's
      default and, once DEF-02 exists, report v3 at the inferred strand as a
      separate row;
    - the same cell-calling anchor where 2.2.0 allows one.
  - **Where it runs:**
    - **context** (already-seen data, R2.2, never tuning):
      - EBV SRR12682296 (10xv2);
      - HSV-1 SRR8315713 (Drop-seq);
      - HHV-6B SRR20710641 (5′, `-x 10xv2`);
      - evonk's own EBV SRR6825024;
      - covid x213/x216 and tonsil x223.
    - **confirmatory:** the truth-panel holdout, as an SCI-04 row.
  - **Report, per dataset:**
    - total and per-virus molecules;
    - called cells and the infected-cell denominator;
    - virus names and grouping;
    - false calls on the host-only and planted-homology negatives.
  - **Attribute every difference to a known implementation mechanism:**
    - barcode correction (SW-13/SW-20, F-018);
    - multimap allocation and siblings (F-017);
    - poly-G/low-complexity handling (F-019);
    - the strand default (F-020);
    - naming and grouping (MECH-A);
    - the summary headline (SW-16).

    Any unattributed difference is a finding in its own right.
  - **Dependencies:**
    - the context runs can start now as exploratory;
    - the confirmatory run needs G3/G4 and the frozen defaults (WP1E R2.1/R2.4);
    - the budget comes from the SCI-04 comparator allotment (WP1E R3.3).
Primary comparison rules: identical FASTQ hashes and viral sequences, same host
release, outcome-independent cell anchor, audited feature intersection, and
unique-only molecule parity. Recommended ambiguity-aware outputs are a separate
secondary analysis. Run at least EBV, HHV-6B, and HSV-1 for both dedicated
comparators.

`G5b` passes when every preregistered row is complete or transparently failed,
raw outputs and commands are retained, and an independent audit finds no
outcome-selected barcodes or mismatched denominator, annotation, feature, or
count layer.

### WP6C — Dataset run plan (new 2026-10-05)

Plan: `docs/plans/2026-10-05-dataset-run-plan.md` (per-dataset run specs,
blockers, order). User rule (2026-10-05): every dataset is analysed on the
latest ViralScan version, recorded as version + commit. Every existing run is
`3.0.0.dev0`; latest is `3.0.0.dev1 @ 12ceac7`. Experiments stay deferred
behind the package-completion plan unless the user releases a row.

- [x] `DSR-01` — viral-only GTF extracted once from cat42d `combined.gtf`.
  - **Done 2026-10-06:** `viral_ref_cat42d/build/viral_panel.gtf` (sha256
    `d80b1fdb…`; the `combined.gtf` rows whose seqname is not `ENST*`). 2,343
    seqnames, matching `viral.fa`. It is byte-identical to cat42
    `viral_panel.gtf`, so every earlier cat42d run already used it.
  - Against `panel.t2g` it has 1 naming mismatch (`ROTA_A_RVA_s4_gp1` vs
    `ROTA_A_RVA_s4_g p1`, a space in the gene_id), plus
    `HUM_PARVO_unassigned_gene_1`, which has no indexed transcript.
- [x] `CHEM-01` — **closed 2026-10-10** (run_manifest record and `--require-chemistry-sanity` exist; bulk stays advice-only). Library diagnostics (user request 2026-10-06): `viralscan
  check-chemistry` plus a post-`kb count` gate (`src/viralscan/chemistry_check.py`,
  tests in `tests/test_chemistry_check.py`).
  - **Single-cell vs bulk:** barcode repeat rate (best of 16 and 12 bp) and a
    continuation rate. Reads sharing a prefix that also share the next UMI-length
    bases are cDNA, not barcodes. Measured on real R1: EBV 0.31, HHV-6B 0.43,
    Drop-seq 0.02, tonsil 0.24; cDNA mates (R2) 0.93–0.97. A cDNA read passed as
    `-s1` is called "bulk (or -s1/-s2 swapped)".
  - **3′ vs 5′:** R1 TSO/poly-T decides; on trimmed R1 (tonsil 28 bp) the
    strand pilot decides (F-020 shape); a conflict is reported, never resolved.
  - **Post-run gate:** warn-only in `scripts/analysis.py`: p_pseudoaligned
    < 40 % warns, < 10 % flags an error, skipped after a host filter
    (reads are already host-free).
  - Done: the diagnosis is in `run_manifest.json` (`analysis.py:140-159`); the gate is fail-closed behind
    `--require-chemistry-sanity`. Decision: the bulk path only advises (`kb count -x BULK` is not wired).
- [x] `CHEM-02` — TusslChemistry (`tchem`, Apache-2.0 Rust binary, v0.1.0) as a second opinion on the chemistry
  (user request and decisions 2026-10-17/18: second opinion, prebuilt release binary, include in viralscan, use it to
  resolve the Track A chemistry). Branch `codex/tchem-second-opinion` (worktree `vs_work/tchem`, off `codex/drop-py39` 5e551c9).
  - **Code:** `src/viralscan/tchem.py` (binary via `--tchem-bin`, `$VIRALSCAN_TCHEM`, `PATH`; lists via `--tchem-whitelists`,
    `$VIRALSCAN_TCHEM_WHITELISTS`, then tchem's own `$TUSSL_CHEMISTRY_WHITELISTS`; two calls per sample: `--format json` for
    the verdict, `--emit kb` for `-x/-w/--strand`; compared on barcode+UMI geometry), flags `--tchem {auto,require,off}`,
    `--tchem-bin`, `--tchem-whitelists` on the run and on `check-chemistry`, a `tchem` block per sample in the manifest,
    generated CLI reference, changelog and `config/public_ship_scope.json` (tchem.py in both lists) updated. Policy (user): agreement passes; a disagreement with `-x`, or a tchem
    `ambiguous`/`no_match`, stops the run unless `--force-technology`; a missing tchem never stops it in `auto`.
  - **Verified against the real binary (2026-10-18):** release `tchem-v0.1.0-x86_64-unknown-linux-musl.tar.gz` sha256
    `9112e5018c26868317fa6547490094b59915927ef1482bb361f0b60ce5fe488e` matches the published `.sha256`; installed at
    `mdmanurung/tools/tchem_0.1.0`, 23 barcode lists fetched and packed in `tools/tchem_whitelists` (10x Flex skipped: Cell
    Ranger only). The real `--emit kb` text and JSON (schema 1.1) parse as the code expects (the former open item); the
    captured outputs are fixtures in `tests/test_tchem.py`, plus an env-gated real-binary test. One regression found and fixed:
    `menu._tchem_second_opinion` indexed an empty `blocks` when a forced `-x` skipped detection (KeyError).
  - **Track A re-check (read-only, 80 libraries with reads on disk, `-x` = recorded chemistry):** HNSCC 34/34 and hpv16 2/2
    agree (10x 3' v2 / 3' v3), 8 tonsil libraries agree (3' v3 and Multiome GEX). New calls: `188573`, `182656`, `269550`
    are 10x Multiome GEX, HIGH, forward, `-x 0,0,16:0,16,28:1,0,0` (their earlier "unresolved" was the missing multiome on-list in
    the manifest, not the library). `187505`, `167868`, `187508` stay AMBIGUOUS between 10x 3' v2 and 5' v1/v2 even at 1M reads
    (same 16+10 geometry, 2 trailing bases, cDNA strand unclear); `182636` is 5' v1/v2 and `187507` 3' v2, both MEDIUM. HHV-6B and
    EBV FASTQs are deleted, so they could not be re-checked; HSV-1 Drop-seq libraries read `no_call` (tchem names Drop-seq but
    gives no kb technology).
  - **Left:** decide the five non-HIGH libraries (tchem narrows them to the 10x v2 geometry but cannot call 3' vs 5').
- [ ] `DSR-10` — standard arm set on every dataset (user, 2026-10-06; see
  `.living/decisions.md` "Same standard arm set on every dataset"). Each
  dataset × arm cell is run or marked n/a with a reason in the run plan's arm
  matrix. The two-step arm needs `DSR-11`.
- [x] `DSR-11` — virus-only cat42d index for the two-step arm: same `viral.fa`,
  DSR-01 GTF and GRCh38 D-list as `panel.idx`, kb 0.50.1
  (`scripts/slurm_cat42d_virus_only_index.sh`, job 25702322, started 2026-10-06).
- [~] `DSR-12` — SFL tonsil x223 full native runs, standard arm set (`TONSIL-02`).
  - 1.6 B read pairs (cellranger metrics), so each arm is ~42 alloc core-h
    combined and ~180 two-step by F-024.
  - The `artefact` filter is single-threaded, about 37 h wall.
  - Whitelist: cellranger raw barcode universe (3,491,354 barcodes, built as in
    the covid RUNBOOK), at
    `viralscan_work/sfl_tonsil/ref/cellranger_raw_whitelist.txt`.
  - Code is pinned at `vs_pinned/4346dc8`.
  - Pilot first: job 25702321, 4 M pairs, check-whitelist plus 3 strands.
  - Then the full arms, using the strand the pilot measures.
- [~] `DSR-02` — one parameterised read-validation script (EC → bus records →
  reads → competitive minimap2 vs GRCh38 + panel → per-read verdict),
  replacing the hand-edited `run.sh` copies used for F-025.
  - 2026-10-12 (`codex/dsr-evidence-tooling`): `scripts/dsr02_validate_call.py` is the one script. Steps 1-4
    were already `viralscan evidence` (EC → bus → exact reads → competitive minimap2, molecule verdicts); the
    two sbatch files carried the cluster paths. The script takes one call (`--run-dir/--virus/--out`) or row N
    of `calls.tsv` (`--calls --row`), with `--panel-fasta`, `--panel-for DATASET=FASTA` (per-dataset panel),
    `--host-fasta`, `--dry-run`, `--postprocess-only`; skips a finished directory; writes `code_sha.txt`
    (dossier E08) and drops the ~GB competitive FASTA / BLAST db unless `--keep-intermediates`. New: a
    per-read table `read_verdicts.tsv.gz` + `read_verdict_summary.tsv` (virus_best / host_best / tie /
    unaligned x repeat class `tso` / `telomeric` / `low_complexity` / `none`), so the F-025 telomere reads are
    labelled, not just counted. `molecule_verdict` now shares `side_score` / `verdict_of` with it. Both sbatch
    scripts call the driver. Tests `tests/test_dsr02_validate_call.py` (SAM text, no binaries).
  - Decision: the `telomeric` label (>= 50 % of bases in TTAGGG/CCCTAA hexamers, >= 2 hexamers) is an
    annotation only, PROVISIONAL, and gates no verdict; the user owns any threshold that does.
  - **Left:** run it on a real evidence directory (needs samtools + minimap2 output; first use is the
    `DSR-13` evidence regeneration on a pinned commit). `dsr02_verdicts.py` does not read the per-read table yet.
- [ ] `DSR-03` — HPV16 GSE189670 rerun on the target commit; read-validate
  HPV118 / HPV29, which appear in both the infected and parental raft.
- [ ] `DSR-04` — COVID x213/x216 rerun: cat42d, `--strand reverse`, read
  filter off / artefact, emptydrops.
- [x] `DSR-05` — HHV-6A residual in SRR20710641: read validation (job
  25701668, 10 HHV-6A-only + 45 shared 6A/6B records, 95 reads). Result
  F-026: the residual is HHV-6B (EM bleed + single strain SNPs), no
  read-level HHV-6A. F-025 (HHV-6B in EBV = telomere repeats) came from the
  same method.
- [ ] `DSR-06` — GSE164690 HNSCC one-lane pilot, HN18 CD45− vs HN01.
- [ ] `DSR-07` — E-MTAB-13687 tonsil atlas negative-control subset, samples
  chosen from metadata before any viral output.
- [x] `DSR-08` — record the git commit next to `viralscan_version` in
  `run_manifest.json` / `reference_provenance.json`.
  - 2026-10-10 (wave 4): `run_manifest.json` already stamped `software_identity()`
    (cb8f728; commit, dirty flag, build hash). `reference_provenance.json` now
    carries the same `software_identity` block; `tests/test_reference_provenance.py`
    pins it. Nothing left.
- [ ] `DSR-16` — exploratory round 2 (grill 3, 2026-10-12): explore ledger, parameterised round-2 SLURM scripts and `dsr_round2/GUARD.md`, round-2 vs round-1 comparison report; then the runs on the pinned commit; `SENS-07` as round 2b.
- [~] `T2-AMEND-01` — the V03 (GSE336913, influenza A) protocol amendment, added 2026-10-13 (the row did not exist on `main`;
  created here per CLAUDE.md). User Q79/Q80: register V03 as a CONTEXT row with pending assets, in one record. **Done**
  (`codex/sci05-r12-protocol`, `DEV-027`): the row, the harmonization and `frozen_inputs` digests and both schema copies move
  together. **Not done, by design:** scoring V03 needs `PRE-01`, a 5-prime read model (`supported_scope.five_prime`), an
  outcome-independent read-structure-only gene-expression-experiment rule, influenza A in the frozen panel, a whole-tier
  promotion rule (EBV, HHV-6B, HSV-1, iSLK.219, BC-1, not V03 alone) and a scored-T2 role/endpoint/denominator
  (SCI-05 round 12 B1), and counts-only reporting (three clean libraries). V03 FASTQs were being downloaded on 2026-10-13
  (Track B, data acquisition, no quantification; `PFR-12`); no V03 library may be run before the commit that freezes its row.
  Study-repo `CAL-11` must be gated on that commit (owner: the Codex session in the study repo).
  **Q87 (user, 2026-10-13): the scored-T2 class (B1) is deferred; V03 stays a context row** until the user reopens it.
- [ ] `DISC-01` — broad virus discovery study (grill 3, 2026-10-12): private study repo, tracker corrections, fact-finding (V02 replacement, V04 labels, healthy backgrounds), then registry, preflight, adjudicator (on DSR-02), longitudinal compiler; runs gated on WP5, G3 and `LOCK-02`.
  - 2026-10-13: package tools merged: ADJ-01 (#34), LONG-01/02 (#33), CAL-13 (#35), LOCK-02 schema/validator (#36). PRE-01 preflight is in the study repo, owned by the Codex session there (Q45, user 2026-10-13; check `git status` there before any study-repo commit). Left: `LOCK-01` values, WP5, G3, downloads.
  - 2026-10-13 (Q54-Q57, Q59-Q61, Q67, Q68, user-confirmed): before G3 only data prep (downloads, `PRE-01`, pilots, T2T index), no quantification until G3 and WP5.
    Calibration sets now (V01 local, V02 `GSE226438` 3+3, V03 if raw GEX exists, V04 5 CHB + 5 FC, `CAL-08` backgrounds, N01) and a 1-2 run pilot per
    discovery/replication set; bulk is streamed after `LOCK-02`. Data under `viralscan_work/discovery/{fastq,runs,checksums}`, 90 % disk cap. Host-homology
    screen: T2T-CHM13v2.0 (`GCF_009914755.1`) + minimap2 `sr` index; HPRC for L1+ candidates only. Agents run downloads under a GUARD. `LOCK-01` carries the locus /
    positive-cell / sibling-bleed definitions, the rare-super-expressor route, `CAL-13` strata and decoy ratio (calibration data only); the lock cites the
    sha256 of the protocol file bytes; the lock schema is checkout-only; lane = `flowcell:lane` (`PRE-01`, study repo, Codex).
- [ ] `T2-AMEND-01` — (2026-10-13, from the grill 3 note) influenza `GSE336913` (V03) joins the package's scored T2 rows by an amendment to the validation
  protocol, reviewed in `SCI-05` round 12 (Q66). Pending: V03 needs raw GEX to exist (Q55); nothing is amended before the review; no scored run before G3.
  - SCI-05 round 12 (2026-10-13): Blocker B1, the protocol has no scored-T2 class. **Q80 (user): register V03 as a context row with pending assets (like
    the KSHV rows) in one `DEV-027` amendment** after B1/M5/M6 are resolved; scored only after `PRE-01`, a 5' read model and an outcome-independent rule
    for picking the GEX experiment. Track B still uses V03 as a calibration positive.
- [x] `LOCK-02` (software only; Track B, `codex/lock02-schema`, 2026-10-13) — generic schema and validator for the study's
  `protocol-v1.lock.json`. `schemas/v3/protocol_lock.schema.json` (Draft 2020-12; its own container, not an extension of
  `validation_protocol.schema.json`, which it cites by `package_validation_protocol_sha256`) and
  `scripts/validate_protocol_lock.py` (`--template` prints a skeleton with `null` for every value; checks schema, refuses
  `test_only` without `--allow-test-only`, recomputes every `X_sha256` that has an `X_path` sibling; exit 0 valid, 1 invalid,
  2 unreadable, 3 incomplete). The `adjudication` object is nested exactly as `adjudicate_candidates.py` (ADJ-01) reads it
  (`distributed_route.{min_unique_molecules,min_positive_cells,min_independent_loci}`, `max_sibling_allocated_fraction`,
  `min_virus_best_fraction`, `min_complex_body_fraction`, `l3_min_programme_cells`), so one lock feeds `--lock` directly.
  Tests `tests/test_protocol_lock.py` (TEST_ONLY fixtures; no value in code or template). **Left open: `LOCK-01`** (the
  calibration that writes the numbers and the sha256s) and writing the real lock; the template is incomplete by design.
  Decisions: nulls are schema-valid but make a lock `incomplete`, never `valid`; `test_only` must appear at the top level
  and inside `adjudication` (the adjudicator reads the latter); `status` is `exploratory` until repo G5, then `final`; not
  added to `MANIFEST.in` or the public ship scope (study tooling, run from a checkout).
  - 2026-10-13 (Q56, Q59, Q60, user-confirmed): the lock cites the sha256 of the protocol file bytes (not a re-serialisation) and its schema is checkout-only (never
    shipped); `LOCK-01` fills the numbers from calibration data. No code change in this row.
- [~] `CAL-13` — optional decoy layer for false-discovery calibration (Track B generic tool; discovery plan M4 / decision 6).
  **Software done 2026-10-13 (`codex/cal13-decoys`); no index was built and no real library was scored.**
  `scripts/cal13_decoys.py` with two subcommands. `make`: panel FASTA(s) (+ `--unrelated-host-fasta`, registered unshuffled
  as class `unrelated_host`, no download code) -> `decoys.fa`, `decoys.gtf` (one gene/transcript/exon per decoy, the
  `build-ref` whole-genome writer), `decoys_manifest.tsv` (`decoy_id, class, source_id, seed, k, rep, length, sha256,
  source_sha256`). Decoys are `DECOY_<source>_<rep>` k-let-preserving shuffles (default k = 6, `--copies` for the
  decoy:target ratio), exact k-mer counts and end (k-1)-mers kept; Altschul-Erickson Eulerian path with the last-exit
  tree from Wilson's algorithm (uniform over admissible sequences, tested by enumeration); stdlib only (`ushuffle` is not
  installed and was not added); per-decoy RNG stream from sha256(seed, k, id), so record order does not matter.
  `fdr`: `--calls` TSV (`library_id, target_id, is_decoy, stratum, score`, zeros included) -> `target_pvalues.tsv`
  (`p_empirical = (b+1)/(m+1)` against the same library and stratum; `q_td` stratified target-decoy, +1 conservative,
  decoy:target ratio aware; `q_bh` = BH on the p-values) and `library_summary.tsv` (calls at `--threshold`, expected false
  calls = sum over strata of n_targets * D/m with a Jeffreys-Beta Monte Carlo interval, optional comparison with the
  pre-registered `false_discovery.expected_false_L2_calls_per_library` / `interval` read from `--protocol`; `TBD*`
  placeholders are ignored). Nothing is baked in: threshold and seed are required, strata are caller-supplied.
  Tests `tests/test_cal13_decoys.py` (k-let composition, determinism, hand-computed p/q, stratification).
  - **Decoys do not model host homology, TSO chimeras or poly-A sinks** (masks, D-list and 5-10 real negative
    libraries cover those); they model only composition-matched random k-mer sharing and unrelated-host viruses.
  - Decisions: p and q are per library (the null is depth-matched; nothing pooled across libraries); a decoy shuffle
    keeps its source's length and 6-mer composition, so a decoy inherits its source's stratum from the manifest
    (`source_id`); `--manifest` fails a library that lacks any manifest decoy (sparse input would shrink m).
  - **Left (user / protocol):** pre-register the target and interval and the score threshold (new optional protocol key
    would be `false_discovery.optional_decoy_layer.score_threshold`: not read yet), decoy:target ratio, seeds and
    strata definitions; build the separate decoy index with the same host and D-list; run on 5-10 real negatives.
  - 2026-10-13 (Q56, Q57, user-confirmed): strata and decoy:target ratio are pre-registered in `LOCK-01`, set from calibration data only; the decoys run once in
    calibration (V positives + 5-10 real negatives), not per discovery library. The `score_threshold` key stays unread until `LOCK-01`.
  - Q47 (user, 2026-10-13): keep the conservative +1 target-decoy estimate `(D+1)/(r*T)`.
- [~] `DSR-13` — round 1: all-dataset reruns on `bbf1821` (user, 2026-10-06; plan
  `~/.claude/plans/read-handoff-md-and-refactored-waffle.md`).
  - Frozen inputs: `viralscan_work/dsr_round1/manifest.tsv` (19 samples; GSE190558
    lanes and GSE154900 runs concatenated per GSM by `scripts/slurm_dsr_concat_inputs.sh`,
    lane → GSM checked against ENA `experiment_title`), pinned worktree
    `vs_pinned/bbf1821` (clean, `v3.0.0.dev1-41-gbbf1821`).
  - One input gate: `scripts/slurm_dsr_chem_array.sh` (`check-chemistry`, no `-x`
    passed); `-x` and `--strand` for every arm come from its JSON, and the JSON is
    copied into each run dir.
  - One run script, `scripts/slurm_dsr_rerun_array.sh ARM`: arms differ only by flags
    (`combined_off`, `combined_artefact`, `twostep`). Comparability guard
    (`dsr_round1/GUARD.md`): runners may change resources only.
  - HSV-1 runs all arms with the SW-21 bypass (SW-21 closed 2026-10-03); GSE154900 at
    full depth (user, 2026-10-06).
  - Decisions: ERR13027122 is `n/a` (0.1–0.2 % pseudoalignment on every strand, so
    not a GEX library; probably HTO/ADT). The STAR host index is GRCh38-2024-A genome
    while cat42d's host is Ensembl 116 cDNA, so combined vs two-step also differ in
    host annotation (noted, not a blocker).
  - Chemistry results: EBV 10xv2 / forward; HHV-6B 10xv2 5′ / reverse (settles v2 vs v3);
    HSV-1 dropseq / forward; HPV16 rafts, GSE164690, ERR13027010/027 3′ forward;
    COVID x213/x216 and tonsil x223 5′ / reverse.
- [x] `DSR-14` — `viralscan evidence` OOM (3 jobs: 25702323 at 48 G, 25702383 at 128 G,
  25712429 at 256 G, peak 231 GB). Cause was NOT BLAST (the earlier note was wrong): none
  reached BLAST. `coverage_depth_points()` ran `samtools depth -a` on a BAM whose header is
  all of GRCh38, printing about 1-2 G host positions into Python dicts for a plot that draws
  only covered points. Fixed by dropping `-a` (`evidence.py:861`, test in `test_evidence.py`).
  Evidence now runs on its own pin (`vs_pinned/<fixsha>`); quant arms stay on `bbf1821`.
  Open: rerun the F-025 gate at 64 G to read real MaxRSS and BLAST output; BLAST against a
  prebuilt host DB instead of rebuilding per call.
- [x] `DSR-15` — twostep cell calling and denominators. The kb matrix of a `--host-filter
  starsolo` run holds no host UMIs, so emptyDrops failed on 7 of 8 twostep cells (1-2,659
  barcodes) and, where it ran, `n_called_cells` counted only barcodes with non-host reads
  (EBV twostep: 263/263 = 100 % infected vs 2,720 called cells in `combined_off`). Fix:
  emptyDrops on STARsolo's host matrix (`Solo.out/Gene/raw`, `emptydrops.R` now accepts
  genes x barcodes); `n_called_cells` / `n_comparable_cells` come from that set;
  `results/host_called_cells.tsv` records it (`cellcalling.py`, `detection.py`). Check on
  EBV twostep via `rerun-multimap`: 2,764 cells, EBV 91.2 % (combined_off 93.4 %), molecules
  unchanged. Guarded on `host_index`, so combined arms are untouched. Open: redetect the
  other twostep cells; pin the commit for the twostep arm.
- [x] `DSR-16` — one reference cell set per sample for every arm (user, 2026-10-06). Arms call
  cells on different matrices, so their own `pct_infected_called` are not comparable. The
  reference is the `combined_off` emptyDrops set (full kb matrix, no filter);
  `scripts/dsr_common_cells.py` counts each arm's `per_cell_viral.tsv` over exactly those
  barcodes -> `dsr_round1/common_cells_summary.tsv` (matrix cells, no rerun). EBV: 2,720 cells in
  every arm; EBV 93.4 / 93.4 / 94.3 % (off / artefact / twostep). Report this table, not the
  per-arm `pct_infected_called`. `scripts/slurm_dsr_twostep_redetect.sh` (DSR-15) writes `twostep_v2`.
  New datasets: run `combined_off` first.
  Vendor cells (user, 2026-10-06): `--reference-cells DATASET/SAMPLE=<cellranger barcodes>`
  overrides the emptyDrops set for that sample. Use it for `sfl_tonsil/x223` (cellranger
  `sample_filtered_feature_bc_matrix/barcodes.tsv.gz`, 50,000 cells). Atlas ERR13027027 has no
  local cellranger output; its twostep v1 failed in emptyDrops (18 barcodes) and is redetected
  as 25715969.
  Redetect is detection-only (user, 2026-10-06): the first twostep_v2 runs (25715231, 25715437)
  reran STARsolo + kb, because snakemake's default `code` trigger fired on the older quant
  commit, then the kb `mv` collided with the copied `kb-python/counts_unfiltered`. The script now
  runs from a profile with `rerun-triggers: mtime` (dry-run: only detection + umap), and "done" =
  `run_complete.json`, since the copy inherits v1's results. `dsr_common_cells.py` skips an
  unfinished twostep_v2 and drops v1 `twostep` once v2 is complete. Resubmitted as 25720116.
  Result: 10 rows done in 3-6 min each. Twostep agrees with combined_off except hpv16 SRR19537339,
  where twostep_v2 calls HPV77 (263 cells); evidence (25720142) shows these are human chr1 reads that
  STAR missed, no viral reads -> F-028 (`.living/findings/twostep-hpv77-call-is-host-reads-star-missed.md`).
  - 2026-10-13 (round-2 tooling, branch `codex/dsr16-round2-tooling`; software only, no job submitted):
    `scripts/explore_ledger.py` (hash-chained append-only TSV, flock-serialised, `check` fails on an edited or
    removed row or a Track B accession from `config/track_b_denylist.txt`; a removed *tail* needs `--expect-head`);
    round scripts take `ROUND_DIR` / `CODE_DIR` / `REF_BUILD` (defaults = round 1, `tests/test_dsr_round_env.py`),
    append to `$ROUND_DIR/explore_ledger.tsv` after each run outside round 1 and refuse a deny-listed manifest row;
    `docs/plans/dsr_round2_GUARD.md` is the GUARD template; `scripts/dsr_round_compare.py` writes the round 2 vs
    round 1 TSV + markdown (round-1 `twostep_v2` is the `twostep` baseline, `detected` filter). DEF-05 two-step
    cell calling is automatic with `--host-filter starsolo --cell-calling emptydrops`, so the `twostep` arm needs no
    new flag and round 2 has no `twostep_v2` redetect pass. Decision: most round-1 `viral_summary.tsv` files predate
    the `detected` column, so the comparison reads them as all-detected and says so in its header; a `disappeared`
    call can then be a round-1 candidate-only virus. Open: pin `<SHA>`, copy chem JSONs and the manifest minus
    `ERR13027122` / `GSM4682312`, then launch.
  - 2026-10-13 resume review (same tooling branch): recovered Claude's completed validation receipt at
    `66d8d05`: 3,542 passed, 13 skipped, 135 deselected (416.85 s); ruff check/format, mypy (53 files),
    and data governance passed. Review corrections preserve zero-call runs in the comparison (blank `virus`,
    `no_calls` when both rounds completed without detections); TSV separators are refused before ledger writes,
    truncated final lines are refused, missing-ledger checks fail, and checks hold a shared file lock.
    Current focused validation: 30 passed, 1 skipped (shellcheck unavailable), including common-cell callers;
    ruff check/format, mypy (53 files), data governance and `git diff --check` pass. The full suite was not
    repeated after these research-script-only fixes. No scheduler jobs, downloads, shared pins, pushes or PRs.
- [x] `DSR-17` — evidence replay selected reads with the wrong EC numbering (2026-10-07).
  `replay_exact_target_bus` re-runs `kallisto bus -n` (multithreaded) and then ran
  `bustools capture -e <primary matrix.ec>`; kallisto numbers ECs in discovery order, so the
  ids differ between runs (GSM5725695 combined_off: 93,708 to 224,345 differing `matrix.ec`
  lines; captured 80,608 records with the primary EC file, 66,671 with the replay's own, which
  matches the 66,684 target reads in the twostep_v2 evidence; a rerun gave 57,353). Rare
  targets returned mostly host reads. Fixed: capture and `parse_flagged_target_bus` use
  `lineage_bus/matrix.ec` (`evidence.replay_ec_path`, `evidence_run._replay_ec_map`, which
  also refuses a differing `transcripts.txt`); tests in `tests/test_evidence.py`.
  **Left:** regenerate affected evidence (all DSR-02 verdicts for small targets, the F-028
  HPV77 evidence, any F-025 evidence) on a pinned commit with this fix; verdicts in
  `dsr02_verdicts.py` outputs made before 2026-10-07 are suspect for low-abundance calls.
  Primary runs and `viral_summary.tsv` are unaffected.
- [x] `ALIGN-PLOT-01` — read-by-read host-versus-virus alignment figure (user, 2026-10-14: "is there any
  visualization that we can make to demonstrate alignment of reads in human vs virus?").
  `src/viralscan/competitive_plot.py`: `plot_competitive_alignment(evidence_dir, out_dir=None, max_lanes=200)` writes
  `competitive_alignment.png/.pdf` (A: one lane per read, virus row + host row, clipped, poly-A/T >= 10 nt hatched,
  10x TSO outlined, verdict and both `AS` at the right; B: virus-AS vs host-AS scatter; C: human loci in 1 kb bins;
  more than 200 reads are shown at an even stride with a "showing N of M" note). `cross_call_loci(dirs, prefix)` pools
  several evidence dirs into `PREFIX.tsv` (host 1 kb / viral 100 bp bins; reads, calls, libraries) plus `PREFIX.png`.
  CLI: `viralscan evidence-plot --evidence-dir D [D ...] [-o OUT] [--cross-call PREFIX] [--max-lanes N]`; `viralscan
  evidence` also writes the figure (a plotting failure only logs a warning). Reads `samtools view` text like the rest of
  `evidence.py` (no new dependency; `pysam` is not a project dependency). Decisions: `read_verdicts.tsv.gz` is written by
  `scripts/dsr02_validate_call.py`, not by `viralscan evidence`, so the figure falls back to BAM-derived verdicts when it
  is absent; reads are drawn in original orientation (reverse-strand CIGARs flipped); `SA` tags add supplementary
  segments no record carries; poly-A/T runs tolerate one-base interruptions. Tests `tests/test_competitive_plot.py`;
  CLI table regenerated in `docs/cli_reference.md` (`docs/output_reference.md` left alone: PR #53 edits it).
- [~] `VERDICT-01` — evidence verdicts that can support per-dataset conclusions (user, 2026-10-07; plan
  `~/.claude/plans/be-as-unbiased-but-careful-whenever-shimmering-teacup.md`, scope Phases 0-2).
  Interim rule: every unexpected call is reported as **unverified** until the calibration gate passes; the
  earlier "HPV16 positives are host reads" reading is withdrawn as unproven (it used pre-DSR-17 evidence).
  Done: `dsr02_verdicts.py` read `complex_body_fraction` backwards (it is the share of reads WITH a templated
  body, `anello_align.is_complex_body`), so EBV (1.52 M clean reads) was `low_complexity`; now `< 0.5` is low
  complexity, reads-weighted over virus references, blank = not measured. Regression: EBV `combined_off`
  fraction 0.9999 -> `viral_best`; `tests/test_dsr02_verdicts.py`.
  Done (0.2): `src/viralscan/molecule_verdict.py` labels each (cell, UMI) `virus_best` / `host_best` / `tie` /
  `unaligned` from the best host and virus `AS` in `competitive_reads.raw.bam` (minimap2 `-ax sr` already keeps
  secondary hits within 80 % of the primary, so no new alignment pass); `counted` = a read with
  `assigned_weight > 0`. `viralscan evidence` writes `molecule_verdicts.tsv.gz` and
  `molecule_verdict_summary.tsv`; `dsr02_verdicts.py` decides on counted molecules (not the run-level
  `host_homology` flag) with a PROVISIONAL `MIN_VIRUS_BEST_FRACTION = 0.5`, ties count against the virus, older
  evidence dirs fall back to the old rule. Streams `samtools view` (KSHV has 24-37 M reads). Tests
  `tests/test_molecule_verdict.py`; suite 1,888 pass. Smoke on HSV-1 `combined_off` (pre-DSR-17 evidence, so
  not a calibration result): 31,078 of 32,759 molecules `virus_best`.
  Done 2026-10-12 (`codex/dsr-evidence-tooling`, software only): (a) reference-cell restriction.
  `dsr02_verdicts.py --reference-cells <round>/reference_cells` (and `verdict(dir, cells)`) tallies only counted
  molecules whose barcode is in the frozen list, from `molecule_verdicts.tsv.gz`; a new `molecule_scope` column
  says `restricted` / `all_cells` / `no_molecule_table` / `no_frozen_cells`, so an unrestricted tally is never
  passed off as a cell-level one (assumes evidence CBs and the frozen barcodes use the same format). (b) the
  calibration gate tool `scripts/dsr02_calibration.py <evidence_root> <labels.tsv> --out calibration.json`: the
  user's labels (positive / artefact per evidence dir) in, per-row pass/fail and a record out; `gate_passed`
  needs every row right, at least one of each class, and no missing directory; a record is void when
  `MIN_VIRUS_BEST_FRACTION` or the new named `MIN_COMPLEX_BODY_FRACTION` (0.5, the old literal) changes. It picks
  no labels and no thresholds. Tests `tests/test_dsr02_verdicts.py`, `tests/test_dsr02_calibration.py`.
  **Left (user / cluster):** (0.3) regenerate evidence on a new pin (DSR-17 left-over; the evidence needs
  `molecule_verdicts.tsv.gz`); (0.4) choose the labelled set (EBV, HSV-1, KSHV, HHV-6B, HPV16 positives; F-025,
  F-028, F-019 artefacts), freeze the two rule values, run the gate, before any unexpected call is read.
- [x] `CELLS-01` — one emptyDrops call per sample, used for every cell-level analysis (user, 2026-10-07).
  **Amended 2026-10-14 (Q153, `CELLS-02`):** the single call is host-only: STARsolo `GeneFull` via two-step,
  or the host gene columns of the kb matrix in combined mode, never viral reads. The round-1 host+virus
  `combined_off` sets described below are diagnostics, superseded.
  The call is each sample's `combined_off` emptyDrops result (kb host+virus matrix; `--cell-calling emptydrops`,
  FDR 0.01, lower 100, 10,000 iterations, seed 100), already run, so nothing was recomputed.
  `scripts/dsr06_reference_cells.py <round>` freezes it to `<round>/reference_cells/<ds>__<sample>.tsv`
  (+ `index.tsv`: source hash, parameters, barcodes tested, knee, inflection, code SHA; never overwritten unless
  identical or `--force`) and writes `<round>/cell_level/<ds>__<sample>.tsv.gz` (every operative arm, frozen cells,
  viral molecules > 0). `dsr_common_cells.py` now reads the frozen list first and no longer drops samples whose
  run root differs from the inner directory (x223 had been missing from every table). 17 samples frozen
  (n/a: ERR13027122, GSM4682312). The arms' own cell calls (twostep host matrix) stay as diagnostics only.
  **Supersedes the 2026-10-06 vendor-barcode decision for x223** (our 44,936 cells, not Cell Ranger's 50,000;
  `--reference-cells` still accepts a vendor set explicitly). `common_cells_summary.tsv` rebuilt; the previous
  file is kept as `common_cells_summary.pre_frozen_2026-10-07.tsv`. Tests `tests/test_dsr06_reference_cells.py`.
  Not checked: that emptyDrops is well calibrated per sample (hhv6b knee == inflection; barcode-rank plots).
- [x] `CELLS-02` — host-only cell calling in combined mode (user, 2026-10-14, Q153; plan steps B1 + A4,
  branch `codex/host-only-cell-calling`). `cellcalling._host_only_emptydrops` drops the viral columns
  (`viral_gene_ids`: Virus Identity table, else `log/analysis.txt`) from the kb `cells_x_genes` matrix into a
  scratch copy, runs the unchanged `emptydrops.R` on it, and writes the result to the same
  `kb-python/counts_unfiltered/emptydrops_cells.tsv` (barcodes stay aligned; `total` is now host UMIs, knee and
  inflection are host-only). Fails closed (`CellCallingError`): no partition, partition matches no column, or no
  host column left. New receipt `results/cell_calling_input.json` (`count_layer: kb.host_genes`) in combined mode.
  Two-step (STARsolo `GeneFull`) and `--cell-calling external` unchanged. Docs: `faq.md`, `output_reference.md`,
  `cli_reference.md` (and its `--cell-calling` help), `quickstart.md`. Decision: REF-06 caveat recorded, the
  combined-mode host columns depend on the reference's D-list (cat42b to cat42d moved calls by -22 % on the
  HPV16 rafts); two-step `GeneFull` and an external list are reference-independent. Round-1 frozen
  `reference_cells/` (host+virus) are not regenerated here. Tests `tests/test_cellcalling.py::TestHostOnlyCellCalling`.
  G3 review must cover this change.
- [x] `CELLS-03` — host-only cell anchors, one per library (user, 2026-10-14: emptyDrops on the host transcriptome only,
  one call per library, every analysis of a dataset uses the same cells). `dsr06_reference_cells.py` gains
  `--source host_called` (round 2 `twostep` `results/host_called_cells.tsv`; knee / inflection / barcodes tested from the
  host-only STARsolo GeneFull `emptydrops_cells.tsv`; method `emptydrops-host-only-starsolo-genefull`), `--out-dir` (lists,
  `index.tsv` and `cell_level/` outside any round, whose GUARD allows writes only under runs/, status/, logs/) and
  `--provisional DATASET` (index `status` column). `dsr_common_cells.py --reference-cells-dir DIR`; precedence explicit
  override > dir > round's frozen list > combined_off called cells. 17 anchors frozen to
  `viralscan_work/cell_anchors/` (outside the repo; x213/x216/x223 are institutional): ebv 2,781; hhv6b 482
  (PROVISIONAL: library type undecided); hsv1 4,400; hpv16 9,057 / 7,055; kshv 5,587 / 7,748 / 6,597 / 8,940; hnscc
  5,219 / 1,172; tonsil_atlas 10,884 / 12,455; kshv_ebv 11,487; covid 17,810 / 10,333; x223 43,284. Each equals the
  `is_cell` rows of its emptyDrops table. Round 1 and 2 common-cell summaries rebuilt over the anchors in
  `cell_anchors/summaries/` with a diff against the old per-round summaries (round 2 had none: its before is
  combined_off). `twostep` added to the `cell_level` arms (`arm_is_final` already drops a superseded one). The
  anchors are not yet wired into the dossier or the Snakemake pipeline. Tests `tests/test_dsr06_reference_cells.py`,
  `tests/test_dsr_common_cells.py`.
- [x] `CELLS-04` — host-shutoff diagnostic (Q17): `scripts/dsr_host_shutoff_diagnostic.py <round> --anchors DIR --out-dir DIR`
  lists barcodes called on host+virus (`combined_off`) but not in the anchor, ranked by expected-virus molecules, and
  anchor cells absent from the host+virus call; also checks the barcode format of both lists (no mismatch in any of the
  17 libraries). Flags HSV-1, EBV, KSHV (lytic libraries lose hundreds of >=10-molecule cells; GSM5725696 and -98 about
  19% of called-cell KSHV molecules) and HHV-6B (provisional). Whether the dropped cells are real infected cells is not
  decided. Test `tests/test_dsr_host_shutoff_diagnostic.py`.
- [x] `AUD3-09a` — four small input-robustness fixes from audit 3 (P3-06, P4-04, P4-06, P4-03); tests `tests/test_aud3_09a_input_robustness.py`.
  (1) `read_support.molecule_support_rows` skips a molecule whose viral-aligned primary has `SEQ=*` and no queries FASTA, one aggregated WARNING with the count (was `IndexError`); unmapped `SEQ=*` still gets `sequence_unavailable`.
  (2) `whitelist_preflight.load_whitelist` strips a terminal `-[1-9][0-9]*` GEM-group suffix. (3) `chemistry.cb_umi_geometry` raises `ChemistryError` for an explicit `bc:umi:seq` whose barcode/UMI is not on file 0 with CB at base 0 and UMI directly after; `menu` dies on it when `--host-filter`/`--read-filter`/`--auto-evidence` is set. Non-obvious: `chemistry._same` (the `-x` vs detected check) uses a new lengths-only helper, so plain `kb count` with a file-1 geometry still works; the `cli_reference.md` `--technology` row is generated from parser help so it was left alone (documented in `faq.md`). (4) `knee_cells` docstring only.
  Pinned docs touched (re-pin pending): `docs/faq.md` (`faq-doc`), `docs/output_reference.md` (`output-reference-doc`).
- [x] `SE-01` — super-expressor summary per library (user, 2026-10-15: "for each dataset also do the superexpressor
  analysis, like the HHV-6 study"). `scripts/super_expressors.py --per-cell-viral F --cells ANCHOR` follows Lareau et al.
  2023 (Nature): super-expressor = a cell with >=10 viral UMIs (`--min-umi`). Per virus: number and percentage of
  super-expressor cells over the frozen cell set, share of the virus's UMIs in them, maximum viral UMIs per 10,000,
  and cell counts at >=3/5/10/50 UMIs. Exploratory, no detection claim. Interim table over the 17 round-2 libraries
  and their anchors written outside the repo (`cell_anchors/summaries/super_expressors_round2.tsv`). Not yet done:
  viral gene breadth and immediate-early/early/late split (needs the per-gene viral matrix), host-gene association.
  Test `tests/test_super_expressors.py`.
- [x] `CELLS-05` — fixes from the independent review of PR #64 (`CELLS-02`; branch `codex/cell-calling-review-fixes`, Q154).
  (1) BLOCKER: `seeds.cell_calling` 20260727002 exceeds R's integer max, so `emptydrops.R` got NA from
  `as.integer()` and `set.seed` stopped (combined and two-step). Fixed in code, protocol and schema untouched (the G3
  allowlist forbids editing seeds): the seed is parsed as numeric; a seed within R's integer range is used unchanged
  (seed 100 reproduces the pre-fix result exactly), a larger one is reduced modulo 2147483647 (20260727002 -> 933374179).
  `cellcalling.effective_seed` mirrors it; the script logs `seed=... effective_seed=...`; both receipts record `seed`
  and `effective_seed`. (2) Test gaps: the Virus Identity table is tested to win over `analysis.txt` (mutation
  "ignore the table" now fails two tests); a matrix column the table does not list raises `CellCallingError` (real
  round-2 tables cover all 47,259 columns); matrix/gene-list shape mismatch raises; real-R tests run `emptydrops.R` at
  seed 20260727002 and 100 on a tiny matrix. (3) `_run_emptydrops` unlinks the old `emptydrops_cells.tsv` first; a
  receipt is cleared at the start of every cell-calling call, so none/knee/combined-external runs and failed re-runs
  leave none. (4) Error messages name the remedy (re-run the analysis step, `--called-cells-file`, two-step
  `--host-filter starsolo`); the virus-only index message names two-step; CHANGELOG line. (5) Receipt gains
  `partition_source`, `partition_sha256` (combined), `r_version`, `dropletutils_version`, `code_sha` (null outside
  git), `build_sha256`; existing keys unchanged. (6) Stale `.host_only_*` scratch dirs are swept before a combined
  call; `output_reference.md` says the dossier's `median_umi_called`/knee/inflection are host-only now and
  pre-`CELLS-02` called sets must not be reused as anchors. Not changed: protocol, schema, scoring, the two-step
  and external logic.
- [x] `DOSSIER-01` / `DOSSIER-02` — per-dataset evidence dossier (user, 2026-10-07: list the files, tables and
  results needed to support each finding before running more analyses). Spec
  `docs/plans/2026-10-07-dsr-dossier-spec.md` (`analysis/dsr_round1/dossier_spec.tsv`, `dossier_roles.tsv`);
  read-only checker `scripts/dsr05_dossier.py <round>` writes `<round>/dossier/` (`index.tsv`, one Markdown
  checklist per sample, `round_index.md`, derived `cell_calling_summary.tsv`, `recurrence.tsv`,
  `arm_concordance.tsv`). Claim ladder C0-C3, N1; calls are enumerated from the runs. Tests
  `tests/test_dsr05_dossier.py`.
  - 2026-10-12 (`codex/dsr-evidence-tooling`): the new producers change the spec. E08 (evidence code SHA) now
    names `dsr02_validate_call.py`, which writes `code_sha.txt` (it said "none"); new non-gating items E10
    (`read_verdict_summary.tsv`, DSR-02) and K01 (`cards/evidence_cards.tsv`, CARD-01); 46 -> 48 items. The
    dossier already enumerated calls through `dsr02_enumerate_calls`, which skips `detected = false` (AUD-21), so
    no filter change was needed. Tests: `test_every_producer_script_named_in_the_spec_exists`,
    `test_new_producer_outputs_turn_their_items_ok`. Still no producer: A04 (deposit md5), A07 (STAR index
    fingerprint), J02 (claim map), H01 (labels); the ladder rungs are unchanged (user-owned).
- [x] `DOSSIER-03` — first round report on `dsr_round1`: 2,593 rows, 1,480 ok / 887 missing / 194 n_a /
  17 stale / 12 failed; every sample is C1, no informative negative. Gaps are listed in the spec (evidence
  molecule verdicts and verdict rows missing for all 140 calls, 53 calls without evidence, stale
  `common_cells_summary.tsv`, no positive controls or planted recovery, no labels or claim map, twostep host
  cell sets disagree with `combined_off`). Nothing regenerated yet.
- [~] `CARD-01` — `scripts/dsr03_evidence_cards.py`: per (dataset, virus) tiered evidence card (read-level,
  genome structure, ambient Poisson null, cross-dataset recurrence, cell-type permutation); needs VERDICT-01.
  - 2026-10-12 (`codex/dsr-evidence-tooling`, software only): `python scripts/dsr03_evidence_cards.py <round>
    [--calibration cal.json] [--labels DS/SAMPLE=tsv]` writes `<round>/cards/evidence_cards.tsv` + `.md`, one card
    per call (`dsr02_enumerate_calls.calls`, so `detected = false` rows are skipped) on `combined_off`. Tiers:
    read level (verdict over the frozen reference cells, `pass` only if `molecule_scope` = restricted; repeat-class
    read counts from DSR-02), genome structure (breadth, read-start concentration, hotspot flags; numbers, no
    cutoff), ambient (Poisson upper-tail p from the viral share of barcodes outside the frozen cells in
    `adata.h5ad`; conservative, p only), recurrence (`dossier/recurrence.tsv`), cell-type permutation (label
    shuffle, Bonferroni over labels; `n_a` without `--labels`). `card_status` enforces the interim rule: an
    unexpected virus is `unverified` until a calibration record passes under the current rule values. Tests
    `tests/test_dsr03_evidence_cards.py` (tiny synthetic round, incl. unexpected-but-read-level-pass, off-list molecules).
  - Decisions (conservative): card per (dataset, sample, virus), not per dataset (a multi-sample rule is the
    user's); ambient uses unique-assignment `adata.h5ad` counts, not multimap-corrected estimates; labels file is
    `barcode<TAB>label` (H01's format was unspecified); `--n-perm` 1000 / seed 0 are computational defaults.
  - **Left:** run on a round after the VERDICT-01 0.3 regeneration and a calibration record (user); cell-type
    labels (H01); any pass/fail cutoff on the genome-structure or ambient tiers (user decision).
- [~] `ADJ-01` — generic candidate adjudicator (study Track B / M5): deterministic claim level L0-L3 per (library,
  virus), mapped onto the DSR-02 verdict classes and the CARD-01 reference arm; L4/L5 are never emitted.
  - 2026-10-13 (`codex/adj01-adjudicator`, software only; no real round was run):
    `python scripts/adjudicate_candidates.py <round_dir> [--lock L.json] [--metadata samples.tsv]
    [--germline-flags F.tsv]` writes `<round_dir>/adjudication/candidate_levels.tsv` (dossier item K02, non-gating).
    Candidates = `dsr02_enumerate_calls.calls()` on `combined_off`; `library_id` = `<dataset>__<sample>`.
    Columns (stable; append only): `library_id, virus_key, level, level_reason, flags, unique_molecules, n_loci,
    read_virus_best_fraction, thresholds_sha256, subject_id, timepoint, run_id, lane` (LONG-01 pass-through from the
    metadata TSV), then `virus_name, allocated_molecules, positive_cells, max_sibling_allocated_fraction,
    complex_body_fraction, read_verdict, molecule_scope, gene_programme_cells`, then (2026-10-13, `TIER-01`/Q104)
    `package_tier_analogue` (L0 candidate, L1 probable, L2 strong, L3 strong).
  - Rules: unique layer only (`multimap_evidence.tsv`); allocated mass alone is L0 `allocated_only`; L1 needs
    read-level evidence over the frozen reference cells (`molecule_scope` restricted), virus_best fraction and
    complex-body fraction at the lock's minima, no sibling bleed (allocated / (unique + allocated) when the virus has an
    identity-table `sibling_group`); L2 adds the distributed route (unique molecules, positive cells, loci); L3 needs a
    `gene_program_summary.tsv` model (the nine herpesviruses) and `l3_min_programme_cells`, else `cap:no_gene_programme`;
    Anelloviridae (`claim_scope` screening_only, genus fallback without an identity table) cap at L1 and carry
    `screening_only`. Flags: `screening_only`, `sibling_bleed`, `germline_integration_compatible` (passed through from
    `--germline-flags` or a `viral_summary.tsv` column; nothing detects it), repeat classes `tso`/`telomeric`/
    `low_complexity` (labels only), `test_only_thresholds`. Host expression and `hostresponse` outputs are never read.
  - Decision: the package ships NO number. Thresholds come from a lock JSON (keys named as in the study
    `PROTOCOL_TEMPLATE.yaml` `adjudication:` block; plus `min_virus_best_fraction`, `min_complex_body_fraction`,
    `l3_min_programme_cells`); no lock = every row L0 `no_protocol_lock`, an incomplete lock = L0
    `incomplete_protocol_lock:<key>`. The DSR-02 class is shown as `read_verdict` but does not decide the level (its
    0.5 provisionals would be default numbers). Test fixture `tests/data/adjudication_lock_TEST_ONLY.json`.
  - **Left:** LOCK-01 supplies the numbers (2026-10-14, `LOCK01-PRE`: `n_loci` and `positive_cells` now follow LOCK-01 D1/D2, no longer placeholders); the rare-super-expressor route (deferred, Q124)
    needs a per-cell unique count nothing writes yet; gene-programme catalogue names ("Epstein-Barr virus") differ from
    identity-table names, so a name mismatch reads as no model (fails to L2, visible reason); `max_sibling_allocated_fraction`
    is allocated share, not per-sibling attribution; ADJ-02 host-homology screen is separate. Run on a real round.
  - 2026-10-13 (Q58, `codex/trackb-gates-q53-q70`): **L3 eligibility reads a taxid map.** `scripts/adjudicate_candidates.py::PROGRAMME_NAME_BY_KEY` maps the
    nine gene-programme herpesviruses' `virus_key` to their `gene_programs.tsv` / `gene_program_summary.tsv` name, so an identity-table run ("Human gammaherpesvirus 4",
    `taxid:10376`) reaches L3 without the legacy name matching. Taxids from `virus_catalog.tsv`: EBV 10376, HCMV 10359, HSV-1 10298, HSV-2 10310, HHV-6A 32603,
    HHV-6B 32604, HHV-7 10372, KSHV 37296, VZV 10335 (the six in `evidence.SELECTOR_HANDLES` agree). The exact-name match is tried first, so legacy runs are
    unchanged; anything off the map (or a mapped virus with no summary row) still reads `cap:no_gene_programme`. `test_taxid_map_agrees_with_the_catalogues` pins the
    map against both catalogues. This resolves the name-mismatch item under "Left" above.
- [x] `LOCK01-PRE` — **new, 2026-10-14 (Q123, Q128).** Package prerequisites for the LOCK-01 pre-registration (study repo `protocol/LOCK-01_preregistration.md`, U1 and U6). Branch `codex/lock01-prereqs`.
  - **D1 locus** (`src/viralscan/viral_loci.py`): loci are the connected components of annotated gene intervals that share at least one base on either strand (closed intervals; unstranded, DEC-04). All rows of one `gene_id` are one node, so an origin-wrapping gene (two `part` rows) merges across the origin; a gene with no `gene` row takes its CDS/exon span; different sequences never overlap. The `analysis` rule writes `results/viral_loci.tsv` (`gene_id, seqname, start, end, locus_id`; `locus_id` = smallest `gene_id` of the component) from the build manifest's viral GTF (`provenance.viral_gtf`, else `gtf`), or the bundled panel when there is no manifest, plus the custom `--gtf` files. No annotation found = no file.
  - **D2 positive cell**: `summarize_multimap_evidence(..., called_mask, locus_by_gene)` appends `locus_id` and `n_unique_positive_cells_called` to `multimap_evidence.tsv` (append-only; documented in `docs/output_reference.md`). The latter is per virus, repeated on each of its rows: called barcodes (the run's frozen cell-caller mask, `detection` passes it) with at least one unique-layer molecule at any gene of the virus; the union over genes, not the maximum. `n_unique_viral_cells` (per gene, all barcodes) is unchanged.
  - **Adjudicator**: `unique_layer()` now reports `n_loci` = distinct `locus_id` over genes with a unique molecule, and `positive_cells` = `n_unique_positive_cells_called`. Both placeholders (`ponytail`) are gone. A run that predates the columns (or has an unplaced gene with unique signal) gets `None`, so the candidate stops at L1 with `missing:min_independent_loci` / `missing:min_positive_cells`; nothing is guessed. Calibration needs runs made after this change.
  - Tests (TEST_ONLY fixtures only, no real data, LOCK-01 rule A5): `tests/test_viral_loci.py` (+ `tests/data/viral_loci_TEST_ONLY.gtf`), `tests/test_adjudicate_candidates.py`.
  - **U6 / Q128** lives in the study repo: `registries/package_holdout.tsv` (it is not a package file) gets a `consumable` column and a validator; see the study repo commit on `master`.
  - **Decisions Q123-Q133 (user, 2026-10-14):** Q123 package D1/D2 before calibration (this row). Q124 rare-super-expressor route deferred; every rare-cell finding is capped at L1. Q125 a failing LOSO system gets a `calibration_unstable:<system>` flag (L1 cap for its virus) instead of blocking. Q126 procedural constants set in LOCK-01 section 10. Q127 V04 CHB = P, FC = U; V03 HA-sorted is report-only. Q128 `consumable` column. Q129 absolute counts now; depth-normalised thresholds in a later protocol version. Q130 HHV-6A/6B flagged by hand, L4 blocked for them. Q131 the CAL-08 10k 3' and 5' sets are one donor unit. Q132 tonsil-atlas remainder deferred. Q133 agents sign commits with their own model's Co-Authored-By trailer. Also: the download disk cap was raised to 95%.
- [~] `SENS-07` — per-dataset planted-read sensitivity (exact + held-out divergent plants in each dataset's own
  geometry, per-family LOD, measured capture into `--positive-control-*`); exploratory, not VAL-01 (G3 gates it).
  **Software done 2026-10-12 (`codex/sens07-dlist02`); no real dataset was run.** `scripts/sens07_planted_sensitivity.py`
  `plant` / `score`, built on `plant_anello_10x.py` (not a second planter). The planter gained optional
  `read_len`, `umi_len`, `levels`, `windows`, `error_rate`, `barcode_blocks`, `taken` (background `(CB, UMI)` pairs a
  plant must not reuse, as in D1) and `circular=False` (linear genomes do not wrap); defaults are unchanged.
  `plant`: geometry from `--technology` plus the background's own R2 length (R1 checked against barcode + UMI);
  cells = most frequent barcodes (outcome-independent); one disjoint cell block per (genome, level), so one run
  separates conditions; exact genomes (`--exact-fasta`) and held-out genomes (`--heldout-fasta`, ID-disjoint from the
  panel, `--paf` gives the nearest-panel identity band); planted and unplanted twin FASTQs; `truth.tsv`,
  `selection.tsv`, `plant_manifest.json` (input digests, labelled exploratory, `thresholds: null`).
  `--control-genes` (with `--gene-spans`) plants a gene body as an exact control, drops its parent genome from the
  other plants (so the gene count is only the control), and with `--index --reference-manifest --identity` writes a
  checksum-bound receipt and the `--positive-control-*` argument line (`positive_control_args.tsv`);
  `sensitivity.verify_exact_control` accepts it (test). `score`: per-condition recovery from an h5ad layer (the same
  layer `--positive-control-*` measures) or a kb counts dir, minus the unplanted twin; per (family, mode) capture
  with Wilson interval, `failed` / `over_recovered` / `measured` (SENS-04's states),
  `lod95_true_molecules = minimum_molecules_for_detection(0.95, t) / capture`, lowest level detected in every
  genome. Tests `tests/test_sens07_planted_sensitivity.py`.
  - Decisions: held-out means a caller-supplied FASTA (not mutation: that is `SENS-CAL-01`); `--per-family N`
    picks from the sorted accessions with the seed and invents no selection rule; no pass/fail threshold anywhere
    (`--threshold` is only the molecules that count a condition as detected, default 1). A control receipt needs
    `--count-layer`/`--count-strategy` to match the consuming run.
  - Integration review repair 2026-10-13: scoring accepts only `--group-column family`, matching the
    planting truth's family labels; other groups previously yielded false zero capture. Exact-control
    receipts now compare the actual planted parent sequence's SHA-256 and length with the reference
    manifest before certification; an accession label alone is insufficient. Regression tests cover
    mismatched sequences and unsupported group choices. Real-dataset measurement remains pending.
  - **Left:** the real per-dataset runs (planted and twin through `viralscan`, then `score`); the held-out
    partition per family under `VAL-RA-DECISION-01` (whole-taxon exclusion, ANI < 95 %, two sister taxa kept) and
    the user's choice of families and genomes per dataset; a probit LOD (E8, at least seven levels over several
    libraries; this is n = 1); 3' poly-A end reads (only internal fragments are planted); strand (fixed 2026-10-14 below).
  - 2026-10-14 **planter strand fix (`codex/planter-strand`, Q117):** the round 2/2b planter drew forward genome
    slices, but kb counts stranded, so reads of minus-strand genes were lost (42 % of planted herpes reads) on
    3' (`forward`) libraries, and on `--strand reverse` libraries (covid x213/x216, sfl x223) the plus-strand genes
    were lost instead. **Round 2b numbers are lower bounds** (also: genome-uniform plants count intergenic reads
    against capture). Fix: `plant --transcripts-fasta cdna.fa` (the reference's own sense-strand transcripts, so
    every plant is indexed) draws a transcript in proportion to its length, then a uniform start; `--library-strand
    {forward,reverse}` (default forward) plants the reverse complement for a reverse library; `truth.tsv` gains a
    `transcript` column. Without `--transcripts-fasta` the legacy forward-slice planter runs with a warning.
    `read_gene_spans` kept no strand, so a minus-strand control gene was planted antisense: it now plants the
    gene's sense strand. `slurm_dsr_round2b_sens07_array.sh` passes `--transcripts-fasta $B/cdna.fa` and the row's
    strand. `plant_reads` gained optional per-genome `transcripts`/`flip` (anello callers unchanged; the anello
    5'/3' windows are not touched). Tests: minus-strand transcript reads are antisense to the genome (forward
    library) and in the genome (reverse library). Round 2c waits for a fresh pin after round 2b.
  - 2026-10-13 **submitted** after the user approved Q45-Q52: job 25755522 (`--array=0-16%2`, 8 cpu, 64G, 48 h),
    driver `vs_pinned/bb7a34f`, code `vs_pinned/35de249` (Q49), frozen config and filled GUARD in `dsr_round2b/`;
    exact plants, default families + expected family, `--per-family 2`, seed 20261012 (Q48), first 50M pairs
    (Q50), both herpes labels (Q51), blank-family rows plant the default families only (Q52); rows 1
    (SRR20710641) and 15 (GSM4682311) n/a in `status/n_a.tsv`.
  - (history) round 2b scripts ready (`codex/sens07-round2b`), then awaiting the user's Q48/Q49
    answers and the go to submit. `scripts/slurm_dsr_round2b_sens07_array.sh` (one array task per round-2
    manifest row: `plant` + twin, both through `viralscan` with round 2's `combined_off` flags, `score`; own
    ledger `dsr_round2b/explore_ledger.tsv`, label `explore`, arms `combined_off_planted/_unplanted`),
    `config/dsr_round2b_sens07.template.tsv` (per-row `run`, `expected_family`, `n_background`, `note`),
    `docs/plans/dsr_round2b_sens07_GUARD.md`, `tests/test_dsr_round2b_sens07.py`. Defaults = the Q48/Q49
    recommendations as parameters: exact plants only, `PLANT_FAMILIES` Orthoherpesviridae Herpesviridae
    Papillomaviridae Polyomaviridae (the identity table labels EBV `Herpesviridae` and the other human herpesviruses
    `Orthoherpesviridae`) plus the row's expected family, `PER_FAMILY=2`, `SEED=20261012`, code pin
    `vs_pinned/35de249` (deviation from round 2's `e7c1b08`, recorded in the GUARD and each row's `code_sha.txt`).
    Decisions not in Q48/Q49: `n_background` defaults to the first 50M pairs (the planter holds every (CB, UMI)
    in memory, so the 1.6B-pair libraries cannot be planted whole); rows `SRR20710641` (R1 150 nt) and
    `GSM4682311` (R1 28 nt, forced 10xv2) cannot match the plant geometry and are `run=no` (n/a) until decided;
    `hnscc_gse164690`, `tonsil_atlas`, `sfl_tonsil` have no recorded expected family (blank in the template).
    Caveat for the user: `select_genomes` draws all families from one RNG in sorted order, so a row that adds a
    family sorting before `Papillomaviridae` (e.g. `Coronaviridae` on the covid rows) plants different
    Papillomaviridae genomes than the other rows at the same seed; a per-family seed is a new selection rule, left
    to the user. A 2M-pair `plant` smoke on `SRR8315713` (scratchpad, not `viralscan_work/`) selected 7 genomes in
    21 conditions and took about 12 min under load (the manifest also hashes the full background FASTQs).
    Dry run only; nothing was created under `viralscan_work/` or `vs_pinned/`.
  - 2026-10-13 (Q63, Q64, user-confirmed): **round 2c** plants HHV-6B `SRR20710641` with R1 trimmed to barcode + UMI (the 150-nt R1 no longer blocks it); `GSM4682311`
    stays n/a. From round 2c on the genome selection uses a **per-family seed** (resolves the RNG-order caveat above); round 2b keeps its single seed. Implemented in the next bullet.
  - **Round 2c ready, not submitted** (`codex/sens07-round2c`, 2026-10-13; waits for round 2b job 25755522 to
    finish; round 2b, `dsr_round2b/` and its pins untouched). User-approved Q63/Q64. `plant` gained
    `--selection-seed-mode {global,per-family}` (default `global` = round 2b's one-RNG draw, unchanged;
    `per-family` seeds each family's draw from sha256(seed, family), so a family's genomes are identical in every
    dataset whatever other families it plants; fixes the caveat above) and `--trim-r1-to-geometry` (R1 longer than
    barcode + UMI is trimmed to it in BOTH the planted and unplanted twin, sequence and quality, R2 untouched;
    recorded in `plant_manifest.json` as `r1_trim` {original_length, trimmed_length} and `read_changes`, a declared
    read change; R1 shorter than barcode + UMI is still refused, trim or not; without the flag a longer R1 is
    refused as before). `planter.copy_background` gained optional `r1_len`. Driver: the round-2b array script is
    extended, not forked: `SELECTION_SEED_MODE=per-family` (env) and an optional per-row `trim_r1` config column
    (read by header name, so round 2b's 6-column config reads as `no`); the new flags go after `PLANT_EXTRA`, so
    round 2b's command line is byte-identical; the script refuses a `CODE_DIR` whose planter lacks a requested
    flag, and refuses any round-2c setting while `ROUND_DIR` is (or defaults to) `dsr_round2b`, so a forgotten
    `ROUND_DIR` cannot `rm -rf` into round 2b's tree. `config/dsr_round2c_sens07.template.tsv` (SRR20710641 hhv6b `run=yes trim_r1=yes`, all other rows as 2b),
    `docs/plans/dsr_round2c_sens07_GUARD.md` (pin `vs_pinned/PIN` to be filled after merge; Q63/Q64 changes).
    GSM4682311 (R1 28 nt vs 10xv2 26) is the longer, trimmable case by the same rule, but stays `run=no`/n/a per
    Q63 (FORCE_TECH row, chemistry unresolved). Tests in `tests/test_sens07_planted_sensitivity.py` (per-family picks
    independent of the other families; default draw unchanged; trim identical in both twins, short R1 refused) and
    `tests/test_dsr_round2b_sens07.py` (2c template, GUARD, dry-run flags). Dry run of row 1 prints
    `--n-background 50000000 --selection-seed-mode per-family --trim-r1-to-geometry`. **Next:** after round 2b
    ends and this merges, pin a clean commit, fill `PIN`, freeze `dsr_round2c/{config.tsv,GUARD.md}`, user submits.
- [x] `READS-01` — `viralscan.reads`: `extract_virus_reads(run_dir, viruses, out_dir)` returns, per virus, a
  `VirusReads` (`.fasta`, `.sequences()`, `.lineage()`), plus `read_fasta`, `align_reads`
  (= `align_reads_to_viral`) and `blast_reads` (= `blast_identity`) for BLAST / alignment (user,
  2026-10-06; kallisto virus_detection_sc "extract reads for specific virus IDs"). Wraps
  `run_evidence` without a reference, so it replays the exact molecules the call was counted from
  (no new extraction code). Failure raises `RuntimeError` instead of `sys.exit`.
- [~] `DLIST-01` — host genome + cDNA D-list (Pachter `7_virus_host_captured_dlist_cdna_dna`; user,
  2026-10-06). cat42d's virus-only index is D-listed against the GRCh38 genome only; the notebook uses
  host cDNA + DNA. I0 `scripts/dlist_kmer_check.py` (job 25720864): are the HPV77 false-positive k-mers
  (F-028) in the genome or cDNA D-list? I1 `scripts/slurm_cat42d_virus_only_dlist2_index.sh` (job
  25720865): index in new `viral_ref_cat42d_dlist2/`, frozen cat42d untouched. I2 (controls) after I1.
  I0 result: 6/6 shared 31-mers are in `dlist.fa`. I1 built (25720865, 23 min). 2026-10-07 single-read
  `kallisto bus` test of one HPV77 TSO+CAG read (200 copies): pseudoaligns to Y15175 E2/E4 (ECs 3670,3671)
  on BOTH cat42d and dlist2 indices. Segments: TSO alone 0, CAG alone 0, CAG+tail 0, TSO+CAG (58 bp) 100/100.
  So the match needs the TSO-junction `GGGGCAGCAG…` k-mers, which a D-list does not mask (D-list only stops
  extension at distinguishing flanking k-mers, it does not delete target k-mers also in the decoy).
  **A bigger D-list cannot fix F-028; I2 on dlist2 is not expected to help.** Fix is TSO handling in twostep
  (trim / require evidence), a user decision. dlist2 index kept, frozen cat42d untouched.
- [x] `TSO-01` — F-028 fix, both parts (user, 2026-10-07; DLIST-01 showed a D-list cannot mask the junction).
  (a) `--read-filter tso-trim` (`scripts/read_filter.py`): cuts the 30 nt 10x TSO (`TSO_FULL`) from the start of
  R2 (sequence and quality), drops only pairs left under 31 nt (`r2_short_after_trim`), lineage `tso_trimmed` /
  `untouched`, audit row `tso_trimmed`. Opt-in; `artefact` mode output is unchanged. Tests in
  `tests/test_read_filter.py` incl. the HPV77 junction; `docs/cli_reference.md` regenerated.
  (b) `scripts/dsr_common_cells.py --verdicts evidence/verdicts.tsv` adds `evidence_status`: twostep(_v2) calls
  are `verified` (`viral_best`), `rejected_<verdict>` or `needs_evidence`, from the verdict of this arm's
  evidence dir, else the pre-v2 twostep one, else combined_off's; other arms `not_gated`. A first version
  accepted any combined-arm call as corroboration and passed HPV77 (263 cells, SRR19537339); removed. Not used in DSR round 1 (quant stays on `bbf1821`, user 2026-10-06).
  **Not measured yet:** re-quant of a TSO-heavy sample (SRR19537339 twostep) with `tso-trim`; the
  trim's effect on genuine 5' viral reads is untested.
- [~] `DLIST-02` — host-capture split (`kallisto bus -n` + `bustools capture --complement`) as
  `--host-capture`, per-virus `reads_virus_only` / `reads_also_host`; only if DLIST-01 I2 shows a gain.
  **Software done 2026-10-12 (`codex/sens07-dlist02`); the gate is not met and nothing was run on real data.**
  Opt-in `--host-capture` (default `None`, in `run_safety.omit_when_unset`, so default fingerprints are
  unchanged) adds the `host_capture` rule after `kb_count`: `evidence.replay_host_capture_bus` replays the
  quantified FASTQs with `kallisto bus -n` (the bus + correct chain is now shared with `replay_exact_target_bus`
  through `_replay_quant_bus`; file names unchanged), captures every read whose class holds a viral transcript,
  then `bustools capture -s -x` with the host transcripts. `scripts/host_capture.py` writes
  `results/host_capture.tsv` (`virus_name, reads_total, reads_virus_only, reads_also_host,
  host_capture_status`) and `log/host_capture.done`. Reads, not molecules; a read compatible with two viruses
  counts for both; viruses with no read have no row. Tests `tests/test_host_capture.py` (real `bustools` on a
  hand-written BUS file pins that a class is captured when it holds ANY listed transcript, so a host+virus class
  is `reads_also_host`).
  - Decisions: (1) the row's recipe is read as a **combined host+virus index**, the Pachter-notebook form. On
    a virus-only index there is no host to split on, so every count is `not_applicable` (never 0, as for
    `host_viral_unresolved`). The option B(i) route in `docs/plans/competitive-host-virus-reference.md` (a second
    pass of extracted reads against a host-only index, for users on a virus-only prebuilt index) is **not**
    implemented. (2) A separate TSV, not `viral_summary.tsv` columns: `detection.py` and
    `docs/output_reference.md` belong to M03. (3) New file in `config/public_ship_scope.json` (wheel + sdist).
  - **Left:** the gate (DLIST-01 I2 not run, and the TSO finding says a D-list cannot fix F-028, so a gain is not
    expected: your call whether this ships on); a real-run check on a combined index (`kallisto bus` hangs on
    synthetic tiny inputs in this sandbox, so no end-to-end `kallisto` test exists); `docs/output_reference.md`
    rows for `results/host_capture.tsv` and `log/host_capture.done` (docs owner, pinned file); a
    `viral_summary.tsv` merge if wanted (`detection.py`, M03); B(i).
- [x] `EMC-01` — **done 2026-10-10 (code); the cluster re-measure is yours.** `rerun-multimap em-cell` memory: dense per-barcode theta
  (`multimapping.py:679`) peaks at 130.6 GB on 793 k barcodes (64 G: OOM at 23 min;
  480 G: 39.6 min, 5.28 alloc core-h). Sparse theta would remove the need for 480 G.
  Fix: `cell_theta` now stores `(compatible genes, fitted values)` per cell, bit-identical allocation. Not changed:
  the transient per-cell dense row inside `em_cell_abundances` (O(n_genes) per cell, freed each cell). Lower `--mem` in
  `scripts/slurm_rerun_multimap_em_cell_cost.sh` only after a real 793 k-barcode rerun confirms the peak.
  - **Cluster re-measure done 2026-10-13: PASSED.** No bar was recorded in this row, so the fallback was used
    (completes within 12 h and 128 G on the largest library). Library: round-2 `covid/combined_off/LUM-SJ-x213-g`
    (largest completed round-2 BUS, 16,986 MB; 942,700 barcodes, 47,259 genes, 8 cpus; code `88f02a1`; source dir
    read-only, copy under `viralscan_work/emc01_remeasure/`). `rerun-multimap em-cell` (job 25755609): 1:22:58 wall,
    `/usr/bin/time` max RSS 17.6 GB (sacct MaxRSS 16.2 GB), vs 130.6 GB before the fix. `em-global` (job 25755613):
    1:10:16, 16.7 GB (sacct 15.6 GB). Both runs include detection and UMAP. Molecule audit against the host-conservative
    source: matrix shape (942,700 x 47,259) equal, `molecule_audit` (input, resolved, unique, ambiguous, unresolved,
    allocated ambiguous mass 11,640,844) and layer totals (unique 125,417,062; allocated 11,640,844; combined 137,057,906)
    equal for both EM methods (only 18 host-virus ambiguous molecules, so no conservative-vs-EM difference in totals).
    The `--mem` in `slurm_rerun_multimap_em_cell_cost.sh` can now drop (not changed here; its 793 k run was not repeated).
- [x] `DSR-09` — row-2 downloads (user request 2026-10-05) into
  `$SITE_SCRATCH/ViralScan/benchmark_inputs/dsr_2026-10-05/`:
  SRA array 25701724 (GSE190558 all 16 lanes; GSE164690 SRR13419165 HN18 CD45−
  and SRR13418973 HN01 CD45−, one lane each) via `fasterq-dump
  --include-technical`; ENA array 25701725 with md5 check (GSE154900 4 runs;
  E-MTAB-13687 ERR13027010 / ERR13027122 / ERR13027027). Tonsil rule: every
  scRNA library is a hashed pool of ~100 GB, so the pilot takes the
  lowest-accession run of scRNA Samples 1, 2, 3 (metadata only). HSV-1 is
  already on disk. The GSE190558 lane manifest and the GSE154900 read-cap
  decision are still to be frozen before any run.
  - **Done 2026-10-06.** All 32 tasks COMPLETED with no error or md5 mismatch.
    18/18 SRA runs are `DONE`. SRR17180386 (checked on 100k reads) has I1 8 bp,
    R1 28 bp and R2 cDNA **55 bp**, so GSE190558 cDNA is short, not SRA-trimmed.
    14/14 ENA files are md5-ok. Per-run dirs are gitignored (`benchmark_inputs/dsr_*/*/`).

## WP7 — Score, audit, and freeze results

Objective: turn completed runs into a single immutable scientific result bundle.
Estimated effort: 3-5 days after all runs finish.

- [ ] `RES-01` — compute all preregistered molecule/cell metrics, runtime, peak
  RAM, temporary storage, output size, calibration, and failure rates.
- [ ] `RES-02` — estimate uncertainty by bootstrapping biological samples, never
  by treating cells as independent experimental replicates.
- [ ] `RES-03` — calibrate thresholds on training data, freeze them, then evaluate
  the untouched holdout once; record every deviation.
- [ ] `RES-04` — independently audit inputs, references, barcodes, features,
  denominators, layers, manifests, scheduler records, and failed rows.
- [ ] `RES-05` — freeze `analysis/v3_results_bundle/` with a manifest containing
  every file digest, generation command, environment/tool version, Git SHA, and
  validation report.

`G5` passes only with 100% count-invariant compliance, complete ambiguous/
unresolved reporting, no probable/strong call in synthetic host-only or planted
host-homology negatives, positive recovery with sample-level uncertainty, and no
general superiority claim based on one virus or sample.

## WP8 — Regenerate documentation and claim evidence

Objective: make public instructions and claims match the actual v3 CLI, schemas,
and frozen results. Estimated effort: 4-7 days; quantitative pages wait for `G5`.

### WP8A — User documentation

- [~] `DOC-01` — reconcile README, installation, quickstart, CLI, outputs, API,
  - 2026-10-04: the pip/full-tier text has landed (7d4bae0). Still open: Sphinx
    `-W` (needs a user-provided sphinx env) and the REF-11 data-fetch
    wording.
  - 2026-10-12 (`codex/docs-pass`): Sphinx `-W` PASSES (0 warnings, `-E -a` clean build) in the
    `viralscan_test_full` env (sphinx 9.0.4, myst-parser 5.1.0, nbsphinx 0.9.8, sphinx_rtd_theme
    3.1.0 are now installed there; no new env needed). `-W -b linkcheck` also passes after adding
    `linkcheck_ignore` in `docs/conf.py` for the changelog compare/tag URLs (404 until tags are
    pushed) and the ACP DOI (403 to bots), and fixing the codecov and snakemake redirects. Still open:
    the REF-11 data-fetch wording (user decision) and a CI job that runs the two sphinx commands
    (`.github/workflows/` is owned by another agent).
  FAQ, reference-panel, support, security, and migration docs with v3 contracts.
  - 2026-10-02 (`v3/docs`):
    - Added `docs/migration.md` (rebuild-only) to the toctree, and linked
      SUPPORT/SECURITY from the index.
    - Guard tests cover "2.5.0", faq, reference_panel and migration.
    - Removed the covid Alphatorquevirus counts from cli_reference, faq and
      output_reference, using the Q11 wording, and the "~4×" claim in
      `BENCHMARK_COMPARISON.md`.
    - Earlier, under DOC-02: the v2.5 container commands and the stale
      primary-call mode.
    - Still open:
      - Sphinx `-W`: nbsphinx isn't registered in the codex env and
        `sphinx_rtd_theme` is missing (installs not approved). A stripped
        build still has 37 pre-existing warnings: `csv` lexer, PLAN
        cross-refs in docs/plans, header jumps.
      - pip/full tier text (REL-01/02) and the data-fetch wording (REF-11).
- [x] `DOC-02` — generate CLI/default tables from the parser and test exact
  defaults; remove plain `em`, removed primary-call modes, silent knee fallback,
  v2.5 container commands, and raw-UMI language for fractional estimates.
  - Done 2026-10-02 (`v3/cli`, d2703eb).
    - `menu.build_parser()`; `scripts/gen_cli_reference.py [--check]` writes
      10 generated flag tables into `docs/cli_reference.md`.
    - `tests/test_cli_reference_parity.py` fails on any undocumented flag.
      The 11 that were missing are now covered.
    - The malformed 4-column rows are fixed.
    - Removed `--multimap-primary-call confidence` (showcase runbook) and the
      v2.5 container commands, and changed per-10k wording to molecule
      estimates.
    - A search found no "silent knee fallback" text, and plain `em` survives
      only in the notebooks, which belong to DOC-06.
- [x] `DOC-03` — clearly label every output as observation, model estimate,
  evidence tier, diagnostic flag, or biological interpretation.
  - Done 2026-10-02 (`v3/docs`, 33e78de).
    - Every column table in `output_reference.md` has a Kind column, with a
      legend of 5 labels: observation, model estimate, evidence tier,
      diagnostic flag, biological interpretation.
    - Four emitted viral_summary columns that were undocumented are added.
    - Tests check every label, plus viral_summary and
      `MULTIMAP_EVIDENCE_COLUMNS` coverage.
- [x] `DOC-04` — document combined versus two-step information loss, anellovirus
  - **Done 2026-10-04 (7d4bae0):** pip-tier vs full-tier install text in the README,
    installation.md and faq.md.
  screening limits, pip/full-workflow tiers, and legacy rebuild-only migration.
  - 2026-10-02 (`v3/docs`): FAQ sections on what combined vs two-step
    each lose, and on anellovirus screening limits (commensal prior, Kane
    et al.; F-019/F-021/F-022). Still open: the pip/full tier text, which
    needs REL-01.
- [~] `DOC-05` — execute clean-install quickstart commands and validate the
  resulting run using only documented steps.
  - 2026-10-12 (`codex/docs-pass`): executed in a fresh venv (`python -m venv`, no system
    site-packages; `pip install ".[full]"` from this checkout, `doctor --profile pip` exit 0, plus
    `pip install kb-python==0.30.2`; kallisto/bustools from the conda env on PATH) on the tiny fixture
    `tests/data/evidence_tiny`.
    - Quickstart Mode 2 (`--reference -fasta -gtf`) completed with `--cell-calling none`, and Mode 1
      completed reusing its `output/index/{index.idx,t2g.txt}`. All four documented key results exist
      (`viral_summary.tsv`, `per_cell_viral.tsv`, `report.html`, `adata_multimap.h5ad`), and
      `viralscan validate-run output/` returned `"ok": true`, `"issues": []`.
    - Doc gaps found and fixed in `docs/quickstart.md`: Mode 2 writes its index to `output/index/`
      (was undocumented); the documented command without `--cell-calling none` fails on a 2-barcode input
      (emptyDrops error, `CellCallingError`); added the `validate-run` step.
    - NOT executed (user-owned input): Mode 3 (`-acc`) and `build-ref` need an NCBI contact email and
      network downloads (Ensembl host); the Mode 3/`build-ref` commands remain parse-checked only
      (`tests/test_docs_sources.py`). Needs a user run with a real email to close.

### WP8B — Vignettes and claims

- [x] `DOC-06` — eight notebooks exist but contain pre-v3 calls/values; rebuild
  them with negative and ambiguous examples and only v3 APIs/artifacts.
  - Done 2026-10-02 (`v3/vignettes`, c92b768). All 8 notebooks use v3
    APIs and columns:
    - `em-global`/`em-cell` and the `selected-method` primary call;
    - `viral_molecules_total_est`;
    - external cell calling, with knee as sensitivity-only.

    Legacy numbers are removed: 3.64×, 13–18 %/4,414, 28,922/30,849/19,920,
    and the hostresponse 0.87/0.967. Negative examples are in
    `specificity_true_negative`, `cell_calling_denominators` and
    `cell_type_enrichment`. Ambiguous examples are the HHV-6 sibling bleed and
    the rare-host EM case. The notebooks are committed output-stripped.
    - Verification: the 6 CI notebooks pass when their cells are run in a
      single exec namespace. nbclient could not run because no env has
      `ipykernel`, so CI `nbmake` is the kernel-path confirmation.
    - `results/hostresponse_ebv_matched/` is now unreferenced but still
      tracked; removing it is left to the user. DOC-07 is untouched.
- [~] `DOC-07` — execute six lightweight notebooks in CI and the reference/full-
  - 2026-10-04: the "3 of 6 fail on column drift" note is stale. DOC-06 moved the
    notebooks to v3 columns, and all 6 CI notebooks pass cell-by-cell.
    Not run through a real kernel, because nbconvert is not installed.
    Remaining: the scheduled full-notebook workflow and its logs/hashes.
  - 2026-10-12 (`codex/docs-pass`): the 6 CI notebooks now pass through a REAL kernel with the exact CI
    command (`MPLBACKEND=Agg pytest --nbmake` on the six `ci.yml` paths, `viralscan_test_full` env with
    nbmake 1.5.5 + ipykernel, package on `PYTHONPATH`): 6 passed in 76.71 s, notebooks left unmodified.
    The two `[skip-ci]` notebooks (`quickstart_fastq_to_viral_load`, `building_a_reference`) were not run,
    as the row requires. The 2026-09-27 "3 of 6 fail" note below is superseded. Still open: the
    scheduled workflow that runs the two heavy notebooks and saves logs/hashes (`.github/workflows/` is
    owned by another agent; needs a workflow change plus a real reference/index).
  workflow notebooks in the locked scheduled workflow; save logs and hashes.
  - **2026-09-27:** 3 of the 6 CI notebooks fail today, and they fail
    identically at `af4d5b3`. The breakage is notebook/API drift, not WP1C.
    Found while checking WP1C with nbclient, because `nbmake` is not installed
    locally.
    - `cell_calling_denominators`: `KeyError: "['viral_umi', 'total_umi'] not in index"`.
    - `qc_and_read_evidence`: `KeyError: 'viral_molecules_total_est'`.
    - `specificity_true_negative`: `KeyError: 'total_umi'`.
    - **Passing in both trees:** `multimapping_correction`,
      `cell_type_enrichment` (identical `padj` tables) and
      `host_response_depth_control`.
    - Fix the notebooks' column names against the current output schema before
      closing this row.
- [~] `DOC-08` — add a balanced five-workflow pros/cons table generated from the
  harmonized benchmark rather than rhetorical claims.
  - 2026-10-12 (`codex/docs-pass`): software part done. `scripts/gen_workflow_comparison.py [--check]`
    renders `docs/workflow_comparison.md` (in the Sphinx toctree; `-W` clean) from
    `workflow_matrix.workflows` in `analysis/v3_validation/protocol.yaml`; `tests/test_workflow_comparison.py`
    checks freshness, five families covering every non-legacy workflow once, and strengths plus
    limitations in every row.
    - Reading taken (the plan does not define "five"): the protocol's workflows grouped by strategy into
      ViralScan combined, ViralScan exact-fragment two-step, STARsolo, traditional host subtraction, and
      dedicated tools (Venus/Viral-Track/VIRTUS); the 2.2.0 arms are a predecessor comparison and are
      footnoted.
    - The page prints NO performance number and no ranking: the matrix is unexecuted
      (`workflow_matrix.status: pending`, `analysis/dedicated_tool_comparison/README.md` "not yet run"), and
      the only measured data in the repo, the pre-freeze 2x2 in
      `analysis/reference_strategy_benchmark/outputs/harmonized_2x2_unique.tsv`, predates `CMP-00` and is not
      a `validated_v3` artifact. Still open: add measured per-workflow columns once the `CMP-00..05` results
      exist (G5).
- [~] `DOC-09` — add a claim-registry schema, validator, stale-hash detection, and
  - 2026-10-04: the validator already existed. Added: `--coverage` over every
    public doc, and a report-only `--report-unmarked` heuristic for numeric
    or comparative lines without a claim marker. It finds about 30 lines
    across cli_reference, faq, output_reference and CHANGELOG; some are false
    positives such as parameters, and triaging them is a user decision.
    The manuscript is not yet in scope. Run
    `python scripts/validate_claim_registry.py --coverage --report-unmarked`.
  coverage check for README, docs, manuscript, tables, figures, and captions.
- [~] `DOC-10` — require each quantitative, comparative, validated, performance,
  specificity, and installation claim to resolve to a `validated_v3` artifact.
  - 2026-10-12 (`codex/docs-pass`): registered the two claims that an existing artifact backs, each
    with a marker in `docs/faq.md`, a pinned `software_invariant` row in `claims/registry.json` and
    inventory rows (commits 695f7be, d6d384b; registry row re-pinned in the third commit):
    `lod95-poisson-depth-floor` (ln 20 = 2.996 observed / 5.991 true at capture 0.5;
    `tests/test_detection_sensitivity_contract.py`, `tests/test_sensitivity.py`) and
    `capture-cliff-substitution-model` (0.25 / 0.063 / 0.013 at 10 / 15 / 20 % divergence, 90 bp, k=31;
    `tests/test_sensitivity.py`). `--coverage --report-unmarked` flagged 36 lines before, 32 after.
  - Still UNCOVERED (no existing artifact backs them; a user or later task must supply one, or the wording
    must go): `docs/cli_reference.md:210-223` host-filter `p_unique` table and "3.5x", `:479` "lost 96.5%
    of its reads"; `docs/faq.md:195-196` TTV panel k-mer share (1.36 % / 85.8 %), `:433` "~8 % of the human
    genome" (literature, needs a citation entry), `:149` prevalence wording; `CHANGELOG.md:162,213,336`
    "~3x faster / ~6x lower memory" performance claims (`analysis/multimap_profiling` is untriaged);
    `docs/output_reference.md:239,280,318,353,440,469,701,714` (not touched here, the capture numbers at
    `:353` are covered by the new claim but carry no marker; M03/AUD-13 own that file; `CHANGELOG.md:25,61,63`
    likewise repeat the covered numbers). False positives to allow-list (parameters, badge, not claims):
    `README.md:13`, `cli_reference.md:71,187,597,642,710`, `reference_panel.md:60`, `faq.md:327`.
  - Not done: the `--report-unmarked` gate stays report-only, because turning it into a failing check
    needs the triage above (a user decision); and installation claims have no registry artifact: the
    DOC-05 run log is not committed or inventoried.

`G6a` passes when Sphinx warnings are errors, links and shell examples pass,
all notebooks execute in their declared tier, parser/default parity is exact,
and claim coverage contains no missing or stale artifact.

Documentation gate:

```text
python3 -m sphinx -W -b html docs docs/_build/html
python3 -m sphinx -W -b linkcheck docs docs/_build/linkcheck
PYTHONPATH=src python3 -m pytest tests/test_docs_consistency.py -q
```

## WP9 — Regenerate and review the manuscript

Objective: create a submission package exclusively from the frozen v3 bundle.
Estimated effort: 1-2 writing weeks after `G5`, excluding author review.

- [ ] `MS-01` — replace the historical draft and hard-coded legacy figure script
  with generators for registered values, tables, figures, captions, and
  supplement sourced only from `analysis/v3_results_bundle/`.
- [~] `MS-02` — frame molecule-aware host-virus ambiguity, combined/two-step
  evaluation, read-level specificity/QC, reference provenance, and honest
  evidence tiers; remove every legacy count and unsupported superiority claim.
  - Blocked (2026-10-02): the prose waits on G5/MS-01. This pass's ship-doc removals of the covid counts are logged under `DOC-01`/`DOC-02`, not here.
- [ ] `MS-03` — rerun EBV host response using v3 labels with depth and
  mitochondrial controls; retain only if it replicates. Keep TTV solely as a
  host-homology case study absent orthogonal validation.
- [!] `MS-04` — authors must supply author order, affiliations, ORCIDs, CRediT,
  lead contact, funding, conflicts, acknowledgments, and ethics/consent/data-
  access text. Exclude private COVID libraries unless all requirements are met.
- [ ] `MS-05` — verify every citation, run independent methods/statistics and
  claim-to-artifact reviews, resolve all findings, and generate the final
  submission/availability package.

`G6` passes when every number/table/figure resolves to a frozen artifact digest,
there are no placeholders or ineligible private results, citations support their
sentences, and independent review has no unresolved correctness, denominator,
annotation, uncertainty, ethics, or unsupported-claim finding.

## WP10 — Publish and test `3.0.0rc1`

Objective: expose the exact candidate artifacts to real users without promoting
them as stable. Estimated release work: 1-2 days; testing window: 2-4 weeks.

- [ ] `RC-01` — freeze `3.0.0rc1` version, changelog, migration notes, CFF, locks,
  reference manifest, checksums, SBOM, licence report, and release notes on one
  green protected-main SHA.
- [!] `RC-02` — maintainers configure/approve TestPyPI or PyPI trusted publishing,
  GHCR/GitHub release permissions, and public artifact hosting.
- [ ] `RC-03` — tag `v3.0.0rc1` only after exact-SHA gates; publish prerelease
  wheel/sdist, versioned OCI, SIF, locks, checksums, attestations, reference
  archive, and GitHub prerelease. Do not update `latest`.
- [!] `RC-04` — recruit at least three external laboratories covering supported
  chemistries; give each a tester packet for clean install, tiny workflow,
  `validate-run`, and one supported real workflow.
- [ ] `RC-05` — track every correctness, data-loss, installation, documentation,
  and reproducibility failure to closure; rerun all gates after fixes.

`G7` passes when all three testers submit complete manifests, each supported
path works from public artifacts, and no release-blocking issue remains.

## WP11 — Publish stable `3.0.0` and submit

Objective: publish one exact reviewed build everywhere, archive it, then submit
the methods manuscript. Estimated release work: 1-3 days after approvals.

- [ ] `STB-01` — freeze final truth-panel/comparator/manuscript bundles and prove
  the release candidate's fixes did not change scientific outputs unexpectedly.
- [!] `STB-02` — reserve the Zenodo software DOI and benchmark/archive DOI, keep
  them distinct from the reference-data DOI, and provide them for metadata.
- [ ] `STB-03` — update CFF/metadata, pass exact-SHA gates, tag `v3.0.0`, and
  publish wheel, sdist, versioned plus `latest` OCI, SIF, locks, reference archive,
  release notes, migration guide, checksums, SBOM, licences, and attestations.
- [ ] `STB-04` — after the real PyPI sdist exists, insert its SHA-256 into the
  Bioconda recipe, lint/build/test it, submit the PR, and verify the installed
  package on the tiny workflow.
- [ ] `STB-05` — archive software and benchmark artifacts, verify every public
  DOI/link/download, update `CITATION.cff`, then submit the manuscript.

`G8` passes when the tagged SHA is the reviewed green SHA, public distribution
formats reproduce validated counts, all links and DOIs resolve, the software and
benchmark archives are durable, and the manuscript has been submitted. Journal
acceptance timing is external and is not a software completion condition.

## WP12 — Observe and maintain

- [~] `OPS-01` — publish supported-version, patch-release, deprecation, and
  - 2026-10-04: drafted in `docs/release_policy.md` (0b89ad5), with contacts and
    targets as placeholders pending REL-15. Marked as a draft awaiting
    maintainer approval.
  - 2026-10-12 (`codex/docs-pass`): checked the draft against `SECURITY.md`, `SUPPORT.md`,
    `release.yml`, `environment.yml`, `pyproject.toml` and the reference-manifest schema; added a "What the
    repository pins today" table sourced from those files, and tagged every policy choice
    `DRAFT: maintainer decision` (non-GitHub contact routes, previous-minor support window, version-part
    mapping, patch triggers, deprecation window and the absence of a deprecation mechanism in code,
    triage/fix targets). Stays `[~]`: closing needs maintainer approval and `REL-15` contacts.
  security-response policies with named contact routes.
- [~] `OPS-02` — establish reference-update cadence and manifest compatibility
  - 2026-10-04: drafted in `docs/reference_maintenance.md` (0b89ad5): invariants,
    a placeholder cadence, and manifest schema_version reader rules.
    Needs maintainer sign-off.
  - 2026-10-12 (`codex/docs-pass`): verified the manifest-compatibility section against the code
    (`load_build_manifest` MAJOR check in `virus_identity.py`; reference-manifest schema accepts
    `schema_version` 3.0.0 and 3.0.1, now stated). The invariants stay as drafted. The review interval, the
    update triggers and what a new version ships with are tagged `DRAFT: maintainer decision`; nothing in the
    code enforces a cadence. Stays `[~]` pending maintainer sign-off.
  rules; never silently change a frozen reference under an existing version.
- [ ] `OPS-03` — triage false positives, false negatives, data loss, and reference
  drift as scientific incidents with reproducible packets.
- [ ] `OPS-04` — monitor install/usage failures during the first 90 days and ship
  patch releases from the same gate process.
- [ ] `OPS-05` — schedule independent truth-panel refreshes without changing v3
  thresholds retrospectively.

## Gate dashboard

| Gate | State | Required proof |
|---|:---:|---|
| `G0` governance | `[~]` | prescribed PyPA-frontend archive build and member check; governance validators otherwise pass |
| `G1` software | `[x]` | full unit/property/safety/tiny-workflow suite. 2026-10-05: green in the full env `viralscan_test_full` (integration 73/73, no skips; unit 1,715; ruff clean) |
| `G2` distribution | `[~]` | clean installs, locks, Docker/Apptainer parity, supply-chain reports |
| `G3` preregistration | `[ ]` | reviewed, schema-valid, hashed protocol frozen with no recorded pre-freeze outcome on scored or protected pending-row datasets. 2026-10-13: Q83 choices accepted and pasted (`DEV-028`); next is SCI-05 round 13, then the user's signature and the freeze commit |
| `G4` references | `[~]` | byte-rebuild, GRCh38 D-list, calibrated holdout safeguard |
| `G5` science | `[ ]` | truth holdout, public positives, all comparators, uncertainty, frozen bundle |
| `G6a` docs | `[~]` | parser parity, executable docs, complete validated claim graph |
| `G6` manuscript | `[ ]` | artifact-generated manuscript plus independent audits |
| `G7` release candidate | `[ ]` | public RC artifacts and three-laboratory closure |
| `G8` stable/publication | `[ ]` | stable artifacts, archives/DOIs, Bioconda, manuscript submission |

## External inputs and authority

These items cannot be invented or completed by an implementation agent:

| Input | Needed by | Owner action |
|---|---|---|
| Author/ethics/funding/conflict metadata | `REL-15`, `MS-04` | authors approve final factual text |
| Publishing credentials and protected environments | `RC-02` | repository owner configures services |
| Three external laboratories | `RC-04` | maintainers recruit and coordinate testers |
| Zenodo DOI reservations | `STB-02` | archive owner reserves distinct DOIs |
| Orthogonal anellovirus-positive sample | `REF-10` | collaborator supplies lawful validated data, or claim stays screening-only |
| Published viral annotation panel archive and registered DOI | `REF-11` | archive owner publishes the panel and registers `10.5281/zenodo.20112332`, or the pinned record identifier is corrected |

## Stop rules

- Never tag or publish from a dirty, unreviewed, or non-green SHA.
- Never migrate or reinterpret a pre-v3 H5AD value as a v3 molecule count.
- Never select shared barcodes using ViralScan-positive outcomes.
- Never turn ambiguity/QC flags into biological conclusions automatically.
- Never claim real-anellovirus sensitivity, general superiority, or formal
  infection from a nonzero molecule without the required validation evidence.

## Evidence log

Append one line after each completed item:

```text
YYYY-MM-DD ITEM — command/result; artifact path(s); Git SHA; reviewer if required
```

- 2026-07-22 `BASE-03` — EBV v3 baseline conserved exact molecule mass;
  `analysis/v3_ebv_baseline/host_conservative.json`; SLURM 25316336.
- 2026-07-22 `BASE-04` — real-tool evidence integration passed;
  `tests/integration/test_exact_lineage.py` and
  `tests/integration/test_evidence_chain.py`; Git SHA `a146050`.
- 2026-07-22 `BASE-05` — full non-network suite: 720 passed, 21 deselected;
  Ruff check/format, data-governance check, and draft protocol validation passed.
- 2026-07-22 `SW-01` — wheel and sdist contain all six byte-matched v3 schemas;
  fresh wheel install loaded 6/6, missing-schema paths fail closed; 8 focused
  tests passed; Git SHA `6e1fe66`.
- 2026-07-22 `REL-04` partial — isolated wheel/sdist build and fresh-wheel schema
  smoke test passed; `twine check`, sdist clean-install test, and `SHA256SUMS`
  remain.
- 2026-07-22 `DOC-01`–`DOC-04` partial — public docs and tracked manuscript
  placeholder reject legacy quantitative/private claims and use v3 molecule and
  candidate-evidence terminology; 94 docs/CLI tests passed; Git SHA `3bb7d1b`.
- 2026-07-22 `SCI-01` — schema-valid, non-executable protocol draft with separate
  training/holdout gates; 16 focused validation tests and 4 docs-consistency
  tests passed; data-governance check passed; independent implementation review
  resolved four blockers and passed; `analysis/v3_validation/protocol.yaml`,
  `schemas/v3/validation_protocol.schema.json`, and
  `scripts/validate_v3_protocol.py`; Git SHA `d941261`.
- 2026-07-22 `SCI-02` — outcome-independent harmonization frozen with 24
  endpoint-mapped denominator contracts, canonical and dependent-field digests,
  49 focused adversarial tests, schema/CLI fail-closed checks, data-governance
  pass, and three independent reviewer passes; protocol remains non-executable;
  Git SHA `d941261`.
- 2026-07-25 `GOV-06` partial — all 44 retained BUS rows validated and
  aggregated; six downloaded EBV mates passed authoritative size/MD5 and local
  SHA-256 gates; pair/chemistry audits `25330878` and `25330881` completed all
  five controls; `analysis/legacy_v2_v3/control_inputs.tsv` is the sanitized
  shared manifest. Fresh attempt 1 (`25331024`, `25331026`) failed before
  execution because Slurm-spooled scripts could not resolve the packet root;
  attempt 2 uses a tested explicit packet root in arrays `25331035` and
  `25331037`. Full suite: 758 passed, 48 deselected; base Git SHA `a8f3287`.
- 2026-07-26 `GOV-06` partial — attempt-2 arrays `25331035` and `25331037` are
  terminal with all ten rows failed and every outcome retained. Three
  independent causes, each confirmed from primary evidence
  (`fresh_control_packet_attempt2/status/*.json`, `sacct`, on-disk output
  trees): (1) `run_fresh_control.py` `_v2_artifact_errors` checked a flat output
  layout while the legacy 2.2.0 CLI nests one level under the sample
  identifier, so `v2__SRR12682296`, `v2__SRR12682297`, `v2__SRR12682298`, and
  `v2__SRR6825025` recorded exit 65 with six phantom missing artifacts despite
  `workflow_exit_code` 0 and complete non-empty output trees; (2) all five v3
  rows raised `ViralScanDataError` during config creation because the Zenodo
  viral-annotation cache was never populated in the frozen packet environment;
  (3) `v2__SRR6825024` was genuinely out-of-memory killed in the legacy
  `multimap.py` at roughly 121.4 GiB peak resident set against the frozen
  128 GiB tier ceiling after 5 h 52 m. A fourth defect was found while
  verifying: the status payload satisfies none of `stage`, `attempt_id`, or
  `scientific_parameter_hash` from `protocol.yaml` `required_failure_fields`,
  so all ten records are non-compliant as failure records. Attempt 3 was
  prepared as a v3-only five-row packet and deliberately not submitted; the
  diagnostic is outcome-ineligible and yields priority to `SCI-03`.
  **Superseded 2026-07-28**: this entry originally read "Attempt 3 is frozen".
  It is not, and no attempt-3 packet exists on disk — only
  `benchmark_runs/legacy_v2_v3/fresh_control_packet_attempt2/`. The `GOV-06`
  work-package row and `analysis/legacy_v2_v3/TRACKER.md` both correctly record
  it as blocked on `REF-11`. Same drift class as the `SCI-03`/`SCI-04` entries
  corrected on 2026-07-27; found by the 2026-07-28 multi-agent review.
- 2026-07-27 `SCI-03` partial — `partitions` and `calibration` are frozen,
  schema-validated sections of `analysis/v3_validation/protocol.yaml`. Partitions
  allocate whole biological samples by deterministic stratified assignment over
  the declared factors at a 0.3 holdout fraction, with a single permitted holdout
  evaluation and explicit template/locus/molecule/cell-barcode leakage
  prohibitions. Calibration fixes a training-only threshold grid with a
  zero-false-positive constraint on host-only and planted-homology negatives, a
  conservative tie-breaker, nine endpoint-linked metrics, a probit LOD95 with
  extrapolation prohibited, and a 2000-replicate biological-sample bootstrap that
  forbids treating cells, molecules, reads, or technical repeats as independent
  replicates. Canonical digests: partitions
  `44173743fc7b0c3e9a48b165348faca196c535bd82eee64f4ea771724a3beb04`, calibration
  `252fa8ade5381266887527cde3ddc139a289e81d18fec9be1ce99da4fe4d6d7c`. Seeds
  `root`, `split`, `cell_calling`, `evidence_sampling`, and `bootstrap` are
  frozen; `generation` remains `VAL-01`. The `partitions_metrics` training
  blocker is closed; `data_hashes` stays open. Validator gains `_validate_sci03`
  with digest, seed, endpoint, and factor cross-checks. Draft gate valid;
  training gate still correctly blocked on `SCI-04`, `SCI-05`, `VAL-01`, and
  `REF-09`. Full suite 790 passed, 48 deselected.
- 2026-07-27 `SCI-04` partial — `workflow_matrix` and
  `failure_and_deviation_reporting` are frozen, schema-validated sections. The
  matrix enumerates eight workflows over 52 rows: ViralScan combined across all
  twelve datasets, ViralScan exact-fragment STAR two-step, STARsolo combined and
  host-filter two-step, traditional host subtraction with and without CB/UMI
  retention, and Venus and Viral-Track at published defaults in isolated
  digest-pinned environments. Both dedicated comparators cover EBV, HHV-6B, and
  HSV-1 plus the host-only and host-homology negatives. The kallisto two-step
  path is excluded with an explicit revisit condition rather than omitted.
  Seventeen per-row record fields are required, and five primary comparison rules
  fix identical FASTQ hashes, host release, the outcome-independent anchor, the
  audited feature intersection, and unique-only parity. Canonical digests:
  workflow matrix
  `fb4c6768e87a87f23fa554c1376cd6b3ce422d972d3f2e409d83985c154bde4c`, failure
  reporting `ce913b3ea10f23062ab0f86bbf81ce341ea9b02edad594d57cbab363c5e42975`.
  `_validate_sci04` enforces digest match, planned-section pairing, dataset and
  reference resolution, and row-count agreement, and at the training gate refuses
  any workflow lacking a pinned tool version and container digest.
  `environment_pinning` stays `pending` on `REL-03`, so the `workflow_rows`
  blocker is closed but `tool_environments` remains. Full suite 800 passed, 48
  deselected.
- 2026-07-26 `REF-11` opened — the Zenodo record pinned by
  `src/viralscan/data_fetch.py` is unregistered. `zenodo.org/api/records/20112332`
  returns `{"status": 404, "message": "The persistent identifier is not
  registered."}` and `doi.org/10.5281/zenodo.20112332` returns 404, verified from
  a network-capable host on which an unrelated third-party Zenodo DOI resolved
  200. The installed package under
  `benchmark_runs/legacy_v2_v3/env_full/.../viralscan/data/` contains only
  `anellovirus_accessions.tsv` and `__init__.py`, while 195 GTFs remain in
  `src/viralscan/data/` in the source tree. `viralscan data fetch` therefore
  cannot succeed for any user, which blocks the `GOV-06` attempt-3 cache pin and
  the clean-install paths behind `REL-05` and `DOC-05`.
- 2026-07-27 evidence-log correction — the `SCI-03` and `SCI-04` entries above,
  both dated 2026-07-27, state that `partitions`, `calibration`,
  `workflow_matrix`, and `failure_and_deviation_reporting` "are frozen". That was
  true when written and is no longer: `SCI-05` round 1 restored blockers and all
  four sections returned to `status: pending`, where they remain. The work-package
  rows carry the authoritative status. The entries are left in place because this
  log is append-only history, not current state.
- 2026-07-27 six-agent code review of `codex/viralscan-v3` vs `main` —
  `.living/outputs/reviews/2026-07-27-branch-codex-viralscan-v3.md`, with a
  behavioural tripwire audit alongside it. Two Majors, both fixed here. (1) The
  protocol froze `seeds.cell_calling: 20260727002` and named it the seed source
  for the shared cell anchor, but `call_cells` never passed `seed` to
  `emptydrops_cells`, so emptyDrops always ran at the signature default `100`;
  `emptydrops_seed` and `emptydrops_niters` are now declared config fields with
  CLI flags, and `emptydrops_cells` takes keyword-only required parameters so no
  future caller can omit one silently. (2) `check_git_sha_fields` accepted a
  record's before-digest unverified when its scoped section was absent at the base
  commit or when `git show` failed — the open `R11-F2`; both branches now error.
  Minors fixed: frozen FASTQ identity is re-derived from the bytes in
  `verify_frozen_fastq` rather than trusted from a weeks-old audit (kept out of
  `prepare_fresh_controls`, where it would re-read ~250 GB to close no additional
  window); `Rscript` is preflighted when the resolved cell caller is `emptydrops`,
  since `SW-11` made that path fail closed after `kb_count`, `analysis`, and
  `multimap`; `detect_cells` documents `viral_count_matrix`. One review finding
  was reversed on verification: `SIBLING_VIRUS_PAIRS` omits the protocol's
  EBV/KSHV pair *correctly* — the constant is a runtime EM-bleed heuristic
  requiring near-identity, the protocol list is an evaluation population spanning
  a relatedness gradient, and adding EBV/KSHV to the constant would annotate
  genuine co-infection as artifact. Both sides now say so.
- 2026-07-27 `DEV-019` — the protocol edit above changed
  `partitions.contract_sha256` to
  `5fd9b366001c0101dd174c79d2d446cbc957525b0a375e659937d4aaa4eb1a05`. An earlier
  draft of this entry claimed no ledger record was required because `partitions`
  is `pending`. That was wrong and unverified: `partitions` carries an
  eight-link chain that was reconciled at `0dae092a…` before the edit, so the
  edit orphaned it. `DEV-019` records the change, and declaring it re-digests the
  frozen `frozen_inputs` section, carried in the same record's
  `additional_digest_changes`. All five sections with chains now reconcile against
  the live protocol. The gap that let this happen is closed by
  `pending_section_ledger_drift`: `validate_amendment_ledger` binds only *frozen*
  sections, so a pending section that already has records could drift unnoticed
  until the moment it was frozen. That check now blocks the training phase, which
  must pass before any freeze.
- 2026-08-08 `GOV-03` — `analysis/v3_artifact_inventory.tsv` now has 30
  sanitized, stable-identity rows spanning claim-bearing inputs, references,
  intermediates, results, commands, environment, scheduler evidence, failure
  evidence, schemas, and documentation. The fail-closed schema/hash/cross-link
  validator passed; the provenance-incomplete retained EBV baseline remains
  private and ineligible for a public validated-v3 claim. Base Git SHA
  `26260e7cbc18cc0e7777379e1db12cc594778dd4`.
- 2026-08-08 `GOV-04` — `config/public_ship_scope.json` is the single positive
  wheel, sdist, Docker-context, public-documentation, claim-bearing, and
  governance-text allowlist. Packaging no longer recursively includes `docs/`,
  Docker no longer uses `COPY .`, the allowlisted-text institutional-path scan
  passed, and wheel/sdist/context members matched exactly when the distributions
  were built through the installed setuptools backend. Base Git SHA
  `26260e7cbc18cc0e7777379e1db12cc594778dd4`.
- 2026-08-08 `GOV-05` — the eight-record claim graph is schema-valid and covers
  every allowlisted public marker with checked artifact identities and hashes.
  Legacy counts are explicitly ineligible, reference homology is implemented but
  not calibrated, and the retained EBV baseline is provenance-incomplete rather
  than promoted. Focused governance tests: 18 passed. All three governance
  validators and the draft protocol validator passed. Base Git SHA
  `26260e7cbc18cc0e7777379e1db12cc594778dd4`.
- 2026-08-08 `G0` partial reconciliation — the default non-network suite passed
  958 tests with 73 deselected; Ruff check and changed-file format checks passed.
  The prescribed `python -m build --no-isolation` command is unavailable because
  the locked runtime has no PyPA `build` frontend; no network install was
  attempted. Repo-wide Ruff format remains red on ten protected/out-of-scope
  pre-existing files. A separate corrected-PATH integration run produced 19
  passed, 2 failed, and 17 errors: the shared real-workflow fixture requires a
  missing `conda` executable, and this Snakemake emits no DAG rule listing under
  the tests' `--quiet` invocation. `G0` remains `[~]` and no later gate is
  promoted.
- 2026-09-27 housekeeping — the WP4B2–WP4F work was squashed into `3379b7c`
  so every commit passes its own suite (the first half had registered
  `hpv_genes.*` and `anellovirus_genes.tsv` in the ship-scope allowlist before
  those files existed); the result tree is byte-identical to the pre-squash
  history, kept on `backup/wp4-pre-squash`. The schema 1.1.0 governance
  migration was finished in `f447ff3` (config at 1.1.0, redundant `text_files`
  dropped — all 86 of its entries are already covered by
  `ship_scope_source_paths()`) and its three self-referential inventory rows
  pinned in `7b879f5`. `471f889` makes the two bundled-GTF tests skip when the
  gitignored panel is absent. A clean-checkout worktree at `471f889` passed
  **1,260, skipped 6, failed 0**; `check_data_governance.py` and
  `validate_claim_registry.py --coverage` pass. Nothing pushed.
- 2026-09-27 `PROG-07` re-measurement — EBV LCL `SRR12682296` gene programmes
  after removing `BARF1.2`/`BaRF1.1` from the latent set: unique layer 895
  latent / 236 productive / 311 mixed / 3,094 indeterminate (was 2,240 / 102 /
  445 / 2,968); allocated layer 856 latent (was 1,277). Summary sha256
  `48d0a9d12fbc6ca0…` (pre-change `72218893ff71baa3…`, kept as
  `gene_program_summary.pre-PROG-11.tsv`). Git SHA `50253a6`.

- 2026-09-27 `TONSIL-01` — SFL tonsil x223 (5′ v3 GEX, 24 donors), host-subtracted
  screen from the cellranger BAM. Stages: 118,856,604 unmapped GEX reads →
  24,969,278 pass the prefilter → 34,208 viral hits → 1,511 pass the viral
  filters → 23 survive the host re-check. Anelloviridae 0 and HPV 0 in all 24
  donors. Plant recovery: HPV16 99.2 %, TTV 99.4 %. Scripts and outputs are in
  `benchmark_runs/sfl_tonsil_screen_2026-09-26/tonsil01/` (gitignored): SLURM
  25652114 (extract), 25652133 (screen), 25652132 (plant). `calls.tsv` sha256
  `6f853b9f36dbe1e3…`. All 4 FASTQ md5 checks pass (25652115). Opened at Git SHA
  `a86aa6e`.
- 2026-09-27 `TONSIL-01` sensitivity to held-out strains (SLURM 25652148):
  - Held out: MN770908.1 (Beta), MW679005.1 (Alpha), MW455373.1 (Gamma).
  - Recovery, 5′ window / uniform: Beta 87.7 % / 72.9 % (nearest 94.3 %),
    Alpha 15.5 % / 34.7 % (nearest 84.0 %), Gamma 39.5 % / 8.0 % (no asm20 hit).
  - `plant_ho/out/results/reads_final.tsv` sha256 `6dd488146ceceb78…`.
  - The scripts are gitignored, so this commit's SHA does not pin them. Their
    sha256 prefixes: `build_db.py` df022722c37c, `prefilter.py` 421964e9f693,
    `parse_hits.py` fca465b4fe5c, `plant.py` c1d2bd0b07cd, `heldout.py`
    cd98664c5004, `screen.sbatch` ebb86194b0fb, `extract.sbatch` 54a1b2c3d9a3,
    `plant.sbatch` d07115d49705, `heldout.sbatch` 5bac5aea2a8b.
  - `screen.sbatch` gained an optional reference-dir argument after the main
    run. Its default is unchanged.
- 2026-09-27 `SIMP-01` — dead code deleted. Unit suite 1,275 passed, 0 failed.
  A diff of collected test IDs against `af4d5b3` shows exactly one test
  removed, `test_config_value_serializes_none_as_empty_string` (its function
  was deleted). `tests/test_multimap.py` is inventoried, so its docstrings
  that name the removed `normalize_barcodes` were left unchanged; re-pinning
  it for wording alone is not worth it.
- 2026-09-27 `SIMP-02` — `src/` helpers collapsed, 186 lines removed and 84
  added. Unit suite 1,273 passed. The 2 failures are the expected stale
  `build_reference.py` pins, fixed by the follow-up pin commit. Collected test
  IDs are unchanged.

  Integration outcomes match `af4d5b3` exactly: 18 passed, 19 skipped, and 1
  failure that predates this work. The failure,
  `test_build_anellovirus_reference_produces_labelable_gtf`, still expects the
  pre-WP4F `_gene1` placeholder IDs.

  A before/after probe against an `af4d5b3` worktree gave identical results
  for:
  - the packaged-table paths and row counts;
  - sha256 of a `Path`, a `str`, and `hashlib`;
  - the NCBI cache validity check;
  - the gzip opener;
  - the matrix axis sums.

  All 8 `--help` screens exit 0. The `--verbose`/`--quiet` help text is now the
  same on every subcommand.
- 2026-09-27 `SIMP-03` — dependencies: `pyfiglet` and `seaborn` removed, and
  `scipy` pinned `>=1.11`.

  **BH swap (`enrichment._bh_adjust` → `scipy.stats.false_discovery_control`):**
  - Across 7 input sets, including 300 real `fisher_exact` p-values, old and new
    agree to 1.1e-16, which is one floating-point ulp.
  - NaN and p > 1 cannot reach `_bh_adjust`. `fisher_exact` clips to [0, 1],
    and a brute-force run over 1,296 tables with counts 0–5 found none outside
    that range.
  - The new `test_bh_adjust_matches_reference_step_up` pins the equivalence.

  **Plots moved to matplotlib:**
  - Bar heights, tick order and histogram bin counts are identical to
    seaborn's. Seaborn's error bars had zero length.
  - The histogram keeps seaborn's axis labels.

  **Banner:** the welcome banner is plain text unless `pyfiglet` happens to be
  installed. `menu.py` already fell back, which makes the `tests/conftest.py`
  stub redundant, so it was deleted.

  **Gates:**
  - Unit suite: 1,275 passed. The 1 failure is the expected `environment.yml`
    pin, fixed in the pin commit.
  - Collected test IDs: +1, the new BH test.
  - Integration outcomes: identical to `af4d5b3`.

  **Changed files:** `pyproject.toml`, `environment.yml`,
  `conda-recipe/meta.yaml`, `.github/workflows/ci.yml` and `release.yml`.
- 2026-09-27 WP1C regression follow-up. CI `mypy src/viralscan` found 3 new
  errors from SIMP-02/03, all missing or too-narrow annotations with no runtime
  change:
  - the `_bh_adjust` return type;
  - `_open_maybe_gzip` now accepts `str | Path`;
  - `_write_tsv` parameter types.

  After the fix, mypy matches `af4d5b3` exactly: 57 errors, 0 new. The
  untracked scratch file `kmer3.py` was excluded.

  Other checks:
  - **Names:** an AST check of all 606 `viralscan` names imported across 174
    tracked `.py` files found none missing.
  - **Smoke test:** the CLAUDE.md smoke test passes, including with `pyfiglet`
    blocked, where the fallback banner is used.
  - **Containers:** `Dockerfile` and `Singularity.def` build from
    `environment.yml`, so there is no parity drift.
  - **Unit suite at `e0a7b1f`:** 1,276 passed.

  Added `ANELLO-14` for the integration failure that predates this work.
- 2026-09-27 WP1C final gate at `ff86c3d`:
  - Unit suite: 1,276 passed.
  - Integration: identical to `af4d5b3` (18 passed; `ANELLO-14` fails in both).
  - mypy: 57 errors, the same as `af4d5b3`, 0 new.
  - CI vignettes run with nbclient: 3 pass in both trees; 3 fail in both trees
    with identical errors (logged under `DOC-07`).
  - `cell_type_enrichment` `padj` output is identical before and after the BH
    swap.
- 2026-09-27 `CAT-11` — `scripts/measure_kmer_capture.py` reproduces the
  `REF-01` evidence: EBV `NC_007605.1` 144,283 31-mers; bundled 20-genome panel
  median coverage 0.00 %, zero-coverage 85.80 %, Betatorquevirus 98.4 %
  (n=1,542), Gamma/Samek/Het/Gyro/Mem 100 %; expanded 2,042-genome panel
  leave-one-out 20.55 %, zero-coverage 0.15 %, panel k-mers 4,888,291.
  **Correction:** the expanded panel's fragment capture is a measured **0.5097**,
  not the claimed analytic 1.0000 (F-011). Population
  `viralscan_showcase/.../anellovirus.fa` (2,042 genomes); outputs
  `bundled20.{tsv,json}` and `expanded2042.{tsv,json}` under the session
  scratchpad; 12 tests in `tests/test_measure_kmer_capture.py`. Git SHA `db6c699`.
- 2026-09-27 `CAT-02`/`CAT-03` — seed catalogue and catalogue-driven naming.
  `extras/build_virus_catalog.py --from-fasta references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.fa`
  → `src/viralscan/data/virus_catalog.tsv`: 2,215 accessions, 204 species, 30
  families, 0 failures (170 flatfiles fetched, the rest cached). Naming measured
  on the same reference: unnamed gene IDs 173 → 0, distinct groups 182 → 107,
  and influenza A's 8 segments group into 1 virus instead of 8. 1,307 tests
  pass. Git SHA `f39e18d`.
- 2026-09-27 `CAT-09` first batch — 16 previously-absent viruses added to the
  catalogue, every accession verified by live NCBI lookup: SARS-CoV-2, HIV-1/2,
  HTLV-1/2, HCoV-OC43/NL63/HKU1, hMPV, bocavirus, TSPyV, HPyV6/7, simian foamy,
  influenza D (7 segments) and the 16 HPV genotypes. Catalogue 2,249
  accessions / 232 species / 31 families (was 2,215 / 204 / 30); 714 KB.
  1,307 tests pass. Git SHA `af97cf2`.
- 2026-09-27 `CAT-11` correction — **the capture metric measured the forward
  strand only.** kallisto indexes canonical 31-mers; `kmers()` did not fold
  reverse complements, so a panel genome deposited antisense to a target
  counted as a miss. `--strand canonical` is now the default;
  `--strand forward` regenerates the superseded numbers. Recomputed on the
  2,042-genome anellovirus panel:

  | | forward | canonical |
  |---|---|---|
  | median leave-one-out | 0.2055 | **0.2162** |
  | median p_fragment | 0.5097 | **0.5349** |
  | panel distinct 31-mers | 4,888,291 | 4,843,779 |

  Per genus (median p_fragment, forward → canonical): Alpha 0.788 → 0.808,
  Het 0.837 → 0.837, Beta 0.498 → 0.515, Gamma 0.332 → 0.336, Samek 0.255 →
  0.256, Mem 0.247 → 0.247, Gyro 0.063 → 0.063, unclassified 0.334 → 0.397.
  350 of 2,042 genomes improve and **none degrade**; 84 by > 0.10 absolute and
  55 by > 0.25. Worst case `MH649023.1` 0.0329 → 0.9514. Affected accessions
  cluster in the `MH648xxx`/`MH649xxx` block. 17 tests in
  `tests/test_measure_kmer_capture.py` (was 12); `--self-check` now asserts a
  reverse-complement-only panel still captures its target.
- 2026-09-27 `CAT-12` provenance — **the anellovirus panel's sensitivity
  ceiling is set by upstream dereplication, not by how much diversity exists.**
  The 2,022 non-RefSeq panel genomes come from
  `github.com/clareaulab/human_anellovirus_pangenome`, which took 3,545
  complete human-host anellovirus genomes (NCBI Virus taxon 687329, June 2025)
  and **CD-HIT clustered them at 95 % ANI / 85 % coverage down to 2,023
  representatives**, then resolved genus for ~584 NCBI-unclassified genomes by
  ORF1 protein phylogeny (MAFFT → trimAl → IQ-TREE, 40 clusters).

  Three consequences, each measured:

  1. A 31-mer survives only with zero mismatches, so at the 95 % ANI threshold
     the expected shared fraction is `0.95^31 = 0.2039`. Measured median
     canonical leave-one-out is **0.2162**. The agreement to ~0.01 says
     CD-HIT's threshold, not NCBI's holdings, set our capture. The 1,522
     genomes CD-HIT discarded sit in exactly the divergence band where 31-mer
     pseudoalignment fails, so **restoring them is the one lever that raises
     capture** (2,023 → 3,545, +75 % genomes, ~10.5 Mb). Caveat to test before
     acting: the binomial assumes substitutions are independent and uniform,
     while F-011 showed clustering matters; the two must be reconciled.
  2. **The panel's genus labels are phylogeny-derived and are better than
     GenBank's.** The 582/2,042 panel genomes I found with no GenBank genus are
     the ~584 the upstream ORF1 tree resolved. Upstream reports Betatorquevirus
     1,558 by phylogeny vs 1,360 by NCBI label; measured here independently,
     1,542 vs 1,368 — the same correction. Do **not** "fix" our labels against
     GenBank lineages; that would discard the better assignment.
  3. **The absent Anelloviridae genera are correctly absent.** Every genus
     missing from the panel has **zero** complete genomes with
     `"Homo sapiens"[Host]`, verified by live NCBI query 2026-09-27:
     Lambda 408 complete / 0 human, Eta 256/0, Iota 139/0, Kappa 24/0,
     Rho 19/0, Epsilon 18/0, Pi 14/0, Theta 10/0, Sigma 9/0, Upsilon 8/0,
     Mu 5/0, Xi 4/0, Zeta 2/0; Delta, Nu and Tau have none at all. They are
     swine, feline, canine, tupaia and pinniped viruses. Adding them to a human
     panel would add false-positive surface and no sensitivity. This closes the
     "add more genera" line of `CAT-12` and applies to Gyrovirus too.

  Masking checked and excluded as a confound: the shipped panel FASTA carries
  728 N bases in 5,998,625 (0.01 %), 16 genomes affected, worst 11 %.

- [x] `CLI-17` — user-confirmed 2026-10-13 (Q71-Q76) and implemented on `codex/multimap-3-options`.
  - Q72: `--multimap-method` (and `rerun-multimap`) shows 3 options, `conservative` (default), `em`, `equal`; `unique-weighted`, `sibling-weighted`, `em-global` stay accepted, hidden (no warning, protocol grid untouched).
  - Q73: `em` = `em-cell` because the `EMC-01` re-measure passed (17.6 GB, 1 h 23 min on 942,700 barcodes); `em-global` otherwise.
  - Q74: short names are CLI aliases only (`defaults.MULTIMAP_METHOD_ALIASES`, resolved by the parser `type=`); run manifests, configs, h5ad and the protocol record the canonical engine name (`host-conservative`, `em-cell`/`em-global`, `equal`). Tests: `tests/test_multimap_aliases.py`.
  - Q75: after G5 the default follows the engine the frozen protocol rule selects. Documentation only; no code change until then.
  - Q76: no VBEM.
  - Governance: `public-readme` inventory row repinned to the pushed `eeba10f` (README edit); full suite 3798 passed, 1 failed (that stale pin, now fixed), 13 skipped.
  - Q71 (G4): the WP5 build is the G4 candidate; a separate agent rebuilds once to check byte-identity and does `REF-08`/`REF-09`; `REF-11` waived for the study (`DEC-02` hash pin). Tracked there, not implemented here.

## VR pilot — read cleanup and transcript reconstruction (2026-10-14)

**Takeover 2026-10-14:** Claude took over the uncommitted Codex work with user approval
("Take over, commit to a branch") and ported the code, tests and docs onto
`codex/translated-search-vr-pilot` (off `origin/main`). The Codex "no commit/push" line (AA-02) is
superseded for the code only; real-sample runs (AA-02, VR-03/04) are run separately by other
agents and stay `[~]`. Lint cleanup for VR-02 done in the port (E402/SIM115 `noqa` in
`analysis/viral_read_cleanup/pilot.py`; `analysis/` is outside the ruff gate).

- [~] `VR-01..04` — user authorized the reviewed plan with "run it".
  Separate exploratory pilot; do not change frozen DSR rounds or package defaults.
  Tracker: `docs/plans/2026-10-14-viral-read-cleanup-and-transcript-reconstruction.md`.
  Owner: Codex coordinator until 2026-10-14, then Claude for the code port (see takeover note); limited to the new analysis prototype, its tests,
  this tracker row, and the reviewed plan. Original reads and references remain read-only.
  Pinned optional Cutadapt 5.2 and StringTie 2.2.3 are isolated under `/tmp`.
  Real inputs: separate first-1M-pair subsets of EBV SRR12682296 and HPV16 positive /
  parental SRR19537341 / SRR19537339; exact synthetic tail/orientation controls supplement them.
  Nine focused prototype tests pass. Native STAR/StringTie fixture recovered two transcript
  models (three exons, verified junction support). The first real attempt stopped at a repeat
  coordinate audit; corrected v3 is frozen and dry-run only. Real outcomes,
  source/chemistry qualification and the integration verdict remain pending.
  Attempts/receipts: `.tmp/viral-read-cleanup-20261014*`; no scored G5b/VAL-01 or real-sensitivity claim.
- [ ] `VR-05` — candidate host-virus junction design remains deferred until pilot results.
- [x] `AUD3-07` — audit findings P9-01, P3-05, P1-01 (counting): documentation and pinning tests only; no numerical change.
  `docs/v3_counting_contract.md` now states that `em-global` is an expert value whose stop rule (tol 1e-6) can leave a
  residue on a zero-unique-support viral gene at extreme ratios (2.56 molecules at 1e6 unique / 4e6 ambiguous; tol 1e-9
  gives 0.0077; per-cell 5.9e-4, so `detected` stays false) and must be re-evaluated before a G5 pick could make it the
  default; it also states the gene-less-EC skip convention. `docs/output_reference.md` states the 1e-6 slack on tier
  cut-offs; `docs/cli_reference.md` and the `defaults.py` comment point to it.
  Tests: `tests/test_em_global_residue_pinning.py` (6, under 1 s of test time; the 807 s repro is replaced by
  `em_gene_abundances` with float EC counts). No fail-before: they pin current behaviour. EM rule, tolerance,
  pseudocount, tier slack and EC handling are unchanged by decision.
  Pinned-inventory rows to re-pin: `v3-counting-contract-doc` (`docs/v3_counting_contract.md`) and
  `output-reference-doc` (`docs/output_reference.md`).


- [x] `PERF-B16` (2026-10-19) — sequential anonymous SAM handoff accepted against
  d056986. Native minimap2 finishes before unchanged samtools sort/index argv.
  180 focused tests pass, including real child failures/full stderr/partial files;
  actual minimap2 ENOSPC stops before sort and preserves old output. All retained
  native BAM/BAI artifacts compare raw equal in A/A and A/B. Two independent
  seven-pair sessions: median wall1.00805/1.00163x, sampled process+child RSS
  0.51064/0.49989x; traced peak0.00242/0.00197x. Memory benefit only, no runtime
  improvement claim. Temporary SAM disk28,626,864 bytes for100,000 fixture reads,
  scaling with input. Instrumented stdout-hash timing attempts retained but excluded
  from benefit claims. Fixed workers/search/scientific settings. C00 pending.

## Reviewed roadmap implementation — 2026-10-20

Scope: revised roadmap and implementation handoff, independently checked against
`origin/main` `9151da904cb818ad140c3311f51be7dd5803f738`. The supplied two repository
paths resolve to one dirty, stale checkout; implementation uses an isolated clone.
The baseline was confirmed against GitHub main before edits. No unrelated work
was incorporated. Source implementation is frozen at `ecee46c`.

### Acceptance ledger

All PASS entries below are local software/fixture results, not CI or scientific
validation. Receipt filenames refer to the separate local delivery packet.

| ID | Evidence | Local status |
|---|---|---|
| A1 | `test_complete_universe_legacy_and_fractional`; summary before/after TSVs: 4.0% → 0.2% | PASS |
| A2 | `test_auto_resolution_missing_empty_and_suffix`: full empty header | PASS |
| A3 | Same test: missing authoritative cells is a hard error | PASS |
| A4 | Legacy fixture: warning, 4.0%, explicit denominator kind | PASS |
| A5 | Fractional fixture: 9.6 excluded at 10; 10.0 included | PASS |
| A6 | `test_host_only_null_and_called_universe`: null expectations 2.5 and 7.5 | PASS |
| A7 | `test_missing_host_depth`: missing mask and zero host depth; no estimated null | PASS |
| A8 | Non-called barcode with 100,000 viral molecules excluded in plot fixture | PASS |
| A9 | Exact zero/fractional sidecar values; zero count metadata | PASS |
| A10 | `test_figure_only_rank_cap`: 60,000 TSV rows, 50,000 displayed, annotation captured | PASS |
| A11 | Independent baseline/candidate `detection.run` receipts; exact table/H5AD bytes; mutation test | PASS |
| B1 | Zero-gap plot-array fixture; native samtools 1M coverage-map fixture and paired rendering benchmark | PASS |
| B2 | Mixed host/target plot fixture; host is never densified | PASS |
| B3 | Original OOM regression test unchanged and passing | PASS |
| B4 | Empty comparison preserves “No aligned coverage” text | PASS |
| C1 | Native `test_gap_coordinate_parity_missing_annotation` | PASS |
| C2 | Oversized GTF raises error with both lengths and mismatch JSON | PASS |
| C3 | Missing annotation still draws depth, with explicit status | PASS |
| C4 | Missing target reports available references and target status | PASS |
| C5 | Multiple viral references require target selection | PASS |
| C6 | Native depths equal samtools depth fixture values | PASS |
| C7 | Transcript fixture suppresses incompatible genome track and labels status | PASS |
| C8 | Present zero-depth reference has zero-filled table and no_coverage status | PASS |
| D1 | Native STAR BAM fixture reviewed without invoking kallisto replay for that target | PASS |
| D2 | Below-gate candidate retained with detected=False unchanged | PASS |
| D3 | Missing STAR BAM reports none/not_assessed/alignment_bam_missing | PASS |
| D4 | Candidate-review fixture asserts primary table bytes unchanged | PASS |
| E1 | Planted same-target key: nondegenerate reproducible empirical null and p-value | PASS |
| E2 | Unrelated pool/whitelist pairs do not even load BUS | PASS |
| E3 | Same CB–UMI with a different target is not counted | PASS |
| E4 | Missing BUS gives not_assessed; native BUS preference/fallback tested | PASS |
| F1 | Sphinx HTML build; generated CLI parity and docs consistency tests | PASS |
| F2 | Stale design-only status replaced; no “design only” match in ambient plan | PASS |
| G1 | Protocol remains draft, both execution flags false; exact Git diff empty | PASS |
| G2 | Claims registry unchanged; no scientific performance claim added | PASS |

### Validation and qualification

- Baseline focused suite: 157 passed. Final focused suite: 292 passed, one
  harmless Matplotlib warning for the empty called-cell fixture.
- Broad suite before final typing/edge-case edits: 5,917 passed, 15 skipped,
  163 deselected in 739.39 seconds. Final serial suite on frozen code: **5,921 passed, 15 skipped, 163 deselected**
  in 696.51 seconds; no failures.
- Final selected integration suite: 15 passed. This includes existing
  pipeline-count, rerun-artifact, competitor-evidence, and evidence-chain tests.
- Ruff, mypy (62 source files), data governance, artifact inventory and built
  wheel/sdist/Docker-context ship-scope checks pass. Output-reference inventory
  digest and source commit were refreshed without changing claims.
- Sphinx uses the repository's configured theme and extensions, with missing
  docs-only packages installed in a temporary target; project dependencies are
  unchanged. Intersphinx network fetching and notebook execution were disabled.
- Coverage benchmark: one 8 kb viral reference and one million covered host
  positions, one process on an Intel Xeon E5-2690 v3 shared HPC host. Baseline
  peak RSS 795,848 KiB; candidate 766,828 KiB (0.964×; within the 1.10× limit).
  Covered-point data hashes and both coverage TSV byte hashes are identical.
- One 3 kb coverage map rendered in 1.819 seconds on that host. This is a host
  measurement, not qualification on laptop hardware or a pipeline speed claim.
- Existing real evidence BAMs were not found in the bounded results folders
  checked. WP-3's real-evidence-directory qualification remains **pending** even
  though C1–C8 and the CLI path pass on native synthetic BAMs.
- CI, publication, push, merge, training, holdout and scientific calibration were
  not run. WP-6/WP-8 remain deferred; WP-7 remains a read-only gate.

### Resolved specification discrepancies

- GTF features must fit within the BAM length; their maximum feature end need
  not equal the reference length. This follows the handoff's explicit bound
  check and avoids rejecting partial annotations.
- Recognized transcript accessions and whole-reference exon placeholders omit
  genome tracks before incompatible genome-span validation. CDS-only genome
  annotations are not automatically treated as transcript references.
- Zero coverage is the primary status for empty tracks; annotation status is
  retained separately. Invalid target/annotation inputs raise with status JSON.
- A sampled collision null with no variation returns not_assessed, not a
  fabricated p-value. The finite-sample upper tail uses the +1 correction.
- Insert/start–end confirmation is not implemented; it always remains explicit
  not_assessed. No swapped/infection interpretation is inferred from fingerprints.

Next action: review the six-package local delivery and obtain a permitted real
existing evidence directory for WP-3 qualification. This does not authorize any
scientific experiment, scheduler submission or remote publication.

## Kraken2 taxonomic discovery (DISC, new 2026-10-21)

Plan, audit and ADR: `docs/plans/2026-10-21-kraken2-discovery-plan.md` (branch `codex/kraken2-discovery-plan`).
Design: standalone `viralscan discover --run-dir R -o D` on STAR host-unmapped cDNA mates; separate fingerprint;
never writes the primary run, H5AD layers, `detected` or tiers. Snakemake leaf only on request (M3).

- [x] `DISC-00` — repository audit (ViralScan b8e5146; mg2sc a8c01f0 / fork e744c2c; scrna-workflow e0b9533;
  cellsnake 7adf99d; kraken2 2.17.2 8c190b1), critical review, ADR, contracts, gates, tickets.
  - Note: mg2sc/Cellsnake join BAM and Kraken2 output positionally (`k2sc.py:71`), vote with hash-order ties
    (`:119-121`) and count taxid 0 as a feature; none of it is reused.
  - Note: three barcode spaces exist (raw, kb/bustools-corrected, STARsolo with `--soloCBwhitelist None` without `-w`);
    `called_cells.tsv` may be in STARsolo space under host filtering. Spike `DISC-05` decides correction mechanics.
  - Note: the primary manifest does not hash `host_filtered/` outputs; discovery can only record its own hash.
- [x] `DISC-01` optional discovery env (kraken2 2.17.2), `doctor --profile discovery` (2026-10-21): `environment.discovery.yml`,
  pixi feature/env `discovery`, `DISCOVERY_TOOLS` (kraken2, k2, bustools; `FULL_TOOLS` untouched), CI integration job gains kraken2.
- [x] `DISC-02` `discover` parser + skeleton (2026-10-21): `src/viralscan/discovery.py`, `tests/test_discovery.py` (10).
  Refuses an incomplete run, an ambiguous or unknown `--sample`, `-o` inside the run or non-empty, an incomplete DB,
  and missing `kraken2`/`k2`. G2 test: source-run sha256 inventory unchanged. G1 test: the main parser gains no
  discovery option. After validation it exits 1 with "not implemented yet" until DISC-03..09.
  - Note: `--run-dir` is the run root (where `run_manifest.json` lives), plus `--sample`; the plan's example
    `--run-dir out/S1/` was wrong about the layout.
- [x] `DISC-03` input resolution (2026-10-21): reuse the primary `host_filtered/` pair only when the FASTQ pair count,
  the audit's `retained_host_unmapped` and the lineage's retained rows agree (else refuse); otherwise an explicit
  `--host-subtract starsolo --host-index DIR [--host-filter-star-params]` runs `_starsolo_filter` into
  `-o/.partial/host_subtraction/`, after checking the sample FASTQs exist and match the primary manifest sha256.
  The resolved input is staged as `.partial/input.json`; nothing is published yet. Tests: 21 in `tests/test_discovery.py`.
  - Note: `src/viralscan/discovery.py` added to `config/public_ship_scope.json` (wheel + sdist); DISC-02 broke
    `test_tracked_package_files_are_in_every_allowlist` without it.
- [x] `DISC-04` spike: DB-derived taxonomy = `k2 inspect --report-zero-counts` parsed by indentation (2026-10-21).
  - Note: NCBI rsync is refused from the login node; the test DB uses a hand-written 15-taxon `nodes.dmp`/`names.dmp`
    (reusable as the DISC-10 fixture). Kraken2 2.17.2 is in env `conda/envs/viralscan_discovery`.
  - Note: 2.17.2's confidence denominator includes ambiguous k-mers (`562:11 A:45` is unclassified at 0.5),
    contrary to the manual.
- [x] `DISC-05` spike: bustools fromtext→correct keeps flags (**done**, 2026-10-21); reconciliation on EBV
  SRR12682296 (**done**, 2026-10-21, Slurm `25860937`): `_replay_quant_bus` on the kb input pair with the virus-only
  index, 10xv2 on-list; 1,540,822 of 1,540,822 primary BUS records carry the same corrected CB as discovery's
  `barcodes.tsv.gz` (agreement 1.0, gate >= 0.999). Script and result: `viralscan_work/discovery_ebv_2026-10-21/disc05/`.
  - Note: kallisto 0.52.0 `bus` and `bustools fromtext` both encode `N` as `G` (CB and UMI), so discovery passes `N`
    through to match the primary barcode space instead of rejecting it; plan §6.3 updated. Use kb's bundled bustools
    (same 0.45.1 version as conda, different sha256).
- [x] `DISC-06` streaming pass (2026-10-21): `stream_reads` writes `kraken_input.fa.gz` (R2 as `>ordinal`, gzip, which
  Kraken2 2.17.2 reads natively) and `reads.tsv.gz` (`ordinal, read_id, raw_cb, raw_umi, geometry_ok, cb_has_n,
  umi_has_n, primary_read_filter`). It aborts on mate mismatch, truncation and duplicate read IDs (8-byte digest sort).
  The primary read-filter lineage is joined in lockstep (every ID checked) only for a reused primary `host_filtered/`
  pair; it is `not_joined` for a discovery STAR pass and `not_run` without `--read-filter`. 10xv1 and off-R1 geometries
  are refused before any STAR work. Tests: 36 in `tests/test_discovery.py`.
  - Note: every fragment is classified (ordinals cover `[0, n_input)`); bad barcode geometry is excluded later, by
    ordinal. Plan §5.3/§6.2 updated.
- [x] `DISC-07` Kraken2 wrapper + exactly-once ordinal join (2026-10-21): `run_kraken2` (list-form `kraken2 --db
  --threads --confidence --minimum-hit-groups --report --report-minimizer-data --output [--memory-mapping]`; never
  `--use-names/--paired/--quick`) and `join_kraken2`, which aborts on a duplicate, missing, out-of-range or
  non-integer ordinal, `C` with taxid 0 or `U` with taxid > 0, and writes ordinal-sorted `classifications.tsv.gz`
  plus raw `kraken.out.gz`. CLI: `--confidence` (default 0.05 since `DISC-R7b`, [0, 1]), `--minimum-hit-groups` (default 4 since `DISC-R7b`, >= 1),
  `--memory-mapping`. Tests: 49 in `tests/test_discovery.py`.
  - Verified end to end with real Kraken2 2.17.2 on the DISC-04 test DB through `viralscan discover`: 6 fragments
    gave 83333, 10376, 9606, 1280, U, 10376, matching their source genomes.
- [x] `DISC-08` LCA groups, called/off-cell partition, `cell_taxon.h5ad`; oracle + invariance tests.
  - [x] `DISC-08a` (2026-10-21): taxonomy from `k2 inspect --report-zero-counts` (`taxonomy.tsv` with
    `host_relation`); barcode correction with kb's bustools and the primary's on-list (`resolve_onlist`:
    `-w`, else kb's bundled list, else none), giving `barcodes.tsv.gz` (`onlist / corrected / uncorrectable /
    no_onlist / bad_geometry`); `umi_groups.tsv.gz`: exact (corrected CB, primary-encoded UMI) groups, taxon =
    LCA of classified reads, `R*`-rank LCA = ambiguous, CRC32-sharded and globally sorted. Gzip outputs use
    mtime 0, so they are byte-identical (G8). Tests: 100, including a 40-seed independent oracle checked byte-identical
    across shards 1/2/7 and permuted input.
    - Note: in the real-bustools run, a Hamming-1 variant of an on-list barcode with five on-list neighbours was
      dropped (uncorrectable), as bustools should; the designed truth table had expected a correction.
  - [x] `DISC-08b` (2026-10-21): `load_called_cells` reads the primary `results/called_cells.tsv` (suffix-stripped;
    refuses if absent) and records the cell-calling receipt's barcode space (`starsolo` when called on
    `GeneFull.raw`); `build_cell_matrix` writes `cell_taxon.h5ad` (obs = every called cell, zero rows kept; var = resolved
    taxids with name/rank/parent/`host_relation`; X = int32 CSR of resolved called groups; `uns.discovery.unit =
    resolved_umi_group`) and `offcell_taxon.tsv.gz`; matrix total = resolved called groups or abort (G6). Tests: 102.
    Real run (Kraken2 2.17.2, kb bustools, 10x v3 on-list): EBV 2, E. coli 1 (562), human 1 (`host`) in the one
    cell with evidence; ambiguous/unclassified groups excluded; two empty called cells kept.
- [x] `DISC-CMP` mg2sc/Cellsnake comparison (2026-10-21, user request): `analysis/discovery_mg2sc_comparison/`.
  Synthetic truth (1,200 UMIs; family-shared and EVE-like homology, collisions, novel genome); Kraken2 0.0/2 vs
  Cellsnake 0.05/4 vs fallback 0.01/3; ViralScan LCA vs audit deepest-consistent vs mg2sc vote. Results reproduce
  byte for byte. Settings change at most 8 UMIs. mg2sc vote: 11 % of labels change with the hash seed, 118/121 collisions
  are claimed as one species, 87 root features. LCA: 606 exact vs mg2sc 752; deepest-consistent: 799.
  - [x] Decision 4 changed by the user (2026-10-21) to deepest-consistent: see `DISC-R5`.
- [x] `DISC-R5` (2026-10-21, plan step R5; user decision 2026-10-21, replaces approved decision 4): per-UMI taxon =
  deepest classified taxon when all classified taxids lie on one root-to-leaf path, else their LCA; an `R*`-rank result
  stays `ambiguous`. `lineage_consistent` / `deepest_consistent_taxid` / `conflict` columns unchanged (audit).
  - G7 reworded: no group is counted below the LCA of its *conflicting* classified reads. Plan doc §6.4, §6.1/§7 tables,
    §14 (G7, case table) and §16 decision 4 updated; `_oracle` and the designed-case test updated (107 tests pass).
  - `compare.py` rerun (columns renamed `viralscan_deepest_consistent` / `viralscan_lca_only_previous`, the latter
    recomputed from `taxids`): exact species 799 vs 606 (LCA-only) vs 752 (mg2sc); collision overcalls 8 vs 2 vs 118.
    Outputs reproduce byte for byte (`comparison.tsv` `881a98e0…`, `umi_detail.tsv.gz` `7d03ef77…`).
  - Known limit: a species unit can rest on one species-specific read while the others are ancestral.
- [x] `DISC-R1` (2026-10-21, plan step R1): taxonomy path cache (root-first tuples memoised per tree, lazy per queried
  taxid); `lca` is a common-prefix walk, `host_relation`/`_group_shard` read the cache. No behaviour change: 40-seed oracle
  and a new all-subsets LCA test pass. Not re-timed on the EBV grouping step (needs Slurm; R10 records it).
- [x] `DISC-09` (2026-10-21, plan step R2): `provenance.json`, `qc.json` (all counts, called-cell overlap, fractions),
  `discovery_fingerprint` (sha256 of the output-determining subset), `discovery_complete.json` (output sha256s)
  published last; `input.json` removed; `discover` now exits 0. Same fingerprint on rerun is a no-op; a different one
  is refused unless `--force` (removes only what the marker lists). `--threads`/`--memory-mapping` are not
  fingerprinted. `db_identity` is the sha256 of the three `.k2d` files (+ a refused md5-manifest mismatch).
  Decisions: the fingerprint is computed before Kraken2 (so a no-op costs DB hashing only); `--host-subtract` inputs
  are identified by the manifest's sample fingerprints plus the index path (not its bytes). Taxonomy cache skipped
  (cheap after R1). A crash mid-publish leaves files without the marker and needs a manual clean.
  - Follow-up (2026-10-21, review minors): `fingerprint_inputs` is a deep copy, so `provenance.json` re-hashes to the
    marker (tested); resolved paths/counts moved to `provenance.resolved_input`; `provenance.discovery_fingerprint`
    added; `--force` refuses a marker that lists a non-local name (`../x`); `--minimum-hit-groups` covered by the
    rerun test. New `kb_input_check`: when the primary `kb count` read a `host_filtered/` pair at another path (a
    copied run), its sha256 must match the pair discovery classifies, else refused; the outcome is recorded.
- [x] `DISC-R6` (2026-10-21, plan step R6): `--polya-trim/--no-polya-trim` (default on). `trim_r2` cuts R2 from the first
  run of >= 15 A to the end, then a leading run of >= 15 T (constants `POLYA_MIN_RUN`/`POLYT_MIN_RUN`); applied in
  `stream_reads` before `kraken_input.fa.gz`. `reads.tsv.gz` gains `r2_trim` (none/polyA/polyT/both) and
  `r2_len_classified`; `qc.json` stream counts `n_r2_trimmed`/`n_r2_bases_trimmed`. Recorded in provenance and the
  fingerprint (a toggle needs `--force`). Decision: a read trimmed to nothing is sent to Kraken2 as `N` so it keeps its
  ordinal and returns unclassified (an empty FASTA record is not guaranteed to be emitted). R10 measures the effect.
- [x] `DISC-R7b` (2026-10-21, plan step R7b; user decision 2026-10-21, replaces approved decision 5): `--confidence`
  default 0.05 and `--minimum-hit-groups` default 4 (Cellsnake documented defaults); help text says they are not
  calibrated and gives Kraken2's own 0.0 / 2. Design doc §2.3/§4/§16 and `docs/cli_reference.md` updated. R10 still
  reports EBV at both 0.05/4 and 0.0/2. Tests pass explicit values, so none assumed the old defaults.
- [x] `DISC-10` (2026-10-21, plan step R4): `tests/test_discovery_integration.py` (`-m integration`; skips without
  `kraken2-build`, `kraken2`, `k2`, `bustools`). Builds the hand-taxonomy DB at test time (`--no-masking`, seeded, no
  network), runs `viralscan discover` end to end on a fixture run with real 10x v3 on-list barcodes at 0.5 / 2, and
  asserts `umi_groups`, `cell_taxon.h5ad` and `offcell_taxon.tsv.gz` equal designed truth (exact, parent/child ->
  deepest, conflicting -> LCA, root-LCA -> ambiguous, novel -> unclassified). Pins the confidence denominator: a read
  `562:11 A:45` is unclassified at 0.5 and classified at 0.05 (ambiguous k-mers count in the denominator). Passes
  in 27 s with kraken2 2.17.2. A one-substitution barcode was dropped from the fixture: bustools left it
  uncorrectable, which tests bustools rather than discovery.
- [x] `DISC-11` (2026-10-21, plan step R7): `viralscan discover-collapse OUT --rank {species,genus,family} [-o FILE]`
  (default `OUT/cell_taxon_<rank>.h5ad`). Each `cell_taxon.h5ad` column moves to its nearest ancestor-or-self at the
  Kraken2 rank code S/G/F via `taxonomy.tsv` (subspecies collapse to species); taxa with no ancestor at the rank stay
  as themselves with `above_rank = True`. No reclassification; row sums conserved (gate G6, tested and enforced at
  run time). Non-obvious: a taxon already at the rank is `above_rank = False`; the default output sits inside the
  published dir but is not listed in `discovery_complete.json`.
- [x] `DISC-R8` (2026-10-21, plan step R8): `candidates.tsv` (observation state only), written by `discover` before
  provenance and published (listed in `discovery_complete.json`). One row per taxon with a resolved group: `groups_called`,
  `groups_offcell`, `n_cells`, `reads` (Kraken2 report own reads), `distinct_minimizers`, `minimizer_per_read`,
  `host_relation`, `db_has_9606`, `state = observation`, 0/1 flags `low_minimizer_diversity`, `offcell_dominant`,
  `host_lineage`, `viral_only_db`. Thresholds are documented constants, not calibrated: `LOW_MINIMIZER_DISTINCT = 20`,
  `OFFCELL_DOMINANT_FRACTION = 0.5`.
  - Note: the low-diversity flag uses the absolute distinct-minimizer count, not the ratio, because a deeply
    sequenced small genome has a legitimately tiny distinct/reads ratio (EBV: ~600k reads). The ratio is informational.
    `host_lineage` covers `host` and `host_ancestor`; `viral_only_db` means taxid 9606 is absent from the DB taxonomy.
    No `prioritized`/`sequence_supported` states (need the frozen benchmark, `DISC-12`/`DISC-13`).
- [x] `DISC-R9` (2026-10-21, plan step R9): `docs/discovery.md` (scope: classification is not infection; manual database download with md5
  check, no auto-download; command and Slurm example on partition `all`, since `medium` is capped at 4 CPUs / 8 GB by QOS
  `restrictmedium`; outputs; `qc.json` guide; limits), toctree entry in `docs/index.md`, `output_reference.md` section, CHANGELOG entry.
  - Note: `docs/output_reference.md` is pinned by inventory row 17 (`output-reference-doc`); re-pinned in a follow-up
    `chore(inventory)` commit because the pin names the commit that holds the new file.
  - Note: the only real-data numbers cited are labelled "first EBV run, read level" and "not a performance claim"; R10 replaces them with
    the validated arms.
- [x] `DISC-R10` (2026-10-22, plan step R10): EBV SRR12682296, PlusPF-16 2026-06-26, code `c41e330`, Slurm array
  `25862538` (partition `all`, 8 CPUs). Arms: default 0.05/4 + trim; 0.0/2 + trim; 0.05/4 `--no-polya-trim`. All three
  published `discovery_complete.json`; the same-input rerun logged "nothing to do" (7.7 min, 0.2 GB: DB hashing);
  the source run's sha256 inventory (76 files) was unchanged before/after (G2). All 2,524 called cells seen.
  Per arm (default / 0.0-2 / notrim): unclassified 93.4 / 71.3 / 92.9 %; human clade 4.1 / 20.9 / 4.5 %; EBV clade
  2.3 / 6.8 / 2.3 %; EBV groups in called cells 246,233 / 623,096 / 246,882; EBV+ cells 2,224 / 2,471 / 2,227 (primary
  2,524); per-cell Spearman vs primary EBV molecules 0.890 / 0.966 / 0.891; median ratio 0.26 / 0.82 / 0.27;
  candidates 47 / 1,906 / 56 (low-minimizer 26 / 1,680 / 34). Run 1: 35-37 min, 15.2 GB. Trim touched 2.03 M reads.
  Outputs, inventories, `summarize.py` and `r10_summary.tsv`: `viralscan_work/discovery_ebv_2026-10-21/r10/`.
  - Note: the defaults (user decision R7b) trade ~60 % of EBV groups for a 5x smaller human fraction and a 40x
    shorter candidate list; documented in `docs/discovery.md` with the 0.0/2 alternative. The user kept 0.05/4 on 2026-10-22.
  - Note: R1 grouping was not timed separately (the log has no timestamps). Whole runs took 35-37 min (the pre-R1
    run, 26 min, stopped before grouping finished and before publication; ~8 min is now DB hashing). No grouping
    speed-up is claimed beyond the runs completing. Taxonomy cache still skipped.
  - Called-cells discrepancy (2,524 vs 2,764) explained earlier: `viral_summary.tsv` uses `host_called_cells.tsv`;
    discovery uses `called_cells.tsv`, as intended. The Drop-seq `no_onlist` library (R10.5) was skipped: that run
    has no `called_cells.tsv`.
- [x] `DISC-R11` (2026-10-22, plan step R11): merged `origin/main` (CLI-UNIFY-01) at `5ff8011`, inventory re-pinned;
  full suite in 3 foreground chunks 6,037 passed / 0 failed; `pytest -m integration tests/test_discovery_integration.py`
  2 passed (discovery env); `ruff check .` / `ruff format --check .` clean; `gen_cli_reference.py --check` clean;
  `doctor --profile discovery` ok in the discovery env, fails on missing kraken2/k2 in the main env. "(In development)"
  labels dropped. Branch handed over for push / PR.
  - Note: PR #89 CI lint failed on `mypy src/viralscan` (not run locally before): 7 annotations in `discovery.py`
    and an `out` name reused as Path/str in `menu.main`; fixed, `mypy src/viralscan` clean locally.
- [x] `DISC-CELLS-01` (2026-10-23, `codex/discover-called-cells`, decision Q20a): `discover --called-cells FILE` overrides the
  sample's `called_cells.tsv` (two-step runs call 5-455 cells in uninfected/low-viral libraries vs 3,515-30,163 by EmptyDrops on
  the full matrix). Space from the file's own `cell_calling_input.json`; path in provenance, sha256 in the fingerprint; an empty
  override raises. Evidence producers still read the run's own cells. Test `test_load_called_cells_override_uses_its_own_receipt`.
- [x] `DISC-12` frozen benchmark protocol linked to `HOST-05` / `SW-07`. Draft (2026-10-22, `codex/disc12-protocol`):
  `docs/plans/2026-10-22-discovery-benchmark.md`, not frozen, nothing run. Decisions recorded 2026-10-22 (Q1–Q6, Q12):
  `k2_viral_20260626` vs PlusPF-16; HHV-6B library by max primary unique molecules; HSV-1 Drop-seq excluded (`not_assessed`);
  spike-in truth from array 25864985 if per-read; §5 rules filled. 2026-10-23: all inputs decided (Q3', Q16-Q20); input
  manifest pinned (`docs/plans/disc12/input_manifest.tsv`); runs need `discover --called-cells` (PR #112). **Frozen 2026-10-23** (user approved PR #103; `discover --called-cells` merged as #112). Next: the benchmark run.
- [~] `DISC-13` candidate states + evidence adapter (2026-10-22, Q12): `build_candidate_states` / `assign_state` / Route (ii) rejects kb/star sources (2026-10-22 review).
  `read_evidence_table` in `discovery.py` write `candidate_states.tsv` (`observation`, `corroborated` route i,
  `read_supported_observation` route ii; `corroboration` / `read_support` = `not_assessed|negative|positive`; viral-only DB
  and host lineage capped at `observation`). Pure functions, 4 new tests; `docs/discovery.md` "Candidate states".
  - Note: the package has no per-taxon alignment summary, so route (ii) and (i) run over a documented evidence table
    (`EVIDENCE_COLUMNS`); cut-offs 3 units / 80% / 1 shared cell are uncalibrated constants.
  - Done 2026-10-22: `build_kb_evidence` (route i, `source=kb`) from `molecule_assignments.tsv.gz` + `virus_identity.tsv`;
    API only (no `discover-evidence` CLI: needs parser, menu dispatch and help wiring). Tests `tests/test_discovery_kb_evidence.py`.
    Review fixes: empty/no-called-cell molecule file raises; counts are distinct (CB, UMI); above-species ancestor genes compete
    (`n_adjudicated`) but never make a molecule taxon-best.
  - Done 2026-10-22: route (ii) producer `build_read_support_evidence` (`source=minimap2`, `evidence_read_support.tsv`), API only
    (no CLI, plan §7). Per `--candidate-fasta` taxid (cap `max_candidates=3`): exact R2 reads of the taxon's resolved called-cell
    groups -> `HOST|` + `VIRUS|CAND_*` (+ optional `VIRUS|PANEL_*`) -> alignment (injectable aligner/hit reader). Unmapped reads
    count in `n_adjudicated` (tested, no support); taxon-best needs a strictly higher AS on candidate contigs than on any
    host/panel contig (ties, host/panel wins are not best). No resolved called-cell read -> no row (`not_assessed`). Read IDs
    of the FASTQ pair are checked against `reads.tsv.gz`. Tests `tests/test_discovery_read_support.py` (stubbed aligner; no
    real-minimap2 test: not installed in the dev env). `coverage_table`/QC tables of plan §9.3 not emitted (evidence row only).
  - Done 2026-10-22: route (i) STAR producer `build_star_evidence` (`source=star`, `evidence_star.tsv`), API only (no CLI). Reads
    pruned `star_hostvirus/viral.bam` (injectable record reader), `VIRUS|<acc>` -> taxid via `virus_identity.tsv`, kb lineage rule.
    Molecule = distinct (CB, UB), `-N` CB suffix stripped; `-`/missing CB or UB dropped (Unique-mode UB: lower bound). NH rule:
    a read is taxon-best only if visible T records == `NH` (host/unseen loci invisible after pruning); molecule best when all
    its reads are. Raises on missing files, no molecule, or no called-cell overlap. Tests `tests/test_discovery_star_evidence.py`.
  - TODO: STAR route (i) producer (STAR confirm lives on an unmerged branch).
  - Left: no CLI (plan §7 specifies none); `prioritized` / `sequence_supported` of plan §9.2 not built; cut-offs
    need `DISC-12`.
- [ ] `DISC-14` optional Snakemake leaf (only on request).

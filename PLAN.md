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
  - 2026-10-23: tag `v3.0.0rc1` (on `87591c4`) pushed; release run 40000132978 built and `twine check`ed both
    distributions, then failed at `check_ship_scope.py` (`jsonschema` missing in the release job; nothing published).
    Fix: release.yml installs `jsonschema`. Re-tag needed (the release runs only on a tag push of a CI-green commit).

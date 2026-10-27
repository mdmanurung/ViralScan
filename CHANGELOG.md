# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased (after 3.0.0rc1)

- Added `viralscan host-loss`: D15/D16 truth-read loss at the STAR host filter from `fragment_lineage.tsv.gz`, reported as a lower bound under `posthoc/host_loss/` (PLAN `HOST-05a`). `--star-confirm` now keeps `star_hostvirus/star_read_fates.tsv.gz`, per-read STAR fates written before the prune (PLAN `STAR-FATE-01`).
- Added the kallisto translated search as an exploratory post-hoc arm (AA-11 oracle in the reproducibility repo reproduces the frozen arm P summaries for three chemistries): `viralscan translated-index` (annotated viral proteins + PalmDB RdRP, `kb ref --aa`, plain host D-list, kallisto >= 0.51.1 required) and `viralscan translated` (`kb count --aa` on a completed run, per-virus/protein UMI summaries under `posthoc/translated/<sample>/`). Never changes `detected` or tiers (PLAN `AA-08`..`AA-12`).
- `viralscan run --discover --kraken2-db DIR` runs Kraken2 discovery on every sample once the run is complete, into `posthoc/discovery/<sample>/` (PLAN `DISC-14`).
- `viralscan discover` defaults `-o` to `<run-dir>/posthoc/discovery/<sample>`; an explicit output may sit under the run's `posthoc/` folder (PLAN `DISC-16`).
- Added `viralscan discover-states`: candidate states for a published discovery from the run's own STAR or kb evidence (optional minimap2 read support), plus the Kraken2-only `prioritized` state (`--k`, default 3, uncalibrated). Never changes `detected` or tiers (PLAN `DISC-15`).
- CI runs the real-STAR tests (`star` in the integration environment) and a new end-to-end `run --host-filter starsolo --star-confirm` test on the tiny fixture; `star_hostvirus.star_index` accepts `sa_index_nbases` (PLAN `STAR-05`, `STAR-06`).
- Docs: the site is cut to 13 sidebar pages. `docs/results.md` replaces `inspect_results.md` + `assess_evidence.md`; `api_reference.md` is folded into `api.md`; `reference_help.md` and `project_info.md` are removed. `v3_counting_contract.md`, `panel_integrity.md`, `workflow_comparison.md` and `reference_maintenance.md` leave the sidebar but keep their paths (PLAN `DOCS-REV-02`).

- `viralscan discover --called-cells FILE` replaces the run's own `called_cells.tsv` (e.g. the EmptyDrops cells of a full-matrix run when a two-step run called only virus-carrying barcodes). The file's barcode space comes from its own `cell_calling_input.json`; its path and sha256 are recorded in `provenance.json` and the sha256 enters the fingerprint.

## 3.0.0rc1 (release candidate; date set when tagged)

- Added `viralscan discover` and `viralscan discover-collapse` (PLAN `DISC`): Kraken2 classification of a completed run's host-unmapped reads, grouped by cell barcode and UMI into `cell_taxon.h5ad`, `candidates.tsv` (observation state only), `qc.json` and `provenance.json`, published atomically with a fingerprint and a no-op rerun. Classifications are not infection calls; the primary counts, `detected` and tiers do not change, and no discovery performance claim is made. Optional `environment.discovery.yml` and `viralscan doctor --profile discovery`. See `docs/discovery.md`.

- **Breaking (CLI-UNIFY-02, layout version 2):** post-hoc verbs write under `<run-dir>/posthoc/` by default
  and never modify the run's own outputs. `evidence`: `--output` optional (`posthoc/evidence/<virus>`);
  `rerun-multimap`: `-o` optional (`posthoc/multimap_<method>`, still must be empty); `rerun-programs`: writes
  `posthoc/programs/` (no longer in place, no marker re-stamp); `hostresponse`: `posthoc/hostresponse`;
  `evidence --review-candidates`: `evidence_review.tsv` moves from `results/` to the evidence output
  directory. `run_manifest.json` records `layout_version: 2`; `validate-run` accepts 1 and 2 and ignores
  `posthoc/`. See `docs/migration.md`.

- The evidence step writes scPathoQuant-style coverage maps (`<ref>_coveragemap.<layer>.{png,pdf}`, via `bamCoverage`, `gtf2bed` and `pyGenomeTracks`) plus a read-competition track figure and a stacked histogram, replacing `coverage.raw_vs_deduplicated.png`. New tools in `environment.coverage.yml` (optional; `viralscan doctor --profile coverage`): `pygenometracks`, `deeptools`, `bedops`. The read-competition classification also works on STAR/STARsolo BAMs (no `SA` tag): partners come from the other records of `NH`>1 reads (`samtools view -N`).
- **Breaking:** `viralscan hostresponse -o RUN_DIR` is now `--run-dir RUN_DIR` (it was an input run
  directory, unlike `-o/--output` elsewhere). `-o` errors with the replacement; `--output DIR` is optional and
  names where hostresponse writes (see the layout-version-2 entry above for the default).

- Added `viralscan.star_hostvirus` (PLAN STAR-HOSTVIRUS-01): builds one STAR index of the host genome
  (`HOST|` contigs) plus the viral panel (`VIRUS|`), and aligns host-unmapped read pairs with STARsolo
  (multimappers kept, CB/UB tags) to a coordinate-sorted, indexed BAM for coverage maps. Run from the
  workflow with `viralscan run --star-confirm` (PLAN STAR-04, below).
- New `viralscan star-index` verb (PLAN STAR-03) builds the host + viral-panel STAR index and writes
  `star_hostvirus_index.json` (input sha256, STAR version, arguments). `viralscan run --star-confirm
  --star-hostvirus-index DIR` (PLAN STAR-04; needs `--host-filter starsolo`) aligns the host-filtered reads
  to it and adds `star_reads_primary`, `star_reads_any`, `star_molecules`, `star_host_competing`,
  `star_other_virus` to `viral_summary.tsv` plus `results/star_support_by_virus.tsv`. Labels only: no call
  changes. STAR can miss short fragments, drops reads over 100 loci (HHV-6A/6B/7 telomeric repeats, so `star_host_competing` is a lower bound) and does not wrap circular genomes; the
  spike-in detection-limit experiment (SPIKE-01) has not been run.
- `star_hostvirus` hardening (PLAN STAR-02): barcodes are corrected against the chemistry's on-list by
  default (raw barcodes only for list-less chemistries); `--strand forward|reverse|unstranded` maps to
  `--soloStrand` (it was fixed to Unstranded); STAR 2.7.11b's exit-139 zero-mapped crash is recovered through a
  shared `viralscan.star_run` (also used by `anello_align`); new per-virus summary
  (`star_reads_primary/any`, `star_molecules`, `star_host_competing`, `star_other_virus`) and an optional
  prune to a VIRUS|-only BAM.
- `evidence-plot --coverage-map` now accepts a bare GTF seqid `ACC` for a viral BAM
  reference `CLASS|ACC` (host contigs are never prefix-stripped) and records the matched
  seqid and whether it was normalised in `coverage_map.<target>.json`.
- SuperExpressor summaries now require a complete called-cell denominator. In the
  fixture with 1,000 called cells, 50 viral-positive cells and two reaching 10
  allocated molecules, the derived percentage changes from **4.0% to 0.2%**.
  Explicit legacy subset mode remains available with warning/provenance. Empty
  summaries retain the full header. Primary counts and detection defaults do not change.
- High-burden PNGs use called cells and host-only depth, with exact TSV/JSON
  sidecars, explicit missing-null status, display-floor wording, and labelled
  figure-only truncation. Coverage comparisons preserve uncovered viral bases.
- Added opt-in accession coverage maps, source-labelled candidate review, and
  pool-scoped tier-2 molecule fingerprints. These are diagnostics, with no
  scientific performance claim or change to the frozen validation gate.

---

## [Unreleased]

### Added
- **Research-tree size guard (RESTREE-01).** `.gitignore` rules keep bulk analysis outputs untracked, and `scripts/check_tracked_file_sizes.py` (run by `tests/test_tracked_file_sizes.py`) fails if a tracked or staged `analysis/` or `results/` file exceeds 1 MB, unless listed in `config/large_file_allowlist.txt`.
- **Near-duplicate screen and host-homology warning at reference build (AUD3-05, P2-01).** `validate_reference_records` now also flags records that are a near-copy of another record on either strand (under 5% of the record's canonical 31-mers outside the others, the existing `kmer_twin_of` cut, now in `viralscan.kmer_sharing`); a warning by default, a build failure with the new `--strict-reference` on `build-ref` and `scripts/build_bundled_panel_ref.py`. Declared catalogue siblings are exempt. A diverged (about 10%) homologue is not flagged: its threshold is PLAN REF-08. A build without `--genome-dlist` now logs that host homology was not measured. Exit status is unchanged by default.
- **Probable/strong evidence-tier code (TIER-01).** `viralscan.evidence_tiers` is a pure module built to the frozen `calibration.threshold_search` contract of the v3 protocol: it reads the 45-point cutoff grid from the protocol, selects the most conservative eligible point within one standard error of the best training F1, and returns no tier when no grid point has zero host-only and planted-homology false calls. Without a supplied calibration result the tier is `uncalibrated`; nothing assigns `probable` or `strong` otherwise. Not wired into any output yet, and no calibration was run. `scripts/adjudicate_candidates.py` appends a `package_tier_analogue` column (L0 candidate, L1 probable, L2 strong, L3 strong); package tiers and study claim levels remain separate scales.
- **Discovery-study tooling (Track B, ADJ-01, LONG-01/02, CAL-13, LOCK-02).** Four checkout scripts for the broad virus discovery study; none ships a numeric threshold, the values come from a calibrated protocol lock. `scripts/adjudicate_candidates.py` assigns claim levels L0-L3 per library and virus from the DSR-02 read evidence and the unique layer (host expression is never read; anelloviruses cap at L1, non-herpesviruses at L2; without a lock every call is L0). `scripts/longitudinal_compile.py` builds per-subject depth-normalised trajectories and event labels and demotes a new detection or burden increase whose molecules or loci appear in another subject's library in the same run or lane. `scripts/cal13_decoys.py` makes seeded 6-let-preserving shuffled decoys and scores calls against them (stratified empirical p, target-decoy and BH q, expected false calls per library). `schemas/v3/protocol_lock.schema.json` and `scripts/validate_protocol_lock.py` validate `protocol-v1.lock.json` and re-hash every file it cites.
- **TusslChemistry (`tchem`) second opinion on the chemistry (CHEM-02).** `viralscan` and `viralscan check-chemistry` run [tchem](https://github.com/Tussl-Co/TusslChemistry) (v0.1.0, a Rust binary the user installs; barcode lists via its `scripts/fetch_whitelists.sh`) beside the built-in detector and compare barcode+UMI geometry. Flags `--tchem {auto,require,off}` (default `auto`: used when a binary is found, skipped otherwise), `--tchem-bin`, `--tchem-whitelists` (also `$VIRALSCAN_TCHEM`, `$VIRALSCAN_TCHEM_WHITELISTS`, `$TUSSL_CHEMISTRY_WHITELISTS`). A disagreement with `-x`, or a tchem `ambiguous`/`no_match`, stops the run unless `--force-technology`. The call is recorded under a `tchem` key per sample in the manifest. Checked against the real binary on 80 Track A libraries.
- **Truth-panel generator (VAL-01/02/03/05, tiny scale).** `scripts/generate_truth_panel.py` writes seeded paired FASTQs with truth manifests for the v3 validation panel (synthetic background; `--generation-seed` required; `--production` refuses until the protocol is frozen). `v3_partition.unit_id`/`split_libraries` add the donor-library split. Code only: no outcome run.
- **Truth-panel generator takes real references and barcodes (AUD-18).** `--viral-fasta` + `--gene-spans` (and optional `--sibling-pairs`) plant reads from real genomes and strand-aware gene models, with `unique`/`ambiguous` truth following those models; `--onlist` draws all barcodes from a real 10x on-list. Input hashes go under `inputs` in `run_manifest.json`. Without the options the output is byte-identical to before. See the generator design doc, section 8a.
- **Truth-panel generator: CB-preserving hops (M04-GEN, Arm X).** `--hop-preserve-cb` makes a hopped molecule keep the donor `(CB, UMI)` instead of taking a random uninfected barcode of the recipient. `truth_manifest.tsv` then gains `hop_donor_sample_id`, `hop_donor_molecule_id`, `hop_cb_preserved` and `hop_landing` (`cell`, `empty` or `other`, from the recipient's barcode pool), and `run_manifest.json` gains `hop_mode`. Off by default: without the option the output is byte-identical to before. See the generator design doc, section 8b.
- **Truth-panel generator: real-background `D1` mode (M04-GEN, Arm P).** `--background-r1/--background-r2` (+ `--background-chemistry`, `--n-background`) take a real 10x FASTQ pair as the background: each planted `synthetic_factorial` sample is the library's reads verbatim (names rewritten, qualities kept) plus planted molecules, and every one has an unplanted twin made of the same reads with nothing added. Background reads carry truth `unknown`, planted molecules go into the library's most frequent barcodes, and `twin_sample_id` / `background_library_id` pair the arms. `plant_anello_10x.py` gains `iter_background`, which `copy_background` now uses (its output is unchanged). Without the options the output is byte-identical to before. See the generator design doc, section 8c.
- **Scorer: hopped/ambient molecules in a negative reported three ways, truth tiers T1-T5 (M04-SCORE).** `v3_score.score_negative_libraries` gives the unconditional any-virus false call (the headline, per-virus rates alongside), the rate after excluding libraries whose call is fully explained by documented hop/ambient molecules (with their provenance), and the contamination-attributable share, always together. Tiers T1-T5 are validated and never pooled (mixed tiers are refused; T5 is not scored); the scorer takes the tier as input and assigns none. `origin = unknown` background is accepted and never truth, `hop_*` columns are validated, and D1 manifests are exempt from requested == realised. No acceptance threshold.

- **Per-marker polyA capture fields in the programme catalogue (AUD-17, M-02).** `gene_programs.tsv` gains `polya_capture` (`yes`/`no`/`unknown`) and `capture_basis`. `yes` is annotation-implied from RefSeq, `unknown` is the default; HSV-1 `LAT` is `unknown` (PMID 38445887). A virus with such a marker gets a `caveat`; a productive one makes `reactivating`/`latent` read `indeterminate`. No current call changes (EBV and HCMV markers are all `yes`). The competitive host+virus reference (M-05) is a design note only: `docs/plans/competitive-host-virus-reference.md`.
- **Ambient-RNA screen and cross-run barcode-swap check (M-03).** `viral_summary.tsv` gains nine descriptive columns (`ambient_viral_fraction`, `ambient_viral_fraction_ub95`, `n_ambient_barcodes_with_virus`, `soup_expected_cells_max`, `soup_expected_cells_ge2_max`, `soup_tail_p_cells`, `viral_depth_spearman`, `ambient_assessment`, `ambient_basis`): the calibration-free `rho = 1` soup bound with a Poisson upper confidence bound, `not_assessed` (blank, never 0) for virus-only indexes, `--cell-calling none`/`knee` and runs with no host depth. `detected` and every tier are unchanged. Two-step runs carry per-barcode STARsolo host totals in `host_cells`. `scripts/swap_check.py` tests barcode overlap between libraries a user-written pool manifest declares to share a pool and whitelist. The `contamination` row of `interpretation_flags.tsv` shows the screen but stays `not_assessed`. Doublets and uptake are out of scope (claim limit in `docs/output_reference.md`).

### Changed
- **One CLI grammar; quantification is `viralscan run` (CLI-UNIFY-01, breaking).** `-s1/-s2` become `--r1/--r2`; the single-dash long flags (`-ref`, `-gtf`, `-fasta`, `-f1`, `-umap`, `-acc`, `-mm`, `-v`) and the `--visual False` value form are removed; `--data-cache-dir` becomes `--cache-dir`, run `--cell-types` becomes `--cell-type-labels`, `--hostresponse-control-mito` becomes `--hostresponse-mito-control`, `hostresponse --hostresponse-min-umi`/`--detection-threshold` become `--min-umi`, `validate-run` takes `--run-dir` and `--json PATH`, and `doctor --json` takes a path. Removed spellings fail with the replacement named; prefix matching is off on every subcommand (`allow_abbrev=False`), so an abbreviation that used to resolve, such as `--coverage` for `evidence-plot --coverage-map`, is now rejected: spell the full option. A bare `viralscan [OPTIONS]` still runs `run` with a deprecation warning. Config keys and `run_fingerprint` are unchanged. Table: `docs/migration.md`.
- **Detection and multimapping stop when the Virus Identity table does not match the kb matrix (AUD3-01a, audit P9-02).** A matrix column the table does not list (for example the table lists `vB_g1`, the matrix holds `vB_g1.2`) used to make a real detection vanish without a row, line or warning on `--cell-calling none`, an external list or two-step mode; it now raises an error asking for the analysis step to be re-run. Viral table genes that are absent from the matrix log one warning. Thresholds and denominators are unchanged.
- **Detection gate no longer loses a molecule split 6 or 7 ways to float rounding (AUD3-01b, audit P9-03).** `detect_genes` and the called-cell gate compared summed viral molecules with a bare `>=`, so six or seven copies of 1/k that sum to just under 1 in float64 dropped the virus with no row and no CANDIDATE ONLY. Both now use the same 1e-6 relative tolerance as the infected-cell rule. Thresholds, layers and denominators are unchanged; a mass of 0.99 is still below the gate.
- **UMAP neighbour-enrichment test runs on a viral-free graph (AUD3-02, P5-01/P5-07/P8-01).** The permutation p-value in the `--umap` annotation is computed on a host-only kNN graph (host size factors, HVGs without viral genes, PCA), not on the viral-inclusive UMAP whose embedding contains the genes that define the labels (random labels gave p < 0.05 in 15/15 runs). The UMAP figures are unchanged apart from being labelled viral-inclusive; the annotation states what is tested and that it does not distinguish infection from ambient RNA, uptake or doublets, and the report glossary no longer says clustering "suggests true infection". Fewer than two highly variable genes now logs a warning and skips instead of failing in PCA. The test null is still not calibrated (open item `AUD3-02-OPEN`).
- **Host-response default is `--cv auto` (AUD3-06b, audit P5-02).** `viralscan hostresponse`, `run_hostresponse` and the pipeline (`--hostresponse-cv`) now default to `auto`: group-disjoint CV on `--groups`, or on the one donor/sample-like obs column if exactly one exists (logged), otherwise cell CV with a loud donor-leakage warning (also when group CV cannot run on the auto-selected column; reason in the manifest as `cv_auto_fallback`). `--cv cell` is the explicit legacy mode and `--cv group` stays fail-closed; both are byte-identical to before for a fixed seed. Rows keep the concrete `cv_mode`; the manifest adds `cv_mode_requested`. Single-sample inputs (e.g. the EBV-matched LCL analysis) still run cell CV unchanged, now with the warning.
- **Host-response warnings; selection frequency relabelled uncalibrated (AUD3-06a).** Warnings, help text and docs only: no metric, default, file or column name changes. `hostresponse` now warns when `--groups` is given without `--cv group` (it was silently ignored) or a donor/sample-like obs column goes unused under cell CV (donor leakage possible), when a cell-CV test fold holds fewer than 5 positives (one aggregated warning per virus; with about 10 positives and the default top-50%-depth filter each test fold is 1 positive x 1 negative), and logs the skipped host normalisation at WARNING (the docstring threshold now agrees with the code: max > 10). The randomized-Lasso `stab_prob` / `stable` columns are described as an uncalibrated selection frequency (C is fixed, so noise genes pass `--stab-min-prob` as positives grow: 0 / 0 / 1.5 % / 25 % of noise genes at 30 / 100 / 400 / 1000 positives, 1500 features) and a warning is logged above 300 positives. README/PUBLICATION_READINESS host-response numbers are marked cell-level CV.
- **`em-global` residue, tier slack and gene-less ECs documented and pinned (AUD3-07).** `docs/v3_counting_contract.md` records that `em-global` (expert value, not the default) can stop after two iterations at tol 1e-6 and leave about 2.6 allocated molecules on a zero-unique-support viral gene at 1e6 unique / 4e6 ambiguous host molecules (it does not reach `detected=True`), the 1e-6 relative slack on tier cut-offs, and that gene-less ECs are skipped in the molecule intersection. No numerical change; `tests/test_em_global_residue_pinning.py` pins all three.
- **Combined-mode cell calling is host-only and fails closed on a virus-only index (CELLS-02, CELLS-05).** EmptyDrops sees the host gene columns of the kb matrix, never viral reads. A run whose viral gene partition is missing, matches no column, leaves no host column (a virus-only index: use two-step `--host-filter starsolo` or `--called-cells-file`) or whose Virus Identity table omits a matrix column now stops with a `CellCallingError`. The `--emptydrops-seed` value 20260727002 (the frozen `seeds.cell_calling`) used to stop `emptydrops.R` (it overflows an R integer); it is now reduced modulo 2147483647, other seeds are unchanged, and `cell_calling_input.json` records `seed` and `effective_seed`, the R and DropletUtils versions and the partition file digest. Called sets and `knee`/`inflection`/`total` in `emptydrops_cells.tsv` differ from earlier host-plus-viral calls.
- **`--multimap-method` shows three options: `conservative` (default), `em`, `equal` (Q72-Q74).** They are CLI aliases resolved at parse time (`conservative` = `host-conservative`, `em` = `em-cell`); manifests, configs and outputs record only the canonical name. `unique-weighted`, `sibling-weighted`, `em-global` and the canonical names stay accepted, unlisted. Same for `rerun-multimap`. `em` = `em-cell` because the EMC-01 re-measure passed (17.6 GB, 1 h 23 min on 942,700 barcodes).
- **Infected cells and detection use molecule counts and called cells (AUD-06/AUD-07).** A cell is infected for a virus when its summed allocated molecules reach `--infected-min-molecules` (default 1.0); fractional-only cells are counted as sub-threshold. A virus is `detected` only with at least `--detection-min-cells` (default 1) infected called cells and `detection_threshold` molecules in called cells; the all-barcode candidates stay in `viral_summary.tsv` with `detected` False, and `background_umis` reports empty-droplet molecules. `summary.txt` says "candidate". Counts can drop relative to earlier versions. New columns: `infected_subthreshold_cells`, `infected_subthreshold_called`, `background_umis`, `detected`.
- **Host-response per-cell label has its own threshold (B-16).** New `--hostresponse-min-umi` (default 10); the `hostresponse` subcommand's `--detection-threshold` is deprecated.
- **`lod95_molecules` now means required true molecules (SENS-CORR-03).** It was the
  expected *observed* molecules (2.996 at depth 1e6, capture 0.5); it is now the true
  molecules a virus needs for 95 % detection (5.991), equal to `depth * lod95_per_10k / 1e4`.
  The old quantity is the new `lod95_expected_observed_molecules`. `sensitivity.tsv` gains
  `detected`, `capture_status`, `depth_only_sufficient`, `sensitivity_eligible` and
  `negative_blockers`; `informative_negative` is now `sensitivity_eligible and not detected`.
- **Diagnostic molecule representatives (EVID-CORR-01).** `viralscan evidence` keeps one
  deterministic read per corrected (CB, UMI, reference) molecule for BAM deduplication and
  the read-start profile, so one molecule aligned at two starts is no longer counted twice.
  Evidence tables from earlier runs change in the `umi` layer.

### Fixed
- **summary.txt names the two-step barcode set (AUD3-03 follow-up).** In a two-step (`--host-filter starsolo`) run the per-virus line and the headline said "all barcodes" for what is the host-filtered barcode count; they now say "host-filtered barcodes". Combined-mode text is unchanged.
- **Report headline and two-step denominators (AUD3-03, P3-03/P4-02).** `report.html` headlines the called-cell rate (`pct_infected_called`, comparable rate when there are no called cells) in the run summary and main table, with the all-barcode rate secondary and the badge on the primary rate. In two-step (`--host-filter starsolo`) runs `pct_infected` and `viral_molecules_per_10k_est` are NA (they divided by host-filtered barcodes: 100.0 / 10000.0), with a new `denominator_note` column in `viral_summary.tsv`; combined-index values are unchanged.
- **Input robustness (AUD3-09a, P3-06/P4-04/P4-06/P4-03).** `read_support` skips a molecule whose viral-aligned primary record has `SEQ=*` (one aggregated warning) instead of raising `IndexError`. `check-whitelist` strips a Cell Ranger `-<n>` GEM-group suffix before matching, so a suffixed list is no longer reported as a wrong chemistry. An explicit `bc:umi:seq` `-x` with a file index or offset the host-filter/read-filter/evidence paths cannot honour is now an error, not silently read from R1 base 0 (named technologies and plain `kb count` unchanged). The `knee_cells` docstring now states its real failure on realistic curves; behaviour unchanged.
- **Chemistry sanity gate (AUD3-09b).** A virus-only index without `--host-filter` no longer raises a `p_pseudoaligned` error on every normal run (it holds no host transcriptome, so a low rate is expected); combined indexes keep the error. An explicit `--strand forward` on a host-filtered (two-step) or virus-only run now logs a `strand` warning in `chemistry_sanity.json`, since on a 5' library it silently loses 85-93% of reads; it is never an error. `docs/output_reference.md` states that the BLAST/UMI-dedup FASTA, Q137 gates 1/2/4 and the ambiguity flag use the best-viral representative read of each molecule (diagnostic only). The reference-strategy benchmark note now carries a do-not-cite banner on its ~18x/~3.7x figures.
- **`reactivating` is no longer called for viruses whose latency is not observable in RNA (AUD-19, D-01).** Only EBV can reach it; VZV, HHV-6A/6B, HHV-7 and HCMV give `indeterminate`. A programme call is `indeterminate` when any productive catalogue marker is unresolved against the index (D-02); `n_cells_reactivating` reads `not_applicable` where unreachable (D-03); a catalogue whose rows for one virus disagree on `latency_observable_in_rna`/`panel_completeness` is rejected (D-18). Viruses that are only candidates (`detected` False) get no programme row.
- **Programme overlap groups use shared 31-mers, not transitive exon overlap (AUD-01).** A productive marker (e.g. EBV BMRF1) no longer shares a group with latent markers through a chain of overlapping genes; groups are greedy maximal cliques at a 0.5 shared-31-mer threshold. HSV-2 LAT copy `HHV2s01` is merged into `HHV2s02` (catalogue 75 -> 74 rows), like HSV-1. `gene_programs.tsv` `product`/`note` text now follows the RefSeq product annotation (AUD-03). Programme numbers (breadth, independent groups) change for EBV, HCMV, HHV-6A/6B.
- **Detection summary leads with called cells (B-05).** `summary.txt` per-virus lines now start with the called-cell count and percentage; the all-barcode figure is labelled as such. Anything parsing "cells infected (x%)" from `summary.txt` must be updated.
- **`read_start_gini` is computed over all reference positions (B-21)**, so a single-position pile-up no longer scores 0; a blank `viral_minus_host_bitscore` is "not assessed" rather than host-preferred (B-34). `detect_genes`, `super_expressor`, `measure_positive_control` and `below_gate_counts` reject bad input clearly (B-06..B-09); a `-w` list with Cell Ranger `-1` suffixes matches (B-27).
- **`build-ref` rebuilds into a non-empty output directory correctly (C-05)**: `kb ref` always gets `--overwrite`, so the index and t2g are rebuilt instead of silently kept.
- **Host cDNA cache is keyed by Ensembl release; new `build-ref --ensembl-release N` (C-06).** A release that cannot be determined is an error, not a silent fall-back to 116. Existing flat caches are re-downloaded once.
- **NCBI fetch integrity (C-14/C-15):** the FASTA is checked against the GenBank length and version before caching; duplicate accessions are collapsed (versioned form wins); the NCBI email is only required when something must be downloaded; cache write errors are reported cleanly. `data fetch` reports an unwritable cache directory cleanly (C-16); one short t2g row no longer flips the file to the legacy format (A-17).
- **`validate-run` works from any directory (AUD-10)** and resolved input paths are recorded outside the run fingerprint.
- **Competitive alignment emits secondaries, so exact host/virus ties are `tie` (AUD3-04, finding P3-02):** `align_reads_to_viral` ran `minimap2 -ax sr`, whose preset is `--secondary=no`; a read scoring the same on host and virus was reported as a win for whichever side minimap2 picked, and `molecule_verdicts` never held a `tie`. It now runs `-ax sr --secondary=yes -p 0.8 -N 5` (`evidence.MINIMAP2_ALIGN_OPTIONS`), matching the documented policy. The evidence manifest records `minimap2_options` and, for host-verdict runs too, the minimap2 version. Evidence made before this change has no ties; re-run it before using `molecule_verdicts.tsv.gz`. No verdict rule, tier cutoff or threshold changed.
- **Clean errors instead of tracebacks (C-07/08/09/10):** non-interactive `--reference` prompt honours `--yes`; snakemake failures report the exit status and log directory; missing `-w`, `--host-h5ad`, `--cell-types`, `--called-cells-file` and `--positive-control-receipt` are caught before the run; numeric options are range-checked.
- **Run fingerprint covers more inputs (C-13):** the STAR host index directory, `--host-h5ad`, `--cell-types`, `--host-fasta`, `--viral-fasta` and the fetched `-acc` sequences. `--resume` of such runs started on earlier versions is refused; use `--overwrite`. The resume error now lists the differing options.
- **`--gene-programs` no longer fails at DAG build (AUD-08)**: `rule gene_programs` listed `results/viral_summary.tsv`, which no rule produces.
- **Programme per-cell table keeps cells with selected-layer-only evidence (AUD-02)**, so `n_cells_*_selected_layer` is no longer understated; the never-filled `n_cells` column of `gene_program_cells.tsv` is removed (A-14).
- **GTF normaliser hardening (A-06/07/08):** LMP-2A/2B native-exon wrap genes are tagged `wrap_native`; invalid coordinates raise; skipped two-strand genes and missing flatfiles are logged.
- **Host-response now pools viral genes per virus and reports exact statistics (AUD-04, B-12/13/14/17).** Per-virus outputs are named by virus display name; the differential t-test uses df = n-2-k (matches OLS), zero-variance genes get direction `none` and NaN p, BH runs over finite p only, cpm labelling keeps exactly the stated number of positives, and `observed_group_auc` is now `observed_cv_auc`.
- **`--umap` no longer crashes without `plots/` (AUD-05)**, and fails non-zero when QC removes every cell (B-24).
- **`viralscan evidence` picked the wrong reads for small targets (DSR-17).** The
  exact-read replay re-runs `kallisto bus -n` and then selected reads with the
  *primary* run's `matrix.ec`. kallisto numbers equivalence classes in discovery
  order, which differs between multithreaded runs, so the replay's records named
  different classes: rare targets returned mostly host reads (and a different set
  on every rerun). The replay now captures and parses with its own
  `lineage_bus/matrix.ec` and refuses a replay whose `transcripts.txt` differs from
  the primary run's. Evidence produced before this fix for low-abundance targets
  must be regenerated.
- **Capture model is now exact, and a positive control certifies only its own
  scope (SENS-CORR-01/02).** `fragment_capture` treated overlapping k-mer
  windows as independent and overstated capture (0.90 at 10 % divergence, 90 bp,
  k=31); it now delegates to the exact substitution-only DP
  (`fragment_capture_exact`: 0.25 at 10 %, 0.063 at 15 %, 0.013 at 20 %), and the
  implied-divergence inversion is renamed `substitution_model_implied_divergence`
  and labelled a heuristic. A positive control now needs
  `--positive-control-scope exact_sequence --positive-control-virus-key KEY` to
  certify a row; an existing control with no scope is read as `panel_mechanics`
  and certifies no virus (warned once). `positive_control.json` gains `scope`,
  `target`, `certified_targets`. Run fingerprints of runs that do not set the
  new flags are unchanged.
- **Evidence replay now mirrors the primary `kb count` chain.** `viralscan
  evidence` re-pseudoaligns with the primary run's barcode correction
  (`bustools sort → correct → sort` with the run's on-list, gzipped on-lists
  decompressed first) and with the primary run's `--strand` setting, instead of
  capturing from the uncorrected, default-strand BUS. Without both, the
  evidence layer could trace reads kb had discarded or mis-assign strand.
- **Fast multimap rerun swaps the host-viral evidence layer too.**
  `counts_host_viral_selected` is method-dependent; `rerun-multimap` now swaps
  it from the matching per-method layer (`counts_host_viral_selected_equal` /
  `_host_conservative` / `_unique_weighted`, stored on every run) and falls
  back to a full rerun when they are absent, instead of leaving stale
  evidence-tier labels under the new method.
- **`validate-run` accepts `--no-multimapping` runs.** It discovers kb's
  `adata.h5ad` (existence/readability only — not a v3-contract file) instead of
  erroring on the absent `adata_multimap.h5ad`.
- **Evidence manifest records the run fingerprint again.** `run_manifest.json`
  is read from the run root (one level above the per-sample `--run-dir`).
- **HPV high-risk panel carries type 68, not 69.** The IARC Group 1 / clinical
  14-type list is 16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68;
  `hpv_genes.tsv` is regenerated from NCBI with DQ080079.1 (HPV68a) replacing
  AB027020.1 (HPV69), and a test now pins the set against the external list.
- **Merged multi-accession GTFs no longer glue lines.** The whole-genome
  fallback GTF now ends with a newline (GTF format version bumped to 3, so
  cached GTFs are rebuilt), and GTF attributes are spec-formatted with
  per-attribute semicolons.
- **`fragment_capture` returns 0 for fragments shorter than k** (no k-mer
  exists to pseudoalign) instead of an all-or-nothing clean-fragment
  probability, and `poisson_detection_probability` sums its tail in log space
  so expectations above ~745 no longer underflow to a wrong 1.0.
- **NCBI fetching hardening:** proactive ≤3 req/s (≤10 with API key) request
  spacing, a warning when an unparseable GenBank location interval is dropped,
  `--ncbi-api-key` on the run command (env var was already honored), and panel
  re-fetch now removes stale FASTA/TSV files, not just stale GTFs.
- **`--strand` help text** documents the kallisto 0.52.0 behaviour (the legacy
  per-technology forward default no longer engages; no flag effectively means
  unstranded) instead of referring to a kb per-technology default.
- **Packaging/docs:** the `full` extra requires `snakemake>=9,<10` (matching
  the conda pin), and the README notes the PyPI release can lag this repo.
- **Tests:** the evidence-dispatch and same-reads tests now assert behaviour
  (previously one swallowed all exceptions with no assertion and one grepped
  the source), and the Snakefile host-filter input test execs the real
  `_kb_count_inputs` instead of matching substrings.

### Changed
- ViralScan 3.0 uses one corrected CB-UMI molecule as its count unit. BUS read
  multiplicity is audit-only and pre-v3 corrected/combined counts are incompatible.
- The complete selected-method matrix is `X = counts_unique +
  counts_ambiguous_allocated`; fractional outputs are labelled molecule estimates.
- Multimapping methods are `equal`, `host-conservative`, `unique-weighted`,
  `em-global`, and hierarchical `em-cell`; the ambiguous plain `em` name is removed.
- Expanded anellovirus references are explicit opt-in and fail closed on missing
  accessions unless `--allow-partial-panel` is supplied.

### Added
- Exact run fingerprints, atomic manifests, explicit `--resume`/`--overwrite`,
  `viralscan doctor`, and `viralscan validate-run`.
- Streaming corrected BUS molecule resolution, conservation audits, v3 schemas,
  EM convergence diagnostics, and order/buffer-size invariance tests.

## [2.7.0] - 2026-07-21

### Added
- **Primary-call matrix consistency.** `--multimap-primary-call {legacy,unique-only,confidence}`
  now drives the viral *numerator* everywhere — `viral_summary.tsv`, `per_cell_viral.tsv`, plots,
  and cell-type enrichment — while full-library *denominators* (`umi_per_10k`, `viral_fraction`)
  stay on the complete expression matrix. Documented in `docs/output_reference.md` and the report.
- **`viralscan evidence --read-start-profile`** — per-position 5′ read-start distribution along each
  viral reference (`read_start_profile.tsv`), with `--dedup {umi,markdup,none}` and `--bin-size`.
  Exposes 3′ bias, subgenomic-RNA junctions, and EVE/integration hotspots the coverage summary
  cannot.
- **`viralscan evidence --cell-tags`** — writes `viral_reads.tagged.bam` with `CB`/`UB` tags from
  read names so viral reads can be grouped by cell in IGV.
- **`results/reference_provenance.json`** — every run records the viral reference used (index / t2g /
  GTF, technology, multimap settings) and the viral accessions in-reference and detected, so results
  are traceable to their annotation.
- **Host-response gene-symbol annotation + enrichment** — hostresponse outputs map Ensembl IDs to
  HGNC symbols (mygene.info, best-effort) and translate to symbols before `gget.enrichr`, so results
  are interpretable when the host matrix uses Ensembl IDs.
- **EVE (endogenous viral element) annotation toolchain** — `covid_viralscan/scripts/annotate_eve.py`
  + a hardened SLURM pipeline annotating artifact reads against GRCh38 (Phase A/B/C).
- **Non-human host references** — `build-ref --host` supports mouse, rat, macaque, and 13 other
  Ensembl species.
- **Documentation** — an 8-vignette suite under `docs/vignettes/` (each grounded in a manuscript
  result; 6 run in CI), plus `docs/ROADMAP.md` and `docs/IMPLEMENTATION_PLAN.md`.
- **CI** — an nbmake job executes the CI-runnable vignettes so docs cannot silently rot.

### Changed
- The HTML report's "Infected cells (any virus)" count now counts unique virus-positive barcodes
  rather than cell×virus rows.

### Performance
- **Multimapping correction is ~3× faster with ~6× lower peak memory** on deep samples: duplicate
  `(cell, ec)` records are collapsed (their per-gene shares are linear in count) before the hot loop,
  and a per-record `pd.isna` check was hoisted out. Output is numerically identical (all layers,
  rtol=1e-9, golden-tested).

### Fixed
- **EVE annotation correctness** — Phase B no longer treats BLAST subject-local coordinates as
  chromosome coordinates (clone/scaffold subjects are flagged); chromosome names are normalized
  (`chr7`↔`7`, human `NC_000001`–`24`) with a namespace-mismatch warning; legacy `gi|…|ref|NC_…|`
  subject IDs are handled; SLURM read-extraction surfaces real errors instead of swallowing them.
- **Duplicate reference gene IDs** now raise a clear, actionable error instead of a cryptic crash or
  silent undercount.
- `host_viral_ambig_fraction` is clamped to `[0, 1]`.

### Internal
- Repository slimmed for publication: the 195 bundled viral GTFs are untracked (fetched from Zenodo
  via `viralscan data fetch`); the git history was garbage-collected (`.git` 29 GB → 8.9 MB);
  institutional and third-party paths were scrubbed from shipped files. Dead code and a dead
  sparse-matrix special case were removed; primary-call plumbing was deduplicated.

## [2.6.0] - 2026-07-17

### Added
- **EVE (endogenous viral element) artifact flags in `viral_summary.tsv`** — new
  `eve_risk`, `accession_breadth`, and `host_viral_ambig_fraction` columns help
  distinguish genuine viral signal from host-homology artifacts (e.g. anellovirus reads
  that share sequence similarity with GRCh38 non-coding regions). (closes PLAN SH2.6)
- **Sibling cross-mapping warning** — a `sibling_crossmap_note` column flags closely
  related virus pairs (HHV-6A/6B, HSV-1/2) where one member dominating may reflect
  cross-mapping rather than independent detection. (SH2.3)
- **`evalue_flag` column in `<virus>_depth_diagnostics.csv`** — each stable gene is
  now classified as `fragile` (E < 1.5), `moderate` (1.5 ≤ E < 3), or `robust` (E ≥ 3)
  alongside the existing depth-adjusted odds ratio and E-value.
- **`model_auc_depth_adjusted_mean/sd` in `hostresponse_metrics.csv`** — AUC of the
  stable-gene panel logistic with `log1p(host_depth)` added as a covariate (same balanced
  split as the headline model and the depth-alone baseline). Shows how much of the
  stable-gene predictive power survives explicit depth adjustment.
- **`--hostresponse-label`, `--hostresponse-depth-match`, `--hostresponse-control-mito`,
  `--hostresponse-differential`** now exposed on the main `viralscan` pipeline CLI and
  persisted in `RunConfig` / `DEFAULTS` / the Snakemake entry block. Previously these
  knobs were only accessible via `viralscan hostresponse` standalone.

### Changed
- **Default `--multimap-method` is now `host-conservative`** (was `equal`). For viral
  detection on combined host+virus references, specificity is prioritised: host-virus
  ambiguous equivalence-class mass is kept out of primary viral counts by default. Pass
  `--multimap-method equal` for a fast unbiased first pass, or `em` for iterated allocation.

### Performance
- **Faster multimapping correction** — the EM step was vectorised into a sparse
  matrix-vector formulation, per-EC invariants were hoisted out of the inner loop, and
  weight computation now uses direct CSR buffer access (~3x faster on the hot path;
  output is byte-identical and regression-tested).

### Fixed
- **Packaging: complete dependency declarations.** `anndata` is now an explicit
  dependency (previously satisfied only transitively via scanpy, despite an eager import
  in `scripts/multimap.py`); `environment.yml` gained `anndata` and `scikit-learn`; the
  Docker image installs with `--no-deps` against the complete conda environment.
- **Benchmark portability.** The reference-strategy benchmark resolves STAR from `PATH`
  (or `VIRALSCAN_STAR_BIN`) instead of a hardcoded author-specific absolute path.

## [2.5.0] - 2026-07-03

### Added — host-response scientific-hardening (v2.5)

- **Depth-confound diagnostics (always on).** `hostresponse` now reports the AUC
  from sequencing depth *alone* (`depth_alone_auc` in `hostresponse_metrics.csv`)
  next to the model AUC, plus per-gene depth-adjusted odds ratios and Ding &
  VanderWeele **E-values** (`<virus>_depth_diagnostics.csv`). Surfaces the fact
  that the raw `counts >= threshold` label tracks library size.
- **Depth-independent labels / design.** `--label {raw,cpm,fraction}`
  (depth-normalized, prevalence-matched positive call) and `--depth-match`
  (coarsened-exact depth-matched cohort) give a depth-independent host-response
  estimate.
- **`%mito` control** (`--mito-control`, default on) — per-cell mitochondrial
  fraction added as a covariate to the per-gene E-values, with leave-one-out for
  mitochondrial genes.
- **`--gene-symbols`** — annotate host-response CSVs with HGNC symbols (mygene.info).
- **`--differential`** — genome-wide, depth/%mito-adjusted differential-expression
  table (`<virus>_differential.csv`; partial correlation, p-value, BH-FDR, direction).
- **`viralscan check-whitelist`** — barcode/whitelist chemistry-mismatch preflight;
  the main run also warns automatically when an explicit `--whitelist` is supplied.
- **Called-cell denominators** — `viral_summary.tsv` reports `pct_infected_called`
  (over knee/emptyDrops/external called cells) alongside the all-barcode rate.

### Changed

- Corrected the package docstring: ViralScan supports single-cell RNA-seq only
  (bulk was never supported).

## [2.4.0] - 2026-07-02

### Fixed
- Release-readiness pass: `evidence_run.py` migrated to the typed `RunConfig`
  (last script on the legacy raw-dict path); `RunConfig.from_yaml` now normalizes
  the `output` trailing separator; early `kb` preflight in `build-ref`; canonical
  GitHub URLs (`mdmanurung/ViralScan`); single-source package version
  (`viralscan.__version__`, read dynamically by `pyproject.toml`).

### Added
- **`viralscan hostresponse` subcommand** — run host-response analysis on a
  completed viralscan output directory without re-running Snakemake. Trains
  per-virus L2 logistic regression models (Luebbert et al. 2025) and runs
  randomized Lasso stability selection to identify stably virus-associated host
  genes. Optional `--enrichment` flag runs pathway enrichment via
  `gget.enrichr` (requires `pip install "viralscan[enrichment]"`). The same
  analysis also runs inline during a full `viralscan` run when `--host-h5ad`
  is supplied.
- **CLI reference docs** for `viralscan evidence`, `viralscan rerun-multimap`,
  and `viralscan hostresponse` — all three subcommands were shipped but
  previously undocumented.
- **`README.md` Limitations section** — documents five known limitations: FPR/FNR
  not yet characterized, global-pool EM (vs per-cell), untested chemistries beyond
  10xv2/v3/Drop-seq, cross-homology inflation for HHV-6/KDM2A, and ambient RNA
  not corrected.
- **EM global-pool caveat in `em_gene_abundances` docstring** — explains that EC
  counts are pooled across all cells before EM; returned theta is transcriptome-wide,
  not per-cell; differs from alevin-fry/STARsolo; faster but ignores cell-to-cell
  abundance variation.
- **`tests/test_evidence_subcommand.py`** (18 tests) — parser tests for all
  `viralscan evidence` flags and dispatch tests confirming `main()` routes
  `_subcommand="evidence"` to `run_evidence()`.
- **Planted-signal test for `run_hostresponse`** (`TestHostresponsePlantedSignal`)
  — synthetic 200-cell matrix with 40 virus-positive cells and 5 genes amplified
  10× in positive cells; asserts ≥3 planted genes appear in top-10 by stability
  probability.

### Fixed
- sklearn ≥ 1.8 `FutureWarning` in host-response stability selection: switched
  to `solver='saga', l1_ratio=1.0` on sklearn ≥ 1.8 (was `penalty='l1'`, now
  deprecated); older sklearn still uses `penalty='l1', solver='liblinear'`.

### Changed
- Multimapping is now a Snakemake checkpoint with a new `viralscan rerun-multimap`
  subcommand to switch allocation methods on an existing run without re-running
  `kb count`. The default `--multimap-method` is `equal` (fastest; enables instant
  in-place layer swaps). `host-conservative` remains available and is **recommended
  when host-virus cross-homology matters** (e.g. HHV-6 / *KDM2A* / *DR1*), where it
  keeps host-ambiguous equivalence-class mass out of primary viral counts.

---

## [2.3.0] - 2026-05-13

### Added
- Sphinx + MyST documentation skeleton under `docs/`; Read the Docs config at `.readthedocs.yaml`.
- `CITATION.cff` for software citation.
- `environment.yml` for reproducible conda environments.
- `Dockerfile` (mamba-based) and `Singularity.def` for containerised HPC runs.
- `tests/test_errorhandler.py` — direct unit tests for every `errorhandler()` branch.
- `tests/conftest.py` — session-level `pyfiglet` stub so unit tests pass without the optional dep.
- Type annotations on all public functions in `menu.py`, `utils.py`, and `constants.py`.
- Codecov integration: coverage XML uploaded in CI; badge added to README.
- `viralscan data fetch` for downloading the external viral annotation panel from Zenodo.
- Ambiguity-aware multimapper evidence outputs and diagnostic AnnData layers.
- Optional host pre-subtraction support before viral quantification.
- Cell-type enrichment tables and HTML report section.

### Fixed
- Comma-separated custom FASTA/GTF reference inputs are now materialized before `kb ref`.
- Comma-separated custom GTF files are all parsed during viral accession discovery.
- CI dependency installation now matches the modules imported by tests and type checks.
- Release/container install metadata now targets the local `2.3.0` package.

---

## [2.2.0] - 2026-05-08

### Added
- **PR 14 Bug fixes & performance**
  - C1: Fixed double-counting of unique reads in multimapping mode (`multimap.py`).
  - C2: Fixed detection threshold off-by-one (`>` → `>=` in `detection.py`).
  - C3: Replaced O(n) `list.index()` scan in EC parsing with O(1) `dict` lookup (`multimap.py`).
  - C4: Replaced `iterrows()` with `itertuples()` for ~3× faster BUS iteration (`multimap.py`).
  - C5: Changed `2> /dev/null` → `2>&1` in Snakefile so "no reads pseudoaligned" is capturable.
  - C6: Changed `mkdir` → `mkdir -p` in Snakefile to survive re-runs.
  - C7: Removed redundant `f.close()` after `with` block in `analysis.py`.
  - C8: Fixed file-handle leak in `detection.py` by wrapping `found_genes_file` in `with`.
  - C9: Consolidated conflicting `var_names` resolution blocks in `umap.py`; added O(1) index.
  - C10: Added `sc.pp.highly_variable_genes()` before PCA in both UMAP branches.
- **PR 11 Interpretation & Reporting**
  - Structured TSV outputs: `viral_summary.tsv` and `per_cell_viral.tsv`.
  - Self-contained HTML report via Jinja2 template with embedded base64 PNG plots.
  - Normalised metrics: viral prevalence (% cells) and viral load per 10k UMI.
  - Configurable thresholds via `--se-threshold`, `--detection-threshold`, `--min-counts`, `--min-genes`.
- **PR 12 Combined Host+Virus Reference Builder**
  - `viralscan build-ref` subcommand.
  - `src/viralscan/scripts/build_reference.py`: Ensembl cDNA download + NCBI viral fetch + `kb ref`.
  - Supported host species table (`ENSEMBL_SPECIES`) in `constants.py`.
- **PR 10 NCBI accession → reference**
  - `src/viralscan/scripts/ncbi_fetch.py`: fetch FASTA + GTF from NCBI by accession.
  - `--ncbi-accession` / `-acc` and `--ncbi-email` CLI flags.
  - Per-accession cache under `~/.cache/viralscan/ncbi/` with SHA256 validation.
  - Accession validation regex; exponential backoff on 429/5xx.
- **PR 5 Tests**
  - `tests/test_ncbi_fetch.py` (18 unit tests).
  - `tests/test_cli.py` (8 test classes, 144 tests total passing).
  - `tests/test_createconfig.py`, `tests/test_analysis.py`, `tests/test_build_reference.py`.
- **PR 4 Tooling**
  - `.pre-commit-config.yaml` (ruff, ruff-format, end-of-file-fixer, etc.).
  - GitHub Actions CI (matrix Python 3.9–3.12 × ubuntu/macos) and release workflow.
- **PR 3 Code cleanup**
  - Removed ~620 lines of dead/commented-out code from `detection.py` and `umap.py`.
  - Moved `VIRUS_NAME_MAP` to `constants.py`; shared `load_config()` in `utils.py`.
  - Centralised logging via `configure_logging()` / `setup_script_logging()`.
- **PR 2 Correctness & security**
  - Fixed `--reference`, `--visual`, `--multimapping` boolean flag bugs (§1.1).
  - Replaced `os.system` and `shell=True` subprocess calls with safe `subprocess.run([...])`.
  - Replaced all bare `exit()` with `sys.exit()` via `_die()` helper.
  - Added `_check_required_tools()` preflight for `kb` and `snakemake`.
- **PR 1 Hygiene**
  - README fixes, expanded `.gitignore`, `pyproject.toml` classifiers and optional deps.

### Changed
- `createconfig.py`: boolean config values are now native Python `bool`, not strings.
- `umap.py`: removed `== "True"` string comparisons for config booleans.
- `analysis.py`: replaced `!= "None"` guard with `config.get('gtf')`.
- Snakefile: whitelist comparison no longer uses literal `"None"`.

---

## [2.1.0] - 2025-01-01

*(No structured changelog was kept before v2.2.0. See git log for history.)*

[Unreleased]: https://github.com/mdmanurung/ViralScan/compare/v2.7.0...HEAD
[2.7.0]: https://github.com/mdmanurung/ViralScan/compare/v2.6.0...v2.7.0
[2.6.0]: https://github.com/mdmanurung/ViralScan/compare/v2.5.0...v2.6.0
[2.5.0]: https://github.com/mdmanurung/ViralScan/compare/v2.4.0...v2.5.0
[2.4.0]: https://github.com/mdmanurung/ViralScan/compare/v2.3.0...v2.4.0
[2.3.0]: https://github.com/mdmanurung/ViralScan/compare/v2.2.0...v2.3.0
[2.2.0]: https://github.com/mdmanurung/ViralScan/compare/v2.1.0...v2.2.0
[2.1.0]: https://github.com/mdmanurung/ViralScan/releases/tag/v2.1.0

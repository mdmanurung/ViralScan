# ViralScan 3.0 execution tracker

## Next action

Live tracker = open rows (`[ ]`, `[~]`, `[!]`) with their section headers, the shared contract sections
(How to use, Critical path, Gate dashboard, External inputs, Stop rules) and `CAT-31`. Every done row and the
dated "Next action" log up to 42b59f8 are frozen in the reproducibility repo:
[`provenance/PLAN_archive_2026-10.md`](https://github.com/mdmanurung/ViralScan-reproducibility/blob/main/provenance/PLAN_archive_2026-10.md) (REPO-SPLIT-01 batch 4e, D63/D67/D68).

- [x] `BUILDREF-HOSTDUP-01` (2026-11-03, `codex/buildref-hostdup`, found by repro run `cat08_hpv11`): `build-ref --host human` on Ensembl 116 cDNA failed with `Exact duplicate sequences ... ENST00000604066.1 and ENST00000630556.1`. `_check_reference_structure` now takes `viral_identifiers`; with it, exact duplicates among host records pass and are logged as a count (kallisto 0.51.1 indexes identical targets as one equivalence class, checked on a 3-target FASTA). Viral-viral and viral-host exact duplicates still fail. Without the set (viral-only FASTAs, low-complexity gate) behaviour is unchanged. Test: `test_exact_duplicates_among_host_records_pass_but_viral_ones_fail`.
- [x] `TCHEM-DROPSEQ-01` (2026-11-02, `codex/tchem-dropseq`, user decision, found by the HSV-1 SRR8315713 run): `tchem.compare` returned `no_call` for Drop-seq (tchem layout `Drop-seq` has no kb command), so the default `--tchem auto` stopped a run with `-x dropseq`. It now returns `agree` when the layout name and `-x` match after lower-casing and dropping punctuation. Every other layout without a kb technology stays `no_call`. Test: `tests/test_tchem.py::TestCompare::test_layout_without_kb_command_agrees_only_by_name`.
- [x] `PLAN-SPLIT-EXP-01` (2026-11-02, `codex/plan-split-experiments`, user decision: experiment work lives only in the reproducibility repo): 64 open experiment rows
  (64 row blocks, 62 IDs) moved to the repro `PLAN.md`; the list is in "Moved to the reproducibility tracker". `SCAUTOQC-01` split: code stays, its real-data gate is repro `SCAUTOQC-01-RUN`.
  `CLAUDE.md` "Sibling repository" now names these rows. No code or test change.
- [x] `PLAN-RECON-01` (2026-11-02, `codex/plan-recon`, user: finish the package before analysis): open rows checked against the code. Closed: `SENS-CORR-02`, `ANDET-04` (opt-in). Stale notes corrected: `CAT-41`, `DOC-01`, `DEF-08`, `MECH-E`. Package code still open without a user decision: `AUD3-02-OPEN` (closed in `codex/aud3-02-calibration`). Open on a user decision: `REF-13` CDS normaliser, `DEF-04`, `SCAUTOQC-01` defaults, `DLIST-02`, `CAT-30`, `DOC-10`, `PATHS-01`, and `CAT-04` (the `iciHHV6_possible` flag is name-based by design and pinned in `test_panel_integrity_contract.py`; a real test needs a ciHHV-6 carrier library or a new data-driven rule). Every other open row needs real data, the cluster, the network or maintainer infrastructure. Decision batch (user, 2026-11-02): `DEF-04`, `CAT-04`, `PATHS-01` (#164), `REF-13` (#165), `DLIST-02` (#167) and `CAT-30` closed; `DOC-10` docs pass done (#166), gate goes in the repro at the next re-pin; `SCAUTOQC-01` waits on the user's SRR8315713 run.
- [x] `HOSTLOSS-IT-01` (2026-11-02, `codex/hostloss-starsolo-it`, HANDOFF open item from `REVIEW-FIX-E`; test only): `tests/integration/test_host_filter_star.py::test_host_loss_on_a_real_starsolo_run_and_its_primary_tree_guard` runs `viralscan run --host-filter starsolo` on the tiny fixture (STAR genome from `host.fasta`, `--tchem off`), then `host-loss` on that sample: the default output lands in `posthoc/host_loss/` (D15 0/1), and `--output` under `<sample>/results/` or the run root is refused, `--force` included, with nothing created. Disabling `refuse_primary_tree` makes it fail. The STAR genome build moved to a `_star_genome` helper shared with the SW-07 test.
- [x] `SIGTERM-GROUP-01` (2026-11-02, `codex/sigterm-group`, HANDOFF open item from `REVIEW-FIX-D`): the process-group forwarding in `menu._run_snakemake` moved to `run_safety.run_in_group` (own session, SIGINT/SIGTERM forwarded to the group, SIGTERM exits 143, `check`/`capture_output` like `subprocess.run`). `kb ref` (`build_reference._run_kb_ref`) and the `--strand auto` pilot (`strand._kb_pilot`) now use it, so a SIGTERM during the index build or pilot no longer orphans kb/kallisto. STAR and other native tools still use `subprocess.run`. Tests: `tests/test_run_lifecycle.py` (existing SIGTERM test now runs through the helper, plus a capture/raise test); kb-ref test patches moved from `subprocess.run` to `run_in_group`.
- [x] `CORES-AFFINITY-01` (2026-11-02, `codex/cores-affinity`, HANDOFF open item): the `--cores` clamp in `errorhandler` uses `len(os.sched_getaffinity(0))` (the CPUs this process may use, e.g. a Slurm allocation) and falls back to `os.cpu_count()` where the call is missing (macOS). Before, a job given 8 of 24 CPUs kept `--cores 24`. Test: `tests/test_config_zero_and_cores.py::test_cores_clamped_to_affinity_not_node`.
- [x] `REVIEW-FIX-E` (2026-11-02, `codex/posthoc-safety`, review follow-up; no analysis-output change): `rerun-multimap` validates the merged config (`RunConfig.from_snakemake_config`) before copying, exits 2 on a bad override, and `--multimap-em-tol` / `--multimap-em-max-iter` are checked at parse time; right after the copy it restamps the READMEs NOT FINISHED and deletes the copied `artifact_inventory.json` / `combined_summary.tsv`; on a v1/v2 source it drops `log/analysis.txt`, `log/found_genes.txt` and `plots/SuperExpressor_*` from the copy. `run_safety.refuse_primary_tree` (called by `require_fresh_output`, so every post-hoc verb, `--force` included) refuses an `--output` inside a run (a `run_manifest.json`) or sample tree (`config.yaml` plus `kb-python/`, so a bare project `config.yaml` does not count) outside `posthoc/`; `evidence` calls it too, except under `results/evidence/` (auto-evidence's home). `hostresponse --force` copies the old tree to `.<out>.prev` and restores it if the run raises; `too_few_shared_barcodes` now exits 1 (status row kept; the old outputs are not restored in that case). `auto_evidence` removes `results/evidence/<genus>/` folders (those with `evidence_manifest.json`) of genera no longer detected. Tests: `tests/test_posthoc_safety.py` (8, all fail on main); `tests/test_evidence_subcommand.py` outputs moved under `posthoc/`.
- [x] `REVIEW-FIX-D` (2026-11-02, `codex/run-lifecycle`, adversarial-review fixes, no analysis output changed): `run` takes `<out>/.viralscan.lock` (O_EXCL, pid+host; stale same-host pid replaced with a warning; atexit release) before the output dir is touched, so a second `--resume`/`--overwrite` is refused; READMEs are stamped NOT FINISHED right after `prepare_output_directory` (before chemistry blocks, NCBI fetch, index build, strand pilot); `--resume` runs `snakemake --unlock` first; `_run_snakemake` starts snakemake in its own session and forwards SIGINT/SIGTERM to its group (SIGTERM exits 143); `.<name>.tmp` crash litter and the lock do not count as content or inventory; a file as `--output` or an unwritable parent is a `RunSafetyError`; `summary.txt` Runtime is cumulative over finished attempts; `--resume` keeps the first `chemistry_detection`/`tool_binaries` and appends `resume_attempts`; README log hint fixed; unreachable `resume and overwrite` check deleted. Tests: `tests/test_run_lifecycle.py`; `menu.subprocess.run` patches in `test_posthoc_layout.py`/`test_rerun_multimap.py` became `_run_snakemake` patches. Not done: a stale lock from another host needs a manual delete; a symlink-loop `--output` still raises from `Path.resolve`.
- [x] `REVIEW-FIX-G` (2026-11-02, `codex/cli-errors`): `menu.main()` catches `OSError`/`EOFError` (truncated gz, directory given as a file, unwritable path) and prints one `ERROR  <msg>` line, exit 1; the traceback shows only with `--verbose`. `translated-index` checks `--proteins-dir` and `--host-mask` before the PalmDB download. Tests in `tests/test_cli_mature.py` (3 of 5 CLI cases and the translated check failed on main; `-o <file>` and `-t <dir>` already exit cleanly there, covered by another fix).
- [x] `REVIEW-FIX-F` (2026-11-02, `codex/validate-inventory`): `validate-run` reads every `results/artifact_inventory.json` it finds (outside `posthoc/`) and checks each listed file against the sample directory: missing is `inventory_artifact_missing`, a different sha256 is `inventory_hash_mismatch`. Before, a deleted `count_audit.tsv` or a truncated `per_cell_viral.tsv` still gave ok=true. Runs with no inventory (layout v1/v2) are not flagged. Test: `tests/test_validation.py::test_validate_run_checks_inventory_artifacts`.
- [x] `REVIEW-FIX-A` (2026-11-02, `codex/allocation-science`, user approved the output changes): (1) `em-cell` falls back to the equal split per EC: an EC whose genes have no cell-local unique support stays out of the cell fit and is split equally even when the cell has support elsewhere (before, any support elsewhere let the sample prior decide it); a cell with no support anywhere is unchanged. (2) `host_viral_ambig_fraction` = mixed / (unique viral + mixed) from `counts_unique` and `counts_host_viral_ambiguous`, no clamp, NA only at a zero denominator. (3) `_headline_totals` reads `counts_unique` through `get_layer`. (4) `rerun-multimap` swap rewrites `uns["multimap_diagnostics"]` `method`, `model_scope`, `max_iter`, `tolerance` and drops EM-only blocks. (5) `load_barcodes` raises on a `-1` collision; new `count_audit.tsv` column `offlist_barcode_molecules`. Tests: `tests/test_allocation_science.py`.
- [x] `REVIEW-FIX-C` (2026-11-02, `codex/verdict-consumers`): `read_support_summary.tsv` `n_human_equal_or_better` now counts molecules carrying the `human_equal_or_better` flag (was `not gate_competitive`, which also counted host-only and unaligned molecules); `cell_type_enrichment` raises `ValueError` when two distinct obs barcodes collide after the `-1` strip; `support_fraction` builds `set(obs_names)` once. S unchanged by user decision (frozen protocol): host-only and unaligned molecules stay outside S, now stated in the `tier_measures` docstring and `docs/results.md`, pinned by a test (S=0.5, 99 host-only, 1 unaligned). Tests in `test_read_support.py`, `test_detection.py`, `test_tier_measures.py`.
- [x] `REVIEW-FIX-B` (2026-11-02, `codex/config-zero-and-cores`): an explicit 0 for `emptydrops_lower` / `hostresponse_stab_min_prob` is kept (`x or DEFAULT` became an absent-test in `RunConfig` and in `hostresponse`'s menu defaults). `check_numeric_ranges` (CLI `errorhandler` and RunConfig) now rejects `emptydrops_fdr` outside (0,1), `emptydrops_niters`, `hostresponse_n_stab_iter`, `hostresponse_top_n_genes` < 1, and `knee_min_umi` <= 0 (`knee_cells` takes log10 of every barcode at the floor, so 0 would give log10(0)). `emptydrops_lower=0` and `stab_min_prob=0` stay valid. `--cores` above `os.cpu_count()` logs a WARNING and is clamped. Tests: `tests/test_config_zero_and_cores.py`.
- [x] `REVIEW-FIX-H` (2026-11-02, `codex/deps-and-pair-check`): dependency drift. `scanpy>=1.10` (base and `scautoqc` extra; `umap.py` uses `mask_var=`, which needs it) and `use_highly_variable=True` replaced by `mask_var="highly_variable"`; scAutoQC UMAP flag columns fill NaN with `"None"` on any pandas; `scautoqc.py` copies the QC array sctk may mutate; `cellcalling.py` wraps `mmread` in `csc_matrix`/`coo_matrix` (scipy 1.20 returns `coo_array`); `sctk` only for `python_version < '3.12'` (numpy-groupies build fails on 3.12; `--scautoqc` already dies with an install hint when sctk is missing); `uv.lock` regenerated (scanpy 1.9.8 -> 1.11/1.12). Same PR: `fastq_pair_check` rule (`scripts/fastq_pair_check.py`, parallel to `kb_count`, no dependency on it) counts R1 and R2 reads, fails the run with `R1 has N reads, R2 has M` when they differ (kallisto exits 0 and drops reads when R2 is shorter), and writes `log/fastq_pair_check.json`; it is a `rule all` target, not a `multimap` input, so an old run resumes without rerunning multimap; with its raw FASTQs deleted it writes `"pairs": null` instead of failing. One R1/R2 pair per sample, so no multi-file handling. Checks the raw `sample1`/`sample2`.

- [x] `ALIAS-LAYERS-01` (2026-11-01, `codex/drop-alias-layers`, user approved the output change): the four duplicate layers `counts_original` (= `counts_unique`), `counts_corrected` (= `counts_ambiguous_allocated`), `counts_combined` (= `X`) and `counts_host_viral_selected_sibling_weighted` (= `..._equal`) are no longer written by `scripts/multimap.py` or the `rerun-multimap` swap (the swap also deletes stale copies in an old h5ad). One helper, `utils.get_layer`, returns the canonical layer and falls back to the old name; `multimapping.summarize_multimap_evidence` (was the only reader of `counts_corrected`) and `menu._swap_multimap_layer` use it, and `sibling-weighted` now maps to the `_equal` host-viral layer. `X` and every other layer are value-identical. Docs: layer table, `v3_counting_contract.md`, `migration.md`, CHANGELOG. Tests: swap and evidence readers on an old-style h5ad, a new run writes none of the four; `tests/_frozen_oracles/*` untouched. Ten scripts in the reproducibility repo read these names from h5ads that already exist and keep them, so they are unaffected until they read a new run.

- [x] `FINDABILITY-01` (2026-11-01, `codex/output-findability`, on top of `codex/layout-v3`): run `README.txt` and each sample `README.txt` are written at run start (NOT FINISHED) and rewritten at completion, with version, fingerprint, enabled-option files only, how to read a negative and a safe-to-delete list (`outputs.RESULT_FILES`, `safe_to_delete()`); `report.html` file table generated from the same registry; `docs/output_reference.md` has a full multimap-h5ad layer table; `summary.txt` runtime footer idempotent on `--resume`; `run_safety.require_fresh_output` guards `host-loss`, `rerun-programs`, `hostresponse` (new `--force`), `rerun-multimap`, `translated`; `--run-dir` help says run root vs sample dir. Not changed: `discover`/`evidence` (own publication-aware guards). No layer dropped. **Open (user):** alias layers `counts_original` and `counts_combined` are never read and could go; `counts_corrected` is read by `multimapping.py` and the `rerun-multimap` swap (switch them to `counts_ambiguous_allocated` first); `counts_host_viral_selected_sibling_weighted` is a copy of `..._equal`.

- [x] `CLI-EXIT-02` (2026-11-01, `codex/cli-exit-codes`, adversarial-review follow-up to `CLI-MATURE-01`; no analysis-output change):
  malformed `run_manifest.json` (non-object top level, wrong-typed `options`/`input_fingerprints`/`input_paths`, bad
  `field:index` keys) is now a `validate-run` issue, not a traceback; an unwritable `--json PATH` (`doctor`, `validate-run`,
  `check-chemistry`) and a missing `evidence-plot` directory print one error and exit 2. Flag errors that were exit 1 now exit 2
  (`--cores 0` and the other numeric ranges, `--cell-calling external` without its file, the `--star-confirm`/`--discover`/
  `--auto-evidence` combinations, `--r1`/`--r2` count mismatch, `build-ref` without `--host`/`--virus-accessions`, now checked
  before the tool preflight). `run` now applies RunConfig's own validation (`_precheck_run_config`) right after `errorhandler`,
  before the chemistry sniff, the output directory and any NCBI download; it used to raise a traceback after them.
  `--ncbi-accession` with `-i`/`-t`/`-f1` exits 2. New `_warn_ignored_flags`: one WARNING for `--ncbi-email`/`--ncbi-api-key`
  without `--ncbi-accession`, `--host-fasta`/`--viral-fasta` without `--auto-evidence`, `--multimap-*` with
  `--no-multimapping`, `--multimap-em-*` with a non-`em-*` method, `--emptydrops-*` when cell calling is not emptydrops.
  Not flagged: `--umap-n-neighbors`/`--hvg-*` without `--umap` (the QC step uses them regardless). Tests: `tests/test_cli_mature.py`.

- [x] `RUN-CORRECT-01` (2026-11-01, `codex/run-correctness`, adversarial-review fixes, no analysis output changed): cell-type enrichment
  strips a one-sided `-1` barcode suffix instead of dropping cells; `--gene-programs` with `--no-multimapping` is a usage error (exit 2);
  `write_loci` and `write_cell_type_enrichment` delete a stale `viral_loci.tsv` / `cell_type_enrichment.tsv` when a rerun yields none;
  `rerun-multimap` skips `<sample>/posthoc/` and stamps `layout_version` = `LAYOUT_VERSION`; `kb_count` prints the tail of
  `kb_count.log` when kb exits non-zero; the `.snakemake/log` hint names the working directory's `.snakemake` (snakemake runs in the cwd).

- [x] `MOLVERDICT-01` (2026-11-01, `codex/molecule-verdict-unaligned`, user approved the output change): `molecule_verdict._parse_molecules`
  skipped unmapped records before keying the molecule, so `unaligned` was never emitted and unmapped reads were missing from
  `reads`/`counted`. Unmapped records (flag 0x4 or RNAME `*`, not supplementary) now create or join their molecule without a
  score. Output change: `molecule_verdicts.tsv.gz` gains `unaligned` rows and row order follows the first scored-or-unmapped
  record; `tier_measures` `n_excluded_unaligned` becomes non-zero and the support fraction can move (a molecule counted only
  through an unmapped read). The B18 frozen oracle (faf5667) still holds on every case without unmapped reads; the three
  unmapped cases are pinned to explicit rows. Tests `test_a_counted_unmapped_read_is_an_unaligned_molecule`,
  `test_unmapped_cases_keep_unaligned_molecules`. **Open (user):** pinned repro evidence dirs stay self-consistent (per-dir
  manifest hash), but any evidence re-run changes these bytes; re-pin after.

- [x] `LAYOUT-V3-01` (2026-11-01, `codex/layout-v3`, user-approved breaking moves before 3.0.0 stable): `LAYOUT_VERSION` 3
  (supported 1, 2, 3). `log/analysis.txt` -> `results/analysis_genes.txt`, `log/found_genes.txt` -> `results/found_genes.txt`,
  `plots/SuperExpressor_*.{tsv,json}` -> `results/superexpressor/`, in-run `hostresponse/` -> `results/hostresponse/`;
  `discover`/`translated` default to `<sample>/posthoc/<verb>/` (rerun-multimap stays run-level: it copies the whole run).
  Paths are defined once in `outputs.py` (`SAMPLE_LAYOUT`, `V3Outputs`, `analysis_genes_path`/`found_genes_path` with v1/v2
  fallback); `rerun-multimap`/`--resume` copy the old gene list forward. New `<sample>/README.txt`, run `README.txt` and
  `results/artifact_inventory.json` (`artifact_inventory.py`, written by `write_run_complete`). Docs tree synced.
  Decision: `kb-python/`, `read_filtered/`, `host_filtered/`, `star_hostvirus/` stay at the sample root, not
  `intermediate/` (~200 path strings, frozen oracles, kb staging shell block; `kb-python/` holds the matrix users load).
  **Open (user):** the `output-reference-doc` row of `analysis/v3_artifact_inventory.tsv` in the reproducibility repo pins
  `docs/output_reference.md` to a commit; re-pin it there after this merges (three fields: commit, sha256, `git show` command).
  Note: CI integration `TestManifestRunWithoutGtf` still read `log/analysis.txt`; it now uses `analysis_genes_path`.

- [x] `CLI-MATURE-01` (2026-11-01, `codex/cli-mature`, PR A of the CLI maturity batch; user approved scope): `main()` maps Ctrl-C to exit 130;
  declining the `--reference` prompt and post-parse flag-combination errors exit 1 and 2 (usage errors); `--cores/-c` is the
  primary thread flag on `star-index`, `translated-index`, `translated`, `discover` (`--threads` kept as alias, same dest);
  `translated-index --output/-o` (`--out-dir` alias); `--verbose/--quiet` mutually exclusive everywhere; `--version` and the
  Python >= 3.10 guard run before heavy imports (entry point is now `viralscan.__main__:main`); plain usage off a TTY or with
  `NO_COLOR`; bare `viralscan` prints help to stderr, exit 2; `docs/cli_reference.md` subcommand list is generated. Skipped:
  `evidence_run._die` (own copy, exit 1), grouped subcommand help. No analysis output changed.

- [x] `HELPERS-01` (2026-11-01, user-approved scope; `codex/shared-helpers`): one definition of the shared helpers, no output change.
  `utils.open_maybe_gzip` / `to_int` / `to_float` (`evidence` re-exports them); nine inline gzip-or-open copies and the
  `evidence`/`discovery`/`whitelist_preflight`/`read_filter`/`host_filter` imports now use it; `star_run.plain_whitelist` gains
  `name=` and replaces the copies in `discovery` and `scripts/multimap` (same bytes and file names); `hostresponse` uses
  `run_safety.sha256_file`. Deleted: `utils.load_config`, `anello_align.read_accession_tsv`. Left alone: `evidence` replay
  whitelist copy (text mode, so CRLF would change), `evidence._fasta_digests` (a sequence digest, not a file hash),
  `gzip.open(..., "rt")` calls on always-gzip files. The `write_json` duplicate left with `reference_strategy`.
  `reference_strategy.py` and its three test files left the package (user decision) for the reproducibility repo
  (PR `codex/reference-strategy-move`); the package keeps no importer of it.

- [x] `REPO-SPLIT-4E` (2026-10-27, batch 4e of REPO-SPLIT-01, tracker §18 D61-D77): the package tree keeps package work only.
  Removed here after a byte compare with the reproducibility repo (D65): `.living/findings` (31), `.living/log` (92),
  `.living/HANDOFF_*`, `last-session.md`, `outputs/reviews` (11), 12 `docs/plans` (the 4 cited plans stay), `docs/handoffs`,
  `docs/superpowers`, `PUBLICATION_READINESS_EXECUTION.md`, `skillpacks`, `.github/prompts`. `.living/decisions.md` and
  `learnings.md` keep package and shared entries (experiment entries live in the repro copy). This file keeps open rows,
  section headers, the shared sections and `CAT-31`; the full pre-4e file is frozen at
  `provenance/PLAN_archive_2026-10.md` in the repro. Rehearsed with `audit/rehearse.sh 4e` (0 vanished / 0 regressed).

- [x] `DISC-17` (2026-10-27, repro PR #12 `codex/disc17-metrics`): `scripts/disc12_metrics.py` run on the published DISC-12
  outputs, 182/182 arms ok. The user ratified the five docstring ambiguities ("do all above"). Truth came from each primary
  run's `adata_multimap.h5ad` `counts_unique` layer (= the `unique` rows of `molecule_assignments.tsv.gz`), not from
  `rerun-multimap`, which has no `--multimap-molecule-assignments` flag and cannot run on the pre-`multimap.done` primary
  runs (the earlier row text was wrong). Outcome: rule 2 `none` for every DB/configuration (the covid positive is capture-only
  and its panel is never reached), rule 1 `not_assessed`, rule 3 `not_assessed` (spike-in only), rule 4 reported; no default
  changes proposed. DISC-13 constants stay `calibrated: false`. Tables: repro `docs/plans/disc12/{truth,metrics_2026-10-27}/`.
- [x] `DISC-18` dropped (user, 2026-10-27 grilling Q2): re-rendering `report.html` after discovery would break layout v2,
  and `candidates.tsv` / `candidate_states.tsv` already carry the content; no discovery report page.
- [x] `DISC-19` (2026-10-27, user grilling Q1/Q3/Q5): DISC-12 closed under its frozen rules. Rule 2 found no `k` (the
  capture-only covid positive never reaches its panel), so by the protocol `prioritized` stays unbuilt: the state,
  `PRIORITIZED_MIN_GROUPS` and `discover-states --k` are removed (unreleased since rc1), and Kraken2-only taxa stay
  `observation`. Rule 1 `not_assessed`: the 0.05/4 defaults stay, labelled not calibrated. Test
  `test_kraken2_alone_is_always_observation`. Results note with the fenced no-covid preview (k = 1648 at 0.05/4): repro
  `.living/findings/disc12-benchmark-result-2026-10-27.md`.
- [~] `DISC` Kraken2 cell-resolved taxonomic discovery: plan approved 2026-10-21 (all nine §16 decisions).
  `DISC-R1` (grouping perf) done. `DISC-09` (R2: provenance, fingerprint, publish, resume) done. `DISC-R5` (per-UMI rule: deepest-consistent) done. `DISC-R6` (polyA/T trim before Kraken2) done. `DISC-R7b` (Kraken2 defaults 0.05 / 4) done. `DISC-01` (R3: discovery env, doctor profile) done. `DISC-10` (R4: real-Kraken2 integration test) done. `DISC-11` (R7: `discover-collapse` rank views) done. `DISC-R8` (minimal `candidates.tsv`) done. `DISC-R9` (docs: `docs/discovery.md`) done. `DISC-R10` (EBV real-data validation, 3 arms) done. `DISC-R11` (close-out) done. Defaults 0.05 / 4 confirmed by the user after R10 (2026-10-22). Next: push and PR (user); then `DISC-12`.
- [~] REVIEW-C (WP-3): coverage-map code/fixtures pass; real-directory qualification pending.
  - 2026-11-02 split: the package code and fixtures are done; the real-directory qualification moved to the repro tracker as `REVIEW-C-RUN`. Flip this row to `[x]` when that row closes.
  - Note: real Track A BAMs (`VIRUS|U31784.1`) vs a bare-accession GTF gave `annotation_mismatch`;
    `resolve_reference_annotation` now also accepts `ACC` for viral-class references only and
    records `annotation_seqid` / `annotation_seqid_normalised` in the coverage-map JSON.
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
  **Open (package):** layout review by the user, a `environment.coverage.yml` solve in CI.
  **Moved 2026-11-02 to the reproducibility tracker as `COVMAP-01-RUN`:** a pyGenomeTracks run on a real library beyond HHV-6B AF157706.1. Flip this row to `[x]` when the package items and that row close.
- [~] `REF-11-PREP` (2026-10-23, `codex/panel-archive-script`, user: keep the script): `scripts/make_panel_archive.py` builds the
  deterministic viral-annotation panel archive for Zenodo (zip in the `viralscan data fetch` layout + `PANEL_MANIFEST.json`,
  `SHA256SUMS`, DRAFT `zenodo_metadata.draft.json` with TODO values). Preparation only: nothing uploaded, no DOI reserved.
  Checked: two builds of the local 393-GTF panel at 3.0.0rc1 gave byte-identical zips (sha256 5c5b5fea...49d9), corpus
  check passed; `tests/test_make_panel_archive.py`. **Open (user):** Zenodo deposit and the TODO metadata.
- [x] `RC1-01` (2026-10-23, `codex/release-3.0.0rc1`, decision Q15): version `3.0.0.dev1` -> `3.0.0rc1` in every file
  `scripts/check_versions.py` compares (`__init__`, Dockerfile, Singularity.def, conda recipe, CITATION.cff x2) plus the
  cli_reference and release_policy pins; CHANGELOG section renamed. Historical receipts (release_checks, PLAN history) unchanged.
  **Open (user):** merge after #109 + green main, then push tag `v3.0.0rc1` (TestPyPI per release.yml; not PyPI).
  - 2026-10-23: tag `v3.0.0rc1` (on `87591c4`) pushed; release run 40000132978 built and `twine check`ed both
    distributions, then failed at `check_ship_scope.py` (`jsonschema` missing in the release job; nothing published).
    Fix: release.yml installs `jsonschema`. Re-tag needed (the release runs only on a tag push of a CI-green commit).
  - 2026-11-02 (`codex/rc1-freeze`): closed; the freeze continues as `RC-01`. All versions already read `3.0.0rc1` on main (`check_versions.py` clean, also with `--tag v3.0.0rc1`).
- [x] `PATHS-01` (2026-10-22, `codex/strip-institutional-paths`): `PLAN.md` site paths replaced by `$SITE_ROOT` (archive work root), `$SITE_SCRATCH` (scratch checkout root) and `$COLLAB_ROOT` (collaborator root); 10 lines. **Open:** about 190 more hits in `.living/` and about 300 in tracked analysis outputs, manifests, docs and Slurm scripts (`git grep -c -E "/exports/|lumc"`); scope of a repo-wide scrub is the user's call. **Done (2026-10-22, `codex/strip-paths-docs-scripts`):** tracked `docs/`, `scripts/` and `covid_viralscan/` now use `$SITE_ROOT`/`$SITE_SCRATCH`/`$COLLAB_ROOT` placeholders in prose and the fail-fast env vars `VIRALSCAN_SITE_ROOT`, `VIRALSCAN_SITE_SCRATCH`, `VIRALSCAN_COLLAB_ROOT` in scripts (`#SBATCH` log paths are now relative; `--chdir` lines dropped, submit from the repo root). `.living/`, `analysis/`, `results/` and other outputs deliberately untouched; no history rewrite. **Closed 2026-11-02 (user decision):** no repo-wide scrub. Shipped code and docs are already clean; `.living/` and analysis outputs do not ship.
- [~] `DOCS-PAGES-01` (2026-10-20, `codex/docs-pages`): publish the Sphinx docs on GitHub Pages (user choice over
  Read the Docs). `.github/workflows/docs.yml` builds on push to `main` (and by hand), installs pandoc for nbsphinx
  and `.[docs]`, and deploys with `actions/upload-pages-artifact` v5.0.0 / `actions/deploy-pages` v5.0.1 (SHA-pinned,
  `check_actions_pinned.py` clean). Not built locally (no env has `scanpydoc`). Left for the user: Settings -> Pages ->
  Source "GitHub Actions", push, merge; `[x]` once the first deploy is live. `.readthedocs.yaml` is kept.
- [~] `SCAUTOQC-01` *(2026-11-02 fix: sctk 0.2.0 `fit_gaussian` retries a non-converged GMM with `gmm.fit(x, warm_start=True)`, which sklearn rejects (TypeError; main CI red at de55424). `viralscan.scautoqc._patch_sctk_gmm` swaps in a GMM whose `fit` honours it; `test_sctk_gmm_retry_accepts_warm_start` fails without the shim.)* (2026-10-19, branch `codex/scautoqc`, ported from `codex/scautoqc-wip` 92dd88e): optional 2026-10-22: rebased onto origin/main via tmp/wave2-base-2026-10-22; tests adapted to the `run` subcommand.
  `--scautoqc` SCTK cell/cluster QC flags added in place to the count H5ADs and UMAP outputs; cells are flagged,
  never removed. Opt-in rule/flag/config key are absent when off. Design questions decided by the user 2026-11-02:
  (a) the legacy n_counts/n_genes filter in `umap.py` is UNCHANGED with the flag; (b) the rule now runs AFTER
  detection (which calls the cells and writes `results/called_cells.tsv`) and BEFORE umap, and scores only called
  cells (`add_scautoqc(cells=mask)`); other barcodes keep NA flags and NaN coordinates, and the called-cell count
  is recorded in `adata.uns['scautoqc']` and `log/scautoqc.done` (JSON).
  SCTK ships as the optional extra `viralscan[scautoqc]`, not in `full`. Tests: flag-off, all-cells-fail guard,
  flags-only, missing-sctk error, DAG rule graph (detection -> scautoqc -> umap), called-only scoring. **Package code is done. The real-data gate** (FASTQ-to-UMAP on SRR8315713, clusterwise multi-resolution QC, ribo/globin marker sensitivity)
  **moved 2026-11-02 to the reproducibility tracker as `SCAUTOQC-01-RUN`; flip this row to `[x]` when that row closes.**
  Rebased onto origin/main + shifted-log + tchem (wave-2 prep stack, not upstream): `viralscan.scautoqc` is under mypy
  strict and `viralscan.scripts.scautoqc` under the relaxed Snakemake-script override; legacy filter kept.
  Review follow-up: scoring all barcodes (empty droplets included) was a known flaw; fixed 2026-11-02 by moving the
  rule after cell calling (see the decision above). Added BioMart failure handling + low-mapping warning, case-insensitive `MT-`, extra tests;
  virus-only/no-host index still hard-fails (documented). Clustering path needs `scanpy[leiden]` (the extra).
  - Note (2026-11-02, user decision a + b): legacy filter kept; rule moved after detection, called cells only. Only
    the real-data gate (SRR8315713 FASTQ-to-UMAP, user's run) remains. Detection no longer waits for scautoqc.
  - Note (2026-10-22, PR #100 review): globin markers matched human symbols only, so every mouse `--scautoqc` run
    aborted on the mandatory hb check; `GLOBIN_PATTERN` now also covers mouse (Hba-a1/a2, Hba-x, Hbb-bs/bt/b1/b2,
    Hbb-bh1, Hbb-y, Hbq1a/b), with a regex test that needs no SCTK.
  - Note (2026-10-22, grilling Q8): CI now installs `scanpy[leiden]`, `sctk==0.2.0` and `pybiomart` in the py3.11 Ubuntu
    pytest job only, so `tests/test_scautoqc.py` runs in CI instead of skipping. Not py3.12: sctk 0.2.0's numpy-groupies
    pin fails to build there (`configparser.SafeConfigParser` removed), so the `[scautoqc]` extra is py<=3.11 in practice.
- [x] `PERF-EXEC-02` (2026-10-18, branch `codex/perf-wave1-b24-b23-b04-b06`): wave-1 B24
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
  **Closed 2026-11-01:** wave-1 code is on `main` (`cf84b36`, `2f4e2be`); the timing session and whole-run
  measurement are benchmarking, parked with `PERF-EXEC-04`'s 2026-10-20 deferral. Historical: B08 is
  deferred on its gate and B10 plus wave 2 (B25, B27, B18) are in flight in the
  `20261018-112322-wave2-g00` packet, on `codex/perf-wave2-b25-b27-b18`; that session must re-check
  B18's oracle against `main`'s post-merge `molecule_verdict.py` before trusting it.
- [x] `AUD3-FLAKY` (2026-10-27): not reproduced. `tests/test_cellcalling.py::test_real_r_frozen_protocol_seed_runs_and_is_reduced_modulo_int_max`
  ran 40 times, 8 copies at once in 5 rounds (viralscan_test_full, R/DropletUtils): 40/40 passed in 117-138 s each. The test uses
  per-test `tmp_path` outputs and R's per-process `tempdir()`, so no shared state was found; no code change. Reopen with the failing
  output if it recurs.
- [ ] `AUD3-REF08-DIVERGED` — audit P2-01 residual: the build-time near-duplicate screen (`AUD3-05`) reuses the curation twin cut (>= 95% shared canonical 31-mers), so a host homologue diverged by about 8% (about 90% of its k-mers unshared at k=31) is not screened; `host_homology_status=not_measured` stays the honest label. A homology threshold or alignment-identity cut belongs to `REF-08` (CAT-14 deliberately has none); do not invent one here.
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
- [x] `PERF-EXEC-03` (2026-10-18, original branch `codex/perf-wave2-b25-b27-b18`,
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
  **Closed 2026-11-01:** wave-2 code is on `main` (`5dacb4d`); timing and wave-3 cards are parked with `PERF-EXEC-04`.
- [~] W3-11 SW-09: the 13 parser builders moved verbatim to `src/viralscan/cli/parsers.py` (menu.py 2,530 -> 1,196 lines; `gen_cli_reference.py --check` byte-identical; `menu` re-exports every name). Handlers stay in `menu.py`: they call names the tests patch on `menu` (`_check_required_tools`, `_swap_multimap_layer`, ...), so moving them would break those patches for no gain.
  - [~] M1.3 test tool skips: no change needed. The triage showed
    `test_anellovirus_chain` fails on a stale fixture, not a missing
    dustmasker (ANELLO-14). The 22 tool skips wait on the user's full env.

## Software backlog reconciliation (2026-10-09)

- [x] `QA-02` — `ruff format` the 16 files and fix the 166 mypy errors without blanket ignores. **Format done
  2026-10-09** (15 files; `tests/test_dsr02_enumerate_calls.py` is the user's dirty file and is left
  unformatted; 1,930 unit passes after). **Left:** 166 mypy errors (`evidence.py` 81, `anello_align.py` 30,
  `chemistry.py` 25, `multimap.py` 13, ...), assigned per file to the agent that owns that file.
  **Closed:** mypy errors cleared in `9aacbeb`, `a50eb5c`, `2d66e93`; CI gates `ruff format --check .` and
  `mypy src/viralscan` (green on `main` 2d15473, 2026-11-01).
- [~] `CHEM-01` (closed 2026-10-10: the bulk path stays advice-only), `DEF-05` (done: reader's matrix verified, missing host matrix rejected),
  `DEF-06`, `DSR-08`, `SW-07`, `MECH-F` (index-aware denominators, gene-role contract).
- [~] `REL-04`/`REL-05` repeat on the finished source (clean installs, `doctor`, `validate-run`): **done 2026-10-27 on main
  bf8214b** from a `git archive` export (uv build on local disk): wheel 945,340 B, sdist 1,064,397 B, `twine check` passes both;
  clean uv venvs (Python 3.11) for each, `pip check` clean, `viralscan --version` 3.0.0rc1, `scripts/release_package_smoke.py`
  exits 0 from both (doctor `ok: true`, `validate-run` `ok: true`, `issues: []`, `viralscan` imported from site-packages).
  Still open in this row:
  `REL-03`/`DEF-08`/`CAT-25` (`uv.lock` committed and CI-checked since #141; conda lock replay open); `CAT-06` (index carries its GTF and
  catalogue); `REF-11` and `PANEL-01` WP4 bundling (deferred; the DOI stays external).
- [~] `REF-03` (2026-10-10: `prepare_reference_inputs` now rejects a reference with no viral records when masking
  is requested, with a test; the gtf-manifest check in the NCBI-fetching builder run is network-pending),
  `REF-02`/
  `REF-05`/`REF-07` (real provenance, retrieval vs build time), `REF-04` (deterministic rebuild).
- [~] `PANEL-01` WP2 **done 2026-10-10** (SV40 `NC_001669.1`, AAV1/7/8 `NC_002077.1`/`NC_006260.1`/`NC_006261.1` fetched from NCBI,
  catalogued as `shipped` + `contaminant`, GTFs generated locally and listed in `gtf_corpus_manifest.tsv`; 2,272 tests pass) and **WP3 gates**
  (self-consistency matrix, EC-size gate, `CAT-40` `kallisto inspect`); iciHHV-6 flag and CAR-T vector
  filter under `CAT-04`. B19 genotype 2 stays deferred (no complete reference).
- [~] `REF-13` (done 2026-11-02), `CAT-36`, `CAT-38`/`39` (done 2026-10-10), `DEF-04` (done 2026-11-02: kept and warned) (annotation integrity); stale `AB303562.1` rows in
  `anellovirus_genes.tsv` (removed in `f4320d2`, 2026-10-10); GTF-manifest check in the build lifecycle
  (done: `build_bundled_panel_ref.py:_check_gtf_manifest` runs `write_gtf_manifest.py --check` before any fetch;
  only the network-backed builder run is left).
- [ ] Added from the survey, not in the user's list: DSR tooling that is code (`DSR-02`, `DOSSIER-01..03`,
  `VERDICT-01`, `SENS-07`), `ANDET-07`, `TONSIL-02`, `CAT-19`, `DEF-07`,
  `CAT-14`/`15`/`16`. (2026-11-02: the run rows among these moved to the reproducibility tracker, see "Moved to the reproducibility tracker".)
- [ ] `REL-06..14`, `RC-01..03` (one OCI artifact, SIF from its digest, parity, product-environment audit,
  SHA-pinned Actions, exact-revision gating); `DOC-01` (Sphinx `-W`), `DOC-05`, `DOC-07`, `OPS-01`/
  `OPS-02`, `REL-15`, `STB-04`. Maintainer credentials, `REL-15` and `REF-11` stay user-owned (`[!]`).
- [~] `SW-03` (done 2026-10-10), `MECH-C` (done), `SW-09` (parsers done, handlers deferred), `MECH-E` (open), `CAT-18` (closed, not adopted);
  `EMC-01` (done 2026-10-10; the cluster re-measure is the user's); optional chunked BUS ingestion (`QA-03`, **closed 2026-10-12 without a chunked mode**: ingestion is already one streaming pass over `bustools text` rows (SW-03), so no chunk boundary can split a CB-UMI molecule; the acceptance proof is `tests/test_qa03_bus_ingestion.py` (46 seeded property tests: 7 read-buffer sizes incl. mid-row cuts x every multimap method give identical layers and audit, unsorted text is refused, in-memory row order is irrelevant) on top of the DEF-06 properties in `tests/test_multimapping_properties.py`; a real chunked reader is only worth building if a deep-sample RSS measurement shows the stream is not enough).
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
## Moved to the reproducibility tracker (2026-11-02)

User decision 2026-11-02: experiment work (runs on real or public data, planted benchmark panels, comparator tools, results) is tracked only in
the reproducibility repo's [`PLAN.md`](https://github.com/mdmanurung/ViralScan-reproducibility/blob/main/PLAN.md). The open rows below moved there with their history. This file keeps package code, tests, docs,
packaging, release and CI. An ID in this list that appears elsewhere in this file now lives in the repro tracker.

- Public-data and dataset runs: `RUN-01`..`RUN-05`, `DSR-02`, `DSR-03`, `DSR-04`, `DSR-06`, `DSR-07`, `DSR-10`, `DSR-12`, `DSR-13`, `DSR-16`, `TONSIL-02`, `DISC-01`, `T2-AMEND-01`.
- Comparators: `CMP-00`..`CMP-06`.
- Validation and calibration: `AUD3-G3-EVIDENCE`, `REF-08`, `REF-09`, `REF-10`, `VAL-04`, `VAL-08`, `VAL-10`, `HOST-05`, `SENS-06`, `SENS-07`, `SENS-CAL-01`, `CAL-13`, `VERDICT-01`, `CARD-01`, `ADJ-01`, `DLIST-01`.
- Anellovirus detection: `ANELLO-PRIOR` (with `.1`..`.4`), `ANELLO-12`, `ANELLO-13`, `ANDET-07`, `ANDET-08`, `ANDET-09`, `ANDET-09f`, `HPV-11`.
- Catalogue and reference measurements: `CAT-07`, `CAT-08`, `CAT-10`, `CAT-14`, `CAT-15`, `CAT-16`, `CAT-19`, `CAT-40`.
- Gene programmes: `PROG-15`, `PROG-16`, `PROG-18`.
- Translated search (arm P): `AA-02`, `AA-04`, `AA-05`.
- Split rows: `SCAUTOQC-01` stays here (code done); its real-data gate is `SCAUTOQC-01-RUN` in the repro tracker. Likewise `COVMAP-01` (real-library pyGenomeTracks run: `COVMAP-01-RUN`) and `REVIEW-C` (real-directory qualification: `REVIEW-C-RUN`).
- Checked and not split: `AUD3-REF08-DIVERGED` has no run part. It is a deferral of the homology threshold to `REF-08` (already in the repro tracker), so it stays here as written.

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

No open rows; done rows are in the PLAN archive (see Next action).

## WP0 — Governance and legacy freeze

No open rows; done rows are in the PLAN archive (see Next action).

## WP1 — Core software contracts

### WP1B — Workflow safety and architecture

- [~] `SW-09` — **parsers done 2026-10-10 (`src/viralscan/cli/parsers.py`); handlers deliberately left in `menu.py`, see W3-11.** Split the oversized CLI into thin parsers plus importable service
  functions; convert Snakemake scripts to minimal wrappers without changing
  outputs.
### WP1D — Mechanism consolidation (architecture review, new 2026-09-29)

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
  - 2026-11-02 (`PLAN-RECON-01`): the first three "Still open" items are stale. The catalogue now shows Influenza D per segment, NC_001664 "Human herpesvirus 6A" and KC138720 "Human papillomavirus type XS2". Only the smaller exemplar-name mismatches remain unreviewed.
- [~] `MECH-D` — **tracked under `DEF-02` (reconciled 2026-10-04): one module, one row;
  this closes with DEF-02.** (module landed 2026-10-03) Chemistry module: one geometry for kb, STARsolo and the
  preflight, plus an R1-length / polyT check. The 10x v2 run as `-x 10xv3`
  gave 1.78M "cells" with no error.
- [~] `MECH-E` — the reference pipeline: source adapters → one gate module →
  masking → GTF emitter → index with a recorded D-list. The prototype is
  `viral_panel_max_2026-09-28/01–04`.
  Software side done under `AUD-15` (2026-10-11/12): source adapters, `_prepare_viral_fasta` masking gate, `_emit_viral_gtf` emitter, `_finish_reference_build` index boundary with `genome_dlist_applied` in the manifest. Left: governance re-pin and a real (non-mock) index build.
  - 2026-11-02 (`PLAN-RECON-01`): the governance re-pin is done (repro PR #19 pins package 9e31c65, `reference-manifest-schema` included). Only the real (non-mock) index build is left.
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
- [x] `DEF-04` — *(2026-10-10: no shipped GTF carries a `whole_genome` gene any more, `panel_integrity.py` reports 0; the gene-level exclusion code is deferred: it only affects FASTA-only user references where the whole-genome gene is the virus's only unit, so it needs a design call, not plumbing.)* `whole_genome` tag in the catalogue and its exclusions (Q12). **Closed 2026-11-02 (user decision: keep + warn, `codex/def04-warn`):** the whole-genome gene is kept, since it can be the virus's only unit and excluding it would hide the virus. `prepare_reference_inputs` now logs one warning naming accessions that have only a whole-genome unit; their counts are whole-genome, not per-gene. Pinned by `test_whole_genome_only_virus_is_kept_and_warned`.
  - 2026-10-12 (grill, user-confirmed): deferred past the exploratory public-dataset runs; no shipped reference is affected.
- [~] `DEF-08` — lock `env_full`, and add a snakemake 9 dry-run to CI (Q9).
  - 2026-10-04: CI gains a `snakemake9-dryrun` job (`snakemake>=9,<10` +
    `test_snakefile_dag.py`). It has not run on a runner yet; watch for the
    `connection_pool` build issue. Locking env_full is still open.
  - 2026-10-10 (wave 4): the lock target is `viralscan_test_full`, not
    `benchmark_runs/legacy_v2_v3/env_full` (a stale kb 0.28.2 env). The explicit
    lock is `environment.explicit.txt` (see REL-03). `tests/test_snakefile_dag.py`
    passes 13/13 (incl. the 9 integration-marked dry-runs) with that env's
    snakemake 9.23.1. The CI-runner result of `snakemake9-dryrun` stays the user's.
  - 2026-11-02 (`PLAN-RECON-01`): the CI job `Snakemake 9 dry-run (DEF-08)` runs and passes on every PR. Only the `viralscan_test_full` lock is left.
## WP2 — Install, lock, and artifact parity

### WP2A — Pip and locked environment

- [~] `REL-03` — generate and commit reproducible runtime/development lockfiles
  - 2026-10-10: `uv.lock` committed (`c5c4f9c`); `uv lock --check` passes. Runtime conda lock still open.
  - 2026-10-28 (`codex/lock-refresh`): `uv.lock` was stale since the 2026-10-20/22 `pyproject.toml` changes
    (scautoqc extra, taxonomy); CI never checked it. Regenerated with uv 0.11.14 `uv lock --python 3.11` on local
    disk: 64 packages added (docs and scautoqc extras), none removed, no existing pin moved. The lint job now runs
    `uv lock --check`, so a `pyproject.toml` change without a lock refresh fails CI.
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
  - 2026-10-23 (REL-03-PIPCHECK): the first replay was the v3.0.0rc1 release Docker build (run 40007587542), and its
    `pip check` failed. The lock (and `viralscan_test_full`) held `snakemake-interface-common` 1.23.1, which needs
    `packaging>=26.1`, while snakemake 9.23.1 needs `packaging<26`; `pyseq-align` 1.0.2 needs `cython`, which was absent.
    Fix: pin `snakemake-interface-common-1.23.0-pyh84498cf_0` (`packaging>=24,<26`; satisfies every
    snakemake-interface-* dependent) and add `cython-3.3.0-py311h781f220_0` (deps already locked). Lock is now 488 rows;
    both md5s were checked against the downloaded files. No env replay was done here (installs are the user's job); the
    release Docker job is the replay.
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
  - 2026-11-02 (`codex/rel12-13-ci`): re-checked, all 48 `uses:` (incl. the new attest step) in `ci.yml`, `release.yml`, `docs.yml` and
    `lock-audit.yml` are 40-hex SHAs with `# vX.Y.Z` (`check_actions_pinned.py` clean; `ci-gate` also runs it).
    Added GitHub build-provenance attestations for the wheel, sdist and `SHA256SUMS` in the `publish` job
    (`actions/attest-build-provenance` v4.2.2, `attestations: write`; workflow permissions only, no repo setting).
    Still `[~]`: OCI image, SIF and lock attestations are not added (the image digest exists only in the `container`
    job; needs a tag run to verify), and the attest step is untested until a tag run. Maintainer: verify with
    `gh attestation verify` after the first release; the Node 20 action majors stay as noted above.
- [x] `REL-13` — make release publication depend on green CI for the exact tagged
  SHA; build distributions once and make all downstream jobs consume them.
  - 2026-11-02 (`codex/rel12-13-ci`): already in `release.yml` (landed with the release-CI work; row was stale). `ci-gate`
    runs `scripts/check_release_ci.py` (latest `ci.yml` push run on `main` for `git rev-parse HEAD` must be
    completed/success, run id uploaded as `release-ci.json`); `build` builds and tests wheel+sdist once; `container`,
    `publish` and `release-assets` consume the `dist` artifact after `sha256sum -c`. Pinned by
    `tests/test_release_static.py`. Limit: the gate needs a CI run on the tagged SHA, so tag a commit already on `main`.
- [ ] `REL-14` — add protected TestPyPI/PyPI environments, stable approval, GitHub
  release assets, and separate RC/stable container tag behavior.
- [!] `REL-15` — maintainers must provide final author, licence-holder,
  maintainer/contact, ORCID, affiliation, support-window, and trusted-publisher
  metadata.
## WP3 — Preregister scientific validation

No open rows; done rows are in the PLAN archive (see Next action).

## WP4 — Freeze production references

### WP4A — Profiles and provenance

- [~] `REF-02` — fail-closed fetches and manifests exist; complete accession
  version, taxonomy, snapshot, retrieval date, SHA-256, length, licence, cluster,
  representative status, rationale, and missing-accession fields.
  - Blocked (2026-10-02): the licence fields need `REF-05`. Wiring the manifest into the production builder is `REF-03` (`[ ]`).
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
- [x] `REF-13` — **new, 2026-09-26.** Detection-side reference visibility and
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
  - 2026-11-02 (user decision, `codex/ref13-build-reference`): **wired into `build-reference`.** `build_reference._emit_viral_gtf` now runs
    `normalise_gtf_file` on the viral GTF it writes (combined-reference and anellovirus builds), so `kb ref` sees exon rows for CDS-only
    genes. Test: `test_emitted_viral_gtf_is_normalised_for_kb`. Shipped GTFs untouched. `gene_programs.tsv` regenerated from the 393 local
    GTFs with `extras/build_gene_programs.py`: output byte-identical to the committed file (it already carries the REF-13 groups, so the
    10-10 "not regenerated" note was stale).
  - 2026-11-02 (user decision, same PR): **`viralscan --reference <fasta> <gtf>` wired too.** `menu._normalise_kb_ref_gtf` writes
    `index/reference.normalised.gtf` and kb indexes that; the user's file is never modified. Host cDNA rows (seqname == transcript_id !=
    gene_id) are copied verbatim; only the other rows are normalised. `content_files.gtf` still hashes the user's GTF; new
    `content_files.normalised_gtf` hashes the indexed copy (schema allows extra keys, so no schema change; a changed user GTF still changes
    both hashes and `content_sha256`). Test: `test_reference_path_indexes_a_normalised_gtf_copy`; mutation check (normalise call removed)
    fails both new tests. **Remaining:** HBV P/S origin-wrapping first part stays unindexed (ngs_tools limit, no package fix).
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
### WP4B2 — Host-filtering design and reporting (new 2026-09-26)

Open rows moved to the reproducibility tracker (see "Moved to the reproducibility tracker").

## WP4C — Detection sensitivity and negative-result claims (new 2026-09-26)

- [x] `SENS-CORR-02` — exact substitution-only model. **Done (2026-10-06):**
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
  - 2026-11-02 (`PLAN-RECON-01`): done. `docs/output_reference.md:502` gives the exact values (0.708 / 0.254 / 0.063 / 0.013) and `:515-540` documents the `positive_control.json` fields.
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
## WP4D — Gene-programme inference, layer 2 (new 2026-09-26)

- [~] `PROG-09` — **open.** The catalogue is 9 viruses. Polyomaviruses (JC/BK/KI/WU,
  MCPyV), HPV, HBV, HDV, GBV-C, HIV-1 and HTLV-1 are deliberately excluded: they
  have no latency/lytic dichotomy representable from the panel's protein-coding
  genes, and for HIV/HPV/HBV the interesting state is DNA-level latency, which
  scRNA-seq cannot observe at all. Adding them would mean emitting
  `not_applicable` rows, which is what the code does for a detected virus with no
  model — a decision to make explicitly rather than by omission.
## WP4E — Named HPV ORFs for the oncogene-versus-capsid contrast (new 2026-09-26)

Open rows moved to the reproducibility tracker (see "Moved to the reproducibility tracker").

## WP4F — Real gene structure for the Anelloviridae panel (new 2026-09-26)

### Measured CDS coverage (full panel, 2,042 accessions)

Open rows moved to the reproducibility tracker (see "Moved to the reproducibility tracker").

## WP4G — Anellovirus detection you can trust (new 2026-09-27)

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
- [x] `ANDET-04` — evidence replay reads the raw FASTQs (`evidence_run.py:154-155`,
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
  - 2026-11-02 (`PLAN-RECON-01`): done as an opt-in rule (`--auto-evidence`, Snakefile `rule auto_evidence`). Default-on stays the user's call under `ANELLO-PRIOR.3`.
## WP4H — Comprehensive human-virus catalogue (new 2026-09-27)

- [x] `CAT-04` — risk classes: exclude human endogenous retroviruses; flag
  integrated ciHHV-6, `EVE_RISK_GENERA`, and vector/reagent contaminants.
  - 2026-10-10 (wave 4 reconcile): the **flags are done** in `f4320d2`: catalogue `role` (`contaminant` AAV1-8 and SV40,
    `decoy` MLV/XMRV, `endogenous` HERV-K113) and `risk_class` feed `reference_risk_flags`
    (`virus_grouping.py:261-270`: `endogenous_reference`, `vector_reagent_reference`, `endogenous_overlap_possible`,
    `iciHHV6_possible`) into `viral_summary.tsv` and the report. **Open, needs a decision/data:** HERV exclusion (HERV-K113 is
    `broad`, not shipped, so nothing to exclude yet), a CAR-T/lentiviral vector filter (no vector sequences in the repo), and a
    real ciHHV-6 integration test (today the flag is name-based).
  - Closed as deferred (user, 2026-11-02): reopen when a ciHHV-6 carrier library / vector sequences exist. No data-driven rule now; the name-based `iciHHV6_possible` flag stays (pinned in `tests/test_panel_integrity_contract.py`). HERV exclusion and the vector filter wait for data too.
### WP4H expansion — breadth and measured diversity (new 2026-09-27)

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
### WP4H expansion correction and additions (new 2026-09-27)

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
- [x] `CAT-30` — *(2026-11-02, user-approved: retracted in ViralScan-reproducibility PR #20, `covid_viralscan/RETRACTION_NOTICE_CAT-30.md`; registry claim `ttv-biological-signal` was already `rejected`.)* **retract the covid Alphatorquevirus claim.** `main` @ `4fcd748`
  closed F-005 — *"All 8 detected anellovirus accessions appear in Phase A
  (multi-chromosomal GRCh38 alignment) — accession-level artifact confirmation"*,
  *"no genuine viral infection in these COVID PBMC samples"* — but
  `covid_viralscan/results/*/results/viral_summary.tsv` **still publishes**
  `Alphatorquevirus` 1,167,103.0 UMI (x213-g, x216-g.jul2), 1,605,631.0
  (x216-g), 17.63–19.05 % of cells. The science is settled; the retraction never
  happened. **Blocked pending explicit user approval** — do not edit published
  results unilaterally. `CAT-17` makes this more urgent, not less: the shipped
  panel is the one carrying 170 pure-homopolymer k-mers.
## WP4I — Complete latent/lytic state calling (new 2026-09-27)

Open rows moved to the reproducibility tracker (see "Moved to the reproducibility tracker").

## WP4J — SFL tonsil TTV and HPV screen (new 2026-09-27)

Open rows moved to the reproducibility tracker (see "Moved to the reproducibility tracker").

## WP4K — Catalogue↔index reconciliation and targeted reference adds (new 2026-09-27)

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
- [~] `CAT-36` — *(2026-10-10: reverse gap is zero, every GTF seqname is catalogued with `source`/`refseq`; the
  `oncogenic_class` ask is dropped because `hpv_genes.HIGH_RISK_GENOTYPES` already carries it; the 96-row
  count is measured on a built index and is not re-checked here.)* **the catalogue is now behind the panel by 96 rows.** The build
  reports 96 genomes indexed but absent from `virus_catalog.tsv` (the INSDC HPV
  set, plus `AF157706.1` HHV-6B). They have no provenance, family or tier, and
  `CAT-31`'s reverse check will keep reporting them. Extend the catalogue, and
  propagate `source` / `refseq` / `oncogenic_class` for the HPV rows.
## WP4L — Audit 2026-10 findings (new 2026-10-11)

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
- [~] `AUD-15` — `MECH-E` reference pipeline refactor (user-approved for autonomous work 2026-10-11): source adapters, one gate, masking, GTF emitter; touches pinned `build_reference.py`; mock-verified only.
  **Decision (user, 2026-10-11):** refactor plus manifest schema additions allowed (D-list record, source adapters; bump `reference_manifest.schema.json`, governance re-pin, docs); a golden test must prove byte-identical FASTA/GTF/t2g on fixture inputs for everything not intentionally changed; mock-verified only. Also fold in D-13 (the anellovirus origin-wrap gene is transposed in the `kb ref` cDNA).
  Partial 2026-10-11: one `_prepare_viral_fasta` masking gate, `_finish_reference_build` index/receipt boundary, `source_adapters` and `genome_dlist_applied` in the manifest (`reference_manifest.schema.json` now accepts 3.0.0 and 3.0.1; builder writes 3.0.1), golden test `tests/test_reference_pipeline_aud15.py` (13 cases, byte-identical FASTA/GTF/t2g on fixtures). Then ab9d9da: the GTF emitter is one boundary, `build_reference._emit_viral_gtf` (catalogue and combined adapters; `_anellovirus_gtf` and `_gtf_from_merged_fasta` are thin wrappers). **D-13 closed as a documented limitation, 2026-10-12 (`codex/aud-15c`):** `ngs_tools` (the splitter `kb ref` runs) merges exons into a sorted `SegmentCollection`, so `join(2677..2824,1..80)` of KU243129.1 is indexed as `[1:80]+[2677:2824]` whatever order the GTF gives; this cannot be fixed in ViralScan's GTF without rotating the genome or splitting the gene, both of which change coordinates and the transcript model, so it needs an `ngs_tools`/kb-python change. The rows carry `viralscan_norm "wrap_native"` (as EBV LMP-2A/2B do) and `test_origin_wrap_gene_is_tagged_and_kb_splitter_transposes_it` runs the same splitter on a fixture, keeps the GTF in transcript order (ANELLO-05, never re-sorted) and fails if kb-python ever fixes the order. Impact: one gene, junction-spanning reads only. Still open: governance re-pin (user, once at the end).
- [~] `AUD-16` — `ANDET-01/02/04` follow-ups (user-approved 2026-10-11): per-accession genome-coverage breadth, host-homology coordinates kept for later demotion; auto-evidence stays opt-in.
  **Decision (user, 2026-10-11):** ANDET-02 = keep per-locus host-homology coordinates and add an `eve_risk` column, NO demotion until REF-08 calibrates thresholds; ANDET-01 per-accession genome-coverage breadth and ANDET-04 as in the row (auto-evidence stays opt-in).
  Partial 2026-10-11 (`codex/aud-a1`): ANDET-01 per-virus best-genome breadth done (`alignment_breadth_max`, `alignment_breadth_unique_max`, circular-safe; reporting only). Still open: ANDET-02 per-locus host-homology coordinates and `eve_risk` column; ANDET-04 as in the row.
  ANDET-02 coordinates done 2026-10-11: the builder keeps every host-homology locus (`host_homology_loci` in the reference manifest, PAF validated), detection writes `results/host_homology_loci.tsv` with an `eve_risk` column; no demotion (waits for REF-08). Open: ANDET-04 as in the row.
  2026-10-12 (`codex/aud13-andet-followups`): `results/host_homology_loci.tsv` is now documented (pinned `output_reference.md`, re-pin pending) and a test pins that every column is. The row stays `[~]` (ANDET-04 opt-in, demotion on REF-08).
- [x] `AUD3-02-OPEN` — open item: the viral-free test is still not calibrated. On pure-noise host cells with random virus-positive labels the audit saw p < 0.05 in 4/15 runs without viral genes in the embedding (nominal 5%), cause unexplained; shown by `biomni_bundle/plans3/verification/repro/P5-01_viral_in_embedding.py` (the `viral_in_embedding=False` rows). The test is descriptive only until this is understood; do not read p < 0.05 as evidence of anything.
  - 2026-11-02 (`codex/aud3-02-calibration`): cause found and already fixed. The audit's "viral-free" graph (`87c9e35`, PCA patched to drop viral HVGs) still computed size factors over viral counts, so viral load moved positive cells in host PCA space. `host_only_pca` (P5-07) drops viral columns first. Simulation, 200 null seeds per load: old graph rejects at 0.040 / 0.065 / 0.715 (low / moderate / dominant), current graph at 0.050 in every mode. New test `tests/test_umap.py::TestViralFreeTestGraph::test_null_rejection_rate_is_nominal` (60 seeds, <= 8 rejections; 60/60 with the size-factor leak restored). The report caveat "the null is not calibrated" is replaced. A small p stays descriptive about infection (ambient RNA, uptake, doublets).
## WP5 — Build and validate the truth panel

### WP5A — Generator

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
## WP6 — Run public positives and comparators

### WP6A — Public v3 reruns

Open rows moved to the reproducibility tracker (see "Moved to the reproducibility tracker").

### WP6B — Harmonized workflow matrix

Open rows moved to the reproducibility tracker (see "Moved to the reproducibility tracker").

### WP6C — Dataset run plan (new 2026-10-05)

- [x] `DLIST-02` — host-capture split as `--host-capture`: **removed (user decision, 2026-11-02, `codex/dlist02-remove`): gate not met; recover from commit ca2ce7f.** DLIST-01 I2 never ran, the feature never ran end to end on real kallisto, and the TSO finding says a D-list cannot fix F-028. Removed: the option, its `run_safety.omit_when_unset` entry (default fingerprints unchanged), the `host_capture` rule and config key, `scripts/host_capture.py`, `evidence.replay_host_capture_bus`, `tests/test_host_capture.py`, the ship-scope entries and the doc rows. `evidence._replay_quant_bus` stays (`replay_exact_target_bus` uses it).
## WP7 — Score, audit, and freeze results

No open rows; done rows are in the PLAN archive (see Next action).

## WP8 — Regenerate documentation and claim evidence

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
  - 2026-11-02 (`PLAN-RECON-01`): Sphinx `-W` now runs in CI (`.github/workflows/docs.yml:33`). Only the REF-11 data-fetch wording (user) is left.
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
  - 2026-11-02 (`codex/doc10-gate`, user decision: triage, then make the gate fail): the checker moved to
    ViralScan-reproducibility in 40d1703, so the gate itself is a repro change, made at the next re-pin. Docs
    pass done here: `--report-unmarked` flagged 53 lines before, 26 after. Registered markers added
    (`lod95-poisson-depth-floor`, `capture-cliff-substitution-model`) in CHANGELOG and `output_reference.md`;
    host-filter `p_unique` table, "3.5x", "lost 96.5%", "85-93%", TTV k-mer share, "~8 % of the human genome",
    the -22 % shift and the "~3x / ~6x" performance wording were softened (magnitudes not registered).
    The 26 left are parameters, thresholds, definitions and worked examples; the checker has no allow-list,
    so the repro must allow-list them by text anchor (list in the PR) before the gate can fail.
  - Not done: the repro gate is not failing yet; installation claims still have no registry artifact (the
    DOC-05 run log is not committed or inventoried).
  - 2026-11-02 (`codex/changelog-4pct`): a reviewer read the `CHANGELOG.md` "4.0% to 0.2%" line as 50 of 1,000 = 5 %.
    The line is correct, so it is unchanged and the repro allow-list anchor stands. The fixture
    (`test_complete_universe_legacy_and_fractional`, commit 3f4f99b, since moved to the repro) has 50 viral-positive
    cells, of which 2 reach 10 molecules (the other 48 have 9.6). Legacy denominator = 50 viral-positive called
    cells: 2/50 = 4.0 %. New denominator = 1,000 called cells: 2/1,000 = 0.2 %. The 50 is a denominator, not the numerator.
## WP9 — Regenerate and review the manuscript

No open rows; done rows are in the PLAN archive (see Next action).

## WP10 — Publish and test `3.0.0rc1`

- [~] `RC-01` — freeze `3.0.0rc1` version, changelog, migration notes, CFF, locks,
  reference manifest, checksums, SBOM, licence report, and release notes on one
  green protected-main SHA.
  - 2026-11-02 (`codex/rc1-freeze`, PR only, nothing tagged): version strings (`__init__`, Dockerfile, Singularity.def,
    conda recipe, CITATION.cff x2, `cli_reference.md`, `release_policy.md`) all read `3.0.0rc1`. CHANGELOG: the
    "Unreleased (after 3.0.0rc1)" and `[Unreleased]` sections are folded into one `3.0.0rc1` section, because the
    `v3.0.0rc1` tag (now on `4da4674`, 164 commits behind main) never produced a release. `docs/migration.md` already
    covers layout v3. `uv.lock` unchanged (version is dynamic; no bump needed). Checksums, SBOM and licence report are
    produced by `release.yml` on the tag run, not committed.
    **Open (maintainer):** (1) the last release run 40038144113 failed at "Audit the product image before
    publication" (anchore, severity cutoff high), before any push, so nothing is on GHCR/TestPyPI; triage the findings
    or set an exception. (2) The stale remote tag `v3.0.0rc1` must be deleted and re-created on the chosen green main
    SHA (`RC-03`); this PR does not touch tags. (3) CITATION.cff `date-released`, authors and contact (`REL-15`).
    (4) Reference archive manifest (`REF-11`).
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
## WP11 — Publish stable `3.0.0` and submit

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

No open rows; done rows are in the PLAN archive (see Next action).

## VR pilot — read cleanup and transcript reconstruction (2026-10-14)

No open rows; done rows are in the PLAN archive (see Next action).

## Reviewed roadmap implementation — 2026-10-20

No open rows; done rows are in the PLAN archive (see Next action).

## Kraken2 taxonomic discovery (DISC, new 2026-10-21)

No open rows; done rows are in the PLAN archive (see Next action).

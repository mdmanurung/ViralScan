# ViralScan 3.0 execution tracker

## Next action

Live tracker = open rows (`[ ]`, `[~]`, `[!]`) with their section headers, the shared contract sections
(How to use, Critical path, Gate dashboard, External inputs, Stop rules) and `CAT-31`. Every done row and the
dated "Next action" log up to 42b59f8 are frozen in the reproducibility repo:
[`provenance/PLAN_archive_2026-10.md`](https://github.com/mdmanurung/ViralScan-reproducibility/blob/main/provenance/PLAN_archive_2026-10.md) (REPO-SPLIT-01 batch 4e, D63/D67/D68).

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
- [~] `PATHS-01` (2026-10-22, `codex/strip-institutional-paths`): `PLAN.md` site paths replaced by `$SITE_ROOT` (archive work root), `$SITE_SCRATCH` (scratch checkout root) and `$COLLAB_ROOT` (collaborator root); 10 lines. **Open:** about 190 more hits in `.living/` and about 300 in tracked analysis outputs, manifests, docs and Slurm scripts (`git grep -c -E "/exports/|lumc"`); scope of a repo-wide scrub is the user's call. **Done (2026-10-22, `codex/strip-paths-docs-scripts`):** tracked `docs/`, `scripts/` and `covid_viralscan/` now use `$SITE_ROOT`/`$SITE_SCRATCH`/`$COLLAB_ROOT` placeholders in prose and the fail-fast env vars `VIRALSCAN_SITE_ROOT`, `VIRALSCAN_SITE_SCRATCH`, `VIRALSCAN_COLLAB_ROOT` in scripts (`#SBATCH` log paths are now relative; `--chdir` lines dropped, submit from the repo root). `.living/`, `analysis/`, `results/` and other outputs deliberately untouched; no history rewrite.
- [~] `DOCS-PAGES-01` (2026-10-20, `codex/docs-pages`): publish the Sphinx docs on GitHub Pages (user choice over
  Read the Docs). `.github/workflows/docs.yml` builds on push to `main` (and by hand), installs pandoc for nbsphinx
  and `.[docs]`, and deploys with `actions/upload-pages-artifact` v5.0.0 / `actions/deploy-pages` v5.0.1 (SHA-pinned,
  `check_actions_pinned.py` clean). Not built locally (no env has `scanpydoc`). Left for the user: Settings -> Pages ->
  Source "GitHub Actions", push, merge; `[x]` once the first deploy is live. `.readthedocs.yaml` is kept.
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
- [x] `AUD3-FLAKY` (2026-10-27): not reproduced. `tests/test_cellcalling.py::test_real_r_frozen_protocol_seed_runs_and_is_reduced_modulo_int_max`
  ran 40 times, 8 copies at once in 5 rounds (viralscan_test_full, R/DropletUtils): 40/40 passed in 117-138 s each. The test uses
  per-test `tmp_path` outputs and R's per-process `tempdir()`, so no shared state was found; no code change. Reopen with the failing
  output if it recurs.
- [ ] `AUD3-G3-EVIDENCE` — evidence made before `AUD3-04` (the 35 `trackA_round1` evidence directories and the `dsr_round2` calls) used `minimap2 -ax sr` without secondary alignments, so it has no `tie` verdicts and its manifests lack `minimap2_options`. Replaying the 15 saved `dsr_round2` calls showed no change in verdict counts or passing molecules (none holds an exact tie), but re-run the evidence before `molecule_verdicts.tsv.gz` feeds tier measure S at G3 preparation.
- [ ] `AUD3-REF08-DIVERGED` — audit P2-01 residual: the build-time near-duplicate screen (`AUD3-05`) reuses the curation twin cut (>= 95% shared canonical 31-mers), so a host homologue diverged by about 8% (about 90% of its k-mers unshared at k=31) is not screened; `host_homology_status=not_measured` stays the honest label. A homology threshold or alignment-identity cut belongs to `REF-08` (CAT-14 deliberately has none); do not invent one here.
- [!] `AA-02` — **Deferred 2026-10-20 (user parked translated-search work behind Track A on 2026-10-19; a full-depth run is not authorized).** user-authorized bounded exploratory sanity check: build the
  human-genome/cDNA-masked protein index, then search the first 1,000,000 intact
  paired reads from public EBV-positive SRR12682296 (10xv2) and HSV-1-positive
  SRR8315713 (Drop-seq). Keep immutable subsets, input/reference hashes,
  exploration ledger, expected-virus and off-target summaries. No full-depth
  analysis, held-out negatives or calibration is authorized here. (The Codex-era "no
  commit/push" line is superseded for the code only; see the takeover note in the
  "VR pilot" section.) **Note 2026-10-14 (`AA-04`):** the index behind these counts used the plain, untranslated
  host mask, so its host masking was largely inactive; the AA-02 numbers are unmasked.
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
- [~] W3-11 SW-09: the 13 parser builders moved verbatim to `src/viralscan/cli/parsers.py` (menu.py 2,530 -> 1,196 lines; `gen_cli_reference.py --check` byte-identical; `menu` re-exports every name). Handlers stay in `menu.py`: they call names the tests patch on `menu` (`_check_required_tools`, `_swap_multimap_layer`, ...), so moving them would break those patches for no gain.
  - [~] M1.3 test tool skips: no change needed. The triage showed
    `test_anellovirus_chain` fails on a stale fixture, not a missing
    dustmasker (ANELLO-14). The 22 tool skips wait on the user's full env.

## Software backlog reconciliation (2026-10-09)

- [~] `QA-02` — `ruff format` the 16 files and fix the 166 mypy errors without blanket ignores. **Format done
  2026-10-09** (15 files; `tests/test_dsr02_enumerate_calls.py` is the user's dirty file and is left
  unformatted; 1,930 unit passes after). **Left:** 166 mypy errors (`evidence.py` 81, `anello_align.py` 30,
  `chemistry.py` 25, `multimap.py` 13, ...), assigned per file to the agent that owns that file.
- [~] `CHEM-01` (closed 2026-10-10: the bulk path stays advice-only), `DEF-05` (done: reader's matrix verified, missing host matrix rejected),
  `DEF-06`, `DSR-08`, `SW-07`, `MECH-F` (index-aware denominators, gene-role contract).
- [~] `REL-04`/`REL-05` repeat on the finished source (clean installs, `doctor`, `validate-run`): **done 2026-10-27 on main
  bf8214b** from a `git archive` export (uv build on local disk): wheel 945,340 B, sdist 1,064,397 B, `twine check` passes both;
  clean uv venvs (Python 3.11) for each, `pip check` clean, `viralscan --version` 3.0.0rc1, `scripts/release_package_smoke.py`
  exits 0 from both (doctor `ok: true`, `validate-run` `ok: true`, `issues: []`, `viralscan` imported from site-packages).
  Still open in this row:
  `REL-03`/`DEF-08`/`CAT-25` (review the untracked `uv.lock`); `CAT-06` (index carries its GTF and
  catalogue); `REF-11` and `PANEL-01` WP4 bundling (deferred; the DOI stays external).
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
- [ ] Added from the survey, not in the user's list: DSR tooling that is code (`DSR-02`, `DOSSIER-01..03`,
  `VERDICT-01`, `SENS-07`, `DLIST-02`), `ANDET-07`, `TONSIL-02`, `CAT-19`, `DEF-07`,
  `CAT-14`/`15`/`16`.
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
- [~] `DEF-04` — *(2026-10-10: no shipped GTF carries a `whole_genome` gene any more, `panel_integrity.py` reports 0; the gene-level exclusion code is deferred: it only affects FASTA-only user references where the whole-genome gene is the virus's only unit, so it needs a design call, not plumbing.)* `whole_genome` tag in the catalogue and its exclusions (Q12).
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
- [ ] `CMP-06` (WP6B) — same-input comparison with evonk's original
  ViralScan 2.2.0 on the **same latest reference** as v3, so only the
  implementation differs (user request 2026-09-30).
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
- [ ] `REL-13` — make release publication depend on green CI for the exact tagged
  SHA; build distributions once and make all downstream jobs consume them.
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
### WP4B2 — Host-filtering design and reporting (new 2026-09-26)

- [ ] `HOST-05` — **open.** The sensitivity cost of the v3 STARsolo path is still
  unmeasured. It blocks `MS-02`/`CMP-01`–`CMP-03` and the preregistered
  `D15`/`D16` endpoints; no number should be quoted for it until then.
  - [x] `HOST-05a` (2026-10-27, `codex/host-loss`): `viralscan host-loss --run-dir SAMPLE --truth TSV [--viral-label]
    [--mixed-manifest] [-o]` exposes `host_filter.lost_truth_counts` on a run's `host_filtered/fragment_lineage.tsv.gz`,
    writing `posthoc/host_loss/host_loss.json` with `lower_bound: true` and `not_estimable` when truth reads are missing
    from the lineage. Tests `tests/test_host_loss.py` (2).
  - [x] `HOST-05b` (2026-10-27): `tests/integration/test_host_filter_star.py` parametrised over `STAR_FILTER_PARAM_SETS`
    (`pinned`, `star-default`); the spiked viral read survives both (D15 loss 0, nothing unseen). Runs in CI once STAR-05
    lands.
  - [x] `HOST-05d` (2026-10-27): docs state the loss as a lower bound (`output_reference.md`, `cli_reference.md`).
  - Still open (needs real-data compute): `hosthomolog.sbatch` on the Track A libraries (DISC-12-HH), phase-6 EBV
    acceptance, CMP-01..03. HOST-05 itself stays `[ ]`.
## WP4C — Detection sensitivity and negative-result claims (new 2026-09-26)

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
## WP4D — Gene-programme inference, layer 2 (new 2026-09-26)

- [~] `PROG-09` — **open.** The catalogue is 9 viruses. Polyomaviruses (JC/BK/KI/WU,
  MCPyV), HPV, HBV, HDV, GBV-C, HIV-1 and HTLV-1 are deliberately excluded: they
  have no latency/lytic dichotomy representable from the panel's protein-coding
  genes, and for HIV/HPV/HBV the interesting state is DNA-level latency, which
  scRNA-seq cannot observe at all. Adding them would mean emitting
  `not_applicable` rows, which is what the code does for a detected virus with no
  model — a decision to make explicitly rather than by omission.
## WP4E — Named HPV ORFs for the oncogene-versus-capsid contrast (new 2026-09-26)

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

### Measured CDS coverage (full panel, 2,042 accessions)

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

- [~] `CAT-04` — risk classes: exclude human endogenous retroviruses; flag
  integrated ciHHV-6, `EVE_RISK_GENERA`, and vector/reagent contaminants.
  - 2026-10-10 (wave 4 reconcile): the **flags are done** in `f4320d2`: catalogue `role` (`contaminant` AAV1-8 and SV40,
    `decoy` MLV/XMRV, `endogenous` HERV-K113) and `risk_class` feed `reference_risk_flags`
    (`virus_grouping.py:261-270`: `endogenous_reference`, `vector_reagent_reference`, `endogenous_overlap_possible`,
    `iciHHV6_possible`) into `viral_summary.tsv` and the report. **Open, needs a decision/data:** HERV exclusion (HERV-K113 is
    `broad`, not shipped, so nothing to exclude yet), a CAR-T/lentiviral vector filter (no vector sequences in the repo), and a
    real ciHHV-6 integration test (today the flag is name-based).
- [ ] `CAT-07` — diversity-aware representatives for high-diversity families,
  chosen by leave-one-out k-mer capture (the `REF-01` method), not one exemplar.
- [ ] `CAT-08` — build with native `viralscan build-ref` and measure leave-one-out
  confusability per family: where reads from non-indexed isolates land. Also the
  valid `HPV-11` test.
### WP4H expansion — breadth and measured diversity (new 2026-09-27)

- [ ] `CAT-10` — Tier 2 persistence/reactivation set, curated from primary
  literature with a citation per entry and checked with `verify-references`
  before commit. Each row carries a `persistence_class` (`latent-episomal`,
  `latent-integrated`, `chronic-productive`, `persistent-commensal`,
  `recurrent-lytic`). Set: the 9 human herpesviruses; polyomaviruses (BK, JC,
  MC, WU, KI, TSPyV, HPyV6/7); HIV-1/2 and HTLV-1/2; HBV and HDV; HPV; the
  anelloviruses; parvovirus B19; adenovirus; measles (SSPE); enterovirus; HCV;
  pegivirus.
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
### WP4H expansion correction and additions (new 2026-09-27)

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
## WP4I — Complete latent/lytic state calling (new 2026-09-27)

- [ ] `PROG-15` — measure before adding antisense latency transcripts (VLT, LAT
  intron, LUNA): strandedness (`kallisto bus` runs with no strand flag), the share
  of EBV LCL reads outside annotated exons, and which markers are real index
  targets (the KSHV GTF has 26 exon rows for 96 genes).
- [ ] `PROG-18` — (found 2026-10-03 with PROG-08) validate the KSHV latent
  call on a real dataset: a public PEL or KS scRNA-seq set (BCBL-1/BC-3 PEL
  lines are latently infected). Expect mostly latent, with a lytic minority.
  Not started.
- [ ] `PROG-16` — lytic acceptance test: KSHV `GSE190558` (`RUN-04`), induced vs
  uninduced; HSV-1 `SRR8315713` expected productive with no latent calls.
## WP4J — SFL tonsil TTV and HPV screen (new 2026-09-27)

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
- [ ] `AUD3-02-OPEN` — open item: the viral-free test is still not calibrated. On pure-noise host cells with random virus-positive labels the audit saw p < 0.05 in 4/15 runs without viral genes in the embedding (nominal 5%), cause unexplained; shown by `biomni_bundle/plans3/verification/repro/P5-01_viral_in_embedding.py` (the `viral_in_embedding=False` rows). The test is descriptive only until this is understood; do not read p < 0.05 as evidence of anything.
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
## WP6 — Run public positives and comparators

### WP6A — Public v3 reruns

- [ ] `RUN-01` — HHV-6B `SRR20710641` under measured 10x v3 geometry.
- [ ] `RUN-02` — EBV `SRR12682296` under measured 10x v2 geometry.
- [ ] `RUN-03` — HSV-1 `SRR8315713` under measured Drop-seq geometry.
- [ ] `RUN-04` — KSHV latent/lytic series `GSE190558`.
- [ ] `RUN-05` — KSHV plus EBV replicates `GSE154900`.
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
### WP6C — Dataset run plan (new 2026-10-05)

- [ ] `DSR-10` — standard arm set on every dataset (user, 2026-10-06; see
  `.living/decisions.md` "Same standard arm set on every dataset"). Each
  dataset × arm cell is run or marked n/a with a reason in the run plan's arm
  matrix. The two-step arm needs `DSR-11`.
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
- [ ] `DSR-06` — GSE164690 HNSCC one-lane pilot, HN18 CD45− vs HN01.
- [ ] `DSR-07` — E-MTAB-13687 tonsil atlas negative-control subset, samples
  chosen from metadata before any viral output.
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
  - Not done: the `--report-unmarked` gate stays report-only, because turning it into a failing check
    needs the triage above (a user decision); and installation claims have no registry artifact: the
    DOC-05 run log is not committed or inventoried.
## WP9 — Regenerate and review the manuscript

No open rows; done rows are in the PLAN archive (see Next action).

## WP10 — Publish and test `3.0.0rc1`

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

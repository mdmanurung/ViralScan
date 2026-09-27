# ViralScan 3.0 execution tracker

Status: **active**

Branch: `codex/viralscan-v3`

Last reconciled: 2026-08-08

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
- [~] `GOV-06` — execute the outcome-ineligible ViralScan 2.2.0 versus v3
  diagnostic in `analysis/legacy_v2_v3/`. The identical-BUS arm is complete
  (44/44 valid rows) and the five-control input gate is closed. Fresh
  matched-FASTQ attempt-2 arrays `25331035` and `25331037` are **terminal with
  all ten tasks failed** across three independent causes: a nested-output
  validator defect that manufactured exit 65 on four genuinely successful v2
  rows, an unpopulated Zenodo viral-data cache that killed all five v3 rows
  before quantification, and one genuine v2 out-of-memory on `SRR6825024`.
  Attempt 3 is fully wired and tested but **not frozen**, because the cache pin
  it requires is blocked on `REF-11`; `LVC-13`–`LVC-14` remain unstarted and have
  no comparison tooling yet.

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
- [ ] `SW-03` — add optional compressed molecule-assignment evidence containing
  CB, UMI, ECs, distinct genes, ambiguity class, method, weights, and exclusion
  reason without changing default matrix mass.
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

- [~] `SW-06` — `doctor`, `validate-run`, fingerprints, resume/overwrite safety,
  and atomic manifests exist; finish whole-workflow staging and an atomic
  completion marker.
- [~] `SW-07` — exact-fragment STAR filtering and mate synchronization exist;
  emit a reason for every retained/removed fragment plus lost-truth and
  host-virus-ambiguous boundary counts.
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
- [ ] `SW-09` — split the oversized CLI into thin parsers plus importable service
  functions; convert Snakemake scripts to minimal wrappers without changing
  outputs.
- [~] `SW-10` — run one tiny paired-end fixture through documented CLI commands:
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

- [ ] `REL-01` — move full-workflow-only dependencies such as Snakemake out of
  mandatory pip runtime dependencies; make `doctor --profile pip` check only
  Python/API/reporting/validation capabilities.
- [ ] `REL-02` — define Python 3.11 on `linux-64` as the canonical full-workflow
  toolchain; advertise other platforms only after the same workflow passes.
- [ ] `REL-03` — generate and commit reproducible runtime/development lockfiles
  and an exact external-tool version manifest.
- [~] `REL-04` — build wheel and sdist once, run `twine check`, inspect packaged
  assets, install each in a clean environment, and write `SHA256SUMS`.
- [ ] `REL-05` — run `doctor --profile pip` and `validate-run` against a packaged
  v3 fixture from both clean installations.

### WP2B — Containers and parity

- [~] `REL-06` — Miniforge is digest-pinned; rebuild Docker from the committed
  lock and the exact tested wheel rather than from a mutable source install.
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
- [ ] `REL-12` — pin GitHub Actions by commit SHA and generate provenance/
  attestations for wheel, sdist, OCI, SIF, locks, and checksums.
- [ ] `REL-13` — make release publication depend on green CI for the exact tagged
  SHA; build distributions once and make all downstream jobs consume them.
- [ ] `REL-14` — add protected TestPyPI/PyPI environments, stable approval, GitHub
  release assets, and separate RC/stable container tag behavior.
- [!] `REL-15` — maintainers must provide final author, licence-holder,
  maintainer/contact, ORCID, affiliation, support-window, and trusted-publisher
  metadata.

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

- [~] `SCI-04` — enumerate every ViralScan, STARsolo, traditional alignment,
  Venus, Viral-Track, and VIRTUS row with exact environment/reference
  requirements. Twelve workflows over 78 rows. Both fairness defects are fixed:
  the only exact-truth dataset now reaches all three dedicated comparators, and
  each runs a native-published and a matched-accession-index arm so a tool
  difference is separable from a reference difference. Breadth and tuning
  asymmetries are disclosed rather than removed. Sections remain **`pending`**
  until `SCI-05` round 2 passes; `tool_environments` still blocks on `REL-03`.

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

  Standing limitation, established by the 2026-07-27 code review: every one of
  these rounds read `protocol.yaml`. None could have caught a frozen value that
  no call site consumes, because that is invisible from the protocol side and
  from the output side alike. Protocol review does not substitute for verifying
  that the code honours the contract.

`G3` passes when the protocol validates against its schema, has no unresolved
review finding, is hashed, and outcome-generating jobs have not preceded its
freeze commit.

## WP4 — Freeze production references

Objective: make reference contents reproducible and calibrate host-homology
safeguards without holdout leakage. Estimated effort: 3-5 engineering days plus
about 8 cluster hours per full GRCh38 build.

### WP4A — Profiles and provenance

- [x] `REF-01` — profile names and expanded anellovirus behavior exist; the
  expanded panel is now the **default** rather than opt-in, and its gene IDs
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
- [ ] `REF-03` — apply identical masking, duplicate-ID/sequence validation, and
  manifest generation to dedicated and combined build paths.
- [ ] `REF-04` — make frozen inputs rebuild byte-identical panel FASTA/GTF/t2g
  contents and save a reproducibility audit.
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

### WP4B — Full-genome competition and anellovirus

- [~] `REF-06` — the public `--genome-dlist` path and raw annotations exist;
  build the production curated human-plus-virus index with a checksum-pinned
  GRCh38 D-list.
- [ ] `REF-07` — save the reference manifest, host-homology/low-complexity table,
  index/t2g/GTF, commands, versions, checksums, and build resource accounting.
- [ ] `REF-08` — calibrate homology/complexity exclusion thresholds using only
  preregistered training controls and freeze `thresholds.json`.
- [ ] `REF-09` — prove planted human-homology reads cannot reach probable/strong
  evidence on holdout; retain raw measurements and all excluded calls.
- [!] `REF-10` — an orthogonally confirmed anellovirus-positive sample is needed
  for real sensitivity claims; without it, ship screening support only and label
  every result accordingly.
- [~] `REF-13` — **new, 2026-09-26.** Detection-side reference visibility and
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
    row exists for any of them); the 6 genes with neither CDS nor exon; the
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

## WP4D — Gene-programme inference, layer 2 (new 2026-09-26)

Objective: for viruses layer 1 detected, distinguish latent from productive
expression per cell. Added after measuring that a per-gene comparison cannot do
this at all on this data.

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
- [~] `PROG-07` — **under re-verification (2026-09-27); do not cite these numbers.**
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
- [~] `PROG-08` — **open.** Four of nine viruses have a genuine latency and
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
  35, 39, 45, 51, 52, 56, 58, 59, 66, 69) plus HPV1 and HPV2, which are carried
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
- [~] `HPV-10` — **the scientific limit, recorded so it is not over-read later.**
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
- [ ] `ANELLO-14` — the integration test
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

- [ ] `ANDET-01` — `accession_breadth` is always 1.0: it is computed over
  `found_genes`, which are already detected (`detection.py:166`, `:501-506`).
  Compute it over every index gene of the virus and add per-accession
  genome-coverage breadth, F-005's deciding gate (≤3.41 %).
- [ ] `ANDET-02` — read `host_homology_annotations.tsv` (written at
  `build_reference.py:768`, read by nothing) in detection; demote calls
  concentrated in host-homologous regions; surface `eve_risk` in the report.
- [ ] `ANDET-03` — `claim_scope` column (`screening_only` for Anelloviridae) in
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
- [ ] `ANDET-05` — one genus name per genome: bundled `TTVgp1` IDs resolve to
  "Torque teno virus" while genome-scoped `NC_002076.2_TTVgp1` resolves to
  "Alphatorquevirus". Also fix the `UUKU` and `VARV` aliases.
- [ ] `ANDET-06` — correct three docs: `REF-01`'s "now the default" (the CLI
  default is off), the STAR-defaults claim at `host_filter.py:176-192`, and the
  genome-reference ranking in `docs/cli_reference.md:200`.
- [ ] `ANDET-07` — make the expanded panel the default quantification reference,
  **only after `ANDET-01`–`ANDET-04` land**. `protocol.yaml`'s
  `anellovirus_expanded` entry says "never the default" and sits under the frozen
  `frozen_inputs` digest, so this needs a `DEV-0xx` deviation record.
- [ ] `ANDET-08` — genus-level LOD95: plant held-out genomes from all 8 genera in
  10x v3 geometry into a checksum-pinned healthy-PBMC background (`VAL-04`),
  probit fit per the `SCI-03` method; the unplanted background must yield no
  reported call. A minimal `VAL-01` slice; closes `SENS-06` for this family.

## WP4H — Comprehensive human-virus catalogue (new 2026-09-27)

Objective: a generated, frozen human-host viral catalogue with real gene
structure and deterministic names. No comprehensive index exists: the widest
build (`_misc/viralscan_panel_ref_genomic`) holds ~2,216 genomes, ~88 % of them
anelloviruses, so ~99 other species.

- [ ] `CAT-01` — `build-ref` discards the real GTF (`build_reference.py:609`) and
  turns every non-anellovirus accession into one `{acc}_gene1` gene (`:704`,
  `:930-947`), which would silently disable gene programmes on any natively built
  index. Use the `_genbank_to_gtf` output; emit `exon` rows for CDS-only
  features; placeholders only for CDS-less records.
- [ ] `CAT-02` — `extras/build_virus_catalog.py` → `src/viralscan/data/virus_catalog.tsv`,
  cache-first via `ncbi_fetch.fetch_genbank()`: NCBI Virus RefSeq complete
  genomes with human host ∪ bundled panel ∪ 2,042 anelloviruses ∪ 16 HPV
  genotypes ∪ SARS-CoV-2. Becomes the frozen `broad-discovery` list (`REF-13`).
- [ ] `CAT-03` — names from the catalogue (genome-scoped prefix → species/genus);
  segmented viruses grouped; target 0 % unnamed gene IDs (13.4 % today on the
  covid index).
- [ ] `CAT-04` — risk classes: exclude human endogenous retroviruses; flag
  integrated ciHHV-6, `EVE_RISK_GENERA`, and vector/reagent contaminants.
- [ ] `CAT-05` — duplicate guard (`validate_reference_records`) on every build
  path, including `scripts/build_bundled_panel_ref.py` (the `NC_002076.2`
  duplicate broke `kallisto index` once already).
- [ ] `CAT-06` — the index manifest carries its own viral GTF and catalogue and
  `analysis.py` reads them, so `-gtf` no longer silently drops the panel and
  nothing depends on the unregistered Zenodo DOI (`REF-11`).
- [ ] `CAT-07` — diversity-aware representatives for high-diversity families,
  chosen by leave-one-out k-mer capture (the `REF-01` method), not one exemplar.
- [ ] `CAT-08` — build with native `viralscan build-ref` and measure leave-one-out
  confusability per family: where reads from non-indexed isolates land. Also the
  valid `HPV-11` test.

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
- [~] `PROG-11` — catalogue biology, each change verified against primary
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
- [ ] `PROG-12` — states `latent` / `reactivating` (immediate-early only) /
  `productive` / `mixed` / `indeterminate`; symmetric breadth thresholds (latent
  needs 1 group today, productive 2); a per-marker UMI floor. (Not "honour
  `non_overlapping`": `_breadth` ignores it deliberately — the overlap group is
  the unit.) HCMV needs its own handling here: because latency mirrors
  low-level late-lytic expression, a presence-based `productive` call is not
  specific either, so it needs a per-cell quantity threshold or
  `not_applicable`.
- [ ] `PROG-13` — merge 31-mer-identical repeat copies (HSV LAT/ICP0/ICP4 in
  TRL/IRL, VZV ORF62/ORF63 in TRS/IRS) to one gene_id in t2g, so they reach the
  unique layer.
- [ ] `PROG-14` — unresolved markers fail loudly instead of a log line
  (`scripts/gene_programs.py:180-187`); overlap groups from the active index's GTF.
- [ ] `PROG-15` — measure before adding antisense latency transcripts (VLT, LAT
  intron, LUNA): strandedness (`kallisto bus` runs with no strand flag), the share
  of EBV LCL reads outside annotated exons, and which markers are real index
  targets (the KSHV GTF has 26 exon rows for 96 genes).
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
- [ ] `TONSIL-02` — native 5′ support so `viralscan` itself can run this library:
  - a `--strand` option passed to `kb count` (overlaps `PROG-15`);
  - the cellranger cell barcodes as `-w`;
  - `-x 0,0,16:0,16,28:1,0,0`.

  Measure on a 1 M-read subsample first: host mapping rate under forward,
  reverse and unstranded, then `viralscan check-whitelist`. The native counts
  must agree with `TONSIL-01`.

## WP5 — Build and validate the truth panel

Objective: create deterministic read/molecule truth across supported chemistry
and ambiguity regimes. Estimated effort: 1-2 engineering weeks plus compute.

### WP5A — Generator

- [ ] `VAL-01` — implement a seeded generator spanning viral abundance, infected-
  cell fraction, host homology, sibling viruses, low complexity, PCR duplication,
  CB/UMI collisions, ambient/index hopping, and 10x v2/v3/Drop-seq geometry.
  - **2026-09-26:** `VAL-01` also unblocks `SENS-06` (WP4C). A generator that
    plants a target at a known abundance is exactly the positive control
    `SENS-04` consumes, so the two should be built together: it turns
    `informative_negative` from always-false into a measured quantity, and
    supplies the `E8`/`D17` probit input at the same time.
- [ ] `VAL-02` — emit paired FASTQs, `truth_manifest.tsv`, read/molecule truth
  tables, barcode/chemistry metadata, input hashes, and a run manifest.
- [ ] `VAL-03` — add synthetic host-only and adversarial GRCh38-homology
  conditions plus reagent/empty-droplet controls when available.
- [ ] `VAL-04` — add a checksum-pinned 10x healthy-donor PBMC v3 presumed-negative
  observational control; never label it absolute ground truth.
- [ ] `VAL-05` — add planted target and sibling-virus positives for EBV, HHV-6,
  HSV-1/2, KSHV, and a separately scored anellovirus panel.

### WP5B — Scorer and tiny gate

- [ ] `VAL-06` — implement molecule/cell precision, recall, F1, AUPRC, sibling
  confusion, host-homology false positives, burden concordance, calibration, and
  limit-of-detection scoring.
- [ ] `VAL-07` — validate the scorer against hand-computed fixtures and reject
  denominator, feature, barcode, or count-layer mismatch.
- [ ] `VAL-08` — run the golden tiny panel end to end and prove exact planted-
  molecule recovery, count conservation, determinism, row-order invariance, and
  chunk-size invariance.
- [ ] `VAL-09` — freeze generator/scorer version, seeds, manifests, and expected
  tiny outputs before cluster-scale execution.
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

Primary comparison rules: identical FASTQ hashes and viral sequences, same host
release, outcome-independent cell anchor, audited feature intersection, and
unique-only molecule parity. Recommended ambiguity-aware outputs are a separate
secondary analysis. Run at least EBV, HHV-6B, and HSV-1 for both dedicated
comparators.

`G5b` passes when every preregistered row is complete or transparently failed,
raw outputs and commands are retained, and an independent audit finds no
outcome-selected barcodes or mismatched denominator, annotation, feature, or
count layer.

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
  FAQ, reference-panel, support, security, and migration docs with v3 contracts.
- [~] `DOC-02` — generate CLI/default tables from the parser and test exact
  defaults; remove plain `em`, removed primary-call modes, silent knee fallback,
  v2.5 container commands, and raw-UMI language for fractional estimates.
- [~] `DOC-03` — clearly label every output as observation, model estimate,
  evidence tier, diagnostic flag, or biological interpretation.
- [~] `DOC-04` — document combined versus two-step information loss, anellovirus
  screening limits, pip/full-workflow tiers, and legacy rebuild-only migration.
- [ ] `DOC-05` — execute clean-install quickstart commands and validate the
  resulting run using only documented steps.

### WP8B — Vignettes and claims

- [~] `DOC-06` — eight notebooks exist but contain pre-v3 calls/values; rebuild
  them with negative and ambiguous examples and only v3 APIs/artifacts.
- [ ] `DOC-07` — execute six lightweight notebooks in CI and the reference/full-
  workflow notebooks in the locked scheduled workflow; save logs and hashes.
- [ ] `DOC-08` — add a balanced five-workflow pros/cons table generated from the
  harmonized benchmark rather than rhetorical claims.
- [ ] `DOC-09` — add a claim-registry schema, validator, stale-hash detection, and
  coverage check for README, docs, manuscript, tables, figures, and captions.
- [ ] `DOC-10` — require each quantitative, comparative, validated, performance,
  specificity, and installation claim to resolve to a `validated_v3` artifact.

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
  security-response policies with named contact routes.
- [ ] `OPS-02` — establish reference-update cadence and manifest compatibility
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
| `G1` software | `[~]` | full unit/property/safety/tiny-workflow suite |
| `G2` distribution | `[~]` | clean installs, locks, Docker/Apptainer parity, supply-chain reports |
| `G3` preregistration | `[ ]` | reviewed, schema-valid, hashed protocol frozen before outcomes |
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

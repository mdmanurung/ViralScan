# ViralScan 3.0 stable/G8 execution checklist

Plan: `docs/plans/2026-08-08-viralscan-v3-stable-g8.md`

Mark an item complete only after its named evidence passes. Update `PLAN.md` in
the same commit as each completed tracker item. Execute strictly in gate order.

## Start and baseline

- [x] Read the frozen plan, this checklist, the executor handoff, `PLAN.md`, and
  `CONTEXT.md` completely.
- [x] Verify branch `codex/viralscan-v3` and starting HEAD
  `26260e7cbc18cc0e7777379e1db12cc594778dd4`.
- [x] Verify the three persistent dirty-file SHA-256 values from the frozen plan.
- [x] Record whether the untracked `.mycelium` lock still exists without changing
  or staging it.
- [x] Confirm the initial runtime reports Python 3.11 from
  `benchmark_runs/legacy_v2_v3/env_full/bin/python`.

## G0 — governance

- [x] Add failing tests for a positive public ship-scope allowlist.
- [x] Add `config/public_ship_scope.json` with explicit wheel, sdist, Docker-
  context, public-doc, and claim-bearing memberships.
- [x] Replace recursive public-doc packaging with explicit allowlisted entries.
- [x] Make the Docker build consume an allowlisted staged context or the tested
  wheel instead of `COPY .`.
- [x] Extend data-governance tests to reject institutional paths in every
  allowlisted text file.
- [x] Add failing schema-boundary tests for incomplete claim records.
- [x] Add `schemas/v3/claim_registry.schema.json`.
- [x] Implement `scripts/validate_claim_registry.py` with fail-closed schema,
  artifact-existence, SHA-256, status, and `--coverage` checks.
- [x] Add tests rejecting unregistered public claim markers and legacy numerical
  evidence presented as v3.
- [x] Expand `claims/registry.json` with all required provenance fields.
- [x] Add failing tests for incomplete, duplicate, missing, and stale inventory
  rows.
- [x] Add `schemas/v3/artifact_inventory.schema.json`.
- [x] Implement `scripts/validate_artifact_inventory.py`.
- [x] Create `analysis/v3_artifact_inventory.tsv` without institutional absolute
  paths.
- [x] Implement `scripts/check_ship_scope.py` for wheel, sdist, and effective
  Docker-context member lists.
- [x] Build wheel and sdist once into a temporary directory and pass the member-
  list check.
- [x] Run focused governance tests and all three governance validators.
- [x] Run the draft protocol validator without opening outcome data.
- [x] Append dated reconciliation evidence to `PLAN.md`; update its next action
  and `G0` state without rewriting historical log entries.
- [x] Recheck the frozen dirty paths and inspect `git diff --check`.
- [ ] Commit the verified `G0` slice atomically if every exit criterion passes.

## G1 — software contracts

- [ ] Add the optional compressed molecule-assignment evidence schema and tests.
- [ ] Implement assignment evidence without changing default matrix mass.
- [ ] Add tests proving missing material outputs invalidate the correct Snakemake
  rules.
- [ ] Declare all material workflow artifacts and validated completion manifests
  in the DAG.
- [ ] Complete atomic whole-workflow staging, resume, overwrite, and corruption-
  recovery behavior.
- [ ] Add service-boundary tests before splitting CLI parsers from services.
- [ ] Split CLI/service code without changing quantification syntax or outputs.
- [ ] Add the versioned per-fragment host-decision schema.
- [ ] Add adversarial cDNA-mate host-filter classification unit tests.
- [ ] Pin STAR version and every alignment/filter parameter.
- [ ] Implement one decision row per input fragment and exact partition checks.
- [ ] Enforce retained-mate synchronization and fail closed on malformed/truncated
  FASTQ.
- [ ] Add real-STAR fixtures for every frozen alignment class.
- [ ] Run focused unit/property/mutation/rerun/corruption/integration tests.
- [ ] Run the documented tiny workflow through evidence, BAM, BLAST, plots, IGV,
  and `validate-run`.
- [ ] Compare representative downstream artifacts with the pre-change baseline.
- [ ] Append `G1` evidence to `PLAN.md` and commit only if the full gate passes.

## G2 — distribution

- [ ] Add tests for minimal `pip`, canonical `workflow`, and deprecated `full`
  doctor profiles.
- [ ] Move full-workflow-only packages out of mandatory pip dependencies.
- [ ] Implement and document `doctor --profile pip|workflow` plus the `full` alias.
- [ ] Generate canonical Python 3.11 `linux-64` runtime/development locks and an
  external-tool manifest.
- [ ] Build wheel and sdist once and inspect their exact contents.
- [ ] Metadata-check and install wheel and sdist in clean environments with real
  dependency resolution.
- [ ] Build OCI from the tested wheel and lock when an approved local runtime is
  available.
- [ ] Build SIF from the OCI digest when an approved local runtime is available.
- [ ] Run and compare the identical tiny workflow in conda, OCI, and SIF.
- [ ] Generate parity, audit, scan, SBOM, licence, checksum, and attestation
  artifacts.
- [ ] Pin actions by SHA and harden exact-SHA build-once RC/stable release paths.
- [ ] Prove a no-publish dry run and verify RC never receives `latest`.
- [ ] Append `G2` evidence to `PLAN.md` and commit only if the full gate passes.

## G3 — protocol freeze

- [ ] Add a duplicate-key-rejecting YAML loader and regression test.
- [ ] Replace stale protocol blocker/readiness text.
- [ ] Derive one canonical workflow-row count from the manifest rows.
- [ ] Add endpoint-to-row applicability schema and coverage tests.
- [ ] Add mixed host-virus truth to every required paired workflow arm.
- [ ] Freeze the continuous score and nested probable/strong rules.
- [ ] Remove nominal-confidence metrics for non-probabilistic scores.
- [ ] Freeze thresholds, ties, missing/not-estimable behavior, factor levels,
  seeds, partitions, BCa bootstrap, and failure states.
- [ ] Invoke the execution guard in the manifest builder.
- [ ] Invoke the execution guard immediately before each local/scheduled row.
- [ ] Implement an atomic immutable single-use holdout receipt.
- [ ] Test pre-freeze training rejection, pre-threshold holdout rejection, and
  second/concurrent holdout rejection.
- [ ] Obtain and record independent passing review for the exact protocol SHA.
- [ ] Write the freeze sidecar and append `G3` evidence to `PLAN.md`.

## G4 — references and mechanics freeze

- [ ] Freeze source releases, accessions, transformations, D-list, homology
  annotations, licences, duplicates, and comparator linkages.
- [ ] Build references twice in independent roots and compare content hashes.
- [ ] Implement the hand-calculated outcome-ineligible generator/scorer
  microfixture.
- [ ] Freeze generator/scorer versions, factors, seeds, schemas, budgets, and tiny
  expected outputs.
- [ ] Extend the freeze sidecar with every `G4` hash.
- [ ] Append `G4` evidence to `PLAN.md` without generating training/holdout
  outcomes.

## G5 — scientific execution

- [ ] Obtain explicit approval for each needed network, scheduler, restricted-
  data, credential, and outcome-execution action.
- [ ] Generate independent training and holdout panels after both mechanical gates
  pass.
- [ ] Run training through the guard and freeze thresholds/training manifests.
- [ ] Run holdout once through the guard and consume the immutable receipt.
- [ ] Execute every frozen workflow row and retain every attempt state.
- [ ] Score exact truth and contextual public evidence under their distinct
  estimands.
- [ ] Compute biological-sample BCa uncertainty and descriptive resource metrics.
- [ ] Freeze the complete result bundle and independent audit.
- [ ] Append `G5` evidence to `PLAN.md` only when every exit criterion passes.

## G6a and G6 — docs and manuscript

- [ ] Regenerate parser-derived CLI/default documentation and enforce parity.
- [ ] Generate every public quantitative claim and application artifact from the
  frozen G5 bundle.
- [ ] Rebuild notebooks as parameterized read-only consumers and execute their
  proper tiers.
- [ ] Pass Sphinx warnings/link checks, shell examples, quickstart, notebooks, and
  claim coverage.
- [ ] Obtain independent human editorial review and append `G6a` evidence.
- [ ] Generate every manuscript number, table, figure, caption, and supplement
  from the frozen bundle.
- [ ] Rerun host response with v3 labels plus depth/mitochondrial controls and
  include it only on frozen reproduction.
- [ ] Obtain approved author/ethics/funding/conflict/data-access metadata.
- [ ] Pass independent methods/statistics, citation, and claim-artifact audits.
- [ ] Append `G6` evidence only when every manuscript claim resolves to a digest.

## G7 and G8 — publication

- [ ] Obtain explicit RC publication approval before creating or pushing a tag.
- [ ] Publish `3.0.0rc1` artifacts without `latest` from the exact green SHA.
- [ ] Collect three complete independent laboratory tester manifests.
- [ ] Close release blockers and rerun every affected gate.
- [ ] Prove RC fixes preserve frozen scientific outputs or record an amendment.
- [ ] Obtain software and benchmark DOI reservations and stable publication
  approval.
- [ ] Publish `v3.0.0` artifacts from the exact reviewed SHA.
- [ ] Update, build, submit, and verify Bioconda with the real PyPI sdist hash.
- [ ] Verify durable archives, metadata, DOI/public links, and public-format count
  parity.
- [ ] Submit the manuscript and record confirmation before marking `G8` complete.

## Initial G0 validation commands

- [x] Run
  `benchmark_runs/legacy_v2_v3/env_full/bin/python -m pytest tests/test_data_governance.py tests/test_claim_registry.py tests/test_artifact_inventory.py tests/test_ship_scope.py -q`.
- [x] Run
  `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/check_data_governance.py`.
- [x] Run
  `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/validate_claim_registry.py --coverage`.
- [x] Run
  `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/validate_artifact_inventory.py`.
- [x] Run
  `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/validate_v3_protocol.py`.
- [ ] Build distributions with
  `benchmark_runs/legacy_v2_v3/env_full/bin/python -m build --no-isolation --outdir <temporary-directory>`.
- [x] Run `scripts/check_ship_scope.py` against the built wheel, sdist, and
  effective Docker context.
- [x] Run `git diff --check`, `git status --short --untracked-files=all`, and the
  frozen dirty-path SHA-256 checks.
- [ ] Run the full non-network repository gate recorded in the reconciled
  `PLAN.md` before closing `G0`.

`twine` is not installed in the frozen initial environment. Do not download it
without approval. Its metadata check belongs to `G2`; report the missing tool if
reached before an approved locked environment supplies it.

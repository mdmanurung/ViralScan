# ViralScan 3.0 stable/G8 frozen execution plan

Status: **frozen for execution**

Frozen: 2026-08-08

Repository: `/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan`

This plan is the execution contract derived from the user-supplied ViralScan 3.0
Stable/G8 Completion Plan and a repository-grounded requirements, risk, and test
review. `PLAN.md` remains the append-only operational ledger. Historical evidence
in that file must not be rewritten.

## Baseline that must be preserved

- Branch: `codex/viralscan-v3`
- HEAD: `26260e7cbc18cc0e7777379e1db12cc594778dd4`
- Upstream relation at freeze: 42 commits ahead of `origin/codex/viralscan-v3`
- Pre-existing dirty paths, excluded from every edit and commit:
  - `.living/INDEX.md` — SHA-256
    `0e8a386ed2b18f8601a3847520f379b455f02e12109e6c752d48880f92d99b34`
  - `docs/review-clear-execute-plan.md` — SHA-256
    `27e8832e1eb7ff118d771583cb260e6a5a076e7fc386fde4d8333876ae5064dd`
  - `docs/review-clear-execute-tasks.md` — SHA-256
    `5e9b80328fa0cecfee1b15639acfa24204daf5ea39a69c5fd1e0e7f8ce80cad7`
  - `.mycelium/locks/INDEX.md-2828a51c35284fa97c79.lock` — SHA-256
    `00fee6eaaab6d660431a706eda6373c0c03e0d147897290ce96863b9afb60577`

The lock file may be removed by its owning tool. Do not recreate, stage, or alter
it; if it disappears, record that fact and continue protecting the other three
paths.

## Gate law

1. Execute gates in order. Existing later-gate code may be tested or inventoried,
   but it is not promoted while an earlier gate is open.
2. No truth-panel, training, public-positive, comparator, or holdout outcome may be
   generated, scheduled, opened, summarized, or inferred until the mechanical
   `G3` and `G4` gates pass for the exact frozen hashes.
3. A local technical pass is not scientific validity or promotion. Human review
   and owner approval remain separate evidence.
4. Failed, unfavorable, timeout, OOM, invalid, unsupported, incomparable, and not
   applicable attempts remain explicit records. Missing values are never zero.
5. Each completed tracker item updates `PLAN.md` in the same atomic commit. Dated
   corrections are appended; historical evidence is not rewritten.
6. Heavy computation belongs in Python/R scripts and is orchestrated by
   Snakemake. Notebooks are parameterized, read-only consumers of immutable
   outputs.
7. Use the locked Python 3.11 runtime at
   `benchmark_runs/legacy_v2_v3/env_full/bin/python` for the initial local gate.
   Bare `python3` is Python 3.6.8 and is not a valid project runtime.

## Reconciled design decisions

The following repository contradictions are resolved by this plan:

- The stale `PLAN.md` next-action text does not reopen `SW-02`, `SW-04`, or
  `SW-05`; their current checked rows and commits are retained. `G0` is the next
  gate.
- `G4` freezes reference identity, transformations, D-lists, homology annotations,
  licences, duplicate checks, and byte-reproducible builds. Any threshold or tier
  calibrated from training outcomes belongs to `G5-training`; any exact-negative
  proof using holdout outcomes belongs to `G5-holdout`. This removes the current
  `G4`/`G5` dependency cycle.
- `doctor --profile workflow` is the canonical full-workflow spelling. The current
  `full` spelling remains a tested deprecated compatibility alias through 3.x.
- The workflow matrix has one canonical, schema-validated row count derived from
  its rows. Prose counts are never authoritative. Endpoint-to-row applicability is
  checked mechanically, including mixed host-virus truth populations.
- Candidate evidence remains any nonzero selected-method molecule support and is
  never called validated infection. Probable and strong tiers are separate,
  nested, deterministic training-only rules. Non-probabilistic scores do not
  receive nominal confidence or calibration-error claims.

## G0 — Governance and legacy eligibility

Start here. No other gate may be promoted until every exit check below passes.

1. Reconcile `PLAN.md` with current code, tests, commits, artifacts, protocol
   state, and the frozen dirty baseline. Append dated corrections for stale
   next-action and status text.
2. Define `config/public_ship_scope.json` as the single positive allowlist for
   public wheel, sdist, OCI build context, documentation, and claim-bearing files.
   Exclude `.living`, `.mycelium`, operational plans and handoffs, diagnostic and
   historical benchmark archives, private paths, scratch outputs, scheduler
   packets, and pre-v3 numerical evidence through packaging rules, not deletion.
3. Replace recursive documentation packaging with explicit public documentation
   inclusions. Make the Docker build consume a staged allowlisted product context
   or an exact tested wheel; it must not copy the repository wholesale.
4. Add `schemas/v3/claim_registry.schema.json` and
   `schemas/v3/artifact_inventory.schema.json`. Add fail-closed validators at
   `scripts/validate_claim_registry.py`,
   `scripts/validate_artifact_inventory.py`, and `scripts/check_ship_scope.py`.
5. Expand `claims/registry.json`. Every public claim records status, source
   location, artifact identity and SHA-256, Git SHA, input/reference hashes,
   schema, layer, denominator, generation command, validation scope, and any
   legacy/rejection reason required by status.
6. Create `analysis/v3_artifact_inventory.tsv`. Inventory each claim-bearing
   input, reference, intermediate, result, command, environment, scheduler record,
   and failure record with stable identity, public/private scope, software and
   counting version, Git SHA, input/reference hashes, artifact SHA-256, schema,
   layer, denominator, generation command, and rebuild eligibility. Never place
   an institutional absolute path in a tracked row.
7. Reject missing inventory artifacts, stale hashes, duplicate identities,
   incomplete claim metadata, unregistered claim markers in allowlisted public
   prose, and pre-v3 numerical evidence presented as v3.
8. Build wheel and sdist into a temporary directory and compare their member lists
   with the allowlist. Compute and validate the effective Docker context without
   building or publishing an image.

`G0` closes only when governance, inventory, claim coverage, legacy-ineligibility,
and package/context membership checks pass. The incomplete outcome-ineligible
legacy diagnostic may remain recorded as partial if it is not used to support a
public v3 claim.

## G1 — Core software contracts

1. Complete optional compressed molecule-assignment evidence without changing
   default matrix mass or the five selected multimapping methods.
2. Finish whole-workflow atomic staging, material completion manifests, resume and
   overwrite safety, and corruption recovery. Snakemake declares every material
   artifact, not only `.done` sentinels; missing or corrupt outputs invalidate the
   producing rule.
3. Split the oversized CLI into thin parsers and importable services without
   changing established quantification CLI syntax or v3 output semantics.
4. Rebuild exact-fragment host filtering around the cDNA mate. Pin STAR version and
   every alignment/filter parameter. Remove a pair for any accepted primary,
   secondary, supplementary, chimeric, multimapping, or excessive-locus human
   alignment; retain only an explicit no-accepted-host-alignment decision.
5. Emit one compressed decision row per input fragment with canonical read ID,
   CB, UMI, decision, reason, alignment class/count, MAPQ/flags where applicable,
   and retained-output location. Require exact input = retained + removed equality
   and synchronized retained mates.
6. Add real-STAR fixtures for unmapped, unique host, host multimapping, accepted
   plus rejected alignments, secondary, supplementary, chimeric, excessive-locus,
   malformed-pair, and truncated-FASTQ cases.
7. Run unit, property, mutation, schema-boundary, order/chunk invariance, rerun
   invalidation, corruption recovery, integration, and the full documented tiny
   end-to-end path including evidence/BAM/BLAST/plots/IGV and `validate-run`.
8. Compare representative downstream artifacts before and after changes. Require
   exact equality when semantics are unchanged and a documented migration when
   they change.

`G1` closes only with finite, non-negative, exact molecule-partition and matrix-
mass invariants for every selected method and no known correctness or data-loss
defect. Unsafe CB-UMI-wide kallisto subtraction remains unavailable.

## G2 — Install, container, and release engineering

1. Define a minimal pip/API profile for Python 3.9-3.12 and a canonical locked
   Python 3.11 `linux-64` workflow profile. Remove Snakemake and native full-
   workflow dependencies from mandatory pip requirements.
2. Implement `doctor --profile pip|workflow`; retain `full` as a deprecated alias
   for `workflow`. The pip profile checks only installed parsing, validation,
   reporting, and documented Python APIs.
3. Generate reproducible runtime/development locks and an exact external-tool
   manifest.
4. Build wheel and sdist once; inspect, metadata-check, and install each in clean
   environments with real dependency resolution and without `--no-deps`.
5. Build OCI from the tested wheel and lock. Build SIF from the resulting OCI
   digest, not from a second dependency solve.
6. Run the identical tiny workflow through locked conda, OCI, and SIF. Compare
   validated molecule counts and deterministic hashes; document allowed
   nondeterminism.
7. Audit the locked product environment, scan OCI, and generate SBOMs, licences,
   checksums, and attestations.
8. Pin GitHub Actions by commit SHA; make publication depend on green CI for the
   exact tagged SHA; fan out one build; separate RC/stable paths; protect PyPI
   environments; prohibit `latest` for RC; prove a no-publish dry run.

`G2` closes only when clean installs and all three workflow formats reproduce the
validated tiny result. Building containers locally is allowed only when the
runtime is already available and causes no network or registry write; publication
always requires approval.

## G3 — Freeze the scientific protocol

Outcome access remains disabled throughout this gate.

1. Replace stale blocker text and derive one canonical workflow-row count.
2. Add schema-validated endpoint-to-row applicability. Every endpoint population
   has applicable rows; paired arms cover identical datasets; mixed host-virus
   truth is present wherever required.
3. Freeze the continuous score formula, direction, aggregation, missing values,
   ties, factor levels, seeds, partitions, threshold grids, biological-sample BCa
   bootstrap, small-sample/not-estimable behavior, and failure states.
4. Remove nominal-confidence metrics unless a genuinely probabilistic score is
   preregistered. Freeze nested probable and strong training-only rules.
5. Reject duplicate YAML keys. Add tests for the existing duplicate-key class.
6. Implement the execution guard in the validation-manifest builder and
   immediately before every local or scheduled row.
7. Implement an atomically created, immutable single-use holdout receipt binding
   protocol, code, reference, environment, partition, workflow-manifest,
   training-result, and threshold hashes. Reject a second or concurrent receipt.
8. Prove training rejection before freeze, holdout rejection before threshold
   hashes, and second-holdout rejection.
9. Obtain independent protocol reviews of the current SHA until no Major remains;
   store a freeze sidecar for the passing SHA.

`G3` closes only when the protocol and guard are schema-valid, buildable, reviewed,
hash-frozen, and still mechanically prohibit outcome access until `G4` passes.

## G4 — Freeze references and generator/scorer mechanics

Outcome access remains disabled throughout this gate.

1. Freeze checksum-pinned human and viral source releases, exact accessions, GTF
   transformations, D-list, outcome-independent homology annotations, licences,
   duplicate checks, and native/matched comparator linkage.
2. Build references twice from independent clean roots and require byte-equivalent
   content hashes.
3. Implement an outcome-ineligible generator/scorer microfixture with hand-
   calculated truth projection and expected metrics. It cannot substitute for the
   preregistered panel.
4. Freeze generator/scorer code versions, numeric factors, seeds, truth schemas,
   compute/storage budget, and tiny expected outputs without generating training
   or holdout outcomes.
5. Record protocol, code, reference, environment, partition, workflow-manifest,
   generator, and scorer hashes in the freeze sidecar.

`G4` closes on immutable reference and mechanics evidence. Training-derived
homology thresholds and tier calibration are explicitly deferred to
`G5-training`; holdout negative proof is deferred to `G5-holdout`.

## G5 — Execute and freeze the scientific comparison

This gate requires explicit approval for network downloads, scheduler submissions,
restricted data access, or credential use. Do not start it merely because local
code is ready.

1. Generate complete training and holdout panels with separate sample/seed
   lineages and no shared template, locus, molecule, UMI, or barcode.
2. Run training only after the training guard passes. Calibrate thresholds and any
   training-derived homology/tier rules only on training; freeze
   `thresholds.json`, the training-result manifest, and hashes.
3. Open holdout exactly once after its guard passes and atomically consume the
   receipt.
4. Execute every frozen row: ViralScan combined; exact-fragment STAR filter to
   ViralScan; STARsolo combined/two-step; traditional subtraction with/without
   CB/UMI preservation; native/matched Venus, Viral-Track, and VIRTUS; and
   native/matched ViralScan 2.2.0 diagnostic arms.
5. Run comparators at published defaults in digest-pinned environments; never
   silently patch algorithms. Retain every planned row and every attempt state.
6. Keep estimands distinct: pure within-tool reference strategy, hybrid product
   workflow, and cross-tool unique-molecule parity. Report ambiguity-aware model
   outputs separately.
7. Score exact synthetic truth for precision, recall, F1, AUPRC, sibling
   confusion, host-homology false calls, and host-filter loss. Treat real public
   positives/negatives as contextual evidence, not ground truth.
8. Use biological-sample BCa bootstrap. Record hardware-specific wall time, peak
   RSS, scratch, and output size descriptively.
9. Freeze one immutable bundle of raw row/attempt manifests, failures, metrics,
   figures, commands, versions, hashes, deviations, and an independent audit.

`G5` closes only with 100% ViralScan invariant compliance, complete ambiguity and
unresolved reporting, no probable/strong calls in required exact negatives, all
rows retained, and no unsupported general-superiority claim.

## G6a — Public documentation and application

1. Restrict editorial rewriting to the public documentation allowlist.
2. Generate CLI/default tables from the parser and fail CI on drift.
3. Keep README concise: purpose, supported use, recommended combined workflow,
   installation, smallest example, first outputs, combined/two-step distinction,
   limitations, and links.
4. Generate all quantitative prose, tables, figures, captions, and claims from the
   frozen G5 bundle; never transcribe values manually.
5. Build the eight-part combined-versus-two-step application specified in the
   source plan, with contrast type and annotation incompatibility explicit.
6. Rebuild notebooks as parameterized read-only views; execute light notebooks in
   CI and heavy notebooks in the locked tier.
7. Pass Sphinx `-W`, link checking, CLI parity, shell examples, clean quickstart,
   notebooks, claim coverage, and independent human editorial review.

## G6 — Manuscript

1. Generate every manuscript number, table, figure, caption, and supplement from
   the frozen bundle.
2. Report the full matrix, failures, incomparable endpoints, uncertainty,
   chemistry confounding, reference differences, and tuning asymmetry.
3. Rerun host response using v3 labels plus depth and mitochondrial controls;
   include it only if the frozen analysis reproduces.
4. Keep anellovirus screening-only without an orthogonally validated positive;
   keep TTV as a homology case without a biological-positive claim.
5. Obtain approved authorship, affiliation, ORCID, CRediT, funding, conflict,
   ethics, consent, and data-access text.
6. Complete independent methods/statistics, citation, and claim-to-artifact audits
   with no unresolved Major.

`G6` closes only when every manuscript claim resolves to a frozen artifact digest.

## G7 — Release candidate

1. Freeze one green protected-main SHA as `3.0.0rc1`.
2. After explicit publication approval, publish versioned prerelease wheel, sdist,
   OCI, SIF, locks, checksums, SBOMs, licences, attestations, reference archive,
   and GitHub prerelease. Do not update `latest`.
3. Obtain three independent laboratory reports spanning supported chemistries;
   each performs clean install, tiny workflow, `validate-run`, and one real
   supported workflow.
4. Close every correctness, data-loss, install, documentation, and reproducibility
   failure; rebuild and rerun every gate after fixes.

## G8 — Stable publication and submission

1. Prove RC fixes did not change frozen scientific outputs except by documented,
   audited amendment.
2. Reserve distinct software and benchmark/archive DOIs.
3. Freeze `v3.0.0` from the exact reviewed green SHA.
4. After explicit stable-publication approval, publish wheel, sdist, versioned and
   `latest` OCI, SIF, locks, references, checksums, SBOMs, licences, attestations,
   release notes, and migration guide.
5. Insert the real PyPI sdist hash into Bioconda; lint, build, submit, and verify
   the tiny workflow.
6. Archive software and benchmark artifacts, verify DOI/public links, update
   citation metadata, and submit the manuscript.

`G8` passes only when public formats reproduce validated counts, archives are
durable, metadata is final, and manuscript submission is confirmed. WP12 is
post-completion maintenance.

## Public interfaces and artifacts

- Preserve quantification CLI spelling and v3 output semantics.
- Add/finalize `doctor --profile pip|workflow`; retain the `full` alias.
- Version fragment-decision, workflow-row/attempt, truth, claim-registry,
  artifact-inventory, and result-bundle schemas.
- Keep `X` as selected-method molecule estimates, with non-overlapping
  `counts_unique` and `counts_ambiguous_allocated`. Cross-tool parity consumes
  only `counts_unique`.
- Publish no quantitative claim without status, artifact hash, Git SHA,
  input/reference hashes, layer, denominator, command, and validation scope.

## Authorization ledger

Allowed now:

- Read repository and existing local artifact metadata.
- Make reversible edits inside this repository, excluding the frozen dirty paths.
- Run local non-outcome tests, validators, builds, archive inspections, and dry
  runs using already-installed tools.
- Create atomic commits for verified tracker items, with the corresponding
  `PLAN.md` update in the same commit.

Explicit approval required:

- Any network download or live service query.
- Any SLURM or other scheduler submission.
- Container registry, PyPI/TestPyPI, GitHub release/tag/push, DOI/archive,
  Bioconda, external-laboratory, email/message, or manuscript-submission action.
- Credential access beyond ordinary already-configured local read-only tooling.
- Training, public-positive, comparator, or holdout outcome execution after the
  mechanical gates pass.

Forbidden in this packet:

- Editing, staging, cleaning, or committing the frozen dirty paths.
- Deleting or overwriting historical/private/scientific artifacts to make a gate
  pass.
- Inspecting scientific outcomes before their exact gate and explicit execution
  approval.
- Tagging or publishing from the current release workflow.
- Claiming `G7` or `G8` without external evidence and owner promotion.

## Stop conditions

Stop and report the exact blocker when protocol state, hashes, permissions, data,
compute, storage, metadata, dependencies, or required approvals are absent. Stop
before any destructive action. Stop if a frozen dirty path changes. Stop if an
implementation decision would alter public quantification behavior beyond this
plan. A missing optional test tool is reported; it is not installed from the
network without approval.

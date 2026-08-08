# No-Parent-History Executor Handoff: ViralScan stable G8

## Objective

Advance ViralScan strictly from the earliest open gate toward stable `G8`,
implementing and verifying all locally authorized work and stopping at the first
genuine scientific, dependency, permission, or external-approval boundary.
`RC-READY` is not completion; never report `G8` without its external evidence.

## Repository

`/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan`

Revised Plan Path:
`/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan/docs/plans/2026-08-08-viralscan-v3-stable-g8.md`

Task List Path:
`/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan/docs/plans/2026-08-08-viralscan-v3-stable-g8-tasks.md`

## Frozen Plan

1. Verify the frozen branch, HEAD, dirty-path hashes, runtime, and current
   `PLAN.md` state. Preserve all pre-existing changes.
2. Close `G0` first: positive ship-scope allowlist, claim-registry and artifact-
   inventory schemas/validators, complete sanitized inventory, packaging/context
   enforcement, coverage tests, and dated tracker reconciliation.
3. After `G0` genuinely passes, close `G1`: material Snakemake outputs and atomic
   staging, assignment evidence, thin CLI/services, full exact-fragment STAR
   decisions, adversarial real-tool fixtures, end-to-end and invariance proof.
4. Close `G2`: honest minimal pip tier, locked workflow tier, `workflow` profile
   with compatible `full` alias, clean wheel/sdist installs, build-once OCI/SIF
   parity, supply-chain evidence, and safe exact-SHA RC/stable automation.
5. With outcomes still disabled, close `G3` and `G4`: repair/freeze the protocol,
   endpoint coverage, tier/score rules, dual execution guard, single-use holdout
   receipt, reviewed freeze sidecar, byte-reproducible references, and hand-
   calculated outcome-ineligible generator/scorer microfixture. Training-derived
   calibration belongs to `G5-training`, not `G4`.
6. Stop and obtain explicit approval before any network, scheduler, restricted-
   data, credential, training, public-positive, comparator, or holdout action.
   After approval, execute training, freeze thresholds, consume holdout once, run
   every frozen row, preserve all attempts, and freeze the audited `G5` bundle.
7. Generate `G6a` documentation and the `G6` manuscript only from the frozen
   bundle; require parser parity, executable docs/notebooks, complete claim
   coverage, approved factual metadata, and independent audits.
8. Stop and obtain explicit RC publication approval. Publish `3.0.0rc1` without
   `latest`, collect three independent laboratory reports, close blockers, and
   rerun gates before `G7` promotion.
9. Stop and obtain explicit stable-publication/submission approval. Publish the
   exact reviewed `v3.0.0`, verify Bioconda and durable archives/DOIs, and record
   manuscript submission before `G8` promotion.

## Authorization

Allowed:

- Read repository files and existing local artifact metadata.
- Make reversible repository-local edits required by the frozen plan, excluding
  protected dirty paths.
- Run local non-outcome tests, validators, builds, archive inspections, and dry
  runs with already-installed tools.
- Create atomic commits for completed tracker items, including the corresponding
  `PLAN.md` update in the same commit.
- Use bounded implementation/review subagents only when work is partitioned
  safely and the governing skill permits it.

Approval required:

- Network downloads or live-service queries.
- SLURM or any scheduler submission.
- Container registry, TestPyPI/PyPI, GitHub tag/push/release, DOI/archive,
  Bioconda, external-laboratory, message/email, or manuscript-submission actions.
- New credential use or access to restricted data.
- Training, public-positive, comparator, or holdout outcome execution, even after
  its mechanical gate becomes green.

Forbidden or out of scope:

- Editing, staging, cleaning, deleting, or committing the pre-existing dirty
  paths named below.
- Generating, opening, summarizing, or inferring preregistered outcomes before the
  exact `G3`/`G4` freeze and explicit execution approval.
- Destructive history/worktree operations, silent comparator patching, outcome-
  driven exclusion, or conversion of missing/incomparable values to zero.
- Tagging or publication from the current unsafe release workflow.
- Claiming owner/human promotion or external completion from local evidence.

## Constraints and Non-Goals

- Start on branch `codex/viralscan-v3` at HEAD
  `26260e7cbc18cc0e7777379e1db12cc594778dd4`.
- Protect these pre-existing paths and initial hashes:
  - `.living/INDEX.md` —
    `0e8a386ed2b18f8601a3847520f379b455f02e12109e6c752d48880f92d99b34`
  - `docs/review-clear-execute-plan.md` —
    `27e8832e1eb7ff118d771583cb260e6a5a076e7fc386fde4d8333876ae5064dd`
  - `docs/review-clear-execute-tasks.md` —
    `5e9b80328fa0cecfee1b15639acfa24204daf5ea39a69c5fd1e0e7f8ce80cad7`
  - `.mycelium/locks/INDEX.md-2828a51c35284fa97c79.lock` —
    `00fee6eaaab6d660431a706eda6373c0c03e0d147897290ce96863b9afb60577`
- The `.mycelium` lock may disappear through its owner. Do not recreate it.
- Use `benchmark_runs/legacy_v2_v3/env_full/bin/python` for the initial gate; it
  is Python 3.11.15. Bare `python3` is not a valid project runtime here.
- Preserve quantification CLI spelling and v3 output semantics. Make `workflow`
  canonical for doctor while retaining `full` as a deprecated 3.x alias.
- Never rewrite historical `PLAN.md` evidence; append dated reconciliation.
- `twine` is absent from the initial environment. Do not install it from the
  network without approval.
- The incomplete legacy v2/v3 diagnostic is not required to close `G0` when it is
  ineligible for all public v3 claims and its partial state remains explicit.

## Validation Commands

- `benchmark_runs/legacy_v2_v3/env_full/bin/python -m pytest tests/test_data_governance.py tests/test_claim_registry.py tests/test_artifact_inventory.py tests/test_ship_scope.py -q`
- `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/check_data_governance.py`
- `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/validate_claim_registry.py --coverage`
- `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/validate_artifact_inventory.py`
- `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/validate_v3_protocol.py`
- `benchmark_runs/legacy_v2_v3/env_full/bin/python -m build --no-isolation --outdir <temporary-directory>`
- `benchmark_runs/legacy_v2_v3/env_full/bin/python scripts/check_ship_scope.py --wheel <wheel> --sdist <sdist> --docker-context .`
- `git diff --check`
- `git status --short --untracked-files=all`
- Run the full non-network test/lint/format/integration gate recorded in the
  reconciled `PLAN.md` before closing each applicable software slice.
- Re-run SHA-256 checks for every persistent protected dirty file before staging
  and after tests.

## Stop Conditions

- Stop if repository state drift affects the plan or a protected dirty file
  changes.
- Stop if implementation would change public quantification semantics beyond the
  frozen plan.
- Stop before any outcome access unless both mechanical gates pass and the user
  explicitly approves that outcome execution.
- Stop before destructive, scheduler, network, credential, deployment,
  publication, external-message, DOI, Bioconda, or manuscript-submission action
  without explicit approval.
- Stop if required permissions, dependencies, data, compute, storage, hashes,
  metadata, or independent human review are missing.
- Stop if a gate cannot be evidenced; do not weaken or relabel it.

## Required Artifacts

- The three frozen packet files remain tracked and internally consistent.
- Updated checklist boxes only for verified work.
- Reconciled append-only `PLAN.md` evidence for each completed item.
- Focused and full validation output, including downstream invariance evidence
  where semantics are unchanged.
- Atomic commits for completed tracker items only.
- Final report naming changed files, test results, exact gate state, preserved
  dirty state, blockers, approvals still required, and the single next action.

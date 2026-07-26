# ViralScan v3 next-session handoff

Date: 2026-07-22

## Pause instruction

The user explicitly paused analysis. Do not launch benchmarks, downloads, cluster
jobs, comparator runs, releases, or manuscript regeneration until a new session
authorizes the relevant task.

## Repository state

- Branch: `codex/viralscan-v3`
- Local and upstream head: `a8f32876e6ec4a5d8cf29b21aaf325ccaa50e64f`
- The 18-commit v3 series is pushed to `origin/codex/viralscan-v3`.
- Verified gate: 720 passed, 21 deselected; Ruff check/format, data governance,
  draft protocol validation, and the two live evidence integrations passed.
- Pre-existing `.living/**` edits remain deliberately uncommitted. Preserve them.
- Large EBV BUS/text intermediates, scheduler logs, and the local pre-v3 tracker
  archive are intentionally ignored and must not be committed.

The authoritative task state is in [`PLAN.md`](../../PLAN.md). Locked product and
scientific contracts are in
[`docs/plans/2026-07-22-viralscan-v3-correctness-and-publication.md`](../plans/2026-07-22-viralscan-v3-correctness-and-publication.md).
Do not recreate those details here.

## Next sessions

### Session 1 — SCI-03 preregistration freeze

Start with `SCI-03` in `PLAN.md`. Without inspecting holdout outcomes, freeze:

- biological-sample train/holdout partitions and deterministic seeds;
- the continuous cell score and evidence-tier calibration/search rule;
- numeric simulation factors and limit-of-detection estimation;
- molecule/cell/sample metrics and biological-sample bootstrap uncertainty.

Update `analysis/v3_validation/protocol.yaml`, its schema and adversarial tests.
Run the draft, training, and holdout validators; training/holdout should remain
closed until their explicitly documented prerequisites exist. Obtain an
independent contract review before marking the item complete.

### Session 2 — GOV-03 artifact inventory

This can proceed independently of SCI-03. Create
`analysis/v3_artifact_inventory.tsv` for every claim-bearing input, reference,
intermediate, result, and scheduler record. Record location, public/private ship
scope, software/counting version, input/reference hashes, artifact SHA-256, and
rebuild eligibility. Extend claim-registry validation only as needed; never add
private clinical outputs or institutional paths to ship scope.

### Session 3 — REL-01 honest pip tier

Move full-workflow-only dependencies (especially Snakemake/native-tool concerns)
out of mandatory pip runtime requirements. Make `doctor --profile pip` test only
the Python/API/reporting/validation tier, while conda/container profiles retain
the full toolchain contract. Verify wheel and sdist clean installs and packaged
schemas, then update installation docs and the tracker.

### After those sessions

Execute `SCI-04` to register exact comparator rows and environments before any
outcomes are viewed. The current v3 matrix includes STARsolo, traditional
alignment, Venus, and Viral-Track. VIRTUS/VIRTUS2 has not been registered or run;
decide whether to add it during SCI-04, not after benchmark results exist.

## Restart checklist

```bash
git switch codex/viralscan-v3
git status --short --branch
git rev-parse HEAD
sed -n '1,180p' PLAN.md
python3 scripts/check_data_governance.py
python3 scripts/validate_v3_protocol.py
```

Before committing any future slice, run its focused tests and then the repository
gate recorded in `PLAN.md`. Update the task checkbox, evidence log, next-action
pointer, and claim registry in the same atomic commit.

## Suggested skills

- `review-clear-execute` for executing the accepted tracker slice.
- `tdd` for implementation tasks with new behavior.
- `atomic-commits` when checkpointing completed slices.
- `i-have-adhd` for short, next-action-first progress updates.

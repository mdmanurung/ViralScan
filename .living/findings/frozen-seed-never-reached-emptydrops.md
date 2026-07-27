# F-008 — The frozen cell-calling seed never reached emptyDrops

**Date**: 2026-07-27
**Status**: fixed
**Severity**: major
**Found by**: mycelium six-agent code review of `codex/viralscan-v3` vs `main`

## What was wrong

`analysis/v3_validation/protocol.yaml` freezes `seeds.cell_calling: 20260727002`
and declares
`harmonization.cell_universe.host_only_emptydrops.seed_source: seeds.cell_calling`.
The plumbing existed on both ends and was severed in the middle:

- `src/viralscan/scripts/emptydrops.R:65` calls `set.seed(seed)` from `args[[6]]`.
- `emptydrops_cells` appended `str(seed)` as the sixth argv element.
- `call_cells` passed `rscript`, `fdr`, `lower`, and `niters` — and stopped.

So DropletUtils always ran at the signature default `seed=100`. `RunConfig` had
no seed field at all, and `emptydrops_niters` was likewise absent, reaching R only
through a `getattr(config, "emptydrops_niters", 10000)` fallback that could never
resolve to anything but the literal.

## Why it mattered

emptyDrops is a Monte-Carlo test: the seed decides which barcodes land on the FDR
boundary. The shared cell anchor is computed once per dataset and reused across
every workflow row, so it is the denominator under every cell-level metric —
including `counts_unique` molecule precision/recall/F1, the single axis on which
the v3-versus-2.2.0 improvement claim is allowed to rest. The anchor was therefore
not reproducible from the frozen protocol, which is exactly what the protocol
promises it is.

## Why eleven review rounds missed it

They read `protocol.yaml`. A frozen constant that no call site consumes is
invisible from the protocol side — the declaration is present and correct — and
invisible from the output side, because the run succeeds and produces plausible
cells. Only reading the consuming code exposes it.

## Fix

`emptydrops_seed` and `emptydrops_niters` are declared configuration
(`DEFAULTS` + `RunConfig` + `--emptydrops-seed` / `--emptydrops-niters`), and
`emptydrops_cells` takes keyword-only required parameters with no defaults, so
omitting one is a `TypeError` rather than a silent substitution.

Tests: `tests/test_cellcalling.py::TestAutoCellCalling` —
`test_configured_seed_and_niters_reach_emptydrops`,
`test_emptydrops_cells_requires_an_explicit_seed`,
`test_emptydrops_seed_is_forwarded_to_the_r_command` (asserts the argv position
the R script reads).

## Generalisation

Every entry under `seeds:` in the protocol needs a "does this reach its consumer"
test. See `.living/outputs/reviews/2026-07-27-branch-codex-viralscan-v3-tripwires.md`
(tripwire2, slug `frozen-seed-reaches-consumer`).

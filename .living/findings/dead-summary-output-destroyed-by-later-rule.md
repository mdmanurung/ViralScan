# F-009 — summary.txt totals were computed and destroyed on every run

**Date**: 2026-07-28
**Status**: fixed
**Severity**: minor (no wrong number published; a correct one was never published)
**Found by**: SW-04 artifact-by-artifact audit of `rerun-multimap`

## What was wrong

Two Snakemake rules opened the same file with a truncating mode:

- `src/viralscan/scripts/multimap.py` — `open(f"{config.output}/summary.txt", "w")`,
  writing three totals: viral molecules in the unique matrix, total viral
  molecules under the selected method, and cells with viral reads.
- `src/viralscan/scripts/detection.py:791` — `open(f"{config.output}/summary.txt", "w")`,
  writing entirely different content.

`detection` runs after `multimap` in the DAG (`create_config → kb_count →
analysis → multimap → detection → umap`), so it truncated the file and multimap's
three totals never survived to disk. They were computed at `multimap.py:314-315`
and referenced nowhere else.

## Why it survived

Nothing consumed them. A repo-wide search for the literal strings found no test,
no documentation, no notebook, and no downstream parser — only the write site
itself. The numbers looked published to anyone reading `multimap.py`, and their
absence looked like they had never been requested to anyone reading
`summary.txt`. No test could fail, because no test knew they should exist.

## Why it mattered

"Total viral molecules (selected method)" is the headline quantity of the whole
tool, and it is method-dependent. A user switching allocation methods would look
for it in the summary and not find it.

## Fix

`detection` is now the sole writer of `summary.txt` and derives the three totals
from the H5AD via `_headline_totals`, rather than receiving them from multimap.
Deriving them at the later rule has a second benefit: `rerun-multimap`'s fast path
swaps the selected layer *without* re-running multimap, and because detection
always re-runs, the totals track the current layer. That is the `SW-04`
"summaries" clause satisfied by construction rather than by remembering to
regenerate.

`multimap` now logs the same totals instead of writing them.

Tests: `tests/test_detection.py::TestHeadlineTotals` — four cases, including
`test_totals_follow_the_selected_layer`, which is the rerun-staleness guard.

## Generalisation

Grep for every `open(..., "w")` against a shared output path and check the DAG
order of its writers. Two rules writing one file with truncating mode is a
last-writer-wins race decided by rule order, and it is invisible from either
call site alone.

# The pinned viral-panel Zenodo DOI is unregistered, breaking `viralscan data fetch` for all users

**Date**: 2026-07-26
**Branch**: `codex/viralscan-v3`
**Tracked as**: `REF-11` in `PLAN.md`
**Status**: open, owner-gated

## Finding

`src/viralscan/data_fetch.py` pins:

```python
VIRAL_DATA_DOI = "10.5281/zenodo.20112332"
VIRAL_DATA_RECORD_ID = "20112332"
```

That record does not exist.

- `https://zenodo.org/api/records/20112332` returns
  `{"status": 404, "message": "The persistent identifier is not registered."}`
- `https://doi.org/10.5281/zenodo.20112332` returns 404

## Why this is not a network problem

Both probes ran from `res-hpc-exe029`, which has working outbound HTTPS. An unrelated
third-party Zenodo DOI referenced in `scripts/ttv_public_datasets.json`
(`10.5281/zenodo.14408301`) returned 200 from the same host in the same session. The
record is missing; egress is fine.

## Blast radius

The 195 panel GTFs are still present in `src/viralscan/data/` in the source tree, but the
installed package carries only `anellovirus_accessions.tsv` and `__init__.py` — verified
against `benchmark_runs/legacy_v2_v3/env_full/lib/python3.11/site-packages/viralscan/data/`.
The panel was moved out of the wheel on the assumption the archive would be fetched.

Consequently `viralscan data fetch` cannot succeed for any user, and any bundled-panel run
from a clean install fails at config creation with:

```
ViralScanDataError: Viral reference annotations were not found at ~/.cache/viralscan/data.
Run `viralscan data fetch` before using the bundled viral reference panel, or pass custom
annotations with -gtf.
```

This is how it surfaced: all five v3 rows of the `GOV-06` fresh-control attempt-2 arrays
died on it before reaching quantification.

Blocks `REL-05` (clean-install `validate-run`), `DOC-05` (documented quickstart from a clean
install), the `G2` distribution gate, and the `GOV-06` attempt-3 cache pin.

## Resolution paths

1. **Owner action (correct fix):** publish the panel archive and register the DOI, or
   correct the pinned record identifier if the archive lives elsewhere.
2. **Local unblock (needs a decision):** build the cache from the GTFs still in
   `src/viralscan/data/`, pinned by content hash, with provenance recorded as repo source at
   a Git SHA. This requires relaxing `cache_valid()`, which currently hard-requires
   `manifest["doi"] == VIRAL_DATA_DOI`. Writing the unregistered DOI into a locally built
   manifest was rejected: it would fabricate provenance for content that never came from the
   archive.

## Guard to add

A CI check that resolves the pinned DOI. A placeholder identifier committed ahead of
publication is indistinguishable from a working one until a user — or five cluster jobs —
runs it.

# Reference updates and manifest compatibility

> **Draft for maintainer approval (PLAN `OPS-02`).** Cadences are proposals;
> the invariants in the first section are not.

## Invariants

1. **A frozen reference is never changed in place.** Once a reference build,
   panel archive or accession list is released or cited under a version label,
   its bytes and its identifier do not change. Any change, however small
   (one accession added, one sequence re-masked, one GTF corrected), produces a
   new version label and a new checksum.
2. **Identity is the checksum, not the name.** A reference is identified by the
   `fasta_sha256` (and per-sequence `sha256`) in its `reference_manifest.json`,
   and the Zenodo archive by its DOI and SHA-256. A matching name with a
   different checksum is a different reference.
3. **Results are tied to the reference they used.** Each run writes a
   `reference_manifest.json` next to its reference inputs; report the
   manifest's `fasta_sha256` with any result. Re-running against a newer
   reference is a new run, not a replacement.
4. **Superseded versions stay retrievable.** Old reference versions remain
   downloadable (Zenodo keeps each version) and are marked superseded, with a
   reason, rather than deleted.
5. **Pre-registered analyses keep their reference.** Analyses whose protocol
   freezes a reference digest keep it; changing it needs a recorded deviation.
6. **Thresholds are not retuned retrospectively** when a reference is refreshed
   (PLAN `OPS-05`).

## Update cadence (proposed)

| Trigger | Action |
|---|---|
| Scheduled review `<proposed: every 6 months>` | Re-check source accessions in RefSeq/GenBank for withdrawn, suspended or superseded records; compare against the manifest; publish a new version only if something material changed |
| A reported error (wrong annotation, wrong taxon, contaminated record) | Treated as a scientific incident; new reference version plus a note on which prior results are affected |
| New viral taxa of interest | Added as a new version through the profile mechanism (`--profile`); never silently into an existing one |
| Code release | Does not by itself change a reference; panels are reviewed independently of code releases, as `SUPPORT.md` states |

Each new reference version ships with: the manifest, a changelog entry listing
added, removed and changed records, the new checksums, and the previous
version it supersedes.

## Manifest compatibility

Manifests carry a `schema_version` (`MAJOR.MINOR.PATCH`). ViralScan checks it
on read.

| Change in the manifest schema | Version bump | Reader behaviour |
|---|---|---|
| New optional field | MINOR | Older readers ignore it; newer readers treat absence as "not assessed", never as a pass |
| Field renamed, removed, or meaning changed; new required field | MAJOR | A reader refuses a MAJOR it does not know, with an explicit error naming the expected version, instead of guessing |
| Wording or documentation only | PATCH | No effect |

Rules:

- A reader must reject a manifest with an unknown MAJOR and must fail closed
  (error), not fall back to a default or skip the check.
- A missing field that a check depends on is reported as "not assessed".
- Schema JSON files under `schemas/` are versioned with the package; a schema
  change is listed in `CHANGELOG.md` and follows the deprecation policy in
  [release policy](release_policy.md).
- The reference-build manifest and the packaged panel manifest are versioned
  independently of the package version; a package release states which manifest
  schema versions it reads.

Implementation notes (current code): `reference_manifest.json` is written with
`schema_version` `3.0.0` and validated against
`reference_manifest.schema.json`; the index build manifest is read with a
MAJOR check (`load_build_manifest` in `src/viralscan/virus_identity.py`).
Cadence enforcement is a maintainer process, not code.

New reference builds record `reference_manifest.json` and
`reference_reproducibility.json`. The production `--reference` path applies the
same default viral masking and low-complexity gate as `build-ref`: dustmasker
windows 64 and 30 at level 30, followed by the targeted k-mer mask. Host cDNA is
identified from its GTF transcript rows and is preserved; GTF coordinates remain
unchanged. Raw and prepared FASTAs both undergo duplicate ID/sequence checks.

Per-sequence provenance distinguishes a recorded retrieval date from the build
date. New NCBI and Ensembl fetches save a checksum-bound retrieval receipt.
Existing caches without a receipt retain an unknown retrieval date; their file
timestamps are not substituted. NCBI taxonomy comes from the cached GenBank
organism and taxon qualifiers. Local FASTAs retain `local_input` provenance with
unknown taxonomy unless source metadata is supplied. Cluster representatives
and excluded members retain CD-HIT decisions; partial-panel fetch failures retain
their accession and reason.

Source licences are structured records with `status: unreviewed`, empty terms,
and an explicit review requirement. This records the remaining review; it does
not grant redistribution rights. Reviewed source terms remain a maintainer task.

`content_sha256` and the reproducibility audit cover sequence/provenance content
and deterministic FASTA/GTF/t2g checksums. Output paths, build/retrieval dates,
software revision, commands, resource receipts and the binary index checksum
are kept outside that digest. The build receipt records child CPU time and
elapsed time for `kb ref`, the peak RSS (`peak_rss_kib`, the high-water mark over
every child process the build reaped, not `kb` alone) and the path, version and
SHA-256 of `kb`, `kallisto`, `bustools`, `dustmasker` and `cd-hit-est`. Toy tests
establish deterministic preparation and mocked t2g content. They do not prove
binary kallisto index reproducibility or a real full-panel rebuild.

<!-- MAINTAINER: approve the proposed six-month reference review cadence before OPS-02 closes. -->

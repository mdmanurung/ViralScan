# Release, support and security policy

> **Draft for maintainer approval (PLAN `OPS-01`).** Contact routes below are
> placeholders until maintainer metadata is settled (PLAN `REL-15`). Anything
> tagged **DRAFT: maintainer decision** is a policy choice that the repository does not fix;
> untagged statements restate `SECURITY.md`, `SUPPORT.md`, the release workflow or the
> pins listed under "What the repository pins today". The
> root-level [`SUPPORT.md`](https://github.com/mdmanurung/ViralScan/blob/main/SUPPORT.md)
> and [`SECURITY.md`](https://github.com/mdmanurung/ViralScan/blob/main/SECURITY.md)
> remain the short statements; this page is the full policy behind them.

## Contact routes

**DRAFT: maintainer decision** (PLAN `REL-15`): every route other than GitHub Issues and
GitHub Security Advisories.

| Purpose | Route |
|---|---|
| Reproducible bugs, feature requests | GitHub Issues |
| Suspected vulnerabilities | GitHub Security Advisories (private); fallback `<maintainer security contact>` |
| Scientific incidents (false positive/negative, data loss, reference drift) | GitHub Issues with a reproducible packet; restricted data by `<maintainer contact>` only |
| Citation, licence, authorship questions | `<maintainer contact>` |

Never attach human-subject FASTQs, H5AD files, clinical metadata, credentials or
institutional paths to any public route.

## Supported versions

Stated in `SECURITY.md`: 3.0 development receives security fixes during the release-candidate
period, the latest stable 3.x line receives security and correctness fixes once 3.0 is stable,
and pre-v3 scientific output formats are unsupported. The finer split below is a draft.

| Release line | Status | Receives |
|---|---|---|
| 3.0 pre-releases (`3.0.0.dev*`, release candidates) | Active development | All fixes; no stability promise between candidates |
| Latest stable 3.x | Supported | Security and correctness patches |
| Earlier stable 3.x minors | Not supported once a newer 3.x minor is stable (DRAFT: maintainer decision) | Upgrade to the latest 3.x |
| 2.x and earlier | Scientific output unsupported | Nothing; counts must be rebuilt (see [migration](migration.md)) |

Only the latest stable minor line is supported. A support window for the previous minor
line is **DRAFT: maintainer decision** and is not promised here.

## What the repository pins today

These are read from the files named, not chosen here. A support promise for any of them is
**DRAFT: maintainer decision**.

| Item | Pinned value | Source |
|---|---|---|
| Package version | `3.0.0.dev1` | `src/viralscan/__init__.py` |
| Python (pip tier) | `>=3.10`; CI runs 3.10, 3.11 and 3.12 on ubuntu and macOS | `pyproject.toml`, `.github/workflows/ci.yml` |
| Full workflow toolchain | Python 3.11 on `linux-64`; `kb-python` 0.30.2, kallisto 0.52.0, bustools 0.45.1, STAR 2.7.11b, snakemake 9.23.1 | `environment.yml` |
| Release channel | tag `vMAJOR.MINOR.PATCH` publishes to PyPI and may advance the container `latest`; any other tag goes to TestPyPI | `.github/workflows/release.yml` |
| Manifest schema | `reference_manifest.json` `schema_version` accepts `3.0.0` and `3.0.1` | `schemas/v3/reference_manifest.schema.json` |

## Versioning

**DRAFT: maintainer decision**: the mapping of change types to version parts below. Releases use
`MAJOR.MINOR.PATCH` (the release workflow treats only a tag of exactly that shape as stable).

- **PATCH**: bug, security and correctness fixes; documentation. It must not
  change a documented output column, its meaning, or a default that alters
  scientific results.
- **MINOR**: new flags, subcommands or output columns that are additive and
  default-off or default-neutral. A change to a scientific default (for example
  which reference or multimapping method is used by default) is MINOR at the
  least, is recorded in `CHANGELOG.md`, and for pre-registered analyses needs a
  deviation record.
- **MAJOR**: removal of a deprecated interface, or a change to the molecule
  counting contract or an output schema major version.

## Patch releases

The release workflow requires the latest completed, successful `ci.yml` push
run on `main` for the exact tagged commit. It builds distributions once,
checks their checksums, and gives the same artifacts to downstream jobs.
The OCI image is scanned before publication; its digest supplies the generated
Apptainer definition. This workflow preparation does not establish runtime parity.

Maintainers must resolve immutable Action pins and configure Trusted Publishers
and protected `testpypi`/`pypi` environments in repository settings. A workflow
environment declaration does not enable required reviewers. Stable publication
uses PyPI and may advance the container's `latest`; candidates use TestPyPI and
keep `latest` unchanged. Licence and vulnerability reports require review.

The triggers for a patch release are **DRAFT: maintainer decision**. Proposed: a security fix; a correctness bug that changes
reported viral counts, calls or denominators; data loss or a silent wrong-sample
result; an installation break on a supported platform.

Every patch release must:

1. Pass the same gates as a minor release (lint, tests, `viralscan doctor`,
   docs-consistency tests); no gate is skipped for urgency.
2. Add a `CHANGELOG.md` entry stating the affected versions and, for a
   correctness fix, whether previously produced results should be re-run.
3. Be built from a tagged commit with checksums recorded for the release
   artefacts.

If a correctness bug invalidates earlier results, the advisory states which
versions and which output tables are affected.

## Deprecation

**DRAFT: maintainer decision**: the whole section, including the minimum of one MINOR release
before removal. No deprecation mechanism exists in the code yet; this describes the intended one.

- A CLI flag, output column, config key or Python API name is deprecated
  by (a) emitting a warning that names the replacement and the first release in
  which it may be removed, and (b) listing it in `CHANGELOG.md` and in the
  affected reference page.
- A deprecated interface keeps working for at least one MINOR release after
  the release that deprecates it, and is removed only in a release that
  `CHANGELOG.md` flags as removing it.
- Scientific-output changes are not "deprecated quietly": a changed meaning of
  an existing column is a schema-version change (`schema_version`), not a
  deprecation.
- Security or correctness fixes may remove an unsafe behaviour immediately;
  the changelog says so.

## Security response

Scope: code in this repository, the packaged data files, and the release
artefacts. Third-party tools (`kb`, `snakemake`, STAR, kallisto) are reported
upstream, with a note in the advisory if ViralScan needs a version pin.

| Step | Target |
|---|---|
| Acknowledge a private report | within 7 days (as stated in `SECURITY.md`) |
| Triage and confirm or reject, with severity | DRAFT: maintainer decision: triage target |
| Fix, patch release, advisory published | DRAFT: maintainer decision: fix target; coordinated with the reporter |
| Credit | reporter named in the advisory unless they decline |

Reports are handled privately until a fixed release exists. Fixes ship as a
patch release of the latest stable line; older lines are not backported.

<!-- MAINTAINER: approve every item tagged DRAFT: maintainer decision (private contact, support windows, patch triggers, deprecation window, security response targets) before OPS-01 closes. -->

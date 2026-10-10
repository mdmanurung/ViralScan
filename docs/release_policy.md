# Release, support and security policy

> **Draft for maintainer approval (PLAN `OPS-01`).** Contact routes below are
> placeholders until maintainer metadata is settled (PLAN `REL-15`). The
> root-level [`SUPPORT.md`](https://github.com/mdmanurung/ViralScan/blob/main/SUPPORT.md)
> and [`SECURITY.md`](https://github.com/mdmanurung/ViralScan/blob/main/SECURITY.md)
> remain the short statements; this page is the full policy behind them.

## Contact routes

| Purpose | Route |
|---|---|
| Reproducible bugs, feature requests | GitHub Issues |
| Suspected vulnerabilities | GitHub Security Advisories (private); fallback `<maintainer security contact>` |
| Scientific incidents (false positive/negative, data loss, reference drift) | GitHub Issues with a reproducible packet; restricted data by `<maintainer contact>` only |
| Citation, licence, authorship questions | `<maintainer contact>` |

Never attach human-subject FASTQs, H5AD files, clinical metadata, credentials or
institutional paths to any public route.

## Supported versions

| Release line | Status | Receives |
|---|---|---|
| 3.0 pre-releases (`3.0.0.dev*`, release candidates) | Active development | All fixes; no stability promise between candidates |
| Latest stable 3.x | Supported | Security and correctness patches |
| Earlier stable 3.x minors | Not supported once a newer 3.x minor is stable | Upgrade to the latest 3.x |
| 2.x and earlier | Scientific output unsupported | Nothing; counts must be rebuilt (see [migration](migration.md)) |

Only the latest stable minor line is supported. `<support window, if any, for the previous minor>`
is a maintainer decision and is not promised here.

## Versioning

Releases use `MAJOR.MINOR.PATCH`.

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

A patch release is cut for: a security fix; a correctness bug that changes
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
| Triage and confirm or reject, with severity | `<triage target, maintainers to set>` |
| Fix, patch release, advisory published | `<fix target, maintainers to set>`; coordinated with the reporter |
| Credit | reporter named in the advisory unless they decline |

Reports are handled privately until a fixed release exists. Fixes ship as a
patch release of the latest stable line; older lines are not backported.

<!-- MAINTAINER: approve private contact, support windows and security response targets before OPS-01 closes. -->

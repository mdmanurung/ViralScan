# Maintainer readiness checklist

**Draft for maintainer fill-in (PLAN `REL-24`; grill session 2, Q1-Q18).**
Unchecked fields are unresolved. This draft grants no publishing authority,
reserves no identifier and makes no new support or response-time promise.
Record the owner, decision date and evidence link for each completed item.

## Identity and contacts (`REL-15`, `MS-04`)

- [ ] Final software authors and citation order: `<fill in>`
- [ ] Licence holder and approved licence text: `<fill in>`
- [ ] Maintainer name and public contact: `<fill in>`
- [ ] Author ORCIDs, affiliations and CRediT roles: `<fill in>`
- [ ] Funding, conflicts, acknowledgments, ethics/consent and data-access wording: `<fill in>`

Update [CITATION.cff](CITATION.cff), [pyproject.toml](pyproject.toml) and
[LICENSE](LICENSE) only after the factual text is approved. Private dataset
results remain excluded unless their access and publication requirements are met.

## Repository and publishing settings (`RC-02`, `REL-14`)

- [ ] `main` branch protection and required CI checks verified by repository owner: `<evidence>`
- [ ] Protected `testpypi` and `pypi` environments, reviewers and permitted deployment refs: `<evidence>`
- [ ] Trusted Publisher owner/repository/workflow/environment registered with TestPyPI and PyPI: `<evidence>`
- [ ] GHCR, GitHub release and public artifact-hosting permissions approved: `<evidence>`

Record settings and links only; never paste tokens or other secrets here.
The [release workflow](.github/workflows/release.yml) declarations do not configure
service accounts or required reviewers. No `v*` tags before the RC gates pass.

## Support and reference policies (`OPS-01`, `OPS-02`, `REL-11`)

- [ ] Supported versions and support window approved: `<decision>`
- [ ] Private security contact, triage/fix targets and existing response policy reviewed: `<decision>`
- [ ] Reference review cadence and responsible maintainer approved: `<decision>`
- [ ] Weekly lock-audit reviewer and evidence-review route assigned: `<owner and route>`
- [ ] Vulnerability exceptions or supported-Python decision documented per affected pin: `<decision and evidence>`

Review [release policy](docs/release_policy.md),
[reference maintenance](docs/reference_maintenance.md), [SECURITY.md](SECURITY.md)
and [SUPPORT.md](SUPPORT.md) together. Their draft placeholders remain unresolved
until approval. The [weekly audit](.github/workflows/lock-audit.yml) reports findings
without blocking CI or release; failed/skipped scans require review and are not
clean results. It changes no lock pins. Conda OSV matches are approximate;
release OCI scanning remains a separate gate.

## Candidate and archives (`RC-01`–`RC-05`, `STB-02`)

- [ ] Exact reviewed main SHA and RC-gate receipts: `<SHA and evidence>`
- [ ] External laboratories, contacts and tester-packet owner: `<fill in>`
- [ ] Testers receive a pinned-commit installation reference before release artifacts exist: `<SHA and packet>`
- [ ] Reference archive DOI/publication status verified (`REF-11`): `<evidence>`
- [ ] Separate software and benchmark/archive DOI owners assigned: `<fill in>`

Use [PLAN.md](PLAN.md) as the gate and status authority. DOI reservation, archive
publication, credentials, tester recruitment and final sign-off remain maintainer
work. Completing this draft alone does not qualify a release or scientific claim.

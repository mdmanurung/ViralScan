# Support and version policy

| Release line | Status | Scope |
|---|---|---|
| 3.0 development | Active | Correctness, security, installation, and documentation |
| Latest stable 3.x | Supported | Security and correctness patches |
| 2.x and earlier | Scientific-output unsupported | Source remains archived; corrected counts must be rebuilt |

Use GitHub Issues for reproducible bugs and feature requests. Include
`viralscan --version`, `viralscan doctor --profile full`, the redacted
`run_manifest.json`, and `viralscan validate-run` output. Never upload restricted
human-subject data. Reference panels are reviewed independently of code releases;
their manifests and checksums define their support status.

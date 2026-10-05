# Migrating from earlier versions

ViralScan v3 does **not** migrate legacy outputs. There is no converter, and
none is planned: the v3 counting contract, identity table and evidence tiers
are not reconstructable from pre-v3 result files.

- **Rebuild only.** Rebuild the reference with `viralscan build-ref` and rerun
  v3 from the original FASTQ files.
- **Do not compare or carry over old values.** Result values from earlier
  versions (counts, percentages, virus lists) must not be compared with v3
  output or copied into v3 tables, reports or manuscripts. Differences reflect
  changed definitions, not just changed numbers.
- **Old run directories** are not accepted by `viralscan validate-run` as v3
  runs; keep them as archives and write new results to a fresh `--output`.

See the [output reference](output_reference.md) for the v3 schema and the
[FAQ](faq.md) for which comparisons across strategies are valid.

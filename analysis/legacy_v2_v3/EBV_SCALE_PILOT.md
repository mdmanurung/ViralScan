# SRR12682296 scale-pilot comparison

Status: **descriptive comparison complete; causal attribution prohibited**

The frozen scale pilot processed Emma's archived ViralScan 2.2.0
`SRR12682296` BUS with the archived 10x v2 whitelist and the frozen v3
host-conservative molecule allocator. The scientific run completed
successfully in one allocator pass. SLURM job `25329992` was marked failed only
because a post-run shell status-print expression had invalid quoting; the
retained wrapper and v3 status records are both `success`, all input/source
hashes match before and after, and the corrected post-check returns zero.

## Scale-pilot evidence

- Raw BUS records: 97,812,898.
- Whitelist on-list records before molecule resolution: 95,758,526.
- Cells: 414,715; genes: 41,409; equivalence classes: 225,802.
- Input CB-UMI molecules: 53,709,866.
- Unique molecules: 50,952,132.
- Ambiguous molecules: 1,388,191.
- Unresolved molecules: 1,369,543.
- Resolved/selected matrix mass: 52,340,323.
- Ignored repeated-read multiplicity: 29,743,645.
- Allocation time: 249.36 seconds.
- Full wrapper wall time: 6:25.66.
- GNU-time peak RSS: 4,231,476 KiB.
- Scratch written before cleanup: 5,370,306,338 bytes; scratch was cleaned.
- All 99 reconstructed legacy feature totals matched the archived summary.

Artifact SHA-256:

- scale result:
  `6cb8a1228d41616301365c90c5eab4e5c521ae75ed4121c60eb9f4734f08b6e5`;
- scale count audit:
  `42382f00ed43b8f972181cda942b6db0c4a444994e00d004910f16b3b1a6214b`;
- scale v3 status:
  `133bc71ad5bc80f73e2a1f78684351269fe403a7f85a0eb3e434a889a23d9883`.

## Retained v3 baseline comparison

The retained baseline is
`analysis/v3_ebv_baseline/host_conservative.json` at SHA-256
`2863487425742175e73c49709df22f746918631dcb647e36a46a21cb2d6f221e`.
It was recorded by commit `9b9913b` from a different fresh combined run
directory. Its benchmark script explicitly called `prepare_resolved_bus` with
`whitelist=None`.

| Metric | Retained no-whitelist baseline | Frozen archived-input scale pilot |
|---|---:|---:|
| BUS records | 103,145,071 | 97,812,898 |
| Cells | 848,191 | 414,715 |
| Genes | 43,451 | 41,409 |
| Equivalence classes | 388,677 | 225,802 |
| Input molecules | 57,957,364 | 53,709,866 |
| Unique molecules | 51,976,300 | 50,952,132 |
| Ambiguous molecules | 4,220,847 | 1,388,191 |
| Unresolved molecules | 1,760,217 | 1,369,543 |
| Selected matrix mass | 56,197,147 | 52,340,323 |
| Ignored read multiplicity | 31,129,183 | 29,743,645 |

Relative to the retained baseline, the scale pilot has 5.17% fewer BUS
records, 51.11% fewer cells, 7.33% fewer input molecules, 67.11% fewer
ambiguous molecules, and 6.86% lower selected matrix mass.

These differences cannot be attributed to the whitelist alone. The two runs
use different retained BUS inputs, cell and feature universes, equivalence
class sets, and preparation histories. They therefore do not estimate package
superiority or a whitelist effect. The valid conclusion is narrower: the
earlier baseline omitted a whitelist and is not an interchangeable comparator
for the frozen archived-input pilot.

## Retained H5AD is a legacy-count comparator

The earlier retained `adata_multimap.h5ad` from the fresh combined
SRR12682296 run is not the July-22 molecule-unit JSON baseline. Its recorded
configuration has `whitelist: null`, and kb-python records `-w None` with no
`bustools correct` step. At the nearest supported source snapshot
`9dbf4a6557019357e6820df6a5fe2920482c33e9`, multimapping consumes raw BUS rows,
does not read the UMI column, and scales allocations by BUS `count`.
`counts_original` comes from the separate `bustools count --umi-gene` matrix;
`counts_combined` therefore combines layers with legacy/read-multiplicity
semantics rather than the frozen v3 CB-UMI molecule contract.

The July-4 job did not record an execution SHA or dirty-worktree state, so
`9dbf4a...` is the nearest supported snapshot rather than cryptographic proof
of the exact executed source. The retained H5AD remains useful for reconstructing
the original calls, but it is not a molecule-unit truth target and must not be
used to claim v3 numerical parity.

## Resource decision

The below-4-GiB tier remains adequate for this row. The 3.13 GB raw BUS used
4.04 GiB peak RSS, 5.00 GiB of scratch output, and 6.5 minutes wall time under
a 32 GiB, 6-hour request. This supports retaining the prespecified tier for
similar rows, but it does not validate the at-least-4-GiB tier; those rows
remain gated by their own resource limits and the full-array review.

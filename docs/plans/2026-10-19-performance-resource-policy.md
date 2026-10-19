# Resource policy proposal — pending separate authorization

No resource or scheduling change is implemented by the performance finalization.
B11/B30 retain their separate authorization gate. The current request authorizes
byte-preserving optimizations and their local qualification; it does not waive the
plan's prohibition on changing effective workers, scientific parameters or SLURM.

## Verified current requests

Installed help and wrapper source were inspected in the full-test environment.
These are requested worker defaults, not measured concurrent CPU utilization.

| Entry point | Current native worker request | Source of request |
|---|---:|---|
| `kb count`, both whitelist branches | 8 | installed kb default, forwarded to kallisto |
| direct `kallisto bus` with no flag | 1 | native default; distinct from kb count |
| `kb ref`, both package entry points | 8 | installed kb default, forwarded to index |
| reference homology minimap2 | 3 | native default |
| reference CD-HIT | all visible CPUs | explicit `-T 0 -M 8000` |
| multimap bustools sort | configured cores | explicit `-t config.cores` |

Workflow `--cores` defaults to six. `kb_count` and `multimap` lack explicit rule
thread declarations; auto evidence and host capture reserve rule threads but read
`config.cores` directly. Host filtering/anello alignment consume the rule's threads.
Snakemake caps rule threads at global cores and exports OMP/BLAS/NumExpr environment
values from rule threads. Adding a declaration alone can therefore change effective
library workers, even if the native argv still requests eight.

Source reach: `Snakefile`, `cli/parsers.py`, `menu.py`, `runconfig.py`,
`scripts/multimap.py`, `scripts/auto_evidence.py`, `scripts/host_capture.py`,
`scripts/build_reference.py` and their parser/config/DAG/native integration tests.
Installed-source and help receipts are retained in the remaining-wave evidence
packet's `baseline.json` and `cards/b11/tool_defaults.json`.

## Concrete proposed next scope

1. Qualify an explicit legacy profile: `kb count -t 8`, both `kb ref -t 8`,
   reference minimap2 `-t 3`; preserve existing replay/STAR/sort/evidence requests,
   CD-HIT `-T 0 -M 8000`, numerical-library environment and CPU affinity.
2. Compare implicit/explicit real native A/A and A/B with complete raw BUS/EC,
   H5AD, FASTA, cluster, PAF, index/t2g, downstream, error and partial-file receipts.
   Record requested budget, reserved CPUs, native argv, inherited environment,
   affinity and observed workers separately. No resource change is a speedup claim.
3. Propose a total-admission policy only after those receipts. An eight-worker
   request cannot fit a one/six-core allocation without rejection or changing
   workers. Either choice changes current behavior and needs a compatibility
   decision. Bounding CD-HIT likewise changes its existing all-CPU request.

A concrete later policy could define `--cores` as total admitted CPUs and reject
fixed requests exceeding the allocation. This would reject today's accepted
one-core configuration and is **not** an approved implementation. Do not silently
substitute `min(default, cores)` or increase the default core budget.

Required dry-run matrix: 1/6/8/32 cores, both whitelist paths, optional rules,
rerun-multimap and external Snakemake overrides. Mutation checks must detect omitted
worker flags, unscaled config values, numerical-environment drift, changed search/
identity/memory parameters and accidental CD-HIT single-threading.

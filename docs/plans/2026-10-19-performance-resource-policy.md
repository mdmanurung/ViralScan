# Resource policy — `--cores` is the total CPU budget (approved 2026-10-19)

**Approved and implemented 2026-10-19 (`total_budget_v1`, PLAN B11/B30, branch
`codex/cpu-budget`).** `--cores` is the total concurrent CPU budget. Native workers are
bounded by it (never raised): `kb count`/`kb ref` `-t min(8, cores)`, homology minimap2
`-t min(3, cores)`, CD-HIT `-T cores` (`-M 8000` unchanged). `kb_count` and `multimap`
now declare `threads` (kb_count `min(8, cores)`, multimap `cores`), and auto_evidence,
host_capture and multimap read `snakemake.threads`, so `--set-threads` is honoured.
Numerical libraries stay pinned at one worker (kb_count shell exports OMP etc = 1;
multimap.py pins them before numpy loads). `build-ref --cores N` is new; unset keeps the
legacy requests (8 / 3 / all CPUs) so unflagged builds are argv-identical. The workflow
default stays 6, so a default workflow run now uses `kb count -t 6`. Lower worker counts
can change thread-dependent raw bytes; that needs its own compatibility decision, no
comparator was relaxed. Manifest fields for requested/effective workers are not yet
recorded (argv in logs and receipts is the record). The text below is the pre-approval
state.


The user explicitly approved the local qualification phase on 2026-10-19.
The candidate makes installed legacy worker requests explicit: kb count/ref eight,
homology minimap2 three. Existing Snakemake declarations, numerical-library
environment, CPU affinity, CD-HIT settings and scientific options stay fixed.
No changed worker count, scheduling policy, SLURM job or large build is authorized.
See [qualification](2026-10-19-legacy-worker-qualification.md) for current receipts
and their limits. B11/B30 remain partial until the total-admission decision and
its later qualification are resolved.

## Verified current requests

Installed help and wrapper source were inspected in the full-test environment.
These are requested worker defaults, not measured concurrent CPU utilization.

| Entry point | Current native worker request | Source of request |
|---|---:|---|
| `kb count`, both whitelist branches | 8 | explicit `-t 8`, preserving installed kb default |
| direct `kallisto bus` with no flag | 1 | native default; distinct from kb count |
| `kb ref`, both package entry points | 8 | explicit `-t 8`, preserving installed kb default |
| reference homology minimap2 | 3 | explicit `-t 3`, preserving native default |
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

## Approved legacy scope and later policy boundary

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

# G02 representative whole-run proposal — pending approval

> Path placeholders: `$SITE_SCRATCH` = site scratch root (institutional filesystem paths, replaced by placeholders).

No scheduler job is submitted by this document. The authorized code delivery can
finish independently; production stage attribution remains unmeasured until this
separate run is approved and succeeds. No production speedup follows from the local
component ratios.

## Bounded run

Use the retained public development sample SRR12682296, with the exact scientific
configuration in `cost01/rerun_em_cell/SRR12682296/SRR12682296/config.yaml`:
10xv2, knee calling, em-cell, 100 EM iterations, tolerance 1e-6, pseudocount 1,
visualizations enabled, read filter/hostresponse/UMAP disabled. Preserve every
configuration value except a fresh isolated output path. Do not relabel this as a
default-emptyDrops run or as coverage of optional filtering/evidence/model stages.

Read-only inputs, verified present on 2026-10-19:

- R1/R2: `$SITE_SCRATCH/ViralScan/benchmark_inputs/reference_strategy/SRR12682296/SRR12682296_{1,2}.fastq.gz` (2.27GB/8.67GB).
- Index/t2g/GTF: `$SITE_SCRATCH/ViralScan/ebv_latest_ref_2026-09-27/full_panel/{index.idx,t2g.txt,viral_final.gtf}`.

Use pre-performance main `6910c2d` as baseline and the final merged main revision as
candidate. Before execution freeze full commit/runtime/input/binary/config hashes,
confirm the runtime diff contains only the delivered performance changes, and verify
that installed kb/kallisto can read the retained index. If compatibility or another
preflight fails, stop without rebuilding the index or changing tools/parameters.

Proposed maximum: **three sequential SLURM jobs, each 8 CPUs, 32GiB, 6 hours**,
order baseline A, baseline A, candidate B. No retries or additional jobs. Request the
same eligible partition/node constraints as the retained cost run after checking its
receipt; do not choose an unverified partition. Reserve at most 100GiB new output
space; preflight free space and fail before launch if unavailable. Existing data and
outputs are never overwritten or removed. Retained `.start`/`.end` receipts differ by 2377 seconds (39min37sec); this
is historical evidence, not a current run estimate or a measured memory bound.

Run the package's standard Snakemake `all` target with each source's Snakefile,
`--cores 8`, and its copied configuration in a new workspace. Keep native defaults,
CPU affinity and numerical-library environment identical between arms; retain the
actual submitted argv, environment, scheduler accounting, tool paths and requested
versus observed workers. No new reference build, internet call or host-response model fit.
Use writable per-run XDG/Numba/Matplotlib caches; do not flush shared caches.

## Acceptance and limits

Capture per-rule wall/CPU/I/O, peak parent/child RSS, complete logs, partial outputs,
full relative file sets and raw hashes. Establish baseline A/A variability first;
compare raw matrices/tables/plots and compressed/native artifacts. Record truthful
executing-source provenance differences as precise literal spans with original
hashes retained. Unexplained A/A or A/B differences block equivalence; never replace
raw comparison with decompression, sorting, rounding or serialization rewrites.

Three jobs establish bounded whole-run attribution and an A/A check, **not** the
seven-pair/two-session precision required for a new whole-run benefit claim. Any
expanded repetition requires a new bounded approval. Optional stages remain covered
only by their separate local/native qualification. No scientific readiness claim,
held-out data, download, changed worker policy or deletion is authorized by approval
of these three jobs. Resource-policy work has its own linked proposal.

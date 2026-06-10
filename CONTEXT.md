# ViralScan

Domain and architecture vocabulary for ViralScan — a Snakemake-driven CLI that
quantifies viral load from paired-end single-cell FASTQ samples using
`kb-python` (kallisto + bustools). This file pins terms so architecture reviews
and design conversations stay consistent; see `LANGUAGE.md` in the
`improve-codebase-architecture` skill for the architecture vocabulary
(module, interface, seam, depth).

## Language

### A run and its context

**Run**:
A single invocation of viral quantification for one sample pair — from FASTQ
input through the six Snakemake rules to the per-sample output directory.
_Avoid_: job, pipeline run, execution.

**Run Config**:
The validated, typed description of a Run — every parameter the rules need
(paths, thresholds, flags, multimap policy). Coerced and validated exactly once,
at the CLI→workflow edge (`from_cli`); read back by downstream rules as a
trusted typed value (`from_yaml`).
_Avoid_: config dict, settings, params, options.

**Run Context**:
A passive value passed to each worker script's `run(ctx)`, bundling the
**Run Config** with the **Kb Count Outputs** for that Run. Answers "what" and
"where"; it does not perform I/O.
_Avoid_: session, environment, state, handle.

**Kb Count Outputs**:
The module that owns the on-disk layout `kb-python` writes under a Run's output
directory — the named file paths (adata, bus, ec, barcodes, genes…) and the
resolution of which adata is current for a Run.
_Avoid_: paths, file map, output dir.

### Viral evidence

**Viral Accession**:
A viral gene identifier drawn from a reference GTF's `gene_id` field; the unit
ViralScan looks for in a quantified sample.
_Avoid_: viral gene id (when precision matters), virus id.

**Detection**:
Calling a Viral Accession present in a sample because its total UMI count meets
the `detection_threshold` (inclusive `>=`).
_Avoid_: hit, call, positive.

**Multimapper Correction**:
Redistribution of UMI mass from equivalence classes that span multiple genes,
per the configured `multimap_method`. Owned by `multimapping.py`.
_Avoid_: ambiguity resolution, dedup.

## Relationships

- A **Run** has exactly one **Run Config**.
- A **Run Context** bundles one **Run Config** and one **Kb Count Outputs**.
- **Run Config** is validated once at `from_cli`; every later read is `from_yaml`
  against a file only `createconfig` writes.
- **Kb Count Outputs** resolves the current adata for a Run from the
  **Run Config**'s `multimapping` flag — not from file existence. If the flag is
  set the multimap adata must exist; a missing file is an error, not a fallback.
- **Detection** reads counts for each **Viral Accession** from the current adata.

## Example dialogue

> **Dev:** "When `detection` needs the count matrix, does it decide which adata
> to read, or does the **Run Context** hand it the current one?"
> **Maintainer:** "The **Kb Count Outputs** resolves the current adata — there's
> one rule for it, so `detection` and `umap` can't drift. `detection` just asks
> for `ctx.outputs.current_adata`."

## Flagged ambiguities

- "config" was used for three different things: the argparse Namespace, the
  Snakemake `--config` dict, and the YAML on disk. Resolved: the typed value is
  the **Run Config**; the others are representations it is built from, not the
  thing itself.
- "viral gene" and **Viral Accession** were used interchangeably. Resolved: they
  denote the same identifier; prefer **Viral Accession** in interfaces.
- Mapping a **Viral Accession** to a virus name was done two ways — substring
  match (`detection`) and prefix match (`umap`) — which disagreed on 151/2692
  bundled gene IDs and were both buggy. Resolved: one **boundary-aware** rule in
  `virus_grouping.py` (key matches when the gene ID equals it, or starts with it
  and the next char is `_` or a digit; longest key wins).

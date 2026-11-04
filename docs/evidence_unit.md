# Evidence unit for read-level checks (PLAN `EVID-00`)

This note fixes what one unit of read-level evidence is, before the splice check (`EVID-02`), the
gene-strand check (`EVID-03`), the read browser (`EVID-04`) and a later evidence verdict (`EVID-06`)
are built. Each of them counts the same thing, so their numbers share one denominator. The note was
added after the 2026-11-03 plan review (Codex C4, Biomni C3).

## The molecule

A molecule is one corrected (CB, UMI) pair of an on-list barcode, as in `molecule_verdicts.tsv.gz`
of a `viralscan evidence` run. Its **member reads** are the rows of `read_lineage.tsv.gz` with that CB
and UB. `competitive_reads.raw.bam` names each read `<CB>_<UB>_<read_number>|<read_id>`, so a BAM
record joins its lineage row on `read_number` (and `read_id`, the part after `|`).

Every check counts molecules. A read never counts on its own: PCR copies of one molecule are one piece
of support, and a molecule supports a feature (a junction, a strand, a gene) when **any** of its member
reads does. Read counts are reported beside molecule counts, never instead of them.

## Which molecules a check uses

| Category | Definition (from `molecule_verdicts.tsv.gz`) | Used by the checks |
|---|---|---|
| `counted_virus_best` | `counted` is true and `verdict` is `virus_best` | yes: the denominator of every check |
| `counted_host_best` | `counted` is true and `verdict` is `host_best` | no; reported as a count |
| `counted_tie` | `counted` is true and `verdict` is `tie` | no; reported as a count |
| `counted_unaligned` | `counted` is true and `verdict` is `unaligned` | no; reported as a count |
| `not_counted` | `counted` is false | no; reported as a count |

Each check writes the five counts, so a reader can see how many molecules it did not look at.

The single representative read per molecule (`viral_reads.umi_dedup.fasta`) is **not** an input: it
is the read with the best viral alignment, so it favours the virus when a molecule mixes viral-looking
and better-scoring host reads (`docs/output_reference.md`, `viralscan evidence` output).

## When there is no unit

- **Competitor mode** (`--competitor-fasta`): no host-versus-virus verdict exists
  (`molecule_verdict = "not_produced_in_competitor_mode"`). Every check that needs the unit writes
  `not_assessed` with that reason and computes nothing.
- **Manifests without `minimap2_options`** (before AUD3-04): their BAMs hold one record per read, so
  `tie` cannot occur and exact ties fall arbitrarily to one side. Checks run, and their outputs carry
  `verdict_ties_resolvable = false`.
- **No counted virus-best molecule**: the check writes its counts and `not_evaluable`.

## What the unit does not cover

- Molecules the kallisto index never captured (for example the HHV-6B 3′ flank reads of F-027)
  are outside every evidence run. A post-hoc check cannot recover them.
- The unit says which molecules align better to the virus than to the host reference given. It does
  not establish infection, reactivation or cell of origin: ambient RNA and integrated genomes
  (ciHHV-6) give viral molecules too.

## Outputs every check carries

Each check's table and receipt record: the evidence run's `evidence_manifest.json` SHA-256, the five
category counts, the read counts beside them, and `not_assessed` / `not_evaluable` with the reason
when the unit is missing.

## Gene-relative strand (`EVID-03`)

`strand_qc.tsv` uses the `counted_virus_best` member reads as its unit and writes the five category counts
(molecules and reads) in every row. It also reads the member reads of `counted_host_best` molecules, to
calibrate the expected read orientation of the library on host exons that belong to one gene (at least 200
blocks, 80 % agreeing; otherwise `not_assessed`). This is the only check that uses host-best molecules, and
only for calibration. The table cannot carry the SHA-256 of `evidence_manifest.json`, which is written after
it; the manifest lists the hash of `strand_qc.tsv` instead. The default of `--strand` is not decided by this
check.

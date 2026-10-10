# Panel integrity checks

The bundled reference builder checks the local GTF corpus hash manifest before
downloads. The default is
`analysis/panel_expansion/gtf_corpus_manifest.tsv`; use `--gtf-manifest` to select
another reviewed corpus. The manifest records hashes and sizes, with unknown
source dates. It identifies local files and does not imply GTF bundling in the
wheel or validate their biological annotations.

Run the offline structure and catalogue matrix check with:

```bash
PYTHONPATH=src python scripts/panel_integrity.py --output /tmp/panel_integrity.tsv
PYTHONPATH=src python scripts/write_gtf_manifest.py \
  --check analysis/panel_expansion/gtf_corpus_manifest.tsv
```

The matrix includes every shipped accession, deliberate index exclusions, and
reverse GTF accessions absent from the catalogue. Real anellovirus gene-table
exons are included because the builder emits them during assembly. Missing gene
annotations remain visible even when the builder can create a whole-genome
fallback from its FASTA. Duplicate records, duplicate gene/transcript ownership,
malformed rows and transcripts without exons fail the structure check. These
checks are separate from hashing: a matching hash does not repair annotation.

The current local corpus has 2,626 reported structure problems, mostly legacy
CDS transcripts without exons and reused transcript IDs, plus 47 anellovirus
annotation gaps. Its reverse catalogue check has no unlisted GTF accessions.
The four stale gene rows for the excluded `AB303562.1` duplicate have been
removed; the canonical `NC_038359.1` annotations remain. (2,628 became 2,629 on
2026-10-10: Human astrovirus `NC_001943.1` wrote `Non structural gene` as its
column-3 feature type, so two CDS rows were invisible to this check; as `CDS`
they now show the same exon-less shape as the other legacy records. 2,629 became 2,626 when the duplicate `NC_002076.2` GTF
`Torque_teno_virus_NC_002076.gtf` was removed in favour of `Alphatorquevirus homin1.gtf`.)

Only an explicit `gene_biotype` or `transcript_biotype` of `whole_genome` labels a
feature as such; gene names and long spans do not classify it. As of 2026-10-10
the local corpus tags no gene as `whole_genome`: HPV-8 `M12737.1` and HPV-71
`NC_039089.1` were given modelled CDS (`scripts/model_hpv_cds_gtfs.py`, QC in
`analysis/panel_expansion/hpv_modelled_cds_qc.tsv`), and the genome-spanning
`misc_RNA` / `prim_transcript` of HAV `NC_001489.1` and HTLV-2 `NC_001488.1` were
removed (`scripts/fix_thin_gtfs.py`), keeping their real CDS. `whole_genome`
genes still arise for FASTA-only user references; excluding them from
gene-level outputs and programmes (`DEF-04`) is deferred because the
whole-genome gene is such a virus's only unit. This check does not silently change counts or gene membership.

The EC budget helper accepts a normalized JSON object with integer
`max_ec_size` and `discarded_ec_count` values:

```bash
PYTHONPATH=src python scripts/panel_integrity.py --output /tmp/panel_integrity.tsv \
  --ec-summary /tmp/normalized_inspect.json --max-ec-size 20 --max-discarded-ecs 0
```

The numbers above illustrate explicit budgets, not calibrated defaults. The
helper rejects missing/invalid metrics and exceeded budgets. It does not run
or parse `kallisto inspect`; CAT-40 needs the real index experiment and an
audited normalization of that result. `--max-ec-size`/`--max-discarded-ecs`
without `--ec-summary` is a usage error (the budgets would otherwise be ignored),
and a catalogue with no `shipped` rows or an unknown `panel` value is rejected
rather than checked vacuously.

`viral_summary.tsv` and the HTML report share `reference_risk_flags` from the
run Virus Identity table. Endogenous and contaminant roles produce
`endogenous_reference` and `vector_reagent_reference`; the curated `HHV-6`
sibling group carries `iciHHV6_possible`. An EVE-risk annotation additionally
carries `endogenous_overlap_possible`. These labels require interpretation;
they do not identify the source, diagnose integration, suppress counts or
establish infection. Automatic HERV exclusion and a calibrated CAR-T/vector
filter remain separate unfinished decisions.

SV40 (`NC_001669.1`) and AAV1/7/8 (`NC_002077.1`, `NC_006260.1`, `NC_006261.1`) were fetched from NCBI,
catalogued as `shipped` + `contaminant` and given generated GTFs on 2026-10-10 (PLAN `PANEL-01` WP2). The GTFs
stay local and are listed in `analysis/panel_expansion/gtf_corpus_manifest.tsv`.
`analysis/panel_expansion/promotions_pending_records.tsv` is the overlay that recorded the decisions;
`scripts/panel_promote.py ... --check` passes now that all four accessions are in the catalogue. The catalogue
generator writes the requested set rather than merging it, so never run it against the live catalogue without
`--out`.

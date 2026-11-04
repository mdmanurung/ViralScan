---
orphan: true
---

# Reference source licences (DRAFT)

**Status: draft, needs maintainer/legal review.** Nothing here is legal advice. Every row
below was written from the source's published terms and has not been checked by a
maintainer. The manifest keeps `source_licence.status: unreviewed` and
`review_required: true` until someone signs this table off (PLAN `REF-05`).

The machine-readable draft text in each manifest row comes from
`src/viralscan/source_licences.py`. Keep it in step with this table.

## What ViralScan fetches or redistributes

| Source | Used for | Published terms (summary) | Where to check | Redistributed by ViralScan? | Open question for the reviewer |
|---|---|---|---|---|---|
| NCBI GenBank / RefSeq nucleotide records (`NCBI nuccore`) | Viral genome FASTA and CDS features, fetched at build time | NCBI places no restrictions on use of its database content. A submitter can still claim patent or copyright on a record. INSDC policy is no restrictions on use or distribution. | <https://www.ncbi.nlm.nih.gov/home/about/policies/>, <https://www.insdc.org/policy/> | The packaged GTFs are ViralScan-built annotations derived from these records. The archive of GTFs is planned for Zenodo (`REF-11`). Sequences are fetched, not shipped. | Do derived GTFs need a per-record claim check, or is the INSDC statement enough? Which licence do we put on the GTF archive? |
| NCBI Taxonomy (`virus_taxonomy.tsv`, retrieved 2026-10-16) | Species, genus, family lookups | Same NCBI policy as above. | <https://www.ncbi.nlm.nih.gov/home/about/policies/> | Yes, as a packaged table. | Does the table need an NCBI attribution line? |
| Ensembl release 116 human cDNA and GTF (`Ensembl cDNA`) | Competitive host transcriptome in the combined index | Ensembl data are freely available without restriction under the EMBL-EBI terms of use. Cite the Ensembl release. | <https://www.ensembl.org/info/about/legal/disclaimer.html>, <https://www.ebi.ac.uk/about/terms-of-use> | No. Fetched at build time. The index holds its k-mers. | Is an index built from Ensembl cDNA a redistribution of Ensembl data? |
| 10x Genomics refdata-gex-GRCh38-2024-A `genome.fa` | Genome D-list for cat42d (`REF-06`) | Not found in a form we could quote. The sequence is GRCh38 from the Genome Reference Consortium and the annotation is GENCODE/Ensembl. The 10x bundle terms were not checked. | <https://www.10xgenomics.com/support/software/cell-ranger/downloads> | No. Used as a local file. The cat42d index holds D-list k-mers only. | What are the 10x bundle terms, and do the GRC and GENCODE terms cover the genome FASTA? |
| Clareau lab `human_anellovirus_pangenome` (commit `3ed77e1`) | Accession list for the anellovirus panel | The repository publishes no licence, paper or DOI. | <https://github.com/clareaulab/human_anellovirus_pangenome> | No. Only accession identifiers are used. Sequences and CDS are re-fetched from NCBI. | Is an accession list enough of a derivation to need permission? Ask the authors? |
| igv.js 2.15.11 (`src/viralscan/vendor/igv/`) | Offline evidence browser | MIT. Licence file is vendored and shipped. | `src/viralscan/vendor/igv/LICENSE` | Yes, with its licence text. | None. |
| User-supplied FASTA/GTF (`local_input`, `viralscan --reference`) | Custom indexes | Unknown to ViralScan. The user owns the check. | n/a | No. | None. |

## What the manifest records today

- NCBI rows: `source: "NCBI nuccore"`, terms text above, `status: unreviewed`.
- Ensembl host rows: `source: "Ensembl cDNA"`, terms text above, `status: unreviewed`.
- Local input rows: `source: "local_input"`, terms text says the user owns the check.
- The 10x D-list is recorded in the manifest `genome_dlist` block (path and SHA-256), not as a
  sequence row, so it carries no licence field. Whether it should is a reviewer decision.

## Not yet assessed

- Third-party data inside the Zenodo panel archive (`REF-11`) beyond the rows above.
- Any licence for the ViralScan-authored annotations themselves. The maintainer has to pick
  one (the package code is MIT). This is a decision, not something the code can settle.

## Sign-off

Reviewer, date, and the commit of this file at sign-off go here. Until then the manifest
fields stay `unreviewed`.

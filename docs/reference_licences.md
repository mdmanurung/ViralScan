---
orphan: true
---

# Reference source licences (working position)

**Status: working position (user decision 2026-11-04). One item remains: a legal glance before the
Zenodo deposit, so the manifest keeps `source_licence.status: unreviewed` and `review_required: true`
until it is done (PLAN `REF-05`); the table still needs maintainer/legal review.** Nothing here is
legal advice. Every row was written from the source's published terms.

The machine-readable draft text in each manifest row comes from
`src/viralscan/source_licences.py`. Keep it in step with this table.

## What ViralScan fetches or redistributes

| Source | Used for | Published terms (summary) | Where to check | Redistributed by ViralScan? | Working position |
|---|---|---|---|---|---|
| NCBI GenBank / RefSeq nucleotide records (`NCBI nuccore`) | Viral genome FASTA and CDS features, fetched at build time | NCBI places no restrictions on use of its database content. A submitter can still claim patent or copyright on a record. INSDC policy is no restrictions on use or distribution. | <https://www.ncbi.nlm.nih.gov/home/about/policies/>, <https://www.insdc.org/policy/> | The packaged GTFs are ViralScan-built annotations derived from these records. The archive of GTFs is planned for Zenodo (`REF-11`). Sequences are fetched, not shipped. | The ViralScan GTF archive is released under CC BY 4.0. Derived GTFs get no per-record claim check; the archive carries an NCBI attribution and a disclaimer that submitters may hold rights in individual records. |
| NCBI Taxonomy (`virus_taxonomy.tsv`, retrieved 2026-10-16) | Species, genus, family lookups | Same NCBI policy as above. | <https://www.ncbi.nlm.nih.gov/home/about/policies/> | Yes, as a packaged table. | The table carries the same NCBI attribution and disclaimer as the GTF archive. |
| Ensembl release 116 human cDNA and GTF (`Ensembl cDNA`) | Competitive host transcriptome in the combined index | Ensembl data are freely available without restriction under the EMBL-EBI terms of use. Cite the Ensembl release. | <https://www.ensembl.org/info/about/legal/disclaimer.html>, <https://www.ebi.ac.uk/about/terms-of-use> | No. Fetched at build time. The index holds its k-mers. | An Ensembl-derived index is redistributable, with the Ensembl release cited. |
| 10x Genomics refdata-gex-GRCh38-2024-A `genome.fa` | Genome D-list for cat42d (`REF-06`) | Not found in a form we could quote. The sequence is GRCh38 from the Genome Reference Consortium and the annotation is GENCODE/Ensembl. The 10x bundle terms were not checked. | <https://www.10xgenomics.com/support/software/cell-ranger/downloads> | No. Used as a local file. The cat42d index holds D-list k-mers only. | The 10x GRCh38-2024-A bundle is never redistributed. ViralScan uses it as a local file only. |
| Clareau lab `human_anellovirus_pangenome` (commit `3ed77e1`) | Accession list for the anellovirus panel | The repository publishes no licence, paper or DOI. | <https://github.com/clareaulab/human_anellovirus_pangenome> | No. Only accession identifiers are used. Sequences and CDS are re-fetched from NCBI. | Accession identifiers only; sequences are re-fetched from NCBI. Credit the repository. The pangenome release question is `ANELLO-15`. |
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
- Licence of the ViralScan-authored annotations: CC BY 4.0 for the GTF archive (working position above); the package code stays MIT.

## Sign-off

Remaining item: a legal glance before the Zenodo deposit. Reviewer, date, and the commit of this
file at sign-off go here. Until then the manifest fields stay `unreviewed`.

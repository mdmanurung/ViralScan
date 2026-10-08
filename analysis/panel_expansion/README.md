# Panel expansion candidates (PLAN `PANEL-01`)

Review input for promoting human-relevant viruses into the shipped panel. Nothing here changes the panel.

| File | What it is |
|---|---|
| `evonk_candidates.tsv` | 297 RefSeq accessions an older colleague folder annotated that the catalogue does not hold (221 genomes + 76 HPV RefSeq aliases). Accession list only: no GTF or FASTA is taken from that folder. |
| `human_relevant_curated.tsv` | Curated human-relevance rows (`status` = `proposed` or `accepted`). Only `accepted` rows count; `UNCERTAIN` in `basis` marks ones to confirm. |
| `candidates.tsv` | Generated review table (see `scripts/panel_candidates.py` for the column meanings). |

Regenerate (the sequence-twin check needs the two FASTAs and the built panel FASTA; without them
`twin_checked` is `no` and no row is excluded as a twin):

    python scripts/panel_candidates.py \
      --panel-fasta <build>/viral.fa \
      --candidate-fasta <colleague>/virtus2_reference/filtered_new_viruses.fasta \
      --candidate-fasta <colleague>/fasta_viruses/viruses.fasta

A row is promoted only after the user accepts it; `relevance` says why it may belong, never that it should.

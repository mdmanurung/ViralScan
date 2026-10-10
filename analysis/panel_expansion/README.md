# Panel expansion candidates (PLAN `PANEL-01`)

Review input for promoting human-relevant viruses into the shipped panel. Nothing here changes the panel.

| File | What it is |
|---|---|
| `evonk_candidates.tsv` | 297 RefSeq accessions an older colleague folder annotated that the catalogue does not hold (221 genomes + 76 HPV RefSeq aliases). Accession list only: no GTF or FASTA is taken from that folder. |
| `human_relevant_curated.tsv` | Curated human-relevance rows (`status` = `proposed` or `accepted`). Only `accepted` rows count; `UNCERTAIN` in `basis` marks ones to confirm. |
| `census.tsv` | NCBI nuccore census (2026-10-08): 401 RefSeq viral records whose `/host` is Homo sapiens; 123 are not in the catalogue. Raw answer cached in `census_raw.json`; rerun offline, `--refresh` to query again. The host qualifier is submitter-provided (the list includes phage and bacterial records), so it is a lower bound and a review list, not a verdict. |
| `candidates.tsv` | Generated review table (see `scripts/panel_candidates.py` for the column meanings). |
| `kmer_sharing.tsv` | WP1b: per candidate, the fraction of its k=31 k-mers that are unique, outside the panel, shared with its group, shared with another group (`scripts/panel_kmer_sharing.py`; columns in its docstring). Candidate FASTAs are cached in `~/.cache/viralscan/panel_candidates/`. |
| `pool_vhdb.tsv` | WP2: Virus-Host DB viruses with a human host (taxid 9606), marked `panel`/`candidate`/`excluded`/`gap` by accession, with a review-order `tier` (`scripts/panel_pool.py --vhdb <virushostdb.tsv>`; the database table is not committed). |
| `kmer_partners.tsv` | WP1b: top 3 sharing partners of every candidate with under half its k-mers outside the panel, with `length_ratio` and `frac_of_partner` to tell a twin from a fragment. |
| `kmer_panel_baseline.tsv` | WP1b, `--baseline`: each panel genome scored against the rest of the panel, the sharing the shipped panel already accepts. |
| `cat28_panviral_marker_check.tsv` | CAT-28: each `gene_programs.tsv` marker row against the panviral annotation table, `matched` / `unmatched` / `accession_mismatch` with a `detail` (`extras/cat28_panviral_marker_check.py`; the vendored input is gitignored). Validation only; the catalogue is not changed. |
| `roles_round1.tsv` | WP2: accepted role overlays for existing AAV references (`contaminant`) and HERV-K113 (`endogenous`), applied with `python scripts/panel_promote.py analysis/panel_expansion/roles_round1.tsv`. Existing MLV/XMRV decoys keep `decoy`; SV40 and AAV1/7/8 remain candidates pending NCBI records and promotion. Roles label interpretation and do not filter counts. |

Regenerate (the sequence-twin check needs the two FASTAs and the built panel FASTA; without them
`twin_checked` is `no` and no row is excluded as a twin):

    python scripts/panel_census.py            # offline from census_raw.json; --refresh queries NCBI
    python scripts/panel_candidates.py \
      --panel-fasta <build>/viral.fa \
      --candidate-fasta <colleague>/virtus2_reference/filtered_new_viruses.fasta \
      --candidate-fasta <colleague>/fasta_viruses/viruses.fasta

A row is promoted only after the user accepts it; `relevance` says why it may belong, never that it should.

Then WP1b and the k-mer exclusions (rows with under 5 % of k-mers outside the panel become `kmer_twin_of:<panel id>`;
the three steps are stable under rerun, and the first command above must carry the same `--candidate-fasta` arguments again,
otherwise the sequence-twin exclusions disappear):

    NCBI_EMAIL=<you> python scripts/panel_kmer_sharing.py --fetch --panel-fasta <build>/viral.fa --baseline --vmr <VMR_MSL41.xlsx>   # numpy; VMR from https://ictv.global/vmr (not committed)
    python scripts/panel_candidates.py --panel-fasta <build>/viral.fa --candidate-fasta ... --candidate-fasta ...   # as above

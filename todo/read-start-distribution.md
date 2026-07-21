# Read-start distribution + PCR-duplicate handling in `viralscan evidence`

**Status:** planned (2026-07-21) · **Priority:** medium · **Effort:** ~1–2 days · **Owner:** unassigned

## Motivation

Profile the *positional* distribution of viral reads along the viral genome, as in
Chen et al., *Nature* 2024 (s41586-024-07575-x, Fig. 14): "we removed PCR duplicates
with picard MarkDuplicates and tallied the location within the SARS-CoV-2 genome
using the start of each sequencing read."

A per-position read-start histogram exposes structure that the current per-reference
coverage summary (breadth / mean depth / #reads) cannot:

- **3′ bias** of 10x 3′ chemistry (start pile-up near the 3′ end).
- **Subgenomic-RNA junctions** for coronaviruses (leader–body junction pile-ups) —
  a genuine-transcription signature.
- **EVE / integration hotspots** — artifact reads concentrate at 1–2 loci
  (complements `accession_breadth` and the EVE pipeline's `samtools depth`).
- **Coverage evenness** — genuine infection spreads across the genome; contamination
  or mismapping piles at few positions.

## What exists today (do not duplicate)

- `viralscan evidence` (`src/viralscan/scripts/evidence_run.py` → `src/viralscan/evidence.py`)
  extracts the reads behind viral calls (`extract_viral_reads`), aligns them to the
  viral FASTA with minimap2 (`align_reads_to_viral` → BAM+BAI), and scores per-reference
  coverage (`coverage_table` = `samtools coverage`) and per-read BLAST identity.
- The COVID EVE pipeline (`covid_viralscan/scripts/slurm_eve_analysis.sh`) already runs
  `samtools depth -a` for per-base **depth** (`depth/per_base_*.tsv`) — full-length
  coverage, not read-start, and only in that ad-hoc script.

**Gaps vs the paper:** (a) no PCR-duplicate handling on the evidence BAM
(`extract_viral_reads` keeps *every* read per `(CB,UMI)`); (b) no read-**start**
tally — coverage/depth count every base a read spans, not the 5′ leftmost position.

## Design

Add to `src/viralscan/evidence.py`:

```python
def read_start_distribution(
    bam: str,
    *,
    dedup: str = "umi",          # "umi" | "markdup" | "none"
    strand_aware: bool = True,   # tally the 5' end (POS on +, end on -)
    bin_size: int = 1,
) -> list[dict]:                 # rows: {reference, position, n_read_starts, n_reads}
    ...
```

**Dedup modes:**
- `umi` (default; correct for droplet scRNA-seq): collapse to one alignment per
  `(CB, UMI, reference)`. The extracted reads are named `<CB>_<UMI>_<n>`, so parse the
  read name and keep the first alignment per `(CB, UMI, ref)`. This is the
  scRNA-appropriate dedup and matches how the count matrix is built (bustools UMI
  collapse) — unlike picard's position-based marking.
- `markdup` (to reproduce the paper literally): shell out to `samtools markdup`
  (preferred, no Java dep) or picard `MarkDuplicates`, then read `-F 0x400`.
- `none`: raw reads (current behaviour).

**Tally:** for each surviving primary alignment (`-F 0x904`), take the 5′ start —
leftmost `POS` for `+` strand, `POS + reference_length` for `−` strand when
`strand_aware` — bin by `bin_size`, count per reference. Prefer `pysam`
(`aln.reference_start`, `aln.is_reverse`, `aln.reference_length`) if importable; else
`samtools view -F 0x904` + awk on cols 2 (FLAG), 3 (RNAME), 4 (POS) with CIGAR for
the reverse-strand end.

**Outputs** (under the evidence output dir):
- `read_start_profile.tsv` — `reference, position, n_read_starts, n_reads_dedup`.
- optional `plots/read_start_<virus>.png` — position (x) vs read-start count (y),
  reusing the existing seaborn/matplotlib plotting style.

**CLI** (`evidence_run.py`): `--read-start-profile` flag (+ `--dedup {umi,markdup,none}`,
`--bin-size N`). Runs after the BAM is built; requires the same `--viral-fasta`/minimap2
path already used for coverage.

## Tests

- Synthetic SAM/BAM (built with pysam or a hand-written SAM) with known read starts and
  deliberate duplicate `(CB,UMI)` reads → assert (a) the per-position histogram, (b)
  `dedup="umi"` collapses duplicates to one, (c) `dedup="none"` keeps them, (d)
  strand-aware 5′-end tally for a reverse read.
- Keep the pure tally/dedup logic in a small unit-testable function (mirror the existing
  `_parse_coverage_output` / `_parse_blast_output` split so no live samtools is needed).

## Caveats to document

- 10x 3′ vs 5′ chemistry changes where the start pile-up sits; document that the profile
  reflects capture chemistry, not only biology.
- Strandedness: a true 5′-start profile needs strand handling; `--strand-aware` default on,
  with a note that minimap2 `-ax sr` orientation must be trusted.
- UMI dedup ≠ picard MarkDuplicates; they answer different questions (molecule-level vs
  position-level). Report which was used in the TSV header.

## Related

- Enables a downstream **subgenomic-RNA junction** detector (coronaviruses) — see ROADMAP.
- Complements `accession_breadth` and the EVE `samtools depth` loci analysis.
- See `.living/decisions.md` 2026-07-21 "Viral-read positional profiling" for the assessment.

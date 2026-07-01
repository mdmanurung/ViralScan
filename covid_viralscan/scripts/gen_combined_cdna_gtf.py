#!/usr/bin/env python3
"""Generate combined_cdna.gtf for the covid_viralscan kb ref repair.

Root cause of the hung build (job 25137178): the original combined.gtf used the
Ensembl *chromosomal* GTF (seqnames = 1, 2, X, …) which does not match the cDNA
FASTA headers (ENST transcript IDs).  ngs_tools's genome-splitting step could not
find any human sequence and hung indefinitely.

Fix: generate a cDNA-level GTF where seqname = ENST transcript ID (matching the
FASTA header), spanning coordinates 1–length (the full transcript).  Viral GTFs
(viral_whole_genome.gtf, viral_genome.gtf) are appended unchanged but with the
same Step-2.5 sanitization that slurm_build_ref.sh applied.

Usage::

    python gen_combined_cdna_gtf.py \\
        viralscan_ref/combined.fa \\
        viralscan_ref/viral/viral_whole_genome.gtf \\
        references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf \\
        viralscan_ref/combined_cdna.gtf
"""
from __future__ import annotations

import sys
from pathlib import Path


def _parse_gene_id(header_parts: list[str]) -> str:
    """Extract ENSG gene ID from Ensembl cDNA FASTA header.

    Header example:
        >ENST00000632585.1 cdna scaffold:… gene:ENSG00000282172.1 gene_biotype:…
    Returns "ENSG00000282172.1"; falls back to the transcript ID.
    """
    for part in header_parts:
        if part.startswith("gene:"):
            return part[5:]
    return header_parts[0]


def _write_gtf_record(out, seqname: str, gene_id: str, tx_id: str, length: int) -> None:
    attrs = f'gene_id "{gene_id}"; transcript_id "{tx_id}";'
    for feature in ("gene", "transcript", "exon"):
        out.write(f"{seqname}\tEnsembl_cDNA\t{feature}\t1\t{length}\t.\t+\t.\t{attrs}\n")


def _append_viral_gtf(out, path: Path) -> int:
    """Append a viral GTF with the Step-2.5 sanitization rules.

    1. Skip blank lines and GFF3 '###' directives (ngs_tools raises GtfEntryError).
    2. Prefix 'unassigned_transcript_N' IDs with the seqname to make them
       globally unique across all viral GTFs (avoids t2g collisions).
    """
    count = 0
    with open(path) as fh:
        for raw in fh:
            stripped = raw.strip()
            if not stripped or stripped.startswith("###"):
                continue
            if 'transcript_id "unassigned_transcript' in raw:
                seqname = raw.split("\t", 1)[0]
                raw = raw.replace(
                    'transcript_id "unassigned_transcript',
                    f'transcript_id "{seqname}-unassigned_transcript',
                    1,
                )
            out.write(raw)
            count += 1
    return count


def main() -> None:
    if len(sys.argv) < 4:
        sys.exit(
            "Usage: gen_combined_cdna_gtf.py combined.fa viral1.gtf [viral2.gtf …] output.gtf"
        )

    fasta_path = Path(sys.argv[1])
    output_path = Path(sys.argv[-1])
    viral_gtf_paths = [Path(p) for p in sys.argv[2:-1]]

    print(f"FASTA:      {fasta_path}  ({fasta_path.stat().st_size / 1e9:.2f} GB)", file=sys.stderr)
    print(f"Viral GTFs: {[str(p) for p in viral_gtf_paths]}", file=sys.stderr)
    print(f"Output:     {output_path}", file=sys.stderr)

    enst_count = 0
    non_enst_count = 0

    with open(fasta_path) as fasta, open(output_path, "w") as out:
        current_id: str | None = None
        current_gene: str = ""
        current_len: int = 0

        def flush() -> None:
            nonlocal enst_count, non_enst_count, current_id, current_len
            if current_id is None:
                return
            if current_id.startswith("ENST"):
                enst_count += 1
                _write_gtf_record(out, current_id, current_gene, current_id, current_len)
                if enst_count % 500_000 == 0:
                    print(f"  … {enst_count:,} ENST entries written", file=sys.stderr)
            else:
                non_enst_count += 1

        for raw in fasta:
            raw = raw.rstrip()
            if raw.startswith(">"):
                flush()
                parts = raw[1:].split()
                current_id = parts[0]
                current_gene = _parse_gene_id(parts)
                current_len = 0
            else:
                current_len += len(raw)

        flush()  # last record

    print(f"  ENST transcripts written: {enst_count:,}", file=sys.stderr)
    print(f"  Non-ENST sequences skipped (handled by viral GTFs): {non_enst_count}", file=sys.stderr)

    with open(output_path, "a") as out:
        for gtf_path in viral_gtf_paths:
            n = _append_viral_gtf(out, gtf_path)
            print(f"  Appended {n:,} lines from {gtf_path.name}", file=sys.stderr)

    total_lines = sum(1 for _ in open(output_path))
    size_mb = output_path.stat().st_size / 1e6
    print(f"Done. Output: {total_lines:,} lines, {size_mb:.1f} MB", file=sys.stderr)


if __name__ == "__main__":
    main()

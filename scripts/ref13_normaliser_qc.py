#!/usr/bin/env python3
"""REF-13 QC: kb cDNA of every bundled viral GTF before and after the normaliser.

For each GTF the cDNA is generated with the installed ``ngs_tools`` (the code path of
``kb ref``) from the raw file and from :func:`viralscan.gtf_normalise.normalise_gtf_file`,
on the genomes cached by ``ncbi_fetch``. One row per gene: the rule applied, the transcript
count and cDNA length before/after, and whether the gene's transcript sequences are identical
(transcript names are ignored). Read-only: the GTFs and the cache are never written.

    python scripts/ref13_normaliser_qc.py -o analysis/panel_expansion/ref13_normaliser_qc.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "src"))

from viralscan.gtf_normalise import normalise_gtf_file  # noqa: E402
from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR  # noqa: E402

ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
COLUMNS = (
    "gene_id",
    "file",
    "rule",
    "transcripts_before",
    "transcripts_after",
    "length_before",
    "length_after",
    "bp_removed",
    "identical",
)


def cdna(lines: list[str], fasta: Path, work: Path) -> dict[str, dict[str, str]]:
    """``{gene_id: {transcript_id: sequence}}`` exactly as ``kb ref`` would split it."""
    import ngs_tools as ngs

    gtf, out = work / "in.gtf", work / "cdna.fa"
    gtf.write_text("\n".join(lines) + "\n")
    genes, transcripts = ngs.gtf.genes_and_transcripts_from_gtf(str(gtf), use_version=True)
    ngs.fasta.split_genomic_fasta_to_cdna(str(fasta), str(out), genes, transcripts)
    result: dict[str, dict[str, str]] = defaultdict(dict)
    with ngs.fasta.Fasta(str(out), "r") as handle:
        for entry in handle:
            result[entry.attributes["gene_id"]][entry.name] = entry.sequence
    return result


def genome_fasta(seqnames: set[str], cache: Path, work: Path) -> Path:
    path = work / "genomes.fa"
    with path.open("w") as handle:
        for seq in sorted(seqnames):
            text = (cache / seq / f"{seq}.fasta").read_text()
            body = text.split("\n", 1)[1]
            handle.write(f">{seq}\n{body}")
    return path


def gene_rules(lines: list[str]) -> dict[str, str]:
    rules: dict[str, str] = {}
    for line in lines:
        cols = line.split("\t")
        if len(cols) == 9 and cols[2] == "transcript":
            attrs = dict(ATTR.findall(cols[8]))
            rules.setdefault(attrs["gene_id"], attrs.get("viralscan_norm", "passthrough"))
    return rules


def qc_file(path: Path, cache: Path) -> list[dict[str, object]]:
    raw = path.read_text().splitlines()
    new = normalise_gtf_file(path, cache)
    seqnames = {line.split("\t", 1)[0] for line in raw if line and line[0] != "#"}
    rules = gene_rules(new)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        fasta = genome_fasta(seqnames, cache, work)
        before, after = cdna(raw, fasta, work), cdna(new, fasta, work)
    rows = []
    for gene in sorted(set(before) | set(after)):
        b, a = before.get(gene, {}), after.get(gene, {})
        lb, la = sum(map(len, b.values())), sum(map(len, a.values()))
        rows.append(
            {
                "gene_id": gene,
                "file": path.name,
                "rule": rules.get(gene, "passthrough"),
                "transcripts_before": len(b),
                "transcripts_after": len(a),
                "length_before": lb,
                "length_after": la,
                "bp_removed": lb - la,
                "identical": "yes" if sorted(b.values()) == sorted(a.values()) else "no",
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=root / "src" / "viralscan" / "data")
    parser.add_argument("--cache-dir", type=Path, default=Path(DEFAULT_CACHE_DIR))
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)
    rows: list[dict[str, object]] = []
    for gtf in sorted(args.data_dir.glob("*.gtf")):
        rows += qc_file(gtf, args.cache_dir)
    with args.output.open("w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    changed = [r for r in rows if r["identical"] == "no"]
    print(f"{len(rows)} genes: {len(rows) - len(changed)} identical, {len(changed)} changed")
    print("changed by rule:", dict(Counter(str(r["rule"]) for r in changed)))
    print("bp removed (changed genes):", sum(int(str(r["bp_removed"])) for r in changed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
Annotate EVE analysis results with GRCh38 gene information.

Inputs:
  --phase-a-dir:   directory with {ACC}_{sample}_grch38.depth.txt files
  --phase-b-blast: blastn hits (covered viral seqs vs NT human)
  --phase-c-paf:   minimap2 PAF (full viral panel vs GRCh38, asm20)
  --gtf:           GRCh38 genes.gtf (CellRanger format)
  --outdir:        output directory

Outputs:
  phase_a_loci.tsv      human loci hit by artifact reads, per accession
  phase_b_blast.tsv     BLAST hits annotated with gene context
  phase_c_panel.tsv     panel-wide EVE screen summary per accession
  eve_summary_report.txt  human-readable summary
"""

import argparse
import collections
import gzip
import os
import re
import sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--phase-a-dir", required=True)
    p.add_argument("--phase-b-blast", required=True)
    p.add_argument("--phase-c-paf", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--key-accs", required=True)
    return p.parse_args()


# ── GTF interval lookup (no bedtools needed) ──────────────────────────────────


def load_gtf_genes(gtf_path):
    """
    Returns dict: chrom → sorted list of (start, end, gene_name, gene_type)
    Parses gene-level records only (not transcript/exon).
    """
    genes = collections.defaultdict(list)
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, start, end = parts[0], int(parts[3]), int(parts[4])
            attrs = parts[8]
            gene_name = _attr(attrs, "gene_name") or _attr(attrs, "gene_id") or "?"
            gene_type = _attr(attrs, "gene_type") or _attr(attrs, "gene_biotype") or "?"
            genes[_normalize_chrom(chrom)].append((start, end, gene_name, gene_type))
    # Keep deterministic output and nearest-gene tie handling.
    for chrom in genes:
        genes[chrom].sort()
    return genes


def _attr(attrs, key):
    for token in attrs.split(";"):
        token = token.strip()
        if token.startswith(key + " ") or token.startswith(key + "\t"):
            val = token.split(None, 1)[1].strip('"')
            return val
    return None


def _normalize_chrom(chrom):
    """Canonicalize a chromosome name so 'chr7'/'7' and 'chrM'/'MT'/'M' match.

    UCSC-style ('chr7') and Ensembl-style ('7') references otherwise never match
    the GTF keys, which silently annotates every locus as 'intergenic'. Strips a
    leading 'chr' (any case) and folds mitochondrial aliases to 'MT'.
    """
    if chrom is None:
        return chrom
    c = chrom.strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c.upper() in ("M", "MT", "MTDNA"):
        return "MT"
    return c


def _is_chromosome_subject(sseqid):
    """True if a BLAST subject is a whole human chromosome (coords are genomic).

    Human chromosome RefSeq accessions are NC_000001..NC_000024; mito is
    NC_012920. Clones/scaffolds (AC_*, NT_*, NW_*, ...) carry subject-LOCAL
    coordinates that must NOT be looked up against chromosome-level GTF intervals.
    Accepts bare accessions (modern ``-outfmt 6``) and legacy ``gi|...|ref|NC_...|``.
    """
    # Human chromosomes are NC_000001..NC_000024 (+ mito NC_012920) — not the
    # broader NC_0000\d\d, which would also accept e.g. mouse NC_000067.
    return bool(re.search(r"(?:^|\|)(NC_0000(?:0[1-9]|1\d|2[0-4])|NC_012920)", sseqid or ""))


def _warn_namespace(queried, genes_by_chrom, phase):
    """Warn if no queried chromosome matches the GTF namespace (silent-fail guard)."""
    known = set(genes_by_chrom)
    if queried and not (queried & known):
        sys.stderr.write(
            f"WARNING [{phase}]: none of the {len(queried)} query chromosomes match "
            f"the GTF namespace (query e.g. {sorted(queried)[:3]} vs GTF e.g. "
            f"{sorted(known)[:3]}). All gene annotations will be 'intergenic' — check "
            f"that the alignment reference and GTF share one assembly and naming.\n"
        )


def annotate_locus(chrom, pos, genes_by_chrom, window=50000):
    """Return nearest gene within window bp, or 'intergenic'."""
    chrom = _normalize_chrom(chrom)
    if chrom not in genes_by_chrom:
        return "intergenic", "?", -1
    records = genes_by_chrom[chrom]
    best = ("intergenic", "?", 10**9)
    for s, e, gname, gtype in records:
        if s <= pos <= e:
            return f"in_gene:{gname}", gtype, 0
        dist = min(abs(pos - s), abs(pos - e))
        if dist < best[2]:
            best = (f"near:{gname}", gtype, dist)
    if best[2] > window:
        return "intergenic", "?", best[2]
    return best


def merge_position_depths(positions, max_gap=5):
    """Merge covered 1-based positions into local loci while preserving depth."""
    if not positions:
        return []
    sorted_positions = sorted(positions)
    loci = []
    start = end = sorted_positions[0][0]
    max_depth = sorted_positions[0][1]
    n_covered = 1
    for pos, depth in sorted_positions[1:]:
        if pos <= end + max_gap:
            end = pos
            max_depth = max(max_depth, depth)
            n_covered += 1
        else:
            loci.append((start, end, n_covered, max_depth))
            start = end = pos
            max_depth = depth
            n_covered = 1
    loci.append((start, end, n_covered, max_depth))
    return loci


# ── Phase A annotation ────────────────────────────────────────────────────────


def annotate_phase_a(phase_a_dir, key_accs, genes_by_chrom, outdir):
    key_accs = set(key_accs)
    rows = []
    queried_chroms = set()
    for fname in sorted(os.listdir(phase_a_dir)):
        if not fname.endswith("_grch38.depth.txt"):
            continue
        acc_sample = fname.replace("_grch38.depth.txt", "")
        parts = acc_sample.rsplit("_", 1)
        if len(parts) != 2:
            continue
        acc, sample = parts[0], parts[1]
        if key_accs and acc not in key_accs:
            continue
        fpath = os.path.join(phase_a_dir, fname)
        # Count reads per chromosome from idxstats file
        idxstats_path = fpath.replace(".depth.txt", ".idxstats.txt")
        chrom_reads = {}
        if os.path.exists(idxstats_path):
            with open(idxstats_path) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 3 and int(p[2]) > 0:
                        chrom_reads[p[0]] = int(p[2])
        # Collect high-coverage positions
        loci = collections.defaultdict(list)
        try:
            with open(fpath) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 3:
                        continue
                    chrom, pos, depth = p[0], int(p[1]), int(p[2])
                    if depth >= 2:
                        loci[chrom].append((pos, depth))
        except FileNotFoundError:
            continue
        for chrom, positions in loci.items():
            if chrom == "*":
                continue
            queried_chroms.add(_normalize_chrom(chrom))
            for start, end, n_bases, max_depth in merge_position_depths(positions):
                gene_ann, gene_type, dist = annotate_locus(
                    chrom, (start + end) // 2, genes_by_chrom
                )
                rows.append(
                    {
                        "accession": acc,
                        "sample": sample,
                        "human_chrom": chrom,
                        "human_start": start,
                        "human_end": end,
                        "n_covered_bases": n_bases,
                        "max_depth": max_depth,
                        "chrom_reads": chrom_reads.get(chrom, "?"),
                        "gene_annotation": gene_ann,
                        "gene_type": gene_type,
                        "dist_to_gene": dist,
                    }
                )
    _warn_namespace(queried_chroms, genes_by_chrom, "Phase A")
    rows.sort(
        key=lambda r: (
            -r["max_depth"],
            r["accession"],
            r["sample"],
            r["human_chrom"],
            r["human_start"],
        )
    )
    out_path = os.path.join(outdir, "phase_a_loci.tsv")
    with open(out_path, "w") as out:
        cols = [
            "accession",
            "sample",
            "human_chrom",
            "human_start",
            "human_end",
            "n_covered_bases",
            "max_depth",
            "chrom_reads",
            "gene_annotation",
            "gene_type",
            "dist_to_gene",
        ]
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"Phase A: {len(rows)} human loci written to {out_path}")
    return rows


# ── Phase B annotation ────────────────────────────────────────────────────────


def annotate_phase_b(blast_path, genes_by_chrom, outdir):
    """Annotate NT BLAST hits (human taxid) with gene context."""
    rows = []
    cols_in = [
        "qseqid",
        "sseqid",
        "stitle",
        "pident",
        "length",
        "qlen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    out_path = os.path.join(outdir, "phase_b_blast.tsv")
    out_cols = cols_in + ["human_chrom_parsed", "gene_annotation", "gene_type", "dist_to_gene"]
    if not os.path.exists(blast_path) or os.path.getsize(blast_path) == 0:
        with open(out_path, "w") as out:
            out.write("\t".join(out_cols) + "\n")
        print(f"Phase B: no BLAST hits; wrote header to {out_path}")
        return rows
    with open(blast_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(cols_in):
                continue
            row = dict(zip(cols_in, p))
            # Parse chromosome and position from stitle (NT accession)
            # stitle format: "Homo sapiens chromosome X, GRCh38.p14 ..." or similar
            stitle = row["stitle"]
            chrom = _parse_chrom_from_stitle(stitle)
            if _is_chromosome_subject(row["sseqid"]):
                sstart, send = int(row["sstart"]), int(row["send"])
                mid = (sstart + send) // 2
                gene_ann, gene_type, dist = annotate_locus(chrom, mid, genes_by_chrom)
            else:
                # Clone/scaffold subject: sstart/send are subject-local, not genomic,
                # so a chromosome-interval lookup would report an unrelated gene.
                gene_ann, gene_type, dist = "subject_not_chromosome", "?", -1
            row["human_chrom_parsed"] = chrom
            row["gene_annotation"] = gene_ann
            row["gene_type"] = gene_type
            row["dist_to_gene"] = dist
            rows.append(row)
    rows.sort(key=lambda r: float(r["evalue"]))
    with open(out_path, "w") as out:
        out.write("\t".join(out_cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in out_cols) + "\n")
    print(f"Phase B: {len(rows)} annotated BLAST hits written to {out_path}")
    return rows


def _parse_chrom_from_stitle(stitle):
    """Attempt to extract chromosome name from NT sequence title."""
    # CellRanger-style: "chr1" or "chrX"
    m = re.search(r"\bchr[\dXYMT]+\b", stitle)
    if m:
        return m.group(0)
    # NCBI-style: "chromosome 1", "chromosome X"
    m = re.search(r"chromosome\s+(\w+)", stitle, re.IGNORECASE)
    if m:
        c = m.group(1)
        return f"chr{c}" if not c.startswith("chr") else c
    return "unknown"


# ── Phase C annotation ────────────────────────────────────────────────────────


def annotate_phase_c(paf_path, genes_by_chrom, outdir, min_mapq=10, min_aln_len=50):
    """Summarise minimap2 PAF: which viral accessions hit GRCh38, and where?"""
    acc_hits = collections.defaultdict(list)
    out_path = os.path.join(outdir, "phase_c_panel.tsv")
    cols = [
        "accession",
        "n_alignments",
        "total_aln_bases",
        "best_aln_len",
        "best_pident",
        "best_chrom",
        "best_coords",
        "best_gene_annotation",
        "best_gene_type",
    ]
    if not os.path.exists(paf_path):
        with open(out_path, "w") as out:
            out.write("\t".join(cols) + "\n")
        print(f"Phase C: PAF file not found; wrote header to {out_path}")
        return []
    queried_chroms = set()
    with open(paf_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            q_name = p[0]
            aln_len = int(p[10])
            mapq = int(p[11])
            if aln_len < min_aln_len or mapq < min_mapq:
                continue
            t_name = p[5]
            t_start = int(p[7]) + 1
            t_end = int(p[8])
            queried_chroms.add(_normalize_chrom(t_name))
            pident = int(p[9]) / aln_len * 100 if aln_len > 0 else 0
            gene_ann, gene_type, dist = annotate_locus(
                t_name, (t_start + t_end) // 2, genes_by_chrom
            )
            acc_hits[q_name].append(
                {
                    "human_chrom": t_name,
                    "human_start": t_start,
                    "human_end": t_end,
                    "aln_len": aln_len,
                    "mapq": mapq,
                    "pident": round(pident, 1),
                    "gene_annotation": gene_ann,
                    "gene_type": gene_type,
                }
            )
    rows = []
    for acc, hits in sorted(acc_hits.items()):
        total_aln = sum(h["aln_len"] for h in hits)
        best = max(hits, key=lambda h: h["aln_len"])
        rows.append(
            {
                "accession": acc,
                "n_alignments": len(hits),
                "total_aln_bases": total_aln,
                "best_aln_len": best["aln_len"],
                "best_pident": best["pident"],
                "best_chrom": best["human_chrom"],
                "best_coords": f"{best['human_start']}-{best['human_end']}",
                "best_gene_annotation": best["gene_annotation"],
                "best_gene_type": best["gene_type"],
            }
        )
    _warn_namespace(queried_chroms, genes_by_chrom, "Phase C")
    rows.sort(key=lambda r: -r["total_aln_bases"])
    with open(out_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"Phase C: {len(rows)} accessions with GRCh38 hits written to {out_path}")
    return rows


# ── Summary report ────────────────────────────────────────────────────────────


def write_summary(phase_a, phase_b, phase_c, outdir):
    out_path = os.path.join(outdir, "eve_summary_report.txt")
    with open(out_path, "w") as out:
        out.write("EVE Analysis Summary Report\n")
        out.write("=" * 60 + "\n\n")

        out.write("== Phase A: Human loci hit by artifact reads ==\n")
        if phase_a:
            # Group by accession
            by_acc = collections.defaultdict(list)
            for r in phase_a:
                by_acc[r["accession"]].append(r)
            for acc, rows in sorted(by_acc.items()):
                chroms = set(r["human_chrom"] for r in rows)
                max_d = max(r["max_depth"] for r in rows)
                genes = set(r["gene_annotation"] for r in rows)
                out.write(
                    f"  {acc}: chroms={','.join(sorted(chroms))} "
                    f"max_depth={max_d} gene={','.join(sorted(genes))}\n"
                )
        else:
            out.write("  No Phase A results.\n")

        out.write("\n== Phase B: BLAST hits (NT human-only) ==\n")
        if phase_b:
            seen = set()
            for r in phase_b[:50]:
                key = (r["qseqid"], r["human_chrom_parsed"])
                if key in seen:
                    continue
                seen.add(key)
                out.write(
                    f"  {r['qseqid']} → {r['human_chrom_parsed']} "
                    f"({r['pident']}% id, e={r['evalue']}) "
                    f"gene={r['gene_annotation']}\n"
                )
        else:
            out.write("  No BLAST hits.\n")

        out.write("\n== Phase C: Panel-wide EVE screen (viral panel vs GRCh38) ==\n")
        if phase_c:
            out.write(
                f"  Total accessions with GRCh38 homology (alnlen≥50, mapq≥10): {len(phase_c)}\n\n"
            )
            out.write(
                f"  {'Accession':<20} {'N_alns':>6} {'TotalBp':>8} "
                f"{'BestPct':>7} {'Chrom':<8} {'Gene annotation'}\n"
            )
            out.write("  " + "-" * 75 + "\n")
            for r in phase_c[:40]:
                out.write(
                    f"  {r['accession']:<20} {r['n_alignments']:>6} "
                    f"{r['total_aln_bases']:>8} {r['best_pident']:>7} "
                    f"{r['best_chrom']:<8} {r['best_gene_annotation']}\n"
                )
            if len(phase_c) > 40:
                out.write(f"  ... and {len(phase_c) - 40} more (see phase_c_panel.tsv)\n")
        else:
            out.write("  No Phase C results.\n")

    print(f"Summary report written to {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    key_accs = args.key_accs.split()

    print(f"Loading GTF: {args.gtf} ...")
    genes_by_chrom = load_gtf_genes(args.gtf)
    print(
        f"  Loaded {sum(len(v) for v in genes_by_chrom.values())} gene records "
        f"across {len(genes_by_chrom)} chromosomes."
    )

    print("\nAnnotating Phase A (reads → GRCh38 loci)...")
    phase_a = annotate_phase_a(args.phase_a_dir, key_accs, genes_by_chrom, args.outdir)

    print("\nAnnotating Phase B (BLAST hits)...")
    phase_b = annotate_phase_b(args.phase_b_blast, genes_by_chrom, args.outdir)

    print("\nAnnotating Phase C (panel vs GRCh38)...")
    phase_c = annotate_phase_c(args.phase_c_paf, genes_by_chrom, args.outdir)

    print("\nWriting summary report...")
    write_summary(phase_a, phase_b, phase_c or [], args.outdir)


if __name__ == "__main__":
    main()

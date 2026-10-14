"""Independent viral loci (LOCK-01 D1): connected components of overlapping annotated gene intervals.

Two gene intervals are dependent when they share at least one base on either strand (quantification is
unstranded, DEC-04). Rows of one ``gene_id`` are one node, so an origin-wrapping gene (two ``part`` rows at
the two ends of a circular genome) joins whatever overlaps either part: that is the merge across the origin.
A gene with no ``gene`` row takes the span of its CDS/exon rows. Loci are positional only; this module does
not look at counts. ``results/viral_loci.tsv`` is written by the ``analysis`` rule and read by ``detection``.
"""

import gzip
import re
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

LOCI_FILENAME = "viral_loci.tsv"
LOCI_COLUMNS = ["gene_id", "seqname", "start", "end", "locus_id"]
_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_FALLBACK_FEATURES = {"CDS", "exon", "start_codon", "stop_codon"}


def gene_intervals(
    lines: Iterable[str], gene_ids: set[str] | None = None
) -> dict[str, list[tuple[str, int, int]]]:
    """``{gene_id: [(seqname, start, end), ...]}``: ``gene`` rows, else the CDS/exon span of that gene."""
    genes: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    fallback: dict[str, dict[str, list[int]]] = defaultdict(dict)
    for line in lines:
        cols = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(cols) < 9 or cols[2] not in {"gene", *_FALLBACK_FEATURES}:
            continue
        m = _GENE_ID.search(cols[8])
        if not m or (gene_ids is not None and m.group(1) not in gene_ids):
            continue
        gid, start, end = m.group(1), int(cols[3]), int(cols[4])
        if start > end:
            raise ValueError(f"{gid}: start {start} > end {end}")
        if cols[2] == "gene":
            genes[gid].append((cols[0], start, end))
        else:
            span = fallback[gid].setdefault(cols[0], [start, end])
            span[0], span[1] = min(span[0], start), max(span[1], end)
    for gid, spans in fallback.items():
        if gid not in genes:
            genes[gid] = [(seq, s, e) for seq, (s, e) in sorted(spans.items())]
    return dict(genes)


def assign_loci(genes: dict[str, list[tuple[str, int, int]]]) -> dict[str, str]:
    """``{gene_id: locus_id}``; the locus id is the smallest gene_id of its component."""
    parent = {g: g for g in genes}

    def find(g: str) -> str:
        while parent[g] != g:
            parent[g] = parent[parent[g]]
            g = parent[g]
        return g

    by_seq: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for gid, ivs in genes.items():
        for seq, s, e in ivs:
            by_seq[seq].append((s, e, gid))
    for seq_ivs in by_seq.values():
        seq_ivs.sort()
        end, holder = seq_ivs[0][1], seq_ivs[0][2]
        for s, e, gid in seq_ivs[1:]:
            if s <= end:  # closed intervals: sharing one base is an overlap
                parent[find(gid)] = find(holder)
                if e > end:
                    end, holder = e, gid
            else:
                end, holder = e, gid
    members: dict[str, list[str]] = defaultdict(list)
    for gid in genes:
        members[find(gid)].append(gid)
    return {gid: min(group) for group in members.values() for gid in group}


def write_loci(gtf_paths: Iterable[Path], gene_ids: set[str], out: Path) -> Path | None:
    """Write ``out`` (``LOCI_COLUMNS``) for the genes of ``gene_ids`` found in ``gtf_paths``; None if none found.

    Intervals of one gene seen in several files are pooled (identical rows collapse)."""
    genes: dict[str, set[tuple[str, int, int]]] = defaultdict(set)
    for path in gtf_paths:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for gid, ivs in gene_intervals(fh, gene_ids).items():
                genes[gid].update(ivs)
    if not genes:
        return None
    pooled = {gid: sorted(ivs) for gid, ivs in genes.items()}
    locus = assign_loci(pooled)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(LOCI_COLUMNS) + "\n")
        for gid in sorted(pooled):
            for seq, s, e in pooled[gid]:
                fh.write(f"{gid}\t{seq}\t{s}\t{e}\t{locus[gid]}\n")
    return out


def read_loci(results_dir: str | Path) -> dict[str, str]:
    """``{gene_id: locus_id}`` from ``<results_dir>/viral_loci.tsv``; empty when the run wrote none."""
    path = Path(results_dir) / LOCI_FILENAME
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, li = header.index("gene_id"), header.index("locus_id")
        return {c[gi]: c[li] for c in (line.rstrip("\n").split("\t") for line in fh) if len(c) > li}

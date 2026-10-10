"""Build-time GTF normaliser for the bundled viral annotations (PLAN ``REF-13``).

kb-python 0.30.2 (``ngs_tools`` 1.8.6) reads only ``gene``/``transcript``/``exon`` rows and
ignores ``CDS``. Most bundled records are CDS-only, so ``kb ref`` indexes each such gene as one
transcript spanning its gene row (ID = ``gene_id``). That is the right default (the gene span
is what reads from UTRs and spliced-out regions should count towards), but it leaves three
defects this module fixes without touching the shipped GTF files:

* 2,517 CDS-only genes have no ``exon`` row, so ``panel_integrity`` reports them and the
  ``unassigned_transcript_N`` IDs on their CDS rows are shared between genes and files;
* a single-protein spliced CDS keeps its intron in the cDNA (~102 kb over ~102 genes);
* a gene with two gene rows (terminal repeats, circular origin) silently keeps only the last.

:func:`normalise_viral_gtf` is a pure function of one GTF's lines (plus optional flatfile CDS
joins). The rules, in order:

1. A gene that already has ``exon`` rows, or that has no ``gene``/``CDS`` rows at all, passes
   through unchanged. This is also what makes the function idempotent.
2. Otherwise one transcript is emitted per gene-row copy, ID = ``gene_id`` (``-c2``, ``-c3`` ...
   for further copies), and its CDS/start/stop rows are retagged to it. Its exon is the gene
   span, clipped to the gene row: ``kb`` cDNA is byte-identical to today's.
3. Introns are removed only when the gene's CDS rows all belong to one protein: the exons
   are the union of the CDS and stop-codon blocks, the first and last extended to the gene
   row. Gaps shorter than :data:`MIN_INTRON_BP` are programmed frameshifts, not introns.
   Multi-protein genes (adenovirus 3'-coterminal families, HBV S) keep the span: their union
   "gaps" are other proteins' territory, not introns.
4. Genes in :data:`RETAINED_INTRON_GENES` also keep the span transcript (``-span``) beside
   the spliced one, so a retained-intron isoform stays indexable.
5. Repeat copies (two gene rows, no ``part`` attribute: HHV-6B DR/B genes in both terminal
   repeats) never join: one transcript per copy.
6. Origin-wrapping genes (gene rows carrying ``part``: HBV P and S) cannot be represented,
   because ``ngs_tools`` concatenates exons in genomic order. ngs_tools' "last gene row wins"
   behaviour is kept and recorded as ``wrap_last_row``.
7. A CDS-bearing gene with no gene row (B19V 11 kDa) gets one spanning the CDS rows.

Every emitted transcript row carries ``viralscan_norm "<rule>"``.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

#: Gaps between CDS blocks shorter than this are frameshifts (e.g. influenza PA-X), not introns.
MIN_INTRON_BP = 11

#: Genes whose spliced and span-based isoforms are both indexed. HCMV UL111A (``wtgp098``,
#: cmvIL-10) is the latency marker. Its latency isoform LAcmvIL-10 has a different splicing
#: pattern from the lytic transcript (Jenkins et al. 2004, J Virol 78:1440, abstract); that it
#: retains an intron is NOT verified, so the span transcript is kept as a conservative second.
RETAINED_INTRON_GENES = frozenset({"HUM_CYTO_HHV5wtgp098"})

_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
_TX_ATTR = re.compile(r'transcript_id\s+"[^"]*"')
_GENE_ATTR = re.compile(r'gene_id\s+"[^"]*"\s*;?')
_CDS_LIKE = ("CDS", "start_codon", "stop_codon")

Block = tuple[int, int]  # 1-based, inclusive


@dataclass
class _Row:
    index: int
    line: str
    cols: list[str] = field(default_factory=list)
    attrs: dict[str, str] = field(default_factory=dict)

    @property
    def feature(self) -> str:
        return self.cols[2] if self.cols else ""

    @property
    def span(self) -> Block:
        return int(self.cols[3]), int(self.cols[4])


def _parse(index: int, line: str) -> _Row:
    cols = line.split("\t")
    if line.startswith("#") or len(cols) != 9 or not (cols[3].isdigit() and cols[4].isdigit()):
        return _Row(index, line)
    return _Row(index, line, cols, dict(_ATTR.findall(cols[8])))


def _merge(blocks: Iterable[Block], min_gap: int) -> list[Block]:
    """Union of blocks; a gap shorter than ``min_gap`` bp is bridged."""
    merged: list[list[int]] = []
    for start, end in sorted(blocks):
        if merged and start - merged[-1][1] - 1 < min_gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def _retag(row: _Row, transcript_id: str) -> str:
    attrs = row.cols[8]
    if _TX_ATTR.search(attrs):
        attrs = _TX_ATTR.sub(f'transcript_id "{transcript_id}"', attrs, count=1)
    else:
        attrs = _GENE_ATTR.sub(
            lambda m: f'gene_id "{row.attrs["gene_id"]}"; transcript_id "{transcript_id}"; ',
            attrs,
            count=1,
        )
    return "\t".join([*row.cols[:8], attrs.rstrip()])


def _transcript_lines(
    seq: str, strand: str, gene_id: str, transcript_id: str, exons: Sequence[Block], rule: str
) -> list[str]:
    attrs = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; '
    head = f"{seq}\tviralscan\t%s\t%d\t%d\t.\t{strand}\t.\t"
    lines = [head % ("transcript", exons[0][0], exons[-1][1]) + attrs + f'viralscan_norm "{rule}";']
    ordered = exons if strand != "-" else exons[::-1]
    lines += [
        head % ("exon", s, e) + attrs + f'exon_number "{n}";' for n, (s, e) in enumerate(ordered, 1)
    ]
    return lines


def _overlap(block: Block, other: Block) -> int:
    return max(0, min(block[1], other[1]) - max(block[0], other[0]) + 1)


def _rewrite_gene(
    gene_id: str,
    gene_rows: list[_Row],
    cds_rows: list[_Row],
    joins: Mapping[Block, Sequence[Block]],
) -> tuple[list[str], dict[int, str], str | None]:
    """``(new transcript/exon lines, {row index: retagged line}, synthesised gene line)``."""
    anchor = gene_rows[0] if gene_rows else cds_rows[0]
    seq, strand = anchor.cols[0], anchor.cols[6]
    added_gene = None
    rule_prefix = ""
    if not gene_rows:
        spans = [r.span for r in cds_rows]
        lo, hi = min(s for s, _ in spans), max(e for _, e in spans)
        added_gene = (
            f"{seq}\tviralscan\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t"
            f'gene_id "{gene_id}"; gene_biotype "protein_coding"; viralscan_norm "gene_row_added";'
        )
        copies = [(lo, hi)]
        rule_prefix = "gene_row_added+"
    else:
        copies = list(dict.fromkeys(r.span for r in gene_rows))
        if len(gene_rows) > 1 and all("part" in r.attrs for r in gene_rows):
            copies, rule_prefix = [gene_rows[-1].span], "wrap_last_row+"
        elif len(copies) > 1:
            rule_prefix = "repeat_copies+"

    members: dict[int, list[_Row]] = defaultdict(list)
    for row in cds_rows:
        best = max(range(len(copies)), key=lambda k: (_overlap(row.span, copies[k]), -k))
        members[best].append(row)

    lines: list[str] = []
    retag: dict[int, str] = {}
    for k, (gs, ge) in enumerate(copies):
        transcript_id = gene_id if k == 0 else f"{gene_id}-c{k + 1}"
        blocks: list[Block] = []
        proteins: set[str] = set()
        for row in members[k]:
            retag[row.index] = _retag(row, transcript_id)
            if row.feature == "CDS":
                proteins.add(
                    row.attrs.get("protein_id") or row.attrs.get("product") or "unidentified"
                )
                blocks += joins.get(row.span, [row.span])
            elif row.feature == "stop_codon":
                blocks.append(row.span)
        clipped = [(max(s, gs), min(e, ge)) for s, e in blocks if s <= ge and e >= gs]
        merged = _merge(clipped, MIN_INTRON_BP)
        spliced = len(merged) > 1 and len(proteins) == 1 and "unidentified" not in proteins
        rule = rule_prefix + ("spliced" if spliced else "span")
        exons = [(gs, merged[0][1]), *merged[1:-1], (merged[-1][0], ge)] if spliced else [(gs, ge)]
        lines += _transcript_lines(seq, strand, gene_id, transcript_id, exons, rule)
        if spliced and gene_id in RETAINED_INTRON_GENES:
            lines += _transcript_lines(
                seq, strand, gene_id, f"{gene_id}-span", [(gs, ge)], "retained_intron_isoform"
            )
    return lines, retag, added_gene


def _adopt_orphan_cds(by_gene: dict[str, list[_Row]]) -> None:
    """Re-key CDS rows whose gene has no gene row onto the gene-only gene at the same locus.

    B19V's 11 kDa CDS rows say ``HUM_PARVO_unassigned_gene_1`` while the gene row for the same
    4890-5174 span is ``HUM_PARVO_B19V_gp4``, which has no CDS. Giving the orphan its own gene
    row would index the locus twice, as two genes with identical transcripts.
    """
    bare: dict[tuple[str, str, Block], str] = {}
    for gene_id, rows in by_gene.items():
        if {r.feature for r in rows} == {"gene"} and len({r.span for r in rows}) == 1:
            bare[(rows[0].cols[0], rows[0].cols[6], rows[0].span)] = gene_id
    for gene_id, rows in list(by_gene.items()):
        on_one_strand = len({(r.cols[0], r.cols[6]) for r in rows}) == 1
        if {r.feature for r in rows} - set(_CDS_LIKE) or not on_one_strand:
            continue
        span = (min(r.span[0] for r in rows), max(r.span[1] for r in rows))
        target = bare.get((rows[0].cols[0], rows[0].cols[6], span))
        if target:
            for row in rows:
                row.attrs["gene_id"] = target
                row.cols[8] = _GENE_ATTR.sub(f'gene_id "{target}"; ', row.cols[8], count=1)
            by_gene[target] += by_gene.pop(gene_id)


def normalise_viral_gtf(
    lines: Iterable[str], cds_joins: Mapping[str, Iterable[Sequence[Block]]] | None = None
) -> list[str]:
    """Normalise one viral GTF (see the module docstring); returns GTF lines, no newlines.

    ``cds_joins`` maps a seqname to the spliced CDS locations of its flatfile
    (:func:`cds_joins_from_genbank`). A CDS row whose span equals a join's whole span is a
    collapsed join (HHV-6B) and is replaced by the join's blocks.
    """
    rows = [_parse(i, line) for i, line in enumerate(lines)]
    by_gene: dict[str, list[_Row]] = defaultdict(list)
    for row in rows:
        if row.cols and row.attrs.get("gene_id"):
            by_gene[row.attrs["gene_id"]].append(row)

    _adopt_orphan_cds(by_gene)

    insert_after: dict[int, list[str]] = {}
    insert_before: dict[int, list[str]] = {}
    replace: dict[int, str] = {}
    for gene_id, grows in list(by_gene.items()):
        features = {r.feature for r in grows}
        gene_rows = [r for r in grows if r.feature == "gene"]
        cds_rows = [r for r in grows if r.feature in _CDS_LIKE]
        if features & {"exon", "transcript"} or not (gene_rows or cds_rows):
            continue
        if len({(r.cols[0], r.cols[6]) for r in grows}) > 1:
            continue  # rows on two sequences/strands cannot be one gene
        joins = {
            (min(s for s, _ in join), max(e for _, e in join)): sorted(join)
            for join in (cds_joins or {}).get(grows[0].cols[0], ())
            if len(join) > 1
        }
        new_lines, retag, added_gene = _rewrite_gene(gene_id, gene_rows, cds_rows, joins)
        replace.update(retag)
        if added_gene:
            insert_before[min(r.index for r in grows)] = [added_gene]
            insert_after[min(r.index for r in grows)] = new_lines
        else:
            insert_after[gene_rows[-1].index] = new_lines

    out: list[str] = []
    for row in rows:
        out += insert_before.get(row.index, [])
        out.append(replace.get(row.index, row.line))
        out += insert_after.get(row.index, [])
    return out


def cds_joins_from_genbank(genbank_text: str) -> list[list[Block]]:
    """Blocks of every multi-interval ``CDS`` feature of a GenBank flatfile (1-based)."""
    from viralscan.scripts.ncbi_fetch import _parse_location, iter_features

    joins = []
    for key, location, _qualifiers in iter_features(genbank_text):
        if key == "CDS" and "join(" in location:
            blocks = [(s, e) for s, e, _strand in _parse_location(location)]
            if len(blocks) > 1:
                joins.append(blocks)
    return joins


def normalise_gtf_file(path: str | Path, cache_dir: str | Path | None = None) -> list[str]:
    """:func:`normalise_viral_gtf` of a GTF file, with joins from the cached flatfiles.

    The flatfile of each seqname is read from the ``ncbi_fetch`` cache when present; without
    it a collapsed join (HHV-6B) cannot be seen and keeps its span.
    """
    from viralscan.scripts.ncbi_fetch import NCBIFetchError, genbank_cache_path

    lines = Path(path).read_text().splitlines()
    joins: dict[str, list[list[Block]]] = {}
    for seq in sorted({line.split("\t", 1)[0] for line in lines if line and line[0] != "#"}):
        try:
            flatfile = genbank_cache_path(seq, cache_dir)
        except NCBIFetchError:
            continue  # not an accession, so there is no flatfile to consult
        if flatfile.is_file():
            joins[seq] = cds_joins_from_genbank(flatfile.read_text())
    return normalise_viral_gtf(lines, joins)

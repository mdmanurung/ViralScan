"""One self-contained igv.js page for the reads behind a viral call (PLAN ``EVID-04``, ``viralscan evidence --html``).

The page lets a reader look at the alignments that sit under the tables of ``viralscan evidence``. It embeds
igv.js 2.15.11 (``vendor/igv/igv.min.js``, MIT, licence in the page and in the wheel), the viral genome, a BAM
of reads, an optional BAM of splice reads and a BED12 gene track. Nothing is fetched at view time.

The population is the unit of ``docs/evidence_unit.md``: the member reads of counted virus-best molecules. A page
draws at most ``MAX_READS`` of them per viral contig, chosen with a fixed seed, and at most ``MAX_SPLICED``
splice reads per contig from ``splice/splice.bam``. The page states the rule, the population sizes and the genes
left out. **No number on the page is computed from the drawn sample**: counts, junctions and verdicts live in the
TSV tables, which use every molecule.

When the page would exceed ``MAX_HTML_BYTES`` the read sample is halved until it fits, and the page says so. It
is never cut without a statement.
"""

import base64
import csv
import hashlib
import html
import json
import random
import re
import tempfile
from collections import Counter
from collections.abc import Iterable, Sequence
from importlib import resources
from pathlib import Path
from typing import Any

from viralscan.evidence import _run
from viralscan.molecule_verdict import _read_number
from viralscan.run_safety import sha256_file
from viralscan.splice import (
    LINEAGE_FILE,
    MAX_CLIP_FRACTION,
    MAX_INTRON_LENGTH,
    MAX_MISMATCH_FRACTION,
    MIN_ANCHOR,
    MIN_INTRON_LENGTH,
    RAW_BAM,
    SPLICE_DIR,
    VERDICTS_FILE,
    SpliceError,
    _member_reads,
    _molecule_categories,
    _stream_bam,
    junctions_of,
)
from viralscan.utils import open_maybe_gzip

# ---------------------------------------------------------------------------
# Page rules. Fixed before any real data were looked at; change them only with a PLAN row.
# ---------------------------------------------------------------------------

IGV_VERSION = "2.15.11"
#: Reads drawn per viral contig. Keeps the page small and igv.js responsive.
MAX_READS = 4000
#: Splice reads drawn per viral contig (reads with a junction that passed the ``--splice`` rules).
MAX_SPLICED = 1500
#: Fixed seed of every random draw; the seed of one contig is ``f"{SEED}:{contig}"``.
SEED = 20261103
#: Page size cap in bytes. Above it the read sample is halved (down to ``MIN_READS``) and the page says so.
MAX_HTML_BYTES = 25_000_000
#: Smallest read sample per contig the size rule may reach. Below that the page is refused.
MIN_READS = 100
#: igv.js loads at most this many reads per window (``samplingDepth``); 10000 is its maximum.
IGV_SAMPLING_DEPTH = 10_000
#: Gene track: a feature longer than this (EBV latency transcripts span most of the genome) is left out and counted.
MAX_FEATURE_BP = 20_000
#: CDS of one gene further apart than this are separate copies (terminal repeats), one BED row each.
COPY_GAP = 10_000

HTML_FILE = "igv_browser.html"
TEMPLATE = "igv_browser.html"
#: Content-Security-Policy of the page: no network, only inline script/style and data/blob URIs.
CSP = (
    "default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; "
    "img-src data: blob:; font-src data:; connect-src data: blob:; worker-src blob:"
)
STRAND_NOTE = (
    "Reads are drawn as aligned to the viral genome. 10x 3' reads are sense to the gene they cover and "
    "10x 5' reads are antisense, so the read strand against the gene strand depends on the chemistry."
)
NO_METRIC_NOTE = (
    "This page is for looking. It holds a sample of the reads, so count nothing here: counts, junctions and "
    "verdicts are in the TSV tables of the evidence directory, which use every molecule."
)


class IgvError(ValueError):
    """The page cannot be built as specified (missing contig, page too big, unreadable input)."""


class IgvNotEvaluable(IgvError):
    """There is no molecule unit to draw from (competitor mode, no verdicts, no counted virus-best molecule)."""


# ---------------------------------------------------------------------------
# Vendored igv.js
# ---------------------------------------------------------------------------


def vendor_text(name: str) -> str:
    """Text of a file in ``viralscan/vendor/igv/``."""
    return resources.files("viralscan").joinpath(f"vendor/igv/{name}").read_text(encoding="utf-8")


def _template() -> str:
    return (
        resources.files("viralscan").joinpath(f"templates/{TEMPLATE}").read_text(encoding="utf-8")
    )


# ---------------------------------------------------------------------------
# Genome and genes
# ---------------------------------------------------------------------------


def plain_contig(name: str) -> str:
    """Contig name without the competitive BAM's ``VIRUS|`` prefix."""
    return name.removeprefix("VIRUS|")


def read_genome(viral_fasta: str | Path, contigs: Iterable[str]) -> dict[str, str]:
    """Sequences of the named contigs (first header token) of a FASTA; a missing one is an error."""
    wanted = set(contigs)
    genome: dict[str, str] = {}
    name: str | None = None
    with open_maybe_gzip(viral_fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in genome:
                    raise IgvError(f"Duplicate FASTA ID {name!r} in {viral_fasta}")
                if name in wanted:
                    genome[name] = ""
            elif name in genome:
                genome[name] += line.strip()
    absent = wanted - set(genome)
    if absent:
        raise IgvError(f"{sorted(absent)} are in the BAM but not in {viral_fasta}")
    return genome


def fasta_and_index(genome: dict[str, str], width: int = 60) -> tuple[bytes, bytes]:
    """FASTA bytes (fixed line width) and the matching ``.fai`` bytes, contigs in dict order."""
    fasta = bytearray()
    fai: list[str] = []
    for name, seq in genome.items():
        header = f">{name}\n".encode()
        fasta += header
        fai.append(f"{name}\t{len(seq)}\t{len(fasta)}\t{width}\t{width + 1}\n")
        for i in range(0, len(seq), width):
            fasta += seq[i : i + width].encode() + b"\n"
    return bytes(fasta), "".join(fai).encode()


def bed12(
    gtfs: Sequence[Path], contigs: Iterable[str]
) -> tuple[str, list[tuple[str, str, int]], int]:
    """``(BED12 text, [(contig, name, length)] left out as too long, rows read)`` for *contigs* from GTFs.

    One row per transcript (a gene without transcript rows: per gene), blocks from exon rows, else CDS rows,
    else the span. A CDS-only gene reaches out to the UTR ends of its gene row. A CDS-only gene annotated in
    both terminal repeats gives one row per copy. igv.js labels BED12 rows; it did not label the transcripts of
    a combined GTF.
    """
    wanted = set(contigs)
    by_contig: dict[str, list[tuple[list[str], dict[str, str]]]] = {c: [] for c in sorted(wanted)}
    n_rows = 0
    for path in gtfs:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                f = line.rstrip("\n").split("\t")
                if line.startswith("#") or len(f) != 9 or plain_contig(f[0]) not in wanted:
                    continue
                if f[2] not in ("gene", "transcript", "exon", "CDS"):
                    continue
                by_contig[plain_contig(f[0])].append(
                    (f, dict(re.findall(r'(\w+)\s+"([^"]*)"', f[8])))
                )
                n_rows += 1
    lines: list[str] = []
    left_out: list[tuple[str, str, int]] = []
    for contig, rows in by_contig.items():
        with_transcripts = {a.get("gene_id", "") for f, a in rows if f[2] == "transcript"}
        features: dict[str, dict[str, Any]] = {}
        for f, attr in rows:
            gene = attr.get("gene_id", "")
            if f[2] == "gene" and gene in with_transcripts:
                continue  # drawn per transcript instead
            key = (
                attr.get("transcript_id", gene)
                if gene in with_transcripts and f[2] != "gene"
                else gene
            )
            label = attr.get("gene_name") or attr.get("gene") or gene.split("_")[-1] or key
            feat = features.setdefault(
                key, {"name": label, "strand": f[6], "exon": [], "CDS": [], "span": []}
            )
            span = (int(f[3]) - 1, int(f[4]))
            feat["span" if f[2] in ("gene", "transcript") else f[2]].append(span)
        bed: list[tuple[int, int, str, str, int, int, list[tuple[int, int]]]] = []
        for feat in features.values():
            blocks = sorted(feat["exon"] or feat["CDS"] or feat["span"])
            if not blocks:
                continue
            copies = [[blocks[0]]]
            for block in blocks[1:]:
                if not feat["exon"] and block[0] - copies[-1][-1][1] > COPY_GAP:
                    copies.append([block])
                else:
                    copies[-1].append(block)
            for copy in copies:
                if not feat["exon"]:  # CDS-only blocks reach out to the gene row's UTR ends
                    for g0, g1 in feat["span"]:
                        if g0 <= copy[0][0] and copy[-1][1] <= g1 and g1 - g0 < COPY_GAP:
                            copy = [(g0, copy[0][1]), *copy[1:]]
                            copy[-1] = (copy[-1][0], g1)
                start, end = copy[0][0], max(e for _, e in copy)
                if end - start > MAX_FEATURE_BP:
                    left_out.append((contig, feat["name"], end - start))
                    continue
                thick = [c for c in feat["CDS"] if start <= c[0] < end] or copy
                bed.append(
                    (start, end, feat["name"], feat["strand"], min(t[0] for t in thick),
                     max(t[1] for t in thick), copy)
                )  # fmt: skip
        for s, e, name, strand, ts, te, blocks in sorted(bed):
            sizes = ",".join(str(b1 - b0) for b0, b1 in blocks)
            starts = ",".join(str(b0 - s) for b0, _ in blocks)
            lines.append(
                f"{contig}\t{s}\t{e}\t{name}\t0\t{strand}\t{ts}\t{te}\t0\t{len(blocks)}\t{sizes}\t{starts}\n"
            )
    return "".join(lines), left_out, n_rows


# ---------------------------------------------------------------------------
# Reads: seeded draw from the BAM stream
# ---------------------------------------------------------------------------


class Reservoir:
    """Uniform random sample of at most *cap* items from a stream (algorithm R), seeded per contig."""

    def __init__(self, cap: int, seed: str) -> None:
        self.cap = cap
        self.rng = random.Random(seed)
        self.items: list[list[str]] = []
        self.seen = 0

    def add(self, item: list[str]) -> None:
        self.seen += 1
        if len(self.items) < self.cap:
            self.items.append(item)
        else:
            j = self.rng.randrange(self.seen)
            if j < self.cap:
                self.items[j] = item

    def draw(self, cap: int, seed: str) -> list[list[str]]:
        """At most *cap* of the kept items, uniformly (a uniform sample of a uniform sample)."""
        if cap >= len(self.items):
            return list(self.items)
        return random.Random(seed).sample(self.items, cap)


def _clean(fields: list[str], contig: str) -> list[str]:
    """SAM fields renamed to the plain contig, mate fields blank, ``SA`` tags dropped (they name other contigs)."""
    out = list(fields)
    out[2] = contig
    out[6], out[7], out[8] = "*", "0", "0"
    return out[:11] + [t for t in out[11:] if not t.startswith("SA:Z:")]


def _sample_member_reads(
    bam: Path, reads: dict[int, tuple[str, str]], cap: int
) -> tuple[dict[str, Reservoir], Counter[str]]:
    """Reservoirs of viral primary records of member reads, per contig, and counts of what was not drawable."""
    pools: dict[str, Reservoir] = {}
    count: Counter[str] = Counter()
    for fields in _stream_bam(bam):
        number = _read_number(fields[0])
        if number not in reads:
            continue
        flag = int(fields[1])
        if flag & 0x104:  # secondary, or unmapped: the primary record speaks for the read
            continue
        count["member_reads_primary"] += 1
        if not fields[2].startswith("VIRUS|"):
            count["primary_on_host"] += 1
            continue
        contig = plain_contig(fields[2])
        pool = pools.setdefault(contig, Reservoir(cap, f"{SEED}:{contig}"))
        pool.add(_clean(fields, contig))
    return pools, count


def _sample_splice_reads(bam: Path, cap: int) -> dict[str, Reservoir]:
    """Reservoirs of splice-BAM records that carry a junction accepted by the ``--splice`` rules."""
    pools: dict[str, Reservoir] = {}
    for fields in _stream_bam(bam):
        flag = int(fields[1])
        nm = next((int(t[5:]) for t in fields[11:] if t.startswith("NM:i:")), None)
        accepted, _ = junctions_of(flag, int(fields[3]), fields[5], nm, len(fields[9]))
        if not accepted:
            continue
        contig = fields[2]
        pool = pools.setdefault(contig, Reservoir(cap, f"{SEED}:splice:{contig}"))
        pool.add(_clean(fields, contig))
    return pools


def _bam_bytes(genome: dict[str, str], records: list[list[str]]) -> tuple[bytes, bytes]:
    """Sorted BAM and BAI bytes for SAM-field records on *genome*'s contigs (samtools, list-form)."""
    with tempfile.TemporaryDirectory() as tmp:
        sam, bam = Path(tmp) / "reads.sam", Path(tmp) / "reads.bam"
        header = ["@HD\tVN:1.6\tSO:coordinate"] + [
            f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s in genome.items()
        ]
        rank = {name: i for i, name in enumerate(genome)}
        body = ["\t".join(r) for r in sorted(records, key=lambda r: (rank[r[2]], int(r[3]), r[0]))]
        sam.write_text("\n".join(header + body) + "\n", encoding="utf-8")
        _run(["samtools", "view", "-b", "--no-PG", "-o", str(bam), str(sam)])
        _run(["samtools", "index", str(bam)])
        return bam.read_bytes(), Path(f"{bam}.bai").read_bytes()


def _uri(data: bytes) -> str:
    return "data:application/octet-stream;base64," + base64.b64encode(data).decode("ascii")


def _json_in_script(obj: object) -> str:
    """JSON safe inside a ``<script>`` element."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).replace("</", "<\\/")


# ---------------------------------------------------------------------------
# Public entry
# ---------------------------------------------------------------------------


def _page(
    template: str,
    *,
    title: str,
    contigs: list[dict[str, Any]],
    data: dict[str, Any],
    receipt: dict[str, Any],
    gene_note: str,
    size_note: str,
) -> str:
    rows = "".join(
        f"<tr><td>{html.escape(c['name'])}</td><td>{c['length']}</td><td>{c['reads_available']}</td>"
        f"<td>{c['reads_drawn']}</td><td>{c['splice_available']}</td><td>{c['splice_drawn']}</td></tr>"
        for c in contigs
    )
    subs = {
        "__TITLE__": html.escape(title),
        "__CSP__": html.escape(CSP, quote=True),
        "__IGV_VERSION__": IGV_VERSION,
        "__CONTIG_ROWS__": rows,
        "__GENE_NOTE__": html.escape(gene_note),
        "__SIZE_NOTE__": html.escape(size_note),
        "__STRAND_NOTE__": html.escape(STRAND_NOTE),
        "__NO_METRIC_NOTE__": html.escape(NO_METRIC_NOTE),
        "__MAX_READS__": str(receipt["read_cap_used"]),
        "__MAX_SPLICED__": str(receipt["splice_cap_used"]),
        "__SEED__": str(SEED),
        "__POPULATION__": html.escape(receipt["population_text"]),
        "__SPLICE_RULE__": html.escape(receipt["splice_rule_text"]),
        "__LICENSE__": html.escape(vendor_text("LICENSE")),
        "__DATA_JSON__": _json_in_script(data),
        "__RECEIPT_JSON__": _json_in_script(receipt),
    }
    # the igv.js source goes in last: it is one long line and must not be scanned for placeholders
    page = re.sub("|".join(map(re.escape, subs)), lambda m: subs[m.group(0)], template)
    return page.replace("__IGV_JS__", vendor_text("igv.min.js"))


def _no_splice_text(ev: Path) -> str:
    """Why the page has no splice track: --splice never ran, or it ran and no junction passed the rules."""
    summary = ev / SPLICE_DIR / "splice_summary.tsv"
    if summary.is_file():
        with summary.open(newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"), {})
        if row.get("status") == "assessed":
            return "No splice track: --splice ran: 0 spliced reads passed the junction rules."
        if row.get("status"):
            return f"No splice track: --splice ran: {row['status']} ({row.get('reason', '')})."
    return "No splice track: this evidence directory has no splice/splice.bam (run viralscan evidence --splice)."


def write_igv_page(
    evidence_dir: str | Path,
    *,
    viral_fasta: str | Path | None = None,
    annotation_gtfs: Sequence[Path] = (),
    out_dir: str | Path | None = None,
    max_reads: int = MAX_READS,
    max_html_bytes: int = MAX_HTML_BYTES,
) -> Path:
    """Write ``igv_browser.html`` for one evidence directory and return its path.

    Raises ``IgvNotEvaluable`` when the directory has no molecule unit to draw from and ``IgvError`` for
    inconsistent inputs or a page that cannot fit *max_html_bytes* at ``MIN_READS`` reads per contig.
    The page has no timestamp: the same inputs give the same bytes.
    """
    ev = Path(evidence_dir)
    verdicts, lineage, bam = (ev / n for n in (VERDICTS_FILE, LINEAGE_FILE, RAW_BAM))
    if not verdicts.is_file():
        raise IgvNotEvaluable(
            f"{ev} has no {VERDICTS_FILE} (competitor mode or a run from before molecule verdicts)"
        )
    try:
        counts, molecules = _molecule_categories(verdicts)
    except SpliceError as exc:
        raise IgvError(str(exc)) from exc
    if not molecules:
        raise IgvNotEvaluable(f"{ev} has no counted virus-best molecule to draw")
    reads = _member_reads(lineage, molecules)
    pools, tally = _sample_member_reads(bam, reads, max_reads)
    if not pools:
        raise IgvNotEvaluable(f"No member read of {ev} has a primary alignment on a viral contig")

    splice_bam = ev / SPLICE_DIR / "splice.bam"
    splice_pools = (
        _sample_splice_reads(splice_bam, min(MAX_SPLICED, max_reads))
        if splice_bam.is_file()
        else {}
    )

    fasta = viral_fasta
    if fasta is None and (ev / SPLICE_DIR / "target_genome.fa").is_file():
        fasta = ev / SPLICE_DIR / "target_genome.fa"
    if fasta is None:
        raise IgvError(
            "The page needs the viral genome: pass --viral-fasta (or run with --splice first)"
        )
    order = sorted(
        {*pools, *splice_pools}, key=lambda c: (-(pools[c].seen if c in pools else 0), c)
    )
    genome = read_genome(fasta, order)
    genome = {c: genome[c] for c in order}
    fasta_bytes, fai_bytes = fasta_and_index(genome)
    bed_text, left_out, n_gtf_rows = bed12(annotation_gtfs, genome)

    manifest_path = ev / "evidence_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    label = manifest.get("target", {}).get("label", ev.name)
    n_member = len(reads)
    n_drawable = sum(p.seen for p in pools.values())
    population_text = (
        f"Member reads of the {counts['counted_virus_best']} counted virus-best molecules "
        f"({n_member} reads in {LINEAGE_FILE}); {n_drawable} have a primary alignment on a viral contig and can "
        f"be drawn, {tally['primary_on_host']} have their primary alignment on the host and are never drawn."
    )
    splice_rule_text = (
        "Splice track: reads of splice/splice.bam (minimap2 splice:sr on the viral genome alone) with a junction "
        f"that passes the --splice rules (anchor >= {MIN_ANCHOR} bp each side, intron {MIN_INTRON_LENGTH} bp to "
        f"{MAX_INTRON_LENGTH} bp, NM per aligned base <= {MAX_MISMATCH_FRACTION}, soft clip <= {MAX_CLIP_FRACTION} "
        f"of the read). At most {MAX_SPLICED} splice reads per contig, drawn the same way."
        if splice_pools
        else _no_splice_text(ev)
    )
    gene_note = (
        f"Gene track from {len(annotation_gtfs)} GTF file(s); {len(left_out)} feature(s) longer than "
        f"{MAX_FEATURE_BP} bp left out"
        + (
            ": " + ", ".join(f"{n} ({c}, {length} bp)" for c, n, length in sorted(left_out)[:20])
            if left_out
            else ""
        )
        + ("." if len(left_out) <= 20 else f" and {len(left_out) - 20} more.")
        if annotation_gtfs and bed_text
        else "No gene track: no annotation GTF with rows on these contigs was given."
    )

    inputs = {
        n: sha256_file(p)
        for n, p in ((VERDICTS_FILE, verdicts), (LINEAGE_FILE, lineage), (RAW_BAM, bam))
    }
    inputs["viral_fasta"] = sha256_file(fasta)
    if splice_pools:
        inputs[f"{SPLICE_DIR}/splice.bam"] = sha256_file(splice_bam)
    if manifest_path.is_file():
        inputs["evidence_manifest.json"] = sha256_file(manifest_path)
    inputs.update({f"gtf:{p}": sha256_file(p) for p in annotation_gtfs})

    template = _template()
    cap = max_reads
    while True:
        splice_cap = min(MAX_SPLICED, cap)
        drawn = {c: pools[c].draw(cap, f"{SEED}:{c}:{cap}") if c in pools else [] for c in genome}
        drawn_splice = {
            c: splice_pools[c].draw(splice_cap, f"{SEED}:splice:{c}:{splice_cap}")
            if c in splice_pools
            else []
            for c in genome
        }
        contigs = [
            {
                "name": c,
                "length": len(genome[c]),
                "reads_available": pools[c].seen if c in pools else 0,
                "reads_drawn": len(drawn[c]),
                "splice_available": splice_pools[c].seen if c in splice_pools else 0,
                "splice_drawn": len(drawn_splice[c]),
            }
            for c in genome
        ]
        reads_bam, reads_bai = _bam_bytes(genome, [r for c in genome for r in drawn[c]])
        data: dict[str, Any] = {
            "fasta": _uri(fasta_bytes),
            "fai": _uri(fai_bytes),
            "bam": _uri(reads_bam),
            "bai": _uri(reads_bai),
            "contigs": [c["name"] for c in contigs],
            "igv_sampling_depth": IGV_SAMPLING_DEPTH,
            "visibility_window": max(len(g) for g in genome.values()),
        }
        if splice_pools:
            splice_bam_bytes, splice_bai = _bam_bytes(
                genome, [r for c in genome for r in drawn_splice[c]]
            )
            data.update(splice_bam=_uri(splice_bam_bytes), splice_bai=_uri(splice_bai))
        if bed_text:
            data["bed"] = _uri(bed_text.encode())
        size_note = (
            f"The read sample was reduced from {max_reads} to {cap} reads per contig (and splice reads to "
            f"{splice_cap}) to keep the page under {max_html_bytes} bytes: a second uniform draw from the first, "
            f'seeded with the text "<seed>:<contig>:{cap}".'
            if cap < max_reads
            else f"No size reduction: the page fits {max_html_bytes} bytes at {cap} reads per contig."
        )
        receipt: dict[str, Any] = {
            "schema_version": "1",
            "check": "evid_04_html",
            "igv_js": {
                "version": IGV_VERSION,
                "sha256": hashlib.sha256(vendor_text("igv.min.js").encode()).hexdigest(),
            },
            "constants": {
                "MAX_READS": MAX_READS,
                "MAX_SPLICED": MAX_SPLICED,
                "SEED": SEED,
                "MAX_HTML_BYTES": max_html_bytes,
                "MIN_READS": MIN_READS,
                "MAX_FEATURE_BP": MAX_FEATURE_BP,
                "COPY_GAP": COPY_GAP,
            },
            "read_cap_used": cap,
            "splice_cap_used": splice_cap,
            "reduced_for_size": cap < max_reads,
            "category_counts": dict(counts),
            "member_reads": n_member,
            "member_reads_primary_on_host": tally["primary_on_host"],
            "contigs": contigs,
            "genes_left_out": [
                {"contig": c, "name": n, "bp": length} for c, n, length in sorted(left_out)
            ],
            "gtf_rows_read": n_gtf_rows,
            "inputs": inputs,
            "population_text": population_text,
            "splice_rule_text": splice_rule_text,
        }
        page = _page(
            template,
            title=f"viralscan evidence: {label}",
            contigs=contigs,
            data=data,
            receipt=receipt,
            gene_note=gene_note,
            size_note=size_note,
        )
        if len(page.encode("utf-8")) <= max_html_bytes:
            break
        if cap <= MIN_READS:
            raise IgvError(
                f"The page is {len(page.encode('utf-8'))} bytes at {cap} reads per contig, above the "
                f"{max_html_bytes} byte cap: the genome and igv.js alone are too large."
            )
        cap = max(MIN_READS, cap // 2)
    out = Path(out_dir) if out_dir is not None else ev
    out.mkdir(parents=True, exist_ok=True)
    path = out / HTML_FILE
    path.write_text(page, encoding="utf-8")
    return path

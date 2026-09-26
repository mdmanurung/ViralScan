"""One-time generator for ``src/viralscan/data/anellovirus_genes.tsv``.

The bug this fixes
------------------
The 2,042-accession Anelloviridae panel used to reach every built reference as
one placeholder gene per genome — ``{accession}_gene1`` spanning the whole
genome, emitted by
:func:`viralscan.scripts.build_reference._genome_as_transcript_gtf`.  The
GenBank path was never taken for the panel even though NCBI annotates real CDS
features: ``scripts/build_bundled_panel_ref.py`` and
``build_anellovirus_reference`` both discarded the GenBank-derived GTF that
``ncbi_fetch._fetch_one`` had already written and regenerated a placeholder
from the FASTA.

The consequence is not a naming inconvenience.  A whole-genome transcript
shares sequence with every other genome in the panel, so reads cross-map in
proportion to conservation and the panel's entire viral signal is absorbed by
whichever genome is most conserved.  Measured on a COVID scRNA-seq run:
1,169,272 anellovirus UMI in total, of which 1,167,103 (99.8 %) landed in
``MW455439.1_gene1`` and were reported as "Alphatorquevirus" — against 1,241
UMI in Betatorquevirus.  That ratio is a conservation rank, not a measurement.

What this generator does
------------------------
For each accession in the packaged ``anellovirus_accessions.tsv`` it reads the
cached GenBank flatfile (``~/.cache/viralscan/ncbi/<acc>/<acc>.gb``), falling
back to a single ``efetch`` when the flatfile is not cached, and writes one row
per real CDS feature.

It is cache-first and resumable
-------------------------------
2,042 GenBank retrievals at NCBI's 3 requests/second unauthenticated ceiling is
~12 minutes of wall clock, and a re-run must be free.  Two properties give that:

* the raw flatfile is cached by :func:`viralscan.scripts.ncbi_fetch.fetch_genbank`
  with a ``.sha256`` sidecar, so a second run reads from disk;
* nothing is written until every accession has been attempted, so an
  interrupted run loses no work beyond the flatfiles it had already cached.

``--limit`` / ``--accessions`` restrict the run to a subset so the generator can
be exercised without touching the network at all.

Measured coverage (full panel, not a sample)
-------------------------------------------
1,995 of 2,042 accessions (97.7 %) carry at least one CDS feature, totalling
2,515 genes.  The other 47 have a bare ``source`` feature and nothing else, so a
placeholder is the only honest annotation for them.  Per-accession CDS count:
1 CDS in 1,740 accessions, 2 in 71, 3 in 115, 4 in 58, 5 in 10, 6 in 1.

Anelloviridae are circular ssDNA, but NCBI annotates in a *linear*
representation and only 488 of the 2,042 records even declared ``circular`` on
the LOCUS line.  Exactly one CDS in the whole panel is origin-spanning —
``KU243129.1`` ``join(2677..2824,1..80)`` — and no interval anywhere falls
outside ``[1, length]``, so wrap-around is handled and flagged rather than
assumed away; see ``ncbi_fetch._origin_spans``.

Genogroup
---------
Genogroup is **not** derivable, and is not invented.  Measured over all 2,042
panel accessions, no record carries a ``/genogroup`` qualifier and exactly two
carry ``/genotype``: ``NC_014081.1`` (``"6"``, ``/organism`` ``Torque teno virus
3``) and ``NC_014094.1`` (``"28"``, organism ``Torque teno virus 6``).  Those
two contradict their own organism number, which is the direct evidence that the
species and genogroup namespaces are independent — inferring one from the other
would be wrong on the only two records where the question can be checked.

``source_genotype`` therefore carries NCBI's ``/genotype`` verbatim for those two
accessions and is empty for the other 2,040, rather than being back-filled from
free text.  Retained alongside it: ``source/isolate`` (1,946 records, laboratory
sample codes such as ``MDJHem2`` or ``SAfiA-468-6``) and ``source/strain`` (31),
so a downstream classifier can be fitted without re-fetching.

Usage
-----
    python extras/build_anellovirus_genes.py --email you@example.org
    python extras/build_anellovirus_genes.py --limit 50 --out /tmp/subset.tsv
    python extras/build_anellovirus_genes.py --email you@example.org \\
        --api-key $NCBI_API_KEY --min-interval 0.12

Exits non-zero when the annotated fraction of the panel falls below
``--min-coverage`` (default 0.95), so a silent NCBI regression cannot ship a
catalogue that is mostly placeholders.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from viralscan.anellovirus import load_accession_table  # noqa: E402
from viralscan.scripts.ncbi_fetch import (  # noqa: E402
    DEFAULT_CACHE_DIR,
    NCBIFetchError,
    catalogue_rows,
    fetch_genbank,
)

DATA_DIR = REPO_ROOT / "src" / "viralscan" / "data"
DEFAULT_OUT = DATA_DIR / "anellovirus_genes.tsv"
ACCESSIONS_TSV = DATA_DIR / "anellovirus_accessions.tsv"

TSV_COLUMNS = (
    "accession",
    "gene_id",
    "transcript_id",
    "gene",
    "gene_symbol",
    "product",
    "n_exons",
    "exons",
    "strand",
    "genome_length",
    "topology",
    "origin_spanning",
    "source_genotype",
    "source_isolate",
    "source_strain",
)

#: NCBI's unauthenticated ceiling is 3 requests/second. The default leaves a
#: margin under it and only serialises this generator, not other callers.
DEFAULT_MIN_INTERVAL = 0.34


def select_accessions(
    accessions: list[str] | None, limit: int | None, per_genus: int | None
) -> list[dict[str, str]]:
    """Return panel rows to process, in accession order.

    ``per_genus`` takes an even slice of every genus, which is what a coverage
    audit wants: the panel is 75 % Betatorquevirus, so a plain prefix would
    measure one genus and call it the panel.
    """
    rows = load_accession_table()
    if accessions:
        wanted = {a.strip() for a in accessions if a.strip()}
        rows = [row for row in rows if row["accession"].strip() in wanted]
    if per_genus:
        by_genus: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_genus[row["genus"].strip()].append(row)
        rows = [row for genus in sorted(by_genus) for row in by_genus[genus][:per_genus]]
    if limit:
        rows = rows[:limit]
    return rows


def collect(
    rows: list[dict[str, str]],
    cache_dir: Path,
    email: str | None,
    api_key: str | None,
    min_interval: float,
) -> tuple[list[dict[str, object]], list[str], list[dict[str, int]]]:
    """Return ``(gene_rows, failed_accessions, per_genus_counts)``.

    Failures are collected rather than raised so one dead accession does not
    abort a 2,042-accession run; the caller decides whether the loss is
    acceptable.
    """
    gene_rows: list[dict[str, object]] = []
    failures: list[str] = []
    per_genus: dict[str, int] = defaultdict(int)
    last_call = 0.0

    for index, row in enumerate(rows, 1):
        accession = row["accession"].strip()
        genus = row["genus"].strip()
        if index % 250 == 0:
            print(f"  progress: {index} / {len(rows)}", flush=True)
        wait = min_interval - (time.monotonic() - last_call)
        if wait > 0:
            time.sleep(wait)
        last_call = time.monotonic()
        try:
            _path, genbank_text = fetch_genbank(
                accession, email=email, api_key=api_key, cache_dir=cache_dir
            )
        except NCBIFetchError as exc:
            failures.append(f"{accession}: {exc}")
            continue
        catalogue = catalogue_rows(accession, genbank_text)
        if catalogue:
            per_genus[genus] += 1
            gene_rows.extend(catalogue)

    return gene_rows, failures, dict(per_genus)


def audit(gene_rows: list[dict[str, object]]) -> list[str]:
    """Return coordinate / uniqueness violations found in *gene_rows*.

    Every exon must lie inside ``[1, genome_length]``, the strand must be ``+``
    or ``-``, ``n_exons`` must agree with the exon list, and ``gene_id`` must be
    unique across the whole panel — the last is the property the old
    whole-genome placeholder violated by construction, since 150 genomes share
    the bare symbol ``ORF1``.
    """
    problems: list[str] = []
    seen: dict[str, str] = {}
    for row in gene_rows:
        gene_id = str(row["gene_id"])
        if gene_id in seen:
            problems.append(f"duplicate gene_id {gene_id} ({seen[gene_id]} and {row['accession']})")
        seen[gene_id] = str(row["accession"])
        length = int(row["genome_length"])
        if str(row["strand"]) not in ("+", "-"):
            problems.append(f"{gene_id}: strand {row['strand']!r}")
        blocks = [
            tuple(int(v) for v in piece.split(":"))  # type: ignore[misc]
            for piece in str(row["exons"]).split(",")
            if piece
        ]
        if len(blocks) != int(row["n_exons"]):
            problems.append(f"{gene_id}: n_exons {row['n_exons']} but {len(blocks)} exon block(s)")
        for start, end in blocks:
            if length and (start < 1 or end > length):
                problems.append(f"{gene_id}: exon {start}..{end} outside 1..{length}")
            if end < start:
                problems.append(f"{gene_id}: exon {start}..{end} is inverted")
    return problems


def report(
    gene_rows: list[dict[str, object]],
    panel_rows: list[dict[str, str]],
    per_genus: dict[str, int],
    failures: list[str],
) -> float:
    """Print a per-genus coverage table; return the annotated fraction."""
    total = len(panel_rows)
    annotated = sum(per_genus.values())
    by_genus_total: dict[str, int] = defaultdict(int)
    for row in panel_rows:
        by_genus_total[row["genus"].strip()] += 1

    print()
    print(f"{'genus':<22}{'panel':>7}{'annotated':>11}{'genes':>8}{'fraction':>10}")
    genus_of = {row["accession"].strip(): row["genus"].strip() for row in panel_rows}
    genes_by_genus: dict[str, int] = defaultdict(int)
    for row in gene_rows:
        genes_by_genus[genus_of.get(str(row["accession"]), "Anelloviridae")] += 1
    for genus in sorted(by_genus_total, key=lambda g: -by_genus_total[g]):
        count = per_genus.get(genus, 0)
        fraction = count / by_genus_total[genus] if by_genus_total[genus] else 0.0
        print(
            f"{genus:<22}{by_genus_total[genus]:>7}{count:>11}"
            f"{genes_by_genus.get(genus, 0):>8}{fraction:>9.1%}"
        )
    print(
        f"{'TOTAL':<22}{total:>7}{annotated:>11}{len(gene_rows):>8}"
        f"{annotated / total if total else 0.0:>9.1%}"
    )

    spliced = sum(1 for row in gene_rows if int(row["n_exons"]) > 1)
    wrapping = sum(1 for row in gene_rows if str(row["origin_spanning"]) == "true")
    circular = sum(1 for row in gene_rows if str(row["topology"]) == "circular")
    print()
    print(f"  spliced genes (>1 exon) : {spliced} / {len(gene_rows)}")
    print(f"  origin-spanning genes  : {wrapping} / {len(gene_rows)}")
    print(f"  records declaring circular on LOCUS: {circular} / {len(gene_rows)}")
    if failures:
        print(f"  fetch failures         : {len(failures)} / {total}")
        for failure in failures[:10]:
            print(f"    {failure}")
        if len(failures) > 10:
            print(f"    … and {len(failures) - 10} more")
    return annotated / total if total else 0.0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument(
        "--accessions",
        nargs="+",
        help="explicit accession subset (default: the whole packaged panel)",
    )
    parser.add_argument("--limit", type=int, help="process at most N panel rows")
    parser.add_argument("--per-genus", type=int, help="take at most N rows from each genus")
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.95,
        help="exit non-zero if the annotated fraction falls below this (default: 0.95)",
    )
    parser.add_argument(
        "--min-interval",
        type=float,
        default=DEFAULT_MIN_INTERVAL,
        help=f"seconds between efetch calls (default: {DEFAULT_MIN_INTERVAL})",
    )
    parser.add_argument(
        "--no-verify",
        action="store_true",
        help="skip the coordinate/uniqueness audit",
    )
    args = parser.parse_args(argv)

    if not args.email:
        parser.error("NCBI requires an email: pass --email or set NCBI_EMAIL")

    panel_rows = select_accessions(args.accessions, args.limit, args.per_genus)
    if not panel_rows:
        parser.error("no panel rows selected")
    print(f"Processing {len(panel_rows)} accessions (cache: {args.cache_dir}) → {args.out}")

    gene_rows, failures, per_genus = collect(
        panel_rows, args.cache_dir, args.email, args.api_key, args.min_interval
    )
    coverage = report(gene_rows, panel_rows, per_genus, failures)

    if not args.no_verify:
        problems = audit(gene_rows)
        if problems:
            for problem in problems[:20]:
                print(f"AUDIT: {problem}", file=sys.stderr)
            print(f"{len(problems)} audit violation(s); refusing to write.", file=sys.stderr)
            return 1
        print(f"  audit                 : OK ({len(gene_rows)} genes, no violations)")

    if coverage < args.min_coverage:
        print(
            f"ERROR: annotated fraction {coverage:.1%} is below "
            f"--min-coverage {args.min_coverage:.1%}; refusing to write.",
            file=sys.stderr,
        )
        return 1

    gene_rows.sort(
        key=lambda row: (str(row["accession"]), int(row["n_exons"]), str(row["gene_id"]))
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in gene_rows:
            writer.writerow({column: row.get(column, "") for column in TSV_COLUMNS})

    print(f"\nWrote {len(gene_rows)} genes for {len(per_genus)} accessions → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

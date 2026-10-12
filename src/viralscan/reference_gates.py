"""Build-time gates on a viral reference index (PLAN CAT-40, CAT-19, CAT-14, CAT-15, CAT-16).

Small, independent checks. None of them builds an index; each reads a finished artefact (a
kallisto index or its recorded ``kallisto inspect`` text, a host-homology table, a ``t2g``, a
GTF, a ``viral_summary.tsv``) and returns a list of problems. An empty list is a pass. No
threshold lives here: a gate fails on a *regression against a recorded baseline*, on an
*explicit budget the caller passes*, or on a plain structural fact (a marker that does not
resolve, an accession with no measurement).

Command line (``python -m viralscan.reference_gates <gate> --help``)::

    inspect          kallisto inspect -> index_inspect.json; fail on regression (CAT-40, CAT-19)
    host-crosstalk   per-accession host-homologous fraction table (CAT-14)
    host-only        a host-only library must report no viral call (CAT-14)
    programmes       gene-programme markers must resolve in the index t2g/GTF (CAT-15)
    scale            gene-count growth and downstream memory estimate (CAT-16)
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import subprocess
import sys
import time
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Optional

from viralscan.validation import kb_tools, tool_path

log = logging.getLogger("viralscan")

# --- CAT-40 / CAT-19: kallisto inspect --------------------------------------------------

#: ``kallisto inspect`` lines worth recording. The first two are the gate; the rest is context.
_INSPECT_FIELDS = {
    "max_ec_size": r"max EC size\s*=\s*([\d,]+)",
    "discarded_ec_count": r"number of ECs discarded\s*=\s*([\d,]+)",
    "n_unitigs": r"number of unitigs\s*=\s*([\d,]+)",
    "n_targets": r"number of targets:\s*([\d,]+)",
    "n_kmers": r"number of k-mers:\s*([\d,]+)",
    "n_dlist_kmers": r"number of D-list k-mers:\s*([\d,]+)",
    "k": r"k-mer length:\s*([\d,]+)",
    "index_version": r"Index version number\s*=\s*([\d,]+)",
}
INSPECT_GATED = ("max_ec_size", "discarded_ec_count")
INSPECT_RECORD = "index_inspect.json"


def parse_inspect(text: str) -> dict[str, int]:
    """Integers out of ``kallisto inspect`` text (stdout and stderr, in either order).

    Raises ``ValueError`` when a gated metric is absent: a kallisto that does not print
    ``max EC size`` or ``number of ECs discarded`` cannot be gated, and a silent zero would
    read as a pass.
    """
    found = {
        key: int(match.group(1).replace(",", ""))
        for key, pattern in _INSPECT_FIELDS.items()
        if (match := re.search(pattern, text))
    }
    missing = [key for key in INSPECT_GATED if key not in found]
    if missing:
        raise ValueError(f"kallisto inspect output has no {', '.join(missing)}")
    return found


def run_inspect(index: Path) -> dict[str, Any]:
    """``kallisto inspect`` on *index*, with the binary that built it (``tool_path``, REL-16).

    A conda ``kallisto`` of the same version aborts on an index ``kb ref`` built
    (``bad_array_new_length``, measured 2026-10-12), so this follows ``kb``'s own binary.
    """
    kallisto = tool_path("kallisto")
    if kallisto is None:
        raise RuntimeError("kallisto is not on PATH; use --from-text with a recorded inspect")
    done = subprocess.run(  # noqa: S603
        [kallisto, "inspect", str(index)], capture_output=True, text=True, check=False
    )
    if done.returncode != 0:
        raise RuntimeError(f"kallisto inspect failed ({done.returncode}): {done.stderr[-300:]}")
    record: dict[str, Any] = dict(parse_inspect(done.stdout + "\n" + done.stderr))
    record.update(
        kallisto=kallisto,
        kallisto_version=kb_tools().get("kallisto", (None, None))[0],
        index=str(index),
        index_bytes=index.stat().st_size,
    )
    return record


def check_ec_regression(current: dict[str, Any], baseline: dict[str, Any]) -> list[str]:
    """Problems when either gated metric is larger than in *baseline* (CAT-40).

    Larger is worse for both: more genomes sharing one k-mer set, more ECs dropped by
    kallisto's pseudo-inverse limit. Equal passes; a legitimate panel growth that raises
    them is accepted by recording a new baseline, which is a decision, not a default.
    """
    problems = []
    for key in INSPECT_GATED:
        now, before = current.get(key), baseline.get(key)
        if type(now) is not int or type(before) is not int:
            problems.append(f"{key}: current={now!r} baseline={before!r} are not both integers")
        elif now > before:
            problems.append(f"{key} regressed: {before} -> {now}")
    return problems


def record_index_inspect(index: Path, out: Optional[Path] = None) -> Optional[dict[str, Any]]:
    """Write ``index_inspect.json`` beside *index* (every build records it, CAT-40).

    Non-fatal for the caller: an inspect failure after hours of ``kb ref`` must not discard
    the index, so it returns None and logs; the missing file is what the gate then reports.
    """
    try:
        record = run_inspect(index)
    except (RuntimeError, ValueError, OSError) as exc:
        log.error("kallisto inspect not recorded for %s: %s", index, exc)
        return None
    (out or index.with_name(INSPECT_RECORD)).write_text(json.dumps(record, indent=2) + "\n")
    return record


# --- CAT-14: host cross-talk ------------------------------------------------------------


def _union_bases(intervals: Iterable[tuple[int, int]]) -> int:
    covered, end = 0, -1
    for start, stop in sorted(intervals):
        start = max(start, end)
        if stop > start:
            covered += stop - start
            end = stop
    return covered


def fasta_lengths(path: Path) -> dict[str, int]:
    """``{first header token: sequence length}`` for a FASTA."""
    lengths: dict[str, int] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        elif name is not None:
            lengths[name] += len(line.strip())
    return lengths


HOST_CROSSTALK_COLUMNS = (
    "accession_version",
    "status",
    "genome_length",
    "host_homologous_bases",
    "host_homologous_fraction",
    "n_loci",
    "max_identity",
    "max_alignment_coverage",
)


def host_crosstalk_report(
    table: Path, expected: Optional[dict[str, int]] = None
) -> tuple[list[dict[str, Any]], list[str]]:
    """Per-accession host-homologous fraction from ``host_homology_annotations.tsv`` (CAT-14).

    The fraction is the union of the viral intervals of every recorded host locus over the
    genome length, so overlapping loci count once; ``max_alignment_coverage`` is the single
    best alignment's share, which is what the table's own column holds. The length comes
    from *expected* (the indexed FASTA); without it the fraction is blank, never 0.
    An accession with no row, or whose status is not ``measured``, is a problem: an unmeasured
    genome must not look host-free.
    """
    with table.open(newline="") as handle:
        rows = {r["accession_version"]: r for r in csv.DictReader(handle, delimiter="\t")}
    report: list[dict[str, Any]] = []
    problems: list[str] = []
    for acc in sorted(set(rows) | set(expected or {})):
        row = rows.get(acc)
        if row is None:
            problems.append(f"{acc}: no host-homology row")
            continue
        status = row.get("host_homology_status", "")
        if status != "measured":
            problems.append(f"{acc}: host homology is {status or 'not recorded'}, not measured")
        loci = json.loads(row.get("host_homology_loci") or "[]")
        bases = _union_bases((int(x["viral_start"]), int(x["viral_end"])) for x in loci)
        length = (expected or {}).get(acc)
        report.append(
            {
                "accession_version": acc,
                "status": status,
                "genome_length": length if length else "",
                "host_homologous_bases": bases if status == "measured" else "",
                "host_homologous_fraction": (
                    round(bases / length, 6) if length and status == "measured" else ""
                ),
                "n_loci": len(loci),
                "max_identity": row.get("host_homology_max_identity", ""),
                "max_alignment_coverage": row.get("host_homology_max_query_coverage", ""),
            }
        )
    return report, problems


def host_only_calls(viral_summary: Path) -> list[str]:
    """Viruses a host-only library reported as a call: ``detected`` or any called infected cell.

    A host-only negative must produce none (CAT-14). Reads ``results/viral_summary.tsv``.
    """
    with viral_summary.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if rows and "detected" not in rows[0]:
        raise ValueError(f"{viral_summary} has no 'detected' column; not a v3 viral_summary")
    calls = []
    for row in rows:
        called = row.get("infected_called") or "0"
        if row["detected"].strip().lower() in {"true", "1"} or float(called) > 0:
            calls.append(row.get("virus_name", "?"))
    return calls


# --- CAT-15: gene programmes survive the index ------------------------------------------


def t2g_gene_ids(t2g: Path) -> set[str]:
    """Gene IDs (column 2) of a kb ``t2g``."""
    with t2g.open() as handle:
        return {parts[1] for line in handle if len(parts := line.rstrip("\n").split("\t")) > 1}


def check_programme_markers(
    gene_ids: Iterable[str],
    viruses: Optional[Iterable[str]] = None,
    panel_form: str = "bundled",
    catalogue: Optional[list[dict[str, str]]] = None,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Every catalogue marker of each virus must be a gene of the index (CAT-15).

    Resolution is ``gene_programs.resolve_markers`` (exact string match), so this fails
    exactly when a run's programme calling would silently lose markers. *viruses* defaults to
    every virus in the catalogue.
    """
    from viralscan.gene_programs import load_catalogue, resolve_markers

    rows = catalogue if catalogue is not None else load_catalogue()
    ids = set(gene_ids)
    names = sorted(set(viruses) if viruses else {r["virus"] for r in rows})
    report: list[dict[str, Any]] = []
    problems: list[str] = []
    for virus in names:
        resolved, unresolved = resolve_markers(virus, ids, panel_form, rows)
        if not resolved and not unresolved:
            problems.append(f"{virus}: no catalogue rows")
        report.append(
            {
                "virus": virus,
                "resolved": len(resolved),
                "unresolved": [u.gene_id for u in unresolved],
            }
        )
        if unresolved:
            problems.append(
                f"{virus}: {len(unresolved)} of {len(resolved) + len(unresolved)} markers do not "
                f"resolve ({', '.join(u.gene_id for u in unresolved[:5])})"
            )
    return report, problems


# --- CAT-16: scale check ----------------------------------------------------------------

#: The row's own expectation ("gene count grows ~10-50x"); exceeding it is a warning only.
EXPECTED_GROWTH = (10, 50)
FLOAT64_BYTES = 8
_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')


def gtf_scale(gtf: Path) -> dict[str, int]:
    """Gene, transcript, contig and exon-base counts of a GTF."""
    genes: set[str] = set()
    transcripts: set[str] = set()
    contigs: set[str] = set()
    exon_bp = 0
    with gtf.open() as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(cols) < 9:
                continue
            attrs = dict(_ATTR.findall(cols[8]))
            contigs.add(cols[0])
            if "gene_id" in attrs:
                genes.add(attrs["gene_id"])
            if "transcript_id" in attrs:
                transcripts.add(attrs["transcript_id"])
            if cols[2] == "exon":
                exon_bp += int(cols[4]) - int(cols[3]) + 1
    return {
        "n_genes": len(genes),
        "n_transcripts": len(transcripts),
        "n_contigs": len(contigs),
        "exon_bp": exon_bp,
    }


def scale_report(
    gtf: Path,
    baseline_gtf: Optional[Path] = None,
    n_cells: Optional[int] = None,
    max_genes: Optional[int] = None,
    max_dense_gib: Optional[float] = None,
) -> tuple[dict[str, Any], list[str]]:
    """Gene-count growth and the dense-matrix memory a detection run would need (CAT-16).

    ``analysis.py`` is timed on this GTF (its own ``extract_gene_ids``). ``detection.py``
    densifies the cells x index-genes submatrix of one virus at a time
    (``matrix_for_genes(...).toarray()``), so the estimate is
    ``n_cells * max genes in any one virus * 8`` bytes **per dense copy**; how many copies are live at
    once is not assumed. ``accession_breadth`` over all index genes is ANDET-01 (done). Fails
    only on the explicit ``max_genes`` / ``max_dense_gib`` budgets; growth beyond the row's
    10-50x expectation is a warning.
    """
    from viralscan.scripts.analysis import extract_gene_ids
    from viralscan.virus_catalog import merged_name_map
    from viralscan.virus_grouping import group_genes_by_virus
    from viralscan.virus_identity import gtf_gene_ids

    gene_ids = sorted(gtf_gene_ids(gtf))
    report: dict[str, Any] = {"gtf": str(gtf), **gtf_scale(gtf)}
    started = time.perf_counter()
    with gtf.open() as handle:
        extract_gene_ids(handle)
    report["analysis_extract_seconds"] = round(time.perf_counter() - started, 3)
    problems: list[str] = []
    warnings: list[str] = []
    if baseline_gtf is not None:
        before = gtf_scale(baseline_gtf)
        report["baseline"] = {"gtf": str(baseline_gtf), **before}
        if before["n_genes"]:
            growth = report["n_genes"] / before["n_genes"]
            report["gene_growth_factor"] = round(growth, 2)
            if growth > EXPECTED_GROWTH[1]:
                warnings.append(f"gene count grew {growth:.1f}x, above the expected ~10-50x")
    if n_cells is not None:
        # detection.py densifies one virus's index genes at a time, not the whole panel.
        groups, _ = group_genes_by_virus(gene_ids, merged_name_map())
        biggest = max(groups, key=lambda v: len(groups[v]), default="")
        report["n_viruses"] = len(groups)
        report["max_genes_per_virus"] = len(groups.get(biggest, []))
        report["max_genes_virus"] = biggest
        gib = n_cells * report["max_genes_per_virus"] * FLOAT64_BYTES / 2**30
        report["dense_gib_per_copy"] = round(gib, 3)
        if max_dense_gib is not None and gib > max_dense_gib:
            problems.append(f"dense viral matrix {gib:.2f} GiB per copy exceeds {max_dense_gib}")
    if max_genes is not None and report["n_genes"] > max_genes:
        problems.append(f"{report['n_genes']} genes exceeds budget {max_genes}")
    report["warnings"] = warnings
    return report, problems


# --- build hooks ------------------------------------------------------------------------


def run_scale_gate(gtf: Path, baseline_gtf: Optional[Path], out: Path) -> str:
    """Write ``scale_check.json`` before the index build; return a one-line summary (CAT-16)."""
    report, _ = scale_report(gtf, baseline_gtf)
    out.write_text(json.dumps(report, indent=2) + "\n")
    growth = report.get("gene_growth_factor")
    return f"{report['n_genes']} genes" + (f", {growth}x the baseline" if growth else "")


def run_build_gates(
    out_dir: Path, index: Path, t2g: Path, ec_baseline: Optional[Path] = None
) -> list[str]:
    """Post-``kb ref`` gates of ``build_bundled_panel_ref.py``; returns the failures.

    Records ``index_inspect.json`` and ``programme_markers.json`` in *out_dir* either way.
    """
    failures: list[str] = []
    record = record_index_inspect(index, out_dir / INSPECT_RECORD)
    if record is None:
        failures.append("kallisto inspect could not be recorded (CAT-40)")
    elif ec_baseline is not None:
        failures += check_ec_regression(record, json.loads(ec_baseline.read_text()))
    report, problems = check_programme_markers(t2g_gene_ids(t2g))
    (out_dir / "programme_markers.json").write_text(json.dumps(report, indent=2) + "\n")
    return failures + problems


# --- command line -----------------------------------------------------------------------


def _finish(problems: list[str]) -> int:
    for problem in problems:
        print(f"FAIL {problem}", file=sys.stderr)
    return int(bool(problems))


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="gate", required=True)

    p = sub.add_parser("inspect", help="record kallisto inspect; fail on regression")
    p.add_argument("index", type=Path, nargs="?", help="kallisto index (omit with --from-text)")
    p.add_argument("--from-text", type=Path, help="recorded `kallisto inspect` text")
    p.add_argument(
        "--out", type=Path, help=f"record to write (default: <index dir>/{INSPECT_RECORD})"
    )
    p.add_argument("--baseline", type=Path, help="earlier record; any gated metric above it fails")
    p.add_argument("--max-ec-size", type=int)
    p.add_argument("--max-discarded-ecs", type=int)

    p = sub.add_parser("host-crosstalk", help="per-accession host-homologous fraction")
    p.add_argument("--table", type=Path, help="host_homology_annotations.tsv")
    p.add_argument("--viral-fasta", type=Path, help="the indexed viral.fa (lengths, coverage)")
    p.add_argument("--host-genome", type=Path, help="with --viral-fasta and no --table: measure")
    p.add_argument("--out", type=Path, required=True, help="report TSV to write")

    p = sub.add_parser("host-only", help="a host-only library must report no viral call")
    p.add_argument("viral_summary", type=Path)

    p = sub.add_parser("programmes", help="gene-programme markers must resolve")
    p.add_argument("--t2g", type=Path)
    p.add_argument("--gtf", type=Path)
    p.add_argument("--virus", action="append", help="repeatable; default every catalogue virus")
    p.add_argument("--panel-form", default="bundled", choices=["bundled", "merged", "starsolo"])

    p = sub.add_parser("scale", help="gene-count growth and memory estimate")
    p.add_argument("--gtf", type=Path, required=True)
    p.add_argument("--baseline-gtf", type=Path)
    p.add_argument("--n-cells", type=int)
    p.add_argument("--max-genes", type=int)
    p.add_argument("--max-dense-gib", type=float)
    p.add_argument("--out", type=Path)

    args = parser.parse_args(argv)

    if args.gate == "inspect":
        if bool(args.index) == bool(args.from_text):
            parser.error("give exactly one of INDEX or --from-text")
        record = (
            dict(parse_inspect(args.from_text.read_text()))
            if args.from_text
            else run_inspect(args.index)
        )
        out = args.out or (args.index.with_name(INSPECT_RECORD) if args.index else None)
        if out:
            out.write_text(json.dumps(record, indent=2) + "\n")
        print(json.dumps(record, indent=2))
        problems = []
        if args.baseline:
            problems += check_ec_regression(record, json.loads(args.baseline.read_text()))
        for key, limit in (
            ("max_ec_size", args.max_ec_size),
            ("discarded_ec_count", args.max_discarded_ecs),
        ):
            if limit is not None and record[key] > limit:
                problems.append(f"{key}={record[key]} exceeds budget {limit}")
        return _finish(problems)

    if args.gate == "host-crosstalk":
        table = args.table
        if table is None:
            if not (args.viral_fasta and args.host_genome):
                parser.error("give --table, or --viral-fasta and --host-genome to measure")
            from viralscan.scripts.build_reference import measure_host_homology

            table = args.out.with_name("host_homology_annotations.tsv")
            measure_host_homology(args.viral_fasta, args.host_genome, table)
        lengths = fasta_lengths(args.viral_fasta) if args.viral_fasta else None
        report, problems = host_crosstalk_report(table, lengths)
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, HOST_CROSSTALK_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(report)
        print(f"{len(report)} accessions written to {args.out}")
        return _finish(problems)

    if args.gate == "host-only":
        calls = host_only_calls(args.viral_summary)
        return _finish([f"host-only library reported a call for {v}" for v in calls])

    if args.gate == "programmes":
        if bool(args.t2g) == bool(args.gtf):
            parser.error("give exactly one of --t2g or --gtf")
        if args.t2g:
            ids = t2g_gene_ids(args.t2g)
        else:
            from viralscan.virus_identity import gtf_gene_ids

            ids = gtf_gene_ids(args.gtf)
        report, problems = check_programme_markers(ids, args.virus, args.panel_form)
        print(json.dumps(report, indent=2))
        return _finish(problems)

    report, problems = scale_report(
        args.gtf, args.baseline_gtf, args.n_cells, args.max_genes, args.max_dense_gib
    )
    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return _finish(problems)


if __name__ == "__main__":
    raise SystemExit(main())

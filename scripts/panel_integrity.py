#!/usr/bin/env python3
"""Offline catalogue/GTF integrity matrix and normalized EC-summary budget gate.

The EC input is a normalized JSON object with integer ``max_ec_size`` and
``discarded_ec_count`` fields. This helper does not run or parse kallisto inspect.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any

from viralscan.constants import SIBLING_VIRUS_PAIRS
from viralscan.gtf_normalise import normalise_gtf_file
from viralscan.virus_catalog import load_catalogue

PANELS = {"shipped", "max", "broad", "legacy"}
ATTR = re.compile(r'(\w+)\s+"([^\"]*)"\s*;')


def check_ec_budget(summary: dict[str, Any], max_ec_size: int, max_discarded_ecs: int) -> None:
    """Reject absent/invalid metrics or a summary exceeding explicit budgets."""
    if not isinstance(summary, dict):
        raise ValueError("EC summary must be a JSON object")
    for key, limit in (("max_ec_size", max_ec_size), ("discarded_ec_count", max_discarded_ecs)):
        value = summary.get(key)
        if type(value) is not int or value < 0 or type(limit) is not int or limit < 0:
            raise ValueError(f"{key} and its budget must be nonnegative integers")
        if value > limit:
            raise ValueError(f"{key}={value} exceeds budget {limit}")


def inspect_gtfs(
    paths: list[Path], normalise: bool = False
) -> tuple[dict[str, set[str]], set[str], list[str]]:
    """Validate rows, unique gene ownership and consistent transcript exon structure.

    ``normalise=True`` inspects what ``kb ref`` is given (the REF-13 normalised view, see
    :mod:`viralscan.gtf_normalise`) rather than the shipped files as they are.

    Repeated gene IDs within one locus (gene/transcript/exon rows) are expected;
    reuse in another file or contig is an error. Whole-genome IDs require the
    emitter's explicit biotype tag, never a length or name heuristic.
    """
    seq_genes: dict[str, set[str]] = {}
    owners: dict[str, tuple[Path, str, str]] = {}
    transcript_owners: dict[str, tuple[Path, str, str, str]] = {}
    whole_genome: set[str] = set()
    errors = []
    for path in paths:
        seen_rows: set[str] = set()
        needs_exon: set[str] = set()
        has_exon: set[str] = set()
        lines = normalise_gtf_file(path) if normalise else path.read_text().splitlines()
        for number, raw in enumerate(lines, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            location = f"{path.name}:{number}"
            cols = raw.split("\t")
            if len(cols) != 9:
                errors.append(f"{location}: expected nine GTF columns")
                continue
            try:
                start, end = int(cols[3]), int(cols[4])
                if (
                    start < 1
                    or end < start
                    or cols[6] not in {"+", "-", "."}
                    or cols[7] not in {".", "0", "1", "2"}
                ):
                    raise ValueError
            except ValueError:
                errors.append(f"{location}: invalid coordinates or strand")
                continue
            attrs = dict(ATTR.findall(cols[8]))
            gene, tx = attrs.get("gene_id", ""), attrs.get("transcript_id", "")
            if not gene or (cols[2] in {"exon", "CDS", "transcript"} and not tx):
                errors.append(f"{location}: missing gene_id/transcript_id")
                continue
            if raw in seen_rows:
                errors.append(f"{location}: duplicate GTF record")
            seen_rows.add(raw)
            owner = (path, cols[0], cols[6])
            if gene in owners and owners[gene] != owner:
                errors.append(f"{location}: duplicate gene_id {gene} across loci/files")
            owners[gene] = owner
            seq_genes.setdefault(cols[0], set()).add(gene)
            if (
                attrs.get("gene_biotype") == "whole_genome"
                or attrs.get("transcript_biotype") == "whole_genome"
            ):
                whole_genome.add(gene)
            if cols[2] in {"exon", "CDS", "transcript"}:
                identity = (path, cols[0], gene, cols[6])
                if tx in transcript_owners and transcript_owners[tx] != identity:
                    errors.append(f"{location}: inconsistent transcript {tx}")
                transcript_owners[tx] = identity
                needs_exon.add(tx)
                if cols[2] == "exon":
                    has_exon.add(tx)
        errors.extend(
            f"{path.name}: transcript {tx} has no exon" for tx in sorted(needs_exon - has_exon)
        )
    return seq_genes, whole_genome, errors


def catalogue_matrix(
    rows: list[dict[str, Any]], seq_genes: dict[str, set[str]], exclusions: set[str] | None = None
) -> list[dict[str, Any]]:
    """One row per shipped accession, including explicit excluded coverage gaps."""
    if not rows or any(not row.get("accession") for row in rows):
        raise ValueError("Catalogue needs nonempty accession rows")
    unknown = {row.get("panel", "") for row in rows} - PANELS
    if unknown:
        raise ValueError(f"Unknown catalogue panel value(s): {sorted(unknown)}")
    if not any(row["panel"] == "shipped" for row in rows):
        raise ValueError("Catalogue has no shipped rows; nothing would be checked")
    excluded = exclusions or set()
    by_name: dict[str, set[str]] = {}
    for row in rows:
        name = row.get("common_name") or row.get("species", "")
        if row.get("sibling_group"):
            by_name.setdefault(name, set()).add(row["sibling_group"])
    known = {row["accession"] for row in rows}
    extra = {seq.split(".")[0] for seq in seq_genes} - known
    result = []
    seen: set[str] = set()
    for row in rows:
        acc = row["accession"]
        if row.get("panel") != "shipped":
            continue
        errors = []
        if acc in seen:
            errors.append("duplicate_accession")
        seen.add(acc)
        genes = set().union(*(g for s, g in seq_genes.items() if s.split(".")[0] == acc))
        if not genes and acc not in excluded:
            errors.append("missing_gtf_seqname")
        if row.get("role", "") not in {
            "",
            "target",
            "decoy",
            "contaminant",
            "endogenous",
            "whole_genome",
        }:
            errors.append("invalid_role")
        name = row.get("common_name") or row.get("species", "")
        partner = SIBLING_VIRUS_PAIRS.get(name)
        if partner and (
            not row.get("sibling_group") or by_name.get(partner) != {row["sibling_group"]}
        ):
            errors.append("inconsistent_sibling_group")
        result.append(
            {
                "accession": acc,
                "panel": "shipped",
                "role": row.get("role", ""),
                "sibling_group": row.get("sibling_group", ""),
                "n_genes": len(genes),
                "status": "error" if errors else "excluded" if acc in excluded else "covered",
                "errors": "|".join(errors),
            }
        )
    result.extend(
        {
            "accession": acc,
            "panel": "",
            "role": "",
            "sibling_group": "",
            "n_genes": 0,
            "status": "error",
            "errors": "gtf_accession_not_catalogued",
        }
        for acc in sorted(extra)
    )
    return result


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=root / "src" / "viralscan" / "data")
    parser.add_argument("--catalogue", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--raw",
        action="store_true",
        help="gate on the shipped GTF files as they are, not the normalised view kb ref sees",
    )
    parser.add_argument("--ec-summary", type=Path)
    parser.add_argument("--max-ec-size", type=int)
    parser.add_argument("--max-discarded-ecs", type=int)
    args = parser.parse_args(argv)
    if not args.ec_summary and (args.max_ec_size is not None or args.max_discarded_ecs is not None):
        parser.error(
            "--max-ec-size/--max-discarded-ecs need --ec-summary; budgets would be ignored"
        )
    gtfs = sorted(args.data_dir.glob("*.gtf"))
    seq_genes, whole, errors = inspect_gtfs(gtfs, normalise=not args.raw)
    raw_errors = errors if args.raw else inspect_gtfs(gtfs)[2]
    # The builder emits these real gene-table exons at assembly, not as bundled GTFs.
    gene_table = args.data_dir / "anellovirus_genes.tsv"
    if gene_table.is_file():
        with gene_table.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                seq_genes.setdefault(row["accession"], set()).add(row["gene_id"])
    exclusions = args.data_dir / "index_exclusions.tsv"
    from viralscan.scripts.build_reference import load_index_exclusions

    excluded = set(load_index_exclusions(exclusions))
    matrix = catalogue_matrix(load_catalogue(args.catalogue), seq_genes, excluded)
    with args.output.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "accession",
                "panel",
                "role",
                "sibling_group",
                "n_genes",
                "status",
                "errors",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(matrix)
    if args.ec_summary:
        try:
            check_ec_budget(
                json.loads(args.ec_summary.read_text()), args.max_ec_size, args.max_discarded_ecs
            )
        except ValueError as exc:
            errors.append(str(exc))
    for error in errors:
        print(error)
    view = "raw" if args.raw else "normalised"
    print(f"shipped files as they are: {len(raw_errors)} errors (gate uses the {view} view)")
    print(
        f"{len(errors)} GTF/EC errors; {sum(r['status'] == 'error' for r in matrix)} catalogue gaps; "
        f"{len(whole)} explicitly tagged whole-genome genes"
    )
    return int(bool(errors or any(row["status"] == "error" for row in matrix)))


if __name__ == "__main__":
    raise SystemExit(main())

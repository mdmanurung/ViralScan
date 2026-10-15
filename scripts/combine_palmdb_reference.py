#!/usr/bin/env python3
"""Combine PalmDB RdRP proteins with ViralScan's annotated viral proteins.

Builds the single protein reference for the kallisto translated search: the
tutorial's PalmDB files (Luebbert et al., Nature Biotechnology 2025,
doi:10.1038/s41587-025-02614-y; PalmDB: Edgar et al., Nature 2022,
doi:10.1038/s41586-021-04332-2) plus the output of ``prepare_viral_proteins.py``.

Gene names: annotated proteins keep ``taxon_<taxid>``; PalmDB sequences map to
``palmdb_<rep_ID>`` (the tutorial's cluster representative). No sequence is
altered. Sequence identifiers must be unique; a PalmDB identifier that collides
with an annotated one gets the prefix ``palmdb_``; any other duplicate is fatal.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from viralscan.run_safety import sha256_file

ROOT = Path(__file__).resolve().parents[1]
RAW = "https://raw.githubusercontent.com/pachterlab/LSCHWCP_2023/main/PalmDB/"
# sha256 of the files the tutorial downloads (fetched 2026-10-14); see palmdb_provenance.json.
PALMDB_SHA256 = {
    "palmdb_rdrp_seqs.fa": "6a774cad43022f801e5dc1015a66344205d0d4c125c487a762706ea9687dcb09",
    "palmdb_clustered_t2g.txt": "3c4d468dcaacdbbd34bfa3ad502f9aac06def0d41f2b5c4a4f603622d68d5248",
    "ID_to_taxonomy_mapping.csv": "a9366dc0a4c271b59cef08643f2158eaa1fb296b250d9c689b6bcb20018e4ee8",
}
RANKS = ("phylum", "class", "order", "family", "genus", "species")


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in seqs:
                    raise ValueError(f"Duplicate sequence ID {name} in {path.name}")
                seqs[name] = ""
            elif name is None:
                if line.strip():
                    raise ValueError(f"Sequence before header in {path.name}")
            else:
                seqs[name] += line.strip()
    if any(not s for s in seqs.values()):
        raise ValueError(f"Empty sequence in {path.name}")
    return seqs


def read_tsv(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _norm(name: str) -> str:
    return " ".join(name.casefold().split())


def _label(row: dict[str, str], rep: str) -> str:
    """Third t2g column: lowest classified rank ('.' means unclassified), else the gene id."""
    return next((row[k] for k in reversed(RANKS) if row[k] not in ("", ".")), f"palmdb_{rep}")


def catalogue_index(catalogue: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Return (name -> taxids from species/organism, lineage term -> taxids)."""
    by_name: dict[str, set[str]] = defaultdict(set)
    by_lineage: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(catalogue):
        for key in ("species", "organism"):
            if row.get(key):
                by_name[_norm(row[key])].add(row["taxid"])
        for term in row.get("lineage", "").split(";"):
            if term.strip():
                by_lineage[_norm(term)].add(row["taxid"])
    return by_name, by_lineage


def build(
    annotated_dir: Path,
    palmdb_dir: Path,
    catalogue: Path,
    output: Path,
    expected: Path | None = None,
    blast_hits: Path | None = None,
    verify_palmdb: bool = True,
) -> dict:
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output already exists; use a new directory: {output}")
    palm = {n: palmdb_dir / n for n in PALMDB_SHA256}
    for name, path in palm.items():
        if verify_palmdb and sha256_file(path) != PALMDB_SHA256[name]:
            raise ValueError(f"{name} does not match the pinned tutorial file checksum")

    ann_seqs = read_fasta(annotated_dir / "proteins.faa")
    ann_t2g = [
        line.split("\t") for line in (annotated_dir / "viruses.t2g").read_text().splitlines()
    ]
    if {t[0] for t in ann_t2g} != set(ann_seqs) or len(ann_t2g) != len(ann_seqs):
        raise ValueError("Annotated proteins.faa and viruses.t2g do not list the same targets")
    pal_seqs = read_fasta(palm["palmdb_rdrp_seqs.fa"])
    pal_t2g = [
        line.split("\t") for line in palm["palmdb_clustered_t2g.txt"].read_text().splitlines()
    ]
    pal_rep = {t[0]: t[1] for t in pal_t2g}
    if len(pal_rep) != len(pal_t2g) or set(pal_rep) != set(pal_seqs):
        raise ValueError("PalmDB FASTA and clustered t2g do not list the same unique targets")
    tax_rows = read_tsv(palm["ID_to_taxonomy_mapping.csv"], ",")
    if any(pal_rep.get(r["ID"]) != r["rep_ID"] for r in tax_rows) or len(tax_rows) != len(pal_rep):
        raise ValueError("ID_to_taxonomy_mapping.csv disagrees with palmdb_clustered_t2g.txt")

    # Unique IDs across sources: prefix colliding PalmDB IDs, then re-check.
    pal_id = {i: (f"palmdb_{i}" if i in ann_seqs else i) for i in pal_seqs}
    final_ids = [*ann_seqs, *pal_id.values()]
    if len(set(final_ids)) != len(final_ids):
        raise ValueError("Sequence IDs are still not unique after prefixing")
    n_prefixed = sum(1 for i in pal_seqs if pal_id[i] != i)

    tax = {}  # rep_ID -> taxonomy row
    sizes: dict[str, int] = defaultdict(int)
    for r in tax_rows:
        sizes[r["rep_ID"]] += 1
        if r["ID"] == r["rep_ID"]:
            tax[r["rep_ID"]] = r
    if set(tax) != set(sizes):
        raise ValueError("A PalmDB cluster representative is missing from its own mapping")
    for r in tax_rows:  # members of a cluster share the representative's taxonomy
        if any(r[k] != tax[r["rep_ID"]][k] for k in RANKS):
            raise ValueError(f"{r['ID']} taxonomy differs from its representative {r['rep_ID']}")
    by_name, by_lineage = catalogue_index(catalogue)
    ann_taxids = {t[1].removeprefix("taxon_") for t in ann_t2g}
    cat_rows = {r["taxid"]: r for r in read_tsv(catalogue)}

    gene_rows = []
    for rep in sorted(tax, key=lambda x: (len(x), x)):
        row = tax[rep]
        name = _norm(row["species"])
        match, taxids = "none", set()
        if name in by_name:
            match, taxids = "species_exact", by_name[name]
        elif name in by_lineage:
            match, taxids = "lineage_term", by_lineage[name]
        gene_rows.append(
            {
                "gene": f"palmdb_{rep}",
                "rep_ID": rep,
                "n_sequences": sizes[rep],
                **{k: row[k] for k in RANKS},
                "strandedness": row["strandedness"],
                "catalogue_match": match,
                "catalogue_taxids": ";".join(sorted(taxids)),
                "catalogue_taxids_in_annotated": ";".join(sorted(taxids & ann_taxids)),
            }
        )

    overlap = {
        "annotated": {
            "sequences": len(ann_seqs),
            "virus_genes": len(ann_taxids),
        },
        "palmdb": {
            "sequences": len(pal_seqs),
            "clusters": len(tax),
            "species": len({r["species"] for r in tax.values() if r["species"] != "."}),
            "clusters_without_species": sum(1 for r in tax.values() if r["species"] == "."),
            "clusters_with_species": sum(1 for r in tax.values() if r["species"] != "."),
            "id_collisions_prefixed": n_prefixed,
        },
    }
    matched = [g for g in gene_rows if g["catalogue_taxids"]]
    exact = [g for g in gene_rows if g["catalogue_match"] == "species_exact"]
    overlap["palmdb_vs_catalogue"] = {
        "clusters_matching_catalogue_species_name": len(exact),
        "clusters_matching_catalogue_lineage_term_only": len(matched) - len(exact),
        "species_matching_catalogue_species_name": len({g["species"] for g in exact}),
        "catalogue_taxids_matched_by_species_name": len(
            {t for g in exact for t in g["catalogue_taxids"].split(";")}
        ),
        "clusters_matching_an_annotated_virus": sum(
            1 for g in matched if g["catalogue_taxids_in_annotated"]
        ),
    }
    pal_set = set(pal_seqs.values())
    shared_exact = set(ann_seqs.values()) & pal_set
    overlap["exact_sequence_identity"] = {
        "annotated_proteins_identical_to_a_palmdb_sequence": sum(
            1 for s in ann_seqs.values() if s in pal_set
        ),
        "distinct_sequences_in_common": len(shared_exact),
    }
    if blast_hits is not None:
        overlap["blast_duplicates"] = _blast_duplicates(blast_hits, ann_seqs, pal_seqs, pal_rep)
    if expected is not None:
        presence = _expected_presence(expected, cat_rows, ann_taxids, tax)
    else:
        presence = []

    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    try:
        with (stage / "combined_proteins.faa").open("w") as handle:
            for i, s in ann_seqs.items():
                handle.write(f">{i}\n{s}\n")
            for i, s in pal_seqs.items():
                handle.write(f">{pal_id[i]}\n{s}\n")
        with (stage / "combined.t2g").open("w") as handle:
            for t in ann_t2g:
                handle.write("\t".join(t) + "\n")
            for i in pal_seqs:
                rep = pal_rep[i]
                handle.write(f"{pal_id[i]}\tpalmdb_{rep}\t{_label(tax[rep], rep)}\n")
        _write_tsv(stage / "palmdb_taxonomy.tsv", gene_rows)
        # proteins.tsv in the layout summarize_translated.py reads (target, virus_id, species, product, gene)
        ann_cols = list(read_tsv(annotated_dir / "proteins.tsv")[0])
        prot_rows = read_tsv(annotated_dir / "proteins.tsv") + [
            {
                **dict.fromkeys(ann_cols, ""),
                "target": pal_id[i],
                "protein_id": i,
                "virus_id": f"palmdb_{pal_rep[i]}",
                "species": _label(tax[pal_rep[i]], pal_rep[i]),
                "family": tax[pal_rep[i]]["family"],
                "length_aa": len(pal_seqs[i]),
            }
            for i in pal_seqs
        ]
        _write_tsv(stage / "proteins.tsv", prot_rows)
        if presence:
            _write_tsv(stage / "expected_virus_presence.tsv", presence)
        receipt = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "palmdb_files": {
                n: {"url": RAW + n, "sha256": sha256_file(p)} for n, p in palm.items()
            },
            "annotated_sha256": {
                n: sha256_file(annotated_dir / n) for n in ("proteins.faa", "viruses.t2g")
            },
            "catalogue_sha256": sha256_file(catalogue),
            "script_sha256": sha256_file(Path(__file__)),
            "overlap": overlap,
            "expected_virus_presence": presence,
        }
        receipt["output_sha256"] = {p.name: sha256_file(p) for p in sorted(stage.iterdir())}
        (stage / "reference.json").write_text(json.dumps(receipt, indent=2) + "\n")
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return receipt


def _write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _blast_duplicates(hits: Path, ann_seqs, pal_seqs, pal_rep, min_pident=90.0, min_cov=0.8):
    """blastp outfmt ``6 qseqid sseqid pident length qlen`` of PalmDB (query) vs annotated (subject)."""
    ann_hit, pal_hit = set(), set()
    with hits.open() as handle:
        for line in handle:
            q, s, pident, length, qlen = line.split("\t")[:5]
            if float(pident) >= min_pident and int(length) >= min_cov * int(qlen):
                if q not in pal_seqs or s not in ann_seqs:
                    raise ValueError(f"BLAST hit names an unknown sequence: {q} {s}")
                pal_hit.add(q)
                ann_hit.add(s)
    return {
        "criterion": f"pident >= {min_pident}, aligned length >= {min_cov} x PalmDB length",
        "annotated_proteins_with_a_palmdb_match": len(ann_hit),
        "palmdb_sequences_with_an_annotated_match": len(pal_hit),
        "palmdb_clusters_with_an_annotated_match": len({pal_rep[q] for q in pal_hit}),
    }


def _expected_presence(expected: Path, cat_rows, ann_taxids, tax) -> list[dict]:
    """One row per expected virus: listed taxon ids plus catalogue species matching name_pattern."""
    out = []
    for row in read_tsv(expected):
        ids = {t for t in row["taxon_ids"].split(";") if t}
        if row["name_pattern"]:
            cat_pat = re.compile(row["name_pattern"], re.IGNORECASE)
            ids |= {t for t, r in cat_rows.items() if cat_pat.search(r["species"])}
        pal_pat = row["palmdb_pattern"] or row["name_pattern"]
        pattern = re.compile(pal_pat, re.IGNORECASE) if pal_pat else None
        pal = [
            g
            for g, r in tax.items()
            if pattern and any(pattern.search(r[k]) for k in ("species", "genus", "family"))
        ]
        out.append(
            {
                "name": row["name"],
                "catalogue_taxids": ";".join(sorted(ids & set(cat_rows))),
                "annotated_taxids_with_proteins": ";".join(sorted(ids & ann_taxids)),
                "palmdb_clusters_by_name": len(pal),
                "palmdb_genes": ";".join(f"palmdb_{g}" for g in sorted(pal)[:5]),
                "palmdb_species_by_name": ";".join(sorted({tax[g]["species"] for g in pal})[:5]),
            }
        )
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotated-dir", type=Path, required=True)
    ap.add_argument("--palmdb-dir", type=Path, required=True)
    ap.add_argument("--catalogue", type=Path, default=ROOT / "src/viralscan/data/virus_catalog.tsv")
    ap.add_argument(
        "--expected", type=Path, help="TSV: name, taxon_ids, name_pattern, palmdb_pattern"
    )
    ap.add_argument(
        "--blast-hits", type=Path, help="blastp -outfmt '6 qseqid sseqid pident length qlen'"
    )
    ap.add_argument("--output", type=Path, help="new directory for the combined reference")
    ap.add_argument(
        "--blast-only",
        type=Path,
        metavar="JSON",
        help="only write the BLAST duplicate summary of --blast-hits to JSON (no reference)",
    )
    args = ap.parse_args(argv)
    if args.blast_only:
        ann = read_fasta(args.annotated_dir / "proteins.faa")
        pal = read_fasta(args.palmdb_dir / "palmdb_rdrp_seqs.fa")
        rep = dict(
            line.split("\t")
            for line in (args.palmdb_dir / "palmdb_clustered_t2g.txt").read_text().splitlines()
        )
        result = _blast_duplicates(args.blast_hits, ann, pal, rep)
        args.blast_only.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
        return 0
    if args.output is None:
        ap.error("--output is required")
    receipt = build(
        args.annotated_dir,
        args.palmdb_dir,
        args.catalogue,
        args.output,
        args.expected,
        args.blast_hits,
    )
    print(json.dumps(receipt["overlap"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

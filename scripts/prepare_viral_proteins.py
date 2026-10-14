#!/usr/bin/env python3
"""Prepare annotated proteins of local GTF viruses for kallisto translated search.

Method: Laura Luebbert's kallisto single-cell virus tutorial and Luebbert et al.,
Nature Biotechnology (2025), doi:10.1038/s41587-025-02614-y. Uses GenBank CDS
translations instead of PalmDB RdRPs; no alignment algorithm is implemented here.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import Bio
from Bio import SeqIO

from viralscan.run_safety import sha256_file
from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR, genbank_cache_path

ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = "https://kallisto.readthedocs.io/en/latest/translated/notebooks/virus_detection_sc.html"
METHOD = "https://doi.org/10.1038/s41587-025-02614-y"


def prepare(data_dir: Path, catalogue: Path, cache_dir: Path, output: Path) -> dict:
    # Snakemake creates output parent directories before starting the rule.
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output already exists; use a new directory: {output}")
    with catalogue.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_accession = {r["accession"]: r for r in rows}
    if len(by_accession) != len(rows):
        raise ValueError("Duplicate catalogue accessions")
    gtfs = sorted(data_dir.glob("*.gtf"))
    if not gtfs:
        raise ValueError(f"No GTFs in {data_dir}")
    coverage, proteins, sequences = [], [], {}
    for gtf in gtfs:
        with gtf.open() as handle:
            seqnames = {
                line.split("\t")[0] for line in handle if line.strip() and not line.startswith("#")
            }
        if len(seqnames) != 1:
            raise ValueError(f"Expected exactly one genome accession in {gtf}")
        seqname = seqnames.pop()
        accession = re.sub(r"\.\d+$", "", seqname)
        row = by_accession[accession]
        candidates = [
            genbank_cache_path(a, cache_dir)
            for a in dict.fromkeys([seqname, row["accession_version"], accession])
        ]
        gb = next(
            (p for p in candidates if p.is_file() and sha256_file(p) == row["gb_sha256"]), None
        )
        if gb is None:
            raise ValueError(f"No catalogue-hash-matched GenBank cache for {accession}")
        record = SeqIO.read(gb, "genbank")
        if record.id != row["accession_version"] or ("." in seqname and seqname != record.id):
            raise ValueError(f"GTF/catalogue/GenBank version mismatch: {gtf}, {record.id}")
        cds = [feature for feature in record.features if feature.type == "CDS"]
        n_proteins = 0
        for feature in cds:
            qualifiers = feature.qualifiers
            if "translation" not in qualifiers:
                if "pseudo" in qualifiers or "pseudogene" in qualifiers:
                    continue
                raise ValueError(f"Untranslated non-pseudo CDS in {record.id}: {feature.location}")
            aa = qualifiers["translation"][0]
            protein_id = qualifiers.get("protein_id", [""])[0]
            if not protein_id or not aa or set(aa) - set("ACDEFGHIKLMNPQRSTVWYXBZUOJ"):
                raise ValueError(f"Invalid protein identifier/translation in {record.id}")
            target = f"{record.id}__{protein_id}"
            if target in sequences:
                raise ValueError(f"Duplicate protein target: {target}")
            sequences[target] = aa
            proteins.append(
                {
                    "target": target,
                    "protein_id": protein_id,
                    "accession": record.id,
                    "virus_id": f"taxon_{row['taxid']}",
                    "species": row["species"],
                    "family": row["family"],
                    "gene": qualifiers.get("gene", [""])[0],
                    "product": qualifiers.get("product", [""])[0],
                    "length_aa": len(aa),
                    "location": str(feature.location),
                    "genbank_sha256": row["gb_sha256"],
                }
            )
            n_proteins += 1
        coverage.append(
            {
                "gtf": gtf.name,
                "gtf_sha256": sha256_file(gtf),
                "accession": record.id,
                "species": row["species"],
                "taxid": row["taxid"],
                "family": row["family"],
                "virus_genome": "RNA"
                if "RNA" in row["molecule"].upper() or row["family"] == "Retroviridae"
                else "DNA",
                "reference_molecule": row["molecule"],
                "genbank_sha256": row["gb_sha256"],
                "cds_features": len(cds),
                "proteins": n_proteins,
                "status": "annotated_proteins" if n_proteins else "no_annotated_translation",
            }
        )
    if not proteins:
        raise ValueError("No annotated protein translations")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    try:
        with (stage / "proteins.faa").open("w") as handle:
            for target, aa in sequences.items():
                handle.write(f">{target}\n{aa}\n")
        for name, table in [("proteins.tsv", proteins), ("coverage.tsv", coverage)]:
            with (stage / name).open("w") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(table[0]), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(table)
        with (stage / "viruses.t2g").open("w") as handle:
            for protein in proteins:
                handle.write(f"{protein['target']}\t{protein['virus_id']}\t{protein['species']}\n")
        receipt = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "method": METHOD,
            "tutorial": TUTORIAL,
            "tutorial_author": "Laura Luebbert",
            "adaptation": "GenBank annotated CDS proteins of GTF viruses instead of PalmDB RdRPs; no sequence clustering",
            "catalogue_sha256": sha256_file(catalogue),
            "script_sha256": sha256_file(Path(__file__)),
            "biopython_version": Bio.__version__,
            "gtfs": len(gtfs),
            "proteins": len(proteins),
            "virus_labels": len({r["species"] for r in coverage}),
            "covered_virus_labels": len({r["species"] for r in coverage if r["proteins"]}),
            "gaps": [r["accession"] for r in coverage if not r["proteins"]],
            "retrieval": "Existing GenBank cache verified against catalogue hashes; no downloads",
            "output_sha256": {p.name: sha256_file(p) for p in sorted(stage.iterdir())},
        }
        (stage / "reference.json").write_text(json.dumps(receipt, indent=2) + "\n")
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return receipt


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "src/viralscan/data")
    parser.add_argument(
        "--catalogue", type=Path, default=ROOT / "src/viralscan/data/virus_catalog.tsv"
    )
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    receipt = prepare(args.data_dir, args.catalogue, args.cache_dir, args.output)
    print(
        json.dumps(
            {
                k: receipt[k]
                for k in ["gtfs", "proteins", "virus_labels", "covered_virus_labels", "gaps"]
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

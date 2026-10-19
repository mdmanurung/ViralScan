import csv
import importlib.util
from pathlib import Path

import pytest

Bio = pytest.importorskip("Bio")
from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqFeature import FeatureLocation, SeqFeature  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from viralscan.run_safety import sha256_file  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/prepare_viral_proteins.py"
spec = importlib.util.spec_from_file_location("prepare_viral_proteins", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def inputs(tmp_path):
    data, cache = tmp_path / "data", tmp_path / "cache"
    data.mkdir()
    rows = []
    for number, (family, molecule, translated) in enumerate(
        [
            ("Retroviridae", "DNA", True),
            ("Orthoherpesviridae", "DNA", True),
            ("Flaviviridae", "RNA", False),
        ],
        1,
    ):
        accession = f"NC_{number:06d}"
        record = SeqRecord(Seq("ATG" * 90), id=f"{accession}.1", name=accession, description="test")
        record.annotations["molecule_type"] = "DNA"
        if translated:
            record.features = [
                SeqFeature(
                    FeatureLocation(0, 270),
                    type="CDS",
                    qualifiers={
                        "translation": ["M" * 90],
                        "protein_id": [f"NP_{number:06d}.1"],
                        "product": ["test protein"],
                    },
                )
            ]
        directory = cache / record.id
        directory.mkdir(parents=True)
        gb = directory / f"{record.id}.gb"
        SeqIO.write(record, gb, "genbank")
        (data / f"{accession}.gtf").write_text(
            f'{record.id}\tfixture\texon\t1\t270\t.\t+\t.\tgene_id "test";\n'
        )
        rows.append(
            {
                "accession": accession,
                "accession_version": record.id,
                "species": f"Virus {number}",
                "taxid": str(number),
                "family": family,
                "molecule": molecule,
                "gb_sha256": sha256_file(gb),
            }
        )
    catalogue = tmp_path / "catalogue.tsv"
    with catalogue.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return data, catalogue, cache


def test_full_translations_dna_retrovirus_and_annotation_gap(tmp_path):
    data, catalogue, cache = inputs(tmp_path)
    output = tmp_path / "reference"
    output.mkdir()  # Snakemake creates output parents before executing rules.
    receipt = m.prepare(data, catalogue, cache, output)
    records = list(SeqIO.parse(output / "proteins.faa", "fasta"))
    assert [str(r.seq) for r in records] == ["M" * 90, "M" * 90]
    assert receipt["proteins"] == 2 and receipt["covered_virus_labels"] == 2
    assert receipt["gaps"] == ["NC_000003.1"]
    coverage = list(csv.DictReader((output / "coverage.tsv").open(), delimiter="\t"))
    assert [r["virus_genome"] for r in coverage] == ["RNA", "DNA", "RNA"]
    assert coverage[-1]["status"] == "no_annotated_translation"
    targets = [line.split("\t")[0] for line in (output / "viruses.t2g").read_text().splitlines()]
    assert targets == [r.id for r in records]
    assert all(
        sha256_file(output / name) == digest for name, digest in receipt["output_sha256"].items()
    )


def test_cache_corruption_fails_without_partial_output(tmp_path):
    data, catalogue, cache = inputs(tmp_path)
    gb = cache / "NC_000001.1/NC_000001.1.gb"
    gb.write_text(gb.read_text() + "corrupted\n")
    with pytest.raises(ValueError, match="hash-matched"):
        m.prepare(data, catalogue, cache, tmp_path / "reference")
    assert not (tmp_path / "reference").exists()


def test_gtf_version_mismatch_fails(tmp_path):
    data, catalogue, cache = inputs(tmp_path)
    gtf = data / "NC_000001.gtf"
    gtf.write_text(gtf.read_text().replace("NC_000001.1", "NC_000001.2"))
    with pytest.raises(ValueError, match="version mismatch"):
        m.prepare(data, catalogue, cache, tmp_path / "reference")


def test_existing_reference_is_preserved(tmp_path):
    data, catalogue, cache = inputs(tmp_path)
    output = tmp_path / "reference"
    output.mkdir()
    (output / "sentinel").write_text("keep")
    with pytest.raises(ValueError, match="already exists"):
        m.prepare(data, catalogue, cache, output)
    assert (output / "sentinel").read_text() == "keep"

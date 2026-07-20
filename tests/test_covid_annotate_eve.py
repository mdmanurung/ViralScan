"""Regression tests for the COVID EVE annotation helper."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "covid_viralscan/scripts/annotate_eve.py"
SPEC = importlib.util.spec_from_file_location("covid_annotate_eve", SCRIPT_PATH)
assert SPEC is not None
annotate_eve = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(annotate_eve)


def test_annotate_locus_finds_long_gene_with_many_internal_starts() -> None:
    records = [(1, 1000, "LONG", "protein_coding")]
    records.extend((start, start + 5, f"SHORT{start}", "lncRNA") for start in range(100, 500, 10))
    genes_by_chrom = {"chr1": sorted(records)}

    assert annotate_eve.annotate_locus("chr1", 750, genes_by_chrom) == (
        "in_gene:LONG",
        "protein_coding",
        0,
    )


def test_phase_a_reports_distinct_loci_and_filters_key_accessions(tmp_path: Path) -> None:
    phase_a = tmp_path / "reads"
    outdir = tmp_path / "annotation"
    phase_a.mkdir()
    outdir.mkdir()

    (phase_a / "NC_001_x213_grch38.depth.txt").write_text(
        "\n".join(
            [
                "chr1\t120\t3",
                "chr1\t121\t4",
                "chr1\t1005\t5",
                "chr1\t1006\t2",
            ]
        )
        + "\n"
    )
    (phase_a / "NC_001_x213_grch38.idxstats.txt").write_text("chr1\t2000\t9\t0\n")
    (phase_a / "NC_999_x213_grch38.depth.txt").write_text("chr1\t150\t10\n")

    genes_by_chrom = {
        "chr1": [
            (100, 200, "GENEA", "protein_coding"),
            (1000, 1100, "GENEB", "lncRNA"),
        ]
    }

    rows = annotate_eve.annotate_phase_a(str(phase_a), ["NC_001"], genes_by_chrom, str(outdir))

    assert len(rows) == 2
    assert {row["accession"] for row in rows} == {"NC_001"}
    assert {(row["human_start"], row["human_end"]) for row in rows} == {
        (120, 121),
        (1005, 1006),
    }
    assert {row["gene_annotation"] for row in rows} == {"in_gene:GENEA", "in_gene:GENEB"}

    with (outdir / "phase_a_loci.tsv").open() as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["human_start"] for row in written] == ["1005", "120"]


def test_phase_b_empty_blast_still_writes_normalized_tsv(tmp_path: Path) -> None:
    blast = tmp_path / "covered_vs_nt_human.tsv"
    outdir = tmp_path / "annotation"
    blast.write_text("")
    outdir.mkdir()

    rows = annotate_eve.annotate_phase_b(str(blast), {}, str(outdir))

    assert rows == []
    with (outdir / "phase_b_blast.tsv").open() as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert written == []


def test_phase_c_reports_one_based_inclusive_coordinates(tmp_path: Path) -> None:
    paf = tmp_path / "panel_vs_grch38.paf"
    outdir = tmp_path / "annotation"
    paf.write_text("ACC\t200\t0\t50\t+\tchr1\t1000\t100\t150\t50\t50\t20\n")
    outdir.mkdir()

    genes_by_chrom = {"chr1": [(125, 125, "MID", "protein_coding")]}

    rows = annotate_eve.annotate_phase_c(str(paf), genes_by_chrom, str(outdir))

    assert rows[0]["best_coords"] == "101-150"
    assert rows[0]["best_gene_annotation"] == "in_gene:MID"
    with (outdir / "phase_c_panel.tsv").open() as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert written[0]["best_coords"] == "101-150"

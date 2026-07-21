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
    # GTF keys are normalized (load_gtf_genes strips the 'chr' prefix); a UCSC-style
    # query must still match after normalization.
    genes_by_chrom = {"1": sorted(records)}

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
        "1": [
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

    genes_by_chrom = {"1": [(125, 125, "MID", "protein_coding")]}

    rows = annotate_eve.annotate_phase_c(str(paf), genes_by_chrom, str(outdir))

    assert rows[0]["best_coords"] == "101-150"
    assert rows[0]["best_gene_annotation"] == "in_gene:MID"
    with (outdir / "phase_c_panel.tsv").open() as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert written[0]["best_coords"] == "101-150"


def test_normalize_chrom_folds_chr_prefix_and_mito() -> None:
    assert annotate_eve._normalize_chrom("chr7") == annotate_eve._normalize_chrom("7") == "7"
    assert annotate_eve._normalize_chrom("chrM") == annotate_eve._normalize_chrom("MT") == "MT"


def test_is_chromosome_subject_distinguishes_chromosomes_from_scaffolds() -> None:
    assert annotate_eve._is_chromosome_subject("NC_000007.14") is True
    assert annotate_eve._is_chromosome_subject("NC_012920.1") is True  # mito
    assert annotate_eve._is_chromosome_subject("AC_012345.2") is False  # BAC clone
    assert annotate_eve._is_chromosome_subject("NT_187513.1") is False  # scaffold
    assert annotate_eve._is_chromosome_subject("gi|568336|ref|NC_000007.14|") is True  # legacy id
    assert annotate_eve._is_chromosome_subject("NC_045512.2") is False  # viral, not a chromosome
    assert annotate_eve._is_chromosome_subject("NC_000024.10") is True  # chrY (last human)
    assert annotate_eve._is_chromosome_subject("NC_000067.6") is False  # mouse chr1, not human
    assert annotate_eve._is_chromosome_subject("NC_000025.1") is False  # beyond chr24


def test_phase_b_annotates_chromosome_subject_but_flags_clone(tmp_path: Path) -> None:
    """Chromosome subjects get a gene; clone/scaffold subjects (subject-local
    coordinates) must NOT be looked up against chromosome GTF intervals."""
    blast = tmp_path / "covered_vs_nt_human.tsv"
    outdir = tmp_path / "annotation"
    outdir.mkdir()
    # cols: qseqid sseqid stitle pident length qlen qstart qend sstart send evalue bitscore
    blast.write_text(
        "\n".join(
            [
                # whole-chromosome subject; sstart/send are genomic -> mid 150 -> GENEZ
                "EBV_hit\tNC_000001.11\tHomo sapiens chromosome 1, GRCh38.p14\t"
                "95.0\t100\t100\t1\t100\t120\t180\t1e-40\t200",
                # BAC clone subject; sstart/send are clone-local -> must be flagged
                "EBV_clone\tAC_012345.2\tHomo sapiens BAC clone RP11 from chromosome 1\t"
                "92.0\t100\t100\t1\t100\t45000\t45100\t1e-30\t150",
            ]
        )
        + "\n"
    )
    genes_by_chrom = {"1": [(100, 200, "GENEZ", "protein_coding")]}

    rows = annotate_eve.annotate_phase_b(str(blast), genes_by_chrom, str(outdir))

    by_q = {r["qseqid"]: r for r in rows}
    assert by_q["EBV_hit"]["gene_annotation"] == "in_gene:GENEZ"
    assert by_q["EBV_clone"]["gene_annotation"] == "subject_not_chromosome"


def test_namespace_mismatch_warns(tmp_path: Path, capsys) -> None:
    """A GTF whose chromosomes never match the query namespace triggers a warning
    rather than silently annotating everything 'intergenic'."""
    paf = tmp_path / "panel_vs_grch38.paf"
    outdir = tmp_path / "annotation"
    outdir.mkdir()
    # query target 'scaffold_99' never matches the GTF's chromosome '1'
    paf.write_text("ACC\t200\t0\t50\t+\tscaffold_99\t1000\t100\t150\t50\t50\t20\n")

    annotate_eve.annotate_phase_c(str(paf), {"1": [(1, 2, "G", "protein_coding")]}, str(outdir))

    assert "WARNING [Phase C]" in capsys.readouterr().err

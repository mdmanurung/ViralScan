from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from tests._fastq import write_fastq
from viralscan.reference_strategy import (
    PANEL_ID,
    BenchmarkContractError,
    audit_fastqs,
    audit_manifest,
    sha256_file,
    validate_audit_rows,
    validate_commands,
    validate_fastq_audit_rows,
    validate_manifest,
    write_commands,
)


def _manifest(tmp_path: Path) -> dict:
    virus_fa = tmp_path / "viral_panel.fa"
    virus_gtf = tmp_path / "viral_panel.gtf"
    virus_fa.write_text(">EBV\nACGT\n>HHV-6B\nACGT\n>HSV-1\nACGT\n", encoding="utf-8")
    virus_gtf.write_text(
        'EBV\tViralScan\texon\t1\t4\t.\t+\t.\tgene_id "Epstein_gene";\n'
        'HHV6B\tViralScan\texon\t1\t4\t.\t+\t.\tgene_id "human_herpesvirus_6b_gene";\n'
        'HSV1\tViralScan\texon\t1\t4\t.\t+\t.\tgene_id "human_herpesvirus_1_gene";\n',
        encoding="utf-8",
    )
    fastq_root = tmp_path / "fastq"
    fastq_root.mkdir()
    for dataset in ("hhv6_carT_ref", "lcl_5lines", "hsv1_fibroblast"):
        (fastq_root / dataset).mkdir()
    for srr, dataset in (
        ("SRR20710641", "hhv6_carT_ref"),
        ("SRR12682296", "lcl_5lines"),
        ("SRR8315713", "hsv1_fibroblast"),
    ):
        srr_dir = fastq_root / dataset / srr
        srr_dir.mkdir()
        write_fastq(srr_dir / f"{srr}_1.fastq.gz", [("read1", "ACGT")])
        write_fastq(srr_dir / f"{srr}_2.fastq.gz", [("read1", "TGCA")])

    star_human = tmp_path / "star_human"
    star_virus = tmp_path / "star_virus"
    star_combined = tmp_path / "star_combined"
    for directory in (star_human, star_virus, star_combined):
        directory.mkdir()
        for star_file in ("Genome", "SA", "SAindex", "genomeParameters.txt"):
            (directory / star_file).write_text("fixture\n", encoding="utf-8")
    for filename in (
        "human.fa",
        "human.gtf",
        "human_cdna.fa",
        "human.idx",
        "virus.idx",
        "virus.t2g",
        "combined.idx",
        "combined.t2g",
        "combined.gtf",
        "anello.tsv",
    ):
        (tmp_path / filename).write_text("fixture\n", encoding="utf-8")

    return {
        "panel": PANEL_ID,
        "fastq_root": str(fastq_root),
        "fastqs": {
            "SRR20710641": {"source_url": "https://example.org/SRR20710641"},
            "SRR12682296": {"source_url": "https://example.org/SRR12682296"},
            "SRR8315713": {"source_url": "https://example.org/SRR8315713"},
        },
        "build_commands": {
            "starsolo_human_only": "STAR --runMode genomeGenerate ...",
            "starsolo_all_virus": "STAR --runMode genomeGenerate ...",
            "starsolo_combined": "STAR --runMode genomeGenerate ...",
            "kallisto_human_only": "kb ref ...",
            "kallisto_all_virus": "kb ref ...",
            "kallisto_combined": "kb ref ...",
        },
        "human": {
            "source_release": "GRCh38-2024-A",
            "genome_fasta": str(tmp_path / "human.fa"),
            "genes_gtf": str(tmp_path / "human.gtf"),
            "human_cdna": str(tmp_path / "human_cdna.fa"),
        },
        "provenance": {"note": "fixture"},
        "viral_panel": {
            "id": PANEL_ID,
            "anellovirus_accession_table": str(tmp_path / "anello.tsv"),
            "anellovirus_expected_count": 3,
            "anellovirus_fetched_count": 3,
            "anellovirus_missing_count": 0,
        },
        "references": {
            "starsolo": {
                "human_source_release": "GRCh38-2024-A",
                "human_only": {"genome_dir": str(star_human)},
                "all_virus": {
                    "genome_fasta": str(virus_fa),
                    "genome_gtf": str(virus_gtf),
                    "genome_dir": str(star_virus),
                },
                "combined": {"genome_dir": str(star_combined)},
            },
            "viralscan": {
                "human_source_release": "GRCh38-2024-A",
                "human_only": {"kallisto_index": str(tmp_path / "human.idx")},
                "all_virus": {
                    "kallisto_index": str(tmp_path / "virus.idx"),
                    "t2g": str(tmp_path / "virus.t2g"),
                    "gtf": str(virus_gtf),
                },
                "combined": {
                    "kallisto_index": str(tmp_path / "combined.idx"),
                    "t2g": str(tmp_path / "combined.t2g"),
                    "gtf": str(tmp_path / "combined.gtf"),
                },
            },
        },
    }


def test_validate_manifest_rejects_mismatched_human_source(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    manifest["references"]["viralscan"]["human_source_release"] = "GRCh37"
    with pytest.raises(BenchmarkContractError, match="human source releases"):
        validate_manifest(manifest)


def test_validate_manifest_rejects_single_virus_reference(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    manifest["references"]["starsolo"]["combined"]["genome_dir"] = "starsolo_p22_6/GRCh38_EBV"
    with pytest.raises(BenchmarkContractError, match="single-virus"):
        validate_manifest(manifest, fail_on_single_virus=True)


def test_validate_manifest_requires_provenance_fields(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    manifest["viral_panel"].pop("anellovirus_fetched_count")
    with pytest.raises(BenchmarkContractError, match="provenance"):
        validate_manifest(manifest, require_provenance=True)


def test_audit_manifest_covers_schema_path_fields(tmp_path: Path) -> None:
    rows = audit_manifest(_manifest(tmp_path))
    keys = {row.key for row in rows}
    assert "references.viralscan.all_virus.gtf" in keys
    assert "references.viralscan.all_virus.t2g" in keys
    assert "references.viralscan.combined.gtf" in keys
    assert "viral_panel.anellovirus_accession_table" in keys


def test_audit_manifest_records_hash_features_and_all_expected_presence(tmp_path: Path) -> None:
    rows = audit_manifest(_manifest(tmp_path))
    validate_audit_rows(rows)
    gtf_row = next(row for row in rows if row.key.endswith("genome_gtf"))
    assert gtf_row.exists is True
    assert gtf_row.feature_count == 3
    assert set(gtf_row.expected_virus_presence.split(",")) == {"HHV-6B", "EBV", "HSV-1"}
    assert len(gtf_row.sha256) == 64


def test_sha256_file_hashes_gzip_artifact_bytes(tmp_path: Path) -> None:
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@read1\nACGT\n+\nIIII\n")

    assert sha256_file(path) == hashlib.sha256(path.read_bytes()).hexdigest()


def test_audit_validation_rejects_missing_path(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    Path(manifest["references"]["viralscan"]["all_virus"]["t2g"]).unlink()
    rows = audit_manifest(manifest)
    with pytest.raises(BenchmarkContractError, match="manifest paths are missing"):
        validate_audit_rows(rows)


def test_audit_validation_rejects_incomplete_star_genome_dir(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    Path(manifest["references"]["starsolo"]["all_virus"]["genome_dir"], "SA").unlink()
    rows = audit_manifest(manifest)
    with pytest.raises(BenchmarkContractError, match="incomplete STAR genomeGenerate directories"):
        validate_audit_rows(rows)


def test_audit_validation_rejects_missing_target_presence(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    Path(manifest["references"]["starsolo"]["all_virus"]["genome_fasta"]).write_text(
        ">EBV\nACGT\n",
        encoding="utf-8",
    )
    Path(manifest["references"]["starsolo"]["all_virus"]["genome_gtf"]).write_text(
        'EBV\tViralScan\texon\t1\t4\t.\t+\t.\tgene_id "Epstein_gene";\n',
        encoding="utf-8",
    )
    rows = audit_manifest(manifest)
    with pytest.raises(BenchmarkContractError, match="expected viruses absent"):
        validate_audit_rows(rows)


def test_write_commands_has_all_12_rows_and_no_forbidden_paths(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    run_dir = tmp_path / "run"
    commands = write_commands(run_dir, manifest)
    assert len(commands) == 12
    validate_commands(run_dir / "commands.jsonl")
    records = [json.loads(line) for line in (run_dir / "commands.jsonl").read_text().splitlines()]
    assert {record["target_virus"] for record in records} == {"HHV-6B", "EBV", "HSV-1"}
    assert any(
        "--host-filter" in record["command"]
        for record in records
        if record["method"] == "viralscan"
    )
    star_two_step = next(
        record
        for record in records
        if record["method"] == "starsolo" and record["reference_strategy"] == "two_step"
    )
    assert len(star_two_step["commands"]) == 2
    assert "Unmapped.out.mate1" in " ".join(star_two_step["commands"][1])
    assert "--soloCBlen" in star_two_step["commands"][0]


def test_fastq_audit_records_hashes_and_source_urls(tmp_path: Path) -> None:
    rows = audit_fastqs(_manifest(tmp_path))
    validate_fastq_audit_rows(rows)
    assert len(rows) == 6
    assert all(row.sha256 for row in rows)
    assert all(row.source_url.startswith("https://example.org/") for row in rows)


def test_fastq_audit_rejects_malformed_sequence_quality_lengths(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    bad = Path(manifest["fastq_root"]) / "lcl_5lines" / "SRR12682296" / "SRR12682296_1.fastq.gz"
    with gzip.open(bad, "wt") as handle:
        handle.write("@bad\nACGT\n+\nIII\n")
    rows = audit_fastqs(manifest)
    with pytest.raises(BenchmarkContractError, match="quality string length"):
        validate_fastq_audit_rows(rows)

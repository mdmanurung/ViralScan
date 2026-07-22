"""Regression tests for exact fragment identity in v3 host filtering."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from tests._fastq import write_fastq
from viralscan.scripts import host_filter
from viralscan.scripts.host_filter import (
    canonical_read_id,
    check_host_filter_tools,
    filter_fastq_pairs,
    validate_paired_fastq_ids,
)


def test_host_filter_script_compiles_after_snakemake_preamble() -> None:
    script = Path(host_filter.__file__).read_text(encoding="utf-8")
    compile("snakemake = None\n" + script, str(host_filter.__file__), "exec")


class TestExactReadIdentity:
    def test_canonicalizes_common_mate_headers(self) -> None:
        assert canonical_read_id("@read-1/1 extra") == "read-1"
        assert canonical_read_id("@read-1/2 extra") == "read-1"

    def test_validates_plain_and_gzipped_synchronized_pairs(self, tmp_path: Path) -> None:
        barcode_umi = "A" * 28
        r1 = tmp_path / "R1.fastq.gz"
        r2 = tmp_path / "R2.fastq.gz"
        # The two records deliberately share a CB–UMI. V3 treats them as two
        # fragments because exact read IDs, not CB–UMI tuples, define filtering.
        write_fastq(r1, [("fragment-a/1", barcode_umi), ("fragment-b/1", barcode_umi)])
        write_fastq(r2, [("fragment-a/2", "ACGT"), ("fragment-b/2", "TGCA")])
        assert validate_paired_fastq_ids(str(r1), str(r2)) == 2

    def test_rejects_mate_mismatch(self, tmp_path: Path) -> None:
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        write_fastq(r1, [("fragment-a", "AAAA")])
        write_fastq(r2, [("fragment-b", "TTTT")])
        with pytest.raises(ValueError, match="FASTQ mate mismatch"):
            validate_paired_fastq_ids(str(r1), str(r2))

    def test_rejects_different_record_counts(self, tmp_path: Path) -> None:
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        write_fastq(r1, [("a", "AAAA"), ("b", "CCCC")])
        write_fastq(r2, [("a", "TTTT")])
        with pytest.raises(ValueError, match="different record counts"):
            validate_paired_fastq_ids(str(r1), str(r2))

    def test_rejects_truncated_record(self, tmp_path: Path) -> None:
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        write_fastq(r1, [("good", "AAAA")])
        write_fastq(r2, [("good", "TTTT"), ("truncated", "GGGG")])
        with r1.open("a") as handle:
            handle.write("@truncated\n")
        with pytest.raises(ValueError, match="Truncated FASTQ"):
            validate_paired_fastq_ids(str(r1), str(r2))

    def test_cb_umi_wide_api_is_disabled(self) -> None:
        with pytest.raises(RuntimeError, match="removed in ViralScan v3"):
            filter_fastq_pairs()


def test_filter_audit_records_retained_read_ids(tmp_path: Path) -> None:
    r1 = tmp_path / "filtered_R1.fastq.gz"
    r2 = tmp_path / "filtered_R2.fastq.gz"
    write_fastq(r1, [("kept/1", "A" * 28)])
    write_fastq(r2, [("kept/2", "ACGT")])

    host_filter._write_filter_audit(tmp_path, 3, 1, str(r1), str(r2))

    with (tmp_path / "host_filter_audit.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["fragments"] for row in rows] == ["3", "1", "2"]
    with gzip.open(tmp_path / "fragment_lineage.tsv.gz", "rt") as handle:
        lineage = list(csv.DictReader(handle, delimiter="\t"))
    assert lineage == [
        {"read_id": "kept", "filter_decision": "retained", "reason": "host_unmapped"}
    ]


class TestHostFilterToolPreflight:
    def test_starsolo_mode_reports_missing_star(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(host_filter.shutil, "which", lambda _tool: None)
        with pytest.raises(RuntimeError, match=r"--host-filter starsolo requires STAR"):
            check_host_filter_tools("starsolo")

    def test_starsolo_accepts_star(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(host_filter.shutil, "which", lambda tool: f"/bin/{tool}")
        check_host_filter_tools("starsolo")

    def test_kallisto_is_rejected_even_when_installed(self) -> None:
        with pytest.raises(ValueError, match="not available in ViralScan v3"):
            check_host_filter_tools("kallisto")

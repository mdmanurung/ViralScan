"""Tests for viral read-evidence tracing and extraction (pure logic)."""

from __future__ import annotations

import csv
import gzip

import pytest

from tests._fastq import write_fastq
from viralscan.evidence import (
    _alignment_qc_from_text,
    _is_low_complexity,
    _parse_blast_output,
    _parse_competitive_blast_output,
    _parse_coverage_output,
    _per_cell_qc_from_text,
    _run,
    cb_umi_geometry,
    extract_exact_reads_by_number,
    extract_viral_reads,
    have_tools,
    interpretation_flags,
    parse_flagged_target_bus,
    resolve_viral_target,
    sample_fasta_deterministic,
    viral_assigned_keys,
    viral_equivalence_classes,
    write_competitive_fasta,
    write_igv_session,
)


class TestRunSurfacesErrors:
    def test_failing_command_raises_with_stderr(self) -> None:
        # The fail-fast contract: a non-zero exit raises RuntimeError that
        # INCLUDES the tool's stderr, never a bare swallowed CalledProcessError.
        with pytest.raises(RuntimeError) as exc:
            _run(["sh", "-c", "echo boom-message 1>&2; exit 3"])
        assert "exit 3" in str(exc.value)
        assert "boom-message" in str(exc.value)

    def test_capture_returns_stdout(self) -> None:
        assert _run(["printf", "hello"], capture=True) == b"hello"

    def test_have_tools_reports_missing(self) -> None:
        assert have_tools(["definitely-not-a-real-binary-xyz"]) == [
            "definitely-not-a-real-binary-xyz"
        ]
        assert have_tools(["sh"]) == []


class TestGeometry:
    def test_named_chemistries(self) -> None:
        assert cb_umi_geometry("10xv2") == (16, 10)
        assert cb_umi_geometry("10xv3") == (16, 12)
        # Drop-seq is case-insensitive and returns the correct (12, 8) geometry.
        # Regression guard: the old _cb_umi_lengths() silently returned (16,12)
        # for any non-10x chemistry.
        assert cb_umi_geometry("DROPSEQ") == (12, 8)
        assert cb_umi_geometry("dropseq") == (12, 8)
        assert cb_umi_geometry("DropSeq") == (12, 8)

    def test_explicit_geometry_string(self) -> None:
        # Explicit kallisto "bc:umi:seq" triplets bypass the named-chemistry table.
        assert cb_umi_geometry("0,0,16:0,16,28:1,0,0") == (16, 12)  # 10xv3 layout
        assert cb_umi_geometry("0,0,12:0,12,20:1,0,0") == (12, 8)  # Drop-seq layout

    def test_unknown_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown technology"):
            cb_umi_geometry("nanopore-bonkers")


class TestViralEcs:
    def test_selects_ecs_touching_viral_genes(self) -> None:
        ec_map = {0: [0], 1: [1], 2: [0, 1], 3: [3]}
        viral = {1}
        assert viral_equivalence_classes(ec_map, viral) == {1, 2}


class TestResolveViralTarget:
    _GENES = ["EPSTEIN_BNRF1", "EPSTEIN_EBNA1", "HUM_HERP6B_U90"]

    def test_resolves_exact_gene(self) -> None:
        label, genes = resolve_viral_target("EPSTEIN_EBNA1", self._GENES)
        assert label == "EPSTEIN_EBNA1"
        assert genes == ["EPSTEIN_EBNA1"]

    def test_resolves_registered_alias_without_substring_matching(self) -> None:
        label, genes = resolve_viral_target("EBV", self._GENES)
        assert label == "Epstein-Barr virus"
        assert genes == ["EPSTEIN_BNRF1", "EPSTEIN_EBNA1"]

    def test_resolves_detected_canonical_call(self) -> None:
        label, genes = resolve_viral_target(
            "Epstein-Barr virus",
            self._GENES,
            detected_virus_names=["Epstein-Barr virus"],
        )
        assert label == "Epstein-Barr virus"
        assert len(genes) == 2

    def test_rejects_substring_and_empty_selector(self) -> None:
        with pytest.raises(ValueError, match="No exact viral target"):
            resolve_viral_target("EBNA", self._GENES)
        with pytest.raises(ValueError, match="exact --virus selector is required"):
            resolve_viral_target("", self._GENES)


class TestViralKeys:
    def test_collects_barcode_umi_for_viral_ecs(self) -> None:
        lines = [
            "AAAA\tUUUU1\t2\t1",  # viral EC -> kept
            "AAAA\tUUUU2\t0\t1",  # host EC  -> skipped
            "CCCC\tUUUU3\t2\t3",  # viral EC -> kept
            "garbage line",
        ]
        keys = viral_assigned_keys(lines, viral_ecs={2})
        assert keys == {("AAAA", "UUUU1"), ("CCCC", "UUUU3")}


class TestExactReadLineage:
    def test_parses_flagged_bus_and_host_conservative_weight(self) -> None:
        rows = parse_flagged_target_bus(
            ["AAAA\tUUUU\t2\t1\t1\n", "CCCC\tUUU2\t1\t1\t3\n"],
            {1: [1], 2: [0, 1]},
            target_gene_indices={1},
            viral_gene_indices={1},
            method="host-conservative",
        )
        assert set(rows) == {1, 3}
        assert rows[1]["ambiguity_class"] == "host_virus_ambiguous"
        assert rows[1]["assigned_weight"] == 0.0
        assert rows[3]["assigned_weight"] == 1.0

    def test_extracts_exact_flagged_read_and_writes_lineage(self, tmp_path) -> None:
        r1, r2 = tmp_path / "R1.fastq", tmp_path / "R2.fastq"
        write_fastq(r1, [("a/1", "A" * 28), ("b/1", "C" * 28)])
        write_fastq(r2, [("a/2", "AAAA"), ("b/2", "TTTT")])
        metadata = {
            1: {
                "cb": "C" * 16,
                "ub": "C" * 12,
                "ecs": "2",
                "compatible_genes": "0,1",
                "ambiguity_class": "host_virus_ambiguous",
                "assigned_weight": 0.0,
                "method": "host-conservative",
                "evidence_tier": "candidate_host_virus_ambiguous",
                "exclusion_reason": "",
            }
        }
        fasta = tmp_path / "target.fasta"
        lineage = tmp_path / "read_lineage.tsv.gz"
        stats = extract_exact_reads_by_number(str(r1), str(r2), metadata, str(fasta), str(lineage))
        assert stats.total_reads == 2 and stats.viral_reads == 1
        assert "TTTT" in fasta.read_text() and "AAAA" not in fasta.read_text()
        with gzip.open(lineage, "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert rows[0]["read_number"] == "1"
        assert rows[0]["read_id"] == "b"

    def test_rejects_out_of_range_bus_read_number(self, tmp_path) -> None:
        r1, r2 = tmp_path / "R1.fastq", tmp_path / "R2.fastq"
        write_fastq(r1, [("a", "AAAA")])
        write_fastq(r2, [("a", "TTTT")])
        with pytest.raises(ValueError, match="exact lineage failed"):
            extract_exact_reads_by_number(
                str(r1),
                str(r2),
                {2: {}},
                str(tmp_path / "target.fasta"),
                str(tmp_path / "lineage.tsv.gz"),
            )


class TestParseCoverageOutput:
    """_parse_coverage_output parses samtools coverage TSV without the binary."""

    _HEADER = "#rname\tstart\tend\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
    _ROW_HIT = "virus_A\t0\t29903\t412\t27100\t90.6\t4.1\t37\t60"
    _ROW_ZERO = "virus_B\t0\t2000\t0\t0\t0.0\t0.0\t0\t0"

    def test_parses_covered_references(self) -> None:
        text = "\n".join([self._HEADER, self._ROW_HIT, self._ROW_ZERO])
        rows = _parse_coverage_output(text)
        assert len(rows) == 1
        assert rows[0]["rname"] == "virus_A"
        assert rows[0]["numreads"] == "412"

    def test_excludes_zero_read_references(self) -> None:
        text = "\n".join([self._HEADER, self._ROW_ZERO])
        assert _parse_coverage_output(text) == []

    def test_empty_output(self) -> None:
        assert _parse_coverage_output("") == []

    def test_header_hash_stripped(self) -> None:
        text = "\n".join([self._HEADER, self._ROW_HIT])
        rows = _parse_coverage_output(text)
        # Column key must be "rname" not "#rname"
        assert "rname" in rows[0]
        assert "#rname" not in rows[0]


def test_alignment_qc_reports_depth_hotspots_cells_and_host_competition() -> None:
    header = "@SQ\tSN:VIRUS|v\tLN:10\n@SQ\tSN:HOST|h\tLN:10\n"
    sam = (
        "CB1_U1_0\t0\tVIRUS|v\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
        "CB2_U2_1\t16\tVIRUS|v\t2\t40\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:1\n"
        "CB3_U3_2\t0\tHOST|h\t1\t50\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
    )
    depth = "VIRUS|v\t1\t1\nVIRUS|v\t2\t3\nVIRUS|v\t3\t10\nHOST|h\t1\t1\n"
    rows = _alignment_qc_from_text(header, sam, depth)
    virus = next(row for row in rows if row["reference"] == "VIRUS|v")
    assert virus["reads"] == 2 and virus["cells"] == 2 and virus["molecules"] == 2
    assert virus["breadth_1x"] == pytest.approx(0.3)
    assert virus["breadth_3x"] == pytest.approx(0.2)
    assert virus["breadth_10x"] == pytest.approx(0.1)
    assert virus["median_depth"] == 0.0
    assert virus["host_competitive_fraction"] == pytest.approx(1 / 3)


def test_per_cell_qc_separates_host_and_virus_support() -> None:
    sam = (
        "CB1_U1_0\t0\tVIRUS|v\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
        "CB1_U2_1\t16\tHOST|h\t2\t40\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:1\n"
    )
    rows = _per_cell_qc_from_text(sam)
    assert {(row["cell_barcode"], row["reference_class"]) for row in rows} == {
        ("CB1", "host"),
        ("CB1", "virus"),
    }
    assert all(row["molecules"] == 1 for row in rows)


def test_low_complexity_and_interpretation_flags_are_explicit() -> None:
    assert _is_low_complexity("AAAAAAAAAAC") is True
    assert _is_low_complexity("ACGTACGTACGT") is False
    rows = interpretation_flags(
        [
            {
                "reference_class": "virus",
                "host_competitive_fraction": 0.4,
                "max_50bp_window_fraction": 0.8,
                "breadth_1x": 0.05,
            }
        ],
        [{"viral_minus_host_bitscore": "-1", "low_complexity": "true"}],
        [{"ambiguity_class": "host_virus_ambiguous"}],
    )
    status = {row["flag"]: row["status"] for row in rows}
    assert status["host_homology"] == "flagged"
    assert status["low_complexity"] == "flagged"
    assert status["eve_or_integration_like_hotspot"] == "flagged"
    assert status["contamination"] == "not_assessed"


class TestParseBlastOutput:
    """_parse_blast_output parses blastn -outfmt 6 text without the binary."""

    _BLAST_TEXT = (
        "read1\tNC_045512.2\t98.5\t150\t1e-80\n"
        "read1\tNC_002549.1\t72.0\t130\t1e-20\n"  # second hit for read1 — ignored
        "read2\tNC_045512.2\t95.0\t148\t1e-75\n"
        "read3\tNC_002549.1\t60.0\t100\t1e-10\n"
    )

    def test_best_hit_per_read(self) -> None:
        rows = _parse_blast_output(self._BLAST_TEXT)
        # read1 appears twice; only the first (best) hit must be returned
        read_ids = [r["read"] for r in rows]
        assert read_ids.count("read1") == 1

    def test_all_reads_present(self) -> None:
        rows = _parse_blast_output(self._BLAST_TEXT)
        assert {r["read"] for r in rows} == {"read1", "read2", "read3"}

    def test_fields_populated(self) -> None:
        rows = _parse_blast_output(self._BLAST_TEXT)
        r1 = next(r for r in rows if r["read"] == "read1")
        assert r1["subject"] == "NC_045512.2"
        assert r1["pident"] == "98.5"
        assert r1["length"] == "150"

    def test_empty_output(self) -> None:
        assert _parse_blast_output("") == []

    def test_short_lines_skipped(self) -> None:
        # Lines with fewer than 4 columns are silently ignored
        text = "only\ttwo\n" + "read1\tNC_045512.2\t98.5\t150\t1e-80\n"
        rows = _parse_blast_output(text)
        assert len(rows) == 1
        assert rows[0]["read"] == "read1"

    def test_competitive_parser_reports_host_viral_score_difference(self) -> None:
        rows = _parse_competitive_blast_output(
            "r1\tVIRUS|v\t99\t50\t100\t1e-20\t80\nr1\tHOST|h\t97\t48\t96\t1e-18\t70\n"
        )
        assert rows[0]["top_viral_hit"] == "VIRUS|v"
        assert rows[0]["top_host_hit"] == "HOST|h"
        assert rows[0]["viral_minus_host_bitscore"] == "10"


def test_competitive_fasta_and_igv_session(tmp_path) -> None:
    host, virus = tmp_path / "host.fa", tmp_path / "virus.fa"
    host.write_text(">chr1 host\nAAAA\n")
    virus.write_text(">NC_1 virus\nCCCC\n")
    combined = tmp_path / "competitive.fa"
    write_competitive_fasta(str(host), str(virus), str(combined))
    assert ">HOST|chr1 host" in combined.read_text()
    assert ">VIRUS|NC_1 virus" in combined.read_text()
    bam = tmp_path / "raw.bam"
    bam.write_bytes(b"bam")
    session = tmp_path / "session.xml"
    write_igv_session(str(combined), [str(bam)], str(session))
    assert str(bam.resolve()) in session.read_text()


class TestDeterministicFastaSampling:
    def test_sampling_is_seeded_and_order_invariant(self, tmp_path) -> None:
        records = [f">read{i}\n{'ACGT' * (i + 1)}\n" for i in range(10)]
        source_a = tmp_path / "a.fasta"
        source_b = tmp_path / "b.fasta"
        source_a.write_text("".join(records))
        source_b.write_text("".join(reversed(records)))
        out_a, out_b = tmp_path / "sample_a.fasta", tmp_path / "sample_b.fasta"
        assert sample_fasta_deterministic(str(source_a), out_a, 4, 17) == (4, 10)
        assert sample_fasta_deterministic(str(source_b), out_b, 4, 17) == (4, 10)
        assert set(out_a.read_text().splitlines()[::2]) == set(out_b.read_text().splitlines()[::2])

    def test_sampling_rejects_zero_cap(self, tmp_path) -> None:
        source = tmp_path / "source.fasta"
        source.write_text(">r\nACGT\n")
        with pytest.raises(ValueError, match="max_records"):
            sample_fasta_deterministic(str(source), tmp_path / "out.fasta", 0, 1)


class TestExtractReads:
    def test_extracts_only_viral_assigned_cdna_reads(self, tmp_path) -> None:
        # 10xv2 geometry: CB=16, UMI=10 -> first 26 bases of R1
        cb = "A" * 16
        umi_hit = "C" * 10
        umi_miss = "G" * 10
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        write_fastq(
            r1,
            [("r1", cb + umi_hit), ("r2", cb + umi_miss), ("r3", cb + umi_hit)],
        )
        write_fastq(
            r2,
            [("r1", "ACGTACGTAC"), ("r2", "TTTTTTTTTT"), ("r3", "GGGGCCCCGG")],
        )
        out = tmp_path / "evidence.fasta"
        stats = extract_viral_reads(
            str(r1), str(r2), keys={(cb, umi_hit)}, cb_len=16, umi_len=10, out_fasta=str(out)
        )
        assert stats.total_reads == 3
        assert stats.viral_reads == 2  # r1 and r3 (umi_hit), not r2 (umi_miss)
        body = out.read_text()
        assert "ACGTACGTAC" in body and "GGGGCCCCGG" in body
        assert "TTTTTTTTTT" not in body

    def test_handles_gzipped_input(self, tmp_path) -> None:
        cb, umi = "A" * 16, "C" * 12  # 10xv3
        r1 = tmp_path / "R1.fastq.gz"
        r2 = tmp_path / "R2.fastq.gz"
        write_fastq(r1, [("r1", cb + umi)])
        write_fastq(r2, [("r1", "ACGTACGTACGT")])
        out = tmp_path / "ev.fasta"
        stats = extract_viral_reads(
            str(r1), str(r2), keys={(cb, umi)}, cb_len=16, umi_len=12, out_fasta=str(out)
        )
        assert stats.viral_reads == 1
        assert "ACGTACGTACGT" in out.read_text()

"""Tests for ``strand_qc.tsv`` (PLAN EVID-03): gene-relative strand, calibrated on host reads."""

import csv
import gzip
import hashlib
import shutil
import subprocess
from pathlib import Path

import pytest

from viralscan.evidence import alignment_qc_table
from viralscan.strand_qc import (
    ALL_VIRAL,
    FIELDS,
    MIN_HOST_BLOCKS,
    aligned_blocks,
    classify_block,
    load_features,
    strand_qc_rows,
    write_strand_qc,
)

VIRUS = "VIRUS|V1"
HOST_GTF = (
    'chr1\tx\texon\t1000\t2000\t.\t+\t.\tgene_id "A";\n'  # plus gene
    'chr1\tx\texon\t5000\t6000\t.\t-\t.\tgene_id "B";\n'  # minus gene
    'chr1\tx\texon\t8000\t8500\t.\t+\t.\tgene_id "E";\n'  # overlaps F on the other strand
    'chr1\tx\texon\t8200\t8700\t.\t-\t.\tgene_id "F";\n'
    'chr1\tx\texon\t9000\t9500\t.\t+\t.\tgene_id "G";\n'  # two plus genes overlap
    'chr1\tx\texon\t9200\t9800\t.\t+\t.\tgene_id "H";\n'
)
VIRAL_GTF = (
    'V1\tx\tgene\t100\t1000\t.\t+\t.\tgene_id "g1";\n'
    'V1\tx\tgene\t2000\t3000\t.\t-\t.\tgene_id "g2";\n'
    'V1\tx\tgene\t3500\t4000\t.\t+\t.\tgene_id "g3";\n'
    'V1\tx\tgene\t3800\t4200\t.\t-\t.\tgene_id "g4";\n'
)
TECH = "10xv3"


def sam(n: int, ref: str, pos: int, flag: int, cigar: str = "50M", score: int = 90) -> str:
    return f"CB{n}_UMI{n}_{n}|r{n}\t{flag}\t{ref}\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*\tAS:i:{score}"


class Lib:
    """Synthetic library: molecules with one read each; ``add`` returns the next molecule number."""

    def __init__(self) -> None:
        self.sam: list[str] = []
        self.verdicts: list[dict[str, str]] = []
        self.lineage: list[dict[str, str]] = []

    def add(
        self,
        ref: str,
        pos: int,
        flag: int,
        *,
        cigar: str = "50M",
        verdict: str,
        counted: bool = True,
    ) -> None:
        n = len(self.lineage)
        if ref:
            self.sam.append(sam(n, ref, pos, flag, cigar))
        else:
            self.sam.append(f"CB{n}_UMI{n}_{n}|r{n}\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*")
        self.verdicts.append(
            {
                "cb": f"CB{n}",
                "ub": f"UMI{n}",
                "reads": "1",
                "counted": str(counted),
                "virus_as": "80" if verdict in ("virus_best", "tie") else "",
                "host_as": "80" if verdict in ("host_best", "tie") else "",
                "verdict": verdict,
            }
        )
        self.lineage.append({"read_number": str(n), "cb": f"CB{n}", "ub": f"UMI{n}"})

    def host(self, pos: int, flag: int, k: int = 1) -> None:
        for _ in range(k):
            self.add("HOST|chr1", pos, flag, verdict="host_best")

    def virus(self, pos: int, flag: int, k: int = 1, cigar: str = "50M") -> None:
        for _ in range(k):
            self.add(VIRUS, pos, flag, cigar=cigar, verdict="virus_best")


def run(lib: Lib, tmp_path: Path, *, tech: str | None = TECH, host_gtf: bool = True):
    (tmp_path / "host.gtf").write_text(HOST_GTF)
    (tmp_path / "viral.gtf").write_text(VIRAL_GTF)
    rows = strand_qc_rows(
        lib.sam,
        lib.verdicts,
        lib.lineage,
        technology=tech,
        host_features=load_features([tmp_path / "host.gtf"], ("exon",)) if host_gtf else None,
        viral_features=load_features([tmp_path / "viral.gtf"], ("gene", "exon", "CDS")),
    )
    return {r["reference"]: r for r in rows}


def viral_mix(lib: Lib) -> None:
    lib.virus(200, 0, 7)  # inside g1 (+), forward read
    lib.virus(200, 16, 1)  # inside g1, reverse read
    lib.virus(3850, 0, 2)  # g3 (+) and g4 (-) overlap
    lib.virus(4600, 0, 3)  # no gene


class TestBlocks:
    def test_unspliced_and_spliced(self) -> None:
        assert aligned_blocks(100, "50M") == [(100, 149)]
        assert aligned_blocks(100, "10S20M5D20M") == [(100, 144)]
        assert aligned_blocks(100, "30M200N20M") == [(100, 129), (330, 349)]

    def test_classes(self) -> None:
        assert classify_block(False, [("+", "a")]) == "sense"
        assert classify_block(True, [("+", "a")]) == "antisense"
        assert classify_block(True, [("-", "a")]) == "sense"
        assert classify_block(False, [("+", "a"), ("-", "b")]) == "ambiguous_overlap"
        assert classify_block(False, [("+", "a"), ("+", "b")]) == "sense"
        assert (
            classify_block(False, [("+", "a"), ("+", "b")], single_gene=True) == "ambiguous_overlap"
        )
        assert classify_block(False, [(".", "a")]) == "ambiguous_overlap"
        assert classify_block(False, []) == "outside_annotation"


class TestCalibration:
    def test_three_prime_library_expects_sense(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS)  # + gene, forward read: sense
        lib.host(5100, 16, 20)  # - gene, reverse read: sense
        lib.host(8300, 0, 30)  # overlaps E+ and F-: ambiguous, never used
        lib.host(9300, 0, 30)  # two genes: ambiguous under the single-gene rule
        lib.host(7000, 0, 30)  # no exon: outside
        viral_mix(lib)
        r = run(lib, tmp_path)[VIRUS]
        assert r["status"] == "assessed" and r["expected_orientation"] == "sense"
        assert (r["host_blocks_used"], r["host_blocks_sense"]) == (MIN_HOST_BLOCKS + 20,) * 2
        assert (r["host_blocks_ambiguous"], r["host_blocks_outside"]) == (60, 30)
        assert r["blocks_sense"] == 7 and r["blocks_antisense"] == 1
        assert r["blocks_ambiguous_overlap"] == 2 and r["blocks_outside_annotation"] == 3
        assert r["blocks_expected_orientation"] == 7 and r["blocks_unexpected_orientation"] == 1
        assert r["fraction_expected_orientation"] == pytest.approx(7 / 8)
        assert r["molecules_sense"] == 7 and r["reads_antisense"] == 1
        assert run(lib, tmp_path)[ALL_VIRAL]["blocks_sense"] == 7

    def test_five_prime_library_flips_the_expectation(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 16, MIN_HOST_BLOCKS)  # + gene, reverse read: antisense
        lib.host(5100, 0, 20)  # - gene, forward read: antisense
        lib.virus(200, 16, 9)  # antisense to g1, the library's normal orientation
        lib.virus(200, 0, 1)
        r = run(lib, tmp_path)[VIRUS]
        assert r["expected_orientation"] == "antisense"
        assert r["blocks_antisense"] == 9 and r["blocks_expected_orientation"] == 9
        assert r["fraction_expected_orientation"] == pytest.approx(0.9)

    def test_spliced_read_gives_two_blocks_one_read(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS)
        lib.virus(300, 0, 1, cigar="50M200N50M")
        r = run(lib, tmp_path)[VIRUS]
        assert r["blocks_sense"] == 2 and r["reads_sense"] == 1 and r["molecules_sense"] == 1


class TestUnit:
    def test_only_counted_virus_best_molecules_are_classified(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS)
        lib.virus(200, 0, 2)
        lib.add(VIRUS, 200, 16, verdict="virus_best", counted=False)
        lib.add(VIRUS, 200, 16, verdict="tie")
        lib.add("", 0, 0, verdict="unaligned")
        r = run(lib, tmp_path)[VIRUS]
        assert r["blocks_antisense"] == 0 and r["blocks_sense"] == 2
        assert r["n_molecules_counted_virus_best"] == 2
        assert r["n_molecules_counted_host_best"] == MIN_HOST_BLOCKS
        assert (r["n_molecules_counted_tie"], r["n_molecules_counted_unaligned"]) == (1, 1)
        assert r["n_molecules_not_counted"] == 1 and r["n_reads_not_counted"] == 1
        assert r["verdict_ties_resolvable"] is True

    def test_all_member_reads_count_not_the_representative(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS)
        lib.virus(200, 0, 1)
        # a second read of molecule `m`, on the opposite strand
        m = len(lib.lineage) - 1
        lib.sam.append(sam(1000, VIRUS, 200, 16).replace("CB1000_UMI1000", f"CB{m}_UMI{m}"))
        lib.lineage.append({"read_number": "1000", "cb": f"CB{m}", "ub": f"UMI{m}"})
        r = run(lib, tmp_path)[VIRUS]
        assert r["reads_sense"] == 1 and r["reads_antisense"] == 1
        assert r["molecules_sense"] == 1 and r["molecules_antisense"] == 1
        assert r["n_reads_counted_virus_best"] == 2 and r["n_molecules_counted_virus_best"] == 1

    def test_no_counted_virus_best_molecule_is_not_evaluable(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS)
        r = run(lib, tmp_path)[ALL_VIRAL]
        assert r["status"] == "not_evaluable"


class TestNotAssessed:
    @pytest.mark.parametrize("tech", ["custom:0,0,16", None, ""])
    def test_unknown_chemistry(self, tmp_path: Path, tech: str | None) -> None:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS)
        viral_mix(lib)
        r = run(lib, tmp_path, tech=tech)[ALL_VIRAL]
        assert (r["status"], r["reason"]) == ("not_assessed", "chemistry_unknown")
        assert r["blocks_sense"] == ""  # nothing computed

    def test_no_host_gtf_keeps_viral_counts(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS)
        viral_mix(lib)
        r = run(lib, tmp_path, host_gtf=False)[VIRUS]
        assert (r["status"], r["reason"]) == ("not_assessed", "host_gtf_missing")
        assert r["expected_orientation"] == "" and r["fraction_expected_orientation"] == ""
        assert r["blocks_sense"] == 7

    def test_too_few_host_blocks(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS - 1)
        viral_mix(lib)
        r = run(lib, tmp_path)[VIRUS]
        assert (r["status"], r["reason"]) == ("not_assessed", "too_few_host_blocks")
        assert r["host_blocks_used"] == MIN_HOST_BLOCKS - 1

    def test_mixed_host_orientation(self, tmp_path: Path) -> None:
        lib = Lib()
        lib.host(1100, 0, 150)
        lib.host(1100, 16, 150)
        viral_mix(lib)
        r = run(lib, tmp_path)[VIRUS]
        assert (r["status"], r["reason"]) == ("not_assessed", "host_orientation_mixed")
        assert r["host_sense_fraction"] == 0.5


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
class TestWrite:
    def _evidence_dir(self, tmp_path: Path) -> tuple[Path, Path]:
        lib = Lib()
        lib.host(1100, 0, MIN_HOST_BLOCKS)
        viral_mix(lib)
        out = tmp_path / "ev"
        out.mkdir()
        header = [
            "@HD\tVN:1.6\tSO:unsorted",
            "@SQ\tSN:HOST|chr1\tLN:100000",
            f"@SQ\tSN:{VIRUS}\tLN:5000",
        ]
        sam_path = out / "competitive_reads.raw.sam"
        sam_path.write_text("\n".join(header + lib.sam) + "\n")
        bam = out / "competitive_reads.raw.bam"
        subprocess.run(["samtools", "sort", "-o", str(bam), str(sam_path)], check=True)
        sam_path.unlink()
        with gzip.open(out / "molecule_verdicts.tsv.gz", "wt", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(lib.verdicts[0]), delimiter="\t")
            w.writeheader()
            w.writerows(lib.verdicts)
        with gzip.open(out / "read_lineage.tsv.gz", "wt", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(lib.lineage[0]), delimiter="\t")
            w.writeheader()
            w.writerows(lib.lineage)
        (tmp_path / "host.gtf").write_text(HOST_GTF)
        (tmp_path / "viral.gtf").write_text(VIRAL_GTF)
        return out, bam

    @staticmethod
    def _hashes(directory: Path) -> dict[str, str]:
        return {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in directory.iterdir()
            if p.is_file()
        }

    def test_adds_only_strand_qc_and_leaves_alignment_qc_alone(self, tmp_path: Path) -> None:
        out, bam = self._evidence_dir(tmp_path)
        before_files, before_qc = self._hashes(out), alignment_qc_table(str(bam))
        write_strand_qc(
            bam,
            out,
            technology=TECH,
            host_gtf=tmp_path / "host.gtf",
            viral_gtfs=[tmp_path / "viral.gtf"],
        )
        assert alignment_qc_table(str(bam)) == before_qc
        after = self._hashes(out)
        assert set(after) - set(before_files) == {"strand_qc.tsv"}
        assert {k: after[k] for k in before_files} == before_files
        with (out / "strand_qc.tsv").open() as h:
            reader = csv.DictReader(h, delimiter="\t")
            assert tuple(reader.fieldnames or ()) == FIELDS
            rows = {r["reference"]: r for r in reader}
        assert rows[VIRUS]["status"] == "assessed" and rows[VIRUS]["blocks_sense"] == "7"
        assert rows[ALL_VIRAL]["blocks_outside_annotation"] == "3"

    def test_competitor_mode_computes_nothing(self, tmp_path: Path) -> None:
        out, bam = self._evidence_dir(tmp_path)
        write_strand_qc(
            bam, out, technology=TECH, host_gtf=None, viral_gtfs=[], competitor_mode=True
        )
        with (out / "strand_qc.tsv").open() as h:
            (row,) = list(csv.DictReader(h, delimiter="\t"))
        assert row["status"] == "not_assessed"
        assert row["reason"] == "verdicts_not_produced_in_competitor_mode"

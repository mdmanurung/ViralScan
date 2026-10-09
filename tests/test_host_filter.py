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

    in1 = tmp_path / "in_R1.fastq.gz"
    in2 = tmp_path / "in_R2.fastq.gz"
    write_fastq(in1, [("gone-a/1", "A" * 28), ("kept/1", "A" * 28), ("gone-b/1", "A" * 28)])
    write_fastq(in2, [("gone-a/2", "ACGT"), ("kept/2", "ACGT"), ("gone-b/2", "ACGT")])

    host_filter._write_filter_audit(tmp_path, 3, 1, str(r1), str(r2), str(in1), str(in2))

    with (tmp_path / "host_filter_audit.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_category = {row["category"]: row["fragments"] for row in rows}
    assert by_category["input"] == "3"
    assert by_category["retained_host_unmapped"] == "1"
    assert by_category["removed_host_aligned_or_ambiguous"] == "2"
    assert by_category["pct_retained"] == "33.33"
    # Every pinned STAR filter parameter is recorded, so a run that retains a
    # different fraction can be attributed to a parameter change rather than to
    # a STAR version difference.
    pinned = {
        k[len("star_param:") :]: v for k, v in by_category.items() if k.startswith("star_param:")
    }
    assert pinned["outFilterMismatchNmax"] == "4"
    assert pinned["outFilterMatchNminOverLread"] == "0.9"
    assert pinned["outFilterMultimapNmax"] == "20"
    assert set(pinned) == {
        arg.lstrip("-") for i, arg in enumerate(host_filter.STAR_FILTER_ARGS) if i % 2 == 0
    }
    assert by_category["star_param_set"] == "pinned"
    with gzip.open(tmp_path / "fragment_lineage.tsv.gz", "rt") as handle:
        lineage = list(csv.DictReader(handle, delimiter="\t"))
    # SW-07: removed fragments are listed too, in input order.
    assert [(r["read_id"], r["filter_decision"], r["reason"]) for r in lineage] == [
        ("gone-a", "removed", "host_mapped"),
        ("kept", "retained", "host_unmapped"),
        ("gone-b", "removed", "host_mapped"),
    ]


def _lineage(path: Path, rows) -> str:
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["read_id", "filter_decision", "reason"])
        w.writerows(rows)
    return str(path)


def _truth(path: Path, rows, header="read_id\tlabel\tmolecule_id") -> str:
    path.write_text(header + "\n" + "\n".join(rows) + "\n")
    return str(path)


class TestLostTruthCounts:
    def test_d15_d16_counts(self, tmp_path: Path) -> None:
        lineage = _lineage(
            tmp_path / "l.tsv.gz",
            [
                ("v1", "removed", "host_mapped"),
                ("v2", "retained", "host_unmapped"),
                ("v3", "removed", "host_mapped"),
                ("v4", "removed", "host_mapped"),
                ("h1", "removed", "host_mapped"),
            ],
        )
        truth = _truth(
            tmp_path / "t.tsv",
            ["v1\tviral\tm1", "v2\tviral\tm1", "v3\tviral\tm2", "v4\tviral\tm3", "h1\thost\tm9"],
        )
        got = host_filter.lost_truth_counts(truth, lineage)
        assert got["d15_truth_fragments"] == 4 and got["d15_removed_fragments"] == 3
        assert got["d15_loss_fraction"] == 0.75
        # m1 survives through v2; m2 and m3 have no surviving fragment.
        assert got["d16_truth_molecules"] == 3 and got["d16_lost_molecules"] == 2
        assert got["d16_loss_fraction"] == pytest.approx(2 / 3)
        assert got["not_in_lineage"] == 0

    def test_without_molecule_column_d16_is_none(self, tmp_path: Path) -> None:
        lineage = _lineage(tmp_path / "l.tsv.gz", [("v1", "retained", "host_unmapped")])
        truth = _truth(tmp_path / "t.tsv", ["v1\tviral"], header="read_id\tlabel")
        got = host_filter.lost_truth_counts(truth, lineage)
        assert got["d15_loss_fraction"] == 0.0 and got["d16_lost_molecules"] is None

    def test_truth_read_missing_from_lineage_is_not_estimable(self, tmp_path: Path) -> None:
        lineage = _lineage(tmp_path / "l.tsv.gz", [("v1", "retained", "host_unmapped")])
        truth = _truth(tmp_path / "t.tsv", ["v1\tviral\tm1", "v2\tviral\tm2"])
        got = host_filter.lost_truth_counts(truth, lineage)
        assert got["not_in_lineage"] == 1
        assert got["d15_loss_fraction"] is None and got["d16_loss_fraction"] is None


def test_star_filter_args_are_pinned() -> None:
    """The filter must not inherit the installed STAR's defaults.

    The pinned set is stricter than STAR's defaults on mismatches and coverage
    (4 vs 10, 0.9 vs 0.66) and looser on multimapping (20 vs 10 loci).
    """
    args = dict(
        zip(
            host_filter.STAR_FILTER_ARGS[::2],
            host_filter.STAR_FILTER_ARGS[1::2],
        )
    )
    assert args["--outFilterMismatchNmax"] == "4"
    assert args["--outFilterMatchNminOverLread"] == "0.9"
    assert args["--outFilterMultimapNmax"] == "20"
    assert args["--outFilterType"] == "Normal"
    # Even length, so every entry is a (flag, value) pair.
    assert len(host_filter.STAR_FILTER_ARGS) % 2 == 0
    # STAR rejects `--outSAMflag None` as an unknown parameter value (SW-15).
    assert "--outSAMflag" not in args


def test_star_default_set_pins_the_same_flags_at_star_defaults() -> None:
    """``star-default`` is the protocol's second star_alignment grid point.

    Values are STAR 2.7.11b ``--help`` defaults, passed explicitly so another
    STAR install cannot shift them.
    """
    pinned = host_filter.STAR_FILTER_PARAM_SETS["pinned"]
    default = host_filter.STAR_FILTER_PARAM_SETS["star-default"]
    assert pinned == host_filter.STAR_FILTER_ARGS
    assert pinned[::2] == default[::2]
    args = dict(zip(default[::2], default[1::2]))
    assert args["--outFilterMismatchNmax"] == "10"
    assert args["--outFilterMatchNminOverLread"] == "0.66"
    assert args["--outFilterMultimapNmax"] == "10"
    assert args["--outFilterMismatchNoverReadLmax"] == "1.0"


def test_starsolo_filter_passes_the_selected_param_set(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_fastq(tmp_path / "in_R1.fastq.gz", [("kept/1", "A" * 28)])
    write_fastq(tmp_path / "in_R2.fastq.gz", [("kept/2", "ACGT")])
    calls = []

    def fake_star(cmd, check):
        calls.append(cmd)
        star_tmp = tmp_path / "star_tmp"
        (star_tmp / "Unmapped.out.mate1").write_text("@kept/2\nACGT\n+\nIIII\n")
        (star_tmp / "Unmapped.out.mate2").write_text(
            "@kept/1\n" + "A" * 28 + "\n+\n" + "I" * 28 + "\n"
        )

    monkeypatch.setattr(host_filter.subprocess, "run", fake_star)
    host_filter._starsolo_filter(
        str(tmp_path / "in_R1.fastq.gz"),
        str(tmp_path / "in_R2.fastq.gz"),
        "/host",
        "10xv3",
        None,
        tmp_path,
        str(tmp_path / "R1.fastq.gz"),
        str(tmp_path / "R2.fastq.gz"),
        1,
        "star-default",
    )
    cmd = calls[0]
    assert cmd[cmd.index("--outFilterMismatchNmax") + 1] == "10"
    with (tmp_path / "host_filter_audit.tsv").open() as handle:
        rows = {r["category"]: r["fragments"] for r in csv.DictReader(handle, delimiter="\t")}
    assert rows["star_param_set"] == "star-default"
    assert rows["star_param:outFilterMultimapNmax"] == "10"


class TestStarsoloBarcodeArgs:
    """SW-15: one pure builder for the STARsolo barcode/UMI arguments."""

    @staticmethod
    def _args(technology: str, whitelist=None) -> dict[str, str]:
        flat = host_filter.starsolo_barcode_args(technology, whitelist)
        return dict(zip(flat[::2], flat[1::2]))

    @pytest.mark.parametrize(
        ("technology", "cb_len", "umi_len"),
        [("10xv2", 16, 10), ("10xv3", 16, 12), ("dropseq", 12, 8)],
    )
    def test_geometry_follows_technology(self, technology, cb_len, umi_len) -> None:
        args = self._args(technology)
        assert args["--soloCBlen"] == str(cb_len)
        assert args["--soloUMIstart"] == str(cb_len + 1)
        assert args["--soloUMIlen"] == str(umi_len)

    def test_barcode_read_length_check_is_disabled(self) -> None:
        """10x 5' R1 is 150 bp; STAR aborted with 'barcode length 150 != 28'."""
        assert self._args("10xv2")["--soloBarcodeReadLength"] == "0"

    def test_whitelist_is_passed_or_star_none(self) -> None:
        assert self._args("10xv3")["--soloCBwhitelist"] == "None"
        assert self._args("10xv3", "/x/wl.txt")["--soloCBwhitelist"] == "/x/wl.txt"

    def test_gzipped_whitelist_is_decompressed_for_star(self, tmp_path: Path) -> None:
        gz = tmp_path / "wl.txt.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write("AAACCCAAGAAACACT\nAAACCCAAGAAACCAT\n")

        plain = host_filter._plain_whitelist(str(gz), tmp_path)

        assert plain is not None and not plain.endswith(".gz")
        assert Path(plain).read_text() == "AAACCCAAGAAACACT\nAAACCCAAGAAACCAT\n"
        assert host_filter._plain_whitelist("/x/wl.txt", tmp_path) == "/x/wl.txt"
        assert host_filter._plain_whitelist(None, tmp_path) is None


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

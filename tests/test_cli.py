"""Tests for the CLI (menu.py).

Covers argument parsing, boolean flag regressions, and validation helpers.
No network access; no subprocesses that touch the filesystem beyond tmp dirs.
"""

from __future__ import annotations

import argparse
import subprocess
from unittest.mock import patch

import pytest

from viralscan.defaults import DEFAULTS

# ── helpers ──────────────────────────────────────────────────────────────────


def _parse(argv: list[str]) -> argparse.Namespace:
    """Run create_help() with a mocked sys.argv."""
    with patch("sys.argv", ["viralscan"] + argv):
        from viralscan.menu import create_help

        return create_help()


# ── --help smoke test ─────────────────────────────────────────────────────────


class TestHelpFlag:
    def test_help_exits_zero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            _parse(["--help"])
        assert exc.value.code == 0

    def test_data_fetch_help_exits_zero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            _parse(["data", "fetch", "--help"])
        assert exc.value.code == 0

    def test_build_ref_help_exits_zero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            _parse(["build-ref", "--help"])
        assert exc.value.code == 0


# ── boolean flag regression (§1.1) ───────────────────────────────────────────


class TestBooleanFlags:
    """Regression tests for §1.1: --visual / --multimapping must not accept
    bare strings 'True'/'False' as they used to when type=bool was used."""

    def test_visual_defaults_true(self) -> None:
        args = _parse([])
        assert args.visual is True

    def test_no_visual_sets_false(self) -> None:
        args = _parse(["--no-visual"])
        assert args.visual is False

    def test_visual_flag_sets_true(self) -> None:
        args = _parse(["--visual"])
        assert args.visual is True

    def test_multimapping_defaults_true(self) -> None:
        args = _parse([])
        assert args.multimapping is True

    def test_no_multimapping_sets_false(self) -> None:
        args = _parse(["--no-multimapping"])
        assert args.multimapping is False

    def test_multimapping_flag_sets_true(self) -> None:
        args = _parse(["--multimapping"])
        assert args.multimapping is True

    def test_reference_defaults_false(self) -> None:
        args = _parse([])
        assert args.reference is False

    def test_reference_flag_sets_true(self) -> None:
        args = _parse(["--reference"])
        assert args.reference is True

    def test_umap_defaults_false(self) -> None:
        args = _parse([])
        assert args.umap is False

    def test_umap_flag_sets_true(self) -> None:
        args = _parse(["--umap"])
        assert args.umap is True


# ── default values ────────────────────────────────────────────────────────────


class TestDefaults:
    def test_technology_default(self) -> None:
        assert _parse([]).technology == "10xv3"

    def test_cores_default(self) -> None:
        assert _parse([]).cores == 6

    def test_detection_threshold_default(self) -> None:
        assert _parse([]).detection_threshold == DEFAULTS["detection_threshold"]

    def test_min_counts_default(self) -> None:
        assert _parse([]).min_counts == DEFAULTS["min_counts"]

    def test_min_genes_default(self) -> None:
        assert _parse([]).min_genes == DEFAULTS["min_genes"]

    def test_se_threshold_default(self) -> None:
        assert _parse([]).se_threshold == DEFAULTS["se_threshold"]

    def test_hvg_min_mean_default(self) -> None:
        assert _parse([]).hvg_min_mean == DEFAULTS["hvg_min_mean"]

    def test_hvg_max_mean_default(self) -> None:
        assert _parse([]).hvg_max_mean == DEFAULTS["hvg_max_mean"]

    def test_hvg_min_disp_default(self) -> None:
        assert _parse([]).hvg_min_disp == DEFAULTS["hvg_min_disp"]

    def test_umap_n_neighbors_default(self) -> None:
        assert _parse([]).umap_n_neighbors == DEFAULTS["umap_n_neighbors"]

    def test_multimap_method_default(self) -> None:
        assert _parse([]).multimap_method == DEFAULTS["multimap_method"]

    def test_multimap_pseudocount_default(self) -> None:
        assert _parse([]).multimap_pseudocount == DEFAULTS["multimap_pseudocount"]

    def test_multimap_primary_call_default(self) -> None:
        assert _parse([]).multimap_primary_call == DEFAULTS["multimap_primary_call"]

    def test_output_defaults_none(self) -> None:
        assert _parse([]).output is None

    def test_whitelist_defaults_none(self) -> None:
        assert _parse([]).whitelist is None

    def test_ncbi_accession_defaults_none(self) -> None:
        assert _parse([]).ncbi_accession is None

    def test_data_cache_dir_defaults_none(self) -> None:
        assert _parse([]).data_cache_dir is None


# ── explicit flag parsing ──────────────────────────────────────────────────────


class TestFlagParsing:
    def test_cores_parsed(self) -> None:
        assert _parse(["-c", "12"]).cores == 12

    def test_output_parsed(self) -> None:
        assert _parse(["-o", "myout/"]).output == "myout/"

    def test_technology_parsed(self) -> None:
        assert _parse(["-x", "10xv2"]).technology == "10xv2"

    def test_detection_threshold_parsed(self) -> None:
        assert _parse(["--detection-threshold", "5"]).detection_threshold == 5

    def test_ncbi_accession_parsed(self) -> None:
        assert _parse(["-acc", "NC_002021.3"]).ncbi_accession == "NC_002021.3"

    def test_ncbi_email_parsed(self) -> None:
        assert _parse(["--ncbi-email", "a@b.com"]).ncbi_email == "a@b.com"

    def test_data_cache_dir_parsed(self) -> None:
        assert _parse(["--data-cache-dir", "/shared/viralscan-cache"]).data_cache_dir == (
            "/shared/viralscan-cache"
        )

    def test_hvg_min_mean_parsed(self) -> None:
        assert _parse(["--hvg-min-mean", "0.2"]).hvg_min_mean == 0.2

    def test_hvg_max_mean_parsed(self) -> None:
        assert _parse(["--hvg-max-mean", "2.5"]).hvg_max_mean == 2.5

    def test_hvg_min_disp_parsed(self) -> None:
        assert _parse(["--hvg-min-disp", "0.7"]).hvg_min_disp == 0.7

    def test_umap_n_neighbors_parsed(self) -> None:
        assert _parse(["--umap-n-neighbors", "21"]).umap_n_neighbors == 21

    def test_multimap_method_parsed(self) -> None:
        assert _parse(["--multimap-method", "unique-weighted"]).multimap_method == "unique-weighted"

    def test_multimap_primary_call_parsed(self) -> None:
        assert (
            _parse(["--multimap-primary-call", "unique-only"]).multimap_primary_call
            == "unique-only"
        )

    def test_multimap_pseudocount_parsed(self) -> None:
        assert _parse(["--multimap-pseudocount", "0.25"]).multimap_pseudocount == 0.25

    def test_invalid_multimap_method_rejected(self) -> None:
        with pytest.raises(SystemExit):
            _parse(["--multimap-method", "em"])

    def test_verbose_and_quiet_mutually_exclusive(self) -> None:
        with pytest.raises(SystemExit):
            _parse(["--verbose", "--quiet"])


class TestCommaSeparatedPaths:
    def test_split_comma_paths_trims_and_drops_empty_entries(self) -> None:
        from viralscan.utils import split_comma_paths

        assert split_comma_paths(" a.fastq.gz, ,b.fastq.gz,, c.fastq.gz ") == [
            "a.fastq.gz",
            "b.fastq.gz",
            "c.fastq.gz",
        ]

    def test_config_value_serializes_none_as_empty_string(self) -> None:
        from viralscan.menu import _config_value

        assert _config_value(None) == ""
        assert _config_value("custom.gtf") == "custom.gtf"


class TestBuildKbRefInputs:
    def test_multiple_reference_inputs_are_materialized_before_kb_ref(self, tmp_path) -> None:
        from viralscan.menu import _build_kb_ref

        fasta1 = tmp_path / "a.fa"
        fasta2 = tmp_path / "b.fa"
        gtf1 = tmp_path / "a.gtf"
        gtf2 = tmp_path / "b.gtf"
        fasta1.write_text(">A\nAAAA")
        fasta2.write_text(">B\nBBBB\n")
        gtf1.write_text('A\t.\tgene\t1\t4\t.\t+\t.\tgene_id "A";\n\n')
        gtf2.write_text('B\t.\tgene\t1\t4\t.\t+\t.\tgene_id "B";\n')

        calls = []

        def fake_run(cmd, check):
            calls.append(cmd)
            assert check is True
            return subprocess.CompletedProcess(cmd, 0)

        with patch("viralscan.menu.subprocess.run", side_effect=fake_run):
            _build_kb_ref(tmp_path / "out", f"{fasta1},{fasta2}", f"{gtf1},{gtf2}")

        assert len(calls) == 1
        cmd = calls[0]
        materialized_fasta = tmp_path / "out" / "index" / "input.fasta"
        materialized_gtf = tmp_path / "out" / "index" / "input.gtf"
        assert cmd[-2:] == [str(materialized_fasta), str(materialized_gtf)]
        assert materialized_fasta.read_text() == ">A\nAAAA\n>B\nBBBB\n"
        assert materialized_gtf.read_text() == (
            'A\t.\tgene\t1\t4\t.\t+\t.\tgene_id "A";\nB\t.\tgene\t1\t4\t.\t+\t.\tgene_id "B";\n'
        )


# ── build-ref subcommand ───────────────────────────────────────────────────────


class TestBuildRefSubcommand:
    def test_subcommand_detected(self) -> None:
        args = _parse(["build-ref"])
        assert args._subcommand == "build-ref"

    def test_host_parsed(self) -> None:
        args = _parse(["build-ref", "--host", "human"])
        assert args.host == "human"

    def test_virus_accessions_parsed(self) -> None:
        args = _parse(["build-ref", "--virus-accessions", "NC_045512.2", "NC_002021.1"])
        assert args.virus_accessions == ["NC_045512.2", "NC_002021.1"]

    def test_no_kb_ref_flag(self) -> None:
        args = _parse(["build-ref", "--no-kb-ref"])
        assert args.no_kb_ref is True

    def test_no_kb_ref_defaults_false(self) -> None:
        args = _parse(["build-ref"])
        assert args.no_kb_ref is False

    def test_list_species_flag(self) -> None:
        args = _parse(["build-ref", "--list-species"])
        assert args.list_species is True


# ── data subcommand ───────────────────────────────────────────────────────────


class TestDataSubcommand:
    def test_data_group_detected(self) -> None:
        args = _parse(["data"])
        assert args._subcommand == "data"

    def test_fetch_subcommand_detected(self) -> None:
        args = _parse(["data", "fetch"])
        assert args._subcommand == "data-fetch"

    def test_fetch_cache_dir_parsed(self) -> None:
        args = _parse(["data", "fetch", "--cache-dir", "/tmp/viralscan-cache"])
        assert args.cache_dir == "/tmp/viralscan-cache"

    def test_fetch_force_defaults_false(self) -> None:
        args = _parse(["data", "fetch"])
        assert args.force is False

    def test_fetch_force_parsed(self) -> None:
        args = _parse(["data", "fetch", "--force"])
        assert args.force is True

    def test_fetch_sha256_parsed(self) -> None:
        args = _parse(["data", "fetch", "--sha256", "abc123"])
        assert args.sha256 == "abc123"


# ── _has_valid_fastq_suffix ───────────────────────────────────────────────────


class TestHasValidFastqSuffix:
    def setup_method(self) -> None:
        from viralscan.menu import _has_valid_fastq_suffix

        self.fn = _has_valid_fastq_suffix

    @pytest.mark.parametrize(
        "path",
        ["sample.fastq", "sample.fq", "sample.fastq.gz", "sample.fq.gz"],
    )
    def test_valid_suffixes(self, path: str) -> None:
        assert self.fn(path) is True

    @pytest.mark.parametrize(
        "path",
        ["sample.bam", "sample.txt", "sample.fastq.bz2", ""],
    )
    def test_invalid_suffixes(self, path: str) -> None:
        assert self.fn(path) is False


# ── errorhandler validation (unit, no filesystem) ────────────────────────────


class TestErrorhandler:
    """Test errorhandler() branches using mocked path existence checks."""

    def _args(self, **kwargs):
        """Return a minimal Namespace with sensible defaults."""
        defaults = dict(
            output="out/",
            sample1="s1.fastq.gz",
            sample2="s2.fastq.gz",
            reference=False,
            gtf=None,
            fasta=None,
            f1=None,
            index="idx.idx",
            transcripts="t2g.txt",
            ncbi_accession=None,
        )
        defaults.update(kwargs)
        return argparse.Namespace(**defaults)

    def test_missing_index_calls_die(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args(index=None)
        with patch("os.path.exists", return_value=False), pytest.raises(SystemExit):
            errorhandler(args)

    def test_missing_transcripts_calls_die(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args(transcripts=None)
        with (
            patch("os.path.exists", side_effect=lambda p: p != "t2g.txt"),
            pytest.raises(SystemExit),
        ):
            errorhandler(args)

    def test_invalid_fastq_suffix_calls_die(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args(sample1="sample.bam", sample2="sample.bam")
        with patch("os.path.exists", return_value=True), pytest.raises(SystemExit):
            errorhandler(args)

    def test_sample_count_mismatch_calls_die(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args(sample1="a.fastq.gz,b.fastq.gz", sample2="c.fastq.gz")
        with patch("os.path.exists", return_value=True), pytest.raises(SystemExit):
            errorhandler(args)

    def test_ncbi_accession_with_reference_calls_die(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args(ncbi_accession="NC_002021.3", reference=True)
        with patch("os.path.exists", return_value=True), pytest.raises(SystemExit):
            errorhandler(args)

    def test_reference_without_gtf_calls_die(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args(reference=True, gtf=None, fasta="viral.fa")
        with patch("os.path.exists", return_value=True), pytest.raises(SystemExit):
            errorhandler(args)

    def test_reference_without_fasta_calls_die(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args(reference=True, fasta=None, gtf="viral.gtf")
        with patch("os.path.exists", return_value=True), pytest.raises(SystemExit):
            errorhandler(args)

    def test_valid_index_mode_passes(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args()
        with patch("os.path.exists", return_value=True):
            # should not raise
            errorhandler(args)

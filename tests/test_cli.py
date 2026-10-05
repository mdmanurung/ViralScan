"""Tests for the CLI (menu.py).

Covers argument parsing, boolean flag regressions, and validation helpers.
No network access; no subprocesses that touch the filesystem beyond tmp dirs.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
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
    def test_technology_default_is_detected(self) -> None:
        args = _parse([])
        assert args.technology is None
        assert args.force_technology is False

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

    def test_strand_defaults_none(self) -> None:
        assert _parse([]).strand is None

    @pytest.mark.parametrize("value", ["forward", "reverse", "unstranded"])
    def test_strand_valid_choices_accepted(self, value) -> None:
        assert _parse(["--strand", value]).strand == value

    def test_strand_bad_value_rejected(self) -> None:
        with pytest.raises(SystemExit) as exc:
            _parse(["--strand", "both"])
        assert exc.value.code == 2

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
            _parse(["--multimap-primary-call", "selected-method"]).multimap_primary_call
            == "selected-method"
        )

    def test_multimap_pseudocount_parsed(self) -> None:
        assert _parse(["--multimap-pseudocount", "0.25"]).multimap_pseudocount == 0.25

    def test_invalid_multimap_method_rejected(self) -> None:
        with pytest.raises(SystemExit):
            _parse(["--multimap-method", "bogus-method"])

    @pytest.mark.parametrize("method", ["em-global", "em-cell"])
    def test_em_multimap_method_accepted(self, method: str) -> None:
        assert _parse(["--multimap-method", method]).multimap_method == method

    def test_ambiguous_plain_em_name_rejected(self) -> None:
        with pytest.raises(SystemExit):
            _parse(["--multimap-method", "em"])

    def test_verbose_and_quiet_mutually_exclusive(self) -> None:
        with pytest.raises(SystemExit):
            _parse(["--verbose", "--quiet"])

    def test_resume_and_overwrite_are_mutually_exclusive(self) -> None:
        with pytest.raises(SystemExit):
            _parse(["--resume", "--overwrite"])

    def test_yes_does_not_select_an_output_mode(self) -> None:
        args = _parse(["--yes"])
        assert args.yes is True
        assert args.resume is False
        assert args.overwrite is False


class TestCommaSeparatedPaths:
    def test_split_comma_paths_trims_and_drops_empty_entries(self) -> None:
        from viralscan.utils import split_comma_paths

        assert split_comma_paths(" a.fastq.gz, ,b.fastq.gz,, c.fastq.gz ") == [
            "a.fastq.gz",
            "b.fastq.gz",
            "c.fastq.gz",
        ]


class TestBuildConfigArgs:
    """_build_config_args produces a correct Snakemake --config k=v list."""

    def _make_args(self, **overrides) -> argparse.Namespace:
        defaults = dict(
            cores=4,
            gtf=None,
            fasta=None,
            visual=True,
            reference=False,
            umap=False,
            technology="10xv3",
            whitelist=None,
            multimapping=True,
            se_threshold=10,
            detection_threshold=1,
            min_counts=1000,
            min_genes=200,
            hvg_min_mean=0.0125,
            hvg_max_mean=3.0,
            hvg_min_disp=0.5,
            umap_n_neighbors=15,
            multimap_method="equal",
            multimap_pseudocount=1.0,
            multimap_primary_call="selected-method",
            multimap_em_max_iter=100,
            multimap_em_tol=1e-6,
            cell_types=None,
            data_cache_dir=None,
            host_filter=None,
            host_index=None,
            positive_control_gene=None,
            positive_control_molecules=None,
            require_positive_control=False,
            anellovirus_gene_ids=True,
        )
        defaults.update(overrides)
        return argparse.Namespace(**defaults)

    def _as_dict(self, args: argparse.Namespace, **path_overrides) -> dict[str, str]:
        from viralscan.menu import _build_config_args

        kv_list = _build_config_args(
            args,
            outs=path_overrides.get("outs", "/out/sample/"),
            index=path_overrides.get("index", "/ref/index.idx"),
            transcripts=path_overrides.get("transcripts", "/ref/t2g.txt"),
            f1=path_overrides.get("f1"),
            s1=path_overrides.get("s1", "R1.fastq.gz"),
            s2=path_overrides.get("s2", "R2.fastq.gz"),
        )
        return dict(kv.split("=", 1) for kv in kv_list)

    def test_em_keys_present(self) -> None:
        d = self._as_dict(self._make_args())
        assert "multimap_em_max_iter" in d
        assert "multimap_em_tol" in d

    def test_em_keys_carry_cli_values(self) -> None:
        d = self._as_dict(self._make_args(multimap_em_max_iter=50, multimap_em_tol=1e-4))
        assert int(d["multimap_em_max_iter"]) == 50
        assert float(d["multimap_em_tol"]) == pytest.approx(1e-4)

    def test_booleans_are_lowercase(self) -> None:
        d = self._as_dict(self._make_args(visual=True, umap=False))
        assert d["visual"] == "true"
        assert d["umap"] == "false"

    def test_none_fields_emit_empty_value(self) -> None:
        d = self._as_dict(self._make_args(gtf=None, cell_types=None))
        assert d["gtf"] == ""
        assert d["cell_types"] == ""

    def test_host_filter_attr_maps_to_host_filter_aligner_key(self) -> None:
        d = self._as_dict(
            self._make_args(host_filter="starsolo", host_index="/path/to/index"),
        )
        assert d["host_filter_aligner"] == "starsolo"

    @pytest.mark.parametrize(
        ("chosen", "expected"), [("star-default", "star-default"), (None, "pinned")]
    )
    def test_star_params_reach_the_host_filter_config(self, tmp_path, chosen, expected) -> None:
        """CLI args -> snakemake --config -> config.yaml -> RunConfig, as host_filter.main reads it."""
        from viralscan.runconfig import RunConfig

        d = self._as_dict(
            self._make_args(
                host_filter="starsolo", host_index="/path/to/index", host_filter_star_params=chosen
            ),
            outs=f"{tmp_path}/",
        )
        RunConfig.from_snakemake_config(d).to_yaml(tmp_path / "config.yaml")
        assert RunConfig.from_yaml(tmp_path / "config.yaml").host_filter_star_params == expected

    @pytest.mark.parametrize(("chosen", "expected"), [("artefact", "artefact"), (None, "off")])
    def test_read_filter_reaches_the_config(self, tmp_path, chosen, expected) -> None:
        """CLI args -> snakemake --config -> config.yaml -> RunConfig (DEF-01)."""
        from viralscan.runconfig import RunConfig

        d = self._as_dict(self._make_args(read_filter=chosen), outs=f"{tmp_path}/")
        RunConfig.from_snakemake_config(d).to_yaml(tmp_path / "config.yaml")
        loaded = RunConfig.from_yaml(tmp_path / "config.yaml")
        assert loaded.read_filter == expected
        assert ("read_filtered" in loaded.kb_r1) is (expected == "artefact")


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
            # kb ref writes the t2g; the build manifest is derived from it.
            Path(cmd[cmd.index("-g") + 1]).write_text(
                "A\tA\t\t\tA\t1\t4\t+\nB\tB\t\t\tB\t1\t4\t+\n"
            )
            return subprocess.CompletedProcess(cmd, 0)

        with patch("viralscan.menu.subprocess.run", side_effect=fake_run):
            _build_kb_ref(tmp_path / "out", f"{fasta1},{fasta2}", f"{gtf1},{gtf2}")

        assert (tmp_path / "out" / "index" / "index.idx.build_manifest.json").is_file()
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

    def test_expanded_anellovirus_is_opt_in(self) -> None:
        assert _parse(["build-ref"]).anellovirus is False
        assert _parse(["build-ref", "--anellovirus"]).anellovirus is True

    def test_partial_panel_is_opt_in(self) -> None:
        assert _parse(["build-ref"]).allow_partial_panel is False
        assert _parse(["build-ref", "--allow-partial-panel"]).allow_partial_panel is True

    def test_genome_dlist_is_explicit(self) -> None:
        assert _parse(["build-ref"]).genome_dlist is None
        assert _parse(["build-ref", "--genome-dlist", "GRCh38.fa"]).genome_dlist == "GRCh38.fa"


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

    def test_star_params_without_host_filter_calls_die(self) -> None:
        from viralscan.menu import errorhandler

        args = self._args(host_filter_star_params="star-default")
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


# ── cell-caller preflight ─────────────────────────────────────────────────────


class TestCellCallerPreflight:
    """Cell calling fails closed and runs late, so a missing R must abort early.

    Without this, an environment without R completes kb_count, analysis, and
    multimap before dying for a reason that was knowable before it started.
    """

    def test_auto_without_an_external_list_resolves_to_emptydrops(self) -> None:
        from viralscan.menu import _resolve_cell_calling

        assert _resolve_cell_calling("auto", None) == "emptydrops"

    def test_auto_with_an_external_list_resolves_to_external(self) -> None:
        from viralscan.menu import _resolve_cell_calling

        assert _resolve_cell_calling("auto", "cells.txt") == "external"

    def test_missing_rscript_aborts_before_the_workflow(self) -> None:
        from viralscan.menu import _check_cell_caller_tools

        with patch("shutil.which", return_value=None), pytest.raises(SystemExit):
            _check_cell_caller_tools("emptydrops", "Rscript")

    def test_present_rscript_passes(self) -> None:
        from viralscan.menu import _check_cell_caller_tools

        with patch("shutil.which", return_value="/usr/bin/Rscript"):
            _check_cell_caller_tools("emptydrops", "Rscript")

    def test_other_methods_do_not_require_r(self) -> None:
        from viralscan.menu import _check_cell_caller_tools

        with patch("shutil.which", return_value=None):
            for method in ("none", "knee", "external"):
                _check_cell_caller_tools(method, "Rscript")

    def test_emptydrops_seed_default_matches_defaults(self) -> None:
        assert _parse([]).emptydrops_seed == DEFAULTS["emptydrops_seed"]

    def test_emptydrops_seed_is_settable(self) -> None:
        assert _parse(["--emptydrops-seed", "20260727002"]).emptydrops_seed == 20260727002


class TestRerunRunManifest:
    """The manifest sits at the ROOT of the result tree, above the per-sample
    directories: `prepare_output_directory` writes it into `--output`, and
    `createconfig` then creates one subdirectory per sample beneath it.

    An earlier version of these tests asserted the opposite, following a review
    finding that claimed the manifest lived beside `config.yaml`. Running the
    pipeline end to end (SW-10) showed `out/run_manifest.json` alongside
    `out/<sample>/config.yaml`, so the finding and the tests were both wrong.
    """

    @staticmethod
    def _run_root(tmp_path):
        root = tmp_path / "out"
        (root / "sampleA").mkdir(parents=True)
        (root / "sampleA" / "config.yaml").write_text("x\n", encoding="utf-8")
        (root / "run_manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": "3.0.0",
                    "allocation_method": "equal",
                    "run_fingerprint": "old-fingerprint",
                }
            ),
            encoding="utf-8",
        )
        return root

    def test_manifest_is_restamped_for_the_new_method(self, tmp_path):
        from viralscan.menu import _rewrite_run_manifest

        root = self._run_root(tmp_path)

        assert _rewrite_run_manifest(
            root, source_dir=tmp_path / "src", new_method="host-conservative"
        )

        written = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
        assert written["allocation_method"] == "host-conservative"
        assert written["parent_run_fingerprint"] == "old-fingerprint"
        assert written["run_fingerprint"] != "old-fingerprint"
        assert written["derived_from"] == str(tmp_path / "src")

    def test_a_sample_directory_holds_no_manifest(self, tmp_path):
        """Guards the regression: rewriting per sample finds nothing to rewrite."""
        from viralscan.menu import _rewrite_run_manifest

        root = self._run_root(tmp_path)

        assert (
            _rewrite_run_manifest(
                root / "sampleA", source_dir=tmp_path / "src", new_method="host-conservative"
            )
            is False
        )

    def test_no_temp_file_survives(self, tmp_path):
        from viralscan.menu import _rewrite_run_manifest

        root = self._run_root(tmp_path)
        _rewrite_run_manifest(root, source_dir=tmp_path, new_method="unique-weighted")

        assert list(root.glob(".*tmp")) == []


class TestSnakemakeRunCommand:
    """SW-14: one snakemake invocation shared by ``main`` and ``rerun-multimap``."""

    def test_target_precedes_quiet(self) -> None:
        """snakemake 9's ``--quiet [{all,...} ...]`` consumes a following ``all``."""
        from viralscan.menu import _snakemake_run_command

        cmd = _snakemake_run_command("/x/Snakefile", 4, ["output=/o/"])

        assert cmd.index("all") < cmd.index("--quiet")

    def test_does_not_require_conda(self) -> None:
        from viralscan.menu import _snakemake_run_command

        assert "--use-conda" not in _snakemake_run_command("/x/Snakefile", 4, [])

    def test_config_args_come_last(self) -> None:
        from viralscan.menu import _snakemake_run_command

        cmd = _snakemake_run_command("/x/Snakefile", 2, ["a=1", "b=2"])

        assert cmd[-3:] == ["--config", "a=1", "b=2"]
        assert cmd[cmd.index("--cores") + 1] == "2"


# ── chemistry check (PLAN DEF-02, WP1E Q6) ───────────────────────────────────


class TestResolveChemistry:
    def _args(self, **kw) -> argparse.Namespace:
        base = dict(
            sample1="a/S1_R1.fastq.gz", whitelist=None, technology=None, force_technology=False
        )
        return argparse.Namespace(**{**base, **kw})

    def _detection(self, chem, reason="97.3% on the 10xv2 list"):
        from viralscan.chemistry import Detection

        return Detection(chem, "bundled on-list", reason, 26, None, 10, 100, {"10xv2": 0.973})

    def test_auto_sets_the_detected_technology_and_returns_blocks(self, monkeypatch) -> None:
        from viralscan import chemistry, menu

        monkeypatch.setattr(chemistry, "detect", lambda paths, wl: [self._detection("10xv2")])
        args = self._args()
        blocks = menu._resolve_chemistry(args)
        assert args.technology == "10xv2"
        assert blocks["S1"]["chemistry"] == "10xv2"

    def test_contradicted_x_stops_the_run(self, monkeypatch) -> None:
        from viralscan import chemistry, menu

        monkeypatch.setattr(chemistry, "detect", lambda paths, wl: [self._detection("10xv2")])
        with pytest.raises(SystemExit):
            menu._resolve_chemistry(self._args(technology="10xv3"))

    def test_force_keeps_x_even_when_detection_cannot_run(self, monkeypatch) -> None:
        from viralscan import chemistry, menu

        def boom(paths, wl):
            raise chemistry.ChemistryError("ngs_tools missing")

        monkeypatch.setattr(chemistry, "detect", boom)
        args = self._args(technology="10xv3", force_technology=True)
        assert menu._resolve_chemistry(args) == {}
        assert args.technology == "10xv3"


def test_missing_snakemake_error_names_full_extra(caplog: pytest.LogCaptureFixture) -> None:
    from viralscan import menu

    with patch("shutil.which", side_effect=lambda t: None if t == "snakemake" else "/bin/" + t):
        with pytest.raises(SystemExit) as exc:
            menu._check_required_tools()
    assert exc.value.code != 0
    err = caplog.text
    assert "snakemake" in err and "viralscan[full]" in err

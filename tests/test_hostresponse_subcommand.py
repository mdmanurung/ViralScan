"""Tests for the 'viralscan hostresponse' subcommand (parser + dispatch)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


# ── Parser tests (no heavy deps) ──────────────────────────────────────────────


class TestHostresponseParser:
    def test_help_exits_zero(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "--help"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit) as exc:
                create_help()
        assert exc.value.code == 0

    def test_parses_required_args(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "out/", "--host-h5ad", "host.h5ad"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args._subcommand == "hostresponse"
        assert args.output == "out/"
        assert args.host_h5ad == "host.h5ad"

    def test_requires_output(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "--host-h5ad", "host.h5ad"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()

    def test_requires_host_h5ad(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "out/"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()

    def test_no_use_hvg_flag(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "out/", "--host-h5ad", "h.h5ad", "--no-use-hvg"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.use_hvg is False

    def test_use_hvg_default_true(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "out/", "--host-h5ad", "h.h5ad"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.use_hvg is True

    def test_n_seeds_override(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "out/", "--host-h5ad", "h.h5ad", "--n-seeds", "3"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.n_seeds == 3

    def test_n_seeds_default_none(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "out/", "--host-h5ad", "h.h5ad"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.n_seeds is None

    def test_enrichment_flag(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "out/", "--host-h5ad", "h.h5ad", "--enrichment"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.enrichment is True

    def test_enrichment_default_false(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "out/", "--host-h5ad", "h.h5ad"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.enrichment is False

    def test_stab_min_prob_override(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "o/", "--host-h5ad", "h.h5ad", "--stab-min-prob", "0.75"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.stab_min_prob == pytest.approx(0.75)

    def test_detection_threshold_override(self) -> None:
        with patch("sys.argv", ["viralscan", "hostresponse", "-o", "o/", "--host-h5ad", "h.h5ad", "--detection-threshold", "3"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.detection_threshold == 3


# ── Dispatch tests ────────────────────────────────────────────────────────────


def _make_sample_dir(tmp_path: Path, output_suffix: str = "/") -> Path:
    """Create a minimal viralscan output dir with config.yaml and analysis.txt."""
    from viralscan.runconfig import RunConfig

    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()
    log_dir = sample_dir / "log"
    log_dir.mkdir()

    cfg = RunConfig(
        output=str(sample_dir) + "/",
        multimapping=True,
        detection_threshold=1,
    )
    cfg.to_yaml(sample_dir / "config.yaml")
    (log_dir / "analysis.txt").write_text("NC_001806.2\n")
    return sample_dir


def _make_args(sample_dir: Path, **overrides) -> MagicMock:
    """Build a minimal argparse Namespace for _run_hostresponse_subcommand."""
    args = MagicMock()
    args.output = str(sample_dir)
    args.host_h5ad = "host.h5ad"
    args.n_seeds = None
    args.n_stab_iter = None
    args.use_hvg = True
    args.stab_min_prob = None
    args.top_n_genes = None
    args.detection_threshold = None
    args.enrichment = False
    args.enrichment_db = None
    args.verbose = False
    args.quiet = False
    for k, v in overrides.items():
        setattr(args, k, v)
    return args


class TestRunHostresponseSubcommand:
    def test_calls_run_hostresponse(self, tmp_path: Path) -> None:
        sample_dir = _make_sample_dir(tmp_path)
        fake_h5ad = sample_dir / "virus.h5ad"
        fake_h5ad.touch()

        mock_kb = MagicMock()
        mock_kb.current_adata.return_value = fake_h5ad

        with (
            patch("viralscan.kb_outputs.KbCountOutputs.from_config_output", return_value=mock_kb),
            patch("viralscan.scripts.hostresponse.run_hostresponse") as mock_run,
        ):
            from viralscan.menu import _run_hostresponse_subcommand

            _run_hostresponse_subcommand(_make_args(sample_dir))

        mock_run.assert_called_once()

    def test_passes_host_h5ad_to_run_hostresponse(self, tmp_path: Path) -> None:
        sample_dir = _make_sample_dir(tmp_path)
        fake_h5ad = sample_dir / "virus.h5ad"
        fake_h5ad.touch()

        mock_kb = MagicMock()
        mock_kb.current_adata.return_value = fake_h5ad

        with (
            patch("viralscan.kb_outputs.KbCountOutputs.from_config_output", return_value=mock_kb),
            patch("viralscan.scripts.hostresponse.run_hostresponse") as mock_run,
        ):
            from viralscan.menu import _run_hostresponse_subcommand

            _run_hostresponse_subcommand(_make_args(sample_dir, host_h5ad="custom_host.h5ad"))

        _, call_kwargs = mock_run.call_args
        assert call_kwargs["host_h5ad"] == "custom_host.h5ad"

    def test_out_dir_is_under_sample_dir(self, tmp_path: Path) -> None:
        sample_dir = _make_sample_dir(tmp_path)
        fake_h5ad = sample_dir / "virus.h5ad"
        fake_h5ad.touch()

        mock_kb = MagicMock()
        mock_kb.current_adata.return_value = fake_h5ad

        with (
            patch("viralscan.kb_outputs.KbCountOutputs.from_config_output", return_value=mock_kb),
            patch("viralscan.scripts.hostresponse.run_hostresponse") as mock_run,
        ):
            from viralscan.menu import _run_hostresponse_subcommand

            _run_hostresponse_subcommand(_make_args(sample_dir))

        _, call_kwargs = mock_run.call_args
        assert call_kwargs["out_dir"].endswith("hostresponse")

    def test_cli_n_seeds_overrides_config(self, tmp_path: Path) -> None:
        sample_dir = _make_sample_dir(tmp_path)
        fake_h5ad = sample_dir / "virus.h5ad"
        fake_h5ad.touch()

        mock_kb = MagicMock()
        mock_kb.current_adata.return_value = fake_h5ad

        with (
            patch("viralscan.kb_outputs.KbCountOutputs.from_config_output", return_value=mock_kb),
            patch("viralscan.scripts.hostresponse.run_hostresponse") as mock_run,
        ):
            from viralscan.menu import _run_hostresponse_subcommand

            _run_hostresponse_subcommand(_make_args(sample_dir, n_seeds=2))

        _, call_kwargs = mock_run.call_args
        assert len(call_kwargs["seeds"]) == 2

    def test_enrichment_flag_forwarded(self, tmp_path: Path) -> None:
        sample_dir = _make_sample_dir(tmp_path)
        fake_h5ad = sample_dir / "virus.h5ad"
        fake_h5ad.touch()

        mock_kb = MagicMock()
        mock_kb.current_adata.return_value = fake_h5ad

        with (
            patch("viralscan.kb_outputs.KbCountOutputs.from_config_output", return_value=mock_kb),
            patch("viralscan.scripts.hostresponse.run_hostresponse") as mock_run,
        ):
            from viralscan.menu import _run_hostresponse_subcommand

            _run_hostresponse_subcommand(_make_args(sample_dir, enrichment=True))

        _, call_kwargs = mock_run.call_args
        assert call_kwargs["do_enrichment"] is True

    def test_dies_when_config_yaml_missing(self, tmp_path: Path) -> None:
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()

        with pytest.raises(SystemExit) as exc:
            from viralscan.menu import _run_hostresponse_subcommand

            _run_hostresponse_subcommand(_make_args(empty_dir))
        assert exc.value.code == 1

    def test_dies_when_analysis_txt_missing(self, tmp_path: Path) -> None:
        from viralscan.runconfig import RunConfig

        sample_dir = tmp_path / "sample"
        sample_dir.mkdir()
        cfg = RunConfig(output=str(sample_dir) + "/")
        cfg.to_yaml(sample_dir / "config.yaml")
        # log/analysis.txt deliberately absent

        with pytest.raises(SystemExit) as exc:
            from viralscan.menu import _run_hostresponse_subcommand

            _run_hostresponse_subcommand(_make_args(sample_dir))
        assert exc.value.code == 1

    def test_dies_when_virus_h5ad_missing(self, tmp_path: Path) -> None:
        sample_dir = _make_sample_dir(tmp_path)
        # virus h5ad does NOT exist

        mock_kb = MagicMock()
        mock_kb.current_adata.return_value = sample_dir / "nonexistent.h5ad"

        with (
            patch("viralscan.kb_outputs.KbCountOutputs.from_config_output", return_value=mock_kb),
            pytest.raises(SystemExit) as exc,
        ):
            from viralscan.menu import _run_hostresponse_subcommand

            _run_hostresponse_subcommand(_make_args(sample_dir))
        assert exc.value.code == 1

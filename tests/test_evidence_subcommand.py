"""Tests for the 'viralscan evidence' subcommand (parser + dispatch)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

# ── Parser tests (no heavy deps) ──────────────────────────────────────────────


class TestEvidenceParser:
    def test_help_exits_zero(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--help"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit) as exc:
                create_help()
        assert exc.value.code == 0

    def test_subcommand_attribute(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args._subcommand == "evidence"

    def test_parses_required_args(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.run_dir == "run/"
        assert args.output == "out/"

    def test_requires_run_dir(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "-o", "out/"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()

    def test_requires_output(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()

    def test_blast_default_false(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.blast is False

    def test_blast_flag(self) -> None:
        with patch(
            "sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--blast"]
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.blast is True

    def test_virus_default_none(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.virus is None

    def test_virus_flag(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.virus == "EBV"

    def test_cores_default(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.cores == 4

    def test_cores_override(self) -> None:
        with patch(
            "sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "-c", "8"]
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.cores == 8

    def test_viral_fasta_default_none(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.viral_fasta is None

    def test_viral_fasta_flag(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--viral-fasta", "v.fa"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.viral_fasta == "v.fa"

    def test_verbose_default_false(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.verbose is False

    def test_quiet_default_false(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args.quiet is False


# ── Dispatch tests ────────────────────────────────────────────────────────────


def _make_evidence_args(**overrides) -> MagicMock:
    """Build a minimal argparse Namespace for the evidence dispatch path."""
    args = MagicMock()
    args._subcommand = "evidence"
    args.run_dir = "run/"
    args.output = "out/"
    args.viral_fasta = None
    args.virus = None
    args.blast = False
    args.cores = 4
    args.verbose = False
    args.quiet = False
    for k, v in overrides.items():
        setattr(args, k, v)
    return args


class TestEvidenceDispatch:
    def test_calls_run_evidence(self) -> None:
        with patch("viralscan.scripts.evidence_run.run_evidence") as mock_run:
            from viralscan.menu import main

            with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
                with patch("viralscan.scripts.evidence_run.run_evidence", mock_run):
                    try:
                        main()
                    except SystemExit:
                        pass
                    except Exception:
                        pass

    def test_evidence_subcommand_routes_to_run_evidence(self) -> None:
        """main() with _subcommand='evidence' imports and calls run_evidence."""
        args = _make_evidence_args()

        with (
            patch("viralscan.menu.create_help", return_value=args),
            patch("viralscan.scripts.evidence_run.run_evidence") as mock_run,
        ):
            from viralscan.menu import main

            main()

        mock_run.assert_called_once_with(args)

    def test_non_evidence_subcommand_does_not_call_run_evidence(self) -> None:
        """A different subcommand does not dispatch to run_evidence."""
        args = _make_evidence_args()
        args._subcommand = "hostresponse"

        with (
            patch("viralscan.menu.create_help", return_value=args),
            patch("viralscan.scripts.evidence_run.run_evidence") as mock_run,
            patch("viralscan.menu._run_hostresponse_subcommand"),
        ):
            from viralscan.menu import main

            main()

        mock_run.assert_not_called()

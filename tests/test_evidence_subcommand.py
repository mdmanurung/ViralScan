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
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args._subcommand == "evidence"

    def test_parses_required_args(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.run_dir == "run/"
        assert args.output == "out/"

    def test_requires_run_dir(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "-o", "out/", "--virus", "EBV"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()

    def test_requires_output(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "--virus", "EBV"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()

    def test_blast_default_false(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.blast is False

    def test_blast_flag(self) -> None:
        with patch(
            "sys.argv",
            [
                "viralscan",
                "evidence",
                "--run-dir",
                "run/",
                "-o",
                "out/",
                "--virus",
                "EBV",
                "--blast",
            ],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.blast is True

    def test_virus_is_required(self) -> None:
        with patch("sys.argv", ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()

    def test_virus_flag(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.virus == "EBV"

    def test_cores_default(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.cores == 4

    def test_cores_override(self) -> None:
        with patch(
            "sys.argv",
            [
                "viralscan",
                "evidence",
                "--run-dir",
                "run/",
                "-o",
                "out/",
                "--virus",
                "EBV",
                "-c",
                "8",
            ],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.cores == 8

    def test_viral_fasta_default_none(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.viral_fasta is None

    def test_viral_fasta_flag(self) -> None:
        with patch(
            "sys.argv",
            [
                "viralscan",
                "evidence",
                "--run-dir",
                "run/",
                "-o",
                "out/",
                "--virus",
                "EBV",
                "--viral-fasta",
                "v.fa",
            ],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.viral_fasta == "v.fa"

    def test_verbose_default_false(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args.verbose is False

    def test_quiet_default_false(self) -> None:
        with patch(
            "sys.argv",
            ["viralscan", "evidence", "--run-dir", "run/", "-o", "out/", "--virus", "EBV"],
        ):
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
        """The real argv parser must route `evidence` to run_evidence exactly once."""
        from viralscan.menu import main

        with (
            patch(
                "sys.argv",
                [
                    "viralscan",
                    "evidence",
                    "--run-dir",
                    "run/",
                    "-o",
                    "out/",
                    "--virus",
                    "EBV",
                ],
            ),
            patch("viralscan.scripts.evidence_run.run_evidence") as mock_run,
        ):
            main()

        mock_run.assert_called_once()
        dispatched = mock_run.call_args.args[0]
        assert dispatched.run_dir == "run/"
        assert dispatched.virus == "EBV"

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


# ── Replay input: the reads kb count actually quantified ─────────────────────


class TestReplayFastqs:
    """Exact-lineage replay must re-read the FASTQs ``kb count`` consumed.

    With a host filter active those are ``host_filtered/R{1,2}.fastq.gz``, not
    the raw ``sample1``/``sample2``: replaying the raw files re-counts reads the
    filter removed, so the evidence would describe molecules that were never in
    the matrix.
    """

    def test_uses_host_filtered_reads_when_present(self) -> None:
        from types import SimpleNamespace

        from viralscan.scripts.evidence_run import replay_fastqs

        config = SimpleNamespace(
            sample1="raw_R1.fq.gz",
            sample2="raw_R2.fq.gz",
            kb_r1="out/host_filtered/R1.fastq.gz",
            kb_r2="out/host_filtered/R2.fastq.gz",
        )
        assert replay_fastqs(config) == (
            "out/host_filtered/R1.fastq.gz",
            "out/host_filtered/R2.fastq.gz",
        )

    def test_falls_back_to_raw_reads_for_configs_without_kb_paths(self) -> None:
        from types import SimpleNamespace

        from viralscan.scripts.evidence_run import replay_fastqs

        config = SimpleNamespace(sample1="raw_R1.fq.gz", sample2="raw_R2.fq.gz", kb_r1="", kb_r2="")
        assert replay_fastqs(config) == ("raw_R1.fq.gz", "raw_R2.fq.gz")

    def _run_stubbed_evidence(self, tmp_path):
        """Run run_evidence with every filesystem/tool seam stubbed.

        Returns (mock_replay, mock_extract, out_dir). The run layout mirrors a
        real one: --run-dir is the per-sample directory (holds config.yaml and
        kb-python/), and run_manifest.json lives one level up at the run root.
        """
        from types import SimpleNamespace

        from viralscan.scripts import evidence_run

        run_dir = tmp_path / "run"
        run_dir.mkdir()
        (run_dir / "config.yaml").write_text("technology: 10xv3\n")
        kb_dir = run_dir / "kb-python"
        kb_dir.mkdir()
        for name in ("resolved.bus.txt", "genes.txt", "transcripts.txt", "kb.ec"):
            (kb_dir / name).write_text("")
        kb_r1 = tmp_path / "out" / "host_filtered" / "R1.fastq.gz"
        kb_r2 = tmp_path / "out" / "host_filtered" / "R2.fastq.gz"
        kb_r1.parent.mkdir(parents=True)
        kb_r1.write_bytes(b"")
        kb_r2.write_bytes(b"")

        index_path = tmp_path / "index.idx"
        index_path.write_bytes(b"")
        config = SimpleNamespace(
            sample1="raw_R1.fq.gz",
            sample2="raw_R2.fq.gz",
            kb_r1=str(kb_r1),
            kb_r2=str(kb_r2),
            technology="10xv3",
            transcripts=str(kb_dir / "transcripts.txt"),
            index=str(index_path),
            whitelist=None,
            strand=None,
            multimap_method="equal",
        )
        kb = SimpleNamespace(
            resolved_bus_txt=kb_dir / "resolved.bus.txt",
            genes=kb_dir / "genes.txt",
            transcripts_txt=kb_dir / "transcripts.txt",
            ec=kb_dir / "kb.ec",
            root=kb_dir,
        )
        flagged = tmp_path / "flagged.txt"
        flagged.write_text("")
        identity = SimpleNamespace(viral_gene_ids=lambda: {"g1"})
        stats = SimpleNamespace(viral_reads=0, total_reads=0)

        out_dir = tmp_path / "evidence"
        args = SimpleNamespace(
            run_dir=str(run_dir),
            output=str(out_dir),
            virus="EBV",
            cores=1,
            verbose=False,
            quiet=True,
            viral_fasta=None,
        )

        with (
            patch.object(evidence_run.RunConfig, "from_yaml", return_value=config),
            patch.object(evidence_run.KbCountOutputs, "from_config_output", return_value=kb),
            patch.object(evidence_run, "load_transcripts", return_value=(["tx1"], {"tx1": "g1"})),
            patch.object(evidence_run, "_replay_ec_map", return_value={}),
            patch.object(evidence_run, "load_run_identity", return_value=identity),
            patch.object(
                evidence_run, "resolve_viral_target", return_value=("EBV", {"g1"})
            ),
            patch.object(
                evidence_run, "replay_exact_target_bus", return_value=flagged
            ) as mock_replay,
            patch.object(evidence_run, "parse_flagged_target_bus", return_value={}),
            patch.object(
                evidence_run, "extract_exact_reads_by_number", return_value=stats
            ) as mock_extract,
        ):
            evidence_run.run_evidence(args)
        return mock_replay, mock_extract, out_dir, (kb_r1, kb_r2)

    def test_replay_and_extraction_both_use_the_same_reads(self, tmp_path) -> None:
        """Replay and extraction must consume the same FASTQs — the pair
        ``replay_fastqs`` resolves (host-filtered when present), never the raw
        ``config.sample1``/``sample2``."""
        mock_replay, mock_extract, _out, (kb_r1, kb_r2) = self._run_stubbed_evidence(tmp_path)

        replay_reads = (mock_replay.call_args.kwargs["r1_path"], mock_replay.call_args.kwargs["r2_path"])
        extract_reads = (mock_extract.call_args.args[0], mock_extract.call_args.args[1])
        assert replay_reads == (str(kb_r1), str(kb_r2))
        assert extract_reads == replay_reads

    def test_run_fingerprint_is_read_from_the_run_root_manifest(self, tmp_path) -> None:
        """B5: run_manifest.json lives at the run ROOT, one level above the
        per-sample --run-dir; the evidence manifest must still record its
        run_fingerprint instead of null."""
        import json

        (tmp_path / "run_manifest.json").write_text(
            json.dumps({"run_fingerprint": "fp-from-run-root"})
        )
        _replay, _extract, out_dir, _reads = self._run_stubbed_evidence(tmp_path)

        manifest = json.loads((out_dir / "evidence_manifest.json").read_text())
        assert manifest["run_fingerprint"] == "fp-from-run-root"

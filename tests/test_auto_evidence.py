"""ANDET-04: opt-in ``auto_evidence`` rule (config, validation, genus selection)."""

from __future__ import annotations

import argparse
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from viralscan.menu import _write_run_config, errorhandler
from viralscan.runconfig import RunConfig
from viralscan.scripts.auto_evidence import anello_genera
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable


def _cfg(tmp_path: Path, **extra) -> dict:
    return {
        "output": str(tmp_path / "out") + "/",
        "index": str(tmp_path / "ref" / "index.idx"),
        "transcripts": "t",
        "sample1": "s1",
        "sample2": "s2",
        "technology": "10xv3",
        "visual": "True",
        "multimapping": "True",
        "gtf": "None",
        "fasta": "None",
        "f1": "None",
        "reference": "False",
        "umap": "False",
        "whitelist": "None",
        "emptydrops_seed": 100,
        **extra,
    }


@pytest.fixture
def fastas(tmp_path):
    (tmp_path / "ref").mkdir()
    host = tmp_path / "host.fa"
    viral = tmp_path / "ref" / "viral.fa"
    host.write_text(">h\nACGT\n")
    viral.write_text(">v\nACGT\n")
    return host, viral


def test_default_is_off_and_writes_typed_config(tmp_path):
    rc = RunConfig.from_snakemake_config(_cfg(tmp_path))
    assert rc.auto_evidence is False and rc.host_fasta is None and rc.viral_fasta is None


def test_enabled_config_is_typed_and_resolves_viral_fasta_beside_the_index(tmp_path, fastas):
    host, viral = fastas
    rc = RunConfig.from_snakemake_config(_cfg(tmp_path, auto_evidence="true", host_fasta=str(host)))
    assert rc.auto_evidence is True
    assert rc.viral_fasta == str(viral.resolve())
    path = _write_run_config(rc)
    written = yaml.safe_load(path.read_text())
    assert written["auto_evidence"] is True
    assert written["host_fasta"] == str(host)
    assert RunConfig.from_yaml(path).auto_evidence is True


def test_runconfig_fails_closed(tmp_path, fastas):
    host, viral = fastas
    with pytest.raises(ValueError, match="host_fasta"):
        RunConfig.from_snakemake_config(_cfg(tmp_path, auto_evidence="true"))
    with pytest.raises(ValueError, match="multimapping"):
        RunConfig.from_snakemake_config(
            _cfg(tmp_path, auto_evidence="true", host_fasta=str(host), multimapping="False")
        )
    viral.unlink()  # no viral-only FASTA: never fall back to a host-bearing one
    with pytest.raises(ValueError, match="viral-only"):
        RunConfig.from_snakemake_config(_cfg(tmp_path, auto_evidence="true", host_fasta=str(host)))


def _args(tmp_path: Path, **kwargs) -> argparse.Namespace:
    (tmp_path / "R1.fastq.gz").touch()
    (tmp_path / "R2.fastq.gz").touch()
    defaults = {
        "ncbi_accession": None,
        "reference": False,
        "gtf": None,
        "fasta": None,
        "index": str(tmp_path / "ref" / "index.idx"),
        "transcripts": "/fake/t2g.txt",
        "f1": None,
        "host_filter": None,
        "host_index": None,
        "sample1": str(tmp_path / "R1.fastq.gz"),
        "sample2": str(tmp_path / "R2.fastq.gz"),
        "multimapping": True,
        "auto_evidence": True,
        "host_fasta": None,
        "viral_fasta": None,
    }
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def test_validation_requires_host_fasta(tmp_path, fastas):
    with patch("os.path.exists", return_value=True), pytest.raises(SystemExit) as exc:
        errorhandler(_args(tmp_path))
    assert exc.value.code != 0


@pytest.mark.parametrize(
    "bad", [{"multimapping": False}, {"host_fasta": "/nonexistent.fa"}, {"viral_fasta": "/no.fa"}]
)
def test_validation_rejects_unusable_inputs(tmp_path, fastas, bad):
    host, _ = fastas
    args = _args(tmp_path, host_fasta=str(host))
    vars(args).update(bad)
    with (
        patch("os.path.exists", return_value=True),
        patch("shutil.which", return_value="/bin/x"),
        pytest.raises(SystemExit),
    ):
        errorhandler(args)


def test_validation_accepts_a_complete_opt_in(tmp_path, fastas):
    host, _ = fastas
    with patch("os.path.exists", return_value=True), patch("shutil.which", return_value="/bin/x"):
        errorhandler(_args(tmp_path, host_fasta=str(host)))
    with patch("os.path.exists", return_value=True), patch("shutil.which", return_value=None):
        with pytest.raises(SystemExit):
            errorhandler(_args(tmp_path, host_fasta=str(host)))


def test_default_invocation_keeps_its_run_fingerprint():
    from viralscan.menu import create_help

    with patch("sys.argv", ["viralscan"]):
        args = create_help()
    assert args.auto_evidence is None and args.host_fasta is None and args.viral_fasta is None
    from viralscan.run_safety import build_run_manifest

    options = build_run_manifest(args)["options"]
    assert not {"auto_evidence", "host_fasta", "viral_fasta"} & set(options)


def _gene(gene_id, name, family):
    return GeneIdentity(
        gene_id=gene_id,
        genome_accession=gene_id + ".1",
        status="catalogued",
        viral=True,
        virus_key=f"genus:{name}",
        virus_name=name,
        family=family,
    )


def test_anello_genera_keeps_detected_anelloviridae_and_drops_alignment_only(tmp_path):
    identity = VirusIdentityTable(
        genes=(
            _gene("A", "Alphatorquevirus", "Anelloviridae"),
            _gene("B", "Betatorquevirus", "Anelloviridae"),
            _gene("E", "Human gammaherpesvirus 4", "Orthoherpesviridae"),
        )
    )
    summary = tmp_path / "viral_summary.tsv"
    summary.write_text(
        "virus_name\tdetection_source\n"
        "Human gammaherpesvirus 4\tkallisto\n"
        "Alphatorquevirus\tkallisto+alignment\n"
        "Betatorquevirus\talignment_only\n"
    )
    assert anello_genera(summary, identity) == ["Alphatorquevirus"]

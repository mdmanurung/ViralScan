"""CONTRACT-FP-01: an unset option equals its DEFAULTS value; an option with no default adds no key."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

import viralscan.run_safety as rs
from viralscan.defaults import DEFAULTS


@pytest.fixture
def make_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(rs, "__version__", "0.0.test")
    monkeypatch.chdir(tmp_path)
    for name, text in [("R1.fastq", "r1\n"), ("R2.fastq", "r2\n"), ("idx.idx", "x\n")]:
        Path(name).write_text(text)
    Path("t2g.txt").write_text("a\tb\n")

    def _make(**extra: object) -> argparse.Namespace:
        values = dict(
            sample1="R1.fastq", sample2="R2.fastq", index="idx.idx", transcripts="t2g.txt",
            called_cells_file=None, multimap_method="equal", cell_calling="emptydrops",
            technology="10xv3", output="o", resume=False, overwrite=False, yes=False,
            verbose=False, quiet=False, cores=4, ncbi_api_key=None, ncbi_email=None,
        )  # fmt: skip
        values.update(extra)
        return argparse.Namespace(**values)

    return _make


def _fp(args: argparse.Namespace) -> str:
    return rs.build_run_manifest(args)["run_fingerprint"]


def test_unset_equals_the_defaults_value(make_args) -> None:
    key = "programme_productive_min_molecules"
    assert _fp(make_args(**{key: None})) == _fp(make_args(**{key: DEFAULTS[key]}))


def test_productive_min_molecules_above_default_changes_the_fingerprint(make_args) -> None:
    assert _fp(make_args(programme_productive_min_molecules=2)) != _fp(
        make_args(programme_productive_min_molecules=None)
    )


def test_reference_host_species_is_hashed_only_when_set(make_args) -> None:
    base = _fp(make_args(reference_host_species=None))
    human = _fp(make_args(reference_host_species="Homo sapiens"))
    assert human != base
    assert human == _fp(make_args(reference_host_species="Homo sapiens"))
    assert human != _fp(make_args(reference_host_species="Mus musculus"))
    assert (
        "reference_host_species"
        not in rs.build_run_manifest(make_args(reference_host_species=None))["options"]
    )

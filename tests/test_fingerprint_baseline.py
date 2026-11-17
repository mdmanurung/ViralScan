"""CONTRACT-FP-01: one fingerprint rule for `viralscan run`.

An unset option takes its DEFAULTS value, so an explicit default and a missing flag hash the
same. The ViralScan version is hashed, so --resume across versions is refused. Resource and
helper options (`cores`, `kb_params`) are recorded and never hashed.
"""

import argparse
import json
from pathlib import Path

import pytest

import viralscan.run_safety as rs
from viralscan.cli.parsers import _build_run_parser
from viralscan.defaults import DEFAULTS
from viralscan.run_safety import RunSafetyError, prepare_output_directory


def _run_defaults() -> dict:
    sub = argparse.ArgumentParser().add_subparsers()
    _build_run_parser(sub)
    return {a.dest: a.default for a in sub.choices["run"]._actions if a.dest != "help"}


def _args(tmp_path: Path, **extra) -> argparse.Namespace:
    for name, text in [("R1.fastq", "r1\n"), ("R2.fastq", "r2\n"), ("idx.idx", "x\n")]:
        (tmp_path / name).write_text(text)
    (tmp_path / "t2g.txt").write_text("a\tb\n")
    base = _run_defaults()
    base.update(
        sample1=str(tmp_path / "R1.fastq"),
        sample2=str(tmp_path / "R2.fastq"),
        index=str(tmp_path / "idx.idx"),
        transcripts=str(tmp_path / "t2g.txt"),
        output=str(tmp_path / "out"),
        strand="forward",
    )
    base.update(extra)
    return argparse.Namespace(**base)


@pytest.fixture(autouse=True)
def _fixed_version(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(rs, "__version__", "0.0.test")


def test_the_omit_at_default_machinery_is_gone() -> None:
    for name in ("_OMIT_AT_DEFAULT", "_OMIT_AT_LEGACY", "_omittable_values", "_omit_at_default"):
        assert not hasattr(rs, name), name


def test_every_defaulted_option_is_recorded_at_its_resolved_value(tmp_path: Path) -> None:
    options = rs.build_run_manifest(_args(tmp_path))["options"]
    for key, value in DEFAULTS.items():
        if key in _run_defaults():
            assert options[key] == value, key


@pytest.mark.parametrize(
    "key, value",
    [
        ("programme_min_breadth", 3),
        ("programme_min_umi", 5.0),
        ("programme_latent_min_breadth", 2),
        ("emptydrops_niters", 5000),
        ("emptydrops_seed", 7),
        ("multimap_molecule_assignments", True),
        ("anello_align", True),
        ("anellovirus_gene_ids", False),
        ("discover", True),
        ("kraken2_db", "/db"),
        ("star_confirm", True),
        ("star_hostvirus_index", "/idx"),
        ("tchem", "off"),
        ("tchem_bin", "/bin/tchem"),
        ("tchem_whitelists", "/wl"),
        ("gene_programs", True),
        ("positive_control_gene", "g"),
        ("positive_control_molecules", 10),
        ("require_positive_control", True),
        ("programme_productive_min_molecules", 3),
        ("infected_min_molecules", 1.0),
    ],
)
def test_a_non_default_value_changes_the_fingerprint(tmp_path: Path, key, value) -> None:
    assert value != _run_defaults()[key]
    plain = rs.build_run_manifest(_args(tmp_path))
    changed = rs.build_run_manifest(_args(tmp_path, **{key: value}))
    assert changed["run_fingerprint"] != plain["run_fingerprint"]
    assert changed["options"][key] == value


def test_a_run_from_another_version_never_matches(tmp_path: Path, monkeypatch) -> None:
    old = rs.build_run_manifest(_args(tmp_path))
    monkeypatch.setattr(rs, "__version__", "0.0.next")
    new = rs.build_run_manifest(_args(tmp_path))
    assert old["viralscan_version"] != new["viralscan_version"]
    assert not rs._fingerprint_matches(old, new)


def test_a_legacy_hash_is_not_accepted(tmp_path: Path) -> None:
    """The pre-3.0.0rc2 compatibility branches are gone: only an exact match resumes."""
    new = rs.build_run_manifest(_args(tmp_path))
    old = json.loads(json.dumps(new))
    old["run_fingerprint"] = "0" * 64
    assert not rs._fingerprint_matches(old, new)


def test_refused_resume_names_the_option_that_differs(tmp_path: Path) -> None:
    first = rs.build_run_manifest(_args(tmp_path))
    out = tmp_path / "out"
    out.mkdir()
    (out / "run_manifest.json").write_text(json.dumps(first))
    args = _args(tmp_path, resume=True, emptydrops_seed=7)
    with pytest.raises(RunSafetyError, match="--emptydrops-seed"):
        prepare_output_directory(
            out, rs.build_run_manifest(args), resume=True, overwrite=False, yes=False
        )


def test_resume_with_the_same_options_and_version_is_accepted(tmp_path: Path) -> None:
    first = rs.build_run_manifest(_args(tmp_path))
    out = tmp_path / "out"
    out.mkdir()
    (out / "run_manifest.json").write_text(json.dumps(first))
    again = rs.build_run_manifest(_args(tmp_path, resume=True))
    assert prepare_output_directory(out, again, resume=True, overwrite=False, yes=False) == "resume"

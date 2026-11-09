"""FP-BASELINE-01: options at their default stay out of the run fingerprint.

Manifests written before these options existed, and manifests that stored them at the default,
must both resume against an invocation that leaves them alone. A changed value must refuse.
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
        # A 3' kit with a stated strand; the strand check is not under test here.
        strand="forward",
    )
    base.update(extra)
    return argparse.Namespace(**base)


@pytest.fixture(autouse=True)
def _fixed_version(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(rs, "__version__", "0.0.test")


def test_every_listed_option_is_omitted_at_its_parser_default(tmp_path: Path) -> None:
    options = rs.build_run_manifest(_args(tmp_path))["options"]
    defaults = _run_defaults()
    for key in rs._OMIT_AT_DEFAULT:
        assert key in defaults, f"{key} is not an option of `viralscan run`"
        assert key not in options, key
    # The legacy value is not a default: 2.0 is hashed, 1.0 is not.
    assert options["infected_min_molecules"] == DEFAULTS["infected_min_molecules"] == 2.0


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
    ],
)
def test_a_non_default_value_changes_the_fingerprint(tmp_path: Path, key, value) -> None:
    assert value != _run_defaults()[key]
    plain = rs.build_run_manifest(_args(tmp_path))
    changed = rs.build_run_manifest(_args(tmp_path, **{key: value}))
    assert changed["run_fingerprint"] != plain["run_fingerprint"]
    assert changed["options"][key] == value


def test_infected_min_molecules_one_is_omitted_and_two_is_hashed(tmp_path: Path) -> None:
    one = rs.build_run_manifest(_args(tmp_path, infected_min_molecules=1.0))
    two = rs.build_run_manifest(_args(tmp_path, infected_min_molecules=2.0))
    assert "infected_min_molecules" not in one["options"]
    assert two["options"]["infected_min_molecules"] == 2.0
    assert one["run_fingerprint"] != two["run_fingerprint"]


def _old_manifest(tmp_path: Path, *, stored_defaults: bool) -> dict:
    """What an earlier release wrote: hashed over its own options, with or without the keys."""
    new = rs.build_run_manifest(_args(tmp_path, infected_min_molecules=1.0))
    old = json.loads(json.dumps(new))
    if stored_defaults:
        for key, value in _run_defaults().items():
            if key in rs._OMIT_AT_DEFAULT:
                old["options"][key] = value
    old["run_fingerprint"] = rs._fingerprint(old)
    return old


@pytest.mark.parametrize("stored_defaults", [False, True])
def test_old_manifest_resumes_with_infected_min_molecules_one(
    tmp_path: Path, stored_defaults: bool
) -> None:
    old = _old_manifest(tmp_path, stored_defaults=stored_defaults)
    new = rs.build_run_manifest(_args(tmp_path, infected_min_molecules=1.0))
    assert rs._fingerprint_matches(old, new)
    # End to end: the output directory check accepts --resume.
    out = tmp_path / "out"
    out.mkdir()
    (out / "run_manifest.json").write_text(json.dumps(old))
    args = _args(tmp_path, infected_min_molecules=1.0, resume=True)
    assert (
        prepare_output_directory(
            out, rs.build_run_manifest(args), resume=True, overwrite=False, yes=False
        )
        == "resume"
    )


def test_old_manifest_refuses_under_the_new_default_of_two(tmp_path: Path) -> None:
    old = _old_manifest(tmp_path, stored_defaults=True)
    new = rs.build_run_manifest(_args(tmp_path))
    assert not rs._fingerprint_matches(old, new)


def test_a_stored_non_default_value_still_refuses(tmp_path: Path) -> None:
    old = _old_manifest(tmp_path, stored_defaults=True)
    old["options"]["emptydrops_seed"] = 7
    old["run_fingerprint"] = rs._fingerprint(old)
    new = rs.build_run_manifest(_args(tmp_path, infected_min_molecules=1.0))
    assert not rs._fingerprint_matches(old, new)


def test_a_manifest_whose_hash_does_not_reproduce_is_not_trusted(tmp_path: Path) -> None:
    old = _old_manifest(tmp_path, stored_defaults=True)
    old["run_fingerprint"] = "0" * 64
    new = rs.build_run_manifest(_args(tmp_path, infected_min_molecules=1.0))
    assert not rs._fingerprint_matches(old, new)


def test_refused_resume_names_the_option_that_differs(tmp_path: Path) -> None:
    old = _old_manifest(tmp_path, stored_defaults=False)
    out = tmp_path / "out"
    out.mkdir()
    (out / "run_manifest.json").write_text(json.dumps(old))
    args = _args(tmp_path, resume=True, infected_min_molecules=1.0, emptydrops_seed=7)
    with pytest.raises(RunSafetyError, match="--emptydrops-seed"):
        prepare_output_directory(
            out, rs.build_run_manifest(args), resume=True, overwrite=False, yes=False
        )

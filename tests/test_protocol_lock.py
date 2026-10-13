"""LOCK-02 software: schema and validator for the discovery study's protocol-v1.lock.json.

Every number below is a TEST_ONLY fixture value; none is a recommended or calibrated threshold.
"""

import hashlib
import json
import re
from pathlib import Path

import pytest

from scripts import validate_protocol_lock as vpl

# The keys scripts/adjudicate_candidates.py (ADJ-01) reads, dotted = nested objects.
ADJUDICATOR_KEYS = (
    "distributed_route.min_unique_molecules",
    "distributed_route.min_positive_cells",
    "distributed_route.min_independent_loci",
    "max_sibling_allocated_fraction",
    "min_virus_best_fraction",
    "min_complex_body_fraction",
    "l3_min_programme_cells",
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(tmp_path: Path) -> dict:
    """A complete TEST_ONLY lock whose cited files exist under tmp_path."""
    files = {}
    for name in (
        "gtf_manifest",
        "panel_manifest",
        "index",
        "t2g",
        "dlist",
        "human",
        "protocol",
        "holdout",
        "mdb",
        "lib1",
    ):
        path = tmp_path / f"{name}.txt"
        path.write_text(f"TEST_ONLY {name}\n")
        files[name] = path
    doc = vpl.template()
    doc.update(
        test_only=True,
        protocol_version="TEST_ONLY",
        created_at="2000-01-01T00:00:00Z",
        package_validation_protocol_sha256=_sha(files["protocol"]),
        package_validation_protocol_path="protocol.txt",
        holdout_manifest_sha256=_sha(files["holdout"]),
        holdout_manifest_path="holdout.txt",
        calibration_inputs=[
            {"library_id": "TEST_ONLY_lib1", "sha256": _sha(files["lib1"]), "path": "lib1.txt"}
        ],
    )
    doc["code"] = {k: "TEST_ONLY" for k in doc["code"]} | {"repo_sha": "0" * 40}
    doc["panel"] = {
        "gtf_corpus_manifest_sha256": _sha(files["gtf_manifest"]),
        "gtf_corpus_manifest_path": "gtf_manifest.txt",
        "panel_manifest_sha256": _sha(files["panel_manifest"]),
        "index_sha256": _sha(files["index"]),
        "index_path": "index.txt",
        "t2g_sha256": _sha(files["t2g"]),
        "dlist": {"identity": "TEST_ONLY", "sha256": _sha(files["dlist"]), "path": "dlist.txt"},
        "human_reference": {"identity": "TEST_ONLY", "sha256": _sha(files["human"])},
        "ncbi_snapshot_date": "2000-01-01",
        "taxonomy_snapshot_date": "2000-01-01",
    }
    doc["strand_policy"]["calibration_rule"]["min_positive_viral_molecule_fraction"] = 0.5
    doc["strand_policy"]["reverse_chemistries"] = []
    doc["cell_calling"].update(emptydrops_seed=1, emptydrops_niters=2)
    doc["adjudication"] = {
        "test_only": True,
        "distributed_route": {
            "min_unique_molecules": 1,
            "min_positive_cells": 2,
            "min_independent_loci": 3,
        },
        "max_sibling_allocated_fraction": 0.5,
        "min_virus_best_fraction": 0.5,
        "min_complex_body_fraction": 0.5,
        "l3_min_programme_cells": 4,
    }
    doc["false_discovery"] = {
        "expected_false_l2_per_library": {"target": 1, "interval": [0, 2]},
        "decoy": {"enabled": True, "seeds": [1], "strata": ["TEST_ONLY_bin"]},
    }
    doc["sensitivity"] = {
        "minimum_detectable_burden_reference": "TEST_ONLY",
        "minimum_detectable_burden_sha256": _sha(files["mdb"]),
        "minimum_detectable_burden_path": "mdb.txt",
    }
    return doc


def _write(tmp_path: Path, doc: dict) -> Path:
    path = tmp_path / "protocol-v1.lock.json"
    path.write_text(json.dumps(doc))
    return path


def _run(tmp_path: Path, doc: dict, *flags: str) -> int:
    return vpl.main([str(_write(tmp_path, doc)), *flags])


def test_valid_fixture_passes(tmp_path, capsys) -> None:
    assert _run(tmp_path, _fixture(tmp_path), "--allow-test-only") == vpl.EXIT_VALID
    assert capsys.readouterr().out.strip() == "valid"


def test_non_test_only_lock_passes_without_flag(tmp_path) -> None:
    doc = _fixture(tmp_path)
    del doc["test_only"], doc["adjudication"]["test_only"]
    assert _run(tmp_path, doc) == vpl.EXIT_VALID


def test_test_only_refused_without_flag(tmp_path, capsys) -> None:
    assert _run(tmp_path, _fixture(tmp_path)) == vpl.EXIT_INVALID
    assert "test_only" in capsys.readouterr().err


def test_test_only_spellings_must_agree(tmp_path) -> None:
    doc = _fixture(tmp_path)
    del doc["adjudication"]["test_only"]
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_INVALID


def test_missing_key_fails(tmp_path, capsys) -> None:
    doc = _fixture(tmp_path)
    del doc["adjudication"]["l3_min_programme_cells"]
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_INVALID
    assert "l3_min_programme_cells" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("adjudication", "min_virus_best_fraction", "0.5"),  # string, not number
        ("adjudication", "min_virus_best_fraction", 1.5),  # fraction above 1
        ("adjudication", "l3_min_programme_cells", -1),
        ("cell_calling", "emptydrops_seed", 1.5),
        ("code", "repo_sha", "abc"),
    ],
)
def test_wrong_type_or_range_fails(tmp_path, section, key, value) -> None:
    doc = _fixture(tmp_path)
    doc[section][key] = value
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_INVALID


def test_unknown_adjudication_key_fails(tmp_path) -> None:
    doc = _fixture(tmp_path)
    doc["adjudication"]["min_made_up"] = 1
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_INVALID


@pytest.mark.parametrize("where", ["top", "nested", "list"])
def test_sha_mismatch_fails(tmp_path, capsys, where) -> None:
    doc = _fixture(tmp_path)
    if where == "top":
        doc["package_validation_protocol_sha256"] = "1" * 64
    elif where == "nested":
        doc["panel"]["dlist"]["sha256"] = "1" * 64
    else:
        doc["calibration_inputs"][0]["sha256"] = "1" * 64
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_INVALID
    assert "mismatch" in capsys.readouterr().err


def test_sha_check_follows_file_content_and_missing_file(tmp_path) -> None:
    doc = _fixture(tmp_path)
    (tmp_path / "index.txt").write_text("changed\n")
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_INVALID
    doc = _fixture(tmp_path)
    (tmp_path / "t2g.txt").unlink()
    doc["panel"]["t2g_path"] = "t2g.txt"
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_INVALID


def test_root_overrides_lock_directory(tmp_path) -> None:
    doc = _fixture(tmp_path)
    lock_dir = tmp_path / "elsewhere"
    lock_dir.mkdir()
    lock = lock_dir / "protocol-v1.lock.json"
    lock.write_text(json.dumps(doc))
    assert (
        vpl.main([str(lock), "--allow-test-only"]) == vpl.EXIT_INVALID
    )  # files are not next to the lock
    assert vpl.main([str(lock), "--allow-test-only", "--root", str(tmp_path)]) == vpl.EXIT_VALID


def test_paths_are_optional(tmp_path) -> None:
    doc = _fixture(tmp_path)
    for key in [k for k in doc["panel"] if k.endswith("_path")]:
        del doc["panel"][key]
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_VALID


def test_unreadable_or_duplicate_key_json_is_usage_error(tmp_path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text('{"status": "exploratory", "status": "final"}')
    assert vpl.main([str(bad)]) == vpl.EXIT_USAGE
    assert vpl.main([str(tmp_path / "absent.json")]) == vpl.EXIT_USAGE


def test_template_validates_as_incomplete_not_valid(tmp_path, capsys) -> None:
    assert vpl.main(["--template"]) == vpl.EXIT_VALID
    template = json.loads(capsys.readouterr().out)
    assert vpl.json_schema_errors(template, vpl.DEFAULT_SCHEMA) == []
    assert vpl.main([str(_write(tmp_path, template))]) == vpl.EXIT_INCOMPLETE
    out = capsys.readouterr().out
    assert out.startswith("incomplete") and "adjudication.l3_min_programme_cells" in out
    assert "calibration_inputs" in out


def test_template_carries_no_numbers() -> None:
    def numbers(node):
        if isinstance(node, dict):
            return [n for v in node.values() for n in numbers(v)]
        if isinstance(node, list):
            return [n for v in node for n in numbers(v)]
        return [node] if isinstance(node, (int, float)) and not isinstance(node, bool) else []

    assert numbers(vpl.template()) == []


def test_disabled_decoy_is_not_incomplete(tmp_path) -> None:
    doc = _fixture(tmp_path)
    doc["false_discovery"]["decoy"] = {"enabled": False, "seeds": None, "strata": None}
    assert _run(tmp_path, doc, "--allow-test-only") == vpl.EXIT_VALID


def test_adjudicator_reads_the_adjudication_block_directly(tmp_path) -> None:
    """Same traversal as adjudicate_candidates.load_lock: doc.get('adjudication', doc), dotted = nested."""
    doc = _fixture(tmp_path)
    for key in ADJUDICATOR_KEYS:
        node = doc["adjudication"]
        for part in key.split("."):
            node = node[part]
        assert isinstance(node, (int, float)) and not isinstance(node, bool) and node >= 0
    assert set(doc["adjudication"]) - {"test_only"} == {k.split(".")[0] for k in ADJUDICATOR_KEYS}
    assert (vpl.REPO_ROOT / "schemas/v3/protocol_lock.schema.json").is_file()


def test_schema_does_not_extend_package_validation_protocol() -> None:
    schema = json.loads(vpl.DEFAULT_SCHEMA.read_text())
    refs = re.findall(r'"\$ref": "([^"]+)"', json.dumps(schema))
    assert refs and not [r for r in refs if "validation_protocol" in r]
    assert "package_validation_protocol_sha256" in schema["required"]

"""Release preparation checks exercise policy without publishing or running CI."""

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

from scripts.check_actions_pinned import unpinned_actions
from scripts.check_versions import version_errors
from scripts.release_sif_definition import definition
from scripts.write_explicit_lock import explicit_lines

ROOT = Path(__file__).resolve().parents[1]


def test_action_pin_check_rejects_mutable_tags_and_preserves_local_actions():
    assert unpinned_actions("- uses: actions/checkout@v4\n  uses: owner/repo@main\n") == [
        "actions/checkout@v4",
        "owner/repo@main",
    ]
    assert unpinned_actions(f"- uses: actions/checkout@{'a' * 40}\n- uses: ./local\n") == []
    assert unpinned_actions(f"- uses: docker://image@sha256:{'b' * 64}\n") == []


@pytest.mark.parametrize(
    "tag,latest", [("v3.0.0.dev1", False), ("v3.0.0rc1", False), ("v3.0.0", True)]
)
def test_actual_release_tag_step_preserves_latest_for_prereleases(tmp_path, tag, latest):
    workflow = yaml.safe_load((ROOT / ".github/workflows/release.yml").read_text())
    step = next(s for s in workflow["jobs"]["container"]["steps"] if s.get("id") == "meta")
    output = tmp_path / "output"
    subprocess.run(
        ["bash", "-eu", "-c", step["run"]],
        check=True,
        env={**os.environ, "GITHUB_REF": f"refs/tags/{tag}", "GITHUB_OUTPUT": str(output)},
    )
    fields = dict(line.split("=", 1) for line in output.read_text().splitlines())
    assert fields["version"] == tag
    assert (":latest" in fields["tags"]) is latest


def test_versions_agree_and_recipe_drift_is_detected(tmp_path):
    assert version_errors(ROOT) == []
    for rel in [
        "src/viralscan/__init__.py",
        "Dockerfile",
        "Singularity.def",
        "conda-recipe/meta.yaml",
        "CITATION.cff",
    ]:
        target = tmp_path / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, target)
    recipe = tmp_path / "conda-recipe/meta.yaml"
    recipe.write_text(
        re.sub(r'\{% set version = "[^"]+" %\}', '{% set version = "0.0.0" %}', recipe.read_text())
    )
    assert any("conda-recipe/meta.yaml" in e for e in version_errors(tmp_path))


def test_sif_definition_derives_from_digest_without_a_second_solve():
    ref = f"ghcr.io/mdmanurung/viralscan@sha256:{'c' * 64}"
    text = definition(ref)
    assert f"From: {ref}" in text and "%post" not in text
    with pytest.raises(ValueError, match="mutable tags"):
        definition("ghcr.io/mdmanurung/viralscan:latest")


def test_write_explicit_lock_sorted_urls_with_md5_and_no_local_paths(tmp_path):
    meta = tmp_path / "conda-meta"
    meta.mkdir()
    for name, md5 in (("zlib-1.0-0", "b" * 32), ("abc-2.0-1", "a" * 32)):
        url = f"https://conda.anaconda.org/conda-forge/linux-64/{name}.conda"
        rec = {"url": url, "md5": md5, "package_tarball_full_path": "/home/u/pkgs/x"}
        (meta / f"{name}.json").write_text(json.dumps(rec))
    (meta / "history").write_text("ignored")
    assert explicit_lines(tmp_path) == [
        f"https://conda.anaconda.org/conda-forge/linux-64/abc-2.0-1.conda#{'a' * 32}",
        f"https://conda.anaconda.org/conda-forge/linux-64/zlib-1.0-0.conda#{'b' * 32}",
    ]
    (meta / "bad-1-0.json").write_text(
        json.dumps({"url": "file:///home/u/x.conda", "md5": "c" * 32})
    )
    with pytest.raises(SystemExit):
        explicit_lines(tmp_path)

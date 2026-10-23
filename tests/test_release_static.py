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
from scripts.check_release_ci import successful_ci_run
from scripts.check_versions import version_errors
from scripts.release_sif_definition import definition
from scripts.write_explicit_lock import explicit_lines

ROOT = Path(__file__).resolve().parents[1]


def _ci_run(**changes):
    return {
        "id": 1,
        "head_sha": "a" * 40,
        "path": ".github/workflows/ci.yml",
        "event": "push",
        "head_branch": "main",
        "head_repository": {"full_name": "owner/repository"},
        "run_number": 1,
        "run_attempt": 1,
        "status": "completed",
        "conclusion": "success",
        **changes,
    }


def test_release_ci_accepts_only_complete_success_for_same_commit():
    run = _ci_run()
    assert successful_ci_run([run], "a" * 40, "owner/repository") == run


@pytest.mark.parametrize(
    "change",
    [
        {"head_sha": "b" * 40},
        {"path": ".github/workflows/release.yml"},
        {"status": "in_progress", "conclusion": None},
        {"status": "queued", "conclusion": None},
        {"conclusion": "failure"},
        {"conclusion": "stale"},
        {"conclusion": "cancelled"},
        {"event": "pull_request"},
        {"head_branch": "feature"},
        {"head_repository": {"full_name": "fork/repository"}},
    ],
)
def test_release_ci_rejects_unqualified_run(change):
    with pytest.raises(ValueError):
        successful_ci_run([_ci_run(**change)], "a" * 40, "owner/repository")


def test_release_ci_does_not_reuse_old_success_when_latest_run_or_attempt_failed():
    for later in (
        _ci_run(run_number=2, conclusion="failure"),
        _ci_run(run_attempt=2, status="in_progress", conclusion=None),
    ):
        with pytest.raises(ValueError, match="Latest CI run"):
            successful_ci_run([_ci_run(), later], "a" * 40, "owner/repository")


def test_release_publications_require_ci_gate_and_tested_artifacts():
    jobs = yaml.safe_load((ROOT / ".github/workflows/release.yml").read_text())["jobs"]
    assert jobs["build"]["needs"] == "ci-gate"
    for name in ("container", "publish"):
        assert "ci-gate" in jobs[name]["needs"]
        assert "build" in jobs[name]["needs"]
        assert jobs[name]["environment"] == "${{ needs.build.outputs.publish_environment }}"
        assert any(step.get("with", {}).get("name") == "dist" for step in jobs[name]["steps"])
        assert any("sha256sum -c" in step.get("run", "") for step in jobs[name]["steps"])
    builds = [
        s for s in jobs["container"]["steps"] if "docker/build-push-action@" in s.get("uses", "")
    ]
    assert len(builds) == 1 and builds[0]["with"]["push"] is False
    assert builds[0]["with"]["build-contexts"] == "wheel=./dist"
    scans = [s for s in jobs["container"]["steps"] if "anchore/sbom-action@" in s.get("uses", "")]
    assert scans[0]["with"]["upload-release-assets"] is False
    smoke = next(s for s in jobs["build"]["steps"] if "Clean wheel" in s.get("name", ""))["run"]
    assert "--no-deps" not in smoke and "--system-site-packages" not in smoke
    assert "*.whl" in smoke and "*.tar.gz" in smoke and "release_package_smoke.py" in smoke


def test_container_uses_lock_and_wheel_without_an_independent_sif_solve():
    dockerfile = (ROOT / "Dockerfile").read_text()
    assert "--file /tmp/environment.explicit.txt" in dockerfile
    assert "from=wheel" in dockerfile and "pip install --no-deps /tmp/dist/*.whl" in dockerfile
    assert "COPY src" not in dockerfile and "pip install --no-deps ." not in dockerfile
    assert "%post" not in (ROOT / "Singularity.def").read_text()


def test_release_smoke_fixture_passes_production_validation(tmp_path, monkeypatch):
    from scripts.release_package_smoke import smoke

    monkeypatch.setenv("PYTHONPATH", str(ROOT / "src"))
    smoke(tmp_path)
    assert (tmp_path / "run" / "run_complete.json").is_file()


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
    assert version_errors(ROOT, "v3.0.0rc1") == []
    assert any("release tag" in e for e in version_errors(ROOT, "v0.0.0"))
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


def _lock_audit_workflow():
    return yaml.safe_load((ROOT / ".github/workflows/lock-audit.yml").read_text())


def test_weekly_lock_audit_is_independent_and_preserves_evidence():
    workflow = _lock_audit_workflow()
    triggers = workflow.get("on", workflow.get(True))
    assert set(triggers) == {"schedule", "workflow_dispatch"}
    assert triggers["schedule"] == [{"cron": "23 7 * * 1"}]
    assert workflow["permissions"] == {"contents": "read"}
    job = workflow["jobs"]["audit"]
    assert job["continue-on-error"] is True and "needs" not in job
    assert job["strategy"]["fail-fast"] is False
    assert job["strategy"]["matrix"]["lock"] == ["uv", "conda"]
    steps = job["steps"]
    assert next(s for s in steps if s.get("id") == "audit")["continue-on-error"] is True
    assert steps[-2]["if"] == "always()" and steps[-1]["if"] == "always()"
    assert "upload-artifact@" in steps[-1]["uses"]
    assert unpinned_actions((ROOT / ".github/workflows/lock-audit.yml").read_text()) == []


@pytest.mark.parametrize("lock", ["uv", "conda"])
@pytest.mark.parametrize("code", [0, 1, 2])
def test_weekly_audit_keeps_exit_code_and_log(tmp_path, lock, code):
    job = _lock_audit_workflow()["jobs"]["audit"]
    step = next(s for s in job["steps"] if s.get("id") == "audit")
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    (scripts / "audit_locks.py").write_text(
        f"def audit_uv_lock(path):\n    print('uv fixture'); return {code}\n"
        f"def audit_conda_lock(path):\n    print('conda fixture'); return {code}\n"
    )
    for filename in ("uv.lock", "environment.explicit.txt"):
        (tmp_path / filename).write_text("fixture\n")
    receipts = tmp_path / "receipts"
    result = subprocess.run(
        ["bash", "-eu", "-o", "pipefail", "-c", step["run"]],
        cwd=tmp_path,
        env={
            **os.environ,
            "LOCK_AUDIT_DIR": str(receipts),
            "LOCK_KIND": lock,
            "PYTHONPATH": str(tmp_path),
        },
        capture_output=True,
        text=True,
    )
    assert result.returncode == code
    assert (receipts / "exit_code.txt").read_text().strip() == str(code)
    assert (receipts / "audit.log").read_text().strip() == f"{lock} fixture"
    assert (receipts / "locks.sha256").is_file()


def test_weekly_skipped_audit_is_recorded_as_skipped(tmp_path):
    step = _lock_audit_workflow()["jobs"]["audit"]["steps"][-2]
    receipts, summary = tmp_path / "receipts", tmp_path / "summary"
    subprocess.run(
        ["bash", "-eu", "-o", "pipefail", "-c", step["run"]],
        check=True,
        env={
            **os.environ,
            "LOCK_AUDIT_DIR": str(receipts),
            "CHECKOUT_OUTCOME": "failure",
            "PYTHON_OUTCOME": "skipped",
            "AUDIT_OUTCOME": "skipped",
            "AUDITED_SHA": "a" * 40,
            "LOCK_KIND": "conda",
            "GITHUB_STEP_SUMMARY": str(summary),
        },
    )
    assert "audit=skipped" in (receipts / "status.txt").read_text()
    assert "exit_code=not-recorded" in summary.read_text()
    assert "Skipped or failed is not a clean scan." in summary.read_text()

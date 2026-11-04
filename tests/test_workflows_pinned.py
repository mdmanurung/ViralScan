"""Every workflow action is a commit pin with a version comment; releases need the CI gate (REL-12/13)."""

import re
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = sorted((ROOT / ".github/workflows").glob("*.yml"))
PIN = re.compile(r"[\w.-]+/[\w.-]+(/[\w./-]+)?@[0-9a-f]{40}")
DOCKER_PIN = re.compile(r"docker://.+@sha256:[0-9a-f]{64}")
USES_LINE = re.compile(r"^\s*(?:-\s*)?uses:\s*(\S+)(\s+#\s*v\d+(\.\d+)*\S*)?\s*$")


def _uses(workflow):
    jobs = yaml.safe_load(workflow.read_text())["jobs"]
    for job in jobs.values():
        if "uses" in job:  # reusable workflow call
            yield job["uses"]
        for step in job.get("steps", []):
            if "uses" in step:
                yield step["uses"]


def test_workflows_found():
    assert {p.name for p in WORKFLOWS} >= {"ci.yml", "release.yml"}


@pytest.mark.parametrize("workflow", WORKFLOWS, ids=lambda p: p.name)
def test_every_uses_is_a_full_commit_pin(workflow):
    bad = [
        u
        for u in _uses(workflow)
        if not (u.startswith("./") or PIN.fullmatch(u) or DOCKER_PIN.fullmatch(u))
    ]
    assert bad == []


@pytest.mark.parametrize("workflow", WORKFLOWS, ids=lambda p: p.name)
def test_every_remote_pin_carries_a_version_comment(workflow):
    bad = []
    for line in workflow.read_text().splitlines():
        if re.match(r"^\s*(?:-\s*)?uses:", line):
            m = USES_LINE.match(line)
            remote = m and not m.group(1).startswith("./")
            if not m or (remote and not m.group(2)):
                bad.append(line.strip())
    assert bad == []


def _needs(job):
    need = job.get("needs", [])
    return {need} if isinstance(need, str) else set(need)


def test_every_release_job_depends_on_the_ci_gate():
    jobs = yaml.safe_load((ROOT / ".github/workflows/release.yml").read_text())["jobs"]

    def reaches_gate(name, seen=()):
        needs = _needs(jobs[name])
        return "ci-gate" in needs or any(
            reaches_gate(n, (*seen, name)) for n in needs if n not in seen
        )

    for name in jobs:
        if name != "ci-gate":
            assert reaches_gate(name), f"{name} does not depend on ci-gate"
    for name in ("container", "publish"):
        assert "ci-gate" in _needs(jobs[name]), f"{name} must need ci-gate directly"
    gate = " ".join(s.get("run", "") for s in jobs["ci-gate"]["steps"])
    assert "check_release_ci.py" in gate

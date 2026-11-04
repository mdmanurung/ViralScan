"""Offline replay checks for the committed conda locks (PLAN REL-03, DEF-08, REL-11).

A real replay (``conda create --file``) needs the network and an install; the release
``parity`` job and the CI ``environment-lock`` job do that. These tests catch the lock
mistakes that would fail the replay, and keep the reviewed-exception scaffold honest.
"""

import collections
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

from scripts.audit_locks import apply_exceptions

ROOT = Path(__file__).resolve().parents[1]
LOCKS = ["environment.explicit.txt", "environment.full.explicit.txt"]
ROW = re.compile(
    r"https://conda\.anaconda\.org/(conda-forge|bioconda)/(linux-64|noarch)/"
    r"(.+)-([^-]+)-([^-]+)\.(conda|tar\.bz2)#[0-9a-f]{32}"
)


def _rows(name):
    lines = (ROOT / name).read_text().splitlines()
    assert "@EXPLICIT" in lines
    return [line for line in lines if line and not line.startswith(("#", "@"))]


@pytest.mark.parametrize("name", LOCKS)
def test_lock_rows_are_replayable(name):
    rows = _rows(name)
    assert all(ROW.fullmatch(r) for r in rows), [r for r in rows if not ROW.fullmatch(r)]
    names = collections.Counter(ROW.fullmatch(r).group(3) for r in rows)
    assert [n for n, c in names.items() if c > 1] == []
    assert rows == sorted(rows)


@pytest.mark.parametrize("name", LOCKS)
def test_lock_carries_the_pip_check_corrections(name):
    # REL-03-PIPCHECK: interface-common 1.23.1 needs packaging>=26.1, snakemake 9.23.1 needs <26.
    text = "\n".join(_rows(name))
    assert "/snakemake-interface-common-1.23.0-" in text
    assert "/snakemake-interface-common-1.23.1-" not in text
    assert "/cython-" in text


def test_full_lock_adds_the_dev_tools_to_the_runtime_lock():
    def names(n):
        return {ROW.fullmatch(r).group(3) for r in _rows(n)}

    assert {"pytest", "snakemake", "kb-python"} <= names("environment.full.explicit.txt")
    assert "pytest" not in names("environment.explicit.txt")  # runtime-only image lock


@pytest.mark.network
@pytest.mark.skipif(shutil.which("micromamba") is None, reason="micromamba not on PATH")
@pytest.mark.parametrize("name", LOCKS)
def test_lock_dry_run_replay(name, tmp_path):
    subprocess.run(
        [
            "micromamba", "create", "--dry-run", "--yes",
            "--prefix", str(tmp_path / "env"), "--file", str(ROOT / name),
        ],
        check=True,
        capture_output=True,
        timeout=900,
    )  # fmt: skip


def test_exception_scaffold_is_pending_and_silences_nothing(capsys):
    data = yaml.safe_load((ROOT / "security/reviewed_exceptions.yaml").read_text())
    entries = data["exceptions"]
    assert entries and all(e["status"] == "pending maintainer review" for e in entries)
    assert all(e["advisories"] and e["reason"] and e["scope"] for e in entries)
    e = entries[0]
    hits = {(e["package"], e["version"]): list(e["advisories"])}
    assert apply_exceptions(ROOT / e["lock"], hits) == hits
    assert "pending maintainer review" in capsys.readouterr().out


def test_reviewed_entry_silences_only_its_listed_ids(tmp_path):
    entry = {
        "lock": "uv.lock",
        "package": "x",
        "version": "1",
        "advisories": ["A"],
        "status": "reviewed",
    }
    f = tmp_path / "ex.yaml"
    f.write_text(yaml.safe_dump({"exceptions": [entry]}))
    kept = apply_exceptions(Path("uv.lock"), {("x", "1"): ["A", "B"], ("y", "1"): ["C"]}, f)
    assert kept == {("x", "1"): ["B"], ("y", "1"): ["C"]}

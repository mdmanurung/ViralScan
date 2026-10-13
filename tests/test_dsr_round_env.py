"""DSR-16: the round scripts default to round 1's paths and take ROUND_DIR/CODE_DIR/REF_BUILD overrides."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
M = "/exports/archive/hg-funcgenom-research/mdmanurung"
ROUND1 = {
    "R": f"{M}/viralscan_work/dsr_round1",
    "CODE": f"{M}/vs_pinned/bbf1821",
    "B": f"{M}/viral_ref_cat42d/build",
    "REF": f"{M}/viral_ref_cat42d/build",
}
NAMES = {
    "slurm_dsr_rerun_array.sh": ["R", "CODE", "B", "VIRUS_ONLY"],
    "slurm_dsr_chem_array.sh": ["R", "CODE", "REF"],
    "slurm_dsr_collect.sh": ["R"],
    "slurm_dsr_twostep_redetect.sh": ["R"],
}


def _eval(script: str, env: dict[str, str]) -> dict[str, str]:
    text = (SCRIPTS / script).read_text()
    block = re.search(r"# >>> dsr-env.*?\n(.*?)# <<< dsr-env", text, re.S).group(1)
    probe = block + "".join(f'echo "{n}=${n}"\n' for n in NAMES[script])
    out = subprocess.run(
        ["bash", "-c", probe], env=env, capture_output=True, text=True, check=True
    ).stdout
    return dict(line.split("=", 1) for line in out.splitlines())


@pytest.mark.parametrize("script", sorted(NAMES))
def test_defaults_reproduce_round1(script):
    got = _eval(script, {"PATH": "/usr/bin:/bin"})
    for name, value in got.items():
        if name in ROUND1:
            assert value == ROUND1[name], (script, name)
    if "VIRUS_ONLY" in got:
        assert got["VIRUS_ONLY"] == f"{M}/viral_ref_cat42d/build_virus_only"


@pytest.mark.parametrize("script", sorted(NAMES))
def test_overrides_take_effect(script):
    env = {
        "PATH": "/usr/bin:/bin",
        "ROUND_DIR": "/r2",
        "CODE_DIR": "/code/abc",
        "REF_BUILD": "/refs/viral_ref_new/build",
    }
    got = _eval(script, env)
    assert got["R"] == "/r2"
    if "CODE" in got:
        assert got["CODE"] == "/code/abc"
    for name in ("B", "REF"):
        if name in got:
            assert got[name] == "/refs/viral_ref_new/build"
    if "VIRUS_ONLY" in got:
        assert got["VIRUS_ONLY"] == "/refs/viral_ref_new/build_virus_only"


@pytest.mark.parametrize("script", sorted(NAMES))
def test_bash_syntax_and_no_hardcoded_round1_after_header(script):
    path = SCRIPTS / script
    subprocess.run(["bash", "-n", str(path)], check=True)
    body = path.read_text().split("# <<< dsr-env", 1)[1]
    assert "dsr_round1" not in re.sub(r"(?m)^#.*$", "", body)  # only comments/#SBATCH may name it


def test_shellcheck_if_available():
    if shutil.which("shellcheck") is None:
        pytest.skip("shellcheck not installed")
    for script in NAMES:
        subprocess.run(["shellcheck", "-S", "error", str(SCRIPTS / script)], check=True)

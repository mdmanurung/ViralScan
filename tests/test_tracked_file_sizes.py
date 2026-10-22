import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "check_tracked_file_sizes.py"


def _run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def test_repo_has_no_oversized_tracked_research_files():
    r = _run()
    assert r.returncode == 0, r.stderr


def test_guard_flags_and_allowlists(tmp_path):
    def git(*a):
        subprocess.run(["git", "-C", str(tmp_path), *a], check=True, capture_output=True)

    git("init", "-q")
    (tmp_path / "analysis").mkdir()
    (tmp_path / "analysis" / "big.tsv").write_bytes(b"x" * 2000)
    git("add", "analysis/big.tsv")  # staged only, not committed
    assert _run("--repo", str(tmp_path), "--max-bytes", "1000").returncode == 1
    allow = tmp_path / "allow.txt"
    allow.write_text("# c\nanalysis/big.tsv\n")
    assert (
        _run("--repo", str(tmp_path), "--max-bytes", "1000", "--allowlist", str(allow)).returncode
        == 0
    )

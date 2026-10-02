import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import write_gtf_manifest as m  # noqa: E402


def _dir(tmp_path):
    d = tmp_path / "data"
    d.mkdir()
    (d / "b.gtf").write_text("b\n")
    (d / "a.gtf").write_text("a\n")
    return d


def test_write_and_check(tmp_path):
    d, out = _dir(tmp_path), tmp_path / "m.tsv"
    assert m.main(["--data-dir", str(d), "-o", str(out)]) == 0
    lines = out.read_text().splitlines()
    assert lines[0] == m.HEADER and [l.split("\t")[0] for l in lines[1:]] == ["a.gtf", "b.gtf"]
    assert lines[1].endswith("\t2\tunknown")
    assert m.main(["--data-dir", str(d), "--check", str(out)]) == 0
    (d / "a.gtf").write_text("changed\n")
    assert m.main(["--data-dir", str(d), "--check", str(out)]) == 1
    (d / "a.gtf").write_text("a\n")
    (d / "c.gtf").write_text("x\n")
    assert m.main(["--data-dir", str(d), "--check", str(out)]) == 1  # extra
    (d / "c.gtf").unlink()
    (d / "b.gtf").unlink()
    assert m.main(["--data-dir", str(d), "--check", str(out)]) == 1  # missing


def test_real_data(tmp_path):
    if not list(m.DEFAULT_DIR.glob("*.gtf")):
        pytest.skip("bundled GTFs absent")
    out = tmp_path / "m.tsv"
    assert m.main(["-o", str(out)]) == 0
    assert m.main(["--check", str(out)]) == 0

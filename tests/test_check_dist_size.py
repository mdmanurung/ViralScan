import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "check_dist_size", Path(__file__).parent.parent / "scripts" / "check_dist_size.py"
)
cds = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cds)


def test_under_limit_passes(tmp_path):
    (tmp_path / "a.whl").write_bytes(b"x" * 1000)
    (tmp_path / "a.tar.gz").write_bytes(b"x" * 1000)
    assert cds.main([str(tmp_path)]) == 0


def test_over_limit_fails(tmp_path, capsys):
    (tmp_path / "a.tar.gz").write_bytes(b"x" * 2_000_000)
    assert cds.main([str(tmp_path), "--max-mb", "1"]) == 1
    assert "a.tar.gz" in capsys.readouterr().err


def test_empty_dir_is_error(tmp_path):
    assert cds.main([str(tmp_path)]) == 2

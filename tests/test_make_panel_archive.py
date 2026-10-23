"""scripts/make_panel_archive.py (PLAN REF-11): same panel directory gives the same archive bytes."""

import subprocess
import sys
import zipfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "make_panel_archive.py"


def _build(panel: Path, out: Path) -> Path:
    subprocess.run(
        [sys.executable, str(SCRIPT), "--panel-dir", str(panel), "--out-dir", str(out),
         "--version", "0.0.0test", "--no-corpus-check"],
        check=True, capture_output=True,
    )  # fmt: skip
    return out / "viralscan-viral-panel-0.0.0test.zip"


def test_archive_is_deterministic_and_sorted(tmp_path):
    panel = tmp_path / "panel"
    panel.mkdir()
    for name in ("b_virus.gtf", "a_virus.gtf"):
        (panel / name).write_text(f'NC_1\tx\tgene\t1\t9\t.\t+\t.\tgene_id "{name}";\n')
    first = _build(panel, tmp_path / "out1")
    second = _build(panel, tmp_path / "out2")
    assert first.read_bytes() == second.read_bytes()
    names = zipfile.ZipFile(first).namelist()
    assert names == sorted(names) and {"a_virus.gtf", "b_virus.gtf", "PANEL_MANIFEST.json"} <= set(
        names
    )

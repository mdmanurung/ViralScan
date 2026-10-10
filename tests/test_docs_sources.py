"""Offline documentation checks; Sphinx and clean installs remain separate gates."""

import re
import shlex
from pathlib import Path

from viralscan.menu import build_parser

ROOT = Path(__file__).resolve().parents[1]


def test_public_toctree_targets_exist():
    for text in re.findall(r"```\{toctree\}\n(.*?)```", (ROOT / "docs/index.md").read_text(), re.S):
        for row in text.splitlines():
            target = row.strip()
            if target and not target.startswith(":"):
                assert any(
                    (ROOT / "docs" / (target + ext)).is_file() for ext in (".md", ".rst", ".ipynb")
                ), target


def test_quickstart_commands_parse_without_running_workflows():
    text = (ROOT / "docs/quickstart.md").read_text()
    checked = 0
    for block in re.findall(r"```bash\n(.*?)```", text, re.S):
        for command in block.replace("\\\n", " ").splitlines():
            argv = shlex.split(command, comments=True)
            if not argv or argv[0] != "viralscan":
                continue
            if argv[1:] == ["--help"]:
                continue
            build_parser().parse_args(argv[1:])
            checked += 1
    assert checked >= 4


def test_ci_notebooks_have_no_institutional_paths():
    import json

    notebooks = list((ROOT / "docs/vignettes").glob("*.ipynb"))
    assert len(notebooks) == 8
    for path in notebooks:
        doc = json.loads(path.read_text())
        for cell in doc["cells"]:
            if cell["cell_type"] == "code":
                code = "".join(cell["source"])
                assert not re.search(r"/(?:exports|home|lustre|gpfs)/", code), path

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


def test_public_notebook_cells_have_stable_ids_and_valid_python():
    import json

    for path in (ROOT / "docs/vignettes").glob("*.ipynb"):
        doc = json.loads(path.read_text())
        ids = [cell["id"] for cell in doc["cells"]]
        assert len(set(ids)) == len(ids), path
        assert all(re.fullmatch(r"[A-Za-z0-9_-]{1,64}", value) for value in ids), path
        assert doc["nbformat_minor"] >= 5, path
        for cell in doc["cells"]:
            if cell["cell_type"] == "code":
                compile("".join(cell["source"]), str(path), "exec")


def test_literal_notebook_cli_commands_parse_without_running_workflows():
    import ast
    import json

    checked = 0
    for path in (ROOT / "docs/vignettes").glob("*.ipynb"):
        for cell in json.loads(path.read_text())["cells"]:
            if cell["cell_type"] != "code":
                continue
            for node in ast.walk(ast.parse("".join(cell["source"]))):
                if not isinstance(node, ast.List):
                    continue
                if not all(isinstance(value, ast.Constant) for value in node.elts):
                    continue
                argv = ast.literal_eval(node)
                if not argv or argv[0] != "viralscan" or argv[1:] == ["--version"]:
                    continue
                args = build_parser().parse_args(argv[1:])
                if getattr(args, "_subcommand", None) == "build-ref":
                    assert args.host or args.reference_panel == "anellovirus", path
                checked += 1
    assert checked >= 3


def test_happy_path_commands_parse_and_cover_the_path():
    text = (ROOT / "docs/happy_path.md").read_text()
    seen = []
    for block in re.findall(r"```bash\n(.*?)```", text, re.S):
        for command in block.replace("\\\n", " ").splitlines():
            argv = shlex.split(command, comments=True)
            if not argv or argv[0] != "viralscan":
                continue
            build_parser().parse_args(argv[1:])
            seen.append(" ".join(argv[1:3]))
    for step in (
        "check env",
        "build ref",
        "check chemistry",
        "run -t",
        "check run",
        "status --run-dir",
    ):
        assert any(row.startswith(step) for row in seen), step
    for block in re.findall(r"```python\n(.*?)```", text, re.S):
        compile(block, "docs/happy_path.md", "exec")


def test_happy_path_open_first_files_exist_in_the_run_layout():
    from viralscan.outputs import RESULT_FILES, SAMPLE_LAYOUT

    text = (ROOT / "docs/happy_path.md").read_text()
    section = text.split("## 6. What to open first")[1]
    named = re.findall(r"^\d\. \*\*`([^`]+?)\.?`", section, re.M)
    assert named[:3] == ["report.html", "results/viral_summary_core.tsv", "results/sensitivity.tsv"]
    known = {path for path, *_ in SAMPLE_LAYOUT} | {path for path, *_ in RESULT_FILES}
    assert set(named[:3]) <= known
    for page in ("README.md", "docs/results.md", "docs/output_reference.md"):
        assert "happy_path.md" in (ROOT / page).read_text(), page

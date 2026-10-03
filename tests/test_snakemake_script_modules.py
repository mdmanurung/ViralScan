"""Every module the Snakefile runs as a ``script:`` must import cleanly there.

Snakemake writes its preamble (``import sys; ... snakemake = pickle.loads(...)``)
at the **top** of a copy of the module, above the docstring. A ``from __future__``
import then stops being the first statement and the job dies with

    SyntaxError: from __future__ imports must occur at the beginning of the file

which only shows up at runtime, in the rule, after everything upstream has run.
Caught on the ANDET-09 acceptance run (job 25696097_2), where it surfaced after
the host filter and kb count had already completed.
"""

import ast
import re
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parents[1] / "src" / "viralscan" / "scripts"
SNAKEFILE = Path(__file__).parents[1] / "src" / "viralscan" / "Snakefile"


def snakemake_script_modules() -> list[Path]:
    """The ``script:`` targets named in the Snakefile."""
    text = SNAKEFILE.read_text(encoding="utf-8")
    names = re.findall(r'script:\s*\n\s*"scripts/([A-Za-z_0-9]+\.py)"', text)
    assert names, "no script: targets found; the Snakefile layout changed"
    return [SCRIPTS_DIR / n for n in sorted(set(names))]


@pytest.mark.parametrize("path", snakemake_script_modules(), ids=lambda p: p.name)
def test_no_future_import_in_a_snakemake_script(path: Path) -> None:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    offenders = [
        node.lineno
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module == "__future__"
    ]
    assert not offenders, (
        f"{path.name}:{offenders} has a __future__ import. Snakemake inserts its "
        "preamble above the docstring, so this is a SyntaxError inside the rule. "
        "Delete the import; PEP 585 generics evaluate natively on the supported "
        "Python versions."
    )

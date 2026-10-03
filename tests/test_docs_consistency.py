"""Guard rails keeping user-facing docs in sync with runtime defaults.

These tests read the shipped Markdown docs and assert they do not drift back to
stale version examples or the pre-2.5 multimap default. They run from the repo
root (pytest rootdir), so the relative doc paths resolve directly.
"""

from pathlib import Path

from viralscan import __version__, anello_align
from viralscan.defaults import DEFAULT_MULTIMAP_METHOD

DOCS = [
    Path("README.md"),
    Path("docs/installation.md"),
    Path("docs/quickstart.md"),
    Path("docs/cli_reference.md"),
    Path("docs/output_reference.md"),
    Path("docs/api.md"),
]


def test_docs_do_not_reference_stale_container_version():
    for path in DOCS:
        text = path.read_text()
        assert "2.3.0" not in text, f"{path} still references 2.3.0"
        assert "viralscan_2.3.0" not in text, f"{path} still references old SIF name"


def test_runtime_multimap_default_is_documented_consistently():
    # Default reverted to host-conservative (2026-07-05) for viral-detection specificity.
    assert DEFAULT_MULTIMAP_METHOD == "host-conservative"
    stale_phrases = [
        "default is `equal`",
        "default multimapping method is `equal`",
        "The default `--multimap-method equal`",
        "The default multimapping method is `equal`",
        "default `multimap_method` is `equal`",
    ]
    for path in DOCS:
        text = path.read_text()
        for phrase in stale_phrases:
            assert phrase not in text, f"{path} contains stale default phrase: {phrase}"


def test_cli_reference_mentions_current_version():
    text = Path("docs/cli_reference.md").read_text()
    assert f"**{__version__}**" in text


def test_output_reference_documents_called_cell_column():
    text = Path("docs/output_reference.md").read_text()
    assert "is_called_cell" in text


# --- DOC-01 / DOC-03 additions -------------------------------------------------

import re  # noqa: E402

from viralscan.multimapping import MULTIMAP_EVIDENCE_COLUMNS  # noqa: E402

GUARDED = DOCS + [Path("docs/faq.md"), Path("docs/reference_panel.md"), Path("docs/migration.md")]
KINDS = {
    "observation",
    "model estimate",
    "evidence tier",
    "diagnostic flag",
    "biological interpretation",
}
COLUMN_HEADERS = {"column", "field", "output", "layer"}


def test_guarded_docs_do_not_reference_stale_versions():
    for path in GUARDED:
        text = path.read_text()
        for stale in ("2.3.0", "2.5.0"):
            assert stale not in text, f"{path} still references {stale}"


def test_faq_and_panel_docs_have_no_stale_multimap_default():
    for path in (Path("docs/faq.md"), Path("docs/reference_panel.md")):
        assert "default is `equal`" not in path.read_text(), path


def test_covid_alphatorquevirus_counts_removed():
    for path in GUARDED + [Path("BENCHMARK_COMPARISON.md")]:
        text = path.read_text()
        for stale in ("1,167,103", "1,002,218", "57,715", "~4× more sensitive"):
            assert stale not in text, f"{path} still contains {stale}"


def _cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def _tables(text):
    """Yield (header_cells, rows) for each pipe table in Markdown text."""
    lines = text.splitlines()
    i = 0
    while i < len(lines) - 1:
        if lines[i].startswith("|") and re.match(r"\|[-| :]+\|$", lines[i + 1]):
            j = i + 2
            while j < len(lines) and lines[j].startswith("|"):
                j += 1
            yield _cells(lines[i]), [_cells(r) for r in lines[i + 2 : j]]
            i = j
        else:
            i += 1


def test_output_reference_column_tables_have_kind_labels():
    text = Path("docs/output_reference.md").read_text()
    tables = [(h, r) for h, r in _tables(text) if h[0].lower() in COLUMN_HEADERS]
    assert len(tables) >= 11
    for header, rows in tables:
        assert header[-1] == "Kind", f"table starting {header} lacks a Kind column"
        for row in rows:
            assert row[-1] in KINDS, f"{row[0]}: bad Kind label {row[-1]!r}"


def test_emitted_columns_are_documented():
    text = Path("docs/output_reference.md").read_text()
    src = Path("src/viralscan/scripts/detection.py").read_text()
    block = src.split("def write_tsv_outputs(", 1)[1].split("    columns = [", 1)[1].split("]", 1)[0]
    summary_cols = re.findall(r'"(\w+)"', block)
    assert "viral_molecules_total_est" in summary_cols and len(summary_cols) >= 16
    # ANDET-09 alignment columns are appended from anello_align.SUMMARY_COLUMNS.
    summary_cols += list(anello_align.SUMMARY_COLUMNS)
    for col in summary_cols + list(MULTIMAP_EVIDENCE_COLUMNS):
        assert f"`{col}`" in text, f"{col} not documented in output_reference.md"


def test_migration_doc_is_in_index_and_states_rebuild_only():
    index = Path("docs/index.md").read_text()
    assert "migration" in index and "SUPPORT.md" in index and "SECURITY.md" in index
    assert "does **not** migrate" in Path("docs/migration.md").read_text()

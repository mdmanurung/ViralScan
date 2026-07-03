"""Guard rails keeping user-facing docs in sync with runtime defaults.

These tests read the shipped Markdown docs and assert they do not drift back to
stale version examples or the pre-2.5 multimap default. They run from the repo
root (pytest rootdir), so the relative doc paths resolve directly.
"""
from pathlib import Path

from viralscan import __version__
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
    assert DEFAULT_MULTIMAP_METHOD == "equal"
    stale_phrases = [
        "default is `host-conservative`",
        "default multimapping method is `host-conservative`",
        "default `--multimap-method host-conservative`",
        "The default `--multimap-method host-conservative`",
        "`host-conservative` | Multimapper allocation",
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

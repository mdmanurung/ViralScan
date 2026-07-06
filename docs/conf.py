# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

# Make the package importable without installing it.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

# -- Project information -------------------------------------------------------

project = "ViralScan"
copyright = "2026, Emma Vonk (Leiden University Medical Centre)"
author = "Emma Vonk"

# Single-source the version from the package (src/viralscan/__init__.py), so the
# docs never advertise a stale release. Fall back gracefully if the import fails.
try:
    from viralscan import __version__ as release
except Exception:  # noqa: BLE001 — a doc build must not hard-fail on version import
    release = "0.0.0"
version = ".".join(release.split(".")[:2])

# -- General configuration -----------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # NumPy / Google docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",  # render type annotations in docs
    "myst_parser",  # Markdown source files
    "nbsphinx",  # Jupyter notebook vignettes
]

# MyST extensions
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    # Internal planning / manuscript drafts — not part of the public docs build.
    "manuscript_draft.md",
    "review-*.md",
    "write-docs-prompt.md",
    "showcase_runbook.md",
]

# Source file suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
    ".ipynb": "nbsphinx",
}

# The vignettes include real workflow commands and external downloads. Render
# the saved notebook content in docs builds instead of executing them.
nbsphinx_execute = "never"

# The root document
root_doc = "index"

# Intersphinx — link to external docs
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable", None),
}

# -- Options for HTML output ---------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": True,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 4,
}

# -- Autodoc configuration -----------------------------------------------------

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"
add_module_names = False

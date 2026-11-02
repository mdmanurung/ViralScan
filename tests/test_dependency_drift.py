"""REVIEW-FIX-H: guards for scanpy / scipy / pandas version drift."""

import re
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

from viralscan.scripts import cellcalling

ROOT = Path(__file__).parents[1]


def test_no_deprecated_use_highly_variable():
    assert "use_highly_variable" not in (ROOT / "src/viralscan/scripts/umap.py").read_text()


def test_scanpy_lower_bound_in_pyproject():
    text = (ROOT / "pyproject.toml").read_text()
    assert len(re.findall(r'"scanpy(\[leiden\])?>=1\.10"', text)) == 2


def test_host_matrix_totals_accept_coo_array(tmp_path, monkeypatch):
    """scipy 1.20 mmread returns coo_array; the reader must not depend on the class."""
    (tmp_path / "barcodes.tsv").write_text("AAA\nCCC\n")
    (tmp_path / "features.tsv").write_text("g1\ng2\ng3\n")
    (tmp_path / "matrix.mtx").write_text("")
    arr = sp.coo_array(np.array([[1, 0], [0, 2], [3, 0]]))
    monkeypatch.setattr(scipy.io, "mmread", lambda _path: arr)
    barcodes, totals = cellcalling.read_host_matrix_totals(tmp_path)
    assert barcodes == ["AAA", "CCC"] and totals.tolist() == [4, 2]

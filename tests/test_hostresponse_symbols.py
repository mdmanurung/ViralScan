"""C1: host-response Ensembl->symbol annotation (pure, network-free logic)."""

from __future__ import annotations

import pandas as pd

from viralscan.scripts.hostresponse import (
    _add_symbol_column,
    _looks_like_ensembl,
    _map_ensembl_to_symbols,
)


def test_look_like_ensembl():
    assert _looks_like_ensembl(["ENSG00000141510", "ENSG00000012048"]) is True
    assert _looks_like_ensembl(["GAPDH", "MT-CO1", "TP53"]) is False


def test_add_symbol_column_inserts_next_to_gene():
    df = pd.DataFrame({"gene": ["ENSG00000141510.17", "ENSG00000012048.23"], "weight": [1.0, 2.0]})
    mapping = {"ENSG00000141510": "TP53", "ENSG00000012048": "BRCA1"}
    out = _add_symbol_column(df, mapping)
    assert list(out.columns) == ["gene", "symbol", "weight"]  # inserted right after gene
    assert out["symbol"].tolist() == ["TP53", "BRCA1"]  # version suffix stripped for lookup


def test_add_symbol_column_is_noop_without_map_or_gene_col():
    df = pd.DataFrame({"gene": ["ENSG1"], "w": [1.0]})
    assert "symbol" not in _add_symbol_column(df, {}).columns  # empty map -> no column
    assert "symbol" not in _add_symbol_column(pd.DataFrame({"x": [1]}), {"a": "b"}).columns


def test_map_ensembl_to_symbols_empty_input_no_network():
    # empty input returns {} before any network call
    assert _map_ensembl_to_symbols([]) == {}

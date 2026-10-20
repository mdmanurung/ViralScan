"""Exact figure data and universe contracts; no scientific performance claims."""

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan.scripts import detection as d


def draw(tmp_path, x, host, called):
    a = ad.AnnData(sparse.csr_matrix(x, dtype=float))
    a.var_names = ["host", "virus"]
    a.obs_names = [f"C{i}" for i in range(len(x))]
    before = a.X.copy()
    d.super_expressor(a, "V", ["virus"], tmp_path, host_gene_mask=host, called_mask=called)
    assert (before != a.X).nnz == 0
    assert not len(a.obs.columns)
    stem = tmp_path / "plots" / "SuperExpressor_V"
    return pd.read_csv(str(stem) + ".tsv", sep="\t"), json.loads(open(str(stem) + ".json").read())


def test_host_only_null_and_called_universe(tmp_path):
    table, meta = draw(
        tmp_path,
        [[1, 10], [3, 0], [1000, 100000]],
        np.array([True, False]),
        np.array([True, True, False]),
    )
    assert table.viral_molecules.tolist() == [10, 0]
    assert table.null_expected.tolist() == [2.5, 7.5]
    assert meta["n_universe"] == 2 and meta["n_zero"] == 1 and meta["n_high_burden"] == 1
    assert meta["null_status"] == "ok"


@pytest.mark.parametrize("host", [None, np.array([True, False])])
def test_missing_host_depth(tmp_path, host):
    table, meta = draw(tmp_path, [[0, 0.3], [0, 0]], host, np.ones(2, dtype=bool))
    assert meta["null_status"] == "not_estimable"
    assert table.null_expected.isna().all()
    assert table.viral_molecules.tolist() == [0.3, 0]


def test_figure_only_rank_cap(tmp_path, monkeypatch):
    labels = []
    from matplotlib.axes import Axes

    original = Axes.text

    def capture(self, x, y, s, *args, **kwargs):
        labels.append(s)
        return original(self, x, y, s, *args, **kwargs)

    monkeypatch.setattr(Axes, "text", capture)
    table, meta = draw(
        tmp_path, np.tile([1, 0], (60000, 1)), np.array([True, False]), np.ones(60000, dtype=bool)
    )
    assert len(table) == meta["n_universe"] == 60000
    assert meta["n_plotted"] == meta["max_rank_display"] == 50000
    assert any("showing ranks 1–50000 of 60000" in s for s in labels)


def test_empty_called_universe(tmp_path):
    table, meta = draw(tmp_path, [[1, 1]], np.array([True, False]), np.array([False]))
    assert table.empty and meta["n_universe"] == 0


def test_invalid_masks(tmp_path):
    with pytest.raises(ValueError, match="called_mask"):
        draw(tmp_path, [[1, 1]], None, np.array([1]))


def test_plot_preserves_exact_counts_layers_and_serialized_detection(tmp_path):
    a = ad.AnnData(sparse.csr_matrix([[1.0, 10.0], [3.0, 0.3], [100.0, 20.0]]))
    a.var_names = ["host", "virus"]
    a.obs_names = ["A", "B", "C"]
    a.layers["counts_original"] = a.X.copy()
    called = np.array([True, True, False])
    matrices = [a.X, a.layers["counts_original"]]
    before = [
        (m.dtype.str, m.data.tobytes(), m.indices.tobytes(), m.indptr.tobytes()) for m in matrices
    ]
    kwargs = dict(called_mask=called, viral_count_matrix=a.X, index_genes_by_virus={"V": ["virus"]})
    stats, per_cell = d.compute_stats(a, {"virus": 30.3}, {"V": ["virus"]}, {"V"}, **kwargs)
    d.write_tsv_outputs(stats, per_cell, tmp_path / "before")
    d.super_expressor(
        a, "V", ["virus"], tmp_path, host_gene_mask=np.array([True, False]), called_mask=called
    )
    stats_after, per_after = d.compute_stats(a, {"virus": 30.3}, {"V": ["virus"]}, {"V"}, **kwargs)
    d.write_tsv_outputs(stats_after, per_after, tmp_path / "after")
    for name in ["viral_summary.tsv", "per_cell_viral.tsv"]:
        assert (tmp_path / "before" / "results" / name).read_bytes() == (
            tmp_path / "after" / "results" / name
        ).read_bytes()
    assert before == [
        (m.dtype.str, m.data.tobytes(), m.indices.tobytes(), m.indptr.tobytes()) for m in matrices
    ]
    assert type(stats["V"]["detected"]) is type(stats_after["V"]["detected"])
    assert stats["V"]["detected"] == stats_after["V"]["detected"]

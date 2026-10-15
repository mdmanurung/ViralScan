"""AUD3-01b (audit P9-03): float under-shoot of a split molecule must not drop a virus."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.scripts import detection
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable


def _split_row(k, mass=None):
    """One cell, host gene h1, and k viral genes that each hold `1/k` (float64)."""
    share = 1.0 / k if mass is None else mass / k
    return np.array([50.0, *([share] * k)], dtype=np.float64)


@pytest.mark.parametrize("k", [6, 7])
def test_numpy_float64_split_really_undershoots(k):
    """Guard: the premise of the bug (sum of k copies of 1/k is < 1 in numpy)."""
    assert np.full(k, 1.0 / k).sum() < 1.0


@pytest.mark.parametrize("k", [6, 7])
def test_detect_genes_keeps_k_way_split_molecule(k):
    genes = [f"v{i}" for i in range(k)]
    X = sp.csr_matrix(_split_row(k)[1:][None, :])
    found = detection.detect_genes(genes, X, genes, 1, groups={"V": genes})
    assert set(found) == set(genes)


def test_detect_genes_still_rejects_subthreshold_mass():
    genes = [f"v{i}" for i in range(6)]
    X = sp.csr_matrix(_split_row(6, mass=0.99)[1:][None, :])
    assert detection.detect_genes(genes, X, genes, 1, groups={"V": genes}) == {}


@pytest.mark.parametrize("k", [6, 7])
def test_compute_stats_called_gate_keeps_k_way_split_molecule(k):
    genes = [f"v{i}" for i in range(k)]
    a = ad.AnnData(
        X=sp.csr_matrix(_split_row(k)[None, :]),
        obs=pd.DataFrame(index=["C0"]),
        var=pd.DataFrame(index=["h1", *genes]),
    )
    found = {g: 1.0 / k for g in genes}
    stats, _ = detection.compute_stats(a, found, {"V": genes}, genes)
    assert stats["V"]["detected"] is True


def test_compute_stats_called_gate_still_rejects_subthreshold_mass():
    k = 6
    genes = [f"v{i}" for i in range(k)]
    a = ad.AnnData(
        X=sp.csr_matrix(_split_row(k, mass=0.99)[None, :]),
        obs=pd.DataFrame(index=["C0"]),
        var=pd.DataFrame(index=["h1", *genes]),
    )
    found = {g: 0.99 / k for g in genes}
    stats, _ = detection.compute_stats(a, found, {"V": genes}, genes)
    assert stats["V"]["detected"] is False


def _gi(gid, viral, key="", name=""):
    return GeneIdentity(
        gene_id=gid,
        genome_accession=gid,
        status="catalogued" if viral else "host",
        viral=viral,
        virus_key=key,
        virus_name=name,
    )


@pytest.mark.parametrize("k, mass, expect", [(6, None, True), (7, None, True), (6, 0.99, False)])
def test_viral_summary_through_detection_run(tmp_path, k, mass, expect):
    genes = [f"v{i}" for i in range(k)]
    table = [_gi("h1", False), *[_gi(g, True, "V", "Virus V") for g in genes]]
    run = tmp_path
    (run / "results").mkdir(parents=True)
    (run / "log").mkdir()
    VirusIdentityTable(tuple(table)).write_tsv(run / "results" / "virus_identity.tsv")
    (run / "log" / "analysis.txt").write_text("\n".join(genes) + "\n")
    # 20 cells so emptyDrops is not needed (cell_calling none); the molecule sits in cell 0.
    X = np.zeros((20, 1 + k))
    X[:, 0] = 50.0
    X[0] = _split_row(k, mass)
    a = ad.AnnData(X=sp.csr_matrix(X))  # float64 layer, as the pipeline produces
    a.obs_names = [f"CELL{i}" for i in range(20)]
    a.var_names = ["h1", *genes]
    d = run / "kb-python" / "counts_unfiltered"
    d.mkdir(parents=True)
    a.write_h5ad(d / "adata.h5ad")
    (d / "cells_x_genes.genes.txt").write_text("\n".join(["h1", *genes]) + "\n")
    cfg = RunConfig(
        output=str(run) + "/",
        multimapping=False,
        visual=False,
        cell_calling="none",
        detection_threshold=1,
        detection_min_cells=1,
    )
    detection.run(
        RunContext.from_config(cfg), str(run / "log" / "analysis.txt"), str(run / "log" / "done")
    )
    summary = pd.read_csv(run / "results" / "viral_summary.tsv", sep="\t")
    if expect:
        row = summary.loc[summary["virus_name"].astype(str).str.contains("V")]
        assert len(row) == 1 and bool(row["detected"].iloc[0]) is True
    else:
        assert summary.empty or not summary["detected"].astype(bool).any()

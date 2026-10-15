"""AUD3-01a (audit P9-02): the Virus Identity table must match the kb matrix on every detection path."""

import logging

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.scripts import detection
from viralscan.virus_grouping import IdentityMatrixMismatch, reconcile_identity_to_matrix
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable


def _gi(gid, viral, key="", name=""):
    return GeneIdentity(
        gene_id=gid,
        genome_accession=gid,
        status="catalogued" if viral else "host",
        viral=viral,
        virus_key=key,
        virus_name=name,
    )


def _build(run, table_genes, matrix_genes, X):
    (run / "results").mkdir(parents=True)
    (run / "log").mkdir()
    VirusIdentityTable(tuple(table_genes)).write_tsv(run / "results" / "virus_identity.tsv")
    (run / "log" / "analysis.txt").write_text(
        "\n".join(g.gene_id for g in table_genes if g.viral) + "\n"
    )
    a = ad.AnnData(X=sp.csr_matrix(np.asarray(X, dtype=np.float32)))
    a.obs_names = [f"CELL{i}" for i in range(a.n_obs)]
    a.var_names = matrix_genes
    d = run / "kb-python" / "counts_unfiltered"
    d.mkdir(parents=True)
    a.write_h5ad(d / "adata.h5ad")
    (d / "cells_x_genes.genes.txt").write_text("\n".join(matrix_genes) + "\n")


def _run_detection(run):
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


_TABLE = [
    _gi("h1", False),
    _gi("h2", False),
    _gi("vA_g1", True, "A", "Virus A"),
    _gi("vB_g1", True, "B", "Virus B"),
]


def test_renamed_gene_id_is_an_error_not_a_silent_drop(tmp_path):
    """Scenario 1: table lists vB_g1, the matrix holds it renamed (vB_g1.2) in every cell."""
    X = np.zeros((30, 4))
    X[:, :2] = 50
    X[0, 2] = 3
    X[:, 3] = 5
    _build(tmp_path, _TABLE, ["h1", "h2", "vA_g1", "vB_g1.2"], X)
    with pytest.raises(IdentityMatrixMismatch, match=r"vB_g1\.2.*another run"):
        _run_detection(tmp_path)


def test_table_virus_absent_from_matrix_only_warns(tmp_path, caplog):
    """Scenario 2: the table lists a virus with no matrix column; no counts can exist, so no error."""
    table = [*_TABLE, _gi("vC_g1", True, "C", "Virus C")]
    X = np.zeros((30, 4))
    X[:, :2] = 50
    X[:3, 2] = 3
    X[:3, 3] = 4
    _build(tmp_path, table, ["h1", "h2", "vA_g1", "vB_g1"], X)
    with caplog.at_level(logging.WARNING):
        _run_detection(tmp_path)
    warned = [r for r in caplog.records if "vC_g1" in r.getMessage()]
    assert len(warned) == 1
    assert (tmp_path / "results" / "viral_summary.tsv").is_file()


def test_reconcile_helper_exact_match_is_silent(caplog):
    table = VirusIdentityTable(tuple(_TABLE))
    with caplog.at_level(logging.WARNING):
        reconcile_identity_to_matrix(table, ["h1", "h2", "vA_g1", "vB_g1"], "t.tsv")
    assert not caplog.records

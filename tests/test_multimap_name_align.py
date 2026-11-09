"""Name alignment and loud failure on missing genes in scripts/multimap (MULTIMAP-ALIGN-01)."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan.scripts import multimap


def _adata(barcodes, genes):
    # X[i, j] encodes the pair, so a swap shows in the value.
    x = np.array([[10.0 * i + j for j in range(len(genes))] for i in range(len(barcodes))])
    return ad.AnnData(
        X=sparse.csr_matrix(x),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes),
    )


def test_align_reorders_by_name_when_h5ad_order_differs():
    original = _adata(["B-1", "A-1", "C-1"], ["g2", "g1"])
    to_idx = {"A": 0, "B": 1, "C": 2}  # barcodes.txt order after the -1 strip
    aligned = multimap.align_adata_to_text_files(original, to_idx, ["g1", "g2"])
    assert list(aligned.var_names) == ["g1", "g2"]
    assert list(aligned.obs_names) == ["A-1", "B-1", "C-1"]
    assert aligned.X[0, 0] == 11.0  # A / g1 was row 1, col 1
    assert aligned.X[1, 1] == 0.0  # B / g2 was row 0, col 0
    assert aligned.X[2, 0] == 21.0  # C / g1 was row 2, col 1


def test_align_returns_same_object_when_order_matches():
    original = _adata(["A", "B"], ["g1", "g2"])
    assert multimap.align_adata_to_text_files(original, {"A": 0, "B": 1}, ["g1", "g2"]) is original


def test_align_raises_on_different_barcode_set():
    original = _adata(["A", "B"], ["g1"])
    with pytest.raises(ValueError, match="1 only in the h5ad.*1 only in barcodes.txt"):
        multimap.align_adata_to_text_files(original, {"A": 0, "C": 1}, ["g1"])


def test_align_raises_on_different_gene_set():
    original = _adata(["A"], ["g1", "gX"])
    with pytest.raises(ValueError, match=r"genes: 1 only in the h5ad \(e.g. \['gX'\]\)"):
        multimap.align_adata_to_text_files(original, {"A": 0}, ["g1", "g2"])


def _write_ec(tmp_path, lines):
    path = tmp_path / "matrix.ec"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_read_ec_maps_transcripts_to_gene_indices(tmp_path):
    ec = _write_ec(tmp_path, ["0\t0", "1\t0,1"])
    got = multimap.read_ec(ec, ["t0", "t1"], {"t0": "g0", "t1": "g1"}, ["g0", "g1"])
    assert got == {0: [0], 1: [0, 1]}


def test_read_ec_raises_when_a_transcript_is_missing_from_t2g(tmp_path):
    ec = _write_ec(tmp_path, ["0\t0,1"])
    with pytest.raises(ValueError, match=r"1 transcripts missing from the t2g.*\['t1'\]"):
        multimap.read_ec(ec, ["t0", "t1"], {"t0": "g0"}, ["g0", "g1"])


def test_read_ec_raises_when_a_t2g_gene_is_not_in_genes_txt(tmp_path):
    ec = _write_ec(tmp_path, ["0\t0,1"])
    with pytest.raises(ValueError, match=r"1 t2g genes missing from genes.txt.*\['gZ'\]"):
        multimap.read_ec(ec, ["t0", "t1"], {"t0": "g0", "t1": "gZ"}, ["g0", "g1"])


def test_read_ec_raises_on_transcript_index_past_transcripts_txt(tmp_path):
    ec = _write_ec(tmp_path, ["0\t0,5"])
    with pytest.raises(ValueError, match="1 transcript indices beyond transcripts.txt"):
        multimap.read_ec(ec, ["t0"], {"t0": "g0"}, ["g0"])


def test_read_ec_warns_about_ecs_with_no_transcripts(tmp_path, caplog):
    ec = _write_ec(tmp_path, ["0\t0", "1\t,"])
    with caplog.at_level("WARNING", logger="viralscan"):
        got = multimap.read_ec(ec, ["t0"], {"t0": "g0"}, ["g0"])
    assert got == {0: [0]}
    assert "1 of 2 ECs list no transcripts" in caplog.text

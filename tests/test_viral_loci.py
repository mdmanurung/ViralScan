"""LOCK-01 D1: independent loci from overlapping annotated intervals (TEST_ONLY fixture, no real annotation)."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from viralscan.multimapping import summarize_multimap_evidence
from viralscan.runconfig import RunConfig
from viralscan.viral_loci import assign_loci, gene_intervals, read_loci, write_loci

GTF = Path(__file__).parent / "data" / "viral_loci_TEST_ONLY.gtf"


def _loci() -> dict[str, str]:
    return assign_loci(gene_intervals(GTF.read_text().splitlines()))


def test_overlap_is_at_least_one_base_on_either_strand():
    loci = _loci()
    assert loci["A1"] == loci["A2"]  # opposite strands, overlapping
    assert loci["A3"] != loci["A2"]  # 300 | 301: adjacent, no shared base
    assert loci["A3"] == loci["A4"]  # shares base 400
    assert loci["A5"] not in {loci["A1"], loci["A3"]}
    assert loci["B1"] != loci["A1"]  # same coordinates on another sequence


def test_origin_wrapping_gene_merges_across_the_origin():
    loci = _loci()
    assert loci["C1"] == loci["C2"]  # C2 overlaps only the part of C1 at the start of the genome
    assert loci["C3"] != loci["C1"]


def test_cds_only_gene_takes_the_cds_span():
    loci = _loci()
    assert loci["A6"] == loci["A7"]  # A6 spans 1000-1250 and A7 starts at 1240


def test_write_and_read_roundtrip_restricts_to_requested_genes(tmp_path):
    out = write_loci([GTF], {"A1", "A2", "A5"}, tmp_path / "results" / "viral_loci.tsv")
    got = read_loci(tmp_path / "results")
    assert out is not None and set(got) == {"A1", "A2", "A5"} and got["A1"] == got["A2"] != got["A5"]
    assert write_loci([GTF], {"nope"}, tmp_path / "x.tsv") is None
    assert read_loci(tmp_path / "missing") == {}


class TestLockOneLocusAndPositiveCell:
    """LOCK-01 D1/D2 columns on a TEST_ONLY matrix: positive cells are a union over genes, among called barcodes."""

    def _adata(self):
        import anndata as ad

        adata = ad.AnnData(
            X=sparse.csr_matrix(np.zeros((4, 3))),
            obs=pd.DataFrame(index=["BC1", "BC2", "BC3", "BC4"]),
            var=pd.DataFrame(index=["gA", "gB", "gC"]),
        )
        # gA: BC1, BC2; gB: BC2, BC3; gC (other virus): BC4
        adata.layers["counts_unique_viral"] = sparse.csr_matrix(
            np.array([[1.0, 0, 0], [2.0, 1.0, 0], [0, 3.0, 0], [0, 0, 5.0]])
        )
        return adata

    def test_positive_cells_are_a_union_not_a_maximum(self) -> None:
        res = summarize_multimap_evidence(
            self._adata(), {"V": ["gA", "gB"], "W": ["gC"]}, RunConfig(multimap_method="equal")
        )
        assert res.loc[res.virus_name == "V", "n_unique_positive_cells_called"].tolist() == [3, 3]
        assert res.loc[res.gene_id == "gA", "n_unique_viral_cells"].item() == 2  # per-gene count unchanged
        assert res.loc[res.virus_name == "W", "n_unique_positive_cells_called"].item() == 1

    def test_only_called_barcodes_count_and_locus_is_passed_through(self) -> None:
        res = summarize_multimap_evidence(
            self._adata(),
            {"V": ["gA", "gB"]},
            RunConfig(multimap_method="equal"),
            called_mask=np.array([True, True, False, True]),
            locus_by_gene={"gA": "L1", "gB": "L1"},
        )
        assert res["n_unique_positive_cells_called"].tolist() == [2, 2]
        assert res["locus_id"].tolist() == ["L1", "L1"]

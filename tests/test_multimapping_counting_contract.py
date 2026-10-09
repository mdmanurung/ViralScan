"""Independent hand-counted contracts for opt-in sibling allocation (DEF-06)."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan.kb_outputs import KbCountOutputs
from viralscan.multimapping import build_multimap_layers
from viralscan.runconfig import RunConfig
from viralscan.scripts import multimap
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable


def _fixture():
    rows = [("A", f"H{i}", 0, 7) for i in range(6)]
    rows += [("A", f"V{i}", 1, 3) for i in range(3)]
    rows += [("A", "W", 2, 1)]
    rows += [
        ("A", "S", 3, 8),  # sibling-only
        ("A", "I", 4, 2),
        ("A", "I", 3, 9),  # EC intersection -> siblings
        ("A", "M", 5, 4),  # host + two siblings
        ("A", "G", 6, 1),  # cross-group
        ("A", "N", 7, 1),  # unknown sibling
        ("A", "C", 0, 1),
        ("A", "C", 8, 1),  # unresolved collision
        ("B", "S", 3, 10),  # no local unique support
    ]
    frame = pd.DataFrame(rows, columns=["barcode", "umi", "ec", "count"])
    ecs = {
        0: [0],
        1: [1, 1],
        2: [2],
        3: [1, 2],
        4: [3, 2, 1],
        5: [0, 1, 2],
        6: [1, 3],
        7: [1, 4],
        8: [3],
    }
    return frame, ecs


def _build(frame, ecs, **kwargs):
    return build_multimap_layers(
        frame,
        {"A": 0, "B": 1},
        ecs,
        2,
        5,
        {1, 2, 3, 4},
        sparse.csr_matrix((2, 5)),
        **kwargs,
    )


@pytest.mark.parametrize("seed", range(12))
def test_sibling_allocation_hand_counted_mass_membership_and_order(seed):
    frame, ecs = _fixture()
    shuffled_ecs = {ec: list(reversed(genes)) for ec, genes in reversed(list(ecs.items()))}
    result = _build(
        frame.sample(frac=1, random_state=seed),
        shuffled_ecs,
        method="sibling-weighted",
        sibling_groups={1: "S", 2: "S", 3: "T"},
    )
    # Two S-only molecules weight 4:2 from local unique counts 3:1 + pseudocount.
    # Mixed, cross-group and unknown-group molecules keep their equal allocation.
    np.testing.assert_allclose(
        result.corrected.toarray(),
        [
            [1 / 3, 8 / 3, 1, 1 / 2, 1 / 2],
            [0, 1 / 2, 1 / 2, 0, 0],
        ],
    )
    np.testing.assert_allclose(result.unique.toarray(), [[6, 3, 1, 0, 0], [0] * 5])
    assert result.audit.input_molecules == 17
    assert result.audit.unique_molecules == 10
    assert result.audit.ambiguous_molecules == 6
    assert result.audit.unresolved_molecules == 1
    assert result.corrected.sum() == pytest.approx(6)
    assert result.host_virus_ambiguous_molecules == 1
    assert result.host_viral_selected.sum() == pytest.approx(2 / 3)
    for matrix in (result.unique, result.corrected, result.sibling_weighted):
        assert np.isfinite(matrix.data).all() and (matrix.data >= 0).all()


def test_sibling_zero_support_pseudocount_fraction_and_stream_buffers(tmp_path):
    frame, ecs = _fixture()
    ordered = frame.sort_values(["barcode", "umi", "ec"])
    path = tmp_path / "bus.txt"
    ordered.to_csv(path, sep="\t", header=False, index=False)
    expected = _build(
        frame,
        ecs,
        method="sibling-weighted",
        pseudocount=2,
        sibling_groups={1: "S", 2: "S", 3: "T"},
    )
    # Local A unique support gives 5:3; zero-support B still gives exactly 1:1.
    assert expected.corrected[0, 2] == pytest.approx(2 * 3 / 8 + 1 / 3)
    assert expected.corrected[1, 1] == pytest.approx(1 / 2)
    for buffer_size in (2, 8192):
        streamed = _build(
            path,
            ecs,
            method="sibling-weighted",
            pseudocount=2,
            sibling_groups={1: "S", 2: "S", 3: "T"},
            bus_buffer_size=buffer_size,
        )
        np.testing.assert_allclose(streamed.corrected.toarray(), expected.corrected.toarray())
        assert streamed.audit == expected.audit


def test_sibling_missing_identity_fails_but_empty_curated_membership_is_equal():
    frame, ecs = _fixture()
    with pytest.raises(ValueError, match="Virus Identity"):
        _build(frame, ecs, method="sibling-weighted")
    with pytest.raises(ValueError, match="indexed viral"):
        _build(frame, ecs, method="sibling-weighted", sibling_groups={0: "S"})
    result = _build(frame, ecs, method="sibling-weighted", sibling_groups={})
    np.testing.assert_allclose(result.corrected.toarray(), result.equal.toarray())
    assert result.method_diagnostics["sibling_identity_available"] is True


@pytest.mark.parametrize(
    "method", ["equal", "host-conservative", "unique-weighted", "em-global", "em-cell"]
)
def test_sibling_metadata_does_not_change_existing_allocations(method):
    frame, ecs = _fixture()
    old = _build(frame, ecs, method=method)
    metadata = _build(frame, ecs, method=method, sibling_groups={1: "S", 2: "S"})
    np.testing.assert_array_equal(old.corrected.toarray(), metadata.corrected.toarray())
    assert old.audit == metadata.audit


def test_host_virus_boundary_uses_final_intersection_not_observed_ec_union():
    frame = pd.DataFrame(
        [
            ("A", "U", 0, 20),
            ("A", "U", 1, 5),
            ("A", "M", 0, 7),
            ("A", "C", 1, 1),
            ("A", "C", 2, 1),
        ],
        columns=["barcode", "umi", "ec", "count"],
    )
    result = _build(frame, {0: [0, 1, 2], 1: [1], 2: [0]}, method="host-conservative")
    assert result.audit.unique_molecules == 1
    assert result.audit.ambiguous_molecules == 1
    assert result.audit.unresolved_molecules == 1
    assert result.host_virus_ambiguous_molecules == 1
    assert result.host_viral_ambiguous.sum() == pytest.approx(1)
    assert result.host_viral_selected.sum() == 0


def _identity():
    return VirusIdentityTable(
        (
            GeneIdentity("H", "", "host", False),
            GeneIdentity("V1", "A1", "catalogued", True, role="target", sibling_group="S"),
            GeneIdentity("V2", "A2", "catalogued", True, role="decoy", sibling_group="S"),
            GeneIdentity("V3", "A3", "catalogued", True, role="whole_genome"),
            GeneIdentity("V4", "A4", "uncatalogued", True),
        )
    )


def test_gene_role_catalogue_and_count_audit_roundtrip(tmp_path, monkeypatch):
    frame, ecs = _fixture()
    layers = _build(frame, ecs, method="sibling-weighted", sibling_groups={1: "S", 2: "S"})
    ids = [row.gene_id for row in _identity().genes]
    original = ad.AnnData(
        X=layers.unique, obs=pd.DataFrame(index=["A", "B"]), var=pd.DataFrame(index=ids)
    )
    matrix, viral_counts = multimap.create_new_h5ad(
        layers.unique + layers.corrected,
        original,
        ids,
        ids,
        ids,
        {1, 2, 3, 4},
        5,
    )
    multimap.stamp_gene_identity(matrix, ids, _identity())
    monkeypatch.setattr(
        multimap, "config", RunConfig(output=str(tmp_path), multimap_method="sibling-weighted")
    )
    outputs = KbCountOutputs.from_config_output(str(tmp_path))
    outputs.adata_multimap.parent.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(multimap, "kb", outputs)
    multimap.final_results(viral_counts, original, {1, 2, 3, 4}, matrix, 2, layers)
    loaded = ad.read_h5ad(outputs.adata_multimap)
    assert list(loaded.var["gene_role"]) == ["host", "target", "decoy", "whole_genome", "unknown"]
    assert loaded.uns["index_kind"] == "combined"
    assert loaded.uns["molecule_denominator"] == "all_indexed_genes"
    assert loaded.uns["molecule_audit"]["host_virus_ambiguous_molecules"] == 1
    np.testing.assert_allclose(loaded.X.toarray(), (layers.unique + layers.corrected).toarray())
    np.testing.assert_allclose(
        loaded.layers["counts_multimap_sibling_weighted"].toarray(), layers.corrected.toarray()
    )
    audit = pd.read_csv(tmp_path / "count_audit.tsv", sep="\t")
    assert audit.loc[0, "host_virus_ambiguous_molecules"] == 1


def test_gene_identity_requires_every_indexed_gene_and_consistent_partition():
    matrix = ad.AnnData(
        X=sparse.csr_matrix((1, 2)),
        var=pd.DataFrame({"is_viral": [False, True]}, index=["H", "MISSING"]),
    )
    with pytest.raises(ValueError, match="missing"):
        multimap.stamp_gene_identity(matrix, list(matrix.var_names), _identity())
    matrix.var_names = ["H", "V1"]
    matrix.var["is_viral"] = [False, False]
    with pytest.raises(ValueError, match="disagrees"):
        multimap.stamp_gene_identity(matrix, list(matrix.var_names), _identity())


@pytest.mark.parametrize(
    "viral,kind", [([True, True], "virus_only"), ([False, False], "host_only")]
)
def test_index_composition_and_legacy_source_are_explicit(tmp_path, monkeypatch, viral, kind):
    monkeypatch.setattr(multimap, "output", str(tmp_path))
    matrix = ad.AnnData(
        X=sparse.csr_matrix((1, 2)), var=pd.DataFrame({"is_viral": viral}, index=["G1", "G2"])
    )
    multimap.stamp_gene_identity(matrix, list(matrix.var_names), None)
    assert matrix.uns["index_kind"] == kind
    assert matrix.uns["gene_identity_source"] == "unavailable"
    Path(tmp_path, "log").mkdir()
    Path(tmp_path, "log", "analysis.txt").write_text("G1\nG2\n")
    multimap.stamp_gene_identity(matrix, list(matrix.var_names), None)
    assert matrix.uns["gene_identity_source"] == "legacy_analysis"

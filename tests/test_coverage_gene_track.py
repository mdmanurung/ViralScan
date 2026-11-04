"""Pure layout rules of the coverage-map gene track (PLAN EVID-05)."""

import pytest

from viralscan import coverage_gene_track as gt


def row(feature, start, end, strand="+", **attributes):
    return gt.GtfRow(feature, start, end, strand, attributes)


def test_labels_never_overlap_and_priority_wins():
    length = 100_000
    wanted = [(i, 50_000 + 500 * i, 8_000) for i in range(6)]  # all crowd one spot
    placed, unplaced = gt.place_labels(wanted, length, max_rows=3)
    assert [p.index for p in placed] == [0, 1, 2]  # best-covered first, one per row
    assert [p.row for p in placed] == [0, 1, 2]
    assert unplaced == [3, 4, 5]
    for r in {p.row for p in placed}:
        spans = sorted((p.low, p.high) for p in placed if p.row == r)
        assert all(a[1] < b[0] for a, b in zip(spans, spans[1:], strict=False))


def test_labels_share_a_row_when_apart_and_stay_inside_axis():
    placed, unplaced = gt.place_labels(
        [(0, 50, 1000), (1, 50_000, 1000), (2, 99_990, 1000)], 100_000
    )
    assert unplaced == [] and {p.row for p in placed} == {0}
    assert all(p.low >= 1 and p.high <= 100_000 for p in placed)


def test_long_genes_dropped_and_counted():
    rows = [
        row("gene", 100, 1_000, gene_id="A", gene_name="A"),
        row("gene", 1_000, 30_000, gene_id="LONG", gene_name="LONG"),
        row("gene", 5_000, 25_000, gene_id="EDGE", gene_name="EDGE"),  # 20,001 bp
    ]
    models, n_long = gt.build_gene_models(rows, 100_000)
    assert [m.name for m in models] == ["A"] and n_long == 2
    models, n_long = gt.build_gene_models(rows[:1] + rows[2:], 100_000, max_bp=20_001)
    assert [m.name for m in models] == ["A", "EDGE"] and n_long == 0


def test_cds_only_gene_splits_into_copies_and_reaches_utr_ends():
    rows = [
        row("gene", 90, 1_210, gene_id="U14", gene_name="U14"),
        row("CDS", 100, 1_200, gene_id="U14", gene_name="U14"),
        row("gene", 140_090, 141_210, gene_id="U14", gene_name="U14"),
        row("CDS", 140_100, 141_200, gene_id="U14", gene_name="U14"),
    ]
    models, n_long = gt.build_gene_models(rows, 150_000)
    assert n_long == 0
    assert [m.blocks for m in models] == [((90, 1_210),), ((140_090, 141_210),)]


def test_transcripts_replace_gene_row_and_exons_set_blocks():
    rows = [
        row("gene", 100, 5_000, gene_id="G", gene_name="G"),
        row("transcript", 100, 5_000, "-", gene_id="G", transcript_id="T1", gene_name="G"),
        row("exon", 100, 300, "-", gene_id="G", transcript_id="T1", gene_name="G"),
        row("exon", 4_000, 5_000, "-", gene_id="G", transcript_id="T1", gene_name="G"),
    ]
    (model,), _ = gt.build_gene_models(rows, 10_000)
    assert model.blocks == ((100, 300), (4_000, 5_000)) and model.strand == "-"
    assert model.size_bp == 4_901


def test_origin_spanning_row_counts_both_parts():
    models, n_long = gt.build_gene_models([row("gene", 9_900, 100, gene_id="W")], 10_000)
    assert n_long == 0 and models[0].blocks == ((9_900, 10_000), (1, 100))
    assert models[0].size_bp == 201


def test_label_name_falls_back_to_gene_id_token():
    (model,), _ = gt.build_gene_models([row("gene", 1, 10, gene_id="HUM_HERP6B_U14")], 100)
    assert model.name == "U14"


def test_top_n_by_coverage_skips_uncovered():
    assert gt.select_labelled([0, 5, 9, 0, 7], top_n=2) == [2, 4]
    assert gt.select_labelled([0, 0]) == []


def test_mean_block_depth_binned_matches_unbinned():
    depth = [1.0, 1.0, 3.0, 3.0]  # bin_size 2 over 8 bp
    assert gt.mean_block_depth(depth, 2, [(1, 8)]) == pytest.approx(2.0)
    assert gt.mean_block_depth(depth, 2, [(4, 5)]) == pytest.approx(2.0)

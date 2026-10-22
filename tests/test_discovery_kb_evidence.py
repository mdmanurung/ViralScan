"""DISC-13: the route (i) ``source=kb`` evidence producer, on hand-computed synthetic dirs."""

import gzip
import json
from pathlib import Path

import pytest

from viralscan import discovery
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

TAXONOMY = [  # taxid, parent, rank, name
    (1, 1, "R", "root"),
    (5, 1, "G", "genus A"),
    (10, 5, "S", "virus A"),
    (11, 10, "S1", "virus A strain"),
    (20, 1, "S", "virus B"),
    (30, 1, "S", "virus C"),
]


def _gene(gene_id: str, taxid: str, viral: bool = True) -> GeneIdentity:
    return GeneIdentity(gene_id=gene_id, genome_accession="", status="x", viral=viral, taxid=taxid)


def _setup(tmp_path: Path, molecules, space: str = "X", candidates: str = "10\n20\n30\n"):
    disc, sample = tmp_path / "disc", tmp_path / "sample"
    (sample / "results").mkdir(parents=True)
    disc.mkdir()
    with (disc / discovery.TAXONOMY_TSV).open("w") as h:
        h.write("taxid\tparent_taxid\trank\tname\n")
        h.writelines(f"{t}\t{p}\t{r}\t{n}\n" for t, p, r, n in TAXONOMY)
    (disc / discovery.CANDIDATES_TSV).write_text("taxid\n" + candidates)
    header = "\t".join(discovery.UMI_GROUP_COLUMNS)
    groups = [  # cb, umi, taxid, class
        ("C1", "U1", 10, "resolved"),
        ("C2", "U9", 10, "resolved"),  # called, not a kb molecule
        ("C3", "U6", 10, "resolved"),  # kb molecule but C3 is not called
        ("C1", "U2", 10, "ambiguous"),  # not a resolved group
        ("C2", "U5", 20, "resolved"),
    ]
    with gzip.open(disc / discovery.UMI_GROUPS_TSV, "wt") as h:
        h.write(header + "\n")
        for cb, umi, tax, cls in groups:
            h.write(f"{cb}\t{umi}\t1\t1\t0\t{tax}\t{tax}\t{cls}\t1\t\t0\n")
    VirusIdentityTable(
        (
            _gene("gA1", "10"),
            _gene("gA2", "11"),  # descendant of candidate 10
            _gene("gB", "20"),
            _gene("gG", "5"),  # genus-level gene, an ancestor of candidate 10
            _gene("gH", "", viral=False),
        )
    ).write_tsv(sample / "results" / "virus_identity.tsv")
    (sample / "results" / "called_cells.tsv").write_text("barcode\nC1-1\nC2-1\n")
    (sample / "results" / "cell_calling_input.json").write_text(json.dumps({"count_layer": space}))
    if molecules is not None:
        with gzip.open(sample / "results" / "molecule_assignments.tsv.gz", "wt") as h:
            h.write("barcode\tumi\tec_ids\tstatus\tgenes\tweights\thost_virus_boundary\n")
            for cb, umi, genes in molecules:
                h.write(f"{cb}\t{umi}\t0\tx\t{genes}\t\tfalse\n")
    return disc, sample


MOLECULES = [
    ("C1", "U1", "gA1"),  # unique T10: best
    ("C1", "U2", "gA1;gA2"),  # two genes, both T10: best
    ("C1", "U3", "gA1;gB"),  # other virus ties: adjudicated for 10 and 20, best for neither
    ("C2", "U4", "gA2;gH"),  # host ties: adjudicated for 10 only
    ("C2", "U5", "gB"),  # unique T20: best
    ("C3", "U6", "gA1"),  # T10 best, uncalled cell
    ("C2", "U7", "gH"),  # host only: ignored
    ("C1", "U8", ""),  # unresolved: ignored
]


def test_counts_match_hand_oracle_and_roundtrip(tmp_path):
    disc, sample = _setup(tmp_path, MOLECULES)
    out = discovery.build_kb_evidence(disc, sample)
    assert out == disc / discovery.KB_EVIDENCE_TSV
    rows = discovery.read_evidence_table(out)
    # taxid 30 has no panel gene: no row (not_assessed), never a negative row.
    assert set(rows) == {(10, "corroboration"), (20, "corroboration")}
    t10, t20 = rows[(10, "corroboration")], rows[(20, "corroboration")]
    assert t10["source"] == "kb"
    # T10: adjudicated U1 U2 U3 U4 U6 = 5; best U1 U2 U6 = 3; shared with called groups: (C1, U1) only.
    assert [t10[c] for c in discovery.EVIDENCE_COLUMNS[3:]] == [5, 3, 1, 1]
    # T20: adjudicated U3 U5 = 2; best U5 = 1; shared (C2, U5).
    assert [t20[c] for c in discovery.EVIDENCE_COLUMNS[3:]] == [2, 1, 1, 1]


def test_panel_taxon_without_molecules_is_negative_row(tmp_path):
    disc, sample = _setup(tmp_path, [("C1", "U1", "gA1")])
    rows = discovery.read_evidence_table(discovery.build_kb_evidence(disc, sample))
    assert [rows[(20, "corroboration")][c] for c in discovery.EVIDENCE_COLUMNS[3:]] == [0, 0, 0, 0]


def test_refuses_starsolo_barcode_space(tmp_path):
    disc, sample = _setup(tmp_path, MOLECULES, space="GeneFull.raw")
    with pytest.raises(discovery.DiscoveryError, match="barcode space"):
        discovery.build_kb_evidence(disc, sample)


def test_refuses_when_molecule_barcodes_miss_called_cells(tmp_path):
    disc, sample = _setup(tmp_path, [("ZZ", "U1", "gA1")])
    with pytest.raises(discovery.DiscoveryError, match="barcode spaces differ"):
        discovery.build_kb_evidence(disc, sample)


def test_refuses_missing_molecule_assignments(tmp_path):
    disc, sample = _setup(tmp_path, MOLECULES)
    (sample / "results" / "molecule_assignments.tsv.gz").unlink()
    with pytest.raises(discovery.DiscoveryError, match="molecule_assignments"):
        discovery.build_kb_evidence(disc, sample)


def _counts(rows, taxid):
    return [rows[(taxid, "corroboration")][c] for c in discovery.EVIDENCE_COLUMNS[3:]]


def test_empty_molecule_file_raises(tmp_path):
    disc, sample = _setup(tmp_path, [])
    with pytest.raises(discovery.DiscoveryError, match="no molecules"):
        discovery.build_kb_evidence(disc, sample)


def test_no_molecule_in_called_cell_raises_even_without_panel_hits(tmp_path):
    disc, sample = _setup(tmp_path, [("ZZ", "U1", "gH"), ("ZZ", "U2", "")])
    with pytest.raises(discovery.DiscoveryError, match="barcode spaces differ"):
        discovery.build_kb_evidence(disc, sample)


def test_duplicate_cb_umi_rows_count_once(tmp_path):
    mols = [("C1", "U1", "gA1"), ("C1", "U1", "gA1"), ("C1", "U1", "gA1;gB"), ("C2", "U5", "gB")]
    disc, sample = _setup(tmp_path, mols)
    rows = discovery.read_evidence_table(discovery.build_kb_evidence(disc, sample))
    # T10: one distinct (C1, U1), best via its gA1-only row; the gA1;gB row does not remove it.
    assert _counts(rows, 10) == [1, 1, 1, 1]
    assert _counts(rows, 20)[:2] == [2, 1]


def test_ancestor_above_species_competes_but_never_wins(tmp_path):
    mols = [("C1", "U1", "gG"), ("C1", "U2", "gA1;gG"), ("C1", "U3", "gA2"), ("C2", "U4", "gB")]
    disc, sample = _setup(tmp_path, mols, candidates="5\n10\n11\n20\n")
    rows = discovery.read_evidence_table(discovery.build_kb_evidence(disc, sample))
    # Species 10: U1 U2 U3 adjudicated; the genus gene blocks U1 and U2, only U3 (strain gene) is best.
    assert _counts(rows, 10)[:2] == [3, 1]
    # Genus 5: every viral gene is the genus or a descendant; U1 U2 U3 best.
    assert _counts(rows, 5)[:2] == [3, 3]
    # Strain 11: gA1 (species ancestor) is allowed, genus gG is not: U3 is best, U2 is not.
    assert _counts(rows, 11)[:2] == [3, 1]
    # Sibling species B never gets virus A molecules as taxon-best.
    assert _counts(rows, 20)[:2] == [1, 1]


def test_species_ancestor_gene_is_best_for_strain_candidate(tmp_path):
    disc, sample = _setup(tmp_path, [("C1", "U1", "gA1")], candidates="11\n")
    rows = discovery.read_evidence_table(discovery.build_kb_evidence(disc, sample))
    assert _counts(rows, 11)[:2] == [1, 1]


def test_reads_the_real_molecule_assignments_writer(tmp_path):
    from viralscan.multimapping import _write_molecule_assignments

    disc, sample = _setup(tmp_path, None)
    names = ["gA1", "gB", "gH"]
    # (barcode, umi, ec, count); ec 0 -> gA1, 1 -> gA1+gB tie, 2 -> host only; BUS order is sorted.
    records = [("C1", "U1", 0, 1), ("C1", "U2", 1, 2), ("C2", "U3", 2, 1), ("C2", "U5", 1, 1)]
    _write_molecule_assignments(
        sample / "results" / "molecule_assignments.tsv.gz",
        records,
        {"C1": 0, "C2": 1},
        {0: [0], 1: [0, 1], 2: [2]},
        names,
        {0, 1},
        lambda cell, genes: [1.0 / len(genes)] * len(genes),
    )
    rows = discovery.read_evidence_table(discovery.build_kb_evidence(disc, sample))
    assert _counts(rows, 10)[:2] == [3, 1]  # U1 best; U2 and U5 tie with gB
    assert _counts(rows, 20)[:2] == [2, 0]

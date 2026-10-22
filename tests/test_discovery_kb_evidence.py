"""DISC-13: the route (i) ``source=kb`` evidence producer, on hand-computed synthetic dirs."""

import gzip
import json
from pathlib import Path

import pytest

from viralscan import discovery
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

TAXONOMY = [  # taxid, parent, rank, name
    (1, 1, "R", "root"),
    (10, 1, "S", "virus A"),
    (11, 10, "S1", "virus A strain"),
    (20, 1, "S", "virus B"),
    (30, 1, "S", "virus C"),
]


def _gene(gene_id: str, taxid: str, viral: bool = True) -> GeneIdentity:
    return GeneIdentity(gene_id=gene_id, genome_accession="", status="x", viral=viral, taxid=taxid)


def _setup(tmp_path: Path, molecules, space: str = "X"):
    disc, sample = tmp_path / "disc", tmp_path / "sample"
    (sample / "results").mkdir(parents=True)
    disc.mkdir()
    with (disc / discovery.TAXONOMY_TSV).open("w") as h:
        h.write("taxid\tparent_taxid\trank\tname\n")
        h.writelines(f"{t}\t{p}\t{r}\t{n}\n" for t, p, r, n in TAXONOMY)
    (disc / discovery.CANDIDATES_TSV).write_text("taxid\n10\n20\n30\n")
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
            _gene("gH", "", viral=False),
        )
    ).write_tsv(sample / "results" / "virus_identity.tsv")
    (sample / "results" / "called_cells.tsv").write_text("barcode\nC1-1\nC2-1\n")
    (sample / "results" / "cell_calling_input.json").write_text(json.dumps({"count_layer": space}))
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

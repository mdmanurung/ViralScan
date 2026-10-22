"""DISC-13: the route (i) ``source=star`` evidence producer, stubbed BAM records."""

import pytest

from tests.test_discovery_kb_evidence import _counts, _setup
from viralscan import discovery
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable


def _gene(acc: str, taxid: str) -> GeneIdentity:
    return GeneIdentity(
        gene_id=f"g{acc}", genome_accession=acc, status="x", viral=True, taxid=taxid
    )


def _run(tmp_path, records, *, bam=True):
    disc, sample = _setup(tmp_path, None)
    VirusIdentityTable(
        (_gene("AccA1", "10"), _gene("AccA2", "11"), _gene("AccB", "20"), _gene("AccG", "5"))
    ).write_tsv(sample / "results" / "virus_identity.tsv")
    if bam:
        (sample / "star_hostvirus").mkdir()
        (sample / discovery.STAR_VIRAL_BAM).write_bytes(b"")
    return disc, sample, lambda _bam: iter(records)


def V(name, acc, nh, cb, ub):
    return (name, f"VIRUS|{acc}", nh, cb, ub)


RECORDS = [
    V("r1", "AccA1", 1, "C1-1", "U1"),  # unique T10: best
    V("r2", "AccA1", 2, "C1-1", "U2"),  # both placements visible, both T10 (lineage): best
    V("r2", "AccA2", 2, "C1-1", "U2"),
    V("r3", "AccA1", 2, "C1-1", "U3"),  # NH 2, one visible T record: adjudicated only
    V("r4", "AccA1", 2, "C1-1", "U4"),  # other virus: not best for either
    V("r4", "AccB", 2, "C1-1", "U4"),
    V("r5", "AccB", 1, "C2-1", "U5"),  # unique T20: best
    V("r6", "AccA1", 1, "C1-1", "-"),  # no UB: not a molecule
    V("r7", "AccA1", 1, "-", "U7"),  # no CB: not a molecule
    V("r8", "AccG", 1, "C1-1", "U8"),  # genus ancestor: competes, never wins
    V("r9", "AccA1", 1, "C3-1", "U6"),  # best, but C3 not called
]


def test_counts_match_hand_oracle_and_roundtrip(tmp_path):
    disc, sample, reader = _run(tmp_path, RECORDS)
    out = discovery.build_star_evidence(disc, sample, record_reader=reader)
    assert out == disc / discovery.STAR_EVIDENCE_TSV
    rows = discovery.read_evidence_table(out)
    # taxid 30 has no accession: no row.
    assert set(rows) == {(10, "corroboration"), (20, "corroboration")}
    assert rows[(10, "corroboration")]["source"] == "star"
    # T10: adjudicated U1 U2 U3 U4 U8 U6 = 6; best U1 U2 U6 = 3; shared (C1, U1) only.
    assert _counts(rows, 10) == [6, 3, 1, 1]
    # T20: adjudicated U4 U5; best U5; shared (C2, U5).
    assert _counts(rows, 20) == [2, 1, 1, 1]


def test_corroborated_when_counts_pass(tmp_path):
    recs = [r for r in RECORDS if r[0] in ("r1", "r2", "r9")] + [V("r10", "AccA2", 1, "C1", "U10")]
    disc, sample, reader = _run(tmp_path, recs)
    rows = discovery.read_evidence_table(
        discovery.build_star_evidence(disc, sample, record_reader=reader)
    )
    assert _counts(rows, 10) == [4, 4, 1, 1]
    cand = {"taxid": "10", "name": "A", "rank": "S", "host_lineage": "0", "viral_only_db": "0"}
    assert discovery.assign_state(cand, rows)["state"] == "corroborated"


def test_accession_without_molecule_is_zero_row(tmp_path):
    disc, sample, reader = _run(tmp_path, [V("r1", "AccA1", 1, "C1", "U1")])
    rows = discovery.read_evidence_table(
        discovery.build_star_evidence(disc, sample, record_reader=reader)
    )
    assert _counts(rows, 20) == [0, 0, 0, 0]


def test_missing_bam_raises(tmp_path):
    disc, sample, reader = _run(tmp_path, RECORDS, bam=False)
    with pytest.raises(discovery.DiscoveryError, match="viral.bam"):
        discovery.build_star_evidence(disc, sample, record_reader=reader)


def test_no_called_cell_overlap_raises(tmp_path):
    disc, sample, reader = _run(tmp_path, [V("r1", "AccA1", 1, "ZZ-1", "U1")])
    with pytest.raises(discovery.DiscoveryError, match="barcode spaces differ"):
        discovery.build_star_evidence(disc, sample, record_reader=reader)


def test_no_molecules_raises(tmp_path):
    disc, sample, reader = _run(tmp_path, [V("r1", "AccA1", 1, "C1", "-")])
    with pytest.raises(discovery.DiscoveryError, match="no CB/UB molecule"):
        discovery.build_star_evidence(disc, sample, record_reader=reader)

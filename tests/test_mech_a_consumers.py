"""MECH-A step 4: the consumers read the Virus Identity table, not gene-ID strings.

Phase 1 (count-affecting core): multimap's viral/host partition and detection's
grouping, naming, sibling note and EVE flag. The retired name-keyed rules
(``SIBLING_VIRUS_PAIRS``, ``EVE_RISK_GENERA``) are frozen below so the
catalogue-wide regression proves the table-driven logic reproduces them.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from viralscan import virus_catalog
from viralscan.virus_grouping import (
    group_genes_by_identity,
    legacy_eve_risk,
    load_run_identity,
    virus_facts,
)
from viralscan.virus_identity import (
    ANELLOVIRIDAE,
    GeneIdentity,
    VirusIdentityTable,
    build_identity_table,
)

# Frozen copies of the rules this change retired (constants.py kept them only as
# the no-table fallback).
RETIRED_SIBLING_PAIRS = {
    "Human herpesvirus 6": "Human herpesvirus 6b",
    "Human herpesvirus 6b": "Human herpesvirus 6",
    "Human herpesvirus 1": "Human herpesvirus 2",
    "Human herpesvirus 2": "Human herpesvirus 1",
}
RETIRED_EVE_RISK_GENERA = frozenset(
    {
        "Alphatorquevirus",
        "Betatorquevirus",
        "Gammatorquevirus",
        "Samektorquevirus",
        "Memtorquevirus",
        "Hetorquevirus",
        "Gyrovirus",
        "Anelloviridae",
        "Torque teno virus",
    }
)


@pytest.fixture(scope="module")
def catalogue_table(tmp_path_factory) -> VirusIdentityTable:
    """An identity table over every catalogued genome: two genes per accession."""
    rows = virus_catalog.load_catalogue()
    t2g = tmp_path_factory.mktemp("t2g") / "t2g.txt"
    with open(t2g, "w") as fh:
        for row in rows:
            acc = row["accession_version"] or row["accession"]
            for n in (1, 2):
                fh.write(f"{acc}_gene{n}\t{acc}_gene{n}\t\t\t{acc}\n")
    return build_identity_table(t2g, [])


def test_every_retired_sibling_pair_is_a_table_sibling_group(catalogue_table) -> None:
    facts = virus_facts(catalogue_table)
    for a, b in RETIRED_SIBLING_PAIRS.items():
        assert a in facts and b in facts, (a, b)
        assert facts[a].sibling_group, a
        assert facts[a].sibling_group == facts[b].sibling_group, (a, b)


def test_ebv_types_are_siblings_not_one_virus(catalogue_table) -> None:
    facts = virus_facts(catalogue_table)
    ebv1 = facts["Epstein-Barr virus"]
    ebv2 = next(f for f in facts.values() if f.virus_key == "taxid:12509")
    assert ebv1.virus_key != ebv2.virus_key
    assert ebv1.sibling_group == ebv2.sibling_group == "HHV-4"


def test_every_retired_eve_genus_is_table_artifact_risk(catalogue_table) -> None:
    """Every genome the retired substring test flagged now carries artifact_risk, not eve_risk."""
    facts = virus_facts(catalogue_table)
    by_gene = catalogue_table.by_gene()
    checked_genera = set()
    for row in virus_catalog.load_catalogue():
        acc = row["accession_version"] or row["accession"]
        gene = by_gene[f"{acc}_gene1"]
        retired_flag = (
            row["genus"] in RETIRED_EVE_RISK_GENERA
            or row["family"] in RETIRED_EVE_RISK_GENERA
            or any(g in (row["species"] or "") for g in RETIRED_EVE_RISK_GENERA)
        )
        if retired_flag:
            # ANELLO-PRIOR.4 (2026-10-03): the retired EVE genera (all Anelloviridae) now
            # carry artifact_risk=low_complexity and no eve_risk.
            f = facts[gene.virus_name]
            assert f.artifact_risk == "low_complexity" and not f.eve_risk, (acc, gene.virus_name)
            checked_genera.add(row["genus"] or row["family"])
    # Every retired genus (or its family) is represented in the catalogue.
    for genus in RETIRED_EVE_RISK_GENERA - {"Torque teno virus"}:
        assert (
            genus in checked_genera
            or genus == ANELLOVIRIDAE
            or any(r["genus"] == genus for r in virus_catalog.load_catalogue())
        ), genus


def test_non_eve_viruses_are_not_flagged(catalogue_table) -> None:
    facts = virus_facts(catalogue_table)
    for name in ("Epstein-Barr virus", "Human herpesvirus 1", "Human herpesvirus 6b"):
        assert not facts[name].eve_risk, name
        assert facts[name].artifact_risk == "", name


def _table(*genes: GeneIdentity) -> VirusIdentityTable:
    return VirusIdentityTable(tuple(genes))


def test_uncatalogued_virus_fails_closed_to_the_genus_test() -> None:
    """No risk_class to read: an anellovirus-genus name keeps artifact_risk, never silently empty."""
    table = _table(
        GeneIdentity("g1", "X1", "uncatalogued", True, "accession:X1", "Betatorquevirus"),
        GeneIdentity("g2", "X2", "uncatalogued", True, "accession:X2", "Some other virus"),
    )
    facts = virus_facts(table)
    assert facts["Betatorquevirus"].artifact_risk == "low_complexity"
    assert facts["Betatorquevirus"].eve_risk is False
    assert facts["Some other virus"].artifact_risk == ""
    assert facts["Some other virus"].eve_risk is False
    assert not legacy_eve_risk("Betatorquevirus")


def test_one_eve_gene_flags_the_whole_virus() -> None:
    table = _table(
        GeneIdentity("g1", "A", "catalogued", True, "taxid:1", "V", risk_class=""),
        GeneIdentity("g2", "B", "catalogued", True, "taxid:1", "V", risk_class="eve"),
    )
    assert virus_facts(table)["V"].eve_risk is True


def test_one_low_complexity_gene_labels_the_whole_virus() -> None:
    table = _table(
        GeneIdentity("g1", "A", "catalogued", True, "taxid:1", "V", risk_class=""),
        GeneIdentity("g2", "B", "catalogued", True, "taxid:1", "V", risk_class="low_complexity"),
    )
    facts = virus_facts(table)["V"]
    assert facts.artifact_risk == "low_complexity" and facts.eve_risk is False


def test_gene_to_group_conservation(catalogue_table) -> None:
    """Group sums equal gene sums; no gene is in two groups (incl. genus:<Genus> keys)."""
    rng = np.random.default_rng(7)
    genes = catalogue_table.viral_gene_ids()
    counts = dict(zip(genes, rng.integers(0, 50, size=len(genes)).tolist()))
    groups, detected = group_genes_by_identity(genes, catalogue_table)
    members = [g for gl in groups.values() for g in gl]
    assert len(members) == len(set(members)) == len(genes)
    assert sum(sum(counts[g] for g in gl) for gl in groups.values()) == sum(counts.values())
    # names map 1:1 to keys
    keys = catalogue_table.groups()
    assert len(groups) == len(keys)
    assert any(k.startswith("genus:") for k in keys), "anello genus keys present"
    by_gene = catalogue_table.by_gene()
    for name, gl in groups.items():
        assert {by_gene[g].virus_name for g in gl} == {name}
        assert len({by_gene[g].virus_key for g in gl}) == 1
    assert detected == set(groups)


def test_grouping_skips_host_and_unknown_genes() -> None:
    table = _table(
        GeneIdentity("v1", "A", "catalogued", True, "taxid:1", "V"),
        GeneIdentity("h1", "", "host", False),
    )
    groups, detected = group_genes_by_identity(["v1", "h1", "absent"], table)
    assert groups == {"V": ["v1"]} and detected == {"V"}


def test_load_run_identity_roundtrip_and_absence(tmp_path: Path) -> None:
    assert load_run_identity(tmp_path) is None
    table = _table(GeneIdentity("v1", "A", "catalogued", True, "taxid:1", "V"))
    table.write_tsv(tmp_path / "results" / "virus_identity.tsv")
    assert load_run_identity(tmp_path).viral_gene_ids() == ["v1"]


def test_unparseable_table_raises(tmp_path: Path) -> None:
    path = tmp_path / "results" / "virus_identity.tsv"
    path.parent.mkdir()
    with open(path, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerow(["gene_id"])
    with pytest.raises(ValueError):
        load_run_identity(tmp_path)


# --------------------------------------------------------------------- multimap


def test_multimap_partition_reads_the_table_not_analysis_txt(tmp_path: Path) -> None:
    from viralscan.scripts.multimap import viral_gene_partition

    (tmp_path / "log").mkdir()
    (tmp_path / "log" / "analysis.txt").write_text("hostgene\nv1\n")
    _table(
        GeneIdentity("v1", "A", "catalogued", True, "taxid:1", "V"),
        GeneIdentity("hostgene", "", "host", False),
    ).write_tsv(tmp_path / "results" / "virus_identity.tsv")
    assert viral_gene_partition(str(tmp_path)) == {"v1"}


def test_multimap_partition_legacy_fallback(tmp_path: Path) -> None:
    from viralscan.scripts.multimap import viral_gene_partition

    (tmp_path / "log").mkdir()
    (tmp_path / "log" / "analysis.txt").write_text("v1\nv2\n")
    assert viral_gene_partition(str(tmp_path)) == {"v1", "v2"}
    assert viral_gene_partition(str(tmp_path / "missing")) == set()


# -------------------------------------------------------------------- detection


def test_sibling_note_uses_table_groups() -> None:
    from viralscan.scripts.detection import check_sibling_crossmapping

    stats = {
        "Human herpesvirus 6b": {"viral_molecules_total_est": 5000},
        "Human herpesvirus 6": {"viral_molecules_total_est": 10},
        "Other": {"viral_molecules_total_est": 1},
    }
    groups = {"Human herpesvirus 6b": "HHV-6", "Human herpesvirus 6": "HHV-6"}
    notes = check_sibling_crossmapping(stats, groups)
    assert set(notes) == {"Human herpesvirus 6"}
    assert "possible_em_bleed: 500:1" in notes["Human herpesvirus 6"]


def test_sibling_note_needs_no_name_pair() -> None:
    """EBV-2 is a sibling through the table even though no name pair lists it."""
    from viralscan.scripts.detection import check_sibling_crossmapping

    stats = {
        "Epstein-Barr virus": {"viral_molecules_total_est": 10000},
        "Human gammaherpesvirus 4 type 2": {"viral_molecules_total_est": 100},
    }
    groups = {k: "HHV-4" for k in stats}
    assert "Human gammaherpesvirus 4 type 2" in check_sibling_crossmapping(stats, groups)
    assert check_sibling_crossmapping(stats) == {}  # legacy pairs know nothing of it


def test_write_tsv_outputs_eve_risk_comes_from_facts(tmp_path: Path) -> None:
    import pandas as pd

    from viralscan.scripts.detection import write_tsv_outputs
    from viralscan.virus_grouping import VirusFacts

    stats = {
        "Risky": {
            "viral_molecules_total_est": 5,
            "infected_cells": 1,
            "total_cells": 10,
            "pct_infected": 10.0,
            "viral_molecules_per_10k_est": 1.0,
        }
    }
    facts = {"Risky": VirusFacts("genus:X", "Risky", "", True)}
    write_tsv_outputs(stats, pd.DataFrame(), str(tmp_path), facts=facts)
    with open(tmp_path / "results" / "viral_summary.tsv") as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    assert row["eve_risk"] == "True"

"""MECH-A step 4, phase 2: umap, hostresponse, evidence and gene programmes read the table."""

from __future__ import annotations

from pathlib import Path

import pytest

from viralscan import virus_catalog
from viralscan.evidence import resolve_viral_target
from viralscan.gene_programs import load_catalogue, resolve_markers
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable, build_identity_table


@pytest.fixture(scope="module")
def table(tmp_path_factory) -> VirusIdentityTable:
    """Identity table over every catalogued genome (two genes per accession)."""
    t2g = tmp_path_factory.mktemp("t2g") / "t2g.txt"
    with open(t2g, "w") as fh:
        for row in virus_catalog.load_catalogue():
            acc = row["accession_version"] or row["accession"]
            for n in (1, 2):
                fh.write(f"{acc}_g{n}\t{acc}_g{n}\t\t\t{acc}\n")
    return build_identity_table(t2g, [])


def _genes(table: VirusIdentityTable) -> list[str]:
    return table.viral_gene_ids()


# --------------------------------------------------------------------- evidence


@pytest.mark.parametrize(
    "selector, name",
    [
        ("taxid:10376", "Epstein-Barr virus"),
        ("10376", "Epstein-Barr virus"),
        ("ebv", "Epstein-Barr virus"),
        ("EBV", "Epstein-Barr virus"),
        ("Epstein-Barr virus", "Epstein-Barr virus"),
        ("Human gammaherpesvirus 4", "Epstein-Barr virus"),
        ("hhv6b", "Human herpesvirus 6b"),
        ("hhv6a", "Human herpesvirus 6"),
        ("Human herpesvirus 6b", "Human herpesvirus 6b"),
        ("hsv1", "Human herpesvirus 1"),
        ("kshv", "Human herpesvirus 8"),
        ("10298", "Human herpesvirus 1"),
    ],
)
def test_selector_resolves_through_the_table(table, selector, name) -> None:
    label, genes = resolve_viral_target(selector, _genes(table), identity=table)
    assert label == name
    by_gene = table.by_gene()
    assert genes and {by_gene[g].virus_name for g in genes} == {name}


def test_ebv_type2_is_its_own_selector(table) -> None:
    label, genes = resolve_viral_target("taxid:12509", _genes(table), identity=table)
    assert label == "Human gammaherpesvirus 4 type 2"
    ebv1 = resolve_viral_target("ebv", _genes(table), identity=table)[1]
    assert not set(genes) & set(ebv1)


def test_ttv_selects_the_whole_family(table) -> None:
    label, genes = resolve_viral_target("ttv", _genes(table), identity=table)
    assert label == "Anelloviridae"
    by_gene = table.by_gene()
    assert {by_gene[g].family for g in genes} == {"Anelloviridae"}
    assert len({by_gene[g].virus_key for g in genes}) > 1


def test_genus_key_selects_one_genus(table) -> None:
    label, genes = resolve_viral_target("genus:Betatorquevirus", _genes(table), identity=table)
    assert label == "Betatorquevirus"
    assert {table.by_gene()[g].virus_key for g in genes} == {"genus:Betatorquevirus"}


def test_unknown_selector_and_empty_selector_raise(table) -> None:
    with pytest.raises(ValueError, match="No exact viral target"):
        resolve_viral_target("not-a-virus", _genes(table), identity=table)
    with pytest.raises(ValueError, match="exact --virus selector"):
        resolve_viral_target("", _genes(table), identity=table)


def test_ambiguous_selector_refuses() -> None:
    t = VirusIdentityTable(
        (
            GeneIdentity("a", "A1", "catalogued", True, "taxid:1", "V1", organism="Same"),
            GeneIdentity("b", "B1", "catalogued", True, "taxid:2", "V2", organism="Same"),
        )
    )
    with pytest.raises(ValueError, match="ambiguous"):
        resolve_viral_target("same", ["a", "b"], identity=t)


def test_exact_gene_id_still_wins(table) -> None:
    gene = _genes(table)[0]
    assert resolve_viral_target(gene, _genes(table), identity=table)[0] == gene


# -------------------------------------------------------------- gene programmes


def test_every_programme_virus_is_one_table_virus(table) -> None:
    """The programme catalogue's virus names each map to exactly one virus_key."""
    names = table.virus_names()
    key_of_name = {name: key for key, name in names.items()}
    assert len(key_of_name) == len(names), "display names are unique across keys"
    programme_viruses = {row["virus"] for row in load_catalogue()}
    missing = programme_viruses - set(key_of_name)
    assert not missing, f"programme catalogue names absent from the identity table: {missing}"


def test_resolve_markers_by_virus_key_equals_by_name(table) -> None:
    catalogue = load_catalogue()
    var_names = [row["gene_id_bundled"] for row in catalogue] + ["x"]
    key_of_name = {name: key for key, name in table.virus_names().items()}
    for name in sorted({row["virus"] for row in catalogue}):
        by_name, _ = resolve_markers(name, var_names, catalogue=catalogue)
        by_key, _ = resolve_markers(
            key_of_name[name], var_names, catalogue=catalogue, identity=table
        )
        by_name_ident, _ = resolve_markers(name, var_names, catalogue=catalogue, identity=table)
        assert [m.var_name for m in by_key] == [m.var_name for m in by_name], name
        assert [m.var_name for m in by_name_ident] == [m.var_name for m in by_name], name


# ------------------------------------------------------------------ umap / host


def test_umap_names_come_from_the_table() -> None:
    from viralscan.scripts.umap import _gene_to_virus

    t = VirusIdentityTable(
        (
            GeneIdentity("TTVgp1", "NC_002076.2", "catalogued", True, "genus:X", "Genus X"),
            GeneIdentity("h", "", "host", False),
        )
    )
    # The legacy prefix rule named this gene by prefix; the table names it by identity.
    assert _gene_to_virus(["TTVgp1"], t) == {"TTVgp1": "Genus X"}
    assert _gene_to_virus(["TTVgp1"], None) == {"TTVgp1": "Torque teno virus"}


def test_hostresponse_viral_set_from_table_and_legacy(tmp_path: Path) -> None:
    from viralscan.scripts.hostresponse import _load_viral_accessions

    t = VirusIdentityTable(
        (
            GeneIdentity("v1", "A", "catalogued", True, "taxid:1", "V"),
            GeneIdentity("ENSG1", "", "host", False),
        )
    )
    path = t.write_tsv(tmp_path / "virus_identity.tsv")
    assert _load_viral_accessions(str(path)) == {"v1"}
    legacy = tmp_path / "analysis.txt"
    legacy.write_text("v1\nENSG1\n")
    assert _load_viral_accessions(str(legacy)) == {"v1", "ENSG1"}

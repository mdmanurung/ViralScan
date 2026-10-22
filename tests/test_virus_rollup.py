"""PLAN TAXO-01: NCBI lineage per virus and the species/genus aggregate-up counts."""

from collections import Counter

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan.multimapping import viral_group_matrices
from viralscan.scripts.multimap import stamp_viral_groups
from viralscan.scripts.ncbi_fetch import fetch_taxonomy, parse_taxonomy_xml
from viralscan.virus_identity import (
    COLUMNS,
    TAXONOMY_COLUMNS,
    GeneIdentity,
    VirusIdentityTable,
    _ncbi_lineage,
)
from viralscan.virus_rollup import assignment_columns, group_summaries

TAX_XML = """<?xml version="1.0" ?><TaxaSet>
<Taxon><TaxId>12509</TaxId><ScientificName>Human herpesvirus 4 type 2</ScientificName><Rank>no rank</Rank>
<AkaTaxIds><TaxId>99999</TaxId></AkaTaxIds>
<LineageEx>
<Taxon><TaxId>3044472</TaxId><ScientificName>Orthoherpesviridae</ScientificName><Rank>family</Rank></Taxon>
<Taxon><TaxId>10374</TaxId><ScientificName>Gammaherpesvirinae</ScientificName><Rank>subfamily</Rank></Taxon>
<Taxon><TaxId>10375</TaxId><ScientificName>Lymphocryptovirus</ScientificName><Rank>genus</Rank></Taxon>
<Taxon><TaxId>3050299</TaxId><ScientificName>Lymphocryptovirus humangamma4</ScientificName><Rank>species</Rank></Taxon>
</LineageEx></Taxon>
<Taxon><TaxId>3050299</TaxId><ScientificName>Lymphocryptovirus humangamma4</ScientificName><Rank>species</Rank>
<LineageEx><Taxon><TaxId>10375</TaxId><ScientificName>Lymphocryptovirus</ScientificName><Rank>genus</Rank></Taxon>
</LineageEx></Taxon>
</TaxaSet>"""


def test_parse_taxonomy_xml_lineage_alias_and_self_species():
    t = parse_taxonomy_xml(TAX_XML)
    assert t["12509"]["species_taxid"] == "3050299"
    assert (t["12509"]["genus"], t["12509"]["subfamily"], t["12509"]["family"]) == (
        "Lymphocryptovirus",
        "Gammaherpesvirinae",
        "Orthoherpesviridae",
    )
    assert t["99999"] is t["12509"]  # merged taxid alias
    assert (
        t["3050299"]["species"] == "Lymphocryptovirus humangamma4"
    )  # a species is its own species


def test_packaged_lineage_unites_ebv_types_but_not_binned_species():
    ebv, ebv2 = _ncbi_lineage("10376"), _ncbi_lineage("12509")
    assert ebv["species_taxid"] == ebv2["species_taxid"] == "3050299"
    assert ebv["family_ncbi"] == ebv2["family_ncbi"] == "Orthoherpesviridae"
    binned = _ncbi_lineage("2049444")  # "Gammapapillomavirus sp." is a bin, not a species
    assert binned["species_taxid"] == "" and binned["genus_ncbi"] == "Gammapapillomavirus"
    assert _ncbi_lineage("") == {} and _ncbi_lineage("1") == {}


def _gene(gid, key, name, taxid="", genus="", **kw):
    return GeneIdentity(gid, "", "catalogued", True, key, name, taxid=taxid, genus=genus, **kw)


def test_group_keys_species_genus_and_anellovirus():
    ebv = _gene("e1", "taxid:10376", "EBV", "10376", "Lymphocryptovirus", **_ncbi_lineage("10376"))
    assert ebv.group_key("virus") == "taxid:10376"
    assert ebv.group_key("species") == "species:taxid:3050299"
    assert ebv.group_name("species") == "Lymphocryptovirus humangamma4"
    assert ebv.group_key("genus") == "genus:Lymphocryptovirus"
    anello = _gene(
        "a1",
        "genus:Alphatorquevirus",
        "Alphatorquevirus",
        genus="Alphatorquevirus",
        **_ncbi_lineage("2809145"),
    )
    assert anello.group_key("species") == "genus:Alphatorquevirus"  # ANDET-05 key, never the bin
    assert anello.group_key("genus") == "genus:Alphatorquevirus"
    bare = _gene("u1", "accession:XY1", "XY1")
    assert {bare.group_key(lvl) for lvl in ("virus", "species", "genus")} == {"accession:XY1"}


def test_old_identity_table_without_taxonomy_columns_still_reads(tmp_path):
    p = tmp_path / "virus_identity.tsv"
    cols = [c for c in COLUMNS if c not in TAXONOMY_COLUMNS]
    pd.DataFrame([{c: "" for c in cols} | {"gene_id": "g", "viral": "true"}], columns=cols).to_csv(
        p, sep="\t", index=False
    )
    (g,) = VirusIdentityTable.read_tsv(p).genes
    assert g.viral and g.species_taxid == ""


def test_viral_group_matrices_unique_compatible_host():
    # genes 0,1 -> group 0 (EBV), 2 -> group 1 (EBV type 2); gene 3 is host
    keys = Counter({(0, (0,)): 2, (0, (0, 1)): 1, (1, (1, 2)): 3, (1, (2, 3)): 4})
    u, c, h = viral_group_matrices(keys, {0: 0, 1: 0, 2: 1}, n_cells=2, n_groups=2)
    assert u.toarray().tolist() == [[3, 0], [0, 0]]  # (0,1) is two genes of one virus: still unique
    assert c.toarray().tolist() == [[3, 0], [3, 7]]
    assert h.toarray().tolist() == [[0, 0], [0, 4]]


def _toy_run():
    """Two EBV types (one species) and a host gene; X = host-conservative-like equal split."""
    lin1, lin2 = _ncbi_lineage("10376"), _ncbi_lineage("12509")
    table = VirusIdentityTable(
        (
            _gene("EBV1_a", "taxid:10376", "EBV", "10376", "Lymphocryptovirus", **lin1),
            _gene("EBV2_a", "taxid:12509", "EBV type 2", "12509", "Lymphocryptovirus", **lin2),
            GeneIdentity("HOST", "", "host", False),
        )
    )
    genes = ["EBV1_a", "EBV2_a", "HOST"]
    # cell 0: 10 unique EBV1, 90 shared EBV1/EBV2; cell 1: 1 unique EBV2, 5 shared with host
    keys = Counter({(0, (0,)): 10, (0, (0, 1)): 90, (1, (1,)): 1, (1, (1, 2)): 5})
    X = np.array([[10 + 45, 45, 0], [0, 1, 5]], dtype=float)
    adata = ad.AnnData(X=sparse.csr_matrix(X), var=pd.DataFrame(index=genes))
    adata.uns["multimap_method"] = "host-conservative"

    class Layers:
        viral_molecule_keys = keys

    stamp_viral_groups(adata, genes, table, Layers)
    return adata, table


def test_group_summaries_and_assignment_note(tmp_path):
    adata, table = _toy_run()
    adata.write_h5ad(tmp_path / "m.h5ad")  # detection reads the stamped groups back from the H5AD
    frames = group_summaries(ad.read_h5ad(tmp_path / "m.h5ad"), table, called_mask=[True, True])
    v = frames["virus"].set_index("group_key")
    assert v.loc["taxid:12509", "molecules_unique"] == 1
    assert v.loc["taxid:12509", "molecules_compatible"] == 96  # 1 + 90 shared + 5 host-shared
    assert v.loc["taxid:12509", "molecules_host_shared"] == 5
    assert v.loc["taxid:12509", "molecules_package_total"] == 46
    sp = frames["species"].iloc[0]
    assert sp["members"] == "EBV|EBV type 2" and sp["n_members"] == 2
    assert sp["molecules_unique"] == 101  # every EBV1/EBV2 molecule is unique to the species
    assert sp["molecules_compatible"] == 106
    assert sp["molecules_unique"] <= sp["molecules_package_total"] <= sp["molecules_compatible"]
    assert sp["n_called_cells_unique"] == 2
    cols = assignment_columns(frames, table)
    assert cols["EBV type 2"]["unique_fraction"] == pytest.approx(1 / 96, abs=1e-6)
    assert "species Lymphocryptovirus humangamma4" in cols["EBV type 2"]["aggregate_note"]
    assert cols["EBV"]["aggregate_note"] == ""  # 10/100 unique clears 0.05


def test_no_identity_no_groups():
    adata, _ = _toy_run()
    assert group_summaries(adata, None) == {}
    assert assignment_columns({}, None) == {}


@pytest.mark.network
def test_fetch_taxonomy_live():
    t = fetch_taxonomy(["12509"])
    assert t["12509"]["species_taxid"] == "3050299"

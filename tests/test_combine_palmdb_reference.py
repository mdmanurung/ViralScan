import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/combine_palmdb_reference.py"
spec = importlib.util.spec_from_file_location("combine_palmdb_reference", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

CAT_HEADER = "taxid\tspecies\torganism\tlineage\n"


def inputs(tmp_path, palm_ids=("u1", "u2", "u3"), ann_ids=("NC_1.1__YP_1.1", "NC_2.1__YP_2.1")):
    ann, palm = tmp_path / "ann", tmp_path / "palm"
    ann.mkdir()
    palm.mkdir()
    (ann / "proteins.faa").write_text("".join(f">{i}\nMKV{n}\n" for n, i in enumerate(ann_ids)))
    (ann / "proteins.tsv").write_text(
        "target\tprotein_id\tvirus_id\tspecies\tfamily\tgene\tproduct\tlength_aa\n"
        + "".join(
            f"{i}\tP{n}\ttaxon_{n + 10}\tVirus {n}\tF\tg\tprod\t4\n" for n, i in enumerate(ann_ids)
        )
    )
    (ann / "viruses.t2g").write_text(
        "".join(f"{i}\ttaxon_{n + 10}\tVirus {n}\n" for n, i in enumerate(ann_ids))
    )
    (palm / "palmdb_rdrp_seqs.fa").write_text(
        "".join(f">{i}\nAAAA{n}\n" for n, i in enumerate(palm_ids))
    )
    # u1 is the representative of {u1, u3}; u2 is unclassified
    (palm / "palmdb_clustered_t2g.txt").write_text("u1\tu1\nu3\tu1\nu2\tu2\n")
    (palm / "ID_to_taxonomy_mapping.csv").write_text(
        "ID,rep_ID,phylum,class,order,family,genus,species,strandedness\n"
        "u1,u1,P,C,O,Coronaviridae,Betacoronavirus,Zika virus,+ssRNA\n"
        "u3,u1,P,C,O,Coronaviridae,Betacoronavirus,Zika virus,+ssRNA\n"
        "u2,u2,.,.,.,.,.,.,.\n"
    )
    cat = tmp_path / "cat.tsv"
    cat.write_text(
        CAT_HEADER + "64320\tZika virus\tZika virus\tViruses; Flavivirus\n10\tVirus 0\tVirus 0\t\n"
    )
    return ann, palm, cat


def run(tmp_path, ann, palm, cat, **kw):
    return m.build(ann, palm, cat, tmp_path / "out", verify_palmdb=False, **kw)


def test_combined_files_and_catalogue_match(tmp_path):
    ann, palm, cat = inputs(tmp_path)
    receipt = run(tmp_path, ann, palm, cat)
    out = tmp_path / "out"
    assert (out / "combined_proteins.faa").read_text().count(">") == 5
    t2g = [line.split("\t") for line in (out / "combined.t2g").read_text().splitlines()]
    assert ("u3", "palmdb_u1", "Zika virus") in {tuple(t) for t in t2g}
    assert ("u2", "palmdb_u2", "palmdb_u2") in {
        tuple(t) for t in t2g
    }  # unclassified: gene id label
    assert ("NC_1.1__YP_1.1", "taxon_10", "Virus 0") in {tuple(t) for t in t2g}
    tax = {r["gene"]: r for r in m.read_tsv(out / "palmdb_taxonomy.tsv")}
    assert tax["palmdb_u1"]["n_sequences"] == "2"
    assert tax["palmdb_u1"]["catalogue_taxids"] == "64320"
    assert tax["palmdb_u2"]["catalogue_match"] == "none"
    prot = {r["target"]: r for r in m.read_tsv(out / "proteins.tsv")}
    assert prot["u3"]["virus_id"] == "palmdb_u1" and prot["NC_1.1__YP_1.1"]["product"] == "prod"
    assert receipt["overlap"]["palmdb"]["clusters"] == 2
    assert receipt["overlap"]["palmdb"]["clusters_without_species"] == 1


def test_colliding_id_is_prefixed(tmp_path):
    ann, palm, cat = inputs(tmp_path, ann_ids=("u1", "x"))
    receipt = run(tmp_path, ann, palm, cat)
    assert receipt["overlap"]["palmdb"]["id_collisions_prefixed"] == 1
    assert ">palmdb_u1\n" in (tmp_path / "out/combined_proteins.faa").read_text()
    assert "palmdb_u1\tpalmdb_u1\t" in (tmp_path / "out/combined.t2g").read_text()


def test_duplicate_ids_fail(tmp_path):
    ann, palm, cat = inputs(tmp_path, palm_ids=("u1", "u2", "u2"))
    with pytest.raises(ValueError, match="Duplicate sequence ID u2"):
        run(tmp_path, ann, palm, cat)
    (tmp_path / "b").mkdir()
    ann, palm, cat = inputs(tmp_path / "b", ann_ids=("a", "a"))
    with pytest.raises(ValueError, match="Duplicate sequence ID a"):
        run(tmp_path / "b", ann, palm, cat)


def test_checksum_pin_and_existing_output(tmp_path):
    ann, palm, cat = inputs(tmp_path)
    with pytest.raises(ValueError, match="pinned"):
        m.build(ann, palm, cat, tmp_path / "out")
    (tmp_path / "out").mkdir()
    (tmp_path / "out/x").write_text("")
    with pytest.raises(ValueError, match="already exists"):
        run(tmp_path, ann, palm, cat)


def test_blast_duplicates_and_expected(tmp_path):
    ann, palm, cat = inputs(tmp_path)
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "u1\tNC_1.1__YP_1.1\t95.0\t4\t4\nu3\tNC_1.1__YP_1.1\t80.0\t4\t4\nu2\tNC_2.1__YP_2.1\t99\t1\t4\n"
    )
    expected = tmp_path / "exp.tsv"
    expected.write_text(
        "name\ttaxon_ids\tname_pattern\tpalmdb_pattern\nZIKV\t\t^Zika virus$\t\nVirus0\t10\t\t\n"
    )
    receipt = run(tmp_path, ann, palm, cat, blast_hits=hits, expected=expected)
    dup = receipt["overlap"]["blast_duplicates"]
    assert dup["annotated_proteins_with_a_palmdb_match"] == 1
    assert dup["palmdb_clusters_with_an_annotated_match"] == 1
    rows = {r["name"]: r for r in receipt["expected_virus_presence"]}
    assert rows["ZIKV"]["catalogue_taxids"] == "64320"
    assert rows["ZIKV"]["palmdb_clusters_by_name"] == 1
    assert rows["Virus0"]["annotated_taxids_with_proteins"] == "10"

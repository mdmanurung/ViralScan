"""AA-08..AA-10: translated search ported into the package (synthetic data, no network, no kb runs)."""

import argparse
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan import translated
from viralscan.translated import TranslatedError

# ── summarize / barcodes (from the repro tests) ───────────────────────────────


def toy_matrix(barcodes):
    # taxon_1: 3 UMIs in cell 0, 12 in cell 1; taxon_2: 5 in cell 1
    return ad.AnnData(
        X=sp.csr_matrix(np.array([[3, 0], [12, 5], [0, 0]], dtype=np.float32)),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=["taxon_1", "taxon_2"]),
    )


PROTEINS = {
    "g__p1": ("taxon_1", "Virus one", "prot one", ""),
    "g__p2": ("taxon_1", "Virus one", "", "gene2"),
    "g__q1": ("taxon_2", "Virus two", "prot q", ""),
}
EC = ["0\t0", "1\t1", "2\t2", "3\t0,2"]
BUS = [
    ["AAAA", "u1", "0", "1"],
    ["CCCC", "u1", "0", "1"],
    ["CCCC", "u2", "1", "1"],
    ["CCCC", "u3", "2", "1"],
    ["GGGG", "u4", "3", "1"],
]


def run_summary(cells=None, barcodes=("AAAA-1", "CCCC-1", "GGGG-1")):
    return translated.summarize(
        toy_matrix(list(barcodes)), ["g__p1", "g__p2", "g__q1"], EC, BUS, PROTEINS,
        translated.taxon_ids("1"), cells,
    )  # fmt: skip


def test_summarize_all_barcodes():
    result = run_summary()
    exp = result["expected"]["taxon_1"]
    assert (exp["umis"], exp["barcodes"], exp["barcodes_ge10_umis"]) == (15, 2, 1)
    assert exp["distinct_proteins"] == 2 and exp["label"] == "Virus one"
    assert result["total_virus_umis"] == 20 and result["barcodes_any_virus"] == 2
    assert result["bus_records"] == 5 and result["ambiguous_multi_virus_records"] == 1
    assert result["top_viruses"][0][0] == "taxon_1"


def test_summarize_restricted_to_cells_ignores_suffix():
    result = run_summary({"CCCC"})
    assert result["restricted_to_cells"] and result["barcodes_in_matrix"] == 1
    assert result["expected"]["taxon_1"]["umis"] == 12
    assert result["expected"]["taxon_1"]["barcodes"] == 1
    assert result["bus_records"] == 3 and result["ambiguous_multi_virus_records"] == 0


def test_taxon_ids_accepts_lists_and_prefixes():
    assert translated.taxon_ids("37296;taxon_10376") == ["taxon_37296", "taxon_10376"]
    assert translated.taxon_ids("") == []


def test_check_barcodes_reports_suffix_difference():
    observed = ["AACC-1", "GGTT-1", "TTAA-1", "CCGG-1"]
    report = translated.compare(observed, ["AACC", "GGTT", "TTAA", "ACGT"])
    assert report["raw_overlap"] == 0 and report["normalized_overlap"] == 3
    assert report["format_differs"] and report["frozen_found_fraction"] == 0.75
    assert report["observed"]["with_suffix"] == 4 and report["frozen"]["with_suffix"] == 0
    assert translated.compare(["AACC"], ["GGTT"])["normalized_overlap"] == 0
    assert translated.compare(["AAAC"], ["GTTT"])["revcomp_overlap"] == 1


# ── combine_palmdb (from the repro tests) ─────────────────────────────────────


def _inputs(tmp_path, palm_ids=("u1", "u2", "u3"), ann_ids=("NC_1.1__YP_1.1", "NC_2.1__YP_2.1")):
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
    (palm / "palmdb_clustered_t2g.txt").write_text("u1\tu1\nu3\tu1\nu2\tu2\n")
    (palm / "ID_to_taxonomy_mapping.csv").write_text(
        "ID,rep_ID,phylum,class,order,family,genus,species,strandedness\n"
        "u1,u1,P,C,O,Coronaviridae,Betacoronavirus,Zika virus,+ssRNA\n"
        "u3,u1,P,C,O,Coronaviridae,Betacoronavirus,Zika virus,+ssRNA\n"
        "u2,u2,.,.,.,.,.,.,.\n"
    )
    return ann, palm


def test_combine_palmdb_files_and_labels(tmp_path):
    ann, palm = _inputs(tmp_path)
    receipt = translated.combine_palmdb(ann, palm, tmp_path / "out")
    out = tmp_path / "out"
    assert (out / "combined_proteins.faa").read_text().count(">") == 5
    t2g = {tuple(line.split("\t")) for line in (out / "combined.t2g").read_text().splitlines()}
    assert ("u3", "palmdb_u1", "Zika virus") in t2g
    assert ("u2", "palmdb_u2", "palmdb_u2") in t2g  # unclassified: gene id label
    assert ("NC_1.1__YP_1.1", "taxon_10", "Virus 0") in t2g
    prot = {r["target"]: r for r in translated._read_tsv(out / "proteins.tsv")}
    assert prot["u3"]["virus_id"] == "palmdb_u1" and prot["NC_1.1__YP_1.1"]["product"] == "prod"
    assert receipt["overlap"]["palmdb"]["clusters"] == 2
    assert receipt["overlap"]["palmdb"]["clusters_without_species"] == 1
    assert (out / "reference.json").is_file()


def test_combine_palmdb_prefixes_collisions_and_rejects_duplicates(tmp_path):
    ann, palm = _inputs(tmp_path, ann_ids=("u1", "x"))
    receipt = translated.combine_palmdb(ann, palm, tmp_path / "out")
    assert receipt["overlap"]["palmdb"]["id_collisions_prefixed"] == 1
    assert ">palmdb_u1\n" in (tmp_path / "out/combined_proteins.faa").read_text()
    (tmp_path / "dup").mkdir()
    ann2, palm2 = _inputs(tmp_path / "dup", palm_ids=("u1", "u2", "u2"))
    with pytest.raises(TranslatedError, match="Duplicate sequence ID"):
        translated.combine_palmdb(ann2, palm2, tmp_path / "dup" / "out")


# ── kb --aa wrappers ──────────────────────────────────────────────────────────


def test_translated_index_refuses_old_kallisto(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "viralscan.validation.kb_tools",
        lambda: {"kallisto": ("0.50.0", "/k"), "bustools": ("0.45.1", "/b")},
    )
    with pytest.raises(TranslatedError, match="0.51.1"):
        translated.build_translated_index(tmp_path / "p.faa", tmp_path / "h.fa", tmp_path / "i.idx")


def test_translated_index_command_shape(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "viralscan.validation.kb_tools",
        lambda: {"kallisto": ("0.52.0", "/k/kallisto"), "bustools": ("0.45.1", "/b")},
    )
    monkeypatch.setattr(
        translated.shutil, "which", lambda name: "/bin/kb" if name == "kb" else None
    )
    seen = {}
    monkeypatch.setattr(subprocess, "run", lambda cmd, check: seen.setdefault("cmd", cmd))
    cmd = translated.build_translated_index(
        tmp_path / "p.faa", tmp_path / "h.fa", tmp_path / "i.idx", 3
    )
    assert cmd[:2] == ["/bin/kb", "ref"] and "--aa" in cmd and "--d-list" in cmd
    assert cmd[cmd.index("--d-list") + 1] == str(tmp_path / "h.fa") and "-k" in cmd
    assert "--strand" not in cmd and seen["cmd"] is cmd


def test_count_translated_command_shape(tmp_path, monkeypatch):
    monkeypatch.setattr(
        translated.shutil, "which", lambda name: "/bin/kb" if name == "kb" else None
    )
    monkeypatch.setattr(subprocess, "run", lambda cmd, check: None)
    cmd = translated.count_translated(
        "a_1.fq", "a_2.fq", "10xv3", None, tmp_path / "i", tmp_path / "t", tmp_path / "o", 2
    )
    assert (
        cmd[1:5] == ["count", "--aa", "--parity", "single"] and "-w" not in cmd
    )  # 10xv3 has an on-list
    cmd = translated.count_translated(
        "a_1.fq", "a_2.fq", "dropseq", None, tmp_path / "i", tmp_path / "t", tmp_path / "o", 2
    )
    assert cmd[cmd.index("-w") + 1] == "None"  # no on-list: kb must not knee-call (SW-21)


# ── the verb on an existing counts dir ────────────────────────────────────────


def _counts_dir(tmp_path, barcodes=("AAAA-1", "CCCC-1", "GGGG-1")) -> Path:
    d = tmp_path / "counts"
    (d / "counts_unfiltered").mkdir(parents=True)
    toy_matrix(list(barcodes)).write_h5ad(d / "counts_unfiltered" / "adata.h5ad")
    (d / "transcripts.txt").write_text("g__p1\ng__p2\ng__q1\n")
    (d / "matrix.ec").write_text("\n".join(EC) + "\n")
    (d / "output.unfiltered.bus").write_bytes(b"")
    return d


def _index_dir(tmp_path) -> Path:
    idx = tmp_path / "index"
    idx.mkdir()
    (idx / "proteins.tsv").write_text(
        "target\tvirus_id\tspecies\tproduct\tgene\n"
        + "".join(f"{t}\t{v}\t{s}\t{p}\t{g}\n" for t, (v, s, p, g) in PROTEINS.items())
    )
    (idx / "translated_index.json").write_text(json.dumps({"schema_version": 1}))
    return idx


def test_translated_verb_counts_dir_publishes_posthoc(tmp_path, monkeypatch):
    run = tmp_path / "run"
    (run / "S1" / "results").mkdir(parents=True)
    (run / "S1" / "results" / "called_cells.tsv").write_text("barcode\nAAAA-1\nCCCC-1\n")
    monkeypatch.setattr(
        translated,
        "_sample_config",
        lambda r, s: (
            run / "S1",
            SimpleNamespace(sample1="", sample2="", technology="10xv3", whitelist=None),
        ),
    )
    monkeypatch.setattr(translated, "bus_text_rows", lambda bus: iter(BUS))
    monkeypatch.setattr(translated, "kallisto_translated_ok", lambda: (True, "0.52.0"))
    args = argparse.Namespace(
        run_dir=run, sample=None, index_dir=_index_dir(tmp_path), counts_dir=_counts_dir(tmp_path),
        cells=None, expected_taxa="1", output=None, threads=1,
    )  # fmt: skip
    out = translated.run_translated(args)
    assert out == run / "posthoc" / "translated" / "S1"
    names = sorted(p.name for p in out.iterdir())
    assert names == [
        "barcode_check.json",
        "summary.all.json",
        "summary.cells.json",
        "translated_complete.json",
        "translated_receipt.json",
    ]
    cells = json.loads((out / "summary.cells.json").read_text())
    assert cells["restricted_to_cells"] and cells["barcodes_in_matrix"] == 2
    assert cells["expected"]["taxon_1"]["umis"] == 15
    receipt = json.loads((out / "translated_receipt.json").read_text())
    assert receipt["command"] is None and receipt["exploratory"] is True
    published = json.loads((out / "translated_complete.json").read_text())["outputs"]
    assert published == [n for n in names if n != "translated_complete.json"]  # marker written last


def test_translated_verb_aborts_on_barcode_format_mismatch(tmp_path, monkeypatch):
    run = tmp_path / "run"
    (run / "S1" / "results").mkdir(parents=True)
    (run / "S1" / "results" / "called_cells.tsv").write_text(
        "barcode\nTTTT\nGGCC\n"
    )  # none observed
    monkeypatch.setattr(
        translated,
        "_sample_config",
        lambda r, s: (
            run / "S1",
            SimpleNamespace(sample1="", sample2="", technology="10xv3", whitelist=None),
        ),
    )
    args = argparse.Namespace(
        run_dir=run, sample=None, index_dir=_index_dir(tmp_path), counts_dir=_counts_dir(tmp_path),
        cells=None, expected_taxa=None, output=None, threads=1,
    )  # fmt: skip
    with pytest.raises(TranslatedError, match="barcode formats differ"):
        translated.run_translated(args)
    assert (run / "posthoc" / "translated" / "S1" / "barcode_check.json").is_file()

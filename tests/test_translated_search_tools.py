import gzip
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis" / "translated_search"))
import check_barcodes  # noqa: E402
import make_configs  # noqa: E402
import summarize_translated  # noqa: E402


def fastq(path: Path, length: int) -> Path:
    with gzip.open(path, "wt") as handle:
        for i in range(5):
            handle.write(f"@r{i}\n{'A' * length}\n+\n{'I' * length}\n")
    return path


def write_chem(
    chem_dir: Path, sample: str, tech: str, r1_length: int, strand: str, basis="bundled on-list"
):
    (chem_dir / f"{sample}.json").write_text(
        json.dumps(
            {
                "chemistry": {"chemistry": tech, "basis": basis, "r1_length": r1_length},
                "strand_pilot": {"choice": strand},
            }
        )
    )


def manifest_row(dataset, sample, r1, r2, whitelist="-"):
    return "\t".join([dataset, sample, "S", str(r1), str(r2), str(whitelist), "1", "hint"])


@pytest.fixture
def library(tmp_path):
    chem = tmp_path / "chem"
    chem.mkdir()
    whitelist = tmp_path / "wl.txt"
    whitelist.write_text("AAAA\n")
    rows = []
    for sample, dataset, tech, bc, strand, wl in [
        ("a", "ebv", "10xv2", 26, "forward", "-"),
        ("b", "kshv_gse1", "10xv3", 28, "forward", "-"),
        ("c", "covid", "10xv3", 28, "reverse", whitelist),
        # chem JSON (10xv2) overrides the manifest hint (10xv3)
        ("d", "other", "10xv2", 28, "forward", "-"),
    ]:
        r1, r2 = (
            fastq(tmp_path / f"{sample}_1.fq.gz", bc),
            fastq(tmp_path / f"{sample}_2.fq.gz", 90),
        )
        write_chem(
            chem, sample, tech, bc, strand, "user on-list" if wl != "-" else "bundled on-list"
        )
        rows.append(manifest_row(dataset, sample, r1, r2, wl))
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "dataset\tsample\tgrp\tr1\tr2\twhitelist\tpairs\thint\n" + "\n".join(rows) + "\n"
    )
    base = tmp_path / "base.yaml"
    base.write_text(yaml.safe_dump({"kb": "kb", "work_dir": "w", "whitelist": "stale"}))
    taxa = tmp_path / "taxa.tsv"
    taxa.write_text("dataset\ttaxon_ids\tnote\nebv\t10376\t\nkshv\t37296\t\n")
    return tmp_path, manifest, chem, base, taxa


def test_make_configs_groups_whitelist_and_expected(library):
    tmp_path, manifest, chem, base, taxa = library
    out = tmp_path / "out"
    summary = make_configs.make_configs(manifest, chem, base, out, taxa)
    assert summary["excluded"] == {}
    assert sorted(summary["groups"]) == ["10xv2_forward", "10xv3_forward", "10xv3_reverse_wl-wl"]
    assert summary["groups"]["10xv2_forward"] == ["a", "d"]
    cfg = yaml.safe_load((out / "10xv3_reverse_wl-wl.yaml").read_text())
    assert (
        cfg["technology"] == "10xv3"
        and cfg["whitelist"].endswith("wl.txt")
        and list(cfg["samples"]) == ["c"]
    )
    # a whitelist key in the base config never leaks into a group without a whitelist
    assert "whitelist" not in yaml.safe_load((out / "10xv2_forward.yaml").read_text())
    expected = dict(
        line.split("\t") for line in (out / "expected_virus.tsv").read_text().splitlines()[1:]
    )
    assert expected == {"a": "10376", "b": "37296", "c": "", "d": ""}


def test_make_configs_swap_and_excludes_bad_pairs(library):
    tmp_path, manifest, chem, base, taxa = library
    # manifest lists the cDNA file as r1 for sample a: excluded unless swapped
    text = manifest.read_text().splitlines()
    cols = text[1].split("\t")
    cols[3], cols[4] = cols[4], cols[3]
    text[1] = "\t".join(cols)
    manifest.write_text("\n".join(text) + "\n")
    bad = make_configs.make_configs(manifest, chem, base, tmp_path / "o1", taxa)
    assert "R1 modal length 90" in bad["excluded"]["a"]
    fixed = make_configs.make_configs(manifest, chem, base, tmp_path / "o2", taxa, swap=("a",))
    assert fixed["excluded"] == {}
    assert fixed["samples"]["a"]["corrections"] == ["r1/r2 swapped relative to manifest"]


def test_make_configs_requires_whitelist_for_user_onlist(library):
    tmp_path, manifest, chem, base, taxa = library
    write_chem(chem, "b", "10xv3", 28, "forward", "user on-list")
    summary = make_configs.make_configs(manifest, chem, base, tmp_path / "o", taxa)
    assert "no whitelist" in summary["excluded"]["b"]


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
    return summarize_translated.summarize(
        toy_matrix(list(barcodes)),
        ["g__p1", "g__p2", "g__q1"],
        EC,
        BUS,
        PROTEINS,
        summarize_translated.taxon_ids("1"),
        cells,
    )


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
    assert (
        result["expected"]["taxon_1"]["umis"] == 12
        and result["expected"]["taxon_1"]["barcodes"] == 1
    )
    assert result["bus_records"] == 3 and result["ambiguous_multi_virus_records"] == 0


def test_taxon_ids_accepts_lists_and_prefixes():
    assert summarize_translated.taxon_ids("37296;taxon_10376") == ["taxon_37296", "taxon_10376"]
    assert summarize_translated.taxon_ids("") == []


def test_check_barcodes_reports_suffix_difference():
    observed = ["AACC-1", "GGTT-1", "TTAA-1", "CCGG-1"]
    report = check_barcodes.compare(observed, ["AACC", "GGTT", "TTAA", "ACGT"])
    assert report["raw_overlap"] == 0 and report["normalized_overlap"] == 3
    assert report["format_differs"] and report["frozen_found_fraction"] == 0.75
    assert report["observed"]["with_suffix"] == 4 and report["frozen"]["with_suffix"] == 0
    assert check_barcodes.compare(["AACC"], ["GGTT"])["normalized_overlap"] == 0
    assert check_barcodes.compare(["AAAC"], ["GTTT"])["revcomp_overlap"] == 1


def test_check_barcodes_cli_abort_threshold(tmp_path, monkeypatch, capsys):
    frozen = tmp_path / "frozen.txt"
    frozen.write_text("\n".join(f"{b}-1" for b in ["AACC", "GGTT", "TTAA", "CCGG", "ACGT"]) + "\n")
    adata = toy_matrix(["AACC", "GGTT", "TTAA"])
    adata.write_h5ad(tmp_path / "obs.h5ad")

    def run(*extra):
        monkeypatch.setattr(
            sys, "argv", ["check_barcodes", str(tmp_path / "obs.h5ad"), str(frozen), *extra]
        )
        return check_barcodes.main()

    assert run() == 2  # 3 of 5 frozen cells found: 60% < 90%
    assert "ABORT" in capsys.readouterr().err
    assert run("--abort-below", "0.5") == 0
    assert run("--abort-below", "0") == 0

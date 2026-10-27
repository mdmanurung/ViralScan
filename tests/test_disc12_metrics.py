"""DISC-12 metrics script (scripts/disc12_metrics.py) on a tiny synthetic run root."""

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

anndata = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")

SPEC = importlib.util.spec_from_file_location(
    "disc12_metrics", Path(__file__).resolve().parents[1] / "scripts" / "disc12_metrics.py"
)
m = importlib.util.module_from_spec(SPEC)
sys.modules["disc12_metrics"] = m
SPEC.loader.exec_module(m)

# root 1 > viruses 10239 > genus 100 > panel species 200 (strain 300) / sibling 201; host 9606 under 2759
TAXONOMY = [(1, 1, "R"), (10239, 1, "D"), (100, 10239, "G"), (200, 100, "S"), (300, 200, "S1"),
            (201, 100, "S"), (2759, 1, "D"), (9606, 2759, "S")]  # fmt: skip
HEADER = ["taxid", "name", "rank", "groups_called", "groups_offcell", "n_cells", "reads",
          "distinct_minimizers", "minimizer_per_read", "host_relation", "db_has_9606", "state",
          "low_minimizer_diversity", "offcell_dominant", "host_lineage", "viral_only_db"]  # fmt: skip


def write_tsv(path, header, rows):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def arm(run, cid, cands, cells):
    """cands: {taxid: (called, off)}; cells: {barcode: {taxid: groups}}."""
    d = run / "out" / cid
    d.mkdir(parents=True)
    write_tsv(d / "taxonomy.tsv", ["taxid", "parent_taxid", "rank", "name", "depth", "host_relation"],
              [(t, p, r, f"t{t}", 0, "host" if t == 9606 else "other") for t, p, r in TAXONOMY])  # fmt: skip
    rows = [(t, f"t{t}", "S", c, o, 1, 1, 30, 1, "host" if t == 9606 else "other", 1, "observation",
             0, 0, int(t == 9606), 0) for t, (c, o) in cands.items()]  # fmt: skip
    write_tsv(d / "candidates.tsv", HEADER, rows)
    taxa = sorted({t for v in cells.values() for t in v})
    x = np.array([[cells[b].get(t, 0) for t in taxa] for b in cells], dtype=np.int32)
    anndata.AnnData(
        X=x.reshape(len(cells), len(taxa)),
        obs=pd.DataFrame(index=pd.Index(list(cells), name="barcode")),
        var=pd.DataFrame(index=pd.Index([str(t) for t in taxa], name="taxid")),
    ).write_h5ad(d / "cell_taxon.h5ad")
    (d / "discovery_complete.json").write_text("{}")
    (run / "receipts").mkdir(exist_ok=True)
    (run / "receipts" / f"{cid}.time.err").write_text(
        "\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02:03\n"
        "\tMaximum resident set size (kbytes): 2048\n\tFile system inputs: 10\n"
    )


def build(tmp_path, better_recovery=True, neg_fp_groups=2, fail=None):
    run = tmp_path / "run"
    run.mkdir(parents=True)
    conds = []
    for conf, hits, trim in ((0.05, 4, 1), (0.1, 4, 1), (0.05, 4, 0)):
        tag = f"c{conf}_h{hits}" + ("" if trim else "_notrim")
        for lib in ("POS", "NEG"):
            cid = f"{lib}__pluspf__{tag}"
            conds.append([cid, lib, "r", "c", "pluspf", "db", conf, hits, trim])
            if cid == fail:
                continue
            if lib == "NEG":
                arm(run, cid, {201: (neg_fp_groups, 0), 9606: (9, 1)}, {"N1": {201: 1}})
            else:
                boost = 2 if (conf == 0.1 and better_recovery) else 1
                arm(run, cid, {200: (3 * boost, 1), 300: (1, 0), 10239: (1, 0), 201: (1, 0)},
                    {"A": {200: 2 * boost}, "B": {300: 1}, "C": {201: 1}})  # fmt: skip
    write_tsv(run / "conditions.tsv", ["cond_id", "library", "run_root", "called_cells", "db_name",
              "db", "confidence", "min_hit_groups", "polya_trim"], conds)  # fmt: skip
    manifest = tmp_path / "manifest.tsv"
    write_tsv(manifest, ["role", "library"], [["positive_x", "POS"], ["negative_x", "NEG"]])
    truth = tmp_path / "truth.tsv"
    write_tsv(truth, ["barcode", "molecules"], [["A-1", 4], ["B-1", 2], ["C-1", 0]])
    out = tmp_path / "out"
    m.main(["--run", str(run), "--manifest", str(manifest), "--panel-taxid", "POS=200",
            "--truth", f"POS={truth}", "--out", str(out)])  # fmt: skip
    rows = {r["cond_id"]: r for r in m.tsv(out / "metrics.tsv")}
    return rows, json.loads((out / "decisions.json").read_text())


def test_metrics_and_rules(tmp_path):
    rows, dec = build(tmp_path)
    pos = rows["POS__pluspf__c0.05_h4"]
    assert pos["panel_reached"] == "species"
    assert (pos["rank_exact"], pos["rank_ancestor"], pos["rank_wrong_branch"]) == ("4", "1", "1")
    assert pos["n_fp_observation"] == "1"  # sibling 201; ancestor 10239 is on the truth lineage
    assert float(pos["truth_cell_fraction"]) == 1.0
    assert float(pos["median_cell_recovery"]) == 0.5  # A 2/4, B 1/2
    assert float(pos["wall_s"]) == 3723 and pos["bytes_read"] == "5120.0"
    neg = rows["NEG__pluspf__c0.05_h4"]
    assert neg["n_fp_observation"] == "1" and neg["max_eligible_groups_nontruth"] == "2"
    assert float(neg["host_leakage"]) == 10 / 12
    assert dec["rule2"]["pluspf"]["0.05/4"] == {"result": "k", "k": 3}
    assert dec["rule1"]["pluspf"]["result"] == "change"
    assert dec["rule1"]["pluspf"]["confidence"] == 0.1
    assert dec["rule3"]["result"] == "not_assessed"


def test_keep_when_no_recovery_gain_and_none_when_panel_below_k(tmp_path):
    _, dec = build(tmp_path, better_recovery=False)
    assert dec["rule1"]["pluspf"]["result"] == "keep"
    _, dec = build(tmp_path / "b", neg_fp_groups=5)
    assert dec["rule2"]["pluspf"]["0.05/4"]["result"] == "none"
    assert dec["rule1"]["pluspf"]["result"] == "not_assessed"


def test_failed_arm_is_not_assessed(tmp_path):
    rows, dec = build(tmp_path, fail="NEG__pluspf__c0.1_h4")
    assert rows["NEG__pluspf__c0.1_h4"]["status"] == "failed"
    assert dec["rule2"]["pluspf"]["0.1/4"]["result"] == "not_assessed"
    assert dec["rule1"]["pluspf"]["result"] == "not_assessed"


def test_positive_without_truth_exits(tmp_path):
    (tmp_path / "run").mkdir()
    manifest = tmp_path / "m.tsv"
    write_tsv(manifest, ["role", "library"], [["positive_x", "POS"]])
    with pytest.raises(SystemExit, match="missing for"):
        m.main(
            ["--run", str(tmp_path / "run"), "--manifest", str(manifest), "--out", str(tmp_path)]
        )

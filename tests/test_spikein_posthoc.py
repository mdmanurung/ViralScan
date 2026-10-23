"""PLAN SPIKE-03: rescore.py and recovery_ratio.py on synthetic bundles."""

import gzip
import importlib
import sys
from pathlib import Path

A = Path(__file__).parent.parent / "analysis" / "spikein_detection_limit"
sys.path.insert(0, str(A))
sys.path.insert(0, str(Path(__file__).parent))
rescore = importlib.import_module("rescore")
rr = importlib.import_module("recovery_ratio")
agg = importlib.import_module("aggregate")
score_condition = importlib.import_module("score_condition")
ta = importlib.import_module("test_spikein_aggregate")


def _full_bundle(root, cid, spiked, host_mols=50):
    """A scorable bundle: truth accession T.1, spike barcode ZZZ (2 truth molecules), plus host gene rows."""
    b = Path(root) / "cond" / cid / "bundle"
    (b / "star").mkdir(parents=True, exist_ok=True)
    (b / "star" / "star_contig_reads.tsv").write_text("contig\tn_any\tn_primary\n")
    (b / "star" / "star_spike_fates.tsv").write_text(
        "name\ton_truth\tfate\n" + ("SPIKE:a:0\t1\tsource\n" if spiked else "")
    )
    if spiked:
        (b / "truth.tsv").write_text("name\tcb\tclass\tcontigs\nSPIKE:a:0\tZZZ\tc\tVIRUS|T.1\n")
    for arm in ("old", "flank"):
        (b / arm).mkdir(exist_ok=True)
        (b / arm / "virus_identity.tsv").write_text(
            "gene_id\tgenome_accession\tviral\trole\nT_g\tT.1\ttrue\t\nENSG1\tENST1.1\tfalse\t\n"
        )
        with gzip.open(b / arm / "kb_molecules.tsv.gz", "wt") as fh:
            fh.write(f"barcode\tgene_id\tunique\tambiguous\nZZZ\tENSG1\t{host_mols}\t0\n")
            if spiked:
                fh.write("ZZZ\tT_g\t2\t0\n")
        (b / arm / "multimap_evidence.tsv").write_text(
            "gene_id\tviral_molecules_total_est\tevidence_tier\nT_g\t2\tcandidate_unique\n"
        )
    return b


def _conds(tmp_path):
    cf = tmp_path / "conditions.tsv"
    cf.write_text("cond_id\tbackground\tsource\tlevel\tseed\tlayout\tcells\tbg_rep\n"
                  "B-L003-s1-c005\tB\tS\t3\t1\tc005\t5\t1\nB-NULL-r1\tB\t\t0\t0\t\t0\t1\n")  # fmt: skip
    return ta.CONDS


def test_rescore_writes_v2_next_to_scores_and_refuses_without_done(tmp_path):
    conds = _conds(tmp_path)
    for c in conds:
        _full_bundle(tmp_path, c["cond_id"], int(c["level"]) > 0)
    (tmp_path / "cond" / "B-L003-s1-c005" / "DONE").write_text("{}")
    old = tmp_path / "cond" / "B-L003-s1-c005" / "bundle" / "scores.tsv"
    old.write_text("keep me\n")
    done, refused = rescore.rescore(tmp_path, conds)
    assert done == ["B-L003-s1-c005"] and len(refused) == 1 and "B-NULL-r1: no DONE" in refused[0]
    assert old.read_text() == "keep me\n"  # never overwritten
    v2 = score_condition.read_tsv(old.with_name("scores_v2.tsv"))
    arm = {r["arm"]: r for r in v2}["old"]
    assert arm["kb_total"] == "2" and arm["kb_total_viral_molecules"] == "2"  # host 50 not counted
    assert not (tmp_path / "cond" / "B-NULL-r1" / "bundle" / "scores_v2.tsv").exists()


def test_rescore_refuses_spiked_without_truth(tmp_path):
    conds = _conds(tmp_path)
    b = _full_bundle(tmp_path, "B-L003-s1-c005", True)
    (b.parent / "DONE").write_text("{}")
    (b / "truth.tsv").unlink()
    done, refused = rescore.rescore(tmp_path, conds[:1])
    assert done == [] and "truth.tsv" in refused[0]


def _row(bg, layout, level, seed, arm, raw, base, n=10):
    net = max(0.0, raw - base)
    return dict(background=bg, layout=layout, level=str(level), seed=str(seed), arm=arm, n_spiked=str(n),
                kb_total=str(net), kb_total_raw=raw, kb_total_in_baseline=base)  # fmt: skip


KEEP = {"per_background": {"B": {"verdict": "keep"}}}


def _scene(old, flank, vd=KEEP, base=0):
    rows = []
    for s, (o, f) in enumerate(zip(old, flank), 1):
        rows += [
            _row("B", "c100", 100, s, "old", o, base),
            _row("B", "c100", 100, s, "flank", f, base),
        ]
    return rr.ratio_rows(rows, vd)


def test_ratio_values_per_seed_range_and_flag():
    (r,) = _scene([10, 10, 10], [10, 5, 6])
    assert r["old_mean"] == 1.0 and abs(r["ratio"] - 0.7) < 1e-9
    assert (r["ratio_seed_min"], r["ratio_seed_max"]) == (0.5, 1.0)
    assert r["needs_user_decision"] is True and r["q15_verdict"] == "keep"


def test_ratio_not_flagged_when_ok_retire_or_low_old_recovery():
    assert _scene([10, 10, 10], [9, 8, 8])[0]["needs_user_decision"] is False  # 0.83 >= 0.80
    retire = {"per_background": {"B": {"verdict": "retire"}}}
    r = _scene([10, 10, 10], [1, 1, 1], retire)[0]
    assert r["ratio"] < rr.RATIO_MIN and r["needs_user_decision"] is False
    assert (
        _scene([0.6, 0.6, 0.6], [0, 0, 0])[0]["needs_user_decision"] is False
    )  # old recovery 0.06 < 0.1


def test_ratio_na_when_old_zero_or_floor_hit():
    r = _scene([0, 0, 0], [3, 3, 3])[0]
    assert (
        r["ratio"] is None
        and r["ratio_na_reason"] == "old_mean is 0"
        and not r["needs_user_decision"]
    )
    r = _scene([2, 10, 10], [10, 10, 10], base=3)[0]  # old seed 1: raw 2 < baseline 3 -> floored
    assert (
        r["ratio"] is None
        and r["ratio_na_reason"] == "net floor hit"
        and not r["needs_user_decision"]
    )


def test_run_refuses_then_writes_files_and_never_changes_verdict(tmp_path, capsys):
    conds = _conds(tmp_path)
    assert rr.run(tmp_path, tmp_path / "out", conds) == 2
    assert not (tmp_path / "out").exists()
    for c in conds:
        ta._write_cond(tmp_path, c)
    assert rr.run(tmp_path, tmp_path / "out", conds) == 0
    md = (tmp_path / "out" / "recovery_ratio.md").read_text()
    assert md.rstrip().endswith(rr.FOOTER) and "uncapped" in md
    head = (tmp_path / "out" / "recovery_ratio.tsv").read_text().splitlines()[0].split("\t")
    assert head == rr.COLS

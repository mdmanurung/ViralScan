import gzip
import importlib
import sys
from pathlib import Path

A = Path(__file__).parent.parent / "analysis" / "spikein_detection_limit"
sys.path.insert(0, str(A))
agg = importlib.import_module("aggregate")
grid = importlib.import_module("grid")
LV = grid.LEVELS


def srow(bg, arm, level, seed, layout, called, unsupported=0, rep=None, n=10):
    """One synthetic scores.tsv row (strings, as read from disk). level 0 = unspiked replicate."""
    r = dict(cond_id=f"{bg}-L{level}-s{seed}-{layout}", background=bg, arm=arm, level=str(level),
             seed=str(seed), layout=layout, bg_rep=str(rep or seed), n_spiked=str(n if level else 0),
             kb_total=str(n if called else 0), star_n_source_contig_any=str(n if called else 0),
             call_kb=str(bool(called)), call_star=str(bool(called)), call_combined=str(bool(called)),
             kb_unsupported_molecules=str(unsupported), n_false_calls="0",
             star_nonspike_reads_total="0", star_contigs_with_reads="0")  # fmt: skip
    return r


def grid_rows(bg, arm, layout, met):
    """met: {level: number of seeds (0-3) with the call true}."""
    return [srow(bg, arm, lv, s, layout, s <= met.get(lv, 0)) for lv in LV for s in (1, 2, 3)]


def lod_of(rows, arm="old", layout="c005", bg="B", call="call"):
    return agg.lod(agg.recovery_rows(rows), bg, arm, layout, call)


def lim(rows, **kw):
    """Comparable number; None = not reached or undetermined."""
    r = lod_of(rows, **kw)
    return None if r["status"] in ("not_reached", "undetermined") else r["num"]


def test_two_of_three_rule():
    assert lim(grid_rows("B", "old", "c005", {lv: 2 for lv in LV})) == 3
    assert lim(grid_rows("B", "old", "c005", {lv: 1 for lv in LV})) is None
    assert lim(grid_rows("B", "old", "c005", {3: 1, 10: 2, 30: 3, 100: 3, 300: 3})) == 10


def test_monotone_rule_and_missing_seed():
    # 10 met but 30 not: the limit is 100, not 10
    m = {3: 0, 10: 3, 30: 1, 100: 2, 300: 3}
    assert lim(grid_rows("B", "old", "c005", m)) == 100
    assert lim(grid_rows("B", "old", "c005", {300: 1})) is None  # top level not met
    rows = [
        r
        for r in grid_rows("B", "old", "c005", {lv: 3 for lv in LV})
        if not (r["level"] == "30" and r["seed"] != "1")
    ]
    # level 30 has one seed only: its call could still be met by the missing seeds -> undetermined, not 100
    assert lod_of(rows)["status"] == "undetermined"


def test_layouts_and_arms_are_separate():
    rows = grid_rows("B", "old", "c005", {lv: 3 for lv in LV}) + grid_rows(
        "B", "old", "c100", {300: 3}
    )
    assert lim(rows, layout="c005") == 3 and lim(rows, layout="c100") == 300


def test_recovery_fraction():
    r = srow("B", "old", 10, 1, "c005", True, n=8)
    r["kb_total"] = "6"
    assert agg.recovery_rows([r])[0]["recovery"] == 0.75


def test_fp_floor_mean_max():
    rows = [
        srow("B", "flank", 0, i, "", False, unsupported=u, rep=i)
        for i, u in ((1, 0), (2, 3), (3, 6))
    ]
    f = next(x for x in agg.fp_floor(rows) if x["metric"] == "kb_unsupported_molecules")
    assert (f["n_libraries"], f["mean"], f["max"]) == (3, 3.0, 6.0)


def _scene(
    old_met, flank_met, old_fp=(0, 0, 0), flank_fp=(0, 0, 0), layouts=("c005", "c100"), bg="B"
):
    rows = []
    for ly in layouts:
        rows += grid_rows(bg, "old", ly, old_met) + grid_rows(bg, "flank", ly, flank_met)
    for i in (1, 2, 3):
        rows.append(srow(bg, "old", 0, 0, "", False, old_fp[i - 1], rep=i))
        rows.append(srow(bg, "flank", 0, 0, "", False, flank_fp[i - 1], rep=i))
    lods = agg.all_lods(agg.recovery_rows(rows))
    return agg.verdict(lods, rows, [bg])


ALL3 = {lv: 3 for lv in LV}


def test_verdict_keep_tie_and_better():
    assert _scene(ALL3, ALL3)["overall"] == "keep"  # tie passes
    assert _scene({30: 3, 100: 3, 300: 3}, ALL3)["overall"] == "keep"  # flank lower limit


def test_verdict_retire_on_limit_or_floor():
    assert _scene(ALL3, {100: 3, 300: 3})["overall"] == "retire"
    assert _scene(ALL3, ALL3, flank_fp=(6, 0, 0))["overall"] == "retire"  # adds 6 > 5
    assert _scene(ALL3, ALL3, flank_fp=(5, 5, 5))["overall"] == "keep"  # adds exactly 5 passes
    assert (
        _scene(ALL3, ALL3, old_fp=(4, 4, 4), flank_fp=(9, 9, 9))["overall"] == "keep"
    )  # adds 5 over old


def test_verdict_old_not_reached_is_undetermined_never_keep():
    v = _scene({}, {})  # vacuous grid: nothing detected by either arm
    assert v["overall"] == "undetermined" and "uninformative" in v["per_background"]["B"]["notes"]
    assert _scene({}, ALL3)["overall"] == "undetermined"
    assert _scene(ALL3, {})["overall"] == "retire"


def test_verdict_missing_data_is_undetermined_and_retire_wins():
    assert _scene(ALL3, ALL3, layouts=("c005",))["overall"] == "undetermined"  # c100 missing
    rows = grid_rows("B", "old", "c005", ALL3) + grid_rows(
        "B", "flank", "c005", ALL3
    )  # no null replicates
    v = agg.verdict(agg.all_lods(agg.recovery_rows(rows)), rows, ["B"], ["c005"])
    assert v["overall"] == "undetermined"
    assert (
        _scene(ALL3, {100: 3}, layouts=("c005",))["overall"] == "retire"
    )  # c100 missing, c005 already fails


def test_overall_across_backgrounds():
    bad = _scene(ALL3, {300: 3}, bg="X")["per_background"]["X"]
    ok = _scene(ALL3, ALL3, bg="Y")["per_background"]["Y"]
    assert bad["verdict"] == "retire" and ok["verdict"] == "keep"


# ------------------------------------------------------------------------------- completeness gate on disk
def _write_cond(root, c, drop=None):
    d = Path(root) / "cond" / c["cond_id"]
    (d / "bundle" / "star").mkdir(parents=True, exist_ok=True)
    lv, spiked = int(c["level"]), int(c["level"]) > 0
    rows = [
        srow(c["background"], a, lv, int(c["seed"]), c["layout"], True, rep=int(c["bg_rep"]))
        for a in agg.ARMS
    ]
    for r in rows:
        r["cond_id"] = c["cond_id"]
    cols = list(rows[0])
    (d / "bundle" / "scores.tsv").write_text(
        "\t".join(cols) + "\n" + "".join("\t".join(r[k] for k in cols) + "\n" for r in rows)
    )
    if spiked:
        (d / "bundle" / "truth.tsv").write_text(
            "name\tcb\tclass\tcontigs\nSPIKE:a:0\tZZZ\tsense_in_gene\tVIRUS|T.1\n"
        )
        (d / "bundle" / "star" / "star_spike_fates.tsv").write_text(
            "name\ton_truth\tfate\nSPIKE:a:0\t1\tsource\n"
        )
    for arm in (
        agg.KB_ARMS
    ):  # kb tables: spike barcode ZZZ carries 2 truth molecules in the spiked run, 1 in the null
        a = d / "bundle" / arm
        a.mkdir(parents=True, exist_ok=True)
        (a / "virus_identity.tsv").write_text("gene_id\tgenome_accession\tviral\nT_g\tT.1\ttrue\n")
        n = 2 if spiked else 1
        with gzip.open(a / "kb_molecules.tsv.gz", "wt") as fh:
            fh.write(f"barcode\tgene_id\tunique\tambiguous\nZZZ\tT_g\t{n}\t0\n")
    if drop != "DONE":
        (d / "DONE").write_text("{}")


CONDS = [dict(cond_id="B-L003-s1-c005", background="B", level="3", seed="1", layout="c005", cells="5", bg_rep="1"),
         dict(cond_id="B-NULL-r1", background="B", level="0", seed="0", layout="", cells="0", bg_rep="1")]  # fmt: skip


def test_check_complete_and_report_runs(tmp_path):
    for c in CONDS:
        _write_cond(tmp_path, c)
    assert agg.check(tmp_path, CONDS) == []
    assert agg.run_report(tmp_path, tmp_path / "out", CONDS) == 0
    assert (tmp_path / "out" / "report.md").exists() and "NOT measured" in (
        tmp_path / "out" / "report.md"
    ).read_text()
    assert "overall" in (tmp_path / "out" / "verdict.tsv").read_text()


def test_check_flags_missing_failed_and_blank(tmp_path):
    _write_cond(tmp_path, CONDS[0], drop="DONE")
    p = agg.check(tmp_path, CONDS)
    assert any("B-L003-s1-c005: no DONE" in x for x in p) and any(
        "B-NULL-r1: no DONE" in x for x in p
    )
    _write_cond(tmp_path, CONDS[1])
    sp = tmp_path / "cond" / "B-NULL-r1" / "bundle" / "scores.tsv"
    lines = sp.read_text().splitlines()
    sp.write_text("\n".join(lines[:2] + lines[3:]) + "\n")  # drop the flank row
    assert any("arms" in x for x in agg.check(tmp_path, CONDS))
    _write_cond(tmp_path, CONDS[1])
    cols = lines[0].split("\t")
    row = lines[1].split("\t")
    row[cols.index("kb_unsupported_molecules")] = ""  # failed score: blank floor
    sp.write_text("\n".join([lines[0], "\t".join(row)] + lines[2:]) + "\n")
    assert any("blank" in x and "kb_unsupported_molecules" in x for x in agg.check(tmp_path, CONDS))


def test_report_refuses_without_writing_when_incomplete(tmp_path, capsys):
    _write_cond(tmp_path, CONDS[0])  # null missing
    assert agg.run_report(tmp_path, tmp_path / "out", CONDS) == 2
    assert not (tmp_path / "out").exists()
    assert "REFUSING" in capsys.readouterr().err
    assert agg.check(tmp_path) != []  # default = the real 66-condition grid: nothing there


def test_baseline_is_subtracted_and_flagged(tmp_path):
    for c in CONDS:
        _write_cond(tmp_path, c)
    rows = agg.apply_baseline(tmp_path, CONDS, agg.load_scores(tmp_path, CONDS))
    old = next(r for r in rows if r["arm"] == "old" and r["level"] == "3")
    # synthetic scores row says raw kb_total 10; the null replicate holds 1 in the spike barcode: net 9
    assert (old["kb_total_raw"], old["kb_total_in_baseline"], old["kb_total"]) == (10.0, 1.0, 9.0)
    assert old["call_kb"] == "True" and old["call_combined"] == "True"
    assert agg.baseline_rows(rows)[0]["max"] == 1.0


def test_background_alone_cannot_pass_the_kb_call(tmp_path):
    for c in CONDS:
        _write_cond(tmp_path, c)
    # null AND spiked both hold 5 molecules in the spike barcode: raw call would be True, net is 0
    for cid in ("B-L003-s1-c005", "B-NULL-r1"):
        with gzip.open(
            tmp_path / "cond" / cid / "bundle" / "old" / "kb_molecules.tsv.gz", "wt"
        ) as fh:
            fh.write("barcode\tgene_id\tunique\tambiguous\nZZZ\tT_g\t5\t0\n")
    sp = tmp_path / "cond" / "B-L003-s1-c005" / "bundle" / "scores.tsv"
    txt = sp.read_text().splitlines()
    cols = txt[0].split("\t")
    out = [txt[0]]
    for ln in txt[1:]:
        f = ln.split("\t")
        if f[cols.index("arm")] == "old":
            f[cols.index("kb_total")] = "5"
        out.append("\t".join(f))
    sp.write_text("\n".join(out) + "\n")
    rows = agg.apply_baseline(tmp_path, CONDS, agg.load_scores(tmp_path, CONDS))
    old = next(r for r in rows if r["arm"] == "old" and r["level"] == "3")
    assert old["kb_total"] == 0.0 and old["call_combined"] == "False"
    flag = next(b for b in agg.baseline_rows(rows) if b["arm"] == "old")
    assert flag["flagged"] is True


def test_check_requires_the_baseline_bundle(tmp_path):
    for c in CONDS:
        _write_cond(tmp_path, c)
    (tmp_path / "cond" / "B-NULL-r1" / "bundle" / "old" / "kb_molecules.tsv.gz").unlink()
    assert any("baseline" in x for x in agg.check(tmp_path, CONDS))


def test_wilson_and_class_table_marks_small_n(tmp_path):
    lo, hi = agg.wilson(3, 3)
    assert 0.4 < lo < 0.5 and hi == 1.0
    assert agg.wilson(0, 0) == (0.0, 1.0)
    for c in CONDS:
        _write_cond(tmp_path, c)
    assert agg.run_report(tmp_path, tmp_path / "out", CONDS) == 0
    md = (tmp_path / "out" / "report.md").read_text()
    assert "insufficient" in md and "Wilson" in md and "pseudo-replicates" in md
    assert "Background baseline" in md


# ----------------------------------------------------------------------------- PLAN SPIKE-03 hardening
def test_threshold_labels_censored_not_reached_undetermined():
    assert lod_of(grid_rows("B", "old", "c005", ALL3))["label"] == "<=3"
    assert lod_of(grid_rows("B", "old", "c005", ALL3))["status"] == "left_censored"
    assert lod_of(grid_rows("B", "old", "c005", {30: 3, 100: 3, 300: 3}))["label"] == "30"
    nr = lod_of(grid_rows("B", "old", "c005", {}))
    assert (nr["label"], nr["status"]) == (">300 (not established)", "not_reached")
    # top level: one seed false, one true, one missing -> could still be met: undetermined, never not_reached
    rows = [
        r
        for r in grid_rows("B", "old", "c005", {300: 1})
        if not (r["level"] == "300" and r["seed"] == "3")
    ]
    assert lod_of(rows)["status"] == "undetermined"
    # two seeds false, one missing: cannot be met even so -> genuinely not met
    rows = [
        r
        for r in grid_rows("B", "old", "c005", {})
        if not (r["level"] == "300" and r["seed"] == "3")
    ]
    assert lod_of(rows)["status"] == "not_reached"
    assert lod_of([])["status"] == "undetermined"  # no rows at all: nothing is known


def test_blank_call_is_missing_not_false():
    rows = grid_rows("B", "old", "c005", ALL3)
    for r in rows:
        if r["level"] == "300":
            r["call_combined"] = ""
    assert lod_of(rows)["status"] == "undetermined"


def test_non_monotone_flag_and_report_warning():
    m = {3: 0, 10: 3, 30: 1, 100: 2, 300: 3}
    t = lod_of(grid_rows("B", "old", "c005", m))
    assert t["non_monotone"] is True and t["label"] == "100"
    assert lod_of(grid_rows("B", "old", "c005", ALL3))["non_monotone"] is False
    rows = grid_rows("B", "old", "c005", m)
    rr = agg.recovery_rows(rows)
    md = agg.report_md(rr, agg.all_lods(rr), [], {"overall": "keep", "per_background": {}}, [], 1)
    assert "non-monotone" in md.lower() and "not a negative result" in md


def test_verdict_missing_condition_is_undetermined_and_backgrounds_from_grid():
    # one background (B) absent from the data entirely: grid.BACKGROUNDS decides what is expected
    v = agg.verdict([], [])
    assert set(v["per_background"]) == set(grid.BACKGROUNDS) and v["overall"] == "undetermined"
    # a missing seed on a decisive level cannot turn into a retire
    rows = grid_rows("B", "old", "c005", ALL3) + [
        r
        for r in grid_rows("B", "flank", "c005", ALL3)
        if not (r["level"] == "300" and r["seed"] != "1")
    ]
    for i in (1, 2, 3):
        rows += [srow("B", a, 0, 0, "", False, 0, rep=i) for a in agg.KB_ARMS]
    lods = agg.all_lods(agg.recovery_rows(rows))
    assert agg.verdict(lods, rows, ["B"], ["c005"])["overall"] == "undetermined"


def test_blank_unsupported_is_not_zero_in_verdict():
    rows = grid_rows("B", "old", "c005", ALL3) + grid_rows("B", "flank", "c005", ALL3)
    for i in (1, 2, 3):
        rows += [srow("B", a, 0, 0, "", False, 0, rep=i) for a in agg.KB_ARMS]
    for r in rows:
        if r["level"] == "0" and r["arm"] == "flank" and r["bg_rep"] == "3":
            r["kb_unsupported_molecules"] = ""
    v = agg.verdict(agg.all_lods(agg.recovery_rows(rows)), rows, ["B"], ["c005"])
    assert v["overall"] == "undetermined" and "paired null" in v["per_background"]["B"]["notes"]


def test_recovery_blank_stays_none_not_zero():
    r = srow("B", "old", 10, 1, "c005", True, n=8)
    r["kb_total"] = ""
    out = agg.recovery_rows([r])[0]
    assert out["n_recovered"] is None and out["recovery"] is None
    r["call_combined"] = ""
    assert agg.recovery_rows([r])[0]["call"] is None


def test_scores_name_flag_reads_scores_v2(tmp_path):
    for c in CONDS:
        _write_cond(tmp_path, c)
    assert agg.check(tmp_path, CONDS, "scores_v2.tsv") != []  # no v2 yet
    for c in CONDS:
        b = tmp_path / "cond" / c["cond_id"] / "bundle"
        (b / "scores_v2.tsv").write_text((b / "scores.tsv").read_text())
    assert agg.check(tmp_path, CONDS, "scores_v2.tsv") == []
    assert agg.run_report(tmp_path, tmp_path / "o", CONDS, "scores_v2.tsv") == 0
    assert "scores_v2.tsv" in (tmp_path / "o" / "report.md").read_text()


def test_q16_wording_in_report_and_tsv(tmp_path):
    for c in CONDS:
        _write_cond(tmp_path, c)
    assert agg.run_report(tmp_path, tmp_path / "o", CONDS) == 0
    md = (tmp_path / "o" / "report.md").read_text().lower()
    assert "computational detection threshold" in md and "limit of detection" not in md
    assert "threshold" in (tmp_path / "o" / "lod.tsv").read_text().splitlines()[0]

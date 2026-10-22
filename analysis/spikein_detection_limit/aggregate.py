"""Aggregate the spike-in grid: limit of detection, null floor, flank keep/retire verdict, STAR gate (PLAN SPIKE-02).

``python aggregate.py check  [--root R]``  completeness gate: every grid condition DONE with a valid scores.tsv
``python aggregate.py report [--root R] [--out DIR]``  runs the gate first and REFUSES (exit 2, nothing written) if
anything is missing; otherwise writes recovery_by_condition.tsv, lod.tsv, fp_floor.tsv, verdict.tsv, star_gate.tsv, report.md.

Rules, all pre-set by the user (2026-10-21/22) except where marked (decision):
 * call: old/flank = call_combined (>= 3 kb molecules AND >= 1 STAR read on the truth contig set; the STAR half comes from
   the old-index run for both kb arms); star = call_star. The kb-only call is reported as a secondary limit.
 * seed rule: a level is MET when the call is true in >= MIN_SEEDS (2) of the 3 seeds, per background x arm x layout.
 * monotone rule: the limit is the lowest level that is met AND every higher level is met; none = not reached.
 * BASELINE (review fix): the kb arm count is NET of the unspiked replicate. Spiked conditions use bg_rep = seed, so
   null replicate r<seed> is the same background subsample without spikes; ``apply_baseline`` subtracts its molecules
   (same spike barcodes, same truth accessions) from kb_total (floored at 0) and recomputes call_kb / call_combined.
   Raw total and baseline stay in ``kb_total_raw`` / ``kb_total_in_baseline``; a background whose baseline reaches
   KB_MIN is flagged in the report (its raw count alone would pass the call).
 * recovery = recovered / spiked: kb = net kb_total (unique + allocated ambiguous, spike barcodes, truth accessions,
   baseline subtracted) / n_spiked; star = reads with any alignment on a truth contig / n_spiked.
 * null floor: unspiked replicates, per arm: unsupported kb molecules per library (mean, max), false calls.
 * an old index that never reaches any level in a background x layout makes the comparison vacuous: undetermined, never keep.
 * keep flank iff, in EVERY background x layout, flank limit <= old limit (a tie passes; both not reached is a tie,
   flagged) AND flank adds <= ADD_MAX unsupported UMIs per null library. (decision) "adds" = per null replicate,
   flank minus old (paired by bg_rep), maximum over replicates; the absolute flank maximum is reported next to it.
   (decision) overall = retire if any background retires, else undetermined if any is, else keep.
"""

import argparse
import collections
import csv
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grid  # noqa: E402
import score_condition  # noqa: E402

MIN_SEEDS = 2
ADD_MAX = 5
MIN_NULL = 3
KB_MIN = score_condition.KB_MIN
MIN_CLASS_N = 30  # per-class STAR recovery below this many spikes is marked insufficient
INF = float("inf")
KB_ARMS = ("old", "flank")
ARMS = ("old", "flank", "star")
CALL_COL = {"old": "call_combined", "flank": "call_combined", "star": "call_star"}
FATES = ("source", "other_virus", "host", "other", "absent")


def truthy(v) -> bool:
    return str(v) == "True"


def read_tsv(path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow([f"{r[h]:.6g}" if isinstance(r[h], float) else r[h] for h in header])


def cond_dir(root, cid) -> Path:
    return Path(root) / "cond" / cid


# ------------------------------------------------------------------------------------------- completeness
def check(root, conds=None) -> list[str]:
    """Problems found in a finished grid (empty list = complete). Never raises on a missing/odd file."""
    conds = conds if conds is not None else grid.build_conditions()
    problems = []
    for c in conds:
        cid, d = c["cond_id"], cond_dir(root, c["cond_id"])
        spiked = int(c["level"]) > 0
        if not (d / "DONE").exists():
            problems.append(f"{cid}: no DONE marker (missing or failed)")
            continue
        sp = d / "bundle" / "scores.tsv"
        if not sp.exists():
            problems.append(f"{cid}: DONE but bundle/scores.tsv missing")
            continue
        try:
            rows = read_tsv(sp)
        except (OSError, csv.Error, UnicodeDecodeError) as e:
            problems.append(f"{cid}: scores.tsv unreadable ({e})")
            continue
        if sorted(r.get("arm", "") for r in rows) != sorted(ARMS):
            problems.append(
                f"{cid}: scores.tsv arms {[r.get('arm') for r in rows]}, need {list(ARMS)}"
            )
            continue
        for r in rows:
            need = ["call_star", "star_n_source_contig_any"] if r["arm"] == "star" else (
                ["call_kb", "call_combined", "kb_total", "kb_unsupported_molecules"])  # fmt: skip
            if not spiked:
                need = ["star_nonspike_reads_total"] if r["arm"] == "star" else (
                    ["kb_unsupported_molecules"])  # fmt: skip
            blank = [k for k in need if r.get(k, "") == ""]
            if blank:
                problems.append(f"{cid}/{r['arm']}: blank {blank}")
            if r.get("cond_id") != cid or r.get("level") != str(c["level"]):
                problems.append(
                    f"{cid}/{r['arm']}: row identity {r.get('cond_id')}/{r.get('level')} mismatch"
                )
            if spiked and int(r.get("n_spiked") or 0) <= 0:
                problems.append(f"{cid}/{r['arm']}: n_spiked is 0 for a spiked condition")
        if spiked:
            for rel in ("truth.tsv", "star/star_spike_fates.tsv"):
                if not (d / "bundle" / rel).exists():
                    problems.append(f"{cid}: bundle/{rel} missing (needed for per-class recovery)")
            nid = null_id(conds, c)
            if nid is None:
                problems.append(
                    f"{cid}: no unspiked replicate with bg_rep {c['bg_rep']} (baseline)"
                )
            else:
                for arm in KB_ARMS:
                    for f in ("kb_molecules.tsv.gz", "virus_identity.tsv"):
                        if not (cond_dir(root, nid) / "bundle" / arm / f).exists():
                            problems.append(f"{cid}: baseline {nid}/bundle/{arm}/{f} missing")
    return problems


def null_id(conds, c) -> str | None:
    """cond_id of the unspiked replicate on the same background subsample (bg_rep) as spiked condition c."""
    return next(
        (
            n["cond_id"]
            for n in conds
            if int(n["level"]) == 0
            and n["background"] == c["background"]
            and str(n["bg_rep"]) == str(c["bg_rep"])
        ),
        None,
    )


def apply_baseline(root, conds, rows) -> list[dict]:
    """Make the kb arms NET of the unspiked replicate (mutates and returns rows; see module doc)."""
    by_id = {c["cond_id"]: c for c in conds}
    for r in rows:
        if int(r["level"]) == 0 or r["arm"] not in KB_ARMS:
            continue
        c = by_id[r["cond_id"]]
        b = cond_dir(root, r["cond_id"]) / "bundle"
        truth = read_tsv(b / "truth.tsv")
        acc = {x.removeprefix("VIRUS|") for x in truth[0]["contigs"].split(",") if x}
        nb = cond_dir(root, null_id(conds, c)) / "bundle" / r["arm"]
        s = score_condition.kb_stats(nb, {t["cb"] for t in truth}, acc, {})
        base, raw = s["uniq"] + s["amb"], float(r["kb_total"])
        net = max(0.0, raw - base)
        r["kb_total_raw"], r["kb_total_in_baseline"], r["kb_total"] = raw, base, net
        r["call_kb"] = str(net >= KB_MIN)
        r["call_combined"] = str(net >= KB_MIN and truthy(r["call_star"]))
    return rows


def baseline_rows(rows) -> list[dict]:
    """Per background x arm: baseline molecules in the spike barcodes, flagged when >= KB_MIN."""
    out = []
    for bg, arm in sorted(
        {(r["background"], r["arm"]) for r in _spiked(rows) if r["arm"] in KB_ARMS}
    ):
        v = [
            float(r["kb_total_in_baseline"])
            for r in _spiked(rows)
            if (r["background"], r["arm"]) == (bg, arm)
            and r.get("kb_total_in_baseline") not in (None, "")
        ]
        if v:
            out.append(dict(background=bg, arm=arm, n=len(v), mean=statistics.fmean(v), max=max(v),
                            flagged=max(v) >= KB_MIN))  # fmt: skip
    return out


def wilson(k, n, z=1.96) -> tuple[float, float]:
    """Wilson score interval for k of n."""
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * (p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5 / d
    return max(0.0, c - h), min(1.0, c + h)


def load_scores(root, conds=None) -> list[dict]:
    conds = conds if conds is not None else grid.build_conditions()
    return [
        r for c in conds for r in read_tsv(cond_dir(root, c["cond_id"]) / "bundle" / "scores.tsv")
    ]


# ---------------------------------------------------------------------------------------------- statistics
def _spiked(rows):
    return [r for r in rows if int(r["level"]) > 0]


def recovery_rows(rows) -> list[dict]:
    out = []
    for r in _spiked(rows):
        n = int(r["n_spiked"])
        got = (
            float(r["kb_total"] or 0)
            if r["arm"] in KB_ARMS
            else float(r["star_n_source_contig_any"] or 0)
        )
        out.append(dict(
            background=r["background"], arm=r["arm"], layout=r["layout"], level=int(r["level"]),
            seed=int(r["seed"]), **{"class": "all"}, n_spiked=n, n_recovered=got,
            recovery=got / n if n else 0.0, call=truthy(r[CALL_COL[r["arm"]]]),
            call_kb=truthy(r["call_kb"]) if r["arm"] in KB_ARMS else "",
        ))  # fmt: skip
    return out


def lod(rows, background, arm, layout, call="call", levels=grid.LEVELS):
    """Limit of detection for one background x arm x layout from recovery_rows(); None = not reached.

    A level is met when >= MIN_SEEDS seeds have the call true (missing seeds count as not met); the limit is the
    lowest level such that it and every higher level are met."""
    k = collections.Counter(
        r["level"] for r in rows
        if (r["background"], r["arm"], r["layout"]) == (background, arm, layout) and r[call] is True
    )  # fmt: skip
    limit = None
    for lv in sorted(levels, reverse=True):
        if k[lv] < MIN_SEEDS:
            break
        limit = lv
    return limit


def all_lods(rrows) -> list[dict]:
    keys = sorted({(r["background"], r["arm"], r["layout"]) for r in rrows})
    out = []
    for bg, arm, layout in keys:
        calls = ["call"] + (["call_kb"] if arm in KB_ARMS else [])
        for col in calls:
            lim = lod(rrows, bg, arm, layout, col)
            out.append(dict(
                background=bg, arm=arm, layout=layout,
                call="combined" if col == "call" and arm in KB_ARMS else "star" if arm == "star" else "kb_only",
                limit="not_reached" if lim is None else lim,
                seeds_called_by_level=";".join(
                    f"{lv}:{sum(1 for r in rrows if (r['background'], r['arm'], r['layout'], r['level']) == (bg, arm, layout, lv) and r[col] is True)}"
                    for lv in grid.LEVELS
                ),
            ))  # fmt: skip
    return out


def fp_floor(rows) -> list[dict]:
    out = []
    nulls = [r for r in rows if int(r["level"]) == 0]
    for bg in sorted({r["background"] for r in nulls}):
        for arm in ARMS:
            sel = [r for r in nulls if r["background"] == bg and r["arm"] == arm]
            if not sel:
                continue
            metrics = (
                ("kb_unsupported_molecules", "n_false_calls") if arm in KB_ARMS
                else ("star_nonspike_reads_total", "star_contigs_with_reads")
            )  # fmt: skip
            for m in metrics:
                v = [float(r[m]) for r in sel]
                out.append(dict(background=bg, arm=arm, metric=m, n_libraries=len(v),
                                mean=statistics.fmean(v), max=max(v)))  # fmt: skip
    return out


def verdict(lods, rows, backgrounds=None, layouts=None) -> dict:
    """Flank keep/retire by the pre-set rule. lods = all_lods() rows (call 'combined'); rows = score rows."""
    backgrounds = backgrounds or sorted({r["background"] for r in lods})
    layouts = layouts or sorted(grid.LAYOUTS)
    lim = {
        (r["background"], r["arm"], r["layout"]): r["limit"]
        for r in lods
        if r["call"] == "combined"
    }
    per = {}
    for bg in backgrounds:
        notes, retire, undet = [], False, False
        for layout in layouts:
            o, f = lim.get((bg, "old", layout)), lim.get((bg, "flank", layout))
            if o is None or f is None:
                undet = True
                notes.append(f"{layout}: limit missing")
                continue
            o_, f_ = (INF if x == "not_reached" else x for x in (o, f))
            if f_ > o_:
                retire = True
                notes.append(f"{layout}: flank limit {f} > old limit {o}")
            elif o_ == INF:
                undet = True
                notes.append(
                    f"{layout}: old index never reached any level (flank {f}): comparison uninformative,"
                    " check the grid wiring"
                )
            else:
                notes.append(f"{layout}: flank {f} <= old {o}")
        by = {a: {int(r["bg_rep"]): float(r["kb_unsupported_molecules"] or 0) for r in rows
                  if r["background"] == bg and int(r["level"]) == 0 and r["arm"] == a} for a in KB_ARMS}  # fmt: skip
        reps = sorted(set(by["old"]) & set(by["flank"]))
        add_max = max((by["flank"][i] - by["old"][i] for i in reps), default=None)
        abs_max = max(by["flank"].values(), default=None)
        if len(reps) < MIN_NULL:
            undet = True
            notes.append(f"only {len(reps)} paired null replicates, need {MIN_NULL}")
        elif add_max > ADD_MAX:
            retire = True
            notes.append(
                f"flank adds up to {add_max:g} unsupported UMIs per null library > {ADD_MAX}"
            )
        else:
            notes.append(
                f"flank adds at most {add_max:g} unsupported UMIs per null library (<= {ADD_MAX})"
            )
        state = "retire" if retire else "undetermined" if undet else "keep"
        per[bg] = dict(
            verdict=state, add_max=add_max, flank_abs_max=abs_max, notes="; ".join(notes)
        )
    vs = [p["verdict"] for p in per.values()]
    overall = (
        "retire" if "retire" in vs else "undetermined" if "undetermined" in vs or not vs else "keep"
    )
    return dict(per_background=per, overall=overall)


def star_by_class(root, conds) -> list[dict]:
    """STAR fate per spike class, summed over seeds: (background, layout, level, class) -> counts (reads from bundles)."""
    acc = collections.defaultdict(collections.Counter)
    for c in conds:
        if int(c["level"]) == 0:
            continue
        b = cond_dir(root, c["cond_id"]) / "bundle"
        fates = {r["name"]: r for r in read_tsv(b / "star" / "star_spike_fates.tsv")}
        for t in read_tsv(b / "truth.tsv"):
            f = fates.get(t["name"])
            for cl in (t["class"], "all"):
                k = acc[(c["background"], c["layout"], int(c["level"]), cl)]
                k["n_spiked"] += 1
                k["on_truth_any"] += bool(f and f["on_truth"] == "1")
                k["fate_" + (f["fate"] if f else "absent")] += 1
    return [dict(background=bg, layout=ly, level=lv, **{"class": cl}, n_spiked=k["n_spiked"],
                 on_truth_any=k["on_truth_any"], **{"fate_" + x: k["fate_" + x] for x in FATES},
                 recovery=k["on_truth_any"] / k["n_spiked"])  # fmt: skip
            for (bg, ly, lv, cl), k in sorted(acc.items())]  # fmt: skip


# ------------------------------------------------------------------------------------------------ report
NOT_MEASURED = """\
- One source virus per background: HHV-6B reads into D03, EBV (type 1 + type 2) reads into R01_SOC4. No other virus, no
  other source library, and no source-background chemistry crossing was tested, so limits do not transfer to other viruses.
- NH policy: the EBV source keeps reads with NH <= 2 whose records all lie on the two EBV genomes (`--max-nh 2`, an
  unvalidated choice); HHV-6B keeps NH = 1 only. The effect of that choice on the limit was not measured.
- No error-free simulation arm in this grid (the `synth` probes were not run here): recovery losses are not separated
  into read-sequence effects versus index effects.
- "First/second source" caveat: spikes come from ONE real library per virus; seeds resample reads and barcodes from that
  library, not new libraries, so seed agreement reflects sampling within one source only.
- Per-class recovery is available for the STAR arm only: kb molecules are summed per barcode and gene, not per spike, so
  the kb arms have class "all" only. The kb arms get their STAR half of the combined call from the old-index run.
- STAR "absent" lumps reads unmapped for any reason (including too short after the 90 nt trim) and reads lost before the
  BAM (host-filter loss); the two are not separated here. The product host-filter loss column was not extracted.
- Limits are on a 5-level grid (3, 10, 30, 100, 300): a limit is a grid level, not an interpolated threshold. The null floor
  is 3 libraries per background; maxima of 3 values bound nothing beyond those libraries.
- Host and virus used the pinned indexes of this repository; depth-matching to real libraries (~45 M pairs) is by subsample.
"""


def _md(header, rows):
    return "\n".join(["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
                     + ["| " + " | ".join(str(x) for x in r) + " |" for r in rows])  # fmt: skip


def report_md(rrows, lods, fp, vd, sc, n_conds, base=()) -> str:
    out = ["# Spike-in detection limit: results", "",
           f"Grid complete: {n_conds} conditions scored. Rules and definitions: see `aggregate.py` docstring and README.", "",
           f"## Flank-index verdict: **{vd['overall'].upper()}**", "",
           *(["**UNDETERMINED: the comparison is uninformative or incomplete (old index never reached a level, "
              "or data missing); see the notes below. This is not a keep.**", ""] if vd["overall"] == "undetermined" else []),
           "Rule: keep iff flank limit <= old limit in every background x layout AND flank adds <= "
           f"{ADD_MAX} unsupported UMIs per null library (paired by replicate, max).", "",
           _md(["background", "verdict", "max added unsupported UMIs (flank - old)", "flank absolute max", "notes"],
               [[b, p["verdict"], p["add_max"], p["flank_abs_max"], p["notes"]]
                for b, p in vd["per_background"].items()]), "",
           "## Limit of detection (lowest level met in >= 2 of 3 seeds, and every higher level met)", "",
           "Cell: limit; in brackets seeds with the call by level 3/10/30/100/300.", ""]  # fmt: skip
    for bg in sorted({r["background"] for r in lods}):
        out += [f"### {bg}", ""]
        for call in ("combined", "star", "kb_only"):
            sel = [r for r in lods if r["background"] == bg and r["call"] == call]
            if sel:
                out += [f"Call: {call}", "", _md(["arm", "layout", "limit", "seeds called by level"],
                        [[r["arm"], r["layout"], r["limit"], r["seeds_called_by_level"]] for r in sel]), ""]  # fmt: skip
    out += ["## Mean recovery (recovered / spiked) by level", ""]
    grp = collections.defaultdict(list)
    for r in rrows:
        grp[(r["background"], r["arm"], r["layout"], r["level"])].append(r["recovery"])
    keys = sorted({k[:3] for k in grp})
    out += [_md(["background", "arm", "layout"] + [str(lv) for lv in grid.LEVELS],
                [list(k) + [f"{statistics.fmean(grp[(*k, lv)]):.2f}" if grp.get((*k, lv)) else "NA" for lv in grid.LEVELS]
                 for k in keys]), "",
            "kb arms: NET kb molecules (unique + allocated ambiguous, spike barcodes, truth accessions, minus the unspiked "
            "replicate of the same background subsample) per spiked read, not capped. star: any alignment on a truth contig.", "",
            "## Background baseline of the truth virus in the spike barcodes (unspiked replicate)", "",
            _md(["background", "arm", "spiked conditions", "mean", "max", f"flag (max >= {KB_MIN})"],
                [[b["background"], b["arm"], b["n"], f"{b['mean']:.3g}", f"{b['max']:g}",
                  "BASELINE >= KB_MIN: raw count alone would pass the call" if b["flagged"] else "ok"] for b in base]), "",
            "## Null floor (unspiked replicates)", "",
            _md(["background", "arm", "metric", "libraries", "mean", "max"],
                [[r["background"], r["arm"], r["metric"], r["n_libraries"], f"{r['mean']:.3g}", f"{r['max']:g}"] for r in fp]), "",
            "## STAR gate evidence: fate of spikes by class (all levels, seeds and layouts pooled per background)", "",
            "Pooled over levels, seeds and layouts: seeds are resamples of ONE source library per virus (pseudo-replicates, "
            "not independent libraries), the high levels dominate the pool, and class composition differs between "
            f"backgrounds, so pooled recovery is not comparable across backgrounds. n is shown beside recovery with a 95% "
            f"Wilson interval; classes with n < {MIN_CLASS_N} are marked insufficient.", ""]  # fmt: skip
    pool = collections.defaultdict(collections.Counter)
    for r in sc:
        for k in ("n_spiked", "on_truth_any", *["fate_" + x for x in FATES]):
            pool[(r["background"], r["class"])][k] += r[k]

    def _ci(k):
        lo, hi = wilson(k["on_truth_any"], k["n_spiked"])
        return f"[{lo:.2f}, {hi:.2f}]"

    out += [_md(["background", "class", "spiked (n)", "on truth contig", "recovery", "95% CI", "n check", "host", "other virus", "other", "absent"],
                [[b, cl, k["n_spiked"], k["on_truth_any"], f"{k['on_truth_any'] / k['n_spiked']:.3f}", _ci(k),
                  "insufficient" if k["n_spiked"] < MIN_CLASS_N else "ok", k["fate_host"],
                  k["fate_other_virus"], k["fate_other"], k["fate_absent"]] for (b, cl), k in sorted(pool.items())]), "",
            "Fates are the PRIMARY alignment (host = host competition; other virus = primary on another viral contig); "
            "on-truth counts any alignment on the truth set.", "",
            "## NOT measured", "", NOT_MEASURED]  # fmt: skip
    return "\n".join(out)


def run_report(root, out_dir, conds=None) -> int:
    conds = conds if conds is not None else grid.build_conditions()
    problems = check(root, conds)
    if problems:
        print(f"aggregate: REFUSING to report, {len(problems)} problem(s):", file=sys.stderr)
        print("\n".join(problems), file=sys.stderr)
        return 2
    rows = apply_baseline(root, conds, load_scores(root, conds))
    rrows, lods, fp = recovery_rows(rows), None, fp_floor(rows)
    lods = all_lods(rrows)
    vd, sc = verdict(lods, rows), star_by_class(root, conds)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "recovery_by_condition.tsv", list(rrows[0]), rrows)
    write_tsv(out / "lod.tsv", list(lods[0]), lods)
    write_tsv(out / "fp_floor.tsv", list(fp[0]), fp)
    write_tsv(out / "verdict.tsv", ["background", "verdict", "add_max", "flank_abs_max", "notes"],
              [dict(background=b, **p) for b, p in vd["per_background"].items()]
              + [dict(background="overall", verdict=vd["overall"], add_max="", flank_abs_max="", notes="")])  # fmt: skip
    write_tsv(out / "star_gate.tsv", list(sc[0]), sc)
    (out / "report.md").write_text(
        report_md(rrows, lods, fp, vd, sc, len(conds), baseline_rows(rows))
    )
    print(f"verdict: {vd['overall']}; wrote {out}", file=sys.stderr)
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("cmd", choices=["check", "report"])
    p.add_argument("--root", default=str(grid.ROOT))
    p.add_argument("--out", help="report directory (default ROOT/aggregate)")
    a = p.parse_args(argv)
    if a.cmd == "check":
        problems = check(a.root)
        print(
            "\n".join(problems)
            if problems
            else f"complete: {len(grid.build_conditions())} conditions"
        )
        sys.exit(1 if problems else 0)
    sys.exit(run_report(a.root, a.out or Path(a.root) / "aggregate"))


if __name__ == "__main__":
    main()

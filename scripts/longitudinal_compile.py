#!/usr/bin/env python3
"""LONG-01 / LONG-02: within-subject viral trajectories and the cross-library contamination demotion.

    python scripts/longitudinal_compile.py --candidates adjudicated.tsv --libraries libraries.tsv \
        --config long.yaml --out-dir out/ [--loci loci.tsv]

Inputs (TSV; extra columns are ignored):

  --candidates  one row per (library, virus) from the adjudicator (ADJ-01): library_id, virus_key, level (L0..L5),
                level_reason, flags (';' list), unique_molecules, n_loci, read_virus_best_fraction,
                thresholds_sha256; optional n_positive_cells; pass-through subject_id, timepoint, run_id, lane.
  --libraries   one row per library, the universe of libraries (a library with no candidate row is "virus not
                detected"): library_id, host_umis or host_reads, n_cells; subject_id, timepoint, run_id, lane
                (filled from the candidates when blank). Visit order comes from the first of `timepoint_order`,
                `biopsy_order`, `time_from_transplant_days` that every library of the subject has, else from a numeric
                `timepoint`. A subject with several visits and no order is skipped, never guessed.
  --loci        optional, one row per (library, virus) evidence key: library_id, virus_key and `locus_id` and/or
                `molecule_key`. Absence of a key means no evidence, so the table must cover every library of a run.

Outputs in --out-dir: `longitudinal_trajectories.tsv` (one row per subject, virus, visit), `longitudinal_events.tsv`
(one row per consecutive visit pair, per virus that is detected at either visit) and `longitudinal_manifest.json`.

Event labels: new_detection, burden_increase, burden_decrease, persistent, cleared, not_testable,
demoted_contamination. The primary metric is depth-normalised burden (viral unique molecules / host UMIs, else
host reads); raw counts, positive-cell fraction and locus breadth are reported beside it. Libraries of one
(subject, visit) are pooled. An event never makes an L4 call: `l4_eligible` only says that nothing blocks one, and
is true only for a new_detection or burden_increase, at level >= L2, with no blocking flag, that passed the
contamination check. Every numeric cutoff comes from --config (LOCK-01 calibrates them); with no
`burden_fold_change` a pair detected at both visits is `not_testable`.

LONG-02: a new_detection or burden_increase becomes `demoted_contamination` when a library of another subject in the
same run (or lane) holds the same locus_id or molecule_key of that virus as a library of the later visit. A lane
matches only when the two run ids are not both known and different, so give the lane flowcell-qualified.
"""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import yaml

LEVELS = {f"L{i}": i for i in range(6)}
UNKNOWN = {"", "UNKNOWN", "NA", "NAN", "NONE"}
ORDER_COLS = ("timepoint_order", "biopsy_order", "time_from_transplant_days")
DEFAULTS = {
    "burden_fold_change": None,  # > 1; no default, calibrated at LOCK-01
    "min_host_umis": None,  # minimum host depth for a library to be testable; calibrated at LOCK-01
    "min_detected_level": "L1",  # a visit is "detected" at or above this level
    "not_testable_flags": ["not_testable"],
    "l4_blocking_flags": ["germline_integration_compatible", "screening_only"],
    "contamination_match": ["molecule_key", "locus_id"],
    "contamination_min_shared": 1,  # shared keys of one type needed to demote
}
TRAJ_COLS = [
    "subject_id", "virus_key", "visit_index", "timepoint", "order", "library_ids", "run_ids", "lanes", "level",
    "detected", "testable", "not_testable_reason", "flags", "unique_molecules", "host_depth", "depth_basis",
    "burden", "n_cells", "n_positive_cells", "positive_cell_fraction", "n_loci", "thresholds_sha256",
]
EVENT_COLS = [
    "subject_id", "virus_key", "from_timepoint", "to_timepoint", "from_libraries", "to_libraries", "event",
    "event_reason", "original_event", "from_burden", "to_burden", "fold_change", "depth_basis",
    "from_unique_molecules", "to_unique_molecules", "from_positive_cell_fraction", "to_positive_cell_fraction",
    "from_n_loci", "to_n_loci", "from_level", "to_level", "flags", "contamination_check",
    "contamination_evidence", "l4_eligible", "l4_block_reason",
]


def _die(msg):
    sys.exit(f"longitudinal_compile: {msg}")


def _num(x):
    try:
        v = float(str(x).strip())
    except (TypeError, ValueError):
        return None
    return v if v == v else None


def _known(x):
    return str(x or "").strip().upper() not in UNKNOWN


def _flags(x):
    return {f.strip() for f in str(x or "").split(";") if f.strip()}


def _level(x):
    s = str(x or "").strip().upper()
    if s not in LEVELS:
        _die(f"level {x!r} is not one of {sorted(LEVELS)}")
    return LEVELS[s]


def _fmt(v):
    return "" if v is None else f"{v:.6g}" if isinstance(v, float) else str(v)


def read_tsv(path):
    with Path(path).open(newline="") as fh:
        return [{k.strip(): (v or "").strip() for k, v in r.items() if k} for r in csv.DictReader(fh, delimiter="\t")]


def load_config(path):
    cfg = dict(DEFAULTS)
    user = (yaml.safe_load(Path(path).read_text()) or {}) if path else {}
    if unknown := set(user) - set(DEFAULTS):
        _die(f"unknown config keys {sorted(unknown)}; known: {sorted(DEFAULTS)}")
    cfg.update(user)
    fold = cfg["burden_fold_change"]
    if fold is not None and not (isinstance(fold, (int, float)) and fold > 1):
        _die("burden_fold_change must be a number > 1 (or unset)")
    if cfg["min_detected_level"] not in LEVELS:
        _die(f"min_detected_level must be one of {sorted(LEVELS)}")
    if not set(cfg["contamination_match"]) <= {"molecule_key", "locus_id"}:
        _die("contamination_match must be a subset of [molecule_key, locus_id]")
    return cfg


def build_libraries(lib_rows, cand_rows):
    """library_id -> metadata; the libraries table is authoritative, candidate pass-through fills blanks."""
    libs = {}
    for r in lib_rows:
        if r["library_id"] in libs:
            _die(f"duplicate library_id {r['library_id']!r} in the libraries table")
        libs[r["library_id"]] = dict(r)
    for r in cand_rows:
        lib = libs.setdefault(r["library_id"], {"library_id": r["library_id"]})
        for k in ("subject_id", "timepoint", "run_id", "lane"):
            if not _known(lib.get(k)) and _known(r.get(k)):
                lib[k] = r[k]
    for lib in libs.values():
        umis, reads = _num(lib.get("host_umis")), _num(lib.get("host_reads"))
        lib["_depth"], lib["_basis"] = (umis, "host_umis") if umis else (reads, "host_reads") if reads else (None, "")
    return libs


def subject_order(libs):
    """{library_id: order value} for one subject's libraries, or None when the visits cannot be ordered."""
    if len({lib.get("timepoint", "") for lib in libs}) == 1 and _known(libs[0].get("timepoint")):
        return {lib["library_id"]: 0.0 for lib in libs}
    if any(not _known(lib.get("timepoint")) for lib in libs):
        return None
    for col in (*ORDER_COLS, "timepoint"):
        vals = [_num(lib.get(col)) for lib in libs]
        if all(v is not None for v in vals):
            return {lib["library_id"]: v for lib, v in zip(libs, vals, strict=True)}
    return None


def _visit(lib_ids, libs, cands):
    """Pool the libraries of one (subject, visit) into the per-virus metrics."""
    ls = [libs[i] for i in lib_ids]
    depths = [lib["_depth"] for lib in ls]
    bases = {lib["_basis"] for lib in ls}
    ncells = [_num(lib.get("n_cells")) for lib in ls]
    v = {
        "libs": lib_ids,
        "timepoint": ";".join(sorted({lib.get("timepoint", "") for lib in ls})),
        "run_ids": ";".join(sorted({lib["run_id"] for lib in ls if _known(lib.get("run_id"))})),
        "lanes": ";".join(sorted({lib["lane"] for lib in ls if _known(lib.get("lane"))})),
        "depth": sum(depths) if all(depths) and len(bases) == 1 else None,
        "basis": bases.pop() if len(bases) == 1 else "",
        "n_cells": sum(ncells) if all(n is not None for n in ncells) else None,
        "thr": {r["thresholds_sha256"] for i in lib_ids for rs in cands.get(i, {}).values() for r in rs
                if r.get("thresholds_sha256")},
    }
    return v


def _metrics(visit, virus, cands, cfg):
    rows = [r for i in visit["libs"] for r in cands.get(i, {}).get(virus, [])]
    level = max((_level(r["level"]) for r in rows), default=None)
    flags = set().union(*(_flags(r.get("flags")) for r in rows)) if rows else set()
    mol = [_num(r.get("unique_molecules")) for r in rows]
    umol = sum(mol) if rows and all(m is not None for m in mol) else None
    npos = [_num(r.get("n_positive_cells")) for r in rows]
    npos = sum(npos) if rows and all(n is not None for n in npos) else None
    nloci = [_num(r.get("n_loci")) for r in rows]
    nloci = max(nloci) if rows and all(n is not None for n in nloci) else None  # ponytail: max, not the union
    depth, basis = visit["depth"], visit["basis"]
    why = ""
    if not depth:
        why = "no_host_depth"
    elif cfg["min_host_umis"] is not None and basis == "host_umis" and depth < cfg["min_host_umis"]:
        why = "below_min_host_umis"
    elif cfg["min_host_umis"] is not None and basis != "host_umis":
        why = "min_host_umis_needs_host_umis"
    elif bad := sorted(flags & set(cfg["not_testable_flags"])):
        why = "flag:" + ",".join(bad)
    return {
        "level": level,
        "detected": level is not None and level >= LEVELS[cfg["min_detected_level"]],
        "flags": flags,
        "umol": umol,
        "burden": umol / depth if umol is not None and depth else None,
        "npos": npos,
        "pos_frac": npos / visit["n_cells"] if npos is not None and visit["n_cells"] else None,
        "nloci": nloci,
        "testable": not why,
        "why": why,
    }


def classify(a, b, va, vb, cfg):
    """(event, reason) for a consecutive pair of one virus; None when it is absent at both visits."""
    if va["thr"] and vb["thr"] and va["thr"] != vb["thr"]:
        return "not_testable", "thresholds_sha256_differs_between_visits"
    if not a["detected"] and not b["detected"]:
        return None
    if not a["detected"]:
        if not a["testable"]:
            return "not_testable", f"earlier_visit_not_testable:{a['why']}"
        return "new_detection", "absent at a testable earlier visit"
    if not b["detected"]:
        if not b["testable"]:
            return "not_testable", f"later_visit_not_testable:{b['why']}"
        return "cleared", "not detected at a testable later visit"
    for tag, m in (("earlier", a), ("later", b)):
        if not m["testable"]:
            return "not_testable", f"{tag}_visit_not_testable:{m['why']}"
        if not m["burden"]:
            return "not_testable", f"{tag}_visit_has_no_burden"
    if va["basis"] != vb["basis"]:
        return "not_testable", "depth_basis_differs_between_visits"
    fold = cfg["burden_fold_change"]
    if fold is None:
        return "not_testable", "burden_fold_change_unset"
    ratio = b["burden"] / a["burden"]
    if ratio >= fold:
        return "burden_increase", f"burden fold change {ratio:.4g} >= {fold:g}"
    if ratio <= 1 / fold:
        return "burden_decrease", f"burden fold change {ratio:.4g} <= 1/{fold:g}"
    return "persistent", f"burden fold change {ratio:.4g} within 1/{fold:g}..{fold:g}"


def neighbours(lib, subject, libs):
    """Other-subject libraries in the same run, or the same lane unless the runs are known to differ."""
    out = []
    for other in sorted(libs):
        o = libs[other]
        if o.get("subject_id") == subject:
            continue
        run_known = _known(lib.get("run_id")) and _known(o.get("run_id"))
        same_run = run_known and lib["run_id"] == o["run_id"]
        same_lane = (_known(lib.get("lane")) and _known(o.get("lane")) and lib["lane"] == o["lane"]
                     and not (run_known and not same_run))
        if same_run or same_lane:
            out.append(other)
    return out


def contamination_check(virus, subject, lib_ids, libs, keys, cfg):
    """(status, evidence) for the libraries of the later visit. status: demoted | clear | not_testable:<why>."""
    if keys is None:
        return "not_testable:no_key_table", ""
    results = []
    for lid in lib_ids:
        own = keys.get((lid, virus))
        nb = neighbours(libs[lid], subject, libs)
        if not own:
            results.append(("not_testable:no_keys_for_library", ""))
        elif not (_known(libs[lid].get("run_id")) or _known(libs[lid].get("lane"))):
            results.append(("not_testable:no_run_or_lane", ""))
        else:
            hits = []
            for other in nb:
                for kind in cfg["contamination_match"]:
                    shared = own.get(kind, set()) & keys.get((other, virus), {}).get(kind, set())
                    if len(shared) >= cfg["contamination_min_shared"]:
                        hits.append(f"{lid}~{other}({libs[other].get('subject_id', '')}):{kind}x{len(shared)}")
            results.append(("demoted", ";".join(hits)) if hits else ("clear", ""))
    for want in ("demoted", "clear"):
        for status, ev in results:
            if status == want:
                return status, ev
    return results[0] if results else ("not_testable:no_key_table", "")


def l4_gate(event, to_level, flags, contam, cfg):
    reasons = []
    if event not in ("new_detection", "burden_increase"):
        reasons.append(f"event_is_{event}")
    if to_level is None or to_level < LEVELS["L2"]:
        reasons.append("level_below_L2")
    reasons += [f"flag:{f}" for f in sorted(flags & set(cfg["l4_blocking_flags"]))]
    if contam != "clear" and event in ("new_detection", "burden_increase"):
        reasons.append(f"contamination_check:{contam}")
    return not reasons, ";".join(reasons)


def compile_longitudinal(cand_rows, lib_rows, loci_rows, cfg):
    """Returns (trajectory rows, event rows, manifest dict). Deterministic."""
    libs = build_libraries(lib_rows, cand_rows)
    cands = {}
    for r in cand_rows:
        cands.setdefault(r["library_id"], {}).setdefault(r["virus_key"], []).append(r)
    keys = None
    if loci_rows is not None:
        keys = {}
        for r in loci_rows:
            for kind in ("molecule_key", "locus_id"):
                if r.get(kind):
                    keys.setdefault((r["library_id"], r["virus_key"]), {}).setdefault(kind, set()).add(r[kind])
    by_subject, excluded = {}, []
    for lid in sorted(libs):
        s = libs[lid].get("subject_id", "")
        if not _known(s) or ";" in s:
            excluded.append(lid)
        else:
            by_subject.setdefault(s, []).append(libs[lid])
    traj, events, skipped = [], [], []
    for subject in sorted(by_subject):
        order = subject_order(by_subject[subject])
        if order is None:
            skipped.append({"subject_id": subject, "reason": "timepoint_order_unknown"})
            continue
        groups = {}
        for lid, o in order.items():
            groups.setdefault(o, []).append(lid)
        visits = [(o, _visit(sorted(groups[o]), libs, cands)) for o in sorted(groups)]
        viruses = sorted({v for _, vis in visits for i in vis["libs"] for v in cands.get(i, {})})
        for virus in viruses:
            ms = [_metrics(vis, virus, cands, cfg) for _, vis in visits]
            if not any(m["detected"] for m in ms):
                continue
            for idx, ((o, vis), m) in enumerate(zip(visits, ms, strict=True)):
                traj.append({
                    "subject_id": subject, "virus_key": virus, "visit_index": idx, "timepoint": vis["timepoint"],
                    "order": o, "library_ids": ";".join(vis["libs"]), "run_ids": vis["run_ids"],
                    "lanes": vis["lanes"], "level": "" if m["level"] is None else f"L{m['level']}",
                    "detected": m["detected"], "testable": m["testable"], "not_testable_reason": m["why"],
                    "flags": ";".join(sorted(m["flags"])), "unique_molecules": m["umol"],
                    "host_depth": vis["depth"], "depth_basis": vis["basis"], "burden": m["burden"],
                    "n_cells": vis["n_cells"], "n_positive_cells": m["npos"], "positive_cell_fraction": m["pos_frac"],
                    "n_loci": m["nloci"], "thresholds_sha256": ";".join(sorted(vis["thr"])),
                })
            for k in range(1, len(visits)):
                (_, va), (_, vb), a, b = visits[k - 1], visits[k], ms[k - 1], ms[k]
                res = classify(a, b, va, vb, cfg)
                if res is None:
                    continue
                event, reason = res
                original, contam, evidence = "", "not_applicable", ""
                if event in ("new_detection", "burden_increase"):
                    contam, evidence = contamination_check(virus, subject, vb["libs"], libs, keys, cfg)
                    if contam == "demoted":
                        original, event = event, "demoted_contamination"
                        reason = f"{original}: same loci/molecules in other subjects' libraries ({evidence})"
                flags = a["flags"] | b["flags"]
                eligible, block = l4_gate(original or event, b["level"], flags, contam, cfg)
                events.append({
                    "subject_id": subject, "virus_key": virus, "from_timepoint": va["timepoint"],
                    "to_timepoint": vb["timepoint"], "from_libraries": ";".join(va["libs"]),
                    "to_libraries": ";".join(vb["libs"]), "event": event, "event_reason": reason,
                    "original_event": original, "from_burden": a["burden"], "to_burden": b["burden"],
                    "fold_change": b["burden"] / a["burden"] if a["burden"] and b["burden"] is not None else None,
                    "depth_basis": va["basis"] if va["basis"] == vb["basis"] else "",
                    "from_unique_molecules": a["umol"], "to_unique_molecules": b["umol"],
                    "from_positive_cell_fraction": a["pos_frac"], "to_positive_cell_fraction": b["pos_frac"],
                    "from_n_loci": a["nloci"], "to_n_loci": b["nloci"],
                    "from_level": "" if a["level"] is None else f"L{a['level']}",
                    "to_level": "" if b["level"] is None else f"L{b['level']}", "flags": ";".join(sorted(flags)),
                    "contamination_check": contam, "contamination_evidence": evidence,
                    "l4_eligible": eligible, "l4_block_reason": block,
                })
    manifest = {"n_libraries": len(libs), "excluded_pooled_or_unknown_subject": excluded, "skipped_subjects": skipped}
    return traj, events, manifest


def _write(path, cols, rows):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([_fmt(r[c]) for c in cols])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--libraries", required=True)
    ap.add_argument("--loci")
    ap.add_argument("--config")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    cands, libs = read_tsv(a.candidates), read_tsv(a.libraries)
    for name, rows, need in (("candidates", cands, ("library_id", "virus_key", "level")), ("libraries", libs, ("library_id",))):
        if rows and (miss := [c for c in need if c not in rows[0]]):
            _die(f"{name} table lacks column(s) {miss}")
    loci = read_tsv(a.loci) if a.loci else None
    if loci and not ({"library_id", "virus_key"} <= set(loci[0]) and {"locus_id", "molecule_key"} & set(loci[0])):
        _die("loci table needs library_id, virus_key and locus_id and/or molecule_key")
    traj, events, manifest = compile_longitudinal(cands, libs, loci, cfg)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    _write(out / "longitudinal_trajectories.tsv", TRAJ_COLS, traj)
    _write(out / "longitudinal_events.tsv", EVENT_COLS, events)
    manifest |= {
        "config": cfg,
        "inputs_sha256": {p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                          for p in (a.candidates, a.libraries, a.loci, a.config) if p},
        "n_events": len(events),
        "event_counts": {e: sum(r["event"] == e for r in events) for e in sorted({r["event"] for r in events})},
    }
    (out / "longitudinal_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    for s in manifest["skipped_subjects"]:
        print(f"longitudinal_compile: skipped {s['subject_id']}: {s['reason']}", file=sys.stderr)


if __name__ == "__main__":
    main()

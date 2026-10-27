#!/usr/bin/env python3
"""DISC-12 metrics and decision rules (protocol docs/plans/2026-10-22-discovery-benchmark.md §4, §5).

Committed before any result was read (blind). One row per arm of ``<run>/conditions.tsv`` from ``<run>/out/<id>/``
and ``<run>/receipts/<id>.time.err``; library roles from the input manifest (``positive_*`` / ``negative_*``).

    python scripts/disc12_metrics.py --run <run root> --manifest docs/plans/disc12/input_manifest.tsv \
        --panel-taxid SRR12682296=10376 --truth SRR12682296=<per-cell tsv> ... --out <dir>

Truth is not derivable from the run, so it is an explicit input for every positive library (missing = exit):
``--panel-taxid LIB=TAXID`` (the panel virus species) and ``--truth LIB=TSV`` with columns ``barcode`` and
``molecules`` (the primary run's per-cell UNIQUE molecules of that virus). A positive whose truth has no cell with
molecules > 0 is ``capture_only`` (protocol §2: covid tests capture, not recovery): it is left out of the recovery
criterion of rule 1 and kept in the capture criterion of rule 2.

Definitions fixed here (protocol ambiguities, see PLAN DISC-12-METRICS):
- panel lineage = the panel species, its descendants and its ancestors; panel groups = resolved groups in the
  species subtree.
- ``prioritized(k)`` (plan §9.2) = a ``candidates.tsv`` row with ``groups_called >= k``, ``host_lineage == 0`` and
  ``viral_only_db == 0``. Off-cell enrichment is reported, not gated; ``artefact_supported_only`` does not exist in
  ``candidates.tsv`` and is not applied; the negative-control comparison is ``not_assessed``.
- false-positive taxa exclude the host lineage in negatives too.
- G15 needs no extra code: viral-only rows fail the eligibility gate, so that DB gets rule 2 ``none`` and rule 1
  ``not_assessed``.
- rule 2 is evaluated per DB and configuration; rule 1 uses the rule-2 ``k`` of the shipped default (0.05 / 4, trim
  on) of the same DB; the ``--no-polya-trim`` arm is reported but is not a rule-1 candidate; "both negatives" is read
  as every negative; the final rule-1 tie-break after confidence is higher minimum hit groups.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

DEFAULT = (0.05, 4)
TIME_FIELDS = {
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": "wall_s",
    "Maximum resident set size (kbytes)": "max_rss_kb",
    "File system inputs": "fs_inputs_blocks",
}
METRIC_COLUMNS = [
    "cond_id", "library", "role", "db_name", "confidence", "min_hit_groups", "polya_trim", "status",
    "panel_taxid", "panel_reached", "n_truth_cells", "n_truth_called", "truth_cell_fraction",
    "median_cell_recovery", "spearman", "n_fp_observation", "max_eligible_groups_nontruth",
    "max_eligible_groups_panel", "rank_exact", "rank_ancestor", "rank_wrong_branch", "host_leakage",
    "panel_called_fraction", "wall_s", "max_rss_kb", "bytes_read",
]  # fmt: skip


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pairs(items: list[str], flag: str) -> dict[str, str]:
    out = {}
    for item in items:
        lib, sep, value = item.partition("=")
        if not sep or not lib or not value:
            sys.exit(f"{flag} expects LIB=VALUE, got {item!r}")
        out[lib] = value
    return out


def parse_time(path: Path) -> dict[str, float | str]:
    """GNU ``/usr/bin/time -v`` receipt: wall seconds, peak RSS (kB), bytes read (512-byte input blocks)."""
    if not path.is_file():
        return {"wall_s": "", "max_rss_kb": "", "bytes_read": ""}
    got: dict[str, float] = {}
    for line in path.read_text().splitlines():
        key, _, value = line.strip().rpartition(": ")
        if key in TIME_FIELDS:
            if TIME_FIELDS[key] == "wall_s":
                got["wall_s"] = sum(
                    float(p) * 60**i for i, p in enumerate(reversed(value.split(":")))
                )
            else:
                got[TIME_FIELDS[key]] = float(value)
    blocks = got.pop("fs_inputs_blocks", None)
    return {
        "wall_s": got.get("wall_s", ""),
        "max_rss_kb": got.get("max_rss_kb", ""),
        "bytes_read": "" if blocks is None else blocks * 512,
    }


class Tree:
    def __init__(self, rows: list[dict[str, str]]):
        self.parent = {int(r["taxid"]): int(r["parent_taxid"]) for r in rows}
        self.rank = {int(r["taxid"]): r["rank"] for r in rows}

    def path(self, taxid: int) -> list[int]:
        """``taxid`` up to the root, inclusive."""
        out = [taxid]
        while self.parent.get(out[-1], out[-1]) != out[-1] and len(out) < 200:
            out.append(self.parent[out[-1]])
        return out

    def genus(self, taxid: int) -> int | None:
        return next((t for t in self.path(taxid) if self.rank.get(t) == "G"), None)


def arm_metrics(cond: dict[str, str], run: Path, role: str, panel: int | None, truth) -> dict:
    out_dir = run / "out" / cond["cond_id"]
    row = {
        "cond_id": cond["cond_id"],
        "library": cond["library"],
        "role": role,
        "db_name": cond["db_name"],
        "confidence": float(cond["confidence"]),
        "min_hit_groups": int(cond["min_hit_groups"]),
        "polya_trim": int(cond["polya_trim"]),
        "panel_taxid": panel or "",
        **parse_time(run / "receipts" / f"{cond['cond_id']}.time.err"),
    }
    if not (out_dir / "discovery_complete.json").is_file():
        return {**row, "status": "failed"}
    row["status"] = "ok"
    tree = Tree(tsv(out_dir / "taxonomy.tsv"))
    cands = tsv(out_dir / "candidates.tsv")
    called = {int(c["taxid"]): int(c["groups_called"]) for c in cands}
    total = sum(int(c["groups_called"]) + int(c["groups_offcell"]) for c in cands)
    eligible = [
        c for c in cands if c["host_lineage"] == "0" and c["viral_only_db"] == "0"
    ]  # prioritized gates other than k
    if cands and cands[0]["db_has_9606"] == "1":
        host = sum(
            int(c["groups_called"]) + int(c["groups_offcell"])
            for c in cands
            if c["host_relation"] == "host"
        )
        row["host_leakage"] = host / total if total else ""
    in_panel = (lambda t: panel in tree.path(t)) if panel else (lambda t: False)
    ancestors = set(tree.path(panel)[1:]) if panel else set()
    nontruth = [
        c
        for c in cands
        if c["host_lineage"] == "0"
        and not in_panel(int(c["taxid"]))
        and int(c["taxid"]) not in ancestors
    ]
    row["n_fp_observation"] = sum(int(c["groups_called"]) >= 1 for c in nontruth)
    nontruth_ids = {c["taxid"] for c in nontruth}
    row["max_eligible_groups_nontruth"] = max(
        (int(c["groups_called"]) for c in eligible if c["taxid"] in nontruth_ids), default=0
    )
    if not panel:
        return row
    row["max_eligible_groups_panel"] = max(
        (int(c["groups_called"]) for c in eligible if in_panel(int(c["taxid"]))), default=0
    )
    genus = tree.genus(panel)
    if panel not in tree.parent:
        row["panel_reached"] = "absent_from_db"
    elif any(n for t, n in called.items() if in_panel(t)):
        row["panel_reached"] = "species"
    elif genus and called.get(genus, 0):
        row["panel_reached"] = "genus"
    else:
        row["panel_reached"] = "none"
    ranks = {"exact": 0, "ancestor": 0, "wrong_branch": 0}
    for c in cands:
        t, n = int(c["taxid"]), int(c["groups_called"])
        if c["host_lineage"] == "1" or not n:
            continue
        ranks["exact" if in_panel(t) else "ancestor" if t in ancestors else "wrong_branch"] += n
    row.update({f"rank_{k}": v for k, v in ranks.items()})
    p_called = sum(int(c["groups_called"]) for c in cands if in_panel(int(c["taxid"])))
    p_all = p_called + sum(int(c["groups_offcell"]) for c in cands if in_panel(int(c["taxid"])))
    row["panel_called_fraction"] = p_called / p_all if p_all else ""
    row.update(cell_metrics(out_dir, in_panel, truth))
    return row


def cell_metrics(out_dir: Path, in_panel, truth: dict[str, float]) -> dict:
    """Truth-cell recovery and Spearman over the arm's called cells (cells absent from the truth file = 0)."""
    import anndata as ad
    import numpy as np
    from scipy.stats import spearmanr

    adata = ad.read_h5ad(out_dir / "cell_taxon.h5ad")
    cols = [i for i, t in enumerate(adata.var_names) if in_panel(int(t))]
    groups = np.asarray(adata.X[:, cols].sum(axis=1)).ravel() if cols else np.zeros(adata.n_obs)
    mols = np.array([truth.get(b, 0.0) for b in adata.obs_names])
    hit = mols > 0
    out: dict = {
        "n_truth_cells": sum(v > 0 for v in truth.values()),
        "n_truth_called": int(hit.sum()),
    }
    if hit.any():
        out["truth_cell_fraction"] = float((groups[hit] > 0).mean())
        out["median_cell_recovery"] = float(np.median(groups[hit] / mols[hit]))
    if np.ptp(groups) and np.ptp(mols):
        out["spearman"] = float(spearmanr(groups, mols).statistic)
    return out


def load_truth(path: Path) -> dict[str, float]:
    from viralscan.scripts.cellcalling import _strip_suffix

    rows = tsv(path)
    if not rows or not {"barcode", "molecules"} <= set(rows[0]):
        sys.exit(f"{path}: needs columns barcode, molecules")
    return {_strip_suffix(r["barcode"]): float(r["molecules"]) for r in rows}


def smallest_k(arms: list[dict]) -> dict:
    """Rule 2 for one DB x configuration: the only candidate k is one above every negative's best eligible taxon."""
    neg = [a for a in arms if a["role"].startswith("negative")]
    pos = [a for a in arms if a["role"].startswith("positive")]
    if not neg or not pos or any(a["status"] != "ok" for a in arms):
        return {"result": "not_assessed", "reason": "missing or failed library"}
    k = max([1, *(a["max_eligible_groups_nontruth"] + 1 for a in neg)])
    short = [a["library"] for a in pos if a["max_eligible_groups_panel"] < k]
    if short:
        return {"result": "none", "k_negatives": k, "panel_below_k": short}
    return {"result": "k", "k": k}


def rule1(by_config: dict, k: int, fp_counts) -> dict:
    """Keep 0.05/4 unless a grid configuration ties or beats it on prioritized FP in every negative AND has a higher
    median per-cell recovery in every recovery positive."""
    base = by_config.get((*DEFAULT, 1))
    if not base or any(a["status"] != "ok" for a in base):
        return {"result": "not_assessed", "reason": "default arm missing or failed"}
    recovery = sorted(
        {a["library"] for a in base if a["role"].startswith("positive") and a.get("n_truth_called")}
    )
    capture_only = sorted(
        {a["library"] for a in base if a["role"].startswith("positive")} - set(recovery)
    )
    if not recovery:
        return {"result": "not_assessed", "reason": "no positive with truth cells"}

    def lib(arms, name):
        return next(a for a in arms if a["library"] == name)

    qualifying = []
    for (conf, hits, trim), arms in sorted(by_config.items()):
        if not trim or (conf, hits) == DEFAULT:
            continue
        if any(a["status"] != "ok" for a in arms):
            return {"result": "not_assessed", "reason": f"arm {conf}/{hits} missing or failed"}
        negs = [a["library"] for a in base if a["role"].startswith("negative")]
        if any(fp_counts(lib(arms, n), k) > fp_counts(lib(base, n), k) for n in negs):
            continue
        if all(
            lib(arms, p).get("median_cell_recovery", -1)
            > lib(base, p).get("median_cell_recovery", -1)
            for p in recovery
        ):
            fp = sum(fp_counts(lib(arms, n), k) for n in negs)
            qualifying.append((fp, -conf, -hits, conf, hits))
    result = {"k": k, "recovery_positives": recovery, "capture_only": capture_only}
    if not qualifying:
        return {**result, "result": "keep", "confidence": DEFAULT[0], "min_hit_groups": DEFAULT[1]}
    _, _, _, conf, hits = min(qualifying)
    return {**result, "result": "change", "confidence": conf, "min_hit_groups": hits}


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--panel-taxid", action="append", default=[], metavar="LIB=TAXID")
    p.add_argument("--truth", action="append", default=[], metavar="LIB=TSV")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)

    roles = {r["library"]: r["role"] for r in tsv(a.manifest)}
    panels = {k: int(v) for k, v in pairs(a.panel_taxid, "--panel-taxid").items()}
    truths = pairs(a.truth, "--truth")
    positives = sorted(lib for lib, role in roles.items() if role.startswith("positive"))
    missing = [lib for lib in positives if lib not in panels or lib not in truths]
    if missing:
        sys.exit(f"every positive needs --panel-taxid and --truth; missing for {missing}")
    truth = {lib: load_truth(Path(truths[lib])) for lib in positives}

    conds = tsv(a.run / "conditions.tsv")
    unknown = sorted({c["library"] for c in conds} - set(roles))
    if unknown:
        sys.exit(f"libraries without a manifest role: {unknown}")
    rows = [
        arm_metrics(
            c, a.run, roles[c["library"]], panels.get(c["library"]), truth.get(c["library"], {})
        )
        for c in conds
    ]
    offcell = []
    for c in conds:
        cand = a.run / "out" / c["cond_id"] / "candidates.tsv"
        if (a.run / "out" / c["cond_id"] / "discovery_complete.json").is_file():
            for r in tsv(cand):
                n_c, n_o = int(r["groups_called"]), int(r["groups_offcell"])
                offcell.append([c["cond_id"], r["taxid"], n_c, n_o, n_c / (n_c + n_o)])

    def eligible_fp(arm: dict, k: int) -> int:
        """Prioritized(k) false-positive taxa of one negative arm."""
        path = a.run / "out" / arm["cond_id"] / "candidates.tsv"
        return sum(
            r["host_lineage"] == "0" and r["viral_only_db"] == "0" and int(r["groups_called"]) >= k
            for r in tsv(path)
        )

    decisions: dict = {"rule2": {}, "rule1": {}}
    for db in sorted({r["db_name"] for r in rows}):
        by_config: dict = {}
        for r in rows:
            if r["db_name"] == db:
                key = (r["confidence"], r["min_hit_groups"], r["polya_trim"])
                by_config.setdefault(key, []).append(r)
        r2 = {
            f"{c}/{h}{'' if t else '/notrim'}": smallest_k(arms)
            for (c, h, t), arms in by_config.items()
        }
        decisions["rule2"][db] = r2
        default_k = r2.get(f"{DEFAULT[0]}/{DEFAULT[1]}", {})
        if default_k.get("result") != "k":
            decisions["rule1"][db] = {
                "result": "not_assessed",
                "reason": f"no rule-2 k at the default ({default_k.get('result', 'missing')})",
            }
        else:
            decisions["rule1"][db] = rule1(by_config, default_k["k"], eligible_fp)
    decisions["rule3"] = {"result": "not_assessed", "reason": "spike-in only (protocol §5.3)"}
    decisions["rule4"] = {
        "result": "reported",
        "note": "DBs side by side; viral-only DB capped at observation (G15)",
    }
    decisions["rule5"] = "any rule above that could not be evaluated says not_assessed"

    a.out.mkdir(parents=True, exist_ok=True)
    with (a.out / "metrics.tsv").open("w", newline="") as handle:
        w = csv.DictWriter(handle, METRIC_COLUMNS, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    with (a.out / "offcell.tsv").open("w", newline="") as handle:
        w2 = csv.writer(handle, delimiter="\t")
        w2.writerow(["cond_id", "taxid", "groups_called", "groups_offcell", "called_fraction"])
        w2.writerows(offcell)
    (a.out / "decisions.json").write_text(json.dumps(decisions, indent=1, sort_keys=True) + "\n")
    print(f"{len(rows)} arms -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

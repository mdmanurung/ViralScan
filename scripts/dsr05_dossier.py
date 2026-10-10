#!/usr/bin/env python3
"""DOSSIER-02: check a DSR round against the per-dataset evidence dossier spec (read-only).

    python scripts/dsr05_dossier.py <round_dir> [--out <round_dir>/dossier] [--no-derived]

Reads `analysis/dsr_round1/dossier_spec.tsv` (what must exist) and `dossier_roles.tsv` (role and acceptance per
sample), walks the round and writes, under --out:

  index.tsv              one row per (item, unit): status ok / missing / stale / failed / n_a / not_started
  <dataset>__<sample>.md per-sample checklist, the highest claim rung reachable, and what blocks the next one
  round_index.md         one line per sample
  cell_calling_summary.tsv, recurrence.tsv, arm_concordance.tsv   derived from files already on disk

Nothing is regenerated and no run directory is touched. Claim ladder (see the plan): C0 cells called, C1 virus
UMI counts, C2 read-level support, C3 cell-type claim, N1 "no panel virus above the LOD". An item gates a rung
only on the sample's reference arm, `combined_off`.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SPEC = REPO / "analysis" / "dsr_round1" / "dossier_spec.tsv"
ROLES = REPO / "analysis" / "dsr_round1" / "dossier_roles.tsv"
REFERENCE_ARM = "combined_off"
ARMS = ("combined_off", "combined_artefact", "twostep_v2")
RUNGS = ("C0", "C1", "C2", "C3")
MIN_CALL_MOLECULES = 3
HASH_LIMIT = 20_000_000  # bytes; larger files get a size only
INDEX_COLS = [
    "item_id",
    "section",
    "unit",
    "arm",
    "path",
    "status",
    "sha256",
    "size",
    "mtime",
    "reason",
]


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def slug(virus: str) -> str:
    """Same rule as scripts/slurm_dsr02_evidence_array.sh."""
    return re.sub(r"[^A-Za-z0-9]+", "_", virus).strip("_")


def inner_name(root: Path, sample: str) -> str | None:
    """Sample directory inside a run root: usually the sample, but `x223` holds `LUM-SJ-x223`."""
    if (root / sample).is_dir():
        return sample
    subs = (
        [p.name for p in root.iterdir() if p.is_dir() and (p / "log").is_dir()]
        if root.is_dir()
        else []
    )
    return subs[0] if len(subs) == 1 else None


def not_applicable(round_dir: Path) -> dict[tuple[str, str, str], str]:
    path = round_dir / "status" / "n_a.tsv"
    return (
        {(r["dataset"], r["sample"], r["arm"]): r["reason"] for r in _tsv(path)}
        if path.is_file()
        else {}
    )


def calls_by_unit(round_dir: Path) -> dict[tuple[str, str, str], list[dict[str, str]]]:
    """Calls per (dataset, sample, arm), enumerated from the runs themselves, not from a stale calls.tsv.

    Same rule as `scripts/dsr02_enumerate_calls.py` (>= 3 molecules, or any anellovirus call). The sample is the
    run-root name (`x223`), not the inner directory (`LUM-SJ-x223`).
    """
    out: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for r in load_script("dsr02_enumerate_calls").calls(round_dir, MIN_CALL_MOLECULES):
        sample = Path(r["run_dir"]).parent.name
        if r["arm"] in ARMS:
            out[r["dataset"], sample, r["arm"]].append({k: str(v) for k, v in r.items()})
    return out


class Ctx:
    """Placeholders for one (sample, arm, call) unit."""

    def __init__(
        self,
        round_dir: Path,
        ds: str,
        sample: str,
        arm: str,
        role: dict[str, str],
        call: dict | None = None,
    ):
        self.round_dir, self.ds, self.sample, self.arm, self.role, self.call = (
            round_dir,
            ds,
            sample,
            arm,
            role,
            call,
        )
        self.root = Path("runs") / ds / arm / sample
        inner = inner_name(round_dir / self.root, sample)
        self.inner_name = inner or sample
        self.inner = self.root / self.inner_name

    def expand(self, pattern: str) -> list[str]:
        out = []
        evs = []
        if self.call is not None:
            d = f"{self.ds}__{self.sample}__{self.arm}__{slug(self.call['virus'])}"
            evs = [f"evidence_v2/{d}", f"evidence/{d}"]
        for alt in pattern.split("|"):
            for ev in evs or [""]:
                out.append(
                    alt.replace("{ev}", ev)
                    .replace("{root}", str(self.root))
                    .replace("{inner}", str(self.inner))
                    .replace("{ds}", self.ds)
                    .replace("{sample}", self.sample)
                    .replace("{arm}", self.arm)
                    .replace("{vendor_barcodes}", self.role.get("vendor_barcodes", ""))
                    .replace("{cell_labels}", self.role.get("cell_labels", ""))
                )
        return out

    def first_existing(self, pattern: str) -> tuple[Path | None, Path]:
        paths = [
            self.round_dir / p if not p.startswith("/") else Path(p)
            for p in self.expand(pattern)
            if p
        ]
        for p in paths:
            if p.exists():
                return p, p
        return None, paths[0] if paths else self.round_dir


def _json_get(path: Path, dotted: str):
    node = json.loads(path.read_text())
    for key in dotted.split("."):
        node = node.get(key) if isinstance(node, dict) else None
    return node


def file_check(found: Path, check: str) -> tuple[str, str]:
    """Status and reason for an existing path under `check`."""
    if found.is_file() and found.stat().st_size == 0:
        return "missing", "file is empty"
    if check.startswith("json:"):
        value = _json_get(found, check[5:])
        return ("ok", "") if value not in (None, "", {}, []) else ("missing", f"no {check[5:]}")
    if check.startswith("json_ne:"):
        key, bad = check[8:].split("=", 1)
        value = _json_get(found, key)
        return (
            ("ok", "")
            if value not in (None, bad)
            else ("missing", f"{key} is {value!r}: not certified")
        )
    if check.startswith("contains:"):
        return ("ok", "") if check[9:] in found.read_text() else ("missing", f"no {check[9:]}")
    if check == "anello_align":
        with found.open(newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        filled = any((r.get("alignment_status") or "") not in ("", "not_run") for r in rows)
        return (
            ("ok", "") if filled else ("n_a", "alignment_* columns empty: --anello-align not run")
        )
    return "ok", ""


def newest_mtime(round_dir: Path, dep: str, ctx: Ctx) -> float:
    """Newest mtime among the files an item is derived from (0 when none exist)."""
    if dep == "{arm_run_complete}":
        paths = [round_dir / ctx.root / "run_complete.json"]
    else:  # {twostep_run_complete}
        paths = list((round_dir / "runs").glob("*/twostep_v2/*/run_complete.json"))
    return max((p.stat().st_mtime for p in paths if p.is_file()), default=0.0)


def sha256_of(path: Path) -> str:
    if not path.is_file() or path.stat().st_size > HASH_LIMIT:
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── derived tables (from files already on disk) ───────────────────────────────────────────────────
def _called(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    return {ln.strip() for ln in path.open() if ln.strip() and ln.strip() != "barcode"}


def cell_calling_row(
    round_dir: Path, ds: str, sample: str, arm: str, ctx: Ctx
) -> dict[str, object] | None:
    drops = round_dir / ctx.inner / "kb-python" / "counts_unfiltered" / "emptydrops_cells.tsv"
    called = _called(round_dir / ctx.inner / "results" / "called_cells.tsv")
    host_called = _called(round_dir / ctx.inner / "results" / "host_called_cells.tsv")
    if not drops.is_file() and not called and not host_called:
        return None
    tested = n_cell = 0
    knee = inflection = ""
    totals: list[int] = []
    if drops.is_file():
        with drops.open(newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                tested += 1
                knee, inflection = r.get("knee", ""), r.get("inflection", "")
                if r.get("is_cell") in ("TRUE", "True", "true", "1"):
                    n_cell += 1
                    totals.append(int(float(r["total"])))
    totals.sort()
    return {
        "dataset": ds,
        "sample": sample,
        "arm": arm,
        "barcodes_tested": tested,
        "emptydrops_is_cell": n_cell,
        "called_cells_file": len(called),
        "host_called_cells_file": len(host_called) if host_called else "",
        "knee": knee,
        "inflection": inflection,
        "median_umi_called": totals[len(totals) // 2] if totals else "",
        "_called": called,
        "_host": host_called,
    }


def write_derived(
    round_dir: Path, out: Path, samples: list[tuple[str, str]], common_rows: list[dict[str, str]]
):
    cc_rows: list[dict[str, object]] = []
    for ds, sample in samples:
        per_arm = {}
        for arm in ARMS:
            row = cell_calling_row(round_dir, ds, sample, arm, Ctx(round_dir, ds, sample, arm, {}))
            if row:
                per_arm[arm] = row
        ref = per_arm.get(REFERENCE_ARM)
        for row in per_arm.values():
            other = row["_host"] or row["_called"]
            ref_set = ref["_called"] if ref else set()
            union = len(other | ref_set)
            row["jaccard_with_combined_off"] = (
                f"{len(other & ref_set) / union:.3f}" if union and ref else ""
            )
            cc_rows.append({k: v for k, v in row.items() if not k.startswith("_")})
    cols = [
        "dataset",
        "sample",
        "arm",
        "barcodes_tested",
        "emptydrops_is_cell",
        "called_cells_file",
        "host_called_cells_file",
        "knee",
        "inflection",
        "median_umi_called",
        "jaccard_with_combined_off",
    ]
    _write(out / "cell_calling_summary.tsv", cc_rows, cols)

    ref_called: dict[str, set[str]] = defaultdict(
        set
    )  # virus -> samples where combined_off calls it
    best: dict[str, int] = defaultdict(int)
    for r in common_rows:
        if r["arm"] == REFERENCE_ARM and int(r["infected_in_reference"]) > 0:
            ref_called[r["virus"]].add(f"{r['dataset']}/{r['sample']}")
            best[r["virus"]] = max(best[r["virus"]], int(r["infected_in_reference"]))
    rec = [
        {
            "virus": v,
            "n_samples": len(s),
            "n_datasets": len({x.split("/")[0] for x in s}),
            "max_cells_in_one_sample": best[v],
            "samples": ",".join(sorted(s)),
        }
        for v, s in ref_called.items()
    ]
    rec.sort(key=lambda r: (-int(r["n_datasets"]), -int(r["n_samples"])))
    _write(
        out / "recurrence.tsv",
        rec,
        ["virus", "n_samples", "n_datasets", "max_cells_in_one_sample", "samples"],
    )

    pivot: dict[tuple[str, str, str], dict[str, str]] = defaultdict(dict)
    for r in common_rows:
        pivot[r["dataset"], r["sample"], r["virus"]][r["arm"]] = r["infected_in_reference"]
    conc = [
        {"dataset": d, "sample": s, "virus": v, **{a: arms.get(a, "") for a in ARMS}}
        for (d, s, v), arms in sorted(pivot.items())
    ]
    _write(out / "arm_concordance.tsv", conc, ["dataset", "sample", "virus", *ARMS])
    return {r["sample"] for r in cc_rows}, {r["sample"] for r in conc}


def _write(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


# ── the check ─────────────────────────────────────────────────────────────────────────────────────
def evaluate(round_dir: Path, out: Path, derived: bool = True) -> tuple[list[dict], dict]:
    spec = _tsv(SPEC)
    roles = _tsv(ROLES)
    na = not_applicable(round_dir)
    calls = calls_by_unit(round_dir)
    verdicts: set[str] = set()
    for name in ("evidence_v2", "evidence"):
        p = round_dir / name / "verdicts.tsv"
        if p.is_file():
            verdicts = {r["evidence_dir"] for r in _tsv(p)}
            break

    run_samples = [(r["dataset"], r["sample"]) for r in roles if r["role"] != "not_started"]
    derived_ok: dict[str, set[str]] = {"cell_calling_summary": set(), "arm_concordance": set()}
    if derived:
        common = load_script("dsr_common_cells")
        common_rows = [{k: str(v) for k, v in r.items()} for r in common.rows(round_dir)]
        derived_ok["cell_calling_summary"], derived_ok["arm_concordance"] = write_derived(
            round_dir, out, run_samples, common_rows
        )
        derived_ok["recurrence"] = {"round"} if (out / "recurrence.tsv").is_file() else set()
    else:
        for name in ("cell_calling_summary", "arm_concordance", "recurrence"):
            p = out / f"{name}.tsv"
            derived_ok[name] = (
                ({r.get("sample", "round") for r in _tsv(p)} | {"round"}) if p.is_file() else set()
            )

    index: list[dict] = []
    for role in roles:
        ds, sample = role["dataset"], role["sample"]
        if role["role"] == "not_started":
            index.append(
                {
                    "item_id": "-",
                    "section": "-",
                    "unit": f"{ds}",
                    "arm": "",
                    "path": "",
                    "status": "not_started",
                    "sha256": "",
                    "size": "",
                    "mtime": "",
                    "reason": role["acceptance"],
                }
            )
            continue
        sample_na = na.get((ds, sample, "all"))
        for item in spec:
            scope = item["scope"]
            arms = [a for a in item["arms"].split(",") if a] or [REFERENCE_ARM]
            if scope == "arm":
                units = [(a, None) for a in arms]
            elif scope == "call":
                units = [(a, c) for a in ARMS for c in calls.get((ds, sample, a), [])]
            else:
                units = [(REFERENCE_ARM, None)]
            for arm, call in units:
                ctx = Ctx(round_dir, ds, sample, arm, role, call)
                row = {
                    "item_id": item["item_id"],
                    "section": item["section"],
                    "unit": f"{ds}/{sample}" + (f"/{call['virus']}" if call else ""),
                    "arm": arm if scope in ("arm", "call") else "",
                    "path": "",
                    "status": "",
                    "sha256": "",
                    "size": "",
                    "mtime": "",
                    "reason": item["gap_note"],
                }
                index.append(row)
                arm_na = na.get((ds, sample, arm)) if scope in ("arm", "call") else None
                if sample_na or arm_na:
                    row["status"], row["reason"] = "n_a", sample_na or arm_na
                    continue
                check = item["check"]
                if check.startswith("derived:"):
                    name = check[8:]
                    ok = (sample in derived_ok[name]) or (
                        name == "recurrence" and derived_ok.get(name)
                    )
                    row["path"] = f"dossier/{name}.tsv"
                    row["status"] = "ok" if ok else "missing"
                    continue
                if check == "no_producer":
                    row["status"] = "missing"
                    continue
                if check == "vendor" and not role["vendor_barcodes"]:
                    row["status"], row["reason"] = "n_a", "no vendor barcode set for this sample"
                    continue
                if check == "labels" and not role["cell_labels"]:
                    row["status"], row["reason"] = (
                        "missing",
                        "no cell-type or donor labels recorded for this sample",
                    )
                    continue
                if check == "role_row":
                    row["status"], row["path"] = "ok", "analysis/dsr_round1/dossier_roles.tsv"
                    continue
                if check == "findings":
                    ids = [f for f in role["findings"].split(",") if f]
                    texts = [
                        p.read_text(errors="ignore")
                        for p in (REPO / ".living" / "findings").glob("*.md")
                    ]
                    have = [i for i in ids if any(i in t for t in texts)]
                    row["status"] = "ok" if ids and len(have) == len(ids) else "missing"
                    row["reason"] = (
                        f"findings {','.join(have) or 'none'} of {role['findings'] or 'none recorded'}"
                    )
                    continue
                if check == "verdict_row":
                    d = f"{ds}__{sample}__{arm}__{slug(call['virus'])}"
                    row["path"], row["status"] = (
                        "verdicts.tsv",
                        "ok" if d in verdicts else "missing",
                    )
                    continue
                found, shown = ctx.first_existing(item["path_pattern"])
                row["path"] = (
                    str(shown.relative_to(round_dir))
                    if shown.is_relative_to(round_dir)
                    else str(shown)
                )
                if found is None:
                    arm_dir = round_dir / ctx.root
                    failed = (
                        scope == "arm"
                        and arm_dir.is_dir()
                        and not (arm_dir / "run_complete.json").is_file()
                    )
                    row["status"] = "failed" if failed else "missing"
                    continue
                status, reason = file_check(found, check)
                row["status"], row["reason"] = (
                    status,
                    reason or ("" if status == "ok" else row["reason"]),
                )
                st = found.stat()
                row["size"], row["mtime"], row["sha256"] = (
                    st.st_size,
                    int(st.st_mtime),
                    sha256_of(found),
                )
                dep = item["stale_if_older_than"]
                if status == "ok" and dep and st.st_mtime < newest_mtime(round_dir, dep, ctx):
                    row["status"], row["reason"] = "stale", f"older than {dep.strip('{}')}"
    return index, {"roles": roles, "calls": calls}


def rung_table(
    index: list[dict], spec: list[dict], ds: str, sample: str, calls: dict
) -> tuple[dict[str, bool], dict[str, list[str]], str]:
    gates = {s["item_id"]: set(g for g in s["gates"].split(",") if g) for s in spec}
    unit = f"{ds}/{sample}"
    blockers: dict[str, list[str]] = defaultdict(list)
    for r in index:
        if not r["unit"].startswith(unit) or r["item_id"] == "-":
            continue
        if r["arm"] not in ("", REFERENCE_ARM):
            continue
        if r["status"] in ("ok", "n_a"):
            continue
        for g in gates.get(r["item_id"], ()):
            blockers[g].append(
                f"{r['item_id']} {r['status']}: {r['unit'].removeprefix(unit).lstrip('/') or sample}"
            )
    d01 = [
        r["status"]
        for r in index
        if r["unit"] == unit and r["arm"] == REFERENCE_ARM and r["item_id"] == "D01"
    ]
    if d01 and all(x == "n_a" for x in d01):
        return dict.fromkeys((*RUNGS, "N1"), False), blockers, "n_a"
    has_ref = "ok" in d01
    n_calls = len(
        {
            r["unit"]
            for r in index
            if r["section"] == "E"
            and r["unit"].startswith(unit + "/")
            and r["arm"] == REFERENCE_ARM
        }
    )
    if not n_calls:  # C2 would be vacuously true: nothing was called, so nothing was read-validated
        blockers["C2"].append("no combined_off call to read-validate")
    reached = {"C0": not blockers["C0"], "C1": False, "C2": False, "C3": False, "N1": False}
    reached["C1"] = reached["C0"] and not blockers["C1"]
    reached["C2"] = reached["C1"] and not blockers["C2"]
    reached["C3"] = reached["C2"] and not blockers["C3"]
    reached["N1"] = reached["C1"] and not blockers["N1"]
    top = (
        next((g for g in reversed(RUNGS) if reached[g]), "none")
        if has_ref or reached["C0"]
        else "none"
    )
    return reached, blockers, top


def write_reports(
    round_dir: Path, out: Path, index: list[dict], ctx: dict
) -> list[tuple[str, str, str, str]]:
    spec = _tsv(SPEC)
    sym = {
        "ok": "ok",
        "missing": "MISSING",
        "stale": "STALE",
        "failed": "FAILED",
        "n_a": "n/a",
        "not_started": "not started",
    }
    summary = []
    for role in ctx["roles"]:
        ds, sample = role["dataset"], role["sample"]
        if role["role"] == "not_started":
            summary.append((ds, sample, "not_started", role["acceptance"]))
            continue
        reached, blockers, top = rung_table(index, spec, ds, sample, ctx["calls"])
        rows = [r for r in index if r["unit"].startswith(f"{ds}/{sample}")]
        n = defaultdict(int)
        for r in rows:
            n[r["status"]] += 1
        lines = [
            f"# {ds} / {sample}",
            "",
            f"Role: **{role['role']}**. {role['acceptance']}",
            "",
            f"Highest claim rung reachable now: **{top}**. N1 (no panel virus above the LOD): "
            f"**{'informative' if reached['N1'] else 'uninformative'}**.",
            "",
            "| Rung | Reached | Blocked by |",
            "|---|---|---|",
        ]
        for g in (*RUNGS, "N1"):
            b = blockers[g]
            lines.append(
                f"| {g} | {'yes' if reached[g] else 'no'} | {'; '.join(b[:4])}{' ...' if len(b) > 4 else ''} |"
            )
        lines += [
            "",
            "## Checklist",
            "",
            "| Item | Unit | Arm | Status | Note |",
            "|---|---|---|---|---|",
        ]
        for r in rows:
            note = (r["reason"] or "")[:90]
            lines.append(
                f"| {r['item_id']} | {r['unit'].removeprefix(f'{ds}/{sample}').lstrip('/') or '-'} | {r['arm']} | {sym[r['status']]} | {note} |"
            )
        (out / f"{ds}__{sample}.md").write_text("\n".join(lines) + "\n")
        summary.append((ds, sample, top, ", ".join(f"{k} {v}" for k, v in sorted(n.items()))))
    lines = [
        "# DSR round dossier",
        "",
        "| Dataset | Sample | Highest rung | Item status counts |",
        "|---|---|---|---|",
    ]
    lines += [f"| {d} | {s} | {t} | {c} |" for d, s, t, c in summary]
    (out / "round_index.md").write_text("\n".join(lines) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round_dir", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--no-derived", action="store_true", help="skip rebuilding the derived tables")
    args = ap.parse_args(argv)
    round_dir = args.round_dir.resolve()
    out = (args.out or round_dir / "dossier").resolve()
    out.mkdir(parents=True, exist_ok=True)
    index, ctx = evaluate(round_dir, out, derived=not args.no_derived)
    _write(out / "index.tsv", index, INDEX_COLS)
    summary = write_reports(round_dir, out, index, ctx)
    counts: dict[str, int] = defaultdict(int)
    for r in index:
        counts[r["status"]] += 1
    print(
        f"{len(index)} rows: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())),
        file=sys.stderr,
    )
    for ds, sample, top, _ in summary:
        print(f"{ds}\t{sample}\t{top}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

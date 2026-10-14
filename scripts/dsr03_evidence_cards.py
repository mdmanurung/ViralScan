#!/usr/bin/env python3
"""CARD-01: one tiered evidence card per (dataset, sample, virus) call on the reference arm.

    python scripts/dsr03_evidence_cards.py <round_dir> [--calibration calibration.json] \\
        [--labels DATASET/SAMPLE=labels.tsv ...] [--n-perm 1000] [--seed 0] [--out <round_dir>/cards]

Writes `<out>/evidence_cards.tsv` (one row per call) and `<out>/evidence_cards.md`. A call is what
`scripts/dsr02_enumerate_calls.py` lists (rows with `detected` false are skipped there) on `combined_off`, the
arm that gates claims. Tiers, each reported, none a conclusion by itself:

  read_level        verdict of `scripts/dsr02_verdicts.py` over counted molecules of the FROZEN reference cells
                    (`reference_cells/`), `pass` only for viral_best with `molecule_scope` = restricted
  genome_structure  breadth and read-start concentration of the best viral reference (alignment_qc.tsv) and the
                    coverage_hotspot / eve_or_integration_like_hotspot flags; numbers only, no cutoff
  ambient           Poisson upper-tail p of the viral molecules observed in the reference cells, given the
                    viral share of the barcodes OUTSIDE them (kb `adata.h5ad`, unique-assignment counts). The
                    expectation scales the pool's viral share by the cells' total UMIs, which over-counts
                    ambient (cell UMIs are mostly endogenous), so p is conservative. p only, no pass/fail
  recurrence        datasets / samples that call the virus (`dossier/recurrence.tsv`; run dsr05_dossier first)
  permutation       cell-type enrichment of infected cells: label permutation within the reference cells,
                    needs `--labels` (barcode, label); `n_a` without labels

`card_status` follows the interim rule (user, 2026-10-07): a virus the dataset is not meant to contain
(`expected` false) is `unverified` until a calibration record (`scripts/dsr02_calibration.py --out`) says the
gate passed under the rule values now in force, whatever its read-level tier shows. Otherwise: `no_evidence`,
`not_restricted` (no frozen cells or no molecule table), `read_level_fail` or `read_level_pass`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import random
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "src") not in sys.path:  # run from a checkout without `pip install -e .`
    sys.path.insert(0, str(REPO / "src"))

from viralscan.sensitivity import poisson_detection_probability  # noqa: E402
from viralscan.virus_grouping import group_genes_by_identity, load_run_identity  # noqa: E402

REFERENCE_ARM = "combined_off"
EVIDENCE_DIRS = ("evidence_v2", "evidence")
CARD_COLS = [
    "dataset", "sample", "virus", "expected", "molecules", "card_status", "status_reason",
    "read_level", "read_verdict", "molecule_scope", "virus_best_molecules", "host_best_molecules",
    "tie_molecules", "reads_telomeric", "reads_tso", "reads_low_complexity", "reads_total",
    "molecules_telomeric", "molecules_tso", "molecules_low_complexity", "molecules_total",
    "genome_structure", "n_virus_references", "best_breadth_1x", "max_50bp_window_fraction",
    "coverage_hotspot", "eve_or_integration_like_hotspot",
    "ambient", "observed_in_cells", "expected_ambient_in_cells", "ambient_poisson_p", "ambient_note",
    "recurrence", "recurrence_samples", "recurrence_datasets",
    "permutation", "perm_top_label", "perm_top_label_p_bonferroni", "permutation_note",
]  # fmt: skip


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


verdicts = _load("dsr02_verdicts")
calibration = _load("dsr02_calibration")
enumerate_calls = _load("dsr02_enumerate_calls")
dossier = _load("dsr05_dossier")


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def evidence_dir(round_dir: Path, dataset: str, sample: str, virus: str) -> Path | None:
    name = f"{dataset}__{sample}__{REFERENCE_ARM}__{dossier.slug(virus)}"
    for root in EVIDENCE_DIRS:
        d = round_dir / root / name
        if (d / "interpretation_flags.tsv").is_file():
            return d
    return None


def read_level(ev: Path | None, cells: set[str] | None) -> dict[str, object]:
    if ev is None:
        return {"read_level": "no_evidence"}
    v = verdicts.verdict(ev, cells)
    if cells is None:
        v["molecule_scope"] = "no_frozen_cells"
    out: dict[str, object] = {
        "read_verdict": v["verdict"],
        "molecule_scope": v["molecule_scope"],
        "virus_best_molecules": v["virus_best_molecules"],
        "host_best_molecules": v["host_best_molecules"],
        "tie_molecules": v["tie_molecules"],
    }
    summary = ev / "read_verdict_summary.tsv"
    if summary.is_file():
        per: dict[str, int] = defaultdict(int)
        for r in _tsv(summary):
            per[r["repeat_class"]] += int(r["reads"])
        out.update(
            reads_telomeric=per["telomeric"],
            reads_tso=per["tso"],
            reads_low_complexity=per["low_complexity"],
            reads_total=sum(per.values()),
        )
    molecular_summary = ev / "molecule_read_verdict_summary.tsv"
    if molecular_summary.is_file():
        per_molecule: dict[str, int] = defaultdict(int)
        for r in _tsv(molecular_summary):
            per_molecule[r["repeat_class"]] += int(r["molecules"])
        out.update(
            molecules_telomeric=per_molecule["telomeric"],
            molecules_tso=per_molecule["tso"],
            molecules_low_complexity=per_molecule["low_complexity"],
            molecules_total=sum(per_molecule.values()),
        )
    if v["molecule_scope"] != "restricted":
        out["read_level"] = "not_restricted"
    else:
        out["read_level"] = "pass" if v["verdict"] == "viral_best" else "fail"
    return out


def genome_structure(ev: Path | None) -> dict[str, object]:
    if ev is None or not (ev / "alignment_qc.tsv").is_file():
        return {"genome_structure": "no_evidence"}
    qc = [
        r
        for r in _tsv(ev / "alignment_qc.tsv")
        if r.get("count_layer") == "deduplicated" and r.get("reference_class") == "virus"
    ]
    if not qc:
        return {"genome_structure": "no_viral_alignment"}
    flags = {r["flag"]: r["status"] for r in _tsv(ev / "interpretation_flags.tsv")}
    return {
        "genome_structure": "reported",
        "n_virus_references": len(qc),
        "best_breadth_1x": max(float(r.get("breadth_1x") or 0) for r in qc),
        "max_50bp_window_fraction": max(float(r.get("max_50bp_window_fraction") or 0) for r in qc),
        "coverage_hotspot": flags.get("coverage_hotspot", ""),
        "eve_or_integration_like_hotspot": flags.get("eve_or_integration_like_hotspot", ""),
    }


def ambient(inner: Path, cells: set[str] | None, virus: str) -> dict[str, object]:
    """Poisson upper-tail p of the observed viral molecules in the cells under the outside-barcode rate."""
    if cells is None:
        return {"ambient": "n_a", "ambient_note": "no frozen reference cells"}
    adata_path = inner / "kb-python" / "counts_unfiltered" / "adata.h5ad"
    identity = load_run_identity(inner)
    if not adata_path.is_file() or identity is None:
        return {"ambient": "n_a", "ambient_note": "no adata.h5ad or virus_identity.tsv in the run"}
    import anndata  # heavy; only the ambient tier needs it
    import numpy as np

    adata = anndata.read_h5ad(adata_path)
    genes, _ = group_genes_by_identity(list(adata.var_names), identity)
    cols = np.flatnonzero(adata.var_names.isin(genes.get(virus, [])))
    if not len(cols):
        return {"ambient": "n_a", "ambient_note": "no gene of this virus in adata.h5ad"}
    x = adata.X.tocsr() if hasattr(adata.X, "tocsr") else np.asarray(adata.X)
    total = np.asarray(x.sum(axis=1)).ravel()
    viral = np.asarray(x[:, cols].sum(axis=1)).ravel()
    in_cells = np.asarray(adata.obs_names.isin(cells))
    pool_umis = float(total[~in_cells].sum())
    if not in_cells.any() or pool_umis <= 0:
        return {"ambient": "n_a", "ambient_note": "no cells or no outside-barcode UMIs"}
    expected = float(viral[~in_cells].sum()) / pool_umis * float(total[in_cells].sum())
    observed = float(viral[in_cells].sum())
    return {
        "ambient": "reported",
        "observed_in_cells": observed,
        "expected_ambient_in_cells": round(expected, 4),
        "ambient_poisson_p": f"{poisson_detection_probability(expected, int(observed)):.3g}",
    }


def recurrence(round_dir: Path, virus: str) -> dict[str, object]:
    path = round_dir / "dossier" / "recurrence.tsv"
    if not path.is_file():
        return {"recurrence": "n_a"}
    for r in _tsv(path):
        if r["virus"] == virus:
            return {
                "recurrence": "reported",
                "recurrence_samples": r["n_samples"],
                "recurrence_datasets": r["n_datasets"],
            }
    return {"recurrence": "reported", "recurrence_samples": 0, "recurrence_datasets": 0}


def permutation(
    round_dir: Path,
    dataset: str,
    sample: str,
    virus: str,
    cells: set[str] | None,
    labels_path: Path | None,
    n_perm: int,
    seed: int,
) -> dict[str, object]:
    """Which label holds more infected cells than shuffled labels would; Bonferroni over the labels."""
    if labels_path is None:
        return {"permutation": "n_a", "permutation_note": "no cell-type labels given"}
    cell_level = round_dir / "cell_level" / f"{dataset}__{sample}.tsv.gz"
    if cells is None or not cell_level.is_file():
        return {"permutation": "n_a", "permutation_note": "no frozen cells or cell_level table"}
    labels = {r["barcode"]: r["label"] for r in _tsv(labels_path)}
    universe = sorted(b for b in cells if b in labels)
    with gzip.open(cell_level, "rt", newline="") as fh:
        infected = {
            r["barcode"]
            for r in csv.DictReader(fh, delimiter="\t")
            if r["arm"] == REFERENCE_ARM and r["virus_name"] == virus
        }
    infected &= set(universe)
    names = sorted({labels[b] for b in universe})
    if not infected or len(names) < 2:
        return {"permutation": "n_a", "permutation_note": "no labelled infected cells or one label"}
    pool = [labels[b] for b in universe]
    obs = {n: sum(labels[b] == n for b in infected) for n in names}
    rng = random.Random(seed)
    ge = dict.fromkeys(names, 0)
    for _ in range(n_perm):
        shuffled = rng.sample(pool, len(pool))
        hit = dict.fromkeys(names, 0)
        for lab in shuffled[: len(infected)]:  # a random same-size set of cells
            hit[lab] += 1
        for n in names:
            ge[n] += hit[n] >= obs[n]
    p = {n: (1 + ge[n]) / (1 + n_perm) for n in names}
    top = min(names, key=lambda n: (p[n], n))
    return {
        "permutation": "reported",
        "perm_top_label": top,
        "perm_top_label_p_bonferroni": f"{min(1.0, p[top] * len(names)):.3g}",
        "permutation_note": f"{len(infected)} infected of {len(universe)} labelled cells, {n_perm} permutations",
    }


def card_status(call: dict, tiers: dict[str, object], gate: tuple[bool, str]) -> tuple[str, str]:
    """Interim rule first: an unexpected virus stays unverified until the calibration gate has passed."""
    if call["expected"] != "true" and not gate[0]:
        return "unverified", f"not an expected virus and {gate[1]}"
    level = tiers["read_level"]
    if level == "no_evidence":
        return "no_evidence", "no evidence directory for this call"
    if level == "not_restricted":
        return "not_restricted", f"molecule_scope={tiers['molecule_scope']}"
    return ("read_level_pass" if level == "pass" else "read_level_fail"), ""


def build_cards(
    round_dir: Path,
    *,
    calibration_record: Path | None = None,
    labels: dict[str, Path] | None = None,
    n_perm: int = 1000,
    seed: int = 0,
) -> list[dict[str, object]]:
    gate = calibration.gate_status(calibration_record)
    cards = []
    for call in enumerate_calls.calls(round_dir):  # skips detected = false
        if call["arm"] != REFERENCE_ARM:
            continue
        inner = Path(call["run_dir"])
        sample = inner.parent.name  # run-root name; the inner directory can differ (x223)
        cells = _frozen(round_dir, call["dataset"], sample)
        ev = evidence_dir(round_dir, call["dataset"], sample, call["virus"])
        tiers: dict[str, object] = {
            **read_level(ev, cells),
            **genome_structure(ev),
            **ambient(inner, cells, call["virus"]),
            **recurrence(round_dir, call["virus"]),
            **permutation(
                round_dir,
                call["dataset"],
                sample,
                call["virus"],
                cells,
                (labels or {}).get(f"{call['dataset']}/{sample}"),
                n_perm,
                seed,
            ),
        }
        status, reason = card_status(call, tiers, gate)
        cards.append(
            {
                "dataset": call["dataset"],
                "sample": sample,
                "virus": call["virus"],
                "expected": call["expected"],
                "molecules": call["molecules"],
                "card_status": status,
                "status_reason": reason,
                **tiers,
            }
        )
    return cards


def _frozen(round_dir: Path, dataset: str, sample: str) -> set[str] | None:
    path = round_dir / "reference_cells" / f"{dataset}__{sample}.tsv"
    if not path.is_file():
        return None
    return {ln.strip() for ln in path.read_text().splitlines()[1:] if ln.strip()}


def write_cards(cards: list[dict[str, object]], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with (out / "evidence_cards.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CARD_COLS, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(cards)
    lines = ["# Evidence cards", ""]
    tiers = ("read_level", "genome_structure", "ambient", "recurrence", "permutation")
    for c in cards:
        lines += [f"## {c['dataset']} / {c['sample']} / {c['virus']}", ""]
        lines += [f"**{c['card_status']}**" + (f": {c['status_reason']}" if c["status_reason"] else ""), ""]
        lines += [f"- {t}: {c.get(t, '')}" for t in tiers] + [""]
    (out / "evidence_cards.md").write_text("\n".join(lines))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round_dir", type=Path)
    ap.add_argument("--calibration", type=Path, help="record from dsr02_calibration.py --out")
    ap.add_argument("--labels", action="append", default=[], metavar="DATASET/SAMPLE=TSV")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args(argv)
    round_dir = args.round_dir.resolve()
    labels = {k: Path(v) for k, v in (a.split("=", 1) for a in args.labels)}
    cards = build_cards(
        round_dir,
        calibration_record=args.calibration,
        labels=labels,
        n_perm=args.n_perm,
        seed=args.seed,
    )
    write_cards(cards, args.out or round_dir / "cards")
    counts: dict[str, int] = defaultdict(int)
    for c in cards:
        counts[str(c["card_status"])] += 1
    print(f"{len(cards)} cards: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

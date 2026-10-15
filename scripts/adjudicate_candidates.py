#!/usr/bin/env python3
"""ADJ-01: deterministic claim level L0-L3 per (library, virus), thresholds from a protocol lock.

    python scripts/adjudicate_candidates.py <round_dir> [--lock protocol-v1.lock.json] \\
        [--metadata samples.tsv] [--germline-flags flags.tsv] [--out <round_dir>/adjudication/candidate_levels.tsv]

A candidate is what `scripts/dsr02_enumerate_calls.py` lists on `combined_off` (the arm `dsr03_evidence_cards.py`
uses; `detected = false` rows are skipped there). `library_id` is `<dataset>__<sample>`, the key of
`reference_cells/`. The package ships NO pass/fail number: every threshold comes from `--lock`. Without a lock,
or with a lock missing a key, every row is L0 (`no_protocol_lock` / `incomplete_protocol_lock:<key>`) and the
metric columns are still filled. L4 and L5 (longitudinal, orthogonal) are never emitted.

Ladder (first failing gate names `level_reason`; the DSR-02 class is shown as `read_verdict` and decides nothing):

  L0  candidate. Needs nothing.
  L1  unique_molecules > 0 (unique layer only: allocated mass never promotes); read-level evidence over the frozen
      reference cells (`molecule_scope` restricted) with virus_best fraction >= min_virus_best_fraction and
      complex-body fraction >= min_complex_body_fraction; no sibling bleed
      (allocated / (unique + allocated) <= max_sibling_allocated_fraction when the virus has a sibling group).
  L2  L1 plus unique molecules, positive cells and loci at or above the distributed-route minima. Loci and
      positive cells are LOCK-01 D1/D2 (independent overlap components; union of called barcodes), read from
      `multimap_evidence.tsv`; a run that predates those columns stops at L1 (`missing:...`).
  L3  L2 plus a gene-programme model for the virus (the nine herpesviruses, matched by taxid map) and >= l3_min_programme_cells cells
      in a programme state. Any other virus stops at L2 (`cap:no_gene_programme`).
  Anelloviridae (`claim_scope` screening_only) stop at L1 (`cap:anellovirus_screening_only`).

`package_tier_analogue` (last column) is the package evidence tier (`viralscan.evidence_tiers`) this level is the
analogue of: L0 candidate, L1 probable, L2 strong, L3 strong. The package tiers are calibrated molecule-count rules;
the study levels add read-level evidence on top, so it is an analogue, not an equivalence and not the package-assigned
tier of the row (nothing here runs the package tier rule).

Host expression is never read. Repeat classes of `read_verdict_summary.tsv` (tso, telomeric, low_complexity) are
flags only. `germline_integration_compatible` is passed through from `--germline-flags` or a `viral_summary.tsv`
column when one exists; nothing here detects it.

Lock keys (top level, or under an `adjudication` object): `distributed_route.min_unique_molecules`,
`distributed_route.min_positive_cells`, `distributed_route.min_independent_loci`, `max_sibling_allocated_fraction`,
`min_virus_best_fraction`, `min_complex_body_fraction`, `l3_min_programme_cells` (all numbers). Optional
`test_only: true` adds the flag `test_only_thresholds`.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "src") not in sys.path:  # run from a checkout without `pip install -e .`
    sys.path.insert(0, str(REPO / "src"))

from viralscan.evidence_tiers import package_tier_analogue  # noqa: E402
from viralscan.virus_grouping import (  # noqa: E402
    CLAIM_SCOPE_SCREENING,
    legacy_claim_scope,
    load_run_identity,
    virus_facts,
)

#: Column contract. LONG-01 reads this: append new columns at the end, never reorder or rename.
COLS = [
    "library_id", "virus_key", "level", "level_reason", "flags", "unique_molecules", "n_loci",
    "read_virus_best_fraction", "thresholds_sha256",
    "subject_id", "timepoint", "run_id", "lane",
    "virus_name", "allocated_molecules", "positive_cells", "max_sibling_allocated_fraction",
    "complex_body_fraction", "read_verdict", "molecule_scope", "gene_programme_cells",
    "package_tier_analogue",
]  # fmt: skip
METADATA_COLS = ("subject_id", "timepoint", "run_id", "lane")
REQUIRED_KEYS = (
    "distributed_route.min_unique_molecules",
    "distributed_route.min_positive_cells",
    "distributed_route.min_independent_loci",
    "max_sibling_allocated_fraction",
    "min_virus_best_fraction",
    "min_complex_body_fraction",
    "l3_min_programme_cells",
)
ARTEFACT_CLASSES = ("tso", "telomeric", "low_complexity")
#: virus_key -> the `virus` name in `gene_programs.tsv` (and so in `gene_program_summary.tsv`), for the nine
#: gene-programme herpesviruses (Q58). An identity-table run names these differently ("Human gammaherpesvirus 4"),
#: which used to read as "no model". Taxids: `virus_catalog.tsv`; names: `gene_programs.tsv`. Off the map = no model.
PROGRAMME_NAME_BY_KEY = {
    "taxid:10376": "Epstein-Barr virus",
    "taxid:10359": "Human cytomegalovirus",
    "taxid:10298": "Human herpesvirus 1",
    "taxid:10310": "Human herpesvirus 2",
    "taxid:32603": "Human herpesvirus 6A",
    "taxid:32604": "Human herpesvirus 6b",
    "taxid:10372": "Human herpesvirus 7",
    "taxid:37296": "Human herpesvirus 8",
    "taxid:10335": "Varicella-zoster virus",
}
PROGRAMME_STATE_COLS = ("n_cells_productive", "n_cells_latent", "n_cells_mixed", "n_cells_reactivating")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cards = _load("dsr03_evidence_cards")  # evidence_dir, read_level, frozen cells, enumerate_calls
verdicts = cards.verdicts


def _die(msg: str) -> None:
    print(f"adjudicate_candidates: {msg}", file=sys.stderr)
    sys.exit(2)


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------ lock


def load_lock(path: Path | None) -> dict:
    """``{"sha256", "values", "missing", "test_only"}``; ``sha256`` is ``""`` without a lock file."""
    if path is None:
        return {"sha256": "", "values": {}, "missing": None, "test_only": False}
    raw = path.read_bytes()
    try:
        doc = json.loads(raw)
    except ValueError as exc:
        _die(f"{path} is not valid JSON: {exc}")
    if not isinstance(doc, dict):
        _die(f"{path} must hold a JSON object")
    doc = doc.get("adjudication", doc)
    values: dict[str, float] = {}
    missing = None
    for key in REQUIRED_KEYS:
        node = doc
        for part in key.split("."):
            node = node.get(part) if isinstance(node, dict) else None
        if isinstance(node, bool) or not isinstance(node, (int, float)) or node < 0:
            missing = missing or key  # first incomplete key, in REQUIRED_KEYS order
        else:
            values[key] = node
    return {
        "sha256": hashlib.sha256(raw).hexdigest(),
        "values": values,
        "missing": missing,
        "test_only": doc.get("test_only") is True,
    }


# ------------------------------------------------------------------ the ladder


def adjudicate(m: dict, lock: dict) -> tuple[int, str, set[str]]:
    """``(level, level_reason, flags)`` from one candidate's metrics ``m``; pure, no I/O.

    Metrics: unique_molecules, allocated_molecules, positive_cells, n_loci (None = not measured),
    read_evidence (bool), molecule_scope, read_virus_best_fraction, complex_body_fraction,
    sibling_fraction (None = no sibling group), claim_scope, programme_model (bool), programme_cells,
    artefact_flags, germline (bool).
    """
    flags = set(m["artefact_flags"])
    if m["claim_scope"] == CLAIM_SCOPE_SCREENING:
        flags.add("screening_only")
    if m["germline"]:
        flags.add("germline_integration_compatible")
    if lock["test_only"]:
        flags.add("test_only_thresholds")
    if not lock["sha256"]:
        return 0, "no_protocol_lock", flags
    if lock["missing"]:
        return 0, f"incomplete_protocol_lock:{lock['missing']}", flags
    t = lock["values"]

    uniq = m["unique_molecules"]
    if uniq is None:
        return 0, "no_unique_layer", flags
    if uniq <= 0:  # unique-layer primacy: allocated mass alone never leaves L0
        return 0, ("allocated_only" if (m["allocated_molecules"] or 0) > 0 else "no_unique_molecules"), flags
    if not m["read_evidence"]:
        return 0, "no_read_level_evidence", flags
    if m["molecule_scope"] != "restricted":
        return 0, "read_level_not_restricted", flags
    vb, cb, sib = m["read_virus_best_fraction"], m["complex_body_fraction"], m["sibling_fraction"]
    if vb is None:
        return 0, "no_aligned_molecules", flags
    if vb < t["min_virus_best_fraction"]:
        return 0, "below_l1:min_virus_best_fraction", flags
    if cb is None:
        return 0, "missing:complex_body_fraction", flags
    if cb < t["min_complex_body_fraction"]:
        return 0, "below_l1:min_complex_body_fraction", flags
    if sib is not None and sib > t["max_sibling_allocated_fraction"]:
        flags.add("sibling_bleed")
        return 0, "below_l1:max_sibling_allocated_fraction", flags

    # L2: the distributed route (the rare-super-expressor route needs a per-cell unique count nothing writes yet)
    for key, value in (
        ("min_unique_molecules", uniq),
        ("min_positive_cells", m["positive_cells"]),
        ("min_independent_loci", m["n_loci"]),
    ):
        if value is None:
            return 1, f"missing:{key}", flags
        if value < t[f"distributed_route.{key}"]:
            return _cap(1, f"below_l2:{key}", m, flags)

    # L3
    if not m["programme_model"]:
        return _cap(2, "cap:no_gene_programme", m, flags)
    if m["programme_cells"] is None:
        return _cap(2, "missing:gene_programme_cells", m, flags)
    if m["programme_cells"] < t["l3_min_programme_cells"]:
        return _cap(2, "below_l3:l3_min_programme_cells", m, flags)
    return _cap(3, "meets_l3", m, flags)


def _cap(level: int, reason: str, m: dict, flags: set[str]) -> tuple[int, str, set[str]]:
    if m["claim_scope"] == CLAIM_SCOPE_SCREENING and level > 1:
        return 1, "cap:anellovirus_screening_only", flags
    return level, reason, flags


# ------------------------------------------------------------------ metrics from a round


def unique_layer(inner: Path, virus: str) -> dict | None:
    """Unique-layer metrics for ``virus`` from ``results/multimap_evidence.tsv``; None when absent.

    ``n_loci`` = distinct ``locus_id`` over genes with a unique molecule (D1); ``positive_cells`` = the union of
    called barcodes with a unique molecule at any gene (D2)."""
    path = inner / "results" / "multimap_evidence.tsv"
    if not path.is_file():
        return None
    rows = [r for r in _tsv(path) if r["virus_name"] == virus]
    if not rows:
        return None
    # LOCK-01 D1/D2. Both come from the run (`viral_loci`, `summarize_multimap_evidence`); a run that predates
    # them has no such column and the value stays None (`missing:` at L1), never a stand-in.
    unique_rows = [r for r in rows if float(r["viral_molecules_unique"] or 0) > 0]
    placed = "locus_id" in rows[0] and all(r["locus_id"] for r in unique_rows)
    positive = rows[0].get("n_unique_positive_cells_called", "")
    return {
        "unique_molecules": sum(float(r["viral_molecules_unique"] or 0) for r in rows),
        "allocated_molecules": sum(float(r["viral_molecules_ambiguous_allocated"] or 0) for r in rows),
        "n_loci": len({r["locus_id"] for r in unique_rows}) if placed else None,
        "positive_cells": int(float(positive)) if positive != "" else None,
    }


def programme(inner: Path, virus: str, virus_key: str = "") -> tuple[bool, int | None]:
    """``(has a gene-programme model, cells in a programme state)`` from ``results/gene_program_summary.tsv``.

    Rows match on the virus name, else (the nine herpesviruses only) on the catalogue name of ``virus_key``."""
    path = inner / "results" / "gene_program_summary.tsv"
    if not path.is_file():
        return False, None
    names = {virus, PROGRAMME_NAME_BY_KEY.get(virus_key, virus)}
    for r in _tsv(path):
        if r["virus_name"] in names:
            if r.get("panel_completeness", "") in ("", "not_applicable"):
                return False, None
            return True, sum(int(r[c]) for c in PROGRAMME_STATE_COLS if (r.get(c) or "").isdigit())
    return False, None


def sibling_fraction(unique: float | None, allocated: float | None, sibling_group: str) -> float | None:
    """allocated / (unique + allocated) when the virus has a sibling group, else None (not applicable)."""
    if not sibling_group or unique is None:
        return None
    total = unique + (allocated or 0.0)
    return (allocated or 0.0) / total if total else 0.0


def complex_fraction(ev: Path | None) -> float | None:
    if ev is None or not (ev / "alignment_qc.tsv").is_file():
        return None
    virus = [
        r
        for r in _tsv(ev / "alignment_qc.tsv")
        if r["count_layer"] == "deduplicated" and r["reference_class"] == "virus"
    ]
    return verdicts.read_weighted_complex_fraction(virus)


def read_virus_best_fraction(tiers: dict) -> float | None:
    counts = [tiers.get(k) for k in ("virus_best_molecules", "host_best_molecules", "tie_molecules")]
    if any(c in (None, "") for c in counts):
        return None
    vb, hb, tie = (int(c) for c in counts)
    return vb / (vb + hb + tie) if vb + hb + tie else None  # ties count against the virus, as in DSR-02


def load_flags(path: Path | None) -> dict[tuple[str, str], bool]:
    if path is None:
        return {}
    return {
        (r["library_id"], r["virus_key"]): r["germline_integration_compatible"].strip().lower() == "true"
        for r in _tsv(path)
    }


def load_metadata(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    out: dict[str, dict[str, str]] = {}
    for r in _tsv(path):
        if "library_id" not in r:
            _die(f"{path} has no library_id column")
        if r["library_id"] in out:
            _die(f"{path}: duplicate library_id {r['library_id']}")
        out[r["library_id"]] = {c: r.get(c, "") for c in METADATA_COLS}
    return out


def adjudicate_round(
    round_dir: Path,
    lock: dict,
    *,
    metadata: dict[str, dict[str, str]] | None = None,
    germline: dict[tuple[str, str], bool] | None = None,
) -> list[dict[str, object]]:
    rows = []
    for call in cards.enumerate_calls.calls(round_dir):
        if call["arm"] != cards.REFERENCE_ARM:
            continue
        inner = Path(call["run_dir"])
        dataset, sample, virus = call["dataset"], inner.parent.name, call["virus"]
        library_id = f"{dataset}__{sample}"
        identity = load_run_identity(inner)
        facts = virus_facts(identity).get(virus) if identity is not None else None
        virus_key = facts.virus_key if facts else f"name:{virus}"
        sibling_group = facts.sibling_group if facts else ""
        claim_scope = facts.claim_scope if facts else legacy_claim_scope(virus)
        ev = cards.evidence_dir(round_dir, dataset, sample, virus)
        tiers = cards.read_level(ev, cards._frozen(round_dir, dataset, sample))
        layer = unique_layer(inner, virus) or {}
        has_model, programme_cells = programme(inner, virus, virus_key)
        summary = {r["virus_name"]: r for r in _tsv(inner / "results" / "viral_summary.tsv")}.get(virus, {})
        g = (library_id, virus_key)
        flag_col = (summary.get("germline_integration_compatible") or "").strip().lower() == "true"
        m = {
            "unique_molecules": layer.get("unique_molecules"),
            "allocated_molecules": layer.get("allocated_molecules"),
            "positive_cells": layer.get("positive_cells"),
            "n_loci": layer.get("n_loci"),
            "read_evidence": ev is not None,
            "molecule_scope": tiers.get("molecule_scope", ""),
            "read_virus_best_fraction": read_virus_best_fraction(tiers),
            "complex_body_fraction": complex_fraction(ev),
            "sibling_fraction": sibling_fraction(
                layer.get("unique_molecules"), layer.get("allocated_molecules"), sibling_group
            ),
            "claim_scope": claim_scope,
            "programme_model": has_model,
            "programme_cells": programme_cells,
            "artefact_flags": {
                c
                for c in ARTEFACT_CLASSES
                if int(tiers.get(f"molecules_{c}", tiers.get(f"reads_{c}")) or 0) > 0
                # molecules (UMI representatives); legacy evidence without them falls back to the read label
            },
            "germline": flag_col or (germline or {}).get(g, False),
        }
        level, reason, flags = adjudicate(m, lock)
        rows.append(
            {
                "library_id": library_id,
                "virus_key": virus_key,
                "level": f"L{level}",
                "level_reason": reason,
                "flags": ";".join(sorted(flags)),
                "unique_molecules": _num(m["unique_molecules"]),
                "n_loci": _num(m["n_loci"]),
                "read_virus_best_fraction": _num(m["read_virus_best_fraction"]),
                "thresholds_sha256": lock["sha256"],
                **(metadata or {}).get(library_id, {}),
                "virus_name": virus,
                "allocated_molecules": _num(m["allocated_molecules"]),
                "positive_cells": _num(m["positive_cells"]),
                "max_sibling_allocated_fraction": _num(m["sibling_fraction"]),
                "complex_body_fraction": _num(m["complex_body_fraction"]),
                "read_verdict": tiers.get("read_verdict", tiers.get("read_level", "")),
                "molecule_scope": m["molecule_scope"],
                "gene_programme_cells": _num(programme_cells),
                "package_tier_analogue": package_tier_analogue(f"L{level}"),
            }
        )
    return sorted(rows, key=lambda r: (str(r["library_id"]), str(r["virus_key"])))


def _num(x: object) -> object:
    return "" if x is None else (round(x, 6) if isinstance(x, float) else x)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round_dir", type=Path)
    ap.add_argument("--lock", type=Path, help="protocol lock JSON holding the thresholds (none: all L0)")
    ap.add_argument("--metadata", type=Path, help="TSV keyed by library_id: subject_id, timepoint, run_id, lane")
    ap.add_argument("--germline-flags", type=Path, help="TSV: library_id, virus_key, germline_integration_compatible")
    ap.add_argument("--out", type=Path, help="default <round_dir>/adjudication/candidate_levels.tsv")
    args = ap.parse_args(argv)
    round_dir = args.round_dir.resolve()
    rows = adjudicate_round(
        round_dir,
        load_lock(args.lock),
        metadata=load_metadata(args.metadata),
        germline=load_flags(args.germline_flags),
    )
    out = args.out or round_dir / "adjudication" / "candidate_levels.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)
    tally: dict[str, int] = {}
    for r in rows:
        tally[str(r["level"])] = tally.get(str(r["level"]), 0) + 1
    print(f"{len(rows)} candidates: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

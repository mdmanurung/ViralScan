#!/usr/bin/env python3
"""Scorer for the v3 truth panel: precision/recall/F1, AUPRC, sibling and host-virus error (PLAN VAL-06).

Pure functions over sets and mappings, so each metric is checked against hand-computed fixtures
(tests/test_v3_scorer.py, VAL-07). Denominators follow ``analysis/v3_validation/protocol.yaml``
(metrics M2-M6, D6): recall divides by planted truth, precision by called items, and an undefined
ratio is ``None``, never 0, so a run with nothing called cannot read as perfectly wrong or right.

No function launches a workflow or estimates a fitted LOD. Pure metrics accept declared truth,
and the native adapter reads pinned completed output with its count-layer contract. Outcome
executors must enforce the protocol gate before calling this scorer (``score_production`` is that
gated entry). The truth-manifest columns
come from ``protocol.yaml`` so the scorer and generator share one schema.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections.abc import Hashable, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

try:
    from scripts.validate_v3_protocol import require_execution_allowed
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution or ``scripts/`` on sys.path.
    from validate_v3_protocol import require_execution_allowed  # type: ignore[no-redef]

PROTOCOL =Path(__file__).resolve().parents[1] / "analysis" / "v3_validation" / "protocol.yaml"
TRUE, FALSE = {"true", "1", "yes"}, {"false", "0", "no", ""}

#: Truth-tier vocabulary of docs/plans/real-data-sensitivity-specificity.md section 2. Tiers are
#: reported separately and never pooled; which dataset belongs to which tier is the user's decision
#: (G3 round), so the scorer takes the tier as input and only validates it. No cut-off lives here.
TIERS = {
    "T1": "exact: planted molecules, recorded per read or molecule",
    "T2": "design-known: infected vs mock, a cell line that carries the virus",
    "T3": "orthogonal assay on the same sample",
    "T4": "presumed negative",
    "T5": "contextual: published results, not scored",
}
#: truth_manifest origins that are contamination, not infection (generator ``hopped``/``ambient``).
CONTAMINATION_ORIGINS = ("hopped", "ambient")
HOP_COLUMNS = ("hop_donor_sample_id", "hop_donor_molecule_id", "hop_cb_preserved", "hop_landing")


def ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def precision_recall_f1(truth: Iterable[Hashable], called: Iterable[Hashable]) -> dict[str, Any]:
    """M3/M4: recall over truth items, precision over called items; ``None`` when undefined."""
    truth, called = set(truth), set(called)
    tp = len(truth & called)
    precision, recall = ratio(tp, len(called)), ratio(tp, len(truth))
    f1 = ratio(2 * tp, len(called) + len(truth))  # = 2PR/(P+R), defined whenever anything exists
    return dict(
        tp=tp,
        fp=len(called - truth),
        fn=len(truth - called),
        precision=precision,
        recall=recall,
        f1=f1,
    )


def auprc(labels: Sequence[int | bool], scores: Sequence[float]) -> float | None:
    """D6: average precision, sum over distinct score thresholds of (R_k - R_{k-1}) * P_k.

    Items with equal scores enter together, so the result does not depend on input order. ``None``
    when there is no positive, where recall is undefined.
    """
    if len(labels) != len(scores):
        raise ValueError("labels and scores differ in length")
    if any(label not in (0, 1) for label in labels):
        raise ValueError("labels must be binary")
    if any(not math.isfinite(value) for value in scores):
        raise ValueError("scores must be finite")
    positives = sum(1 for label in labels if label)
    if not positives:
        return None
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    tp = seen = 0
    previous_recall = area = 0.0
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and scores[order[j]] == scores[order[i]]:
            tp += 1 if labels[order[j]] else 0
            seen += 1
            j += 1
        recall = tp / positives
        area += (recall - previous_recall) * (tp / seen)
        previous_recall, i = recall, j
    return area


def negative_false_call_rate(called: Mapping[str, bool], complete: Iterable[str]) -> dict[str, Any]:
    """M2: fraction of complete negative samples with at least one probable-or-strong call.

    ``complete`` names the rows that finished and are comparable. A row that failed or timed out is
    not estimable here and must be left out of it (and counted as technical non-completion
    elsewhere), not scored as a correct negative.
    """
    rows = list(complete)
    if len(set(rows)) != len(rows):
        raise ValueError("complete sample denominator contains duplicates")
    if any(not isinstance(called[sample], bool) for sample in rows if sample in called):
        raise ValueError("call statuses must be booleans")
    unknown = [sample for sample in rows if sample not in called]
    if unknown:
        raise ValueError(f"complete samples without a call status: {sorted(unknown)[:3]}")
    false_calls = sum(1 for sample in rows if called[sample])
    return dict(n=len(rows), false_calls=false_calls, rate=ratio(false_calls, len(rows)))


def sibling_confusion(
    planted: Mapping[Hashable, str],
    assigned: Mapping[Hashable, str | None],
    target: str,
    sibling: str,
) -> dict[str, Any]:
    """M5: planted molecules of ``target`` that were assigned to ``sibling``, over all of them.

    ``planted`` maps molecule -> planted virus; ``assigned`` maps molecule -> the virus it was
    assigned to (``None``/absent: not assigned).
    """
    mine = [molecule for molecule, virus in planted.items() if virus == target]
    confused = sum(1 for molecule in mine if assigned.get(molecule) == sibling)
    return dict(n_target=len(mine), n_to_sibling=confused, rate=ratio(confused, len(mine)))


def host_virus_allocation_error(
    truth_host: Iterable[Hashable],
    truth_viral: Iterable[Hashable],
    assigned_viral: Iterable[Hashable],
    assigned_host: Iterable[Hashable],
) -> dict[str, Any]:
    """M6, both directions: host molecules given a viral feature, and viral given a host feature."""
    truth_host, truth_viral = set(truth_host), set(truth_viral)
    host_to_viral = len(truth_host & set(assigned_viral))
    viral_to_host = len(truth_viral & set(assigned_host))
    return dict(
        host_to_viral=ratio(host_to_viral, len(truth_host)),
        viral_to_host=ratio(viral_to_host, len(truth_viral)),
        n_host=len(truth_host),
        n_viral=len(truth_viral),
    )


def negative_cell_false_call_rate(
    called: Iterable[Hashable], negative_anchor: Iterable[Hashable]
) -> dict[str, Any]:
    """D7: false-positive cells over an exact host-only/host-homology negative anchor.

    Call once per tier and negative class. Outside-anchor calls are integrity failures,
    not a reason to widen the prespecified denominator.
    """
    anchor, called = set(negative_anchor), set(called)
    if called - anchor:
        raise ValueError("called cells outside the declared negative anchor")
    return dict(n=len(anchor), false_calls=len(called), rate=ratio(len(called), len(anchor)))


def tier_specificity_ppv(
    truth_positive: Iterable[Hashable],
    called_by_tier: Mapping[str, Iterable[Hashable]],
    anchor: Iterable[Hashable],
    *,
    twin_excluded: Iterable[Hashable] = (),
) -> dict[str, dict[str, Any]]:
    """E9/D18/D19: specificity and PPV, sharing the twin exclusion across all tiers.

    These are tier-specific classification proportions, not nominal-confidence calibration.
    Uncertainty and technical noncompletion must be reported separately by the caller.
    """
    anchor, truth, excluded = set(anchor), set(truth_positive), set(twin_excluded)
    if (truth | excluded) - anchor:
        raise ValueError("truth or twin exclusions outside the declared anchor")
    result = {}
    for tier, calls in called_by_tier.items():
        if tier not in {"candidate", "probable", "strong"}:
            raise ValueError(f"unknown evidence tier {tier!r}")
        calls = set(calls)
        if calls - anchor:
            raise ValueError("called cells outside the declared anchor")
        eligible = anchor - excluded
        truth_eligible, calls = truth & eligible, calls - excluded
        negatives = eligible - truth_eligible
        tp, fp = len(truth_eligible & calls), len(negatives & calls)
        tn = len(negatives) - fp
        result[tier] = dict(
            tp=tp,
            fp=fp,
            tn=tn,
            n_negative=len(negatives),
            n_called=len(calls),
            n_twin_excluded=len(excluded),
            specificity=ratio(tn, len(negatives)),
            ppv=ratio(tp, len(calls)),
        )
    return result


def cross_tool_relative_difference(viralscan: float, comparator: float) -> float | None:
    """M9: signed difference divided by mean unique counts; both zero is not applicable."""
    if any(not math.isfinite(x) or x < 0 or x != int(x) for x in (viralscan, comparator)):
        raise ValueError("unique-only counts must be finite nonnegative integers")
    return ratio(viralscan - comparator, (viralscan + comparator) / 2)


def read_run_count_layer(
    sample_dir: str | Path,
    *,
    barcodes: Sequence[str],
    features: Sequence[str],
    layer: str,
    denominator: float,
    expected_run_fingerprint: str,
    expected_counts_sha256: str,
) -> tuple[Any, dict[str, Any]]:
    """Read audited native unique counts for VAL-08, restricted once to declared cells/genes.

    The caller supplies outcome-independent, exact native barcode and feature identifiers,
    the declared unique layer/denominator, and pinned run/artifact identities. No barcode
    normalization, alias rescue, inferred structural zeros, or workflow execution occurs.
    This matrix adapter cannot establish exact molecule recovery without CB-UMI lineage.
    """
    import anndata as ad
    import numpy as np
    from scipy import sparse

    from viralscan.kb_outputs import KbCountOutputs
    from viralscan.run_safety import sha256_file
    from viralscan.validation import (
        _completion_issues,
        _matrix_issues,
        h5ad_contract,
        packaged_schema_resource,
        validate_json_schema,
    )

    if layer != "counts_unique":
        raise ValueError("unique-only scoring requires the declared counts_unique layer")
    if any(
        not isinstance(value, str) or not re.fullmatch(r"[a-fA-F0-9]{64}", value)
        for value in (expected_run_fingerprint, expected_counts_sha256)
    ):
        raise ValueError(
            "run fingerprint and source artifact digest must be declared SHA-256 values"
        )
    if len(barcodes) == 0 or len(features) == 0:
        raise ValueError("declared barcode and feature universes must be nonempty")
    if len(set(barcodes)) != len(barcodes) or len(set(features)) != len(features):
        raise ValueError("declared barcode or feature universe contains duplicates")
    sample_dir = Path(sample_dir)
    manifest_path = sample_dir.parent / "run_manifest.json"
    if not manifest_path.is_file():
        manifest_path = sample_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    issues = validate_json_schema(manifest, packaged_schema_resource("run_manifest.schema.json"))
    if issues:
        raise ValueError("invalid run manifest: " + "; ".join(i.message for i in issues))
    if manifest.get("run_fingerprint") != expected_run_fingerprint:
        raise ValueError("run fingerprint mismatch")
    path = KbCountOutputs(sample_dir).adata_multimap
    digest = sha256_file(path)
    if digest != expected_counts_sha256:
        raise ValueError("count artifact digest mismatch")
    issues = _completion_issues(manifest_path.parent.resolve(), manifest)
    if issues:
        raise ValueError("invalid completion receipt: " + "; ".join(i.code for i in issues))
    marker_path = manifest_path.parent / "run_complete.json"
    marker = json.loads(marker_path.read_text())
    relative_artifact = str(path.resolve().relative_to(manifest_path.parent.resolve()))
    if marker["artifacts"].get(relative_artifact) != digest:
        raise ValueError("completion receipt does not pin the selected count artifact")
    adata = ad.read_h5ad(path)
    issues = _matrix_issues(
        adata, path, h5ad_contract(), packaged_schema_resource("count_audit.schema.json")
    )
    if issues:
        raise ValueError("invalid count-layer contract: " + "; ".join(i.code for i in issues))
    if adata.uns["count_schema_version"] != "3.0.0":
        raise ValueError("count schema version mismatch")
    if adata.uns["quantification_unit"] != manifest["quantification_unit"]:
        raise ValueError("quantification unit mismatch")
    if adata.uns["multimap_method"] != manifest["allocation_method"]:
        raise ValueError("allocation method mismatch")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("output barcode or feature universe contains duplicates")
    if set(barcodes) - set(adata.obs_names):
        raise ValueError("declared barcode universe is absent from the raw output")
    if set(features) - set(adata.var_names):
        raise ValueError("declared feature universe is absent from the raw output")
    unique = adata.layers[layer]
    values = unique.data if sparse.issparse(unique) else np.asarray(unique)
    if not np.equal(values, np.floor(values)).all():
        raise ValueError("counts_unique contains fractional counts")
    matrix = adata[list(barcodes), list(features)].layers[layer].copy()
    measured = float(matrix.sum())
    if not math.isfinite(denominator) or denominator < 0 or denominator != measured:
        raise ValueError("denominator mismatch for the declared cells/features/count layer")
    provenance = dict(
        run_fingerprint=expected_run_fingerprint,
        run_manifest=str(manifest_path.resolve()),
        run_manifest_sha256=sha256_file(manifest_path),
        run_complete_sha256=sha256_file(marker_path),
        source_artifact=str(path.resolve()),
        source_artifact_sha256=digest,
        source_layer=layer,
        unit="integer-unique-molecule",
        integer_status=True,
        umi_collapse_semantics=manifest["quantification_unit"],
        ambiguity_rule="gene-unique-only",
        allocation_method=manifest["allocation_method"],
        barcodes=list(barcodes),
        features=list(features),
        denominator=measured,
        excluded_output_barcodes=sorted(set(adata.obs_names) - set(barcodes)),
        excluded_output_features=sorted(set(adata.var_names) - set(features)),
    )
    return matrix, provenance


# ── truth tiers and negative libraries (PLAN M04-SCORE) ─────────────────────


def check_tier(tier: str) -> str:
    """Return ``tier`` if it is one of T1-T5, else raise. T5 is context and is never scored."""
    if tier not in TIERS:
        raise ValueError(f"unknown truth tier {tier!r}; expected one of {sorted(TIERS)}")
    if tier == "T5":
        raise ValueError("T5 is contextual (published results) and is not scored")
    return tier


def realised_shortfalls(detail: Mapping[str, Any]) -> list[tuple[str, int, int]]:
    """Constructions a sample asked for but did not realise, as (name, requested, realised).

    ``detail`` is one entry of ``run_manifest.json`` ``samples``. A real-background (D1) sample is
    exempt: it carries ``background_library_id`` and has ``realised.hopped == 0`` by design (one
    library, nobody to hop from) beside ``requested_hopped > 0``.
    """
    if "background_library_id" in detail:
        return []
    done = detail.get("realised", {})
    pairs = {
        "collision": (
            done.get("requested_collision", 0),
            done.get("collision_same", 0) + done.get("collision_disjoint", 0),
        ),
        "ambient": (done.get("requested_ambient", 0), done.get("ambient", 0)),
        "hopped": (done.get("requested_hopped", 0), done.get("hopped", 0)),
    }
    return [(name, want, got) for name, (want, got) in pairs.items() if want != got]


def truth_positive_molecules(rows: Iterable[Mapping[str, str]]) -> set[tuple[str, str, str]]:
    """(cell_barcode, umi, virus) of planted molecules. Background (``unknown``) is never one."""
    return {
        (r["cell_barcode"], r["umi"], r["planted_virus_id"])
        for r in rows
        if r["is_planted_molecule"].strip().lower() in TRUE and r["origin"] != "unknown"
    }


def _contamination(
    rows: Iterable[Mapping[str, str]],
) -> dict[tuple[str, str, str], Mapping[str, str]]:
    """Hopped/ambient viral molecules of one library, keyed (cb, umi, virus), with provenance."""
    found: dict[tuple[str, str, str], Mapping[str, str]] = {}
    for r in rows:
        if r["origin"] in CONTAMINATION_ORIGINS and r["planted_virus_id"]:
            found.setdefault((r["cell_barcode"], r["umi"], r["planted_virus_id"]), r)
    return found


def _evidence(rows: Iterable[Mapping[str, str]]) -> dict[str, Any]:
    rows = list(rows)
    landing: dict[str, int] = {}
    for r in rows:
        if r.get("hop_landing"):
            landing[r["hop_landing"]] = landing.get(r["hop_landing"], 0) + 1
    return dict(
        n_molecules=len(rows),
        n_hopped=sum(1 for r in rows if r["origin"] == "hopped"),
        n_ambient=sum(1 for r in rows if r["origin"] == "ambient"),
        n_with_hop_provenance=sum(1 for r in rows if r.get("hop_donor_sample_id")),
        n_cb_preserved=sum(1 for r in rows if r.get("hop_cb_preserved", "").lower() in TRUE),
        donors=sorted({r["hop_donor_sample_id"] for r in rows if r.get("hop_donor_sample_id")}),
        landing=dict(sorted(landing.items())),
    )


def score_negative_libraries(
    libraries: Sequence[Mapping[str, Any]], viruses: Iterable[str] = ()
) -> dict[str, Any]:
    """Negative-library specificity, hopped/ambient molecules reported THREE ways, together.

    Each library: ``sample_id``, ``tier`` (T1-T4), ``complete`` (default True; an incomplete row is
    not estimable and is left out), ``called`` (virus -> (cb, umi) molecules supporting that call;
    a virus absent from the mapping is not called) and ``truth_rows`` (its rows from
    :func:`read_truth_manifest`). Libraries of different tiers are refused, never pooled.

    1. ``unconditional`` (the HEADLINE): any-virus sample-level false calls over complete libraries,
       contamination counted as false calls, with per-virus rates alongside.
    2. ``excluded_with_evidence``: libraries with documented contamination (a hopped or ambient
       viral molecule in the truth manifest) whose call, if any, is fully explained by it (every
       called virus has supporting molecules and all are hopped/ambient) leave the denominator,
       with their provenance listed. A call with any other supporting molecule stays a false call.
    3. ``contamination_attributable``: the false calls removed by (2), over the same denominator as
       (1), beside the unattributable residual; ``unconditional = attributable + residual``.

    No pass/fail cut-off is applied.
    """
    tiers = {check_tier(lib["tier"]) for lib in libraries}
    if len(tiers) > 1:
        raise ValueError(
            f"refusing to pool truth tiers {sorted(tiers)}; score each tier separately"
        )
    done = [lib for lib in libraries if lib.get("complete", True)]
    by_id = {lib["sample_id"]: lib for lib in done}
    found, excluded, explained = {}, {}, set()
    for sid, lib in by_id.items():
        if any(r["origin"] == "cell" and r["planted_virus_id"] for r in lib["truth_rows"]):
            raise ValueError(f"{sid}: a negative library holds cell-origin planted molecules")
        found[sid] = _contamination(lib["truth_rows"])
        called = {v: set(m) for v, m in lib["called"].items()}
        fully = all(
            m and all((cb, umi, v) in found[sid] for cb, umi in m) for v, m in called.items()
        )
        if found[sid] and fully:
            excluded[sid] = found[sid]
        if called and fully:
            explained.add(sid)
    flags = {sid: bool(lib["called"]) for sid, lib in by_id.items()}
    headline = negative_false_call_rate(flags, list(flags))
    rest = [sid for sid in flags if sid not in excluded]
    names = sorted({v for lib in done for v in lib["called"]} | set(viruses))

    def per_virus(ids: list[str]) -> dict[str, dict[str, Any]]:
        counts = {v: sum(1 for i in ids if v in by_id[i]["called"]) for v in names}
        return {v: dict(false_calls=k, rate=ratio(k, len(ids))) for v, k in counts.items()}

    residual = headline["false_calls"] - len(explained)
    return dict(
        tier=next(iter(tiers), None),
        headline="unconditional",
        n_incomplete=len(libraries) - len(done),
        unconditional=headline | dict(per_virus=per_virus(list(flags))),
        excluded_with_evidence=negative_false_call_rate(flags, rest)
        | dict(
            n_excluded=len(excluded),
            excluded=[
                dict(sample_id=s, evidence=_evidence(excluded[s].values()))
                for s in sorted(excluded)
            ],
            per_virus=per_virus(rest),
        ),
        contamination_attributable=dict(
            n=headline["n"],
            false_calls=len(explained),
            rate=ratio(len(explained), headline["n"]),
            residual_false_calls=residual,
            residual_rate=ratio(residual, headline["n"]),
        ),
    )


def score_negative_libraries_by_tier(
    libraries: Sequence[Mapping[str, Any]], viruses: Iterable[str] = ()
) -> dict[str, dict[str, Any]]:
    """One :func:`score_negative_libraries` report per tier; there is no pooled entry."""
    groups: dict[str, list[Mapping[str, Any]]] = {}
    for lib in libraries:
        groups.setdefault(check_tier(lib["tier"]), []).append(lib)
    return {t: score_negative_libraries(g, viruses) for t, g in sorted(groups.items())}


def score_production(
    phase: str, libraries: Sequence[Mapping[str, Any]], viruses: Iterable[str] = ()
) -> dict[str, dict[str, Any]]:
    """Gated scorer entry for a protocol-governed run (SCI-05 r12 M4).

    Fails closed unless the protocol permits outcome-generating runs at ``phase``; then one
    report per tier. The pure functions above stay ungated so exploratory runs and tests can use
    them; an executor scoring real outcomes must come through here.
    """
    require_execution_allowed(phase)
    return score_negative_libraries_by_tier(libraries, viruses)


def format_negative_report(report: Mapping[str, Any]) -> str:
    """Plain text of one tier's report: all three numbers, the unconditional one labelled HEADLINE."""

    def line(label: str, d: Mapping[str, Any]) -> str:
        rate = "n/a" if d["rate"] is None else f"{d['rate']:.4f}"
        return f"{label}: {d['false_calls']}/{d['n']} = {rate}"

    u, e, c = (
        report[k] for k in ("unconditional", "excluded_with_evidence", "contamination_attributable")
    )
    return "\n".join(
        [
            f"tier {report['tier']}; {report['n_incomplete']} incomplete librar(ies) left out",
            line("(i) HEADLINE unconditional any-virus false call", u),
            line(f"(ii) excluded with evidence ({e['n_excluded']} libraries excluded)", e),
            line("(iii) contamination-attributable", c),
            f"     residual (not attributable): {c['residual_false_calls']}/{c['n']}",
        ]
    )


# ── truth manifest ───────────────────────────────────────────────────────────


def required_columns(artifact: str, protocol: Path = PROTOCOL) -> list[str]:
    """Column contract of ``analysis/v3_validation/generated/<artifact>.tsv`` from the protocol."""
    import yaml

    document = yaml.safe_load(protocol.read_text())

    def walk(node: Any):
        if isinstance(node, dict):
            if (
                str(node.get("path", "")).endswith(f"/{artifact}.tsv")
                and "required_columns" in node
            ):
                yield node["required_columns"]
            for value in node.values():
                yield from walk(value)
        elif isinstance(node, list):
            for value in node:
                yield from walk(value)

    found = list(walk(document))
    if not found:
        raise KeyError(f"no required_columns for {artifact}.tsv in {protocol}")
    return list(found[0])


def _flag(value: str, where: str) -> bool:
    if not isinstance(value, str):
        raise ValueError(f"{where}: missing boolean value")
    text = value.strip().lower()
    if text in TRUE:
        return True
    if text in FALSE:
        return False
    raise ValueError(f"{where}: {value!r} is not a boolean")


def read_truth_manifest(path: str | Path, protocol: Path = PROTOCOL) -> list[dict[str, str]]:
    """Rows of ``truth_manifest.tsv``, rejecting a malformed one rather than scoring it.

    Checks the protocol's column contract, booleans, and the two rules that keep truth honest: a
    planted molecule names its virus and gene, and an ``artefact:*`` row is never planted truth.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:  # type: ignore[call-overload]
        reader = csv.DictReader(handle, delimiter="\t")
        if len(set(reader.fieldnames or [])) != len(reader.fieldnames or []):
            raise ValueError(f"{path}: duplicate truth manifest columns")
        missing = [
            c
            for c in required_columns("truth_manifest", protocol)
            if c not in (reader.fieldnames or [])
        ]
        if missing:
            raise ValueError(f"{path}: truth manifest lacks columns {missing}")
        rows = list(reader)
    for number, row in enumerate(rows, start=2):
        where = f"{path}:{number}"
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"{where}: malformed truth manifest row")
        planted = _flag(row["is_planted_molecule"], where)
        _flag(row["is_infected_cell"], where)
        if planted and not (row["planted_virus_id"] and row["planted_gene_id"]):
            raise ValueError(f"{where}: planted molecule without planted_virus_id/planted_gene_id")
        if row["origin"].startswith("artefact:") and (planted or row["planted_virus_id"]):
            raise ValueError(f"{where}: artefact row labelled as planted viral truth")
        if row["partition"] not in {"training", "holdout"}:
            raise ValueError(f"{where}: partition {row['partition']!r} is not training/holdout")
        if row["origin"] == "unknown" and (planted or row["planted_virus_id"]):
            raise ValueError(f"{where}: unknown-origin background labelled as planted viral truth")
        hop = {c: row[c] for c in HOP_COLUMNS if c in row}  # only in --hop-preserve-cb output
        if any(hop.values()) and row["origin"] != "hopped":
            raise ValueError(f"{where}: hop provenance on a row whose origin is not hopped")
        if hop.get("hop_landing", "") not in {"", "cell", "empty", "other"}:
            raise ValueError(f"{where}: hop_landing {hop['hop_landing']!r} is not cell/empty/other")
        if "hop_cb_preserved" in hop:
            _flag(hop["hop_cb_preserved"], where)
    return rows

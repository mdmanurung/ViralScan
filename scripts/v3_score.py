#!/usr/bin/env python3
"""Scorer for the v3 truth panel: precision/recall/F1, AUPRC, sibling and host-virus error (PLAN VAL-06).

Pure functions over sets and mappings, so each metric is checked against hand-computed fixtures
(tests/test_v3_scorer.py, VAL-07). Denominators follow ``analysis/v3_validation/protocol.yaml``
(metrics M2-M6, D6): recall divides by planted truth, precision by called items, and an undefined
ratio is ``None``, never 0, so a run with nothing called cannot read as perfectly wrong or right.

No function launches a workflow or estimates a fitted LOD. Pure metrics accept declared truth,
and the native adapter reads pinned completed output with its count-layer contract. Outcome
executors must enforce the protocol gate before calling this scorer. The truth-manifest columns
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

PROTOCOL = Path(__file__).resolve().parents[1] / "analysis" / "v3_validation" / "protocol.yaml"
TRUE, FALSE = {"true", "1", "yes"}, {"false", "0", "no", ""}


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
    return rows

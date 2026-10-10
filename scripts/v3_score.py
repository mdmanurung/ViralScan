#!/usr/bin/env python3
"""Scorer for the v3 truth panel: precision/recall/F1, AUPRC, sibling and host-virus error (PLAN VAL-06).

Pure functions over sets and mappings, so each metric is checked against hand-computed fixtures
(tests/test_v3_scorer.py, VAL-07). Denominators follow ``analysis/v3_validation/protocol.yaml``
(metrics M2-M6, D6): recall divides by planted truth, precision by called items, and an undefined
ratio is ``None``, never 0, so a run with nothing called cannot read as perfectly wrong or right.

Nothing here runs an outcome: it does not call ``require_execution_allowed`` and reads no protocol
result. The truth-manifest column contract is read from ``protocol.yaml`` so the scorer and the
generator share one schema.
"""

from __future__ import annotations

import csv
import gzip
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
    rows = [sample for sample in complete]
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
        planted = _flag(row["is_planted_molecule"], where)
        _flag(row["is_infected_cell"], where)
        if planted and not (row["planted_virus_id"] and row["planted_gene_id"]):
            raise ValueError(f"{where}: planted molecule without planted_virus_id/planted_gene_id")
        if row["origin"].startswith("artefact:") and (planted or row["planted_virus_id"]):
            raise ValueError(f"{where}: artefact row labelled as planted viral truth")
        if row["partition"] not in {"training", "holdout"}:
            raise ValueError(f"{where}: partition {row['partition']!r} is not training/holdout")
    return rows

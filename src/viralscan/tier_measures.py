"""The three tier-cutoff measures m, F and S behind the frozen ``probable`` grid (PLAN ``TIER-02``, Q135, Q141).

The contract is ``calibration.freeze_scope`` in ``analysis/v3_validation/protocol.yaml``. Each measure is
computed from the product's own outputs, with no truth label, over the run's **called cells**:

* ``m(c, v)`` (:func:`molecules_per_cell`): summed host-conservative ``adata.X`` over the index genes of ``v`` in
  cell ``c``, fractional weights kept, no denominator; a called cell with no viral molecule has ``m = 0``.
* ``F(s, v)`` (:func:`features_supported`): the number of index genes of ``v`` whose summed ``X`` over the sample's
  called cells is at least one molecule, reported beside ``n_index_genes`` (a virus with fewer index genes than
  ``f`` cannot reach ``probable`` at that grid point, :func:`cutoff_unattainable`).
* ``S(s, v)`` (:func:`support_fraction`): ``virus_best / (virus_best + host_best + tie)`` over the counted,
  UMI-deduplicated molecules of the sample's called cells that have a viral alignment, read from an **evidence
  directory** written by ``viralscan evidence`` (``molecule_verdicts.tsv.gz``). Missing, competitor-mode or
  zero-denominator evidence is **not estimable**, reason-coded, never a guess and never zero.

A cell is ``probable`` for ``v`` when ``m >= M and F >= f and S >= s*`` (:func:`probable_cells`); F and S are per
sample and virus and are broadcast to every cell of the sample.

The comparisons on molecule weight (``m >= M`` and the per-gene ``>= 1`` of F) carry the same relative slack as
:func:`viralscan.utils.infected_masks` (``INFECTED_RELATIVE_TOLERANCE``), so one molecule split
``1/3 + 1/3 + 1/3`` is not demoted by float32 rounding. Sums are float64.

S reads ``molecule_verdicts.tsv.gz`` as written by :mod:`viralscan.molecule_verdict`. Its ``(cb, ub)`` key is the
corrected CB-UMI (PLAN ``Q137`` audits that); this module still collapses any repeated key defensively, so a
molecule is one however many reads or rows carry it. The evidence run extracts only the reads that hit the target
virus's index genes, and the molecule verdicts carry no per-gene coordinates, so "aligned to an index gene of v"
is the evidence run's own restriction plus a present viral alignment score here.
"""

import csv
import gzip
import json
from collections.abc import Collection, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np

from viralscan.molecule_verdict import verdict_of
from viralscan.run_safety import sha256_file
from viralscan.utils import INFECTED_RELATIVE_TOLERANCE, matrix_for_genes, virus_molecules_per_cell

VERDICTS_FILE = "molecule_verdicts.tsv.gz"
MANIFEST_FILE = "evidence_manifest.json"

#: Reason codes of a not-estimable S.
REASON_NO_COUNTED = "no_counted_molecules"  # zero denominator (frozen in freeze_scope)
REASON_NO_EVIDENCE = "evidence_directory_missing"
REASON_NO_VERDICTS = "molecule_verdicts_missing"
REASON_COMPETITOR_MODE = "verdicts_not_produced_in_competitor_mode"


class EvidenceMismatch(ValueError):
    """The evidence directory does not belong to this virus or this matrix; S must not be guessed."""


class FeatureSupport(NamedTuple):
    n_supported: int  # F
    n_index_genes: int


@dataclass(frozen=True)
class SupportFraction:
    """S for one sample and virus; ``value`` is None (never 0) when not estimable, with a ``reason``."""

    value: float | None
    reason: str | None
    n_virus_best: int = 0
    n_host_best: int = 0
    n_tie: int = 0
    n_excluded_host_only: int = 0  # counted, called, aligned only to host: outside the population
    n_excluded_unaligned: int = (
        0  # counted, called, aligned to neither side: outside the population
    )


def _unique(genes: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(genes))


def _called(adata: Any, called_mask: Any) -> np.ndarray:
    mask = np.asarray(called_mask, dtype=bool)
    if mask.shape != (adata.n_obs,):
        raise ValueError(f"called_mask has shape {mask.shape}, expected ({adata.n_obs},)")
    return mask


def molecules_per_cell(adata: Any, index_genes: Iterable[str]) -> np.ndarray:
    """``m(c, v)`` over every barcode: summed ``adata.X`` on *index_genes*, float64, fractions kept."""
    present = [g for g in _unique(index_genes) if g in adata.var_names]
    return np.asarray(virus_molecules_per_cell(adata, adata.X, present), dtype=float)


def features_supported(adata: Any, index_genes: Iterable[str], called_mask: Any) -> FeatureSupport:
    """``F(s, v)``: index genes whose summed ``X`` over the called cells is >= 1 molecule (not merely > 0)."""
    genes = _unique(index_genes)
    present = [g for g in genes if g in adata.var_names]
    called = _called(adata, called_mask)
    if not present or not called.any():
        return FeatureSupport(0, len(genes))
    sub = matrix_for_genes(adata, adata.X, present)[called]
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    per_gene = np.asarray(sub, dtype=float).sum(axis=0)
    return FeatureSupport(int((per_gene >= 1.0 - INFECTED_RELATIVE_TOLERANCE).sum()), len(genes))


def cutoff_unattainable(feature_support: FeatureSupport, f: int) -> bool:
    """True when the virus has fewer index genes than *f*: reported as unattainable, never as zero."""
    return feature_support.n_index_genes < f


def _read_manifest(evidence_dir: Path, virus: str) -> dict[str, Any]:
    path = evidence_dir / MANIFEST_FILE
    if not path.is_file():
        return {}
    manifest = json.loads(path.read_text(encoding="utf-8"))
    label = (manifest.get("target") or {}).get("label")
    if label != virus:
        raise EvidenceMismatch(f"{path} is evidence for target {label!r}, not {virus!r}")
    return dict(manifest)


def _molecules(path: Path) -> dict[tuple[str, str], tuple[bool, int | None, int | None]]:
    """One entry per corrected ``(cb, ub)``: ``(counted, best viral AS, best host AS)``, repeated rows merged."""
    merged: dict[tuple[str, str], tuple[bool, int | None, int | None]] = {}

    def best(a: int | None, b: int | None) -> int | None:
        return b if a is None else a if b is None else max(a, b)

    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["cb"], row["ub"])
            counted, virus_as, host_as = merged.get(key, (False, None, None))
            merged[key] = (
                counted or row["counted"] == "True",
                best(virus_as, int(row["virus_as"]) if row["virus_as"] else None),
                best(host_as, int(row["host_as"]) if row["host_as"] else None),
            )
    return merged


def _not_estimable(reason: str) -> SupportFraction:
    return SupportFraction(None, reason)


def support_fraction(
    evidence_dir: str | Path | None,
    virus: str,
    called_barcodes: Collection[str],
    *,
    obs_names: Collection[str] | None = None,
) -> SupportFraction:
    """``S(s, v)`` from the evidence directory of *virus* for one sample.

    *called_barcodes* are the sample's called cells. The population is every counted molecule (>= 1 read with
    ``assigned_weight > 0``) of a called cell that has a viral alignment; a molecule aligned only to host is
    outside it. ``S = virus_best / (virus_best + host_best + tie)``, pooled over cells, a molecule counting once
    whatever its weight or read count.

    Raises :class:`EvidenceMismatch` when the manifest names another target, or when counted molecules exist and
    none of their barcodes is in *obs_names* (a barcode-format mismatch would otherwise read as no evidence).
    """
    if evidence_dir is None or not Path(evidence_dir).is_dir():
        return _not_estimable(REASON_NO_EVIDENCE)
    directory = Path(evidence_dir)
    manifest = _read_manifest(directory, virus)
    if (manifest.get("competition") or {}).get(
        "molecule_verdict"
    ) == "not_produced_in_competitor_mode":
        return _not_estimable(REASON_COMPETITOR_MODE)
    path = directory / VERDICTS_FILE
    if not path.is_file():
        return _not_estimable(REASON_NO_VERDICTS)
    expected = (manifest.get("outputs") or {}).get(VERDICTS_FILE)
    if expected is not None and sha256_file(path) != expected:
        raise EvidenceMismatch(f"{path} does not match the hash in {MANIFEST_FILE}")

    called = set(called_barcodes)
    counted = {k: v for k, v in _molecules(path).items() if v[0]}
    if obs_names is not None and counted and not any(cb in set(obs_names) for cb, _ in counted):
        raise EvidenceMismatch(
            f"none of the {len(counted)} counted molecules in {path} has a barcode in the matrix; "
            "the evidence belongs to another sample or uses another barcode format"
        )
    n = {"virus_best": 0, "host_best": 0, "tie": 0}
    host_only = unaligned = 0
    for (cb, _), (_, virus_as, host_as) in counted.items():
        if cb not in called:
            continue
        if virus_as is None:
            host_only += host_as is not None
            unaligned += host_as is None
            continue
        n[verdict_of(virus_as, host_as)] += 1
    denominator = sum(n.values())
    return SupportFraction(
        n["virus_best"] / denominator if denominator else None,
        None if denominator else REASON_NO_COUNTED,
        n["virus_best"],
        n["host_best"],
        n["tie"],
        host_only,
        unaligned,
    )


@dataclass(frozen=True)
class VirusMeasures:
    """m (every barcode), F with its index-gene count, and S for one sample and virus."""

    m: np.ndarray
    features: FeatureSupport
    support: SupportFraction
    called: np.ndarray


def measure_virus(
    adata: Any,
    virus: str,
    index_genes: Iterable[str],
    called_mask: Any,
    evidence_dir: str | Path | None,
) -> VirusMeasures:
    """The three measures for one sample and virus; S comes only from *evidence_dir* (None = not estimable)."""
    called = _called(adata, called_mask)
    genes = _unique(index_genes)
    obs = np.asarray(adata.obs_names, dtype=str)
    return VirusMeasures(
        molecules_per_cell(adata, genes),
        features_supported(adata, genes, called),
        support_fraction(evidence_dir, virus, set(obs[called]), obs_names=set(obs)),
        called,
    )


def probable_cells(measures: VirusMeasures, cutoffs: Any) -> np.ndarray:
    """Per-cell ``probable`` call: called and ``m >= M and F >= f and S >= s*`` (inclusive).

    *cutoffs* is an :class:`viralscan.evidence_tiers.Cutoffs` grid point. F and S are the sample's values,
    broadcast to every cell; a not-estimable S, or fewer index genes than ``f``, never reaches ``probable``.
    A sample is called when this holds for at least one cell (D8), which is ``assign_tier`` on
    ``(max m over called cells, F, S)``.
    """
    support = measures.support.value
    if (
        support is None
        or measures.features.n_supported < cutoffs.features
        or support < cutoffs.support
    ):
        return np.zeros(measures.m.shape, dtype=bool)
    return measures.called & (
        measures.m >= float(cutoffs.molecules) * (1.0 - INFECTED_RELATIVE_TOLERANCE)
    )

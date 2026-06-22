"""Typed, validated description of a single ViralScan Run.

``RunConfig`` is the one place where the loose, stringly-typed values that arrive
from the CLI (via Snakemake's ``--config``) are coerced and validated. It is the
single write-side checkpoint: ``createconfig`` — the sole writer of
``config.yaml`` — builds a ``RunConfig`` with :meth:`from_snakemake_config` and
serialises it with :meth:`to_yaml`. Every downstream rule reads the file back
with :meth:`from_yaml`, a *trusted* typed load that performs no re-validation,
because nothing but ``createconfig`` ever writes that file.

See ``CONTEXT.md`` ("Run Config") for the vocabulary.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, Union

import yaml

from viralscan.defaults import DEFAULTS

# Strings that represent "unset" once a value has been round-tripped through
# Snakemake's ``--config`` serialisation. ``menu.py`` emits ``key=`` for optional
# paths, which arrives here as an empty string.
_UNSET_STRINGS = {"", "none", "null"}

_TRUE_STRINGS = {"true", "1", "yes", "y"}
_FALSE_STRINGS = {"false", "0", "no", "n", "", "none", "null"}


def _coerce_bool(value: Any) -> bool:
    """Interpret a config value as a bool, robustly.

    The hazard this exists to kill: ``menu.py`` serialises ``--no-visual`` as the
    literal string ``"False"``, and ``bool("False")`` is ``True``. Coerce by
    meaning, never by Python truthiness of a non-empty string.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    s = str(value).strip().lower()
    if s in _TRUE_STRINGS:
        return True
    if s in _FALSE_STRINGS:
        return False
    raise ValueError(f"Cannot interpret {value!r} as a boolean.")


def _opt(value: Any) -> Union[str, None]:
    """Normalise an optional path/string value: unset sentinels become ``None``."""
    if value is None:
        return None
    s = str(value)
    if s.strip().lower() in _UNSET_STRINGS:
        return None
    return s


@dataclass(frozen=True)
class RunConfig:
    """Every parameter the Snakemake rules need to execute one Run.

    All fields carry defaults so the dataclass can be constructed positionally
    free; the :meth:`from_snakemake_config` constructor always sets them
    explicitly. Path-like fields are kept as ``str`` (not ``Path``) so the YAML
    on disk is identical to what the rules read today.
    """

    output: str = ""
    index: str = ""
    transcripts: str = ""
    sample1: str = ""
    sample2: str = ""
    cores: int = 6
    overwrite: str = "yes"
    gtf: Union[str, None] = None
    fasta: Union[str, None] = None
    visual: bool = True
    f1: Union[str, None] = None
    reference: bool = False
    umap: bool = False
    technology: str = "10xv3"
    whitelist: Union[str, None] = None
    multimapping: bool = True
    se_threshold: int = DEFAULTS["se_threshold"]
    detection_threshold: int = DEFAULTS["detection_threshold"]
    min_counts: int = DEFAULTS["min_counts"]
    min_genes: int = DEFAULTS["min_genes"]
    hvg_min_mean: float = DEFAULTS["hvg_min_mean"]
    hvg_max_mean: float = DEFAULTS["hvg_max_mean"]
    hvg_min_disp: float = DEFAULTS["hvg_min_disp"]
    umap_n_neighbors: int = DEFAULTS["umap_n_neighbors"]
    multimap_method: str = DEFAULTS["multimap_method"]
    multimap_pseudocount: float = DEFAULTS["multimap_pseudocount"]
    multimap_primary_call: str = DEFAULTS["multimap_primary_call"]
    multimap_em_max_iter: int = DEFAULTS["multimap_em_max_iter"]
    multimap_em_tol: float = DEFAULTS["multimap_em_tol"]
    cell_types: Union[str, None] = None
    data_cache_dir: Union[str, None] = None
    host_index: Union[str, None] = None
    host_filter_aligner: Union[str, None] = None
    kb_r1: str = ""
    kb_r2: str = ""

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def from_snakemake_config(cls, cfg_in: dict[str, Any]) -> "RunConfig":
        """Build and validate a RunConfig from Snakemake's ``--config`` mapping.

        This is the single coercion + validation checkpoint. Required keys are
        accessed with ``[]`` (a missing one is a programming error worth raising);
        defaulted reporting parameters use ``.get`` against :data:`DEFAULTS`.
        """
        detection_threshold = int(
            cfg_in.get("detection_threshold", DEFAULTS["detection_threshold"])
        )
        if detection_threshold < 1:
            raise ValueError(
                f"detection_threshold must be >= 1, got {detection_threshold}. "
                "A threshold of 0 or below would flag every viral accession as detected."
            )
        multimap_pseudocount = float(
            cfg_in.get("multimap_pseudocount", DEFAULTS["multimap_pseudocount"])
        )
        if multimap_pseudocount <= 0:
            raise ValueError(
                f"multimap_pseudocount must be > 0, got {multimap_pseudocount}."
            )

        host_index = _opt(cfg_in.get("host_index"))
        # Precompute the FASTQ paths kb_count consumes so the Snakefile shell
        # block never needs a conditional. When host subtraction is active the
        # host_filter rule writes filtered FASTQs to a predictable location.
        if host_index:
            out = cfg_in["output"]
            kb_r1 = cfg_in.get("kb_r1") or os.path.join(out, "host_filtered", "R1.fastq.gz")
            kb_r2 = cfg_in.get("kb_r2") or os.path.join(out, "host_filtered", "R2.fastq.gz")
        else:
            kb_r1 = cfg_in.get("kb_r1") or cfg_in["sample1"]
            kb_r2 = cfg_in.get("kb_r2") or cfg_in["sample2"]

        return cls(
            output=cfg_in["output"],
            index=cfg_in["index"],
            transcripts=cfg_in["transcripts"],
            sample1=cfg_in["sample1"],
            sample2=cfg_in["sample2"],
            cores=int(cfg_in.get("cores", 6)),
            overwrite="yes",
            gtf=_opt(cfg_in["gtf"]),
            fasta=_opt(cfg_in["fasta"]),
            visual=_coerce_bool(cfg_in["visual"]),
            f1=_opt(cfg_in["f1"]),
            reference=_coerce_bool(cfg_in["reference"]),
            umap=_coerce_bool(cfg_in["umap"]),
            technology=cfg_in["technology"],
            whitelist=_opt(cfg_in["whitelist"]),
            multimapping=_coerce_bool(cfg_in["multimapping"]),
            se_threshold=int(cfg_in.get("se_threshold", DEFAULTS["se_threshold"])),
            detection_threshold=detection_threshold,
            min_counts=int(cfg_in.get("min_counts", DEFAULTS["min_counts"])),
            min_genes=int(cfg_in.get("min_genes", DEFAULTS["min_genes"])),
            hvg_min_mean=float(cfg_in.get("hvg_min_mean", DEFAULTS["hvg_min_mean"])),
            hvg_max_mean=float(cfg_in.get("hvg_max_mean", DEFAULTS["hvg_max_mean"])),
            hvg_min_disp=float(cfg_in.get("hvg_min_disp", DEFAULTS["hvg_min_disp"])),
            umap_n_neighbors=int(
                cfg_in.get("umap_n_neighbors", DEFAULTS["umap_n_neighbors"])
            ),
            multimap_method=cfg_in.get("multimap_method", DEFAULTS["multimap_method"]),
            multimap_pseudocount=multimap_pseudocount,
            multimap_primary_call=cfg_in.get(
                "multimap_primary_call", DEFAULTS["multimap_primary_call"]
            ),
            multimap_em_max_iter=int(
                cfg_in.get("multimap_em_max_iter", DEFAULTS["multimap_em_max_iter"])
            ),
            multimap_em_tol=float(
                cfg_in.get("multimap_em_tol", DEFAULTS["multimap_em_tol"])
            ),
            cell_types=_opt(cfg_in.get("cell_types")),
            data_cache_dir=_opt(cfg_in.get("data_cache_dir")),
            host_index=host_index,
            host_filter_aligner=_opt(cfg_in.get("host_filter_aligner")),
            kb_r1=kb_r1,
            kb_r2=kb_r2,
        )

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "RunConfig":
        """Trusted typed load of a ``config.yaml`` written by :meth:`to_yaml`.

        Performs no semantic re-validation: the file is only ever produced by
        ``createconfig`` after ``from_snakemake_config`` already validated it.
        Unknown keys (e.g. from a hand-edited file) are ignored.
        """
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        if not isinstance(data, dict):
            raise ValueError(f"Config file {path} did not contain a YAML mapping.")
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in known})

    # ── serialisation ─────────────────────────────────────────────────────
    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_yaml(self, path: Union[str, Path]) -> None:
        with open(path, "w", encoding="utf-8") as out:
            yaml.dump(self.to_dict(), out)

    def to_snakemake_config_args(self) -> "list[str]":
        """Return a ``k=v`` list for Snakemake's ``--config`` derived from all fields.

        Booleans become ``"true"``/``"false"`` (lowercase); ``None`` becomes
        ``""`` (the unset sentinel understood by :func:`_opt`); everything else is
        ``str(v)``.  This is the single authoritative serialisation of
        ``RunConfig`` → Snakemake wire format, eliminating the hand-maintained
        parallel list in ``menu.py``.
        """
        result = []
        for f in fields(self):
            v = getattr(self, f.name)
            if isinstance(v, bool):
                result.append(f"{f.name}={'true' if v else 'false'}")
            elif v is None:
                result.append(f"{f.name}=")
            else:
                result.append(f"{f.name}={v}")
        return result

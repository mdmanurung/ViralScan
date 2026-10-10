"""Typed, validated description of a single ViralScan Run.

``RunConfig`` is the one place where the loose values that arrive from the CLI are
coerced and validated. It is the single write-side checkpoint: ``menu`` builds a
``RunConfig`` with :meth:`from_snakemake_config` and ``menu._write_run_config`` — the
sole writer of ``config.yaml`` — serialises it with :meth:`to_yaml`. Snakemake reads
that typed file (``--configfile``) and every downstream rule reads it back with
:meth:`from_yaml`, a *trusted* typed load that performs no re-validation, because
nothing but ``_write_run_config`` ever writes that file (PLAN ``MECH-C``).

See ``CONTEXT.md`` ("Run Config") for the vocabulary.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, Union

import yaml

from viralscan.defaults import DEFAULTS, HOSTRESPONSE_MAX_SEEDS
from viralscan.sensitivity import CaptureScope

#: Valid ``positive_control_scope`` values, from the one enum that defines them.
CAPTURE_SCOPES: tuple[str, ...] = tuple(scope.value for scope in CaptureScope)

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
    require_chemistry_sanity: bool = False
    whitelist: Union[str, None] = None
    strand: Union[str, None] = None
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
    multimap_molecule_assignments: bool = DEFAULTS["multimap_molecule_assignments"]
    cell_types: Union[str, None] = None
    data_cache_dir: Union[str, None] = None
    host_index: Union[str, None] = None
    host_filter_aligner: Union[str, None] = None
    host_filter_star_params: str = "pinned"
    read_filter: str = "off"
    kb_r1: str = ""
    kb_r2: str = ""
    host_h5ad: Union[str, None] = None
    hostresponse_n_seeds: int = DEFAULTS["hostresponse_n_seeds"]
    hostresponse_n_stab_iter: int = DEFAULTS["hostresponse_n_stab_iter"]
    hostresponse_use_hvg: bool = True
    hostresponse_stab_min_prob: float = DEFAULTS["hostresponse_stab_min_prob"]
    hostresponse_top_n_genes: int = DEFAULTS["hostresponse_top_n_genes"]
    hostresponse_enrichment: bool = False
    hostresponse_enrichment_db: str = "GO_Biological_Process_2023"
    hostresponse_label: str = DEFAULTS["hostresponse_label"]
    hostresponse_depth_match: bool = DEFAULTS["hostresponse_depth_match"]
    hostresponse_control_mito: bool = DEFAULTS["hostresponse_control_mito"]
    hostresponse_differential: bool = DEFAULTS["hostresponse_differential"]
    hostresponse_cv: str = "cell"
    hostresponse_groups: Union[str, None] = None
    hostresponse_cv_folds: int = 5
    hostresponse_cell_type_column: Union[str, None] = None
    hostresponse_cell_types: Union[list[str], None] = None
    hostresponse_panel_in_fold: bool = False
    hostresponse_permutations: int = 0
    hostresponse_permutation_unit: Union[str, None] = None
    hostresponse_permutation_block: Union[str, None] = None
    hostresponse_min_negative_cells: int = 10
    hostresponse_min_groups: int = 2
    # Cell-calling: report viral rates over called cells (primary) + all barcodes
    cell_calling: str = DEFAULTS["cell_calling"]
    called_cells_file: Union[str, None] = None
    emptydrops_fdr: float = DEFAULTS["emptydrops_fdr"]
    emptydrops_lower: int = DEFAULTS["emptydrops_lower"]
    emptydrops_niters: int = DEFAULTS["emptydrops_niters"]
    emptydrops_seed: int = DEFAULTS["emptydrops_seed"]
    knee_min_umi: float = DEFAULTS["knee_min_umi"]
    cell_caller_rscript: str = DEFAULTS["cell_caller_rscript"]
    # Positive control: a spike-in at a known molecule count. Its presence is
    # the only thing that turns a negative into a certifiable negative, because
    # it measures the k-mer capture term that depth alone cannot supply.
    positive_control_gene: Union[str, None] = None
    positive_control_expected_molecules: Union[float, None] = None
    positive_control_receipt: Union[str, None] = None
    # What the control may certify (SENS-CORR-01). None on a configured control
    # means the legacy behaviour, panel_mechanics: it certifies no virus.
    positive_control_scope: Union[str, None] = None
    positive_control_virus_key: Union[str, None] = None
    require_positive_control: bool = DEFAULTS["require_positive_control"]
    anellovirus_gene_ids: bool = DEFAULTS["anellovirus_gene_ids"]
    # Layer 2: gene-programme inference for viruses layer 1 detected
    gene_programs: bool = DEFAULTS["gene_programs"]
    programme_min_breadth: int = DEFAULTS["programme_min_breadth"]
    programme_latent_min_breadth: int = DEFAULTS["programme_latent_min_breadth"]
    programme_min_umi: float = DEFAULTS["programme_min_umi"]
    # Anellovirus alignment branch (ANDET-09). anello_index is resolved from the
    # kb index's directory (anello_star/) unless given; None means no index.
    anello_align: bool = DEFAULTS["anello_align"]
    anello_index: Union[str, None] = None
    # Read-level evidence with host confirmation for detected anellovirus genera
    # (ANDET-04). viral_fasta is resolved to <index dir>/viral.fa unless given.
    auto_evidence: bool = DEFAULTS["auto_evidence"]
    host_fasta: Union[str, None] = None
    viral_fasta: Union[str, None] = None

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def from_snakemake_config(cls, cfg_in: dict[str, Any]) -> RunConfig:
        """Build and validate a RunConfig from Snakemake's ``--config`` mapping.

        This is the single coercion + validation checkpoint. Required keys are
        accessed with ``[]`` (a missing one is a programming error worth raising);
        defaulted reporting parameters use ``.get`` against :data:`DEFAULTS`.
        """

        def hr_int(name: str, default: int) -> int:
            value = cfg_in.get(f"hostresponse_{name}")
            return default if value is None else int(value)

        hr_cv = _opt(cfg_in.get("hostresponse_cv")) or "cell"
        hr_groups = _opt(cfg_in.get("hostresponse_groups"))
        hr_folds = hr_int("cv_folds", 5)
        hr_permutations = hr_int("permutations", 0)
        hr_min_negative = hr_int("min_negative_cells", 10)
        hr_min_groups = hr_int("min_groups", 2)
        hr_n_seeds = hr_int("n_seeds", DEFAULTS["hostresponse_n_seeds"])
        if not 1 <= hr_n_seeds <= HOSTRESPONSE_MAX_SEEDS:
            raise ValueError(
                f"hostresponse_n_seeds must be between 1 and {HOSTRESPONSE_MAX_SEEDS} "
                f"(the fixed seed list), got {hr_n_seeds}."
            )
        hr_unit = _opt(cfg_in.get("hostresponse_permutation_unit"))
        hr_types = cfg_in.get("hostresponse_cell_types")
        if isinstance(hr_types, str):
            hr_types = json.loads(hr_types) if _opt(hr_types) else None
        if hr_types is not None and (
            not isinstance(hr_types, list)
            or not hr_types
            or any(not isinstance(v, str) or not v.strip() for v in hr_types)
            or len(set(hr_types)) != len(hr_types)
        ):
            raise ValueError("hostresponse_cell_types must be a nonempty list of distinct names.")
        if hr_cv not in {"cell", "group"} or (hr_cv == "group" and not hr_groups):
            raise ValueError("hostresponse_cv must be cell, or group with hostresponse_groups.")
        if hr_folds < 2 or hr_permutations < 0 or min(hr_min_negative, hr_min_groups) < 1:
            raise ValueError(
                "Host-response folds/support must be positive (folds >= 2), permutations >= 0."
            )
        if hr_unit not in {None, "cell_within_block", "group"}:
            raise ValueError("Invalid hostresponse_permutation_unit.")
        if hr_permutations and (hr_cv != "group" or hr_unit is None):
            raise ValueError("Structured permutations require group CV and an explicit unit.")
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
            raise ValueError(f"multimap_pseudocount must be > 0, got {multimap_pseudocount}.")
        # An EM iteration budget of zero makes `range(1, max_iter + 1)` empty, so
        # em_gene_abundances/em_cell_abundances return their unconverged seed
        # weights while the H5AD still records multimap_method as an EM method.
        # Mass conservation still holds, so MoleculeAudit.validate cannot catch it.
        multimap_em_max_iter = int(
            cfg_in.get("multimap_em_max_iter", DEFAULTS["multimap_em_max_iter"])
        )
        if multimap_em_max_iter < 1:
            raise ValueError(
                f"multimap_em_max_iter must be >= 1, got {multimap_em_max_iter}. "
                "A budget below one returns unconverged seed weights labelled as an EM result."
            )
        multimap_em_tol = float(cfg_in.get("multimap_em_tol", DEFAULTS["multimap_em_tol"]))
        if not multimap_em_tol > 0:
            raise ValueError(f"multimap_em_tol must be > 0, got {multimap_em_tol}.")
        multimap_primary_call = cfg_in.get(
            "multimap_primary_call", DEFAULTS["multimap_primary_call"]
        )
        if multimap_primary_call != "selected-method":
            raise ValueError(
                "Pre-v3 multimap primary-call modes are scientifically incompatible with "
                "the v3 count contract. Rebuild with multimap_primary_call=selected-method."
            )

        # Positive control. The gene and its expected molecule count are only
        # meaningful together: a gene with no expected count cannot establish a
        # capture ratio, and an expected count with no gene names nothing to
        # measure. Accepting one without the other would let a run look
        # controlled while measuring nothing.
        positive_control_gene = _opt(cfg_in.get("positive_control_gene"))
        positive_control_expected = cfg_in.get("positive_control_expected_molecules")
        if positive_control_expected is not None and str(
            positive_control_expected
        ).strip().lower() not in {"", "none", "null"}:
            positive_control_expected = float(positive_control_expected)
            if positive_control_expected <= 0:
                raise ValueError(
                    "positive_control_expected_molecules must be > 0, got "
                    f"{positive_control_expected}. A control at zero abundance "
                    "cannot demonstrate recovery."
                )
        else:
            positive_control_expected = None
        if bool(positive_control_gene) != (positive_control_expected is not None):
            raise ValueError(
                "positive_control_gene and positive_control_expected_molecules "
                "must be supplied together. Got "
                f"gene={positive_control_gene!r}, "
                f"expected={positive_control_expected!r}. A control with no known "
                "abundance cannot establish the k-mer capture term, which is the "
                "only thing that makes a negative certifiable "
                "(see viralscan.sensitivity)."
            )
        positive_control_scope = _opt(cfg_in.get("positive_control_scope"))
        positive_control_virus_key = _opt(cfg_in.get("positive_control_virus_key"))
        if positive_control_scope is not None and positive_control_scope not in CAPTURE_SCOPES:
            raise ValueError(
                f"positive_control_scope must be one of {', '.join(CAPTURE_SCOPES)}, "
                f"got {positive_control_scope!r}."
            )
        if (positive_control_scope or positive_control_virus_key) and not positive_control_gene:
            raise ValueError(
                "positive_control_scope / positive_control_virus_key require a positive "
                "control: supply --positive-control-gene and --positive-control-molecules."
            )
        if positive_control_virus_key and positive_control_scope is None:
            raise ValueError(
                "positive_control_virus_key requires positive_control_scope: a target "
                "with no declared scope would be read as the legacy panel_mechanics "
                "control, which certifies no virus."
            )
        if positive_control_scope == "virus_key":
            raise ValueError(
                "positive_control_scope=virus_key needs an approved transfer calibration "
                "(VAL-RA-CAL), which this release does not provide. Use exact_sequence "
                "for one declared target, or panel_mechanics."
            )
        if positive_control_scope == "exact_sequence" and not positive_control_virus_key:
            raise ValueError(
                "positive_control_scope=exact_sequence requires positive_control_virus_key "
                "naming the virus row the control was measured on."
            )
        if positive_control_scope == "panel_mechanics" and positive_control_virus_key:
            raise ValueError(
                "positive_control_virus_key is meaningless with panel_mechanics, which "
                "certifies no virus."
            )
        programme_min_breadth = int(
            cfg_in.get("programme_min_breadth", DEFAULTS["programme_min_breadth"])
        )
        if programme_min_breadth < 1:
            raise ValueError(
                f"programme_min_breadth must be >= 1, got {programme_min_breadth}. "
                "Breadth is counted in distinct non-overlapping overlap groups, "
                "not genes, because EBV's latent and lytic ORFs cross-map through "
                "shared exonic sequence; a breadth of 0 would call every cell "
                "productive on no evidence."
            )
        programme_latent_min_breadth = int(
            cfg_in.get("programme_latent_min_breadth", DEFAULTS["programme_latent_min_breadth"])
        )
        if programme_latent_min_breadth < 1:
            raise ValueError(
                f"programme_latent_min_breadth must be >= 1, got {programme_latent_min_breadth}."
            )
        programme_min_umi = float(cfg_in.get("programme_min_umi", DEFAULTS["programme_min_umi"]))
        if programme_min_umi < 0:
            raise ValueError(f"programme_min_umi must be >= 0, got {programme_min_umi}.")
        gene_programs = _coerce_bool(
            cfg_in.get("gene_programs", DEFAULTS["gene_programs"])
            if cfg_in.get("gene_programs") is not None
            else DEFAULTS["gene_programs"]
        )
        require_positive_control = _coerce_bool(
            cfg_in.get("require_positive_control", DEFAULTS["require_positive_control"])
            if cfg_in.get("require_positive_control") is not None
            else DEFAULTS["require_positive_control"]
        )
        if require_positive_control and not positive_control_gene:
            raise ValueError(
                "require_positive_control is set but no positive control was "
                "supplied. Supply --positive-control-gene and "
                "--positive-control-molecules, or unset the requirement. Refusing "
                "to run is deliberate: a negative result with no control is not "
                "evidence of absence."
            )

        anello_align = _coerce_bool(
            cfg_in.get("anello_align")
            if cfg_in.get("anello_align") is not None
            else DEFAULTS["anello_align"]
        )
        anello_index = None
        if anello_align:
            from viralscan.anello_align import resolve_index

            anello_index = _opt(cfg_in.get("anello_index")) or resolve_index(cfg_in["index"])

        auto_evidence = _coerce_bool(
            cfg_in.get("auto_evidence")
            if cfg_in.get("auto_evidence") is not None
            else DEFAULTS["auto_evidence"]
        )
        host_fasta = _opt(cfg_in.get("host_fasta"))
        viral_fasta = _opt(cfg_in.get("viral_fasta"))
        if auto_evidence:
            # Fail closed (ANDET-04): without both FASTAs `viralscan evidence` cannot
            # confirm reads against the host. Only <index dir>/viral.fa is viral-only;
            # never fall back to combined.fa / reference.prepared.fa, which hold host
            # sequence that evidence would label VIRUS.
            viral_fasta = viral_fasta or str(Path(cfg_in["index"]).resolve().parent / "viral.fa")
            if not host_fasta:
                raise ValueError("auto_evidence requires host_fasta (--host-fasta).")
            if not Path(viral_fasta).is_file():
                raise ValueError(
                    f"auto_evidence needs a viral-only FASTA: {viral_fasta} does not exist. "
                    "Pass --viral-fasta, or use a reference built by "
                    "scripts/build_bundled_panel_ref.py, which writes viral.fa next to the index."
                )
            if not _coerce_bool(cfg_in["multimapping"]):
                raise ValueError(
                    "auto_evidence needs the resolved BUS text that multimapping writes; "
                    "drop --no-multimapping."
                )

        host_index = _opt(cfg_in.get("host_index"))
        read_filter = _opt(cfg_in.get("read_filter")) or "off"
        # Precompute the FASTQ paths kb_count consumes so the Snakefile shell
        # block never needs a conditional. The read filter (DEF-01) runs last,
        # after any host filter, so its output wins.
        if read_filter != "off":
            out = cfg_in["output"]
            kb_r1 = cfg_in.get("kb_r1") or os.path.join(out, "read_filtered", "R1.fastq.gz")
            kb_r2 = cfg_in.get("kb_r2") or os.path.join(out, "read_filtered", "R2.fastq.gz")
        elif host_index:
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
            require_chemistry_sanity=_coerce_bool(cfg_in.get("require_chemistry_sanity", False)),
            whitelist=_opt(cfg_in["whitelist"]),
            strand=_opt(cfg_in.get("strand")),
            multimapping=_coerce_bool(cfg_in["multimapping"]),
            se_threshold=int(cfg_in.get("se_threshold", DEFAULTS["se_threshold"])),
            detection_threshold=detection_threshold,
            min_counts=int(cfg_in.get("min_counts", DEFAULTS["min_counts"])),
            min_genes=int(cfg_in.get("min_genes", DEFAULTS["min_genes"])),
            hvg_min_mean=float(cfg_in.get("hvg_min_mean", DEFAULTS["hvg_min_mean"])),
            hvg_max_mean=float(cfg_in.get("hvg_max_mean", DEFAULTS["hvg_max_mean"])),
            hvg_min_disp=float(cfg_in.get("hvg_min_disp", DEFAULTS["hvg_min_disp"])),
            umap_n_neighbors=int(cfg_in.get("umap_n_neighbors", DEFAULTS["umap_n_neighbors"])),
            multimap_method=cfg_in.get("multimap_method", DEFAULTS["multimap_method"]),
            multimap_pseudocount=multimap_pseudocount,
            multimap_primary_call=multimap_primary_call,
            multimap_em_max_iter=multimap_em_max_iter,
            multimap_em_tol=multimap_em_tol,
            multimap_molecule_assignments=_coerce_bool(
                cfg_in.get("multimap_molecule_assignments")
                if cfg_in.get("multimap_molecule_assignments") is not None
                else DEFAULTS["multimap_molecule_assignments"]
            ),
            cell_types=_opt(cfg_in.get("cell_types")),
            data_cache_dir=_opt(cfg_in.get("data_cache_dir")),
            host_index=host_index,
            host_filter_aligner=_opt(cfg_in.get("host_filter_aligner")),
            host_filter_star_params=_opt(cfg_in.get("host_filter_star_params")) or "pinned",
            read_filter=read_filter,
            kb_r1=kb_r1,
            kb_r2=kb_r2,
            host_h5ad=_opt(cfg_in.get("host_h5ad")),
            hostresponse_cv=hr_cv,
            hostresponse_groups=hr_groups,
            hostresponse_cv_folds=hr_folds,
            hostresponse_cell_type_column=_opt(cfg_in.get("hostresponse_cell_type_column")),
            hostresponse_cell_types=hr_types,
            hostresponse_panel_in_fold=_coerce_bool(
                cfg_in.get("hostresponse_panel_in_fold") or False
            ),
            hostresponse_permutations=hr_permutations,
            hostresponse_permutation_unit=hr_unit,
            hostresponse_permutation_block=_opt(cfg_in.get("hostresponse_permutation_block")),
            hostresponse_min_negative_cells=hr_min_negative,
            hostresponse_min_groups=hr_min_groups,
            hostresponse_n_seeds=hr_n_seeds,
            hostresponse_n_stab_iter=int(
                cfg_in.get("hostresponse_n_stab_iter") or DEFAULTS["hostresponse_n_stab_iter"]
            ),
            hostresponse_use_hvg=_coerce_bool(
                cfg_in.get("hostresponse_use_hvg", True)
                if cfg_in.get("hostresponse_use_hvg") is not None
                else True
            ),
            hostresponse_stab_min_prob=float(
                cfg_in.get("hostresponse_stab_min_prob") or DEFAULTS["hostresponse_stab_min_prob"]
            ),
            hostresponse_top_n_genes=int(
                cfg_in.get("hostresponse_top_n_genes") or DEFAULTS["hostresponse_top_n_genes"]
            ),
            hostresponse_enrichment=_coerce_bool(cfg_in.get("hostresponse_enrichment", False)),
            hostresponse_enrichment_db=cfg_in.get("hostresponse_enrichment_db")
            or "GO_Biological_Process_2023",
            hostresponse_label=cfg_in.get("hostresponse_label") or DEFAULTS["hostresponse_label"],
            hostresponse_depth_match=_coerce_bool(
                cfg_in.get("hostresponse_depth_match", DEFAULTS["hostresponse_depth_match"])
            ),
            hostresponse_control_mito=_coerce_bool(
                cfg_in.get("hostresponse_control_mito", DEFAULTS["hostresponse_control_mito"])
            ),
            hostresponse_differential=_coerce_bool(
                cfg_in.get("hostresponse_differential", DEFAULTS["hostresponse_differential"])
            ),
            cell_calling=cfg_in.get("cell_calling") or DEFAULTS["cell_calling"],
            called_cells_file=_opt(cfg_in.get("called_cells_file")),
            emptydrops_fdr=float(cfg_in.get("emptydrops_fdr") or DEFAULTS["emptydrops_fdr"]),
            emptydrops_lower=int(cfg_in.get("emptydrops_lower") or DEFAULTS["emptydrops_lower"]),
            emptydrops_niters=int(cfg_in.get("emptydrops_niters") or DEFAULTS["emptydrops_niters"]),
            # `or` would swallow a legitimate seed of 0, so test for absence.
            emptydrops_seed=int(
                DEFAULTS["emptydrops_seed"]
                if cfg_in.get("emptydrops_seed") is None
                else cfg_in["emptydrops_seed"]
            ),
            knee_min_umi=float(cfg_in.get("knee_min_umi") or DEFAULTS["knee_min_umi"]),
            cell_caller_rscript=cfg_in.get("cell_caller_rscript")
            or DEFAULTS["cell_caller_rscript"],
            positive_control_gene=positive_control_gene,
            positive_control_receipt=_opt(cfg_in.get("positive_control_receipt")),
            positive_control_expected_molecules=positive_control_expected,
            positive_control_scope=positive_control_scope,
            positive_control_virus_key=positive_control_virus_key,
            require_positive_control=require_positive_control,
            gene_programs=gene_programs,
            programme_min_breadth=programme_min_breadth,
            programme_latent_min_breadth=programme_latent_min_breadth,
            programme_min_umi=programme_min_umi,
            anello_align=anello_align,
            anello_index=anello_index,
            auto_evidence=auto_evidence,
            host_fasta=host_fasta,
            viral_fasta=viral_fasta,
            anellovirus_gene_ids=_coerce_bool(
                cfg_in.get("anellovirus_gene_ids", DEFAULTS["anellovirus_gene_ids"])
                if cfg_in.get("anellovirus_gene_ids") is not None
                else DEFAULTS["anellovirus_gene_ids"]
            ),
        )

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> RunConfig:
        """Trusted typed load of a ``config.yaml`` written by :meth:`to_yaml`.

        Performs no semantic re-validation: the file is only ever produced by
        ``menu._write_run_config`` after ``from_snakemake_config`` already validated it.
        Unknown keys (e.g. from a hand-edited file) are ignored.
        """
        with open(path, encoding="utf-8") as f:
            data = yaml.safe_load(f)
        if not isinstance(data, dict):
            raise ValueError(f"Config file {path} did not contain a YAML mapping.")
        known = {f.name for f in fields(cls)}
        kwargs = {k: v for k, v in data.items() if k in known}
        # `output` is an implicit "must end with a separator" contract (downstream
        # f-strings do f"{output}log/…"). A hand-edited or subcommand-loaded config
        # without the trailing slash would silently yield `…outputlog/…`; normalize
        # before constructing (RunConfig is frozen).
        out = kwargs.get("output")
        if isinstance(out, str) and out and not out.endswith(os.sep):
            kwargs["output"] = out + os.sep
        return cls(**kwargs)

    # ── serialisation ─────────────────────────────────────────────────────
    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_yaml(self, path: Union[str, Path]) -> None:
        with open(path, "w", encoding="utf-8") as out:
            yaml.dump(self.to_dict(), out)

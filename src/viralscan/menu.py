"""
This file is the backbone of the framework. It checks users' input and calls
the snakemake workflow. It handles the Argument Parser, showing the help function.
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, NoReturn, Optional

from viralscan.cli.parsers import (  # noqa: F401  (private builders re-exported for tests and callers)
    _add_hostresponse_evaluation_args,
    _add_verbosity_args,
    _build_check_chemistry_parser,
    _build_check_whitelist_parser,
    _build_data_parser,
    _build_doctor_parser,
    _build_evidence_parser,
    _build_hostresponse_parser,
    _build_ref_parser,
    _build_rerun_multimap_parser,
    _build_rerun_programs_parser,
    _build_validate_run_parser,
    build_parser,
)
from viralscan.defaults import (
    DEFAULTS,
)
from viralscan.run_safety import (
    RUN_MANIFEST,
    clear_run_complete,
    record_strand_inference,
    recorded_strand_inference,
    restamp_run_complete,
    write_run_complete,
)
from viralscan.runconfig import RunConfig
from viralscan.utils import configure_logging, split_comma_paths

log = logging.getLogger("viralscan")


REQUIRED_TOOLS = ("kb", "snakemake")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def _swap_multimap_layer(adata_path: Path, new_method: str) -> bool:
    """Swap counts_corrected in a multimap h5ad to a pre-stored layer for new_method.

    All non-EM layers are computed on every multimap run and stored as named layers.
    This updates the complete v3 matrix and its two compositional layers without
    touching the BUS file. Pre-v3 H5AD files are never numerically migrated.

    Returns True on success, False when the target layer is absent (e.g. h5ad was
    produced by an older ViralScan version without pre-stored layers).  Callers
    should fall back to a full multimap rerun when False is returned.
    """
    import anndata as _ad  # heavy import; keep lazy

    layer_name = {
        "equal": "counts_multimap_equal",
        "host-conservative": "counts_multimap_host_conservative",
        "unique-weighted": "counts_multimap_unique_weighted",
        "sibling-weighted": "counts_multimap_sibling_weighted",
    }[new_method]
    # counts_host_viral_selected is method-dependent (0 for host-conservative,
    # equal/weighted shares otherwise), so it must be swapped from the matching
    # per-method layer too; leaving it behind mislabels the evidence tiers.
    host_viral_layer_name = {
        "equal": "counts_host_viral_selected_equal",
        "host-conservative": "counts_host_viral_selected_host_conservative",
        "unique-weighted": "counts_host_viral_selected_unique_weighted",
        "sibling-weighted": "counts_host_viral_selected_sibling_weighted",
    }[new_method]
    adata = _ad.read_h5ad(str(adata_path))
    if new_method == "sibling-weighted" and not (
        adata.uns.get("multimap_diagnostics", {}).get("sibling_identity_available", False)
    ):
        return False
    if (
        adata.uns.get("count_schema_version") != "3.0.0"
        or "counts_unique" not in adata.layers
        or layer_name not in adata.layers
        or host_viral_layer_name not in adata.layers
    ):
        return False
    selected = adata.layers[layer_name]
    adata.layers["counts_ambiguous_allocated"] = selected
    adata.layers["counts_corrected"] = selected
    adata.X = adata.layers["counts_unique"] + selected
    adata.layers["counts_combined"] = adata.X
    adata.layers["counts_host_viral_selected"] = adata.layers[host_viral_layer_name]
    adata.uns["multimap_method"] = new_method
    adata.write_h5ad(str(adata_path))
    return True


def _run_rerun_programs(args: argparse.Namespace) -> None:
    """Run layer 2 in place on a completed run directory."""
    import yaml as _yaml

    from viralscan.kb_outputs import KbCountOutputs
    from viralscan.runconfig import RunConfig
    from viralscan.scripts import gene_programs as gp_script

    configure_logging(verbose=args.verbose, quiet=args.quiet)

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.is_dir():
        _die(f"Run directory does not exist: {run_dir}")
    config_path = run_dir / "config.yaml"
    if not config_path.is_file():
        _die(
            f"No config.yaml in {run_dir}. Point --run-dir at a viralscan run "
            "directory, not its parent."
        )
    if not (run_dir / "log" / "detection.done").exists():
        _die(
            f"No log/detection.done in {run_dir}, so layer 1 has not finished and "
            "its viral_summary.tsv does not exist yet. Run viralscan first, or "
            "enable --gene-programs on the original run."
        )
    if args.programme_min_breadth is not None and args.programme_min_breadth < 1:
        _die("--programme-min-breadth must be >= 1.")
    if args.programme_latent_min_breadth is not None and args.programme_latent_min_breadth < 1:
        _die("--programme-latent-min-breadth must be >= 1.")
    if args.programme_min_umi is not None and args.programme_min_umi < 0:
        _die("--programme-min-umi must be >= 0.")

    with config_path.open() as handle:
        loaded = _yaml.safe_load(handle)
    payload = dict(loaded or {})
    # Layer 2 reads the H5AD and layer 1's summary. It never re-derives counts
    # from the BUS, so the pre-v3 primary-call guard — which exists to stop
    # multimap re-derivation producing a matrix the v3 count contract cannot
    # describe — does not apply here. Override it explicitly, with a log line,
    # rather than letting a legacy run directory block a read-only analysis or
    # silently bypassing the guard for the paths where it does matter.
    if payload.get("multimap_primary_call") not in (None, "selected-method"):
        log.info(
            "Run directory declares multimap_primary_call=%r. Overriding to "
            "'selected-method' for gene-programme inference only: layer 2 reads "
            "the existing count matrix and does not re-derive counts from the BUS.",
            payload["multimap_primary_call"],
        )
        payload["multimap_primary_call"] = "selected-method"
    # An explicit flag overrides the run's own config.yaml; absent, the run's value stands.
    for key in ("programme_min_breadth", "programme_latent_min_breadth", "programme_min_umi"):
        if getattr(args, key) is not None:
            payload[key] = getattr(args, key)
    config = RunConfig.from_snakemake_config(payload)
    outputs = KbCountOutputs.from_config_output(config.output)
    adata_path = Path(str(outputs.adata_multimap))
    if not adata_path.is_file():
        _die(
            f"Multimap H5AD not found at {adata_path}. Layer 2 reads the "
            "counts_unique_viral layer, which only exists when multimapping ran."
        )

    gp_script.config = config
    gp_script.output = str(run_dir)
    gp_script.main(
        str(adata_path),
        str(run_dir / "results" / "viral_summary.tsv"),
        str(run_dir / "log" / "gene_programs.done"),
    )
    restamp_run_complete(run_dir.parent)
    log.info("rerun-programs complete in %s", run_dir)


def _write_run_config(run_config: RunConfig, *, keep_mtime: bool = False) -> Path:
    """The one writer of a sample's ``config.yaml`` (PLAN ``MECH-C``).

    ``main`` and ``rerun-multimap`` both end here: the validated :class:`RunConfig` is written to
    ``<output>config.yaml`` and Snakemake reads that file (``--configfile``), so there is no
    ``k=v`` wire to coerce back from strings. The file is rewritten only when its content
    changes, and ``keep_mtime`` restores its old mtime: ``config.yaml`` is an input of
    ``kb_count``, ``host_filter`` and ``analysis``, so a fresh mtime would re-run kallisto and
    STAR on the FASTQs for what should be a layer swap (SW-22).
    """
    import yaml as _yaml

    path = Path(run_config.output) / "config.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.parent / "log").mkdir(exist_ok=True)
    if path.is_file() and path.read_text(encoding="utf-8") == _yaml.dump(run_config.to_dict()):
        return path
    stat = path.stat() if path.is_file() else None
    run_config.to_yaml(path)
    if keep_mtime and stat is not None:
        os.utime(path, (stat.st_atime, stat.st_mtime))
    return path


def _snakemake_run_command(snakefile_path: str, cores: int, config_file: str | Path) -> list[str]:
    """The one ``snakemake`` invocation used by ``main`` and ``rerun-multimap`` (SW-14).

    The ``all`` target precedes ``--quiet``: snakemake 9 declares
    ``--quiet [{all,host,progress,reason,rules} ...]``, which consumes a
    following ``all``. No ``--use-conda``: the rules need no conda environments,
    and the flag made a run fail on hosts without ``conda`` on PATH.
    """
    return [
        "snakemake",
        "--snakefile",
        snakefile_path,
        "all",
        "--cores",
        str(cores),
        "--quiet",
        "--configfile",
        str(config_file),
    ]


def _snakemake_unlock_command(snakefile_path: str, config_file: str | Path) -> list[str]:
    """``snakemake --unlock`` for the same Snakefile and config as the run just finished."""
    return [
        "snakemake",
        "--snakefile",
        snakefile_path,
        "--unlock",
        "--configfile",
        str(config_file),
    ]


def _run_rerun_multimap(args: argparse.Namespace) -> None:
    """Re-run multimap → detection → umap with a different method, skipping kb_count."""
    import yaml as _yaml

    from viralscan.kb_outputs import KbCountOutputs
    from viralscan.runconfig import RunConfig

    configure_logging(verbose=args.verbose, quiet=args.quiet)
    _check_required_tools()

    source_dir = Path(args.run_dir).resolve()
    output_dir = Path(args.output).resolve()
    if not source_dir.is_dir():
        _die(f"Source run directory does not exist: {source_dir}")
    if output_dir.exists() and any(output_dir.iterdir()):
        _die(f"New rerun output must be empty or absent: {output_dir}")
    if output_dir == source_dir or source_dir in output_dir.parents:
        _die("Rerun output must be separate from and outside the source run directory.")

    source_config_rels = sorted(
        path.relative_to(source_dir)
        for path in source_dir.glob("*/config.yaml")
        if (path.parent / "log" / "multimap.done").exists()
    )
    if not source_config_rels:
        _die(
            f"No samples with a completed multimap step found under {source_dir}. "
            "Run viralscan first so the multimap checkpoint exists."
        )

    if output_dir.exists():
        output_dir.rmdir()
    shutil.copytree(source_dir, output_dir)
    # The copy is a new, unfinished run; the source's marker must not vouch for it.
    clear_run_complete(output_dir)

    new_method = args.multimap_method
    snakefile_path = os.path.join(os.path.dirname(__file__), "Snakefile")

    # Find all sample subdirs with a completed multimap checkpoint.
    sample_configs = [output_dir / rel for rel in source_config_rels]

    log.info(
        "Switching %d sample(s) to --multimap-method %s",
        len(sample_configs),
        new_method,
    )

    for config_yaml_path in sample_configs:
        sample_dir = config_yaml_path.parent
        rel = config_yaml_path.parent.relative_to(output_dir)
        log.info("[%s] Re-running multimap …", rel)

        with open(config_yaml_path) as f:
            cfg = _yaml.safe_load(f)

        # Rerun re-runs detection, which calls the same fail-closed cell caller,
        # so it needs the same preflight the quant path gets. The method comes
        # from the source run's config rather than from argv.
        _check_cell_caller_tools(
            _resolve_cell_calling(
                str(cfg.get("cell_calling") or DEFAULTS["cell_calling"]).lower(),
                cfg.get("called_cells_file"),
            ),
            cfg.get("cell_caller_rscript") or DEFAULTS["cell_caller_rscript"],
        )

        # Update multimap parameters in the working copy. The trailing separator
        # matters: the Snakefile builds paths as f"{config['output']}log/...".
        cfg["output"] = str(output_dir / rel) + os.sep
        cfg["multimap_method"] = new_method
        cfg["cores"] = args.cores
        if args.multimap_em_max_iter is not None:
            cfg["multimap_em_max_iter"] = args.multimap_em_max_iter
        if args.multimap_em_tol is not None:
            cfg["multimap_em_tol"] = args.multimap_em_tol

        use_em = new_method in {"em-global", "em-cell"}
        swapped = False
        if not use_em and (sample_dir / "results" / "molecule_assignments.tsv.gz").exists():
            # Its per-molecule weights belong to the old method and a layer swap never re-reads
            # the BUS file, so only a full multimap pass can rewrite them (SW-03).
            log.info("[%s] molecule_assignments.tsv.gz present: full multimap rerun", rel)
            use_em = True

        if not use_em:
            # Fast path: layers pre-stored; just overwrite counts_corrected in h5ad.
            # Resolve the h5ad inside the *copy*. The copied config.yaml still
            # names the source run as `output` until it is rewritten below, so
            # reading the path from it swapped the SOURCE run's h5ad (SW-19).
            adata_path = KbCountOutputs(sample_dir).adata_multimap
            if not adata_path.exists():
                log.warning(
                    "[%s] No multimap h5ad at %s — falling back to full multimap rerun",
                    rel,
                    adata_path,
                )
                use_em = True  # trigger the full-rerun path below
            else:
                swapped = _swap_multimap_layer(adata_path, new_method)
                if swapped:
                    log.info("[%s] Layer swapped (instant — no bus-file reprocessing)", rel)
                else:
                    log.warning(
                        "[%s] Pre-stored layer not found in h5ad (older run?) — full multimap rerun",
                        rel,
                    )
                    use_em = True

        # Persist updated config (after any use_em adjustment), keeping its
        # mtime (SW-22): config.yaml is an input of kb_count, host_filter and
        # analysis, so a fresh mtime made snakemake re-run kallisto and STAR on
        # the FASTQs. The sentinels below choose what reruns.
        run_config = RunConfig.from_snakemake_config(cfg)
        _write_run_config(run_config, keep_mtime=True)
        _backfill_identity_table(run_config)

        # Drop sentinels so snakemake re-runs the right rules.
        # hostresponse is method-dependent (its per-virus set depends on the
        # allocation layer) but was omitted here, so it re-ran only if snakemake
        # happened to judge it stale by mtime.
        sentinels = ["log/detection.done", "log/umap.done", "log/hostresponse.done"]
        if not swapped:
            # EM or fallback: re-run multimap itself too.
            sentinels.insert(0, "log/multimap.done")
        for sentinel in sentinels:
            p = sample_dir / sentinel
            if p.exists():
                p.unlink()

        subprocess.run(
            _snakemake_run_command(snakefile_path, args.cores, config_yaml_path), check=True
        )
        subprocess.run(_snakemake_unlock_command(snakefile_path, config_yaml_path), check=True)

    if not _rewrite_run_manifest(output_dir, source_dir=source_dir, new_method=new_method):
        log.warning(
            "no %s at %s to re-stamp; provenance still names the source run's allocation method",
            RUN_MANIFEST,
            output_dir,
        )
    else:
        write_run_complete(output_dir)
    log.info("rerun-multimap complete in new result directory %s.", output_dir)


def _backfill_identity_table(run_config: RunConfig) -> None:
    """Give a pre-35940ec Run its ``virus_identity.tsv`` so ``analysis`` does not rerun."""
    from viralscan.virus_identity import backfill_identity_table

    try:
        backfill_identity_table(run_config)
    except ValueError as exc:
        _die(f"Cannot build the Virus Identity table for {run_config.output}: {exc}")


def _rewrite_run_manifest(run_root: Path, *, source_dir: Path, new_method: str) -> bool:
    """Re-stamp the run manifest for the method the rerun just applied.

    The manifest sits at the **root** of the result tree, one level above the
    per-sample directories: ``prepare_output_directory`` writes it into
    ``--output``, and ``_write_run_config`` then creates a subdirectory per sample
    beneath it. A completed run of ``viralscan -o out`` yields
    ``out/run_manifest.json`` alongside ``out/<sample>/config.yaml``.

    This was briefly changed to rewrite per sample, on a review finding that
    claimed the manifest lived beside ``config.yaml``. It does not. The finding
    was wrong and the change was a regression: it moved the rewrite to a path
    that never exists, so the manifest stopped being updated at all. Confirmed by
    running the pipeline end to end (SW-10) and listing the result tree.
    """
    import hashlib as _hashlib
    import json as _json

    manifest_path = run_root / RUN_MANIFEST
    if not manifest_path.is_file():
        return False
    manifest = _json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["derived_from"] = str(source_dir)
    manifest["parent_run_fingerprint"] = manifest.get("run_fingerprint")
    manifest["allocation_method"] = new_method
    manifest["completion_marker"] = True
    payload = {key: value for key, value in manifest.items() if key != "run_fingerprint"}
    manifest["run_fingerprint"] = _hashlib.sha256(
        _json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    staging = run_root / f".{RUN_MANIFEST}.tmp"
    staging.write_text(_json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    staging.replace(manifest_path)
    return True


def _hostresponse_evaluation_kwargs(args: argparse.Namespace, cfg: RunConfig) -> dict[str, Any]:
    mapping = {
        "cv": "cv_mode",
        "groups": "groups_column",
        "cv_folds": "cv_folds",
        "cell_type_column": "cell_type_column",
        "cell_types": "cell_types",
        "panel_in_fold": "panel_in_fold",
        "permutations": "permutations",
        "permutation_unit": "permutation_unit",
        "permutation_block": "permutation_block",
        "min_negative_cells": "min_negative_cells",
        "min_groups": "min_groups",
    }
    result = {}
    for name, keyword in mapping.items():
        value = getattr(args, name, None)
        result[keyword] = getattr(cfg, f"hostresponse_{name}") if value is None else value
    return result


def _run_hostresponse_subcommand(args: argparse.Namespace) -> None:
    """Run host-response analysis on an existing viralscan output directory."""
    from viralscan.kb_outputs import KbCountOutputs
    from viralscan.runconfig import RunConfig
    from viralscan.scripts.hostresponse import resolve_seeds, run_hostresponse

    configure_logging(verbose=args.verbose, quiet=args.quiet)

    output_dir = Path(args.output).resolve()
    config_yaml = output_dir / "config.yaml"
    if not config_yaml.exists():
        _die(f"config.yaml not found in {output_dir}. Is this a viralscan output directory?")

    analysis_txt = output_dir / "log" / "analysis.txt"
    if not analysis_txt.exists():
        _die(
            f"log/analysis.txt not found in {output_dir}. "
            "The viralscan run must be complete (at least through the analysis step)."
        )

    cfg = RunConfig.from_yaml(str(config_yaml))
    kb = KbCountOutputs.from_config_output(cfg.output)
    virus_h5ad = str(kb.current_adata(multimapping=cfg.multimapping))
    if not Path(virus_h5ad).exists():
        _die(f"Virus h5ad not found at {virus_h5ad}. Ensure the multimap step has completed.")

    from viralscan.defaults import DEFAULTS

    n_seeds = args.n_seeds if args.n_seeds is not None else cfg.hostresponse_n_seeds
    try:
        seeds = resolve_seeds(n_seeds)
    except ValueError as exc:
        _die(str(exc))
    n_stab_iter = (
        args.n_stab_iter
        if args.n_stab_iter is not None
        else (cfg.hostresponse_n_stab_iter or DEFAULTS["hostresponse_n_stab_iter"])
    )
    stab_min_prob = (
        args.stab_min_prob
        if args.stab_min_prob is not None
        else (cfg.hostresponse_stab_min_prob or DEFAULTS["hostresponse_stab_min_prob"])
    )
    top_n_genes = (
        args.top_n_genes
        if args.top_n_genes is not None
        else (cfg.hostresponse_top_n_genes or DEFAULTS["hostresponse_top_n_genes"])
    )
    detection_threshold = (
        args.detection_threshold
        if args.detection_threshold is not None
        else cfg.detection_threshold
    )
    enrichment_db = (
        args.enrichment_db or cfg.hostresponse_enrichment_db or "GO_Biological_Process_2023"
    )

    out_dir = str(output_dir / "hostresponse")
    run_hostresponse(
        virus_h5ad=virus_h5ad,
        host_h5ad=args.host_h5ad,
        viral_accessions_file=str(analysis_txt),
        out_dir=out_dir,
        use_hvg=args.use_hvg,
        seeds=seeds,
        n_stab_iter=n_stab_iter,
        stab_min_prob=stab_min_prob,
        top_n_genes=top_n_genes,
        detection_threshold=detection_threshold,
        do_enrichment=args.enrichment,
        enrichment_db=enrichment_db,
        label=args.label,
        depth_match=args.depth_match,
        control_mito=args.mito_control,
        annotate_symbols=args.gene_symbols,
        differential=args.differential,
        **_hostresponse_evaluation_kwargs(args, cfg),
    )
    restamp_run_complete(output_dir.parent)
    log.info("hostresponse complete. Results in %s", out_dir)


def _run_check_whitelist_subcommand(args: argparse.Namespace) -> None:
    """Run the barcode/whitelist match-rate preflight and exit non-zero on mismatch."""
    from viralscan.whitelist_preflight import check_whitelist

    configure_logging(verbose=args.verbose, quiet=args.quiet)
    if not os.path.exists(args.sample1):
        _die(f"R1 FASTQ not found: {args.sample1}")
    if not os.path.exists(args.whitelist):
        _die(f"Whitelist not found: {args.whitelist}")
    try:
        result = check_whitelist(
            args.sample1,
            args.whitelist,
            args.technology,
            min_match_rate=args.min_match_rate,
            n_sample=args.n_sample,
        )
    except ValueError as exc:
        _die(str(exc))
    if result.ok:
        log.info("%s", result.message)
    else:
        log.error("%s", result.message)
        sys.exit(1)


def _run_check_chemistry_subcommand(args: argparse.Namespace) -> None:
    """Print the library diagnosis; exit 1 when it is unresolved or conflicting."""
    from viralscan import chemistry, chemistry_check

    configure_logging(verbose=args.verbose, quiet=args.quiet)
    if not os.path.exists(args.sample1):
        _die(f"R1 FASTQ not found: {args.sample1}")
    report = chemistry_check.diagnose(
        args.sample1,
        args.sample2,
        whitelist=args.whitelist,
        index=args.index,
        t2g=args.transcripts,
        technology=args.technology,
        cores=args.cores,
    )
    ok = report["library_kind"]["call"] != "unclear" and "conflict" not in report["end"]["basis"]
    if args.technology and report["chemistry"]["chemistry"]:
        try:
            chemistry.resolve(args.technology, [chemistry.Detection(**report["chemistry"])])
        except chemistry.ChemistryError as exc:
            report["advice"].insert(0, f"-x {args.technology} is contradicted: {exc}")
            ok = False
    if report["library_kind"]["call"] == "single-cell" and not report["chemistry"]["chemistry"]:
        ok = False
    print(f"library:   {report['library_kind']['call']} ({report['library_kind']['reason']})")
    print(f"chemistry: {report['chemistry']['chemistry']} ({report['chemistry']['reason']})")
    print(f"end:       {report['end']['call']} ({report['end']['basis']})")
    for line in report["advice"]:
        print(f"advice:    {line}")
    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=2))
    sys.exit(0 if ok else 1)


def create_help() -> argparse.Namespace:
    """Parse the command line. Returns the namespace of all user-given arguments."""
    return build_parser().parse_args()


def _die(message: str, code: int = 1) -> NoReturn:
    log.error(message)
    sys.exit(code)


def _has_valid_fastq_suffix(path: str) -> bool:
    return any(path.endswith(suf) for suf in FASTQ_SUFFIXES)


def errorhandler(args: argparse.Namespace) -> None:
    """
    Validate user input and abort with a clear message if anything is wrong.
    """
    if args.ncbi_accession:
        if args.reference or args.gtf or args.fasta:
            _die("--ncbi-accession is mutually exclusive with --reference / -gtf / -fasta.")
    elif args.reference:
        if args.gtf is None:
            _die("--reference requires -gtf. The code has been terminated.")
        if args.fasta is None:
            _die("--reference requires -fasta. The code has been terminated.")
        for fasta in split_comma_paths(args.fasta):
            if not os.path.exists(fasta):
                _die(f"FASTA path does not exist: {fasta}.")
        for gtf in split_comma_paths(args.gtf):
            if not os.path.exists(gtf):
                _die(f"GTF path does not exist: {gtf}.")
        if args.transcripts is not None or args.index is not None or args.f1 is not None:
            answer = (
                input(
                    "\033[33mYou have provided a path to the index or transcripts file but want to "
                    "create a reference. The reference will be written to the output directory. "
                    "Continue? (yes/y/no/n): \033[0m"
                )
                .strip()
                .lower()
            )
            if answer in ("no", "n"):
                print("You have chosen not to continue. The code has been aborted.")
                sys.exit(0)
            if answer not in ("yes", "y"):
                _die(f"This is not a valid answer: {answer}.")
    else:
        if args.index is None or not os.path.exists(args.index):
            _die(f"Path to index does not exist: {args.index}.")
        if args.transcripts is None or not os.path.exists(args.transcripts):
            _die(f"Path to transcripts does not exist: {args.transcripts}.")

    samples1 = split_comma_paths(args.sample1)
    samples2 = split_comma_paths(args.sample2)
    if len(samples1) != len(samples2):
        _die("--sample1 and --sample2 must have the same number of comma-separated entries.")
    for s1, s2 in zip(samples1, samples2):
        if not os.path.exists(s1) or not os.path.exists(s2):
            _die(f"Sample path does not exist: {s1} or {s2}.")
        if not _has_valid_fastq_suffix(s1):
            _die(f"The forward sample is not in FASTQ format: {s1}.")
        if not _has_valid_fastq_suffix(s2):
            _die(f"The backward sample is not in FASTQ format: {s2}.")

    host_filter = getattr(args, "host_filter", None)
    host_index = getattr(args, "host_index", None)
    if host_index and not host_filter:
        _die("--host-index requires --host-filter.")
    if host_filter and not host_index:
        _die("--host-filter requires --host-index.")
    if getattr(args, "host_filter_star_params", None) and not host_filter:
        _die("--host-filter-star-params requires --host-filter.")
    if host_filter:
        host_index_path = str(host_index)
        if not os.path.exists(host_index_path):
            _die(f"Host index path does not exist: {host_index_path}.")

    if getattr(args, "auto_evidence", None):
        host_fasta = getattr(args, "host_fasta", None)
        viral_fasta = getattr(args, "viral_fasta", None)
        if not host_fasta:
            _die("--auto-evidence requires --host-fasta.")
        if not os.path.isfile(host_fasta):
            _die(f"Host FASTA does not exist: {host_fasta}.")
        if viral_fasta and not os.path.isfile(viral_fasta):
            _die(f"Viral FASTA does not exist: {viral_fasta}.")
        # Only <index dir>/viral.fa is viral-only; a freshly built reference has none.
        built_here = args.reference or args.ncbi_accession
        beside_index = not built_here and os.path.isfile(
            os.path.join(os.path.dirname(os.path.abspath(args.index)), "viral.fa")
        )
        if not viral_fasta and not beside_index:
            _die("--auto-evidence requires --viral-fasta (no viral.fa next to the index).")
        if not args.multimapping:
            _die("--auto-evidence needs multimapping; drop --no-multimapping.")
        missing = [t for t in ("minimap2", "samtools") if shutil.which(t) is None]
        if missing:
            _die(f"--auto-evidence needs {', '.join(missing)} on PATH.")

    log.info("All input data has been checked and is correct.")


def _resolve_chemistry(args: argparse.Namespace) -> dict[str, dict[str, Any]]:
    """Detect each sample's chemistry and set ``args.technology`` (WP1E Q6, F-005).

    Returns the per-sample detection blocks for ``run_manifest.json``. Fails the
    run on ambiguity or a mismatch with ``-x`` unless ``--force-technology``.
    """
    from viralscan import chemistry

    samples = split_comma_paths(args.sample1)
    log.info("Checking the chemistry of %d sample(s) from their R1 reads...", len(samples))
    try:
        detections = chemistry.detect(samples, args.whitelist)
    except chemistry.ChemistryError as exc:
        if not args.force_technology:
            _die(str(exc))
        log.warning("Chemistry check skipped: %s", exc)
        detections = []
    try:
        args.technology = chemistry.resolve(args.technology, detections, args.force_technology)
    except chemistry.ChemistryError as exc:
        _die(str(exc))
    for s1, d in zip(samples, detections):
        level = logging.INFO if d.chemistry == args.technology else logging.WARNING
        log.log(level, "Chemistry of %s: %s (%s)", _sample_id(s1), d.chemistry, d.reason)
    log.info("Running with -x %s", args.technology)
    return {_sample_id(s1): d.as_block() for s1, d in zip(samples, detections)}


def _check_required_tools() -> None:
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        _die(
            "The following required external tools are not on PATH: "
            f"{', '.join(missing)}. The full workflow needs kb-python (provides 'kb') and "
            "snakemake: install the 'full' extra (pip install \"viralscan[full]\") or use "
            "environment.yml (Python 3.11, linux-64)."
        )


def _resolve_cell_calling(method: str, called_cells_file: Optional[str]) -> str:
    """Resolve ``auto`` the same way :func:`viralscan.scripts.cellcalling.call_cells` does."""
    if method == "auto":
        return "external" if called_cells_file else "emptydrops"
    return method


def _check_cell_caller_tools(method: str, rscript: str) -> None:
    """Verify the resolved cell caller's tools are on PATH.

    Cell calling now fails closed, and it runs after kb_count, analysis, and
    multimap. Without this check a missing R turns a several-hour run into a
    late abort for a reason that was knowable before it started.
    """
    if method != "emptydrops":
        return
    if shutil.which(rscript) is None:
        _die(
            f"--cell-calling {method} requires '{rscript}' on PATH (it runs "
            "DropletUtils::emptyDrops). Install R plus DropletUtils, point at "
            "another interpreter, or choose a different --cell-calling method."
        )


def _check_host_filter_tools(aligner: str) -> None:
    """Verify that tools required by the chosen host-filter aligner are on PATH."""
    if aligner == "starsolo":
        if shutil.which("STAR") is None:
            _die(
                "--host-filter starsolo requires 'STAR' on PATH. "
                "Install STARsolo: https://github.com/alexdobin/STAR"
            )
    else:
        _die("ViralScan v3 supports only --host-filter starsolo.")


def _count_lines(path: str) -> int:
    with open(path, encoding="utf-8", errors="replace") as f:
        return sum(1 for _ in f)


def _count_unique_genes(t2g_path: str) -> int:
    """Count unique (col2, col3) pairs in a t2g.txt file."""
    seen: set[tuple[str, str]] = set()
    with open(t2g_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                seen.add((cols[1], cols[2]))
    return len(seen)


def _materialize_reference_input(paths: list[str], destination: Path) -> str:
    """Concatenate reference input files with one newline between files."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open(destination, "wb") as out:
        for path in paths:
            data = Path(path).read_bytes().rstrip(b"\r\n")
            out.write(data)
            if data:
                out.write(b"\n")
    return str(destination)


def _prepare_kb_ref_inputs(output_dir: Path, fasta_arg: str, gtf_arg: str) -> tuple[str, str]:
    """Return FASTA/GTF paths safe to pass to ``kb ref``."""
    fasta_paths = split_comma_paths(fasta_arg)
    gtf_paths = split_comma_paths(gtf_arg)
    if not fasta_paths:
        raise ValueError("At least one FASTA path is required to build a reference.")
    if not gtf_paths:
        raise ValueError("At least one GTF path is required to build a reference.")

    index_dir = output_dir / "index"
    fasta = (
        fasta_paths[0]
        if len(fasta_paths) == 1
        else _materialize_reference_input(fasta_paths, index_dir / "input.fasta")
    )
    gtf = (
        gtf_paths[0]
        if len(gtf_paths) == 1
        else _materialize_reference_input(gtf_paths, index_dir / "input.gtf")
    )
    return fasta, gtf


def _build_run_config(
    args: argparse.Namespace,
    outs: str,
    index: str,
    transcripts: str,
    f1: Optional[str],
    s1: str,
    s2: str,
) -> RunConfig:
    """Build the validated :class:`RunConfig` for one sample.

    Constructs a :class:`~viralscan.runconfig.RunConfig` via
    :meth:`~viralscan.runconfig.RunConfig.from_snakemake_config` (the single
    validation checkpoint); :func:`_write_run_config` writes it as ``config.yaml``, which
    Snakemake reads. There is no hand-maintained key list, so CLI flags like
    ``--multimap-em-max-iter`` cannot be omitted from the Snakemake invocation.
    """
    return RunConfig.from_snakemake_config(
        {
            "output": outs,
            "index": index,
            "transcripts": transcripts,
            "sample1": s1,
            "sample2": s2,
            "cores": args.cores,
            "gtf": args.gtf,
            "fasta": args.fasta,
            "visual": args.visual,
            "f1": f1,
            "reference": args.reference,
            "umap": args.umap,
            "technology": args.technology,
            "require_chemistry_sanity": getattr(args, "require_chemistry_sanity", False),
            "whitelist": args.whitelist,
            "strand": getattr(args, "strand", None),
            "multimapping": args.multimapping,
            "se_threshold": args.se_threshold,
            "detection_threshold": args.detection_threshold,
            "positive_control_gene": getattr(args, "positive_control_gene", None),
            "positive_control_receipt": getattr(args, "positive_control_receipt", None),
            "positive_control_expected_molecules": getattr(
                args, "positive_control_molecules", None
            ),
            "positive_control_scope": getattr(args, "positive_control_scope", None),
            "positive_control_virus_key": getattr(args, "positive_control_virus_key", None),
            "require_positive_control": getattr(
                args, "require_positive_control", DEFAULTS["require_positive_control"]
            ),
            "gene_programs": getattr(args, "gene_programs", DEFAULTS["gene_programs"]),
            "anello_align": getattr(args, "anello_align", DEFAULTS["anello_align"]),
            "auto_evidence": getattr(args, "auto_evidence", None),
            "host_fasta": getattr(args, "host_fasta", None),
            "viral_fasta": getattr(args, "viral_fasta", None),
            "programme_min_breadth": getattr(
                args, "programme_min_breadth", DEFAULTS["programme_min_breadth"]
            ),
            "programme_latent_min_breadth": getattr(
                args, "programme_latent_min_breadth", DEFAULTS["programme_latent_min_breadth"]
            ),
            "programme_min_umi": getattr(args, "programme_min_umi", DEFAULTS["programme_min_umi"]),
            # getattr matches the existing convention for args that older
            # callers may not set (see the host_filter handling below).
            "anellovirus_gene_ids": getattr(
                args, "anellovirus_gene_ids", DEFAULTS["anellovirus_gene_ids"]
            ),
            "min_counts": args.min_counts,
            "min_genes": args.min_genes,
            "hvg_min_mean": args.hvg_min_mean,
            "hvg_max_mean": args.hvg_max_mean,
            "hvg_min_disp": args.hvg_min_disp,
            "umap_n_neighbors": args.umap_n_neighbors,
            "multimap_method": args.multimap_method,
            "multimap_pseudocount": args.multimap_pseudocount,
            "multimap_primary_call": args.multimap_primary_call,
            "multimap_em_max_iter": args.multimap_em_max_iter,
            "multimap_em_tol": args.multimap_em_tol,
            "multimap_molecule_assignments": getattr(
                args, "multimap_molecule_assignments", DEFAULTS["multimap_molecule_assignments"]
            ),
            "cell_types": args.cell_types,
            "data_cache_dir": args.data_cache_dir,
            "host_filter_aligner": getattr(args, "host_filter", None),
            "host_filter_star_params": getattr(args, "host_filter_star_params", None),
            "read_filter": getattr(args, "read_filter", None),
            "host_index": getattr(args, "host_index", None),
            "host_h5ad": getattr(args, "host_h5ad", None),
            "hostresponse_cv": getattr(args, "hostresponse_cv", None),
            "hostresponse_groups": getattr(args, "hostresponse_groups", None),
            "hostresponse_cv_folds": getattr(args, "hostresponse_cv_folds", None),
            "hostresponse_cell_type_column": getattr(args, "hostresponse_cell_type_column", None),
            "hostresponse_cell_types": getattr(args, "hostresponse_cell_types", None),
            "hostresponse_panel_in_fold": getattr(args, "hostresponse_panel_in_fold", None),
            "hostresponse_permutations": getattr(args, "hostresponse_permutations", None),
            "hostresponse_permutation_unit": getattr(args, "hostresponse_permutation_unit", None),
            "hostresponse_permutation_block": getattr(args, "hostresponse_permutation_block", None),
            "hostresponse_min_negative_cells": getattr(
                args, "hostresponse_min_negative_cells", None
            ),
            "hostresponse_min_groups": getattr(args, "hostresponse_min_groups", None),
            "hostresponse_n_seeds": getattr(args, "hostresponse_n_seeds", None),
            "hostresponse_n_stab_iter": getattr(args, "hostresponse_n_stab_iter", None),
            "hostresponse_use_hvg": getattr(args, "hostresponse_use_hvg", True),
            "hostresponse_stab_min_prob": getattr(args, "hostresponse_stab_min_prob", None),
            "hostresponse_top_n_genes": getattr(args, "hostresponse_top_n_genes", None),
            "hostresponse_enrichment": getattr(args, "enrichment", False),
            "hostresponse_enrichment_db": getattr(args, "enrichment_db", None),
            "hostresponse_label": getattr(args, "hostresponse_label", None),
            "hostresponse_depth_match": getattr(args, "hostresponse_depth_match", False),
            "hostresponse_control_mito": getattr(args, "hostresponse_control_mito", True),
            "hostresponse_differential": getattr(args, "hostresponse_differential", False),
            "cell_calling": getattr(args, "cell_calling", None),
            "called_cells_file": getattr(args, "called_cells_file", None),
            "emptydrops_seed": getattr(args, "emptydrops_seed", None),
            "emptydrops_niters": getattr(args, "emptydrops_niters", None),
        }
    )


def _write_sample_summary(
    outs: str,
    elapsed: float,
    n_transcripts: int,
    n_genes: int,
) -> None:
    """Append per-sample runtime and reference stats to the sample's summary.txt."""
    summary_path = os.path.join(outs, "summary.txt")
    os.makedirs(outs, exist_ok=True)
    with open(summary_path, "a", encoding="utf-8") as f:
        f.write(f"\nRuntime: {elapsed:.4f} seconds.\n\n")
        f.write(f"Amount of transcripts in data: {n_transcripts}\n")
        f.write(f"Amount of genes in data: {n_genes}\n")


def _sample_id(s1_path: str) -> str:
    """Derive a per-sample output-directory name from the forward FASTQ path.

    Convention: the stem before the first ``_`` in the filename.
    Example: ``/data/SRR123_R1.fastq.gz`` → ``SRR123``.
    """
    return Path(s1_path).name.split("_")[0]


def _write_reference_manifest(index: str, t2g: str, fasta: str, gtf: str) -> None:
    """Record the index's host/viral gene sets next to it (PLAN ``DEF-03``)."""
    from viralscan.run_safety import sha256_file, software_identity
    from viralscan.virus_identity import gtf_gene_ids, write_build_manifest_from_t2g

    path = write_build_manifest_from_t2g(
        index,
        t2g,
        gtf_gene_ids(gtf),
        provenance={
            "builder": "viralscan --reference",
            "software_identity": software_identity(),
            "fasta": {"path": str(Path(fasta).resolve()), "sha256": sha256_file(Path(fasta))},
            "gtf": {"path": str(Path(gtf).resolve()), "sha256": sha256_file(Path(gtf))},
        },
    )
    log.info("Index build manifest: %s", path)


def _build_kb_ref(output_dir: Path, fasta: str, gtf: str) -> tuple[str, str, str]:
    """Run ``kb ref`` to build an index. Returns (transcripts, index, f1) paths."""
    from viralscan.scripts.build_reference import (
        _run_kb_ref,
        prepare_reference_inputs,
        record_reference_build,
    )
    from viralscan.virus_identity import gtf_gene_ids

    index_dir = output_dir / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    fasta_input, gtf_input = _prepare_kb_ref_inputs(output_dir, fasta, gtf)
    prepared, sequence_manifest = prepare_reference_inputs(
        Path(fasta_input),
        Path(gtf_input),
        index_dir,
        viral_gene_ids=gtf_gene_ids(gtf_input),
        mask=True,
    )
    transcripts = str(index_dir / "t2g.txt")
    index = str(index_dir / "index.idx")
    f1 = str(index_dir / "cdna.fa")
    log.info("Building kb ref index. Depending on the genome this can take a while...")
    command = [
        "kb",
        "ref",
        "-i",
        index,
        "-g",
        transcripts,
        "-f1",
        f1,
        "--overwrite",
        str(prepared),
        gtf_input,
    ]
    resources = _run_kb_ref(command, index_dir)
    record_reference_build(
        sequence_manifest,
        fasta=prepared,
        gtf=Path(gtf_input),
        t2g=Path(transcripts),
        command=command,
        resources=resources,
        index=Path(index),
    )
    log.info("Reference index is done!")
    _write_reference_manifest(index, transcripts, str(prepared), gtf_input)
    return transcripts, index, f1


def _resolve_auto_strand(
    args: argparse.Namespace, output_dir: Path, sample: str, s1: str, s2: str, index: str, t2g: str
) -> str:
    """Strand for ``--strand auto``: reuse the recorded choice on resume, else pilot."""
    from viralscan import strand as _strand

    block = recorded_strand_inference(output_dir, sample) if args.resume else None
    if block is None:
        log.info("Piloting strandedness for %s (first %d pairs)...", sample, _strand.PILOT_READS)
        rates = _strand.run_pilot(s1, s2, index, t2g, args.technology, args.whitelist, args.cores)
        block = _strand.inference_block(rates)
        record_strand_inference(output_dir, sample, block)
    choice = str(block["choice"])
    log.info("Strand for %s: %s (pilot rates %s)", sample, choice, block["rates"])
    return choice


def main() -> None:
    args = create_help()

    # Dispatch to build-ref subcommand if requested.
    if getattr(args, "_subcommand", None) == "data-fetch":
        configure_logging(verbose=args.verbose, quiet=args.quiet)
        from viralscan.data_fetch import ViralScanDataError, fetch_viral_data

        try:
            data_dir = fetch_viral_data(
                cache_dir=args.cache_dir,
                archive_url=args.url,
                expected_sha256=args.sha256,
                force=args.force,
            )
        except ViralScanDataError as exc:
            _die(str(exc))
        log.info("Viral annotation data is available at %s", data_dir)
        return

    if getattr(args, "_subcommand", None) == "data":
        _die("Missing data command. Use `viralscan data fetch`.")

    if getattr(args, "_subcommand", None) == "build-ref":
        from viralscan.scripts.build_reference import build_ref_main

        build_ref_main(args)
        return

    if getattr(args, "_subcommand", None) == "evidence":
        from viralscan.scripts.evidence_run import run_evidence

        run_evidence(args)
        return

    if getattr(args, "_subcommand", None) == "rerun-programs":
        _run_rerun_programs(args)
        return

    if getattr(args, "_subcommand", None) == "rerun-multimap":
        _run_rerun_multimap(args)
        return

    if getattr(args, "_subcommand", None) == "doctor":
        import json as _json

        from viralscan.validation import doctor_report

        report = doctor_report(args.profile)
        print(_json.dumps(report, indent=2, sort_keys=True))
        if not report["ok"]:
            sys.exit(1)
        return

    if getattr(args, "_subcommand", None) == "validate-run":
        import json as _json

        from viralscan.validation import validate_run

        report = validate_run(Path(args.run_dir), verify_inputs=not args.no_verify_inputs)
        rendered = _json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.json_output:
            Path(args.json_output).write_text(rendered, encoding="utf-8")
        print(rendered, end="")
        if not report["ok"]:
            sys.exit(1)
        return

    if getattr(args, "_subcommand", None) == "hostresponse":
        _run_hostresponse_subcommand(args)
        return

    if getattr(args, "_subcommand", None) == "check-whitelist":
        _run_check_whitelist_subcommand(args)
        return

    if getattr(args, "_subcommand", None) == "check-chemistry":
        _run_check_chemistry_subcommand(args)
        return

    configure_logging(verbose=args.verbose, quiet=args.quiet)

    # Validate that required run-mode args are present (they are optional in
    # the argparse definition to allow build-ref to coexist).
    if args.output is None:
        _die("--output / -o is required for viral quantification.")
    if args.sample1 is None:
        _die("--sample1 / -s1 is required for viral quantification.")
    if args.sample2 is None:
        _die("--sample2 / -s2 is required for viral quantification.")

    _check_required_tools()
    errorhandler(args)

    if args.host_filter:
        _check_host_filter_tools(args.host_filter)

    _check_cell_caller_tools(
        _resolve_cell_calling(
            getattr(args, "cell_calling", DEFAULTS["cell_calling"]),
            getattr(args, "called_cells_file", None),
        ),
        getattr(args, "cell_caller_rscript", DEFAULTS["cell_caller_rscript"]),
    )

    from viralscan.run_safety import (
        RunSafetyError,
        build_run_manifest,
        prepare_output_directory,
        record_manifest_block,
    )

    chemistry_blocks = _resolve_chemistry(args)
    output_dir = Path(args.output).resolve()
    try:
        prepare_output_directory(
            output_dir,
            build_run_manifest(args),
            resume=args.resume,
            overwrite=args.overwrite,
            yes=args.yes,
        )
    except RunSafetyError as exc:
        _die(str(exc))
    for sample, block in chemistry_blocks.items():
        record_manifest_block(output_dir, "chemistry_detection", sample, block)
    # REL-16: same-version kallisto/bustools binaries differ, so record which ran.
    from viralscan.validation import tool_provenance

    record_manifest_block(output_dir, "tool_binaries", "run", tool_provenance())

    if args.ncbi_accession:
        from viralscan.scripts.ncbi_fetch import NCBIFetchError, fetch_reference

        accessions = split_comma_paths(args.ncbi_accession)
        ref_dir = output_dir / "ncbi_reference"
        try:
            fasta_path, gtf_path = fetch_reference(
                accessions=accessions,
                out_dir=ref_dir,
                email=args.ncbi_email,
                api_key=args.ncbi_api_key,
            )
        except NCBIFetchError as exc:
            _die(f"NCBI download failed: {exc}")
        args.fasta = str(fasta_path)
        args.gtf = str(gtf_path)
        args.reference = True
        log.info("Fetched NCBI reference for: %s", ", ".join(accessions))

    if args.reference:
        transcripts, index, f1 = _build_kb_ref(output_dir, args.fasta, args.gtf)
    else:
        transcripts = args.transcripts
        index = args.index
        f1 = args.f1

    snakefile_path = os.path.join(os.path.dirname(__file__), "Snakefile")
    samples1 = split_comma_paths(args.sample1)
    samples2 = split_comma_paths(args.sample2)
    output = str(output_dir)

    # Fail fast if two inputs share the same derived sample ID before running anything.
    seen_ids: set[str] = set()
    for s1 in samples1:
        sid = _sample_id(s1)
        if sid in seen_ids:
            _die(
                f"Duplicate derived sample ID '{sid}'. Two --sample1 paths share the same "
                f"filename prefix before the first '_'. Rename your input files or supply "
                f"unique prefixes so each sample gets its own output directory."
            )
        seen_ids.add(sid)

    # Hoist loop-invariant transcript counts (they depend only on the reference).
    n_transcripts = _count_lines(transcripts)
    n_genes = _count_unique_genes(transcripts)

    for s1, s2 in zip(samples1, samples2):
        sample_start = time.time()
        out = _sample_id(s1)
        outs = os.path.join(output, out) + os.sep
        sample_args = args
        if args.strand == "auto":
            chosen = _resolve_auto_strand(args, output_dir, out, s1, s2, index, transcripts)
            sample_args = argparse.Namespace(**{**vars(args), "strand": chosen})
        run_config = _build_run_config(sample_args, outs, index, transcripts, f1, s1, s2)
        config_file = _write_run_config(run_config)
        if args.resume:
            _backfill_identity_table(run_config)
        subprocess.run(_snakemake_run_command(snakefile_path, args.cores, config_file), check=True)

        _write_sample_summary(outs, time.time() - sample_start, n_transcripts, n_genes)

        subprocess.run(_snakemake_unlock_command(snakefile_path, config_file), check=True)

    # Every sample finished (check=True above): the run is complete.
    write_run_complete(output_dir)


if __name__ == "__main__":
    main()

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

from viralscan.defaults import (
    CELL_CALLING_METHODS,
    DEFAULTS,
    MULTIMAP_METHODS,
    MULTIMAP_PRIMARY_CALLS,
)
from viralscan.run_safety import (
    RUN_MANIFEST,
    clear_run_complete,
    record_strand_inference,
    recorded_strand_inference,
    restamp_run_complete,
    write_run_complete,
)
from viralscan.runconfig import CAPTURE_SCOPES, RunConfig
from viralscan.utils import configure_logging, split_comma_paths

try:
    from pyfiglet import figlet_format as _figlet_format
except ImportError:  # pyfiglet is optional

    def _figlet_format(text: str, font: str = "standard", **kwargs: Any) -> Any:
        return text


figlet_format = _figlet_format

log = logging.getLogger("viralscan")


REQUIRED_TOOLS = ("kb", "snakemake")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def _add_verbosity_args(parser: Any) -> None:
    """Add the --verbose/--quiet pair shared by every subcommand."""
    parser.add_argument(
        "--verbose", action="store_true", default=False, help="Enable DEBUG-level logging."
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress INFO messages."
    )


def _build_data_parser(subparsers: Any) -> None:
    """Register the 'data' subcommand group."""
    data = subparsers.add_parser(
        "data",
        help="Manage ViralScan reference data downloads.",
        description="Manage ViralScan's external viral annotation data.",
    )
    data.set_defaults(_subcommand="data")
    data_subparsers = data.add_subparsers(dest="_data_command")
    fetch = data_subparsers.add_parser(
        "fetch",
        help="Download bundled viral GTF annotations from Zenodo.",
        description=(
            "Download the ViralScan viral annotation panel from Zenodo, verify the "
            "archive checksum, and unpack GTF files into ~/.cache/viralscan/data/."
        ),
    )
    fetch.add_argument(
        "--cache-dir",
        default=None,
        help="Root cache directory. Default: ~/.cache/viralscan/",
    )
    fetch.add_argument(
        "--url",
        default=None,
        help="Override archive URL. Intended for tests or mirrors; defaults to Zenodo.",
    )
    fetch.add_argument(
        "--sha256",
        default=None,
        help="Optional expected SHA-256 digest for the downloaded archive.",
    )
    fetch.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Re-download and replace cached GTF files even if data already exists.",
    )
    _add_verbosity_args(fetch)
    fetch.set_defaults(_subcommand="data-fetch")


def _build_ref_parser(subparsers: Any) -> None:
    """Register the 'build-ref' subcommand."""
    p = subparsers.add_parser(
        "build-ref",
        help="Build a combined host + virus kallisto reference from Ensembl + NCBI.",
        description=(
            "Download a host cDNA FASTA + GTF from Ensembl and viral sequences from NCBI, "
            "concatenate them, and optionally run 'kb ref' to build a kallisto index. "
            "This is the recommended host-aware reference setup for ViralScan.\n\n"
            "Example:\n"
            "  viralscan build-ref --host human \\\n"
            "      --virus-accessions NC_045512.2 NC_002021.1 \\\n"
            "      --output ref_human_covid/ --ncbi-email you@example.org\n\n"
            "  viralscan build-ref --list-species"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--host",
        default=None,
        help="Host species, e.g. 'human', 'mouse'. Run --list-species for all options.",
    )
    p.add_argument(
        "--virus-accessions",
        nargs="+",
        metavar="ACCESSION",
        default=None,
        help="One or more NCBI nucleotide accessions, e.g. NC_045512.2.",
    )
    p.add_argument(
        "--profile",
        choices=["curated", "broad-discovery"],
        default="curated",
        help="Frozen combined-reference profile. Default: curated.",
    )
    p.add_argument(
        "--output",
        "-o",
        default="viralscan_ref",
        help="Output directory for reference files. Default: viralscan_ref/",
    )
    p.add_argument(
        "--ncbi-email",
        default=None,
        help="Contact e-mail for NCBI E-utilities (avoids throttling).",
    )
    p.add_argument(
        "--ncbi-api-key",
        default=None,
        help="NCBI API key for higher request rates.",
    )
    p.add_argument(
        "--cache-dir",
        default=None,
        help="Root directory for download cache. Default: ~/.cache/viralscan/",
    )
    p.add_argument(
        "--no-kb-ref",
        action="store_true",
        default=False,
        help="Skip running 'kb ref'; only produce concatenated FASTA and GTF.",
    )
    p.add_argument(
        "--genome-dlist",
        default=None,
        metavar="FASTA",
        help=(
            "Full host-genome FASTA passed to kallisto as a D-list and used for raw "
            "viral host-homology annotation. Requires minimap2 and the full tool profile."
        ),
    )
    p.add_argument(
        "--anellovirus",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Include the full packaged Anelloviridae accession table (2,041 accessions) "
            "in the combined host+viral reference (explicit opt-in; default: off). "
            "When --reference-panel anellovirus is used instead, builds an "
            "Anelloviridae-only reference without a host transcriptome; combine with "
            "--no-mask / --cluster for masking/clustering options."
        ),
    )
    p.add_argument(
        "--allow-partial-panel",
        action="store_true",
        default=False,
        help="Allow an incomplete expanded panel and write missing_accessions.tsv; default fails closed.",
    )
    p.add_argument(
        "--no-mask",
        action="store_true",
        default=False,
        help="Skip dustmasker hard-masking of low-complexity viral sequence; the k-mer gate then allows a few.",
    )
    p.add_argument(
        "--cluster",
        action="store_true",
        default=False,
        help="(--anellovirus) Run cd-hit-est clustering at 95%% identity after masking.",
    )
    p.add_argument(
        "--reference-panel",
        choices=["anellovirus"],
        default=None,
        metavar="PANEL",
        help=(
            "Build a pre-defined reference panel. Currently supported: 'anellovirus'. "
            "Uses the bundled FASTA from `viralscan data fetch` when available, "
            "otherwise falls back to NCBI accession download (same as --anellovirus). "
            "Combine with --no-mask / --cluster for masking/clustering options."
        ),
    )
    p.add_argument(
        "--list-species",
        action="store_true",
        default=False,
        help="Print all supported host species and exit.",
    )
    _add_verbosity_args(p)
    p.set_defaults(_subcommand="build-ref")


def _build_evidence_parser(subparsers: Any) -> None:
    """Register the 'evidence' subcommand."""
    p = subparsers.add_parser(
        "evidence",
        help="Extract, visualize (IGV) and score the reads behind viral calls.",
        description=(
            "Trace reads supporting selected viral CB-UMI molecules in a COMPLETED ViralScan "
            "run, extract them, competitively align against host plus target for IGV, and score "
            "coverage and BLAST evidence. These outputs are diagnostic and do not by themselves "
            "confirm infection.\n\n"
            "Example:\n"
            "  viralscan evidence --run-dir output/sample/ --viral-fasta viruses.fa \\\n"
            "      --virus EBV --blast -o output/sample/evidence/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--run-dir", required=True, help="A completed ViralScan run output directory.")
    p.add_argument("--output", "-o", required=True, help="Directory for evidence outputs.")
    p.add_argument(
        "--viral-fasta",
        default=None,
        help="Viral genome FASTA to align extracted reads against (enables BAM/coverage/BLAST). "
        "Omit for read-extraction only.",
    )
    p.add_argument(
        "--host-fasta",
        default=None,
        help="Full host-genome FASTA required with --viral-fasta for competitive v3 "
        "alignment and BLAST.",
    )
    p.add_argument(
        "--virus",
        required=True,
        help="Required exact target: accession/gene ID, canonical detected virus label, or "
        "registered alias (for example EBV or HHV6B). Substring matching is not used.",
    )
    p.add_argument(
        "--blast",
        action="store_true",
        default=False,
        help="BLAST a sample of extracted reads against the viral reference (requires blast+).",
    )
    p.add_argument(
        "--read-start-profile",
        action="store_true",
        default=False,
        help="Write a per-position 5' read-start distribution along the viral genome "
        "(read_start_profile.tsv). Requires --viral-fasta.",
    )
    p.add_argument(
        "--cell-tags",
        action="store_true",
        default=False,
        help="Write viral_reads.tagged.bam with CB/UB cell-barcode tags (from read names) "
        "for per-cell IGV inspection (group by tag CB). Requires --viral-fasta.",
    )
    p.add_argument(
        "--dedup",
        choices=("umi", "markdup", "none"),
        default="umi",
        help="PCR-duplicate handling for --read-start-profile: umi (collapse per CB+UMI; "
        "default), markdup (samtools markdup), or none. Default: umi.",
    )
    p.add_argument(
        "--bin-size",
        type=int,
        default=1,
        help="Bin width (bp) for the read-start profile. Default: 1.",
    )
    p.add_argument(
        "--sampling-seed",
        type=int,
        default=42,
        help="Seed for deterministic BLAST read sampling. Default: 42.",
    )
    p.add_argument(
        "--cores", "-c", type=int, default=4, help="Threads for minimap2/samtools/blast."
    )
    _add_verbosity_args(p)
    p.set_defaults(_subcommand="evidence")


def _build_rerun_multimap_parser(subparsers: Any) -> None:
    """Register the 'rerun-multimap' subcommand."""
    p = subparsers.add_parser(
        "rerun-multimap",
        help="Switch multimapping method on a completed run without redoing kb count.",
        description=(
            "Re-run the multimapping step with a different algorithm for an existing run,\n"
            "skipping the expensive pseudoalignment (kb count).\n\n"
            "For equal / host-conservative / unique-weighted: all three layers are\n"
            "pre-stored in every multimap h5ad, so the switch is instant (no bus-file\n"
            "reprocessing). For EM methods, the BUS file is reprocessed (slower, still skips\n"
            "kb count). Detection and UMAP are regenerated in the new result directory; "
            "host-response analyses must be rerun separately.\n\n"
            "Tip: the default (host-conservative) is the safe choice; rerun with\n"
            "--multimap-method equal for an equal-allocation pass or em-global for refinement.\n\n"
            "Examples:\n"
            "  viralscan rerun-multimap --run-dir out/ -o out_hc/ "
            "--multimap-method host-conservative\n"
            "  viralscan rerun-multimap --run-dir out/ -o out_em/ "
            "--multimap-method em-global --cores 8"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--run-dir",
        required=True,
        help="Completed source run. It is never modified.",
    )
    p.add_argument(
        "--output",
        "-o",
        required=True,
        help="New result directory; must not already be non-empty.",
    )
    p.add_argument(
        "--multimap-method",
        required=True,
        choices=MULTIMAP_METHODS,
        help="Multimapping resolution method to apply.",
    )
    p.add_argument(
        "--cores",
        "-c",
        default=6,
        type=int,
        help="Number of cores for snakemake workers. Default: 6.",
    )
    p.add_argument(
        "--multimap-em-max-iter",
        type=int,
        default=None,
        metavar="N",
        help="(em only) Maximum EM iterations. Default: preserves existing config value.",
    )
    p.add_argument(
        "--multimap-em-tol",
        type=float,
        default=None,
        metavar="TOL",
        help="(em only) EM convergence tolerance. Default: preserves existing config value.",
    )
    _add_verbosity_args(p)
    p.set_defaults(_subcommand="rerun-multimap")


def _build_rerun_programs_parser(subparsers: Any) -> None:
    """Register the 'rerun-programs' subcommand."""
    p = subparsers.add_parser(
        "rerun-programs",
        help="Infer viral gene programmes on a completed run, in place.",
        description=(
            "Run layer 2 (latent vs productive gene programme) over an existing run,\n"
            "writing results/gene_program_summary.tsv and results/gene_program_cells.tsv.\n\n"
            "Unlike rerun-multimap this operates IN PLACE: layer 2 only reads the count\n"
            "matrices and layer 1's viral_summary.tsv, and adds two files. It changes no\n"
            "counts, so there is no reason to copy the run and no risk of the two\n"
            "directories disagreeing.\n\n"
            "Scope is the viruses layer 1 detected. A virus layer 1 did not call gets\n"
            "no programme row -- that is layer 1's detection-limit problem, not one this\n"
            "command can repair."
        ),
    )
    p.add_argument(
        "--run-dir",
        required=True,
        metavar="DIR",
        help="A completed viralscan run directory (the one holding kb-python/ and results/).",
    )
    p.add_argument(
        "--programme-min-breadth",
        type=int,
        default=DEFAULTS["programme_min_breadth"],
        metavar="N",
        help=(
            "Distinct non-overlapping overlap groups required before a programme is "
            f"called. Must be >= 1. Default: {DEFAULTS['programme_min_breadth']}."
        ),
    )
    _add_verbosity_args(p)
    p.set_defaults(_subcommand="rerun-programs")


def _build_doctor_parser(subparsers: Any) -> None:
    p = subparsers.add_parser("doctor", help="Check Python and full-workflow dependencies.")
    p.add_argument("--profile", choices=("pip", "full"), default="full")
    p.add_argument("--json", action="store_true", default=False)
    p.set_defaults(_subcommand="doctor")


def _build_validate_run_parser(subparsers: Any) -> None:
    p = subparsers.add_parser("validate-run", help="Validate v3 schemas and count invariants.")
    p.add_argument("run_dir", help="ViralScan output directory containing run_manifest.json.")
    p.add_argument("--no-verify-inputs", action="store_true", default=False)
    p.add_argument("--json-output", default=None, metavar="PATH")
    p.set_defaults(_subcommand="validate-run")


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
    }[new_method]
    # counts_host_viral_selected is method-dependent (0 for host-conservative,
    # equal/weighted shares otherwise), so it must be swapped from the matching
    # per-method layer too; leaving it behind mislabels the evidence tiers.
    host_viral_layer_name = {
        "equal": "counts_host_viral_selected_equal",
        "host-conservative": "counts_host_viral_selected_host_conservative",
        "unique-weighted": "counts_host_viral_selected_unique_weighted",
    }[new_method]
    adata = _ad.read_h5ad(str(adata_path))
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
    if args.programme_min_breadth < 1:
        _die("--programme-min-breadth must be >= 1.")

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


def _snakemake_run_command(snakefile_path: str, cores: int, config_args: list[str]) -> list[str]:
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
        "--config",
        *config_args,
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
        stat = config_yaml_path.stat()
        run_config.to_yaml(config_yaml_path)
        os.utime(config_yaml_path, (stat.st_atime, stat.st_mtime))
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

        # The one RunConfig -> `--config` serialisation, not a private copy.
        config_args = run_config.to_snakemake_config_args()

        subprocess.run(_snakemake_run_command(snakefile_path, args.cores, config_args), check=True)

        unlock_cmd = [
            "snakemake",
            "--snakefile",
            snakefile_path,
            "--unlock",
            "--config",
            *config_args,
        ]
        subprocess.run(unlock_cmd, check=True)

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
    ``--output``, and ``createconfig`` then creates a subdirectory per sample
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


def _build_hostresponse_parser(subparsers: Any) -> None:
    """Register the 'hostresponse' subcommand."""
    p = subparsers.add_parser(
        "hostresponse",
        help="Run host-response analysis on an existing viralscan output directory.",
        description=(
            "Associate viral presence with host gene expression via logistic regression\n"
            "(Luebbert et al. 2026 approach) on a completed viralscan run.\n\n"
            "The viralscan output directory must already contain a config.yaml and a\n"
            "completed multimap step (log/multimap.done).  Results are written to\n"
            "<output>/hostresponse/.\n\n"
            "Examples:\n"
            "  viralscan hostresponse -o out/sample/ --host-h5ad host.h5ad\n"
            "  viralscan hostresponse -o out/sample/ --host-h5ad host.h5ad \\\n"
            "    --n-stab-iter 200 --enrichment"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--output",
        "-o",
        required=True,
        help="Existing viralscan sample output directory (contains config.yaml).",
    )
    p.add_argument(
        "--host-h5ad",
        required=True,
        metavar="PATH",
        help="Host gene-expression h5ad (cells × genes, matched to the viralscan run).",
    )
    p.add_argument(
        "--n-seeds",
        type=int,
        default=None,
        metavar="N",
        help="Random seeds for multi-seed L2 regression (default: from config or 6).",
    )
    p.add_argument(
        "--n-stab-iter",
        type=int,
        default=None,
        metavar="N",
        help="Stability-selection iterations (default: from config or 100).",
    )
    p.add_argument(
        "--no-use-hvg",
        dest="use_hvg",
        action="store_false",
        default=True,
        help="Use all genes instead of highly variable genes as features.",
    )
    p.add_argument(
        "--stab-min-prob",
        type=float,
        default=None,
        metavar="P",
        help="Min selection probability to call a gene stably associated (default: 0.6).",
    )
    p.add_argument(
        "--top-n-genes",
        type=int,
        default=None,
        metavar="N",
        help="Top N stable genes to pass to pathway enrichment (default: 50).",
    )
    p.add_argument(
        "--detection-threshold",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Minimum selected-method viral molecule estimate for a candidate-support label "
            "(default: from config or 1)."
        ),
    )
    p.add_argument(
        "--label",
        choices=("raw", "cpm", "fraction"),
        default="raw",
        help=(
            "Candidate-support label: 'raw' (default, depth-confounded estimate>=threshold) or "
            "depth-normalized 'cpm'/'fraction' (prevalence-matched burden per host molecule)."
        ),
    )
    p.add_argument(
        "--depth-match",
        action="store_true",
        default=False,
        help=(
            "Restrict analysis to a coarsened-exact depth-matched cohort; this reduces measured "
            "total-depth imbalance but does not establish independence from depth confounding."
        ),
    )
    p.add_argument(
        "--mito-control",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add %%mito as a covariate to the per-gene E-values (default: on; --no-mito-control to disable).",
    )
    p.add_argument(
        "--gene-symbols",
        action="store_true",
        default=False,
        help="Annotate output CSVs with HGNC symbols from Ensembl IDs via mygene.info (network).",
    )
    p.add_argument(
        "--differential",
        action="store_true",
        default=False,
        help="Write a genome-wide depth/%%mito-adjusted differential table (<virus>_differential.csv).",
    )
    p.add_argument(
        "--enrichment",
        action="store_true",
        default=False,
        help="Run pathway enrichment via gget (requires viralscan[enrichment]).",
    )
    p.add_argument(
        "--enrichment-db",
        default=None,
        metavar="DB",
        help="gget.enrichr database (default: GO_Biological_Process_2023).",
    )
    _add_verbosity_args(p)
    p.set_defaults(_subcommand="hostresponse")


def _run_hostresponse_subcommand(args: argparse.Namespace) -> None:
    """Run host-response analysis on an existing viralscan output directory."""
    from viralscan.kb_outputs import KbCountOutputs
    from viralscan.runconfig import RunConfig
    from viralscan.scripts.hostresponse import DEFAULT_SEEDS, run_hostresponse

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

    n_seeds = args.n_seeds if args.n_seeds is not None else (cfg.hostresponse_n_seeds or 6)
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
        seeds=DEFAULT_SEEDS[:n_seeds],
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
    )
    restamp_run_complete(output_dir.parent)
    log.info("hostresponse complete. Results in %s", out_dir)


def _build_check_whitelist_parser(subparsers: Any) -> None:
    """Register the 'check-whitelist' diagnostic subcommand."""
    p = subparsers.add_parser(
        "check-whitelist",
        help="Check that R1 barcodes match a whitelist (chemistry-mismatch preflight).",
        description=(
            "Sample the first reads of R1, extract the cell barcode with the given\n"
            "technology's geometry, and report the fraction that match the whitelist.\n"
            "A low match rate means --technology/--whitelist do not match the library\n"
            "chemistry, which makes bustools silently discard most reads (finding F-005).\n\n"
            "Example:\n"
            "  viralscan check-whitelist -s1 R1.fastq.gz -w 3M-february-2018.txt -x 10xv3"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--sample1", "-s1", required=True, metavar="R1", help="R1 FASTQ (barcode read).")
    p.add_argument(
        "--whitelist",
        "-w",
        required=True,
        metavar="PATH",
        help="Barcode whitelist (optionally .gz).",
    )
    p.add_argument(
        "--technology", "-x", default="10xv3", help="Single-cell technology (default: 10xv3)."
    )
    p.add_argument(
        "--min-match-rate",
        type=float,
        default=0.5,
        metavar="F",
        help="Minimum acceptable barcode match rate (default: 0.5).",
    )
    p.add_argument(
        "--n-sample",
        type=int,
        default=100_000,
        metavar="N",
        help="Number of R1 reads to sample (default: 100000).",
    )
    _add_verbosity_args(p)
    p.set_defaults(_subcommand="check-whitelist")


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


def _build_check_chemistry_parser(subparsers: Any) -> None:
    """Register the 'check-chemistry' diagnostic subcommand (CHEM-01)."""
    p = subparsers.add_parser(
        "check-chemistry",
        help="Diagnose single-cell vs bulk, chemistry, 3'/5' end and strand before a run.",
        description=(
            "Read-only library check. From the first R1 reads it reports bulk vs single-cell,\n"
            "the on-list/chemistry call and UMI length. With -i/-t/-s2 it also runs the\n"
            "1M-pair strand pilot, which tells 3' from 5' on trimmed R1 and picks --strand.\n"
            "Prints the -x/-w/--strand to use. Exit 1 when the chemistry is unresolved or\n"
            "the evidence conflicts.\n\n"
            "Example:\n"
            "  viralscan check-chemistry -s1 R1.fq.gz -s2 R2.fq.gz -i index.idx -t t2g.txt"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--sample1", "-s1", required=True, metavar="R1", help="R1 FASTQ (barcode read).")
    p.add_argument("--sample2", "-s2", default=None, metavar="R2", help="R2 FASTQ (for the pilot).")
    p.add_argument("--index", "-i", default=None, help="kallisto index (enables the strand pilot).")
    p.add_argument("--transcripts", "-t", default=None, help="t2g file (enables the strand pilot).")
    p.add_argument("--whitelist", "-w", default=None, metavar="PATH", help="Barcode on-list.")
    p.add_argument("--technology", "-x", default=None, help="Technology you plan to use; checked.")
    p.add_argument("--cores", "-c", type=int, default=4, help="Threads for the pilot.")
    p.add_argument("--json", default=None, metavar="PATH", help="Also write the report as JSON.")
    _add_verbosity_args(p)
    p.set_defaults(_subcommand="check-chemistry")


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


def build_parser() -> argparse.ArgumentParser:
    """Build the full argument parser (all subcommands) without parsing."""
    parser = argparse.ArgumentParser(
        usage="\n\033[96m"
        + figlet_format("Welcome to ViralScan", font="big", width=200)
        + "\033[0m",
        prog="viralscan",
        allow_abbrev=False,
        description=(
            "ViralScan — viral load quantification from single-cell RNA-seq.\n\n"
            "Subcommands:\n"
            "  (default)       Quantify viral load from FASTQ samples.\n"
            "  data fetch      Download the viral annotation panel from Zenodo.\n"
            "  build-ref       Build a combined host + virus kallisto reference.\n"
            "  evidence        Trace/extract the reads behind viral calls.\n"
            "  rerun-multimap  Switch multimapping method without redoing kb count.\n"
            "  hostresponse    Run host-response analysis on a completed viralscan run.\n"
            "  check-whitelist Check barcode whitelist/chemistry compatibility.\n"
            "  check-chemistry Diagnose bulk/single-cell, chemistry, 3'/5' and strand.\n\n"
            "Recommended host-aware workflow: run 'viralscan build-ref' once, "
            "then quantify with the generated -i/-t files.\n\n"
            "There are 3 ways to run the default (quantification) mode:\n"
            "  1. Provide a pre-built kallisto index (-t / -i).\n"
            "  2. Provide FASTA + GTF (--reference -fasta ... -gtf ...).\n"
            "  3. Provide NCBI accession numbers (-acc ...) — ViralScan fetches + builds.\n\n"
            "Examples:\n"
            "  viralscan -t t2g.txt -i index.idx -o out/ -s1 R1.fastq.gz -s2 R2.fastq.gz\n"
            "  viralscan build-ref --host human --virus-accessions NC_045512.2 -o ref/\n\n"
            "Run 'viralscan --help' or 'viralscan build-ref --help' for full options."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    from viralscan import __version__

    parser.add_argument("--version", action="version", version=f"viralscan {__version__}")

    subparsers = parser.add_subparsers(dest="_subcommand")
    _build_data_parser(subparsers)
    _build_ref_parser(subparsers)
    _build_evidence_parser(subparsers)
    _build_rerun_multimap_parser(subparsers)
    _build_rerun_programs_parser(subparsers)
    _build_doctor_parser(subparsers)
    _build_validate_run_parser(subparsers)
    _build_hostresponse_parser(subparsers)
    _build_check_whitelist_parser(subparsers)
    _build_check_chemistry_parser(subparsers)

    # ── default (quantification) arguments ────────────────────────────────
    parser.add_argument(
        "--output",
        "-o",
        required=False,
        default=None,
        help="The path to the output directory (required for quantification).",
    )
    parser.add_argument(
        "--sample1",
        "-s1",
        default=None,
        help="The path to the forward FASTQ sample (gunzipped is preferred).",
    )
    parser.add_argument(
        "--sample2",
        "-s2",
        default=None,
        help="The path to the backward FASTQ sample (gunzipped is preferred).",
    )

    parser.add_argument(
        "--transcripts",
        "-t",
        default=None,
        help="The path to the transcripts (t2g) file produced by kb ref.",
    )
    parser.add_argument(
        "--index", "-i", default=None, help="The path to the reference index created by kb ref."
    )
    parser.add_argument(
        "--cores",
        "-c",
        default=6,
        type=int,
        help="The amount of cores the workflow can use. Default: 6.",
    )
    parser.add_argument(
        "--reference",
        "-ref",
        action="store_true",
        help="Build a kb ref index from -fasta and -gtf into the output directory.",
    )
    parser.add_argument(
        "--gtf",
        "-gtf",
        default=None,
        help="Path to GTF files (comma-delimited, without space in-between).",
    )
    parser.add_argument(
        "--fasta",
        "-fasta",
        default=None,
        help="Path to FASTA files (comma-delimited, without space in-between).",
    )
    parser.add_argument(
        "--f1",
        "-f1",
        default=None,
        help="Path to the cDNA FASTA (lamanno, nucleus) or mismatch FASTA (kite) to be generated",
    )
    parser.add_argument(
        "--visual",
        "-v",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add visualizations to the output. Use --no-visual to disable. Default: True.",
    )
    parser.add_argument(
        "--technology",
        "-x",
        default=None,
        help=(
            "Single-cell technology (`kb --list` to view). Default: detected from the "
            "first 100k R1 reads of each sample (on-list match, TSO/poly-T position). "
            "An explicit -x that the reads contradict, or reads that fit no single "
            "chemistry, stop the run. GEM-X 5' needs its on-list via -w."
        ),
    )
    parser.add_argument(
        "--force-technology",
        action="store_true",
        help="Run with the explicit -x even when the chemistry check disagrees or "
        "cannot decide. The detection is still logged.",
    )
    parser.add_argument(
        "--whitelist",
        "-w",
        default=None,
        help="Path to file of whitelisted barcodes. If absent, kb-python's bundled whitelist is used.",
    )
    parser.add_argument(
        "--strand",
        choices=["forward", "reverse", "unstranded", "auto"],
        default=None,
        help=(
            "Read strandedness passed to `kb count --strand`. Default (no flag): "
            "with the pinned kallisto 0.52.0 the effective behaviour is "
            "`unstranded` for every chemistry — kallisto 0.52.0 marks the 10x "
            "technologies strand-specific yet leaves the strand unset, so the "
            "legacy per-technology forward default no longer engages (upstream "
            "regression; kallisto <=0.51.1 defaulted 10x to forward). Set the "
            "flag explicitly for stranded quantification. 10x 5' libraries need "
            "`reverse` or `unstranded`: forward pseudoaligns only 6.5-8.9%% of "
            "reads (F-020). `auto` (opt-in) pilots forward/reverse/unstranded on the first "
            "1M read pairs of each sample and picks one; the choice is recorded in "
            "run_manifest.json."
        ),
    )
    parser.add_argument(
        "--multimapping",
        "-mm",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Take multimapping into account. Use --no-multimapping to disable. Default: True.",
    )
    parser.add_argument(
        "--umap",
        "-umap",
        action="store_true",
        help="Generate a UMAP plot. Significantly increases runtime. Default: off.",
    )

    parser.add_argument(
        "--ncbi-accession",
        "-acc",
        default=None,
        help="One or more NCBI nucleotide accessions (e.g. 'NC_002021.3'), comma-separated. "
        "ViralScan will download FASTA + GTF for each and build the index. "
        "Mutually exclusive with --reference / -fasta / -gtf.",
    )
    parser.add_argument(
        "--ncbi-email",
        default=None,
        help="Contact email for NCBI E-utilities. Falls back to $NCBI_EMAIL.",
    )
    parser.add_argument(
        "--ncbi-api-key",
        default=None,
        help="NCBI API key for higher E-utilities request rates with "
        "--ncbi-accession. Falls back to $NCBI_API_KEY.",
    )
    parser.add_argument(
        "--data-cache-dir",
        default=None,
        metavar="PATH",
        help=(
            "Root cache directory for ViralScan's fetched viral annotation panel. "
            "Default: $VIRALSCAN_CACHE or ~/.cache/viralscan/."
        ),
    )

    # Reporting / threshold parameters
    parser.add_argument(
        "--se-threshold",
        type=int,
        default=DEFAULTS["se_threshold"],
        help=(
            "Selected-method viral molecule estimate above which a cell is flagged as a "
            "'super-expressor'. "
            f"Default: {DEFAULTS['se_threshold']}."
        ),
    )
    parser.add_argument(
        "--detection-threshold",
        type=int,
        default=DEFAULTS["detection_threshold"],
        help=(
            "Minimum total selected-method viral molecule estimate required for candidate "
            "detection support. "
            f"Default: {DEFAULTS['detection_threshold']}."
        ),
    )
    parser.add_argument(
        "--positive-control-gene",
        default=None,
        metavar="GENE_ID",
        help=(
            "Gene ID of a spike-in planted at a known molecule count. It is the only "
            "way to measure the k-mer capture term, and therefore the only way to turn "
            "'no virus detected' into a certifiable negative rather than a sampling "
            "statement. Must be given together with --positive-control-molecules."
        ),
    )
    parser.add_argument(
        "--positive-control-molecules",
        type=float,
        default=None,
        metavar="N",
        help=(
            "Molecules of --positive-control-gene planted in the library. Capture is "
            "measured as observed/N and reported in results/positive_control.json; "
            "capture=1.0 means no loss was measurable and implies nothing about "
            "sequence divergence."
        ),
    )
    parser.add_argument(
        "--positive-control-scope",
        choices=CAPTURE_SCOPES,
        default=None,
        help=(
            "What the positive control may certify. exact_sequence: only the virus row "
            "named by --positive-control-virus-key, and only that exact sequence. "
            "panel_mechanics: the pipeline recovers a planted molecule; certifies no "
            "virus. virus_key: needs an approved transfer calibration and is currently "
            "rejected. Unset on a configured control behaves as panel_mechanics and "
            "warns once."
        ),
    )
    parser.add_argument(
        "--positive-control-virus-key",
        default=None,
        metavar="VIRUS_KEY",
        help=(
            "Virus row (the name in results/sensitivity.tsv) an exact_sequence control "
            "was measured on. Required with --positive-control-scope exact_sequence."
        ),
    )
    parser.add_argument(
        "--require-positive-control",
        action=argparse.BooleanOptionalAction,
        default=DEFAULTS["require_positive_control"],
        help=(
            "Fail the run when nothing is detected and no positive control could "
            "measure a capture term. A mechanics/recovery gate: passing it does not "
            "certify any virus outside the control's declared scope. "
            f"Default: {DEFAULTS['require_positive_control']}."
        ),
    )
    parser.add_argument(
        "--anellovirus-gene-ids",
        action=argparse.BooleanOptionalAction,
        default=DEFAULTS["anellovirus_gene_ids"],
        help=(
            "Treat the expanded anellovirus panel's {accession}_geneN IDs as viral. "
            "Required for any reference built with `viralscan build-ref "
            "--reference-panel anellovirus` (or the bundled-panel builder), because "
            "those GTFs are materialized into the index rather than the panel "
            "directory. Off means 2,021 of 2,041 anellovirus genomes are counted but "
            "never reported. "
            f"Default: {DEFAULTS['anellovirus_gene_ids']}."
        ),
    )
    parser.add_argument(
        "--gene-programs",
        action=argparse.BooleanOptionalAction,
        default=DEFAULTS["gene_programs"],
        help=(
            "Second layer: for viruses the detection rule already called, infer the "
            "viral gene programme (latent vs productive) per cell from uniquely-placing "
            "molecules, and write results/gene_program_summary.tsv and "
            "results/gene_program_cells.tsv. Only nine viruses have a programme model; "
            "for the rest a 'not_applicable' row is emitted so silence is not read as "
            "'programme not detected'. Off by default. "
            f"Default: {DEFAULTS['gene_programs']}."
        ),
    )
    parser.add_argument(
        "--anello-align",
        action=argparse.BooleanOptionalAction,
        default=DEFAULTS["anello_align"],
        help=(
            "Anellovirus alignment branch (PLAN ANDET-09): align every host-unmapped "
            "read to the panel's anellovirus genomes with STARsolo, independent of "
            "kallisto, and add alignment_* evidence columns plus detection_source to "
            "viral_summary.tsv. Runs only with --host-filter starsolo and an "
            "anello_star/ index next to the kb index (built by "
            "scripts/build_bundled_panel_ref.py); otherwise alignment_status records "
            "why it was skipped. The columns are labels, never filters. "
            f"Default: {DEFAULTS['anello_align']}."
        ),
    )
    parser.add_argument(
        "--programme-min-breadth",
        type=int,
        default=DEFAULTS["programme_min_breadth"],
        metavar="N",
        help=(
            "Distinct non-overlapping overlap groups required before a programme is "
            "called. Counted in overlap groups rather than genes because EBV's latent "
            "and lytic ORFs share exonic sequence: on the EBV LCL run a naive per-gene "
            "comparison gives a latent:lytic ratio of 1.15 in a cell line defined by "
            "latency. Must be >= 1. "
            f"Default: {DEFAULTS['programme_min_breadth']}."
        ),
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=DEFAULTS["min_counts"],
        help=(
            f"Minimum total molecule count per cell (for UMAP QC). Default: {DEFAULTS['min_counts']}."
        ),
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=DEFAULTS["min_genes"],
        help=(
            "Minimum detected genes per cell (for UMAP QC filter). "
            f"Default: {DEFAULTS['min_genes']}."
        ),
    )
    parser.add_argument(
        "--hvg-min-mean",
        type=float,
        default=DEFAULTS["hvg_min_mean"],
        help=(
            f"Scanpy highly-variable-gene min_mean parameter. Default: {DEFAULTS['hvg_min_mean']}."
        ),
    )
    parser.add_argument(
        "--hvg-max-mean",
        type=float,
        default=DEFAULTS["hvg_max_mean"],
        help=(
            f"Scanpy highly-variable-gene max_mean parameter. Default: {DEFAULTS['hvg_max_mean']}."
        ),
    )
    parser.add_argument(
        "--hvg-min-disp",
        type=float,
        default=DEFAULTS["hvg_min_disp"],
        help=(
            f"Scanpy highly-variable-gene min_disp parameter. Default: {DEFAULTS['hvg_min_disp']}."
        ),
    )
    parser.add_argument(
        "--umap-n-neighbors",
        type=int,
        default=DEFAULTS["umap_n_neighbors"],
        help=(
            "Number of neighbors for Scanpy graph construction before UMAP. "
            f"Default: {DEFAULTS['umap_n_neighbors']}."
        ),
    )
    parser.add_argument(
        "--cell-calling",
        choices=CELL_CALLING_METHODS,
        default=DEFAULTS["cell_calling"],
        help=(
            "How to identify real (non-empty-droplet) cells so viral rates are reported "
            "over called cells (primary) as well as all barcodes (secondary). "
            "'auto' uses --called-cells-file when supplied, otherwise EmptyDrops; "
            "'external' requires that called-cell list; "
            "'emptydrops' runs DropletUtils::emptyDrops (needs R); "
            "'knee' is a sensitivity-only barcode-rank approximation, not for reported "
            "numbers (it calls empty droplets as cells on real libraries); 'none' = all barcodes. "
            f"Default: {DEFAULTS['cell_calling']}."
        ),
    )
    parser.add_argument(
        "--called-cells-file",
        default=None,
        help=(
            "Path to an external called-cell barcode list (one per line, optional -1 suffix) "
            "used when --cell-calling external. Typically a CellRanger/STARsolo "
            "filtered barcodes.tsv(.gz)."
        ),
    )
    parser.add_argument(
        "--emptydrops-seed",
        type=int,
        default=DEFAULTS["emptydrops_seed"],
        help=(
            "Random seed for DropletUtils::emptyDrops. It is a Monte-Carlo test, so "
            "this decides which borderline barcodes are called. Set it from the "
            "preregistered seed when running under a frozen protocol. "
            f"Default: {DEFAULTS['emptydrops_seed']}."
        ),
    )
    parser.add_argument(
        "--emptydrops-niters",
        type=int,
        default=DEFAULTS["emptydrops_niters"],
        help=(
            "Monte-Carlo iterations for DropletUtils::emptyDrops. "
            f"Default: {DEFAULTS['emptydrops_niters']}."
        ),
    )
    parser.add_argument(
        "--multimap-method",
        choices=MULTIMAP_METHODS,
        default=DEFAULTS["multimap_method"],
        help=(
            "How to allocate multi-gene EC counts. "
            "'equal' splits reads equally (fast, good first pass); "
            "'host-conservative' excludes host-virus ambiguous EC mass from viral genes "
            "(recommended for combined host+virus references); "
            "'unique-weighted' weights by unique-gene evidence; 'em-global' fits a "
            "sample-wide model; 'em-cell' uses cell-local evidence with a global prior. "
            "Use 'viralscan rerun-multimap' to switch methods after the run without "
            "redoing the pseudoalignment. "
            f"Default: {DEFAULTS['multimap_method']}."
        ),
    )
    parser.add_argument(
        "--multimap-pseudocount",
        type=float,
        default=DEFAULTS["multimap_pseudocount"],
        help=(
            "Positive pseudocount used by --multimap-method unique-weighted. "
            f"Default: {DEFAULTS['multimap_pseudocount']}."
        ),
    )
    parser.add_argument(
        "--multimap-primary-call",
        choices=MULTIMAP_PRIMARY_CALLS,
        default=DEFAULTS["multimap_primary_call"],
        help=(
            "V3 detection-count contract. The only supported value is 'selected-method': "
            "summaries use the complete molecule matrix for --multimap-method, while unique "
            "and ambiguous evidence remain separately reported. "
            f"Default: {DEFAULTS['multimap_primary_call']}."
        ),
    )
    parser.add_argument(
        "--multimap-em-max-iter",
        type=int,
        default=DEFAULTS["multimap_em_max_iter"],
        help=(
            "Maximum EM iterations for --multimap-method em-global or em-cell. "
            f"Default: {DEFAULTS['multimap_em_max_iter']}."
        ),
    )
    parser.add_argument(
        "--multimap-em-tol",
        type=float,
        default=DEFAULTS["multimap_em_tol"],
        help=(
            "EM convergence tolerance for --multimap-method em-global or em-cell. "
            f"Default: {DEFAULTS['multimap_em_tol']}."
        ),
    )
    parser.add_argument(
        "--cell-types",
        default=None,
        help="Path to a CSV (barcode,cell_type) providing cell-type labels for per-type viral "
        "enrichment in the HTML report. Optional.",
    )
    parser.add_argument(
        "--host-filter",
        choices=["starsolo"],
        default=None,
        metavar="ALIGNER",
        help=(
            "Optional advanced host-subtraction pre-step before viral quantification. "
            "Removes reads that align to the host genome, reducing false positives. "
            "V3 supports 'starsolo' (full-genome STAR alignment) because it preserves exact "
            "fragment identity. Requires --host-index. "
            "Usually not needed when using a combined host+virus reference."
        ),
    )
    parser.add_argument(
        "--read-filter",
        choices=["off", "artefact", "tso-trim"],
        default=None,
        help=(
            "Read-artefact filter before kb count (after --host-filter, if any). 'off' "
            "(the default) counts every read. 'artefact' drops pairs whose R1 carries the "
            "10x TSO in the barcode/UMI span, whose R2 is reagent (TSO|poly-T, TruSeq "
            "chimera), or whose R2 has no complex body before its first >=15 nt "
            "homopolymer in either orientation. It also drops genuine low-complexity "
            "viral reads. 'tso-trim' keeps every pair but cuts the 10x TSO (30 nt) from "
            "the start of R2, because TSO-led host reads that STAR's host filter misses "
            "pseudoalign to viral CAG tracts (F-028). Writes "
            "read_filtered/read_filter_audit.tsv and fragment_lineage.tsv.gz."
        ),
    )
    parser.add_argument(
        "--host-filter-star-params",
        choices=["pinned", "star-default"],
        default=None,
        help=(
            "STAR filter parameter set for --host-filter starsolo. 'pinned' (the default) "
            "removes only near-identical, near-full-length host alignments (at most 4 "
            "mismatches, 90%% of the read aligned, up to 20 loci). 'star-default' uses "
            "STAR's own defaults (10 mismatches, 66%% aligned, up to 10 loci). Both are "
            "passed explicitly and recorded in host_filter_audit.tsv."
        ),
    )
    parser.add_argument(
        "--host-index",
        default=None,
        metavar="PATH",
        help=(
            "Path to the STAR host-genome directory required by --host-filter, built with "
            "STAR --runMode genomeGenerate."
        ),
    )

    # ── Host-response module (optional) ──────────────────────────────────────
    parser.add_argument(
        "--host-h5ad",
        default=None,
        metavar="PATH",
        help=(
            "Path to a host gene-expression h5ad file (cells × host genes, log-normalised or raw). "
            "When provided, ViralScan trains per-virus logistic regression models predicting virus "
            "presence from host gene expression (Luebbert et al. 2026 approach) and writes results "
            "to <output>/hostresponse/."
        ),
    )
    parser.add_argument(
        "--hostresponse-n-seeds",
        type=int,
        default=None,
        metavar="N",
        help="Number of random seeds for the multi-seed L2 logistic regression (default: 6).",
    )
    parser.add_argument(
        "--hostresponse-n-stab-iter",
        type=int,
        default=None,
        metavar="N",
        help="Iterations for randomized Lasso stability selection (default: 100).",
    )
    parser.add_argument(
        "--hostresponse-use-hvg",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use highly variable genes as features (default: on). --no-hostresponse-use-hvg uses all genes.",
    )
    parser.add_argument(
        "--hostresponse-stab-min-prob",
        type=float,
        default=None,
        metavar="PROB",
        help="Minimum stability probability to call a gene stably selected (default: 0.6).",
    )
    parser.add_argument(
        "--hostresponse-top-n-genes",
        type=int,
        default=None,
        metavar="N",
        help="Top N stable genes to pass to pathway enrichment (default: 50).",
    )
    parser.add_argument(
        "--hostresponse-label",
        choices=("raw", "cpm", "fraction"),
        default=None,
        help=(
            "Virus-support labeling strategy for hostresponse (default: raw = selected-method "
            "molecule estimate >= detection_threshold). 'cpm'/'fraction' normalize by measured "
            "host molecule depth but do not remove all depth confounding."
        ),
    )
    parser.add_argument(
        "--hostresponse-depth-match",
        action="store_true",
        default=False,
        help=(
            "Restrict hostresponse to a depth-matched cohort (removes depth as a design-level "
            "confounder). Recommended when depth_alone_auc is close to model_auc."
        ),
    )
    parser.add_argument(
        "--hostresponse-control-mito",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Include %%mito as a covariate in hostresponse depth-adjusted E-values (default: on). "
            "Disable with --no-hostresponse-control-mito if the host h5ad has no mitochondrial genes."
        ),
    )
    parser.add_argument(
        "--hostresponse-differential",
        action="store_true",
        default=False,
        help=(
            "Run a genome-wide depth-and-mito-adjusted differential expression test alongside "
            "the stability-selection model (written to <output>/hostresponse/<virus>_differential.csv)."
        ),
    )
    parser.add_argument(
        "--enrichment",
        action="store_true",
        default=False,
        help="Run pathway enrichment on stable host genes via gget.enrichr (requires gget; install with pip install 'ViralScan[enrichment]').",
    )
    parser.add_argument(
        "--enrichment-db",
        default=None,
        metavar="DB",
        help="Enrichment database for gget.enrichr (default: GO_Biological_Process_2023).",
    )

    output_mode = parser.add_mutually_exclusive_group()
    output_mode.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Resume only when run_manifest.json exactly matches this invocation.",
    )
    output_mode.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Explicitly replace a non-empty output directory after confirmation.",
    )
    parser.add_argument(
        "--yes",
        "-y",
        action="store_true",
        default=False,
        help="Answer the --overwrite confirmation; does not imply overwrite or resume.",
    )

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Enable DEBUG-level logging.",
    )
    verbosity.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress INFO messages; only show warnings and errors.",
    )

    return parser


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
    validation checkpoint); :func:`_build_config_args` serialises it via
    :meth:`~viralscan.runconfig.RunConfig.to_snakemake_config_args`.
    This eliminates the previously hand-maintained parallel key list and
    ensures CLI flags like ``--multimap-em-max-iter`` are never accidentally
    omitted from the Snakemake invocation.
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
            "whitelist": args.whitelist,
            "strand": getattr(args, "strand", None),
            "multimapping": args.multimapping,
            "se_threshold": args.se_threshold,
            "detection_threshold": args.detection_threshold,
            "positive_control_gene": getattr(args, "positive_control_gene", None),
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
            "programme_min_breadth": getattr(
                args, "programme_min_breadth", DEFAULTS["programme_min_breadth"]
            ),
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
            "cell_types": args.cell_types,
            "data_cache_dir": args.data_cache_dir,
            "host_filter_aligner": getattr(args, "host_filter", None),
            "host_filter_star_params": getattr(args, "host_filter_star_params", None),
            "read_filter": getattr(args, "read_filter", None),
            "host_index": getattr(args, "host_index", None),
            "host_h5ad": getattr(args, "host_h5ad", None),
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


def _build_config_args(
    args: argparse.Namespace,
    outs: str,
    index: str,
    transcripts: str,
    f1: Optional[str],
    s1: str,
    s2: str,
) -> list[str]:
    """Build the Snakemake ``--config k=v`` list for one sample."""
    return _build_run_config(args, outs, index, transcripts, f1, s1, s2).to_snakemake_config_args()


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
    from viralscan.run_safety import sha256_file
    from viralscan.virus_identity import gtf_gene_ids, write_build_manifest_from_t2g

    path = write_build_manifest_from_t2g(
        index,
        t2g,
        gtf_gene_ids(gtf),
        provenance={
            "builder": "viralscan --reference",
            "fasta": {"path": str(Path(fasta).resolve()), "sha256": sha256_file(Path(fasta))},
            "gtf": {"path": str(Path(gtf).resolve()), "sha256": sha256_file(Path(gtf))},
        },
    )
    log.info("Index build manifest: %s", path)


def _build_kb_ref(output_dir: Path, fasta: str, gtf: str) -> tuple[str, str, str]:
    """Run ``kb ref`` to build an index. Returns (transcripts, index, f1) paths."""
    from viralscan.scripts.build_reference import _run_kb_ref

    index_dir = output_dir / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    fasta_input, gtf_input = _prepare_kb_ref_inputs(output_dir, fasta, gtf)
    transcripts = str(index_dir / "t2g.txt")
    index = str(index_dir / "index.idx")
    f1 = str(index_dir / "cdna.fa")
    log.info("Building kb ref index. Depending on the genome this can take a while...")
    _run_kb_ref(
        [
            "kb",
            "ref",
            "-i",
            index,
            "-g",
            transcripts,
            "-f1",
            f1,
            "--overwrite",
            fasta_input,
            gtf_input,
        ],
        index_dir,
    )
    log.info("Reference index is done!")
    _write_reference_manifest(index, transcripts, fasta_input, gtf_input)
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
        config_args = run_config.to_snakemake_config_args()
        if args.resume:
            _backfill_identity_table(run_config)
        subprocess.run(_snakemake_run_command(snakefile_path, args.cores, config_args), check=True)

        _write_sample_summary(outs, time.time() - sample_start, n_transcripts, n_genes)

        unlock_cmd = [
            "snakemake",
            "--snakefile",
            snakefile_path,
            "--unlock",
            "--config",
            *config_args,
        ]
        subprocess.run(unlock_cmd, check=True)

    # Every sample finished (check=True above): the run is complete.
    write_run_complete(output_dir)


if __name__ == "__main__":
    main()

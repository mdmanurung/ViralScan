"""
This file is the backbone of the framework. It checks users' input and calls
the snakemake workflow. It handles the Argument Parser, showing the help function.
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, NoReturn

from viralscan.defaults import DEFAULTS, MULTIMAP_METHODS, MULTIMAP_PRIMARY_CALLS
from viralscan.utils import configure_logging, split_comma_paths

try:
    from pyfiglet import figlet_format as _figlet_format
except ImportError:  # pyfiglet is optional

    def _figlet_format(text: str, **kwargs: object) -> str:
        return text


figlet_format = _figlet_format

log = logging.getLogger("viralscan")


REQUIRED_TOOLS = ("kb", "snakemake")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


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
    fetch.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Enable DEBUG-level logging.",
    )
    fetch.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress INFO messages.",
    )
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
        "--list-species",
        action="store_true",
        default=False,
        help="Print all supported host species and exit.",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Enable DEBUG-level logging.",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress INFO messages.",
    )
    p.set_defaults(_subcommand="build-ref")


def create_help() -> argparse.Namespace:
    """
    This function creates the help function and handles the Argument Parser.
    ---------------------------------------------------------------------
    Returns:
        args (argparse.Namespace): All arguments given by the user to process
    """
    parser = argparse.ArgumentParser(
        usage="\n\033[96m"
        + figlet_format("Welcome to ViralScan", font="big", width=200)
        + "\033[0m",
        prog="viralscan",
        allow_abbrev=False,
        description=(
            "ViralScan — viral load quantification from single-cell RNA-seq.\n\n"
            "Subcommands:\n"
            "  (default)  Quantify viral load from FASTQ samples.\n"
            "  data fetch Download the viral annotation panel from Zenodo.\n"
            "  build-ref  Build a combined host + virus kallisto reference.\n\n"
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

    subparsers = parser.add_subparsers(dest="_subcommand")
    _build_data_parser(subparsers)
    _build_ref_parser(subparsers)

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
        default="10xv3",
        help="Single-cell technology used (`kb --list` to view). Default: 10xv3.",
    )
    parser.add_argument(
        "--whitelist",
        "-w",
        default=None,
        help="Path to file of whitelisted barcodes. If absent, kb-python's bundled whitelist is used.",
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
            "UMI count above which a cell is called a 'super-expressor'. "
            f"Default: {DEFAULTS['se_threshold']}."
        ),
    )
    parser.add_argument(
        "--detection-threshold",
        type=int,
        default=DEFAULTS["detection_threshold"],
        help=(
            "Minimum total viral UMI required to call a virus detected. "
            f"Default: {DEFAULTS['detection_threshold']}."
        ),
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=DEFAULTS["min_counts"],
        help=(
            f"Minimum total UMI per cell (for UMAP QC filter). Default: {DEFAULTS['min_counts']}."
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
        "--multimap-method",
        choices=MULTIMAP_METHODS,
        default=DEFAULTS["multimap_method"],
        help=(
            "How to allocate multi-gene EC counts. "
            "'host-conservative' excludes host-virus ambiguous EC mass from viral genes "
            "and is the recommended default for combined host+virus references; "
            "'equal' preserves legacy equal splitting; 'unique-weighted' weights by "
            "unique-gene evidence. "
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
            "Detection policy for multimapped viral evidence. 'legacy' preserves current "
            "calls, 'unique-only' calls from unambiguous viral signal, and 'confidence' "
            "keeps legacy calls while reporting confidence tiers. "
            f"Default: {DEFAULTS['multimap_primary_call']}."
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
        choices=["starsolo", "kallisto"],
        default=None,
        metavar="ALIGNER",
        help=(
            "Optional advanced host-subtraction pre-step before viral quantification. "
            "Removes reads that align to the host genome/transcriptome, reducing false positives. "
            "Choices: 'starsolo' (STAR genome alignment) or 'kallisto' (pseudo-alignment "
            "against a host cDNA index). Requires --host-index. "
            "Usually not needed when using a combined host+virus reference."
        ),
    )
    parser.add_argument(
        "--host-index",
        default=None,
        metavar="PATH",
        help=(
            "Path to the host genome/index directory required by --host-filter. "
            "For 'starsolo': path to a STAR genome directory (built with STAR --runMode genomeGenerate). "
            "For 'kallisto': path to a host cDNA kallisto index file (.idx), e.g. built with "
            "'kallisto index -i host.idx host_cdna.fa'."
        ),
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

    return parser.parse_args()


def _die(message: str, code: int = 1) -> NoReturn:
    log.error(message)
    sys.exit(code)


def _has_valid_fastq_suffix(path: str) -> bool:
    return any(path.endswith(suf) for suf in FASTQ_SUFFIXES)


def check_output(args: argparse.Namespace) -> None:
    """
    This function checks whether the given output directory already
    exists and shows options to the user.
    """
    path = args.output
    if not os.path.isdir(path):
        return
    if not os.listdir(path):
        return
    answer = (
        input(
            "\033[33mThe output directory already exists and contains files. "
            "Do you want to overwrite this? (yes/y/no/n): \033[0m"
        )
        .strip()
        .lower()
    )
    if answer in ("no", "n"):
        print("You have chosen not to continue. The code has been terminated.")
        sys.exit(0)
    if answer not in ("yes", "y"):
        _die(f"This is not a valid answer: {answer}. The code has been terminated.")
    print(
        "You have chosen to continue. The contents of the existing output directory will be overwritten."
    )
    for filename in os.listdir(path):
        file_path = os.path.join(path, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)
        else:
            shutil.rmtree(file_path)


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
    if host_filter:
        host_index_path = str(host_index)
        if not os.path.exists(host_index_path):
            _die(f"Host index path does not exist: {host_index_path}.")

    log.info("All input data has been checked and is correct.")


def _check_required_tools() -> None:
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        _die(
            "The following required external tools are not on PATH: "
            f"{', '.join(missing)}. Please install kb-python (provides 'kb') and snakemake."
        )


def _check_host_filter_tools(aligner: str) -> None:
    """Verify that tools required by the chosen host-filter aligner are on PATH."""
    if aligner == "starsolo":
        if shutil.which("STAR") is None:
            _die(
                "--host-filter starsolo requires 'STAR' on PATH. "
                "Install STARsolo: https://github.com/alexdobin/STAR"
            )
    elif aligner == "kallisto":
        missing = [t for t in ("kallisto", "bustools") if shutil.which(t) is None]
        if missing:
            _die(
                f"--host-filter kallisto requires {', '.join(missing)} on PATH. "
                "Install kb-python: https://github.com/pachterlab/kb-python"
            )


def _count_lines(path: str) -> int:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return sum(1 for _ in f)


def _count_unique_genes(t2g_path: str) -> int:
    """Count unique (col2, col3) pairs in a t2g.txt file."""
    seen: set[tuple[str, str]] = set()
    with open(t2g_path, "r", encoding="utf-8", errors="replace") as f:
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


def _config_value(value: object) -> str:
    """Serialize optional Snakemake config values without literal ``None`` sentinels."""
    return "" if value is None else str(value)


def _build_kb_ref(output_dir: Path, fasta: str, gtf: str) -> tuple[str, str, str]:
    """Run ``kb ref`` to build an index. Returns (transcripts, index, f1) paths."""
    index_dir = output_dir / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    fasta_input, gtf_input = _prepare_kb_ref_inputs(output_dir, fasta, gtf)
    transcripts = str(index_dir / "t2g.txt")
    index = str(index_dir / "index.idx")
    f1 = str(index_dir / "cdna.fa")
    log.info("Building kb ref index. Depending on the genome this can take a while...")
    subprocess.run(
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
        check=True,
    )
    log.info("Reference index is done!")
    return transcripts, index, f1


def main() -> None:
    start = time.time()
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
    check_output(args)
    errorhandler(args)

    if args.host_filter:
        _check_host_filter_tools(args.host_filter)

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.ncbi_accession:
        from viralscan.scripts.ncbi_fetch import fetch_reference, NCBIFetchError

        accessions = split_comma_paths(args.ncbi_accession)
        ref_dir = output_dir / "ncbi_reference"
        try:
            fasta_path, gtf_path = fetch_reference(
                accessions=accessions,
                out_dir=ref_dir,
                email=args.ncbi_email,
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

    for s1, s2 in zip(samples1, samples2):
        out = Path(s1).name.split("_")[0]
        outs = os.path.join(output, out) + os.sep
        if args.host_index:
            kb_r1 = os.path.join(outs, "host_filtered", "R1.fastq.gz")
            kb_r2 = os.path.join(outs, "host_filtered", "R2.fastq.gz")
        else:
            kb_r1 = s1
            kb_r2 = s2
        config_args = [
            f"output={outs}",
            f"index={index}",
            f"transcripts={transcripts}",
            f"sample1={s1}",
            f"sample2={s2}",
            f"kb_r1={kb_r1}",
            f"kb_r2={kb_r2}",
            f"cores={args.cores}",
            f"gtf={_config_value(args.gtf)}",
            f"fasta={_config_value(args.fasta)}",
            f"visual={args.visual}",
            f"f1={_config_value(f1)}",
            f"reference={args.reference}",
            f"umap={args.umap}",
            f"technology={args.technology}",
            f"whitelist={_config_value(args.whitelist)}",
            f"multimapping={args.multimapping}",
            f"se_threshold={args.se_threshold}",
            f"detection_threshold={args.detection_threshold}",
            f"min_counts={args.min_counts}",
            f"min_genes={args.min_genes}",
            f"hvg_min_mean={args.hvg_min_mean}",
            f"hvg_max_mean={args.hvg_max_mean}",
            f"hvg_min_disp={args.hvg_min_disp}",
            f"umap_n_neighbors={args.umap_n_neighbors}",
            f"multimap_method={args.multimap_method}",
            f"multimap_pseudocount={args.multimap_pseudocount}",
            f"multimap_primary_call={args.multimap_primary_call}",
            f"cell_types={_config_value(args.cell_types)}",
            f"data_cache_dir={_config_value(args.data_cache_dir)}",
            f"host_filter_aligner={args.host_filter or ''}",
            f"host_index={args.host_index or ''}",
        ]
        cmd = [
            "snakemake",
            "--snakefile",
            snakefile_path,
            "--cores",
            str(args.cores),
            "--use-conda",
            "--quiet",
            "all",
            "--config",
            *config_args,
        ]
        subprocess.run(cmd, check=True)

        end = time.time()
        summary_path = os.path.join(outs, "summary.txt")
        os.makedirs(outs, exist_ok=True)
        n_transcripts = _count_lines(transcripts)
        n_genes = _count_unique_genes(transcripts)
        with open(summary_path, "a", encoding="utf-8") as f:
            f.write(f"\nRuntime: {end - start:.4f} seconds.\n\n")
            f.write(f"Amount of transcripts in data: {n_transcripts}\n")
            f.write(f"Amount of genes in data: {n_genes}\n")

        unlock_cmd = [
            "snakemake",
            "--snakefile",
            snakefile_path,
            "--unlock",
            "--config",
            *config_args,
        ]
        subprocess.run(unlock_cmd, check=True)


if __name__ == "__main__":
    main()

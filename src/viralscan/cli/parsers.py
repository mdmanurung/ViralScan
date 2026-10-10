"""Argument parsers for the ``viralscan`` CLI (PLAN ``SW-09``).

Moved verbatim out of ``menu.py``: ``build_parser`` and one ``_build_*_parser`` per subcommand.
Nothing here runs a command. ``menu`` re-exports every name so ``from viralscan.menu import
build_parser`` (and the tests that import the private builders) keep working.
"""

import argparse
from typing import Any

from viralscan.defaults import (
    CELL_CALLING_METHODS,
    DEFAULTS,
    MULTIMAP_METHODS,
    MULTIMAP_PRIMARY_CALLS,
)
from viralscan.runconfig import CAPTURE_SCOPES

try:
    from pyfiglet import figlet_format as _figlet_format
except ImportError:  # pyfiglet is optional

    def _figlet_format(text: str, font: str = "standard", **kwargs: Any) -> Any:
        return text


figlet_format = _figlet_format


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
        "--competitor-fasta",
        default=None,
        help="Explicit alternate/related/decoy FASTA; requires --competitor-manifest.",
    )
    p.add_argument(
        "--competitor-manifest",
        default=None,
        help="Checksum-pinned TSV covering target, host and competitor FASTA records.",
    )
    p.add_argument(
        "--blast-tie-delta",
        type=float,
        default=0.0,
        help="Nonnegative bitscore distance retaining top BLAST ties (default: 0).",
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


def _add_hostresponse_evaluation_args(parser: Any, prefix: str = "") -> None:
    """Keep rerun and initial-run evaluation options on the same public contract."""
    for name, kwargs in (
        ("cv", dict(choices=("cell", "group"), help="Cell split or group-disjoint evaluation.")),
        ("groups", dict(help="obs column declaring donor/sample groups.")),
        ("cv-folds", dict(type=int, help="Number of group-disjoint folds (default: 5).")),
        ("cell-type-column", dict(help="obs column used for cell-type baselines and strata.")),
        ("cell-types", dict(nargs="+", help="Declared cell types to analyze separately.")),
        (
            "panel-in-fold",
            dict(
                action=argparse.BooleanOptionalAction,
                help="Fit the stability-selected panel within each training fold.",
            ),
        ),
        ("permutations", dict(type=int, help="Structured-null replicates (default: 0).")),
        (
            "permutation-unit",
            dict(
                choices=("cell_within_block", "group"),
                help="Explicit exchangeability unit for the null.",
            ),
        ),
        ("permutation-block", dict(help="obs column restricting label exchangeability.")),
        (
            "min-negative-cells",
            dict(type=int, help="Minimum negative cells per stratum (default: 10)."),
        ),
        ("min-groups", dict(type=int, help="Minimum groups per stratum (default: 2).")),
    ):
        parser.add_argument(f"--{prefix}{name}", default=None, **kwargs)


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
        help="Random seeds for multi-seed L2 regression, 1-6 (default: from config or 6).",
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
    _add_hostresponse_evaluation_args(p)
    p.set_defaults(_subcommand="hostresponse")


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
        "--require-chemistry-sanity",
        action="store_true",
        help="Stop after kb count if library sanity reports an error or cannot be checked. "
        "Diagnostics are retained; warnings alone do not stop the run.",
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
        "--positive-control-receipt",
        default=None,
        metavar="JSON",
        help="Pinned control identity/measurement receipt required for exact-sequence certification.",
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
        "--multimap-molecule-assignments",
        action=argparse.BooleanOptionalAction,
        default=DEFAULTS["multimap_molecule_assignments"],
        help=(
            "Write results/molecule_assignments.tsv.gz: one row per corrected CB-UMI molecule "
            "with its EC ids, status (unique/ambiguous/unresolved), compatible genes and the "
            "selected --multimap-method's per-gene weights. Optional, off by default, and it "
            "changes no count. Needs one extra pass over the BUS text; the file is large on "
            "deep samples. "
            f"Default: {DEFAULTS['multimap_molecule_assignments']}."
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
    _add_hostresponse_evaluation_args(parser, "hostresponse-")
    parser.add_argument(
        "--hostresponse-n-seeds",
        type=int,
        default=None,
        metavar="N",
        help="Number of random seeds for the multi-seed L2 logistic regression, 1-6 (default: 6).",
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

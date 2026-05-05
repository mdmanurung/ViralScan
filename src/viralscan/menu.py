"""
This file is the backbone of the framework. It checks users' input and calls
the snakemake workflow. It handles the Argument Parser, showing the help function.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from pyfiglet import figlet_format


REQUIRED_TOOLS = ("kb", "snakemake")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def create_help():
    """
    This function creates the help function and handles the Argument Parser.
    ---------------------------------------------------------------------
    Returns:
        args (argparse.Namespace): All arguments given by the user to process
    """
    parser = argparse.ArgumentParser(
        usage='\n\033[96m' + figlet_format("Welcome to ViralScan", font="big", width=200) + '\033[0m',
        prog="ViralScan",
        description="""
        ViralScan is a computational framework which predicts viral counts.

        There are 3 different ways to run ViralScan:
            1. You already have a viral reference index built with kb-python.
            2. You have FASTA + GTF files and want ViralScan to build the index for you (--reference).
            3. You only have NCBI accession numbers; ViralScan will fetch FASTA + GTF and build the index (--ncbi-accession).

        Option 1:
            viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz

        Option 2:
            viralscan -o output/ --reference -fasta fasta.fasta -gtf gtf.gtf -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz

        Option 3:
            viralscan -o output/ -acc NC_002021.3 -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
            viralscan -o output/ -acc NC_002021.3,NC_001512.1 --ncbi-email you@example.org -s1 ... -s2 ...

        For multiple samples, GTFs, or FASTAs: comma-separate the values without spaces.

        Run 'viralscan --help' for the full list of parameters.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--output', "-o", required=True, help="The path to the output directory (required)")
    parser.add_argument('--sample1', "-s1", required=True, help="The path to the forward FASTQ sample (gunzipped is preferred) (required)")
    parser.add_argument('--sample2', "-s2", required=True, help="The path to the backward FASTQ sample (gunzipped is preferred) (required)")

    parser.add_argument('--transcripts', "-t", default=None, help="The path to the transcripts (t2g) file produced by kb ref.")
    parser.add_argument('--index', "-i", default=None, help="The path to the reference index created by kb ref.")
    parser.add_argument('--cores', '-c', default=6, type=int, help="The amount of cores the workflow can use. Default: 6.")
    parser.add_argument('--reference', '-ref', action='store_true',
                        help="Build a kb ref index from -fasta and -gtf into the output directory.")
    parser.add_argument('--gtf', '-gtf', default=None, help="Path to GTF files (comma-delimited, without space in-between).")
    parser.add_argument('--fasta', '-fasta', default=None, help="Path to FASTA files (comma-delimited, without space in-between).")
    parser.add_argument('--f1', '-f1', default=None, help="Path to the cDNA FASTA (lamanno, nucleus) or mismatch FASTA (kite) to be generated")
    parser.add_argument('--visual', '-v', action=argparse.BooleanOptionalAction, default=True,
                        help="Add visualizations to the output. Use --no-visual to disable. Default: True.")
    parser.add_argument('--technology', '-x', default='10xv3', help="Single-cell technology used (`kb --list` to view). Default: 10xv3.")
    parser.add_argument('--whitelist', '-w', default=None,
                        help="Path to file of whitelisted barcodes. If absent, kb-python's bundled whitelist is used.")
    parser.add_argument('--multimapping', '-mm', action=argparse.BooleanOptionalAction, default=True,
                        help="Take multimapping into account. Use --no-multimapping to disable. Default: True.")
    parser.add_argument('--umap', '-umap', action='store_true',
                        help="Generate a UMAP plot. Significantly increases runtime. Default: off.")

    parser.add_argument('--ncbi-accession', '-acc', default=None,
                        help="One or more NCBI nucleotide accessions (e.g. 'NC_002021.3'), comma-separated. "
                             "ViralScan will download FASTA + GTF for each and build the index. "
                             "Mutually exclusive with --reference / -fasta / -gtf.")
    parser.add_argument('--ncbi-email', default=None,
                        help="Contact email for NCBI E-utilities. Falls back to $NCBI_EMAIL.")

    return parser.parse_args()


def _die(message: str, code: int = 1) -> None:
    print(f"\033[31m{message}\033[0m", file=sys.stderr)
    sys.exit(code)


def _has_valid_fastq_suffix(path: str) -> bool:
    return any(path.endswith(suf) for suf in FASTQ_SUFFIXES)


def check_output(args):
    """
    This function checks whether the given output directory already
    exists and shows options to the user.
    """
    path = args.output
    if not os.path.isdir(path):
        return
    if not os.listdir(path):
        return
    answer = input(
        "\033[33mThe output directory already exists and contains files. "
        "Do you want to overwrite this? (yes/y/no/n): \033[0m"
    ).strip().lower()
    if answer in ("no", "n"):
        print("You have chosen not to continue. The code has been terminated.")
        sys.exit(0)
    if answer not in ("yes", "y"):
        _die(f"This is not a valid answer: {answer}. The code has been terminated.")
    print("You have chosen to continue. The contents of the existing output directory will be overwritten.")
    for filename in os.listdir(path):
        file_path = os.path.join(path, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)
        else:
            shutil.rmtree(file_path)


def errorhandler(args):
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
        for fasta in args.fasta.split(','):
            if not os.path.exists(fasta):
                _die(f"FASTA path does not exist: {fasta}.")
        for gtf in args.gtf.split(','):
            if not os.path.exists(gtf):
                _die(f"GTF path does not exist: {gtf}.")
        if args.transcripts is not None or args.index is not None or args.f1 is not None:
            answer = input(
                "\033[33mYou have provided a path to the index or transcripts file but want to "
                "create a reference. The reference will be written to the output directory. "
                "Continue? (yes/y/no/n): \033[0m"
            ).strip().lower()
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

    samples1 = args.sample1.split(',')
    samples2 = args.sample2.split(',')
    if len(samples1) != len(samples2):
        _die("--sample1 and --sample2 must have the same number of comma-separated entries.")
    for s1, s2 in zip(samples1, samples2):
        if not os.path.exists(s1) or not os.path.exists(s2):
            _die(f"Sample path does not exist: {s1} or {s2}.")
        if not _has_valid_fastq_suffix(s1):
            _die(f"The forward sample is not in FASTQ format: {s1}.")
        if not _has_valid_fastq_suffix(s2):
            _die(f"The backward sample is not in FASTQ format: {s2}.")

    print("\033[32mAll input data has been checked and is correct.\033[0m")


def _check_required_tools() -> None:
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        _die(
            "The following required external tools are not on PATH: "
            f"{', '.join(missing)}. Please install kb-python (provides 'kb') and snakemake."
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


def _build_kb_ref(output_dir: Path, fasta: str, gtf: str) -> tuple[str, str, str]:
    """Run ``kb ref`` to build an index. Returns (transcripts, index, f1) paths."""
    index_dir = output_dir / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    transcripts = str(index_dir / "t2g.txt")
    index = str(index_dir / "index.idx")
    f1 = str(index_dir / "cdna.fa")
    print(
        "\033[32mBuilding kb ref index. Depending on the genome this can take a while...\033[0m"
    )
    subprocess.run(
        ["kb", "ref", "-i", index, "-g", transcripts, "-f1", f1, "--overwrite", fasta, gtf],
        check=True,
    )
    print("\033[32mReference index is done!\033[0m")
    return transcripts, index, f1


def main():
    start = time.time()
    args = create_help()
    _check_required_tools()
    check_output(args)
    errorhandler(args)

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.ncbi_accession:
        from viralscan.scripts.ncbi_fetch import fetch_reference, NCBIFetchError

        accessions = [a.strip() for a in args.ncbi_accession.split(',') if a.strip()]
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
        print(f"\033[32mFetched NCBI reference for: {', '.join(accessions)}\033[0m")

    if args.reference:
        transcripts, index, f1 = _build_kb_ref(output_dir, args.fasta, args.gtf)
    else:
        transcripts = args.transcripts
        index = args.index
        f1 = args.f1

    snakefile_path = os.path.join(os.path.dirname(__file__), "Snakefile")
    samples1 = args.sample1.split(',')
    samples2 = args.sample2.split(',')
    output = str(output_dir)

    for s1, s2 in zip(samples1, samples2):
        out = Path(s1).name.split('_')[0]
        outs = os.path.join(output, out) + os.sep
        config_args = [
            f"output={outs}",
            f"index={index}",
            f"transcripts={transcripts}",
            f"sample1={s1}",
            f"sample2={s2}",
            f"gtf={args.gtf}",
            f"fasta={args.fasta}",
            f"visual={args.visual}",
            f"f1={f1}",
            f"reference={args.reference}",
            f"umap={args.umap}",
            f"technology={args.technology}",
            f"whitelist={args.whitelist}",
            f"multimapping={args.multimapping}",
        ]
        cmd = [
            "snakemake",
            "--snakefile", snakefile_path,
            "--cores", str(args.cores),
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
            "--snakefile", snakefile_path,
            "--unlock",
            "--config",
            *config_args,
        ]
        subprocess.run(unlock_cmd, check=True)


if __name__ == "__main__":
    main()

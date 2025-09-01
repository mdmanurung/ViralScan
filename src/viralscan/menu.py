import argparse
import os
import shutil
import subprocess
import time
from pyfiglet import figlet_format

# This is the main function which calls the Snakefile and handles the Argument Parser
def create_help():
    parser = argparse.ArgumentParser(
        usage='\n\033[96m' + figlet_format("Welcome to ViralScan", font="big", width=200) + '\033[0m',
        prog="ViralScan",
        description="""
        ViralScan is a computational framework which predicts viral load in patients.

        Please note that you need a reference index made with kb ref (kb-python) to run this workflow correctly. 
        If you haven't created a reference index add -ref or --reference to the command and the pipeline will create
        one for you. You must also provide a list with gtf files and fasta files.

        Please note if you have not create a reference index yet, it can be created by ViralScan. The name of the new 
        index and new transcripts should be put in at --index and --transcripts.

        Run 'viralscan --help' to show the help function and all the possibilities within the framework. 

        Example of running:
        ViralScan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # viralscan -t /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/transcriptome/t2g_serratus.txt -i /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/transcriptome/index_serratus.idx -o output/ -s1 /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/ebola_samples_test/SRR10307460_1.fastq.gz -s2 /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/ebola_samples_test/SRR10307460_2.fastq.gz 


    # Required arguments
    parser.add_argument('--transcripts', "-t", required=True, help="The path to the transcripts made with kb ref.")
    parser.add_argument('--index', "-i", required=True, help="The path to the reference index created by kb ref.")
    parser.add_argument('--output', "-o", required=True, help="The path to the output directory.")
    parser.add_argument('--sample1', "-s1", required=True, help="The path to the forward FASTQ sample (gunzipped is preferred).")
    parser.add_argument('--sample2', "-s2", required=True, help="The path to the backward sample (gunzipped is preferred).")

    # Optional arguments
    parser.add_argument('--cores', '-c', default=6, type=int, help="The amount of cores the workflow can use (default: 6).")
    parser.add_argument('--reference', '-ref', type=bool, default=False, help="Should kb ref run? If yes, also provide the path to te gtf, genome, fasta and path to the cDNA fasta files (output). [True/False] (default: False).")
    parser.add_argument('--gtf', '-gtf', default=None, help="Path to GTF files (comma-delimited, without space in-between).")
    parser.add_argument('--fasta', '-fasta', default=None, help="Path to FASTA files (comma-delimited, without space in-between).")
    parser.add_argument('--f1', '-f1', default=None, help="Path to the cDNA FASTA (lamanno, nucleus) or mismatch FASTA (kite) to be generated")
    parser.add_argument('--visual', '-v', default=True, type=bool, help="Add visualizations to the output. [True/False] (default: True)")
    parser.add_argument('--umap', '-umap', default=False, type=bool, help="Do you want to create a umap? Please note the running time will increase significantly. [True/False] (default: False)")
    args = parser.parse_args()
    return args


def check_output(args):
    path = args.output
    # check if the path to the output directory exists
    if os.path.isdir(path):
        if os.listdir(path):
            continue_ = input("\033[33mThe output directory already exists and contains files. Do you want to overwrite this? (yes/y/no/n): \033[0m")
            if continue_.lower() == "no" or continue_.lower() == "n":
                print("You have chosen not to continue. The code has been aborted.")
                exit()
            elif continue_.lower() == "yes" or continue_.lower() == "y":
                print(f"You have chosen to continue. The contents of the existing output directory will be overwritten.")
                for filename in os.listdir(path):
                    file_path = os.path.join(path, filename)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                    else: # this is a directory which is going to be deleted
                        shutil.rmtree(file_path)
            else:
                print(f"This is not a valid answer: {continue_}. Aborting code...")
                exit()
    else:
        # the output directory does not exist yet, continue as normal
        pass

def errorhandler(args):
    # Error handling
    if not os.path.exists(args.index):
        print(f"\033[31mPath to index does not exist: {args.index}. Aborting code.\033[0m")
        exit()
    if not os.path.exists(args.transcripts):
        print(f"\033[31mPath to transcripts does not exist: {args.transcripts}. Aborting code.\033[0m")
        exit()
    if not os.path.exists(args.sample1) or not os.path.exists(args.sample2):
        print(f"\033[31mPath to samples is not correct: {args.sample1} or {args.sample2}. Aborting code.\033[0m")
        exit()
    if args.reference:
        if args.gtf is None:
            print(f"\033[31mThe path to the GTF files is not filled in. If you want to create a reference genome, you need to provide a path to the gtf files. Aborting code.\033[0m")
            exit()
        if args.fasta is None:
            print(f"\033[31mThe path to the GTF files is not filled in. If you want to create a reference genome, you need to provide a path to the gtf files. Aborting code.\033[0m")
            exit()
        if not os.path.exists(args.gtf):
            print(f"\033[31mThe path to the GTF files is incorrect. If you want to create a reference genome, you need to provide a path to the gtf files. Given path: {args.gtf}. Aborting code.\033[0m")
            exit()
        if not os.path.exists(args.fasta):
            print(f"\033[31mThe path to the fasta files is incorrect. If you want to create a reference genome, you need to provide a path to the fasta files. Given path: {args.fasta}. Aborting code.\033[0m")
            exit()

    print("\033[32mAll input data has been checked and is correct.\033[0m")


def main():
    start = time.time()
    args = create_help()
    check_output(args)
    errorhandler(args)
    # Check if creating reference is needed
    if args.reference:
         print(f"You have decided a reference still needs to be created. This will be created in the provided output directory. Please note that creating a reference (depending on the genome) can take a while.")
         os.system(f'kb ref -i {args.index} -g {args.transcripts} -f1 {args.f1} -k 11 --overwrite {args.fasta} {args.gtf}')


    # Call snakemake
    snakefile_path = os.path.join(os.path.dirname(__file__), "Snakefile")
    cmd = [
        "snakemake",
        "--snakefile", snakefile_path,
        "--cores", str(args.cores),
        "--use-conda",
        "--config",
        f"output={args.output}",
        f"index={args.index}",
        f"transcripts={args.transcripts}",
        f"sample1={args.sample1}",
        f"sample2={args.sample2}",
        f"gtf={args.gtf}",
        f"fasta={args.fasta}",
        f"visual={args.visual}",
        f"f1={args.f1}",
        f"reference={args.reference}",
        f"umap={args.umap}"
    ]

    subprocess.run(cmd, check=True)
    end = time.time()
    with open(f"{args.output}/summary.txt", "a") as f:
        f.write(f"Runtime: {end - start:.4f} seconds.")


if __name__ == "__main__":
    main()
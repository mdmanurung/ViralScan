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

        Please note there are 2 different ways to run ViralScan: 
            1. You already have a viral reference index
            2. You don't have a viral reference index and you want ViralScan to create one for you.

        If you select option 1 (so you already have a reference index built with kb-python) then you can run the workflow as follows:
        viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
        All the parameters shown above are mandatory for running ViralScan with option 1.

        If you want option 2 (so you don't have a reference index yet), the workflow can be run as follows:
        viralscan -o output/ -reference True -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz -fasta fasta.fasta -gtf gtf.gtf
        All the parameters shown above are mandatory for running ViralScan with option 2.

        If you want to analyze multiple samples at once or you have multiple gtf and fasta files for creating the reference, 
        these can be added to the parameter which is comma-separated. Please note that there should NOT be a space between
        the file names.

        Example:
        viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample1_1.fastq.gz,sample2_1.fastq.gz -s2 sample1_2.fastq.gz sample2_2.fastq.gz

        Run 'viralscan --help' to show the help function and all other possibilities within the framework. 
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Do NOT create a reference
    # viralscan -t /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/transcriptome/t2g_serratus.txt -i /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/transcriptome/index_serratus.idx -o output_example/ -s1 /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/ebola_samples_test/SRR10307460_1.fastq.gz -s2 /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/ebola_samples_test/SRR10307460_2.fastq.gz 

    # Do create a reference
    # python viralscan_v2/src/viralscan/menu.py -o output -s1 /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/ebola_samples_test/SRR10307460_1.fastq.gz -s2 /exports/archive/hg-funcgenom-research/evonk/data/fasta_viruses/Serratus/ebola_samples_test/SRR10307460_2.fastq.gz -ref True -gtf human_reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf -fasta human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa

    # Required arguments
    parser.add_argument('--output', "-o", required=True, help="The path to the output directory (required)")
    parser.add_argument('--sample1', "-s1", required=True, help="The path to the forward FASTQ sample (gunzipped is preferred) (required)")
    parser.add_argument('--sample2', "-s2", required=True, help="The path to the backward sample (gunzipped is preferred) (required)")

    # Optional arguments
    parser.add_argument('--transcripts', "-t", default=None, help="The path to the transcripts made with kb ref.")
    parser.add_argument('--index', "-i", default=None, help="The path to the reference index created by kb ref.")
    parser.add_argument('--cores', '-c', default=6, type=int, help="The amount of cores the workflow can use. Default: 6.")
    parser.add_argument('--reference', '-ref', type=bool, default=False, help="Should kb ref run? If yes, also provide the path to te gtf, genome, fasta and path to the cDNA fasta files (output). [True/False] (default: False).")
    parser.add_argument('--gtf', '-gtf', default=None, help="Path to GTF files (comma-delimited, without space in-between).")
    parser.add_argument('--fasta', '-fasta', default=None, help="Path to FASTA files (comma-delimited, without space in-between).")
    parser.add_argument('--f1', '-f1', default=None, help="Path to the cDNA FASTA (lamanno, nucleus) or mismatch FASTA (kite) to be generated")
    parser.add_argument('--visual', '-v', default=True, type=bool, help="Add visualizations to the output. [True/False]. Default: True.")
    parser.add_argument('--technology', '-x', default='10xv3', help="Single-cell technology used (`kb --list` to view). Default: 10xv3.")
    parser.add_argument('--whitelist', '-w', default=None, help="Path to file of whitelisted barcodes to correct to. If not " \
                            "provided and bustools supports the technology, a pre-packaged whitelist is used. If not, the " \
                            "bustools whitelist command is used. (`kb --list` to view whitelists)"
                        )
    parser.add_argument('--umap', '-umap', default=False, type=bool, help="Do you want to create a umap? Please note the running time will increase significantly. [True/False]. Default: False.")
    args = parser.parse_args()
    return args


def check_output(args):
    path = args.output
    # check if the path to the output directory exists
    if os.path.isdir(path):
        if os.listdir(path):
            continue_ = input("\033[33mThe output directory already exists and contains files. Do you want to overwrite this? (yes/y/no/n): \033[0m")
            if continue_.lower() == "no" or continue_.lower() == "n":
                print("You have chosen not to continue. The code has been terminated.")
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
                print(f"This is not a valid answer: {continue_}. The code has been terminated.")
                exit()
    else:
        # the output directory does not exist yet, continue as normal
        pass

def errorhandler(args):
    # Error handling
    if not args.reference:
        if not os.path.exists(args.index):
            print(f"\033[31mPath to index does not exist: {args.index}. The code has been terminated.\033[0m")
            exit()
        if not os.path.exists(args.transcripts):
            print(f"\033[31mPath to transcripts does not exist: {args.transcripts}. The code has been terminated.\033[0m")
            exit()  
    # Check sample input    
    if ',' not in args.sample1:
        if not os.path.exists(args.sample1) or not os.path.exists(args.sample2):
            print(f"\033[31mPath to samples is not correct: {args.sample1} or {args.sample2}. The code has been terminated.\033[0m")
            exit()
        valid_sample1 = (args.sample1.endswith(".fastq") or args.sample1.endswith(".fq") or
                 args.sample1.endswith(".fastq.gz") or args.sample1.endswith(".fq.gz"))
        valid_sample2 = (args.sample2.endswith(".fastq") or args.sample2.endswith(".fq") or
                 args.sample2.endswith(".fastq.gz") or args.sample2.endswith(".fq.gz"))
        if not valid_sample1:
            print(f"\033[31mThe forward sample is not in FASTQ format: {args.sample1}. The code has been terminated.\033[0m")
        if not valid_sample2:
            print(f"\033[31mThe backward sample is not in FASTQ format: {args.sample2}. The code has been terminated.\033[0m")

    else:
        samples1 = args.sample1.split(',')
        samples2 = args.sample2.split(',')
        for s1, s2 in zip(samples1, samples2):
            if not os.path.exists(s1) or not os.path.exists(s2):
                print(f"\033[31mPath to samples is not correct: {s1} or {s2}.The code has been terminated.\033[0m")
                exit()
            valid_sample1 = (s1.endswith(".fastq") or s1.endswith(".fq") or
                 s1.endswith(".fastq.gz") or s1.endswith(".fq.gz"))
            valid_sample2 = (s2.endswith(".fastq") or s2.endswith(".fq") or
                    s2.endswith(".fastq.gz") or s2.endswith(".fq.gz"))
            if not valid_sample1:
                print(f"\033[31mThe forward sample is not in FASTQ format: {s1}. The code has been terminated.\033[0m")
            if not valid_sample2:
                print(f"\033[31mThe backward sample is not in FASTQ format: {s2}. The code has been terminated.\033[0m")

    if args.reference:
        if args.transcripts is not None or args.index is not None or args.f1 is not None:
            continue_ = input(f"\033[33mYou have provided a path to the index, or transcripts file but want to create a reference. The reference index will be written to the output directory (-o). Is that what you want? \033[0m")
            if continue_.lower() == "no" or continue_.lower() == "n":
                print("You have chosen not to continue. The code has been aborted.")
                exit()
            elif continue_.lower() == "yes" or continue_.lower() == "y":
                print(f"You have chosen to continue. The given path to the transcipts or index will be ignored.")
            else:
                print(f"This is not a valid answer: {continue_}. The code has been terminated.")
                exit()

        if args.gtf is None:
            print(f"\033[31mThe path to the GTF files is not filled in. If you want to create a reference genome, you need to provide a path to the gtf files. The code has been terminated.\033[0m")
            exit()
        if args.fasta is None:
            print(f"\033[31mThe path to the GTF files is not filled in. If you want to create a reference genome, you need to provide a path to the gtf files. The code has been terminated.\033[0m")
            exit()
        if not os.path.exists(args.gtf):
            print(f"\033[31mThe path to the GTF files is incorrect. If you want to create a reference genome, you need to provide a path to the gtf files. Given path: {args.gtf}. The code has been terminated.\033[0m")
            exit()
        if not os.path.exists(args.fasta):
            print(f"\033[31mThe path to the fasta files is incorrect. If you want to create a reference genome, you need to provide a path to the fasta files. Given path: {args.fasta}. The code has been terminated.\033[0m")
            exit()


    print("\033[32mAll input data has been checked and is correct.\033[0m")


def main():
    start = time.time()
    args = create_help()
    check_output(args)
    errorhandler(args)

    # Check if creating reference is needed and change parameters if necessary
    if args.reference:
        transcripts = f"{args.output}/index/t2g.txt"
        index = f"{args.output}/index/index.idx"
        f1 = f"{args.output}/index/cdna.fa"
        print(f"\033[32mYou have decided a reference still needs to be created. This will be created in the output directory (-o). Please note that creating a reference (depending on the genome) can take a while. Starting kb ref now...\033[0m")
        os.makedirs(f'{args.output}/index/')
        os.system(f'kb ref -i {index} -g {transcripts} -f1 {f1} --overwrite {args.fasta} {args.gtf}')
        print("\033[32mReference index is done!\033[0m")
    else:
        transcripts = args.transcripts
        index = args.index
        f1 = args.f1


    # Call snakemake
    snakefile_path = os.path.join(os.path.dirname(__file__), "Snakefile")
    samples1 = args.sample1.split(',')
    samples2 = args.sample2.split(',')
    if args.output[-1] == "/":
        output = args.output[:-1]
    else:
        output = args.output

    for s1, s2 in zip(samples1, samples2):
        out = f"{s1.split('/')[-1].split('_')[0]}"
        outs = f"{output}/{out}/"
        cmd = [
            "snakemake",
            "--snakefile", snakefile_path,
            "--cores", str(args.cores),
            "--use-conda",
            "--quiet", 
            "all", 
            "--config",
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
            f"whitelist={args.whitelist}"
        ]

        # Start running the snakefile
        subprocess.run(cmd, check=True)

        # Write last things to summary file
        end = time.time()
        with open(f"{outs}/summary.txt", "a") as f:
            f.write(f"\nRuntime: {end - start:.4f} seconds.\n\n")
            transcripts_ = subprocess.check_output(
                "wc -l {}".format(transcripts),
                shell=True, text=True
            ).strip()
            transcripts_ = transcripts_.split()[0]
            genes = subprocess.check_output(
            "cut -f2,3 {} | sort -u | wc -l".format(transcripts),
            shell=True, text=True
            ).strip()
            f.write(f"Amount of transcripts in data: {transcripts_}\n")
            f.write(f"Amount of genes in data: {genes}\n")
            # f.write(f"Number of barcodes: {}")
        
        # Unlock snakemake file
        unlock_cmd = [
            "snakemake",
            "--snakefile", snakefile_path,
            "--unlock",
            "--config",
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
            f"umap={args.umap}"
        ]
        subprocess.run(unlock_cmd, check=True)


if __name__ == "__main__":
    main()
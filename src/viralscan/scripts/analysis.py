# Import packages
import glob
import yaml
from pathlib import Path
import os

configfile = snakemake.params.configfile
with open(configfile, 'r') as f:
    config = yaml.safe_load(f)

def obtain_gtf():
    # Define user input out of CONFIG file.
    viral_accessions = set()

    # get the standard gtf files from this package
    project_root = Path(__file__).resolve().parent.parent

    # Point to folder and get the gtf files
    data_dir = project_root / "data"
    gtf_files = list(data_dir.glob("*.gtf"))

    # Serratus viruses
    for file in gtf_files:
        with open(file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    info = line.split("\t")[8]
                    gene_id = info.split('"')[1]
                    viral_accessions.add(gene_id)
                    
    # check if GTF has been added by user. If so, add them to the viral list
    if config['gtf'] != "None":
        if os.path.exists(config["gtf"]):
            gtf_files = config["gtf"].split(',')
            for file in gtf_files:
                with open(file, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            info = line.split("\t")[8]
                            gene_id = info.split('"')[1]
                            viral_accessions.add(gene_id)

    # write list to file
    with open(f"{config['output']}log/analysis.txt", 'w') as f:
        for v in viral_accessions:
            f.write(v+"\n")
    f.close()
    return viral_accessions


def main():
    obtain_gtf()

main()

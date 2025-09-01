# ViralScan

## Introduction
**ViralScan** is a computational bioinformatics framework to detect viral load in peripheral blood samples. It is designed for the Leiden University Medical Centre (LUMC) to enable the detecttion of (possibly) multiple viruses in one sample.

The tool can use both kb ref and kb count to perform scalable detection of the viral load.

---

## Installation
ViralScan can be installed using `pip` in a Python environment:
```
pip install viralscan (**not working yet**)
```

A conda environment is recommended. You can create one as follows:
```
conda create -n bioenv -c conda-forge -c bioconda snakemake kb-python
```

ViralScan is command-line based. To view available commands and options, run the following command:
```
viralscan --help
```
This will display all available options regarding the software, including an example.

---
## Input Data
### Reference Index
ViralScan requires an index created with kb ref (from thhe kb-python package).
There are 2 options regarding the index:
- provide your own pre-built index, or;
- let ViralScan create its own index for you.

Reference data to create the index can be found on Zenodo (https://zenodo.org/uploads/16792022).

### Samples
ViralScan expects paired-end FASTQ files (potentially gunzipped). Both sample files (forward and backward) must be provided for the analysis.

---
## User Guide
Please Note: all the output from ViralScan (including the logs and plots) will be created in the output folder defined by the user. To run ViralScan as is (without any extra's), run the following:
```
ViralScan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
```

Important: if you decide to create a reference index, make sure the path to the reference index does NOT exist yet. Otherwise it will overwrite the existing index.

If no reference index has been created uet and you want ViralScan to do it for you, run this:
```
# If you have 1 gtf and 1 fasta file
! viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz -ref True -gtf path_to_gtf/gtf.gtf -fasta path_to_fasta/fasta.fasta

# If you have multiple gtfs and fasta files (comma-separated without space)
! Viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz -ref True -gtf path_to_gtf/gtf1.gtf,path_to_gtf/gtf2.gtf -fasta path_to_fasta/fasta1.fasta,path_to_fasta/fasta2.fasta
```

If you only want to create a reference index (with kb ref) and don't want to run ViralScan (so you don't want to quantify the data), do this:
```
# If you have 1 gtf and 1 fasta file
! kb ref -i index.idx -g transcripts.txt -f1 cdna.fasta --overwrite fasta.fasta gtf.gtf

# If you have multiple gtfs and fasta files (comma-separated without space)
! kb ref -i index.idx -g transcripts.txt -f1 cdna.fasta --overwrite fasta1.fasta,fasta2.fasta gtf1.gtf,gtf2.gtf
```

## License
This project is licensed under the MIT License. See the LICENSE file for details.
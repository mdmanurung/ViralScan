# ViralScan

## Introduction
**ViralScan** is a computational bioinformatics framework to detect viral load in samples. It is designed for the Leiden University Medical Centre (LUMC) to enable the detection of multiple viruses.

The tool can use both kb ref and kb count to perform scalable detection of the viral load.

---

## Installation
ViralScan can be installed using `pip` in a Python environment:
```
pip install -i https://test.pypi.org/simple/ ViralScan
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
ViralScan requires an index created with kb ref (from the kb-python package).
There are 2 options regarding the index:
- provide your own pre-built index, or;
- let ViralScan create its own index for you.

Reference data to create the index can be found on Zenodo (https://zenodo.org/uploads/16792022).

### Samples
ViralScan expects paired-end FASTQ files (potentially gunzipped). Both sample files (forward and backward) must be provided for the analysis.
Samples which can be used for testing: SRR20710651, SRR20710645 and SRR10307460

---
## User Guide
There are 3 ways to run ViralScan:
1. You already have a reference index built with kb-python
2. You don't have a reference index and want ViralScan to build one from FASTA + GTF files you provide
3. You don't have a reference index and want ViralScan to download FASTA + GTF for one or more NCBI accessions and build the index for you

Please Note: all the output from ViralScan (including the logs and plots) will be created in the output folder defined by the user. To run ViralScan with option 1, run the following:
```
viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
```

If you want ViralScan to create a reference index for you, you have to provide the GTF and FASTA file(s). Running ViralScan to create a reference index and to perform the quantification can be done as follows:
```
viralscan -o output/ --reference -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz -fasta fasta.fasta -gtf gtf.gtf
```
The index will be placed in the output directory, in a subfolder called `index`.

If you only have NCBI accession numbers (e.g. RefSeq IDs like `NC_002021.3`), ViralScan can fetch the FASTA + GTF for you and build the index. Provide one or more accessions, comma-separated:
```
viralscan -o output/ -acc NC_002021.3 -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
viralscan -o output/ -acc NC_002021.3,NC_001512.1 -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
```
NCBI requires an email address. Either pass `--ncbi-email you@example.org` or set the `NCBI_EMAIL` environment variable. An optional `NCBI_API_KEY` env var is honoured for higher rate limits. Downloaded references are cached under `~/.cache/viralscan/ncbi/` so re-runs do not re-download.

If you have multiple samples, ViralScan can analyze them with 1 command. Just split the names of the samples with a comma (without a space in-between). The same is when you have multiple GTF and FASTA files. For example:
```
# If you have multiple samples
viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample1_1.fastq.gz,sample2_1.fastq.gz -s2 sample1_2.fastq.gz,sample2_2.fastq.gz

# If you have multiple gtfs and fasta files
viralscan -o output/ --reference -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz -fasta fasta1.fasta,fasta2.fasta -gtf gtf1.gtf,gtf2.gtf
```

For information about other parameters or possibilities in ViralScan, call the help function:
```
viralscan --help
```

## License
This project is licensed under the MIT License. See the LICENSE file for details.

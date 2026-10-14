OUT = config["output"]
SCRIPT = config["assembly_script"]
PYTHON = config["python"]
SAMPLES = [*config["samples"], "fixture"]

wildcard_constraints:
    sample="|".join(SAMPLES)

rule all:
    input:
        expand(OUT + "/{sample}/assembly/assembly_summary.json", sample=SAMPLES)

rule reference:
    input:
        script=SCRIPT,
        cfg=config["assembly_config"],
        candidates=lambda w: [] if w.sample == "fixture" else [OUT + "/" + w.sample + "/off/candidate_reads.fastq"],
        fasta=config["viral_fasta"],
        gtf=config["viral_gtf"]
    output:
        fasta=OUT + "/{sample}/assembly/reference.fa",
        gtf=OUT + "/{sample}/assembly/reference.gtf",
        reads=OUT + "/{sample}/assembly/reads.fastq",
        receipt=OUT + "/{sample}/assembly/reference.json"
    shell:
        "{PYTHON:q} {SCRIPT:q} reference {input.cfg:q} {wildcards.sample:q}"

rule index:
    input:
        script=SCRIPT,
        cfg=config["assembly_config"],
        fasta=OUT + "/{sample}/assembly/reference.fa",
        gtf=OUT + "/{sample}/assembly/reference.gtf"
    output:
        OUT + "/{sample}/assembly/star_index/Genome"
    threads: 2
    params:
        cfg=config["assembly_config"]
    shell:
        "{PYTHON:q} {SCRIPT:q} index {params.cfg:q} {wildcards.sample:q}"

rule align:
    input:
        script=SCRIPT,
        cfg=config["assembly_config"],
        index=OUT + "/{sample}/assembly/star_index/Genome",
        reads=OUT + "/{sample}/assembly/reads.fastq"
    output:
        OUT + "/{sample}/assembly/Aligned.sortedByCoord.out.bam"
    threads: 2
    params:
        cfg=config["assembly_config"]
    shell:
        "{PYTHON:q} {SCRIPT:q} align {params.cfg:q} {wildcards.sample:q}"

rule assemble:
    input:
        script=SCRIPT,
        cfg=config["assembly_config"],
        bam=OUT + "/{sample}/assembly/Aligned.sortedByCoord.out.bam"
    output:
        gtf=OUT + "/{sample}/assembly/transcripts.gtf",
        receipt=OUT + "/{sample}/assembly/assembly_summary.json"
    threads: 2
    params:
        cfg=config["assembly_config"]
    shell:
        "{PYTHON:q} {SCRIPT:q} assemble {params.cfg:q} {wildcards.sample:q}"

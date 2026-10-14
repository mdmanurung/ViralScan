"""Optional reference-based reconstruction from the pilot's selected candidate reads."""

from __future__ import annotations

import argparse
import json
import math
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from viralscan.anello_align import iter_fasta
from viralscan.run_safety import sha256_file


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def run(command, log):
    with open(log, "w") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)


def reference(config_path, sample):
    cfg = json.loads(Path(config_path).read_text())
    out = Path(cfg["output"]) / sample / "assembly"
    out.mkdir(parents=True, exist_ok=True)
    if sample == "fixture":
        rng = random.Random(20261014)
        seq = "".join(rng.choice("ACGT") for _ in range(5000))
        seq = seq[:800] + "GT" + seq[802:1198] + "AG" + seq[1200:]
        (out / "reference.fa").write_text(">fixture\n" + seq + "\n")
        exons = [(501, 800, "spliced"), (1201, 1500, "spliced"), (3001, 3400, "single")]
        with open(out / "reference.gtf", "w") as f:
            for start, end, gene in exons:
                f.write(
                    f'fixture\tpilot\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{gene}.1";\n'
                )
        with open(out / "reads.fastq", "w") as f:
            for gene, transcript in [
                ("spliced", seq[500:800] + seq[1200:1500]),
                ("single", seq[3000:3400]),
            ]:
                for start in range(0, len(transcript) - 99, 5):
                    read = transcript[start : start + 100]
                    f.write(f"@fixture|{gene}_{start}\n{read}\n+\n{'I' * len(read)}\n")
    else:
        accessions = cfg["samples"][sample]["accessions"]
        found = {}
        with open(out / "reference.fa", "w") as f:
            for name, seq in iter_fasta(Path(cfg["viral_fasta"])):
                if name in accessions:
                    found[name] = len(seq)
                    f.write(f">{name}\n{seq}\n")
        if set(found) != set(accessions):
            raise ValueError("Missing exact accession in assembly reference")
        with open(out / "reference.gtf", "w") as f:
            for line in Path(cfg["viral_gtf"]).read_text().splitlines():
                if line.split("\t")[0] in accessions:
                    f.write(line + "\n")
        reads = Path(cfg["output"]) / sample / "off/candidate_reads.fastq"
        (out / "reads.fastq").symlink_to(reads)
    write_json(
        out / "reference.json",
        {
            "sample": sample,
            "reference_sha256": sha256_file(out / "reference.fa"),
            "gtf_sha256": sha256_file(out / "reference.gtf"),
            "reads_sha256": sha256_file(out / "reads.fastq"),
            "selection": "fixture truth"
            if sample == "fixture"
            else "combined-index unique-target candidate reads; same detection data",
            "independent_confirmation": False,
            "barcode_correction": False,
            "full_host_genome_competition": False,
        },
    )


def index(config_path, sample):
    cfg = json.loads(Path(config_path).read_text())
    out = Path(cfg["output"]) / sample / "assembly"
    size = sum(len(seq) for _, seq in iter_fasta(out / "reference.fa"))
    index_dir = out / "star_index"
    index_dir.mkdir(exist_ok=True)
    command = [
        cfg["star"],
        "--runMode",
        "genomeGenerate",
        "--runThreadN",
        "2",
        "--genomeDir",
        str(index_dir),
        "--genomeFastaFiles",
        str(out / "reference.fa"),
        "--sjdbGTFfile",
        str(out / "reference.gtf"),
        "--sjdbOverhang",
        "99",
        "--genomeSAindexNbases",
        str(min(7, max(1, int(math.log2(size) / 2) - 1))),
        "--genomeChrBinNbits",
        "12",
        "--outFileNamePrefix",
        str(out / "index_"),
    ]
    write_json(out / "index_command.json", command)
    run(command, out / "index.log")


def align(config_path, sample):
    cfg = json.loads(Path(config_path).read_text())
    out = Path(cfg["output"]) / sample / "assembly"
    bam = out / "Aligned.sortedByCoord.out.bam"
    if (out / "reads.fastq").stat().st_size == 0:
        header = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(
            f"@SQ\tSN:{name}\tLN:{len(seq)}\n" for name, seq in iter_fasta(out / "reference.fa")
        )
        subprocess.run(
            [cfg["samtools"], "view", "-b", "-o", str(bam), "-"], input=header.encode(), check=True
        )
        write_json(out / "align_command.json", {"status": "no_candidate_reads"})
        return
    command = [
        cfg["star"],
        "--runThreadN",
        "2",
        "--genomeDir",
        str(out / "star_index"),
        "--readFilesIn",
        str(out / "reads.fastq"),
        "--outFileNamePrefix",
        str(out) + "/",
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
        "--outSAMstrandField",
        "intronMotif",
        "--outSAMattributes",
        "NH",
        "HI",
        "AS",
        "nM",
        "XS",
        "--alignIntronMax",
        "200000",
        "--outFilterMatchNmin",
        "31",
        "--outFilterScoreMinOverLread",
        "0.6",
        "--outFilterMatchNminOverLread",
        "0.6",
        "--limitBAMsortRAM",
        "1000000000",
    ]
    write_json(out / "align_command.json", command)
    run(command, out / "align.log")
    subprocess.run([cfg["samtools"], "index", str(bam)], check=True)


def assemble(config_path, sample):
    cfg = json.loads(Path(config_path).read_text())
    out = Path(cfg["output"]) / sample / "assembly"
    gtf = out / "transcripts.gtf"
    if (out / "reads.fastq").stat().st_size == 0:
        gtf.write_text("# no candidate reads; assembly not assessable\n")
    else:
        command = [
            cfg["stringtie"],
            str(out / "Aligned.sortedByCoord.out.bam"),
            "-p",
            "2",
            "-G",
            str(out / "reference.gtf"),
            "-o",
            str(gtf),
            "-l",
            sample,
            "-f",
            "0.01",
            "-g",
            "5",
            "--fr",
        ]
        write_json(out / "assembly_command.json", command)
        run(command, out / "assembly.log")
    features = [
        line.split("\t")
        for line in gtf.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    for row in features:
        if len(row) != 9 or int(row[3]) < 1 or int(row[4]) < int(row[3]):
            raise ValueError("Malformed transcript model")
    if sample == "fixture":
        intervals = {(int(row[3]), int(row[4])) for row in features if row[2] == "exon"}
        if not {(501, 800), (1201, 1500), (3001, 3400)} <= intervals:
            raise ValueError("Known spliced/single-exon fixture reconstruction failed")
        junctions = (out / "SJ.out.tab").read_text().splitlines()
        if not any(
            row.split("\t")[:3] == ["fixture", "801", "1200"] and int(row.split("\t")[6]) > 0
            for row in junctions
        ):
            raise ValueError("Fixture has no independent spliced-read junction support")
    n = sum(row[2] == "transcript" for row in features)
    write_json(
        out / "assembly_summary.json",
        {
            "sample": sample,
            "transcript_models": n,
            "exons": sum(row[2] == "exon" for row in features),
            "status": "reconstructed_candidate_fragments" if n else "not_assessable",
            "fixture_truth_verified": sample == "fixture",
            "limitation": "reference-based and detection-selected; no full-transcript, full-genome or infection confirmation claim",
        },
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("step", choices=["reference", "index", "align", "assemble"])
    p.add_argument("config")
    p.add_argument("sample")
    a = p.parse_args()
    globals()[a.step](a.config, a.sample)


if __name__ == "__main__":
    main()

"""STARsolo host-unmapped alignment against a combined HOST|/VIRUS| index (STAR-HOSTVIRUS-01)."""

import random
import shutil
import subprocess

import pytest

from viralscan.star_hostvirus import (
    align_cmd,
    align_host_unmapped,
    build_star_hostvirus_index,
    write_reference,
)


def _seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_write_reference_prefixes(tmp_path):
    (tmp_path / "h.fa").write_text(">chr1 x\nACGT\n")
    (tmp_path / "v.fa").write_text(">AB1.1\nACGT\n")
    (tmp_path / "h.gtf").write_text('chr1\ts\texon\t1\t4\t.\t+\t.\tgene_id "g";\n')
    (tmp_path / "v.gtf").write_text('AB1.1\ts\texon\t1\t4\t.\t+\t.\tgene_id "v";\n')
    fa, gtf = write_reference(
        tmp_path / "h.fa", tmp_path / "v.fa", tmp_path / "h.gtf", tmp_path / "v.gtf", tmp_path / "o"
    )
    assert [x for x in fa.read_text().splitlines() if x[0] == ">"] == [
        ">HOST|chr1 x",
        ">VIRUS|AB1.1",
    ]
    assert [x.split("\t")[0] for x in gtf.read_text().splitlines()] == ["HOST|chr1", "VIRUS|AB1.1"]


def test_align_cmd_geometry_and_multimappers():
    cmd = align_cmd("STAR", "idx", "r1.fq.gz", "r2.fq.gz", "10xv2", "out/", 4)
    assert cmd[cmd.index("--readFilesIn") + 1 : cmd.index("--readFilesIn") + 3] == [
        "r2.fq.gz",
        "r1.fq.gz",
    ]
    assert cmd[cmd.index("--soloCBlen") + 1] == "16" and cmd[cmd.index("--soloUMIlen") + 1] == "10"
    assert cmd[cmd.index("--outFilterMultimapNmax") + 1] == "100"
    assert cmd[cmd.index("--soloCBwhitelist") + 1] == "None"


@pytest.mark.skipif(
    shutil.which("STAR") is None or shutil.which("samtools") is None, reason="needs STAR, samtools"
)
def test_end_to_end_viral_and_host_reads(tmp_path):
    host, virus = _seq(20000, 1), _seq(3000, 2)
    (tmp_path / "h.fa").write_text(f">chr1\n{host}\n")
    (tmp_path / "v.fa").write_text(f">V1.1\n{virus}\n")
    (tmp_path / "h.gtf").write_text(
        'chr1\ts\texon\t100\t900\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    )
    (tmp_path / "v.gtf").write_text(
        'V1.1\ts\texon\t100\t900\t.\t+\t.\tgene_id "v"; transcript_id "tv";\n'
    )
    fa, gtf = write_reference(
        tmp_path / "h.fa",
        tmp_path / "v.fa",
        tmp_path / "h.gtf",
        tmp_path / "v.gtf",
        tmp_path / "ref",
    )
    index = build_star_hostvirus_index(fa, gtf, tmp_path / "idx", threads=1, sa_index_nbases=6)
    cb, umi = "ACGTACGTACGTACGT", "TTGCAATTGC"
    r1, r2 = [], []
    for i, src in enumerate([virus] * 3 + [host] * 2):
        pos = 200 + 10 * i
        r2.append(f"@r{i}\n{src[pos : pos + 90]}\n+\n{'I' * 90}\n")
        r1.append(f"@r{i}\n{cb}{umi}\n+\n{'I' * 26}\n")
    (tmp_path / "R1.fq").write_text("".join(r1))
    (tmp_path / "R2.fq").write_text("".join(r2))
    bam = align_host_unmapped(
        tmp_path / "R1.fq",
        tmp_path / "R2.fq",
        index,
        tmp_path / "aln",
        "10xv2",
        whitelist=None,
        threads=1,
    )
    assert bam.with_name(bam.name + ".bai").exists()
    sam = subprocess.run(["samtools", "view", str(bam)], capture_output=True, text=True, check=True)
    rows = [line.split("\t") for line in sam.stdout.splitlines()]
    assert sorted(r[2] for r in rows) == ["HOST|chr1"] * 2 + ["VIRUS|V1.1"] * 3
    assert all(f"CB:Z:{cb}" in r and f"UB:Z:{umi}" in r for r in rows)

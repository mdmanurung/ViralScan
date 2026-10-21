"""scPathoQuant-style coverage maps (pyGenomeTracks) from a live samtools BAM."""

import shutil
import subprocess

import pytest

from viralscan.evidence import plot_coverage_maps

TOOLS = ("samtools", "bamCoverage", "pyGenomeTracks")


@pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS), reason=f"needs {TOOLS}")
def test_one_map_per_viral_reference_and_layer(tmp_path):
    sam = tmp_path / "reads.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:VIRUS|X1\tLN:1000\n@SQ\tSN:VIRUS|X2\tLN:500\n@SQ\tSN:chr1\tLN:5000\n"
        + "".join(
            f"r{i}\t0\tVIRUS|X1\t{p}\t60\t1M\t*\t0\t0\tA\tI{sa}\n"
            for i, (p, sa) in enumerate(
                [
                    (5, ""),
                    (5, "\tSA:Z:HOST|chr1,10,+,1M,60,0;"),
                    (900, "\tSA:Z:VIRUS|X2,5,+,1M,60,0;"),
                ]
            )
        )
        + "r9\t0\tchr1\t10\t60\t1M\t*\t0\t0\tA\tI\n"
    )
    bam = tmp_path / "reads.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    gtf = tmp_path / "x.gtf"
    gtf.write_text('X1\tt\tgene\t100\t400\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n')
    written = plot_coverage_maps({"raw": str(bam), "deduplicated": str(bam)}, tmp_path, [gtf])
    names = sorted(p.rsplit("/", 1)[1] for p in written)
    kinds = ("coveragemap", "coveragemap_competition", "coverage_stacked")
    assert names == sorted(
        f"VIRUS_X1_{k}.{layer}.{ext}"
        for k in kinds
        for layer in ("deduplicated", "raw")
        for ext in ("pdf", "png")
    )
    assert all((tmp_path / n).stat().st_size > 0 for n in names)


def test_read_category_from_sa_tag():
    from viralscan.evidence import _read_category

    assert _read_category(None) == "unique"
    assert _read_category("HOST|chr1,1,+,1M,60,0;") == "host_competing"
    assert _read_category("VIRUS|X2,1,+,1M,60,0;") == "other_virus"
    assert _read_category("VIRUS|X2,1,+,1M,60,0;HOST|chr1,1,+,1M,60,0;") == "host_competing"

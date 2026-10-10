"""Live explicit-reference competition against toy local tool inputs."""

import json
import random

import pytest

from viralscan import evidence as ev

pytestmark = pytest.mark.integration


def test_live_multiclass_blast_retains_more_than_twenty_ties_and_no_hit(tmp_path):
    missing = ev.have_tools(["blastn", "makeblastdb"])
    if missing:
        pytest.skip(f"missing tools: {missing}")
    rng = random.Random(449)
    sequence = "".join(rng.choice("ACGT") for _ in range(600))
    refs = tmp_path / "refs.fa"
    refs.write_text(
        "".join(f">TARGET|t{i}\n{sequence}\n" for i in range(25))
        + ">HOST|h\n"
        + sequence[::-1]
        + "\n"
    )
    reads = tmp_path / "reads.fa"
    reads.write_text(f">cell_u1_0|r0\n{sequence[120:270]}\n>cell_u2_1|r1\n{'N' * 150}\n")
    sampling = tmp_path / "sampling.json"
    rows = ev.competitive_blast_identity(
        str(reads),
        str(refs),
        str(tmp_path / "blast"),
        threads=1,
        multi_class=True,
        sampling_manifest=str(sampling),
    )
    match = next(r for r in rows if r["read"] == "cell_u1_0|r0")
    absent = next(r for r in rows if r["read"] == "cell_u2_1|r1")
    assert match["n_top_tied_hits"] == 25
    assert match["diagnostic_class"] == "target_specific"
    assert absent["diagnostic_class"] == "no_hits"
    assert len((tmp_path / "blast/blast_raw_hits.tsv").read_text().splitlines()) >= 25
    assert json.loads(sampling.read_text())["max_target_seqs"] == 26
    bam = ev.align_reads_to_viral(str(reads), str(refs), str(tmp_path / "raw.bam"), threads=1)
    classes = {r["reference_class"] for r in ev.alignment_qc_table(bam)}
    assert "target" in classes


def test_live_target_related_host_classes_survive_qc_and_plots(tmp_path):
    missing = ev.have_tools(["minimap2", "samtools", "blastn", "makeblastdb"])
    if missing:
        pytest.skip(f"missing tools: {missing}")
    rng = random.Random(661)
    target = "".join(rng.choice("ACGT") for _ in range(800))
    related_chars = list(target)
    for i in range(0, 800, 23):
        related_chars[i] = "ACGT"[("ACGT".index(related_chars[i]) + 1) % 4]
    related = "".join(related_chars)
    host = "".join(rng.choice("ACGT") for _ in range(800))
    references = tmp_path / "refs.fa"
    references.write_text(f">TARGET|t\n{target}\n>RELATED|r\n{related}\n>HOST|h\n{host}\n")
    reads = tmp_path / "reads.fa"
    reads.write_text(
        f">c_u1_0|r0\n{target[150:350]}\n>c_u2_1|r1\n{related[150:350]}\n>c_u3_2|r2\n{host[150:350]}\n"
    )
    bam = ev.align_reads_to_viral(str(reads), str(references), str(tmp_path / "raw.bam"), threads=1)
    assert {r["reference_class"] for r in ev.alignment_qc_table(bam)} == {
        "target",
        "related_virus",
        "host",
    }
    assert {r["reference_class"] for r in ev.per_cell_alignment_qc(bam)} == {
        "target",
        "related_virus",
        "host",
    }
    rows = ev.competitive_blast_identity(
        str(reads), str(references), str(tmp_path / "blast"), threads=1, multi_class=True
    )
    assert {r["diagnostic_class"] for r in rows} == {
        "target_specific",
        "related_preferred",
        "host_preferred",
    }
    plot = tmp_path / "coverage.png"
    ev.plot_coverage_comparison(bam, bam, str(plot))
    assert plot.stat().st_size > 0

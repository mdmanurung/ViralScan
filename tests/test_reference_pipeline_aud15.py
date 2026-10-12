"""Offline MECH-E artifact invariance; kb indexing itself is mocked."""

import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from viralscan.scripts import build_reference as builder

VIRAL_FASTA = b">FIXTURE_A\nNNNNNNNN\n>FIXTURE_B\nNNNNNNNNNN\n"
HOST_FASTA = b">HOST_TX gene:HOST_GENE\nNNNN\n"
NCBI_GTF = 'FIXTURE_A\tfixture\texon\t2\t7\t.\t+\t.\tgene_id "A_GENE"; transcript_id "A_TX";\n'
PLACEHOLDER_ATTRS = (
    'gene_id "FIXTURE_B_gene1"; transcript_id "FIXTURE_B_tx1"; '
    'gene_name "FIXTURE_B"; gene_biotype "whole_genome";'
)
PLACEHOLDER_GTF = "".join(
    f"FIXTURE_B\tViralScan\t{feature}\t1\t10\t.\t+\t.\t{PLACEHOLDER_ATTRS}\n"
    for feature in ("gene", "transcript", "exon")
)
HOST_GTF = "".join(
    f"HOST_TX\tEnsembl_cDNA\t{feature}\t1\t4\t.\t+\t.\t"
    'gene_id "HOST_GENE"; transcript_id "HOST_TX";\n'
    for feature in ("gene", "transcript", "exon")
)
CATALOGUE_GTF = (
    'FIXTURE_A\tNCBI\texon\t2\t7\t.\t+\t0\tgene_id "A_GENE"; '
    'transcript_id "A_TX"; gene_name "fixture"; gene_biotype "protein_coding"; '
    'product "fixture"; n_exons "1";\n'
)


@pytest.fixture
def offline_inputs(tmp_path, monkeypatch):
    viral = tmp_path / "input.fa"
    viral.write_bytes(VIRAL_FASTA)
    viral_gtf = tmp_path / "input.gtf"
    viral_gtf.write_text(NCBI_GTF)
    host = tmp_path / "host.fa.gz"
    with gzip.open(host, "wb") as handle:
        handle.write(HOST_FASTA)
    host_gtf = tmp_path / "host.gtf.gz"
    with gzip.open(host_gtf, "wb") as handle:
        handle.write(b"unused chromosomal companion annotation\n")
    monkeypatch.setattr(builder, "fetch_host_cdna", lambda *a, **k: (host, host_gtf))
    monkeypatch.setattr(
        "viralscan.scripts.ncbi_fetch.fetch_reference", lambda *a, **k: (viral, viral_gtf)
    )
    monkeypatch.setattr(
        "viralscan.anellovirus.load_gene_table",
        lambda *a, **k: [
            {
                "accession": "FIXTURE_A",
                "gene_id": "A_GENE",
                "transcript_id": "A_TX",
                "product": "fixture",
                "gene_symbol": "",
                "gene": "fixture",
                "n_exons": 1,
                "exons": "2:7",
                "strand": "+",
            }
        ],
    )
    monkeypatch.setattr(
        builder, "mask_low_complexity", lambda src, dst: shutil.copyfile(src, dst) or True
    )
    monkeypatch.setattr(builder.shutil, "which", lambda name: "/mock/kb" if name == "kb" else None)
    monkeypatch.setattr(builder, "reference_tool_versions", lambda: {})
    commands = []

    def fake_kb(cmd, out_dir):
        commands.append(cmd)
        Path(cmd[cmd.index("-i") + 1]).write_bytes(b"mock index")
        # Model only kb's transcript-to-gene mapping boundary, not its sequence splitter.
        pairs = []
        for line in Path(cmd[-1]).read_text().splitlines():
            gene = re.search(r'gene_id "([^"]+)"', line)
            tx = re.search(r'transcript_id "([^"]+)"', line)
            if gene and tx and (tx.group(1), gene.group(1)) not in pairs:
                pairs.append((tx.group(1), gene.group(1)))
        Path(cmd[cmd.index("-g") + 1]).write_text("".join(f"{tx}\t{gene}\n" for tx, gene in pairs))
        return {"elapsed_seconds": 0.0}

    monkeypatch.setattr(builder, "_run_kb_ref", fake_kb)
    return viral, commands


@pytest.mark.parametrize("source", ["combined", "anellovirus-local", "anellovirus-ncbi"])
@pytest.mark.parametrize("mask", [False, True])
@pytest.mark.parametrize("index", [False, True])
@pytest.mark.parametrize("with_dlist", [False, True])
def test_builder_artifacts_match_pre_refactor_golden(
    tmp_path, monkeypatch, offline_inputs, source, mask, index, with_dlist
):
    """Literal outputs frozen before the orchestration refactor, including final newlines."""
    viral, commands = offline_inputs
    dlist = None
    if with_dlist:
        dlist = tmp_path / "genome.fa"
        dlist.write_text(">host\nNNNN\n")
        monkeypatch.setattr(builder, "measure_host_homology", lambda *a, **k: {})
    if source == "combined":
        result = builder.build_combined_reference(
            "human",
            ["FIXTURE_A", "FIXTURE_B"],
            tmp_path / "out",
            mask=mask,
            run_kb_ref=index,
            genome_dlist=dlist,
        )
        want_fasta = HOST_FASTA + VIRAL_FASTA
        want_gtf = HOST_GTF + NCBI_GTF + PLACEHOLDER_GTF
        want_t2g = b"HOST_TX\tHOST_GENE\nA_TX\tA_GENE\nFIXTURE_B_tx1\tFIXTURE_B_gene1\n"
    else:
        result = builder.build_anellovirus_reference(
            tmp_path / "out",
            accessions=["FIXTURE_A", "FIXTURE_B"],
            mask=mask,
            run_kb_ref=index,
            fasta_path=viral if source.endswith("local") else None,
            genome_dlist=dlist,
        )
        want_fasta = VIRAL_FASTA
        want_gtf = CATALOGUE_GTF + PLACEHOLDER_GTF
        want_t2g = b"A_TX\tA_GENE\nFIXTURE_B_tx1\tFIXTURE_B_gene1\n"
    assert result["fasta"].read_bytes() == want_fasta
    assert result["gtf"].read_bytes() == want_gtf.encode()
    if index:
        assert result["t2g"].read_bytes() == want_t2g
        assert len(commands) == 1
        assert commands[0][-2:] == [str(result["fasta"]), str(result["gtf"])]
        assert ("--d-list" in commands[0]) is with_dlist
    else:
        assert result["index"] is result["t2g"] is None
        assert commands == []
    manifest = json.loads(result["manifest"].read_text())
    assert manifest["content_files"]["fasta"] == builder.sha256_file(result["fasta"])
    assert manifest["content_files"]["gtf"] == builder.sha256_file(result["gtf"])
    assert manifest["build_receipt"]["genome_dlist_applied"] is (index and with_dlist)
    assert manifest["source_adapters"] == (
        ["ensembl_cdna", "ncbi_genbank"]
        if source == "combined"
        else ["local_fasta" if source.endswith("local") else "ncbi_genbank"]
    )


@pytest.mark.parametrize("adapter", ["combined", "catalogue", "placeholder"])
def test_shared_gtf_emitter_preserves_mixed_source_order(tmp_path, offline_inputs, adapter):
    """Keep literal pre-refactor separators and catalogue-before-placeholder ordering."""
    fasta = tmp_path / "mixed.fa"
    fasta.write_bytes(b">FIXTURE_B description\nNNNNN\nNNNNN\n>FIXTURE_A\nNNNNNNNN\n")
    gtf = tmp_path / "mixed.gtf"
    if adapter == "catalogue":
        assert builder._anellovirus_gtf(fasta, gtf) == (1, 1)
        expected = CATALOGUE_GTF + PLACEHOLDER_GTF
    elif adapter == "combined":
        sources, placeholders = builder._emit_viral_gtf(
            fasta,
            gtf,
            anello_accessions={"FIXTURE_A"},
            real_gtf_blocks=builder.index_gtf_by_seqname(NCBI_GTF),
        )
        assert sources == {"catalogue": 1, "placeholder": 1}
        assert placeholders == ["FIXTURE_B"]
        expected = PLACEHOLDER_GTF + CATALOGUE_GTF + "\n"
    else:
        builder._gtf_from_merged_fasta(fasta, gtf)
        expected = PLACEHOLDER_GTF + PLACEHOLDER_GTF.replace("FIXTURE_B", "FIXTURE_A").replace(
            "\t10\t", "\t8\t"
        )
    assert gtf.read_bytes() == expected.encode()


@pytest.mark.parametrize("source", ["combined", "anellovirus-local", "anellovirus-ncbi"])
@pytest.mark.parametrize("failure", ["missing-kb", "kb-fails"])
def test_shared_index_boundary_never_records_a_failed_build(
    tmp_path, monkeypatch, offline_inputs, source, failure
):
    viral, _ = offline_inputs
    if failure == "missing-kb":
        monkeypatch.setattr(builder.shutil, "which", lambda name: None)
        expected_error = RuntimeError
    else:

        def fail_kb(command, out_dir):
            raise subprocess.CalledProcessError(1, command)

        monkeypatch.setattr(builder, "_run_kb_ref", fail_kb)
        expected_error = subprocess.CalledProcessError
    out = tmp_path / "out"
    with pytest.raises(expected_error):
        if source == "combined":
            builder.build_combined_reference("human", ["FIXTURE_A", "FIXTURE_B"], out, mask=False)
        else:
            builder.build_anellovirus_reference(
                out,
                accessions=["FIXTURE_A", "FIXTURE_B"],
                mask=False,
                fasta_path=viral if source.endswith("local") else None,
            )
    assert not (out / "reference_reproducibility.json").exists()
    manifest = json.loads((out / "reference_manifest.json").read_text())
    assert "command" not in manifest["build_receipt"]

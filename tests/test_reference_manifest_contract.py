"""Frozen-input reference contracts, independent of kb and remote services."""

import json
import subprocess
from pathlib import Path

import pytest

from viralscan.scripts import build_reference as builder
from viralscan.scripts import ncbi_fetch


def test_local_input_does_not_fabricate_retrieval_or_database(tmp_path):
    fasta = tmp_path / "local.fa"
    fasta.write_text(">V1\nACGTCAGTAC\n")
    manifest = builder.write_reference_manifest(
        fasta,
        tmp_path / "manifest.json",
        profile="curated",
        host_species="none",
        viral_identifiers={"V1"},
    )
    row = json.loads(manifest.read_text())["sequences"][0]
    assert row["retrieved_at"] is None
    assert row["retrieval_status"] == "unknown_local_input"
    assert row["source"] == "local_input"
    assert row["taxonomy"]["organism"] is None
    assert row["source_licence"]["status"] == "unreviewed"


def test_cached_provenance_retains_recorded_retrieval_and_taxonomy(tmp_path):
    acc = "NC_999999.1"
    cache = tmp_path / "cache"
    entry = cache / acc
    entry.mkdir(parents=True)
    fasta = entry / f"{acc}.fasta"
    gb = entry / f"{acc}.gb"
    ncbi_fetch._write_cached(fasta, f">{acc}\nACGTCAGTAC\n")
    ncbi_fetch._write_cached(
        gb,
        f'LOCUS       {acc} 10 bp DNA linear\nVERSION     {acc}\nFEATURES             Location/Qualifiers\n     source          1..10\n                     /organism="Synthetic virus"\n                     /db_xref="taxon:1234"\nORIGIN\n//\n',
    )
    known = ncbi_fetch.cached_accession_provenance(acc, cache)
    assert known["taxonomy"] == {"organism": "Synthetic virus", "taxid": 1234}
    date = known["retrieved_at"]
    assert date is not None
    assert ncbi_fetch.cached_accession_provenance(acc, cache)["retrieved_at"] == date
    fasta.with_suffix(".fasta.retrieval.json").unlink()
    unknown = ncbi_fetch.cached_accession_provenance(acc, cache)
    assert unknown["retrieved_at"] is None
    assert unknown["retrieval_status"] == "unknown_legacy_cache"


def _inputs(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">ENST1\n" + "A" * 60 + "\n>V1\n" + "ACGTCAGTACGT" * 5 + "\n")
    gtf = tmp_path / "input.gtf"
    gtf.write_text(
        'ENST1\ts\texon\t1\t60\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
        'V1\ts\texon\t1\t60\t.\t+\t.\tgene_id "V1_g"; transcript_id "V1_tx";\n'
    )
    return fasta, gtf


def test_production_preparation_masks_only_viral_records_and_preserves_gtf(tmp_path, monkeypatch):
    fasta, gtf = _inputs(tmp_path)
    original_gtf = gtf.read_bytes()
    seen = []

    def fake_mask(source, output):
        records = builder.validate_reference_records(source)
        seen.extend(records)
        builder._write_fasta(output, [(name, "N" * len(seq)) for name, seq in records])
        return True

    monkeypatch.setattr(builder, "mask_low_complexity", fake_mask)
    prepared, manifest = builder.prepare_reference_inputs(
        fasta, gtf, tmp_path / "out", viral_gene_ids={"ENSG1", "V1_g"}
    )
    assert [name for name, _ in seen] == ["V1"]
    assert dict(builder.validate_reference_records(prepared)) == {"ENST1": "A" * 60, "V1": "N" * 60}
    assert gtf.read_bytes() == original_gtf
    data = json.loads(manifest.read_text())
    assert data["input_fasta_sha256"] == builder.sha256_file(fasta)
    assert data["masking"]["windows"] == [64, 30]


@pytest.mark.parametrize("failure", ["duplicate_id", "duplicate_sequence", "mask_unavailable"])
def test_production_preparation_fails_before_indexing(tmp_path, monkeypatch, failure):
    fasta, gtf = _inputs(tmp_path)
    if failure == "duplicate_id":
        fasta.write_text(">V1\nACGT\n>V1\nTGCA\n")
    elif failure == "duplicate_sequence":
        fasta.write_text(">V1\nACGT\n>V2\nACGT\n")
    monkeypatch.setattr(builder, "mask_low_complexity", lambda *args: False)
    with pytest.raises((ValueError, RuntimeError)):
        builder.prepare_reference_inputs(fasta, gtf, tmp_path / "out", viral_gene_ids={"V1_g"})


def test_frozen_reference_content_digest_ignores_build_time_and_output_location(
    tmp_path, monkeypatch
):
    fasta, gtf = _inputs(tmp_path)
    monkeypatch.setattr(
        builder,
        "mask_low_complexity",
        lambda source, output: builder._write_fasta(output, [("V1", "N" * 60)]) or True,
    )
    audits = []
    contents = []
    for name in ("one", "two"):
        prepared, manifest = builder.prepare_reference_inputs(
            fasta, gtf, tmp_path / name, viral_gene_ids={"V1_g"}
        )
        t2g = tmp_path / name / "t2g.txt"
        t2g.write_text("ENST1\tENSG1\nV1_tx\tV1_g\n")
        builder.record_reference_build(manifest, fasta=prepared, gtf=gtf, t2g=t2g)
        audits.append((tmp_path / name / "reference_reproducibility.json").read_bytes())
        contents.append(prepared.read_bytes())
    assert contents[0] == contents[1]
    assert audits[0] == audits[1]


def test_cluster_member_decisions_retain_the_selected_representative(tmp_path):
    report = tmp_path / "panel.fa.clstr"
    report.write_text(">Cluster 0\n0\t60nt, >NC_1.1... *\n1\t60nt, >NC_2.1... at +/96.0%\n")
    data = builder._cluster_provenance(report)
    assert data["NC_1.1"]["representative_status"] == "representative"
    assert data["NC_2.1"] == {
        "cluster": "0",
        "representative_accession": "NC_1.1",
        "representative_status": "excluded_cluster_member",
    }


def test_reference_receipt_resources_and_index_do_not_change_content_hash(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">V1\nACGTCAGTAC\n")
    gtf = tmp_path / "input.gtf"
    gtf.write_text('V1\ts\texon\t1\t10\t.\t+\t.\tgene_id "V1_g"; transcript_id "V1_tx";\n')
    manifest = builder.write_reference_manifest(
        fasta,
        tmp_path / "reference_manifest.json",
        profile="curated",
        host_species="none",
        viral_identifiers={"V1"},
    )
    builder.record_reference_build(manifest, fasta=fasta, gtf=gtf)
    first = json.loads(manifest.read_text())["content_sha256"]
    index = tmp_path / "index.idx"
    index.write_bytes(b"toy binary index")
    builder.record_reference_build(
        manifest,
        fasta=fasta,
        gtf=gtf,
        index=index,
        resources={"elapsed_seconds": 2.5},
        command=["kb", "ref", "toy.fa"],
    )
    data = json.loads(manifest.read_text())
    assert data["content_sha256"] == first
    assert data["build_receipt"]["index"]["sha256"] == builder.sha256_file(index)
    assert data["build_receipt"]["resources"]["elapsed_seconds"] == 2.5


def test_kb_ref_resources_record_child_peak_rss(tmp_path, monkeypatch):
    monkeypatch.setattr(builder.subprocess, "run", lambda *a, **k: None)
    resources = builder._run_kb_ref(["kb", "ref", "x"], tmp_path)
    assert resources["peak_rss_kib"] > 0
    assert resources["peak_rss_scope"] == "max_over_all_child_processes"


def test_build_receipt_records_tool_versions(tmp_path, monkeypatch):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">V1\nACGTCAGTAC\n")
    gtf = tmp_path / "input.gtf"
    gtf.write_text('V1\ts\texon\t1\t10\t.\t+\t.\tgene_id "V1_g"; transcript_id "V1_tx";\n')
    manifest = builder.write_reference_manifest(
        fasta,
        tmp_path / "reference_manifest.json",
        profile="curated",
        host_species="none",
        viral_identifiers={"V1"},
    )
    monkeypatch.setattr(
        builder,
        "tool_provenance",
        lambda names: (
            {n: {"path": f"/bin/{n}", "version": None, "sha256": "0"} for n in names}
            | {"cd-hit-est": None}
        ),
    )

    def fake_run(cmd, **kwargs):
        text = {
            "kb": "kb_python 9.9.9",
            "dustmasker": "dustmasker: 1.2.3\n Package: blast",
        }.get(Path(cmd[0]).name, "")
        return subprocess.CompletedProcess(cmd, 0, stdout=text, stderr="")

    monkeypatch.setattr(builder.subprocess, "run", fake_run)
    builder.record_reference_build(manifest, fasta=fasta, gtf=gtf, command=["kb", "ref"])
    tools = json.loads(manifest.read_text())["build_receipt"]["tools"]
    assert tools["kb"]["version"] == "9.9.9"
    assert tools["dustmasker"]["version"] == "1.2.3"
    assert tools["cd-hit-est"] is None
    assert set(tools) == {"kb", "kallisto", "bustools", "dustmasker", "cd-hit-est"}


def test_production_preparation_rejects_a_reference_with_no_viral_records(tmp_path, monkeypatch):
    fasta, gtf = _inputs(tmp_path)
    monkeypatch.setattr(builder, "mask_low_complexity", lambda *args: True)
    with pytest.raises(ValueError, match="No viral sequence"):
        builder.prepare_reference_inputs(fasta, gtf, tmp_path / "out", viral_gene_ids={"other"})
    assert not (tmp_path / "out" / "reference.prepared.fa").exists()

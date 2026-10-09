"""Independent exact-control receipts: flags and group names are not evidence."""

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from viralscan.run_safety import sha256_file
from viralscan.runconfig import RunConfig
from viralscan.scripts import detection as D
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable


@pytest.fixture
def exact_fixture(tmp_path, monkeypatch):
    accession, gene = "NC_123456.1", "NC_123456.1_control"
    index = tmp_path / "index.idx"
    index.write_bytes(b"toy immutable index")
    planted = tmp_path / "planting_truth.tsv"
    planted.write_text("gene\texpected_molecules\n" + gene + "\t100\n")
    digest = hashlib.sha256(b"ACGTCAGTAC").hexdigest()
    manifest = tmp_path / "reference_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": "3.0.0",
                "profile": "curated",
                "created_at": "2026-01-01T00:00:00Z",
                "fasta_sha256": digest,
                "build_receipt": {
                    "built_at": "2026-01-01T00:00:00Z",
                    "software_identity": {},
                    "index": {"sha256": sha256_file(index)},
                },
                "sequences": [
                    {
                        "accession_version": accession,
                        "sha256": digest,
                        "length": 10,
                        "source": "local_input",
                        "retrieved_at": None,
                        "inclusion_rationale": "control",
                    }
                ],
            }
        )
    )
    receipt = tmp_path / "control.receipt.json"
    data = {
        "schema_version": "3.0.0",
        "gene": gene,
        "accession_version": accession,
        "sequence_sha256": digest,
        "index_sha256": sha256_file(index),
        "reference_fasta_sha256": digest,
        "expected_molecules": 100,
        "count_layer": "X",
        "count_strategy": "kb_count",
        "reference_manifest": {"path": manifest.name, "sha256": sha256_file(manifest)},
        "evidence": {"path": planted.name, "sha256": sha256_file(planted)},
    }
    receipt.write_text(json.dumps(data))
    identity = VirusIdentityTable(
        (
            GeneIdentity(
                gene,
                accession,
                "catalogued",
                True,
                virus_key="taxid:123",
                virus_name="Target virus",
            ),
        )
    )
    monkeypatch.setattr(D, "identity", identity)
    cfg = RunConfig(
        index=str(index),
        positive_control_gene=gene,
        positive_control_expected_molecules=100,
        positive_control_scope="exact_sequence",
        positive_control_virus_key="Target virus",
        positive_control_receipt=str(receipt),
        multimapping=False,
    )
    a = ad.AnnData(X=np.array([[50.0]]))
    a.var_names = [gene]
    return a, cfg, receipt, data


def test_exact_flags_alone_retain_recovery_but_certify_nothing(exact_fixture):
    a, cfg, receipt, _ = exact_fixture
    receipt.unlink()
    capture, detail = D.measure_positive_control(a, cfg)
    assert capture == 0.5 and detail["status"] == "measured"
    assert detail["verification_status"] == "unavailable"
    assert D.positive_control_from(cfg, detail) is None


def test_valid_receipt_certifies_only_its_identity_resolved_exact_sequence(exact_fixture):
    a, cfg, _, _ = exact_fixture
    capture, detail = D.measure_positive_control(a, cfg)
    assert capture == 0.5
    assert detail["verification_status"] == "verified"
    control = D.positive_control_from(cfg, detail)
    assert control is not None
    assert D.certified_viruses(control, ["Target virus", "Unrelated virus"]) == ["Target virus"]


@pytest.mark.parametrize(
    "field",
    [
        "gene",
        "accession_version",
        "sequence_sha256",
        "index_sha256",
        "reference_fasta_sha256",
        "expected_molecules",
        "count_layer",
        "count_strategy",
        "evidence",
    ],
)
def test_mismatched_receipts_retain_measurement_without_certification(exact_fixture, field):
    a, cfg, path, data = exact_fixture
    if field == "expected_molecules":
        data[field] = 200
    elif field == "evidence":
        Path(path.parent / data[field]["path"]).write_text("tampered truth")
    elif field.endswith("sha256"):
        data[field] = "0" * 64
    else:
        data[field] = "mismatch"
    path.write_text(json.dumps(data))
    capture, detail = D.measure_positive_control(a, cfg)
    assert capture == 0.5
    assert detail["verification_status"] == "failed"
    assert detail["verification_detail"]
    assert D.positive_control_from(cfg, detail) is None


def test_one_exact_sequence_does_not_certify_a_multigenome_group(exact_fixture, monkeypatch):
    a, cfg, _, _ = exact_fixture
    first = D.identity.genes[0]
    second = GeneIdentity(
        "OTHER_gene",
        "NC_654321.1",
        "catalogued",
        True,
        virus_key=first.virus_key,
        virus_name=first.virus_name,
    )
    monkeypatch.setattr(D, "identity", VirusIdentityTable((first, second)))
    _, detail = D.measure_positive_control(a, cfg)
    assert detail["verification_status"] == "verified"
    assert detail["certification_status"] == "group_contains_other_accessions"
    assert D.positive_control_from(cfg, detail) is None


def test_report_archives_measurement_without_rewriting_it(exact_fixture, tmp_path):
    a, cfg, _, _ = exact_fixture
    _, detail = D.measure_positive_control(a, cfg)
    path = Path(
        D.write_control_report(detail, str(tmp_path), control=D.positive_control_from(cfg, detail))
    )
    report = json.loads(path.read_text())
    archived = tmp_path / report["measurement_archive"]["path"]
    original = archived.read_bytes()
    assert sha256_file(archived) == report["measurement_archive"]["sha256"]
    D.write_control_report(detail, str(tmp_path), control=D.positive_control_from(cfg, detail))
    assert archived.read_bytes() == original


def test_target_resolves_through_the_identity_table_by_key_or_name(exact_fixture):
    """The flag may name the identity's virus_key or display name; the row is its name."""
    a, cfg, _, _ = exact_fixture
    cfg = RunConfig(**{**vars(cfg), "positive_control_virus_key": "taxid:123"})
    _, detail = D.measure_positive_control(a, cfg)
    assert detail["verification_status"] == "verified"
    control = D.positive_control_from(cfg, detail)
    assert control.target == "Target virus"  # resolved, not the declared string
    assert D.certified_viruses(control, ["taxid:123", "Target virus"]) == ["Target virus"]


def test_a_declared_target_that_the_identity_table_does_not_own_certifies_nothing(exact_fixture):
    a, cfg, _, _ = exact_fixture
    cfg = RunConfig(**{**vars(cfg), "positive_control_virus_key": "Unrelated virus"})
    _, detail = D.measure_positive_control(a, cfg)
    assert detail["verification_status"] == "failed"
    assert "identity" in detail["verification_detail"].lower()
    assert D.positive_control_from(cfg, detail) is None


def _report(detail, cfg, tmp_path):
    control = D.positive_control_from(cfg, detail)
    certified = D.certified_viruses(control, ["Target virus"])
    path = D.write_control_report(
        detail, str(tmp_path), control=control, certified_targets=certified
    )
    return json.loads(Path(path).read_text())


def test_report_schema_covers_verified_failed_and_unavailable_controls(exact_fixture, tmp_path):
    """Each outcome writes a schema-valid report (write_control_report enforces it)."""
    a, cfg, receipt, data = exact_fixture
    _, detail = D.measure_positive_control(a, cfg)
    verified = _report(detail, cfg, tmp_path / "verified")
    assert verified["verification_status"] == "verified"
    assert verified["certified_targets"] == ["Target virus"]
    assert verified["receipt_sha256"] == sha256_file(receipt)

    data["accession_version"] = "NC_000000.1"
    receipt.write_text(json.dumps(data))
    _, detail = D.measure_positive_control(a, cfg)
    failed = _report(detail, cfg, tmp_path / "failed")
    assert failed["verification_status"] == "failed" and failed["verification_detail"]
    assert failed["certified_targets"] == [] and failed["certifies_negatives"] is False

    receipt.unlink()
    _, detail = D.measure_positive_control(a, cfg)
    unavailable = _report(detail, cfg, tmp_path / "unavailable")
    assert unavailable["verification_status"] == "unavailable"
    assert unavailable["certified_targets"] == [] and unavailable["certifies_negatives"] is False


@pytest.mark.parametrize(
    "detail",
    [
        {"status": "not-configured"},
        {"status": "gene-not-in-reference", "gene": "g", "detail": "absent"},
    ],
)
def test_report_schema_covers_controls_that_never_measured(detail, tmp_path):
    path = D.write_control_report(detail, str(tmp_path))
    assert json.loads(Path(path).read_text())["certifies_negatives"] is False


@pytest.mark.parametrize(
    "forged",
    [
        {"verification_status": "unavailable", "verification_detail": "x"},
        {"verification_status": "failed", "verification_detail": "x"},
        {"verification_status": "verified"},  # no receipt, digest or accession recorded
        {"scope": "exact_sequence"},  # exact control certifying with no verification at all
    ],
)
def test_report_schema_refuses_an_unverified_report_that_certifies(forged):
    from viralscan.validation import SchemaContractError, require_schema_valid

    payload = {
        "schema_version": "3.0.0",
        "status": "measured",
        "certified_targets": ["Target virus"],
        "informative_negative_targets": ["Target virus"],
        "certifies_negatives": True,
        "capture_used_for_sensitivity": 0.5,
        **forged,
    }
    with pytest.raises(SchemaContractError):
        require_schema_valid(payload, "positive_control_report.schema.json")


def test_shipped_schema_copies_are_identical():
    root = Path(__file__).resolve().parents[1]
    for name in ("positive_control_receipt.schema.json", "positive_control_report.schema.json"):
        assert (root / "schemas/v3" / name).read_bytes() == (
            root / "src/viralscan/schemas/v3" / name
        ).read_bytes()

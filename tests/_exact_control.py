"""Shared fixture for exact_sequence positive controls (SENS-CORR-01).

An ``exact_sequence`` control certifies a negative only with a verified receipt:
checksum-bound index, sequence manifest and planting evidence, and an indexed viral
accession in the VirusIdentityTable. ``exact_control_kwargs`` builds all of it for a
toy control gene and returns the ``RunConfig`` kwargs that point at it.
"""

from __future__ import annotations

import hashlib
import json

from viralscan.run_safety import sha256_file
from viralscan.scripts import detection as D
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

ACCESSION = "NC_123456.1"


def write_exact_control(
    root, gene, expected=100.0, accession=ACCESSION, count_layer="X", strategy="kb_count"
):
    """Write index, manifest, evidence and receipt under ``root``; return (index, receipt, data)."""
    root.mkdir(parents=True, exist_ok=True)
    index = root / "index.idx"
    index.write_bytes(b"toy immutable index")
    evidence = root / "planting_truth.tsv"
    evidence.write_text(f"gene\texpected_molecules\n{gene}\t{expected:g}\n")
    digest = hashlib.sha256(b"ACGTCAGTAC").hexdigest()
    manifest = root / "reference_manifest.json"
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
    data = {
        "schema_version": "3.0.0",
        "gene": gene,
        "accession_version": accession,
        "sequence_sha256": digest,
        "index_sha256": sha256_file(index),
        "reference_fasta_sha256": digest,
        "expected_molecules": expected,
        "count_layer": count_layer,
        "count_strategy": strategy,
        "reference_manifest": {"path": manifest.name, "sha256": sha256_file(manifest)},
        "evidence": {"path": evidence.name, "sha256": sha256_file(evidence)},
    }
    receipt = root / "control.receipt.json"
    receipt.write_text(json.dumps(data))
    return index, receipt, data


def exact_control_kwargs(tmp_path, monkeypatch, gene, virus_name, expected=100.0):
    """RunConfig kwargs for a verified exact control of ``gene`` -> ``virus_name``.

    Also installs the one-gene VirusIdentityTable the verification resolves through.
    """
    index, receipt, _ = write_exact_control(tmp_path / "exact_control", gene, expected)
    table = VirusIdentityTable(
        (
            GeneIdentity(
                gene,
                ACCESSION,
                "catalogued",
                True,
                virus_key=f"key:{virus_name}",
                virus_name=virus_name,
            ),
        )
    )
    monkeypatch.setattr(D, "identity", table)
    return {
        "index": str(index),
        "positive_control_receipt": str(receipt),
        "positive_control_scope": "exact_sequence",
        "positive_control_virus_key": virus_name,
        "multimapping": False,  # the receipt names the "kb_count" strategy
    }

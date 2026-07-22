"""G1: reference-annotation provenance recorded alongside results."""

from __future__ import annotations

import json
from pathlib import Path

from viralscan.runconfig import RunConfig
from viralscan.scripts.detection import reference_provenance, write_reference_provenance


def _cfg():
    return RunConfig(
        index="ref/index.idx",
        transcripts="ref/t2g.txt",
        gtf="ref/viral.gtf",
        technology="10xv2",
        multimap_method="em-global",
        multimap_primary_call="selected-method",
    )


def test_reference_provenance_records_reference_and_accessions():
    prov = reference_provenance(_cfg(), ["NC_007605.1", "NC_001664.4"], ["Epstein-Barr virus"])
    assert prov["index"] == "ref/index.idx"
    assert prov["multimap_primary_call"] == "selected-method"
    assert prov["n_viral_accessions_in_reference"] == 2
    assert prov["viral_accessions"] == ["NC_001664.4", "NC_007605.1"]  # sorted
    assert prov["viruses_detected"] == ["Epstein-Barr virus"]
    assert prov["viralscan_version"]  # non-empty


def test_write_reference_provenance_emits_valid_json(tmp_path: Path):
    path = write_reference_provenance(_cfg(), ["NC_007605.1"], ["EBV"], str(tmp_path))
    assert Path(path).name == "reference_provenance.json"
    data = json.loads(Path(path).read_text())
    assert data["n_viruses_detected"] == 1 and data["viral_accessions"] == ["NC_007605.1"]

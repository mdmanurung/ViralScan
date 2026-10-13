"""VAL-RA guard: ``virus_key`` positive-control certification stays rejected on every config entry point.

VAL-RA-01..04 and VAL-RA-CAL are conditional on matching evidence (an approved transfer calibration,
SENS-CAL-01 on the cluster). Until it exists nothing may certify a negative for a virus other than the
one the control was measured on: not the CLI/Snakemake config, and not a saved ``config.yaml`` that
detection and the ``rerun-*`` commands read through ``RunConfig.from_yaml``.
"""

from __future__ import annotations

import pytest
import yaml

from tests.test_positive_control import TestScopeConfig
from viralscan.runconfig import VIRUS_KEY_REJECTED, RunConfig

BASE = {
    "positive_control_gene": "SPIKE",
    "positive_control_expected_molecules": "100",
    "positive_control_scope": "virus_key",
    "positive_control_virus_key": "HHV4",
}


def test_snakemake_config_rejects_virus_key():
    with pytest.raises(ValueError, match="transfer calibration"):
        TestScopeConfig()._cfg(**BASE)


def test_saved_config_yaml_cannot_smuggle_virus_key_past_the_guard(tmp_path):
    good = TestScopeConfig()._cfg(
        positive_control_scope="exact_sequence", positive_control_virus_key="HHV4"
    )
    path = tmp_path / "config.yaml"
    good.to_yaml(path)
    assert RunConfig.from_yaml(path).positive_control_scope == "exact_sequence"
    data = yaml.safe_load(path.read_text())
    data["positive_control_scope"] = "virus_key"  # a hand edit
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError, match="transfer calibration"):
        RunConfig.from_yaml(path)


def test_the_rejection_message_names_the_missing_evidence_and_the_alternatives():
    assert "VAL-RA-CAL" in VIRUS_KEY_REJECTED
    assert "exact_sequence" in VIRUS_KEY_REJECTED and "panel_mechanics" in VIRUS_KEY_REJECTED

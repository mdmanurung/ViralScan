"""PANEL-01 WP2: scripts/panel_pool.py marks Virus-Host DB human viruses by where their accessions already are."""

from __future__ import annotations

import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "panel_pool", Path(__file__).resolve().parents[1] / "scripts" / "panel_pool.py"
)
pp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pp)


def test_lineage_names() -> None:
    lineage = "Viruses; Riboviria; Pisuviricota; Picornavirales; Picornaviridae; Paavivirinae; Aalivirus; Aalivirus apekidu"
    assert pp.lineage_names(lineage) == ("Picornaviridae", "Aalivirus")
    assert pp.lineage_names("Viruses; Unclassified") == ("", "")


def test_classify_prefers_panel_then_candidate_and_defaults_to_gap() -> None:
    panel_of = {"A": "shipped", "B": "max", "C": "max"}
    cand = {"C": False, "D": True}
    assert pp.classify(["X", "A"], panel_of, cand) == ("panel", "shipped")
    assert pp.classify(["D", "C"], panel_of, cand) == ("candidate", "max")
    assert pp.classify(["D"], panel_of, cand) == ("excluded", "")
    assert pp.classify(["B"], panel_of, cand) == ("catalogue", "max")
    assert pp.classify(["X"], panel_of, cand) == ("gap", "")
    assert pp.classify([], panel_of, cand) == ("gap", "")


def test_tier_ranks_disease_over_evidence() -> None:
    assert pp.tier({"DISEASE": "Viral wart", "evidence": "Literature"}) == "disease"
    assert pp.tier({"DISEASE": "", "evidence": "Literature, RefSeq"}) == "refseq_evidence"
    assert pp.tier({"DISEASE": "", "evidence": "Literature"}) == "literature_only"

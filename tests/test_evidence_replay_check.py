"""EVID-COUNTED-01: the replay must agree with the primary kb count's barcodes."""

from __future__ import annotations

import pytest

from viralscan.kb_outputs import KbCountOutputs
from viralscan.scripts.evidence_run import check_replay_matches_matrix


def _lineage(*cbs: str, weight: float = 1.0) -> dict[int, dict[str, object]]:
    return {i: {"cb": cb, "assigned_weight": weight} for i, cb in enumerate(cbs, 1)}


def _kb(tmp_path, barcodes: list[str] | None, corrected: bool) -> KbCountOutputs:
    kb = KbCountOutputs(tmp_path)
    kb.counts_dir.mkdir(parents=True)
    if barcodes is not None:
        kb.barcodes.write_text("\n".join(barcodes) + "\n")
    if corrected:
        kb.kb_corrected_bus.write_text("")
    return kb


def test_matching_replay_passes(tmp_path) -> None:
    kb = _kb(tmp_path, ["AAAA", "CCCC"], corrected=True)
    check_replay_matches_matrix(_lineage("AAAA", "CCCC"), kb, "kb-python/whitelist.txt")


def test_primary_corrected_but_replay_has_no_whitelist_fails(tmp_path) -> None:
    kb = _kb(tmp_path, ["AAAA"], corrected=True)
    with pytest.raises(ValueError, match="primary run corrected barcodes"):
        check_replay_matches_matrix(_lineage("AAAA"), kb, None)


def test_off_matrix_barcode_in_a_counted_read_fails(tmp_path) -> None:
    kb = _kb(tmp_path, ["AAAA"], corrected=False)
    with pytest.raises(ValueError, match="1 of 2 replay barcodes"):
        check_replay_matches_matrix(_lineage("AAAA", "GGGG"), kb, None)


def test_weight_zero_reads_are_not_checked(tmp_path) -> None:
    kb = _kb(tmp_path, ["AAAA"], corrected=False)
    check_replay_matches_matrix(_lineage("GGGG", weight=0.0), kb, None)

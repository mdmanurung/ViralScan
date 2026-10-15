"""AUD3-05 (audit P2-01): near-duplicate screen between reference records, and the host-homology warning.

The screen reuses the curation-time twin cut (`viralscan.kmer_sharing.KMER_TWIN_BELOW`, canonical 31-mers), so a
diverged homologue below ~95 % shared k-mers is NOT flagged here: its threshold belongs to PLAN REF-08.
"""

import importlib.util
import json
import logging
import random
from pathlib import Path

import pytest

from viralscan.scripts import build_reference as br

REPO = Path(__file__).resolve().parents[1]


def _dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _write(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in records))
    return path


HOST = _dna(1500, 1)
VIRAL = _dna(1500, 2)


def _one_mismatch(seq: str) -> str:
    sub = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return seq[:700] + sub[seq[700]] + seq[701:]


def _revcomp(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _diverged(seq: str, rate: float = 0.10) -> str:
    rng = random.Random(3)
    return "".join(c if rng.random() > rate else rng.choice("ACGT") for c in seq)


@pytest.mark.parametrize(
    "make", [_one_mismatch, _revcomp], ids=["1-mismatch", "reverse-complement"]
)
def test_near_copy_is_warned_by_default_and_fails_under_strict(tmp_path, caplog, make):
    fa = _write(tmp_path / "r.fa", [("HOST_TX", HOST), ("VIR1.1", VIRAL), ("VIR2.1", make(HOST))])
    with caplog.at_level(logging.WARNING, logger="viralscan"):
        assert [i for i, _ in br.validate_reference_records(fa)] == ["HOST_TX", "VIR1.1", "VIR2.1"]
    assert "near-duplicate" in caplog.text and "VIR2.1" in caplog.text and "HOST_TX" in caplog.text
    with pytest.raises(br.NearDuplicateError, match="VIR2.1"):
        br.validate_reference_records(fa, near_duplicates="error")
    br.validate_reference_records(fa, near_duplicates="off")  # opt-out stays silent


def test_ten_percent_diverged_copy_is_not_flagged_threshold_is_ref08(tmp_path, caplog):
    # about 9 % of its 31-mers survive, far from the 95 % twin cut; any lower cut is PLAN REF-08's to set
    fa = _write(tmp_path / "r.fa", [("HOST_TX", HOST), ("VIR2.1", _diverged(HOST))])
    with caplog.at_level(logging.WARNING, logger="viralscan"):
        br.validate_reference_records(fa, near_duplicates="error")
    assert "near-duplicate" not in caplog.text


def test_exact_duplicate_still_fails_closed(tmp_path):
    fa = _write(tmp_path / "r.fa", [("HOST_TX", HOST), ("VIR1.1", HOST)])
    with pytest.raises(ValueError, match="Exact duplicate sequences"):
        br.validate_reference_records(fa)


def test_declared_siblings_are_not_flagged(tmp_path, caplog):
    fa = _write(tmp_path / "r.fa", [("STRAIN_A.1", HOST), ("STRAIN_B.1", _one_mismatch(HOST))])
    with caplog.at_level(logging.WARNING, logger="viralscan"):
        br.validate_reference_records(
            fa, near_duplicates="error", sibling_groups={"STRAIN_A": "hhv6", "STRAIN_B": "hhv6"}
        )
    assert "near-duplicate" not in caplog.text
    with pytest.raises(br.NearDuplicateError):  # different declared groups are not siblings
        br.validate_reference_records(
            fa, near_duplicates="error", sibling_groups={"STRAIN_A": "x", "STRAIN_B": "y"}
        )


def test_viral_identifiers_scope_the_screen_to_pairs_with_a_viral_record(tmp_path, caplog):
    host_iso = _one_mismatch(HOST)
    fa = _write(tmp_path / "r.fa", [("HOST_TX", HOST), ("HOST_ISO", host_iso), ("VIR1.1", VIRAL)])
    # host-host near-copies (isoforms/paralogs) are normal and ignored; the only viral record is unique
    br.validate_reference_records(fa, near_duplicates="error", viral_identifiers={"VIR1.1"})
    fa2 = _write(tmp_path / "r2.fa", [("HOST_TX", HOST), ("VIR2.1", _revcomp(HOST))])
    with pytest.raises(br.NearDuplicateError, match="VIR2.1"):
        br.validate_reference_records(fa2, near_duplicates="error", viral_identifiers={"VIR2.1"})


def test_manifest_warns_that_host_homology_was_not_measured(tmp_path, caplog):
    fa = _write(tmp_path / "r.fa", [("VIR1.1", VIRAL)])
    with caplog.at_level(logging.WARNING, logger="viralscan"):
        m = br.write_reference_manifest(
            fa,
            tmp_path / "m.json",
            profile="curated",
            host_species="human",
            viral_identifiers={"VIR1.1"},
        )
    assert "Host homology was NOT measured" in caplog.text
    assert (
        json.loads(m.read_text())["sequences"][0]["host_homology_status"]
        == "not_assessed_no_host_genome"
    )


def test_manifest_strict_reference_fails_on_non_sibling_pair(tmp_path):
    fa = _write(tmp_path / "r.fa", [("HOST_TX", HOST), ("VIR2.1", _one_mismatch(HOST))])
    kwargs = dict(profile="curated", host_species="human", viral_identifiers={"VIR2.1"})
    br.write_reference_manifest(fa, tmp_path / "m.json", **kwargs)  # default: warn only
    with pytest.raises(br.NearDuplicateError):
        br.write_reference_manifest(fa, tmp_path / "m2.json", strict_reference=True, **kwargs)


def test_build_ref_cli_flag_defaults_off_and_toggles():
    from unittest.mock import patch

    from viralscan.menu import create_help

    def parse(*extra):
        with patch("sys.argv", ["viralscan", "build-ref", *extra]):
            return create_help()

    assert parse().strict_reference is False
    assert parse("--strict-reference").strict_reference is True
    assert parse("--strict-reference", "--no-strict-reference").strict_reference is False


def test_bundled_panel_builder_has_the_same_flag():
    spec = importlib.util.spec_from_file_location(
        "bbpr", REPO / "scripts" / "build_bundled_panel_ref.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    parser = mod._build_arg_parser()
    default = parser.parse_args([]) if not any(a.required for a in parser._actions) else None
    if default is not None:
        assert default.strict_reference is False
    dests = {a.dest for a in parser._actions}
    assert "strict_reference" in dests

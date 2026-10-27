"""DISC-15: ``discover-states`` picks evidence sources by availability and records provenance."""

import argparse
import csv
import json
from pathlib import Path

import pytest

from viralscan import discovery
from viralscan.discovery import DiscoveryError

FP = "a" * 64


def _candidate(taxid, name, groups_called, host_lineage="0", viral_only="0") -> dict:
    row = {c: "0" for c in discovery.CANDIDATE_COLUMNS}
    row.update(
        taxid=str(taxid),
        name=name,
        rank="S",
        groups_called=str(groups_called),
        state="observation",
        host_relation="unrelated",
        host_lineage=host_lineage,
        viral_only_db=viral_only,
    )
    return row


def _discovery(tmp_path: Path, candidates, run_fp=FP) -> Path:
    out = tmp_path / "disc"
    out.mkdir()
    with (out / discovery.CANDIDATES_TSV).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=discovery.CANDIDATE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidates)
    (out / discovery.PROVENANCE_JSON).write_text(
        json.dumps(
            {
                "primary_run_fingerprint": run_fp,
                "discovery_fingerprint": "d" * 64,
                "resolved_input": {"r1": "R1.fq.gz", "r2": "R2.fq.gz"},
            }
        )
    )
    (out / discovery.COMPLETE_JSON).write_text("{}")
    return out


def _run(tmp_path: Path, fp=FP) -> Path:
    root = tmp_path / "run"
    (root / "S1" / "results").mkdir(parents=True)
    (root / "S1" / "config.yaml").write_text("output: x\n")
    (root / "run_manifest.json").write_text(json.dumps({"run_fingerprint": fp}))
    (root / "run_complete.json").write_text("{}")
    return root


def _evidence(out: Path, name: str, taxid: int, source: str) -> Path:
    path = out / name
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(discovery.EVIDENCE_COLUMNS)
        writer.writerow([taxid, "corroboration", source, 10, 10, 2, 5])
    return path


def _states(out: Path) -> dict[str, dict[str, str]]:
    with (out / discovery.CANDIDATE_STATES_TSV).open() as handle:
        return {r["taxid"]: r for r in csv.DictReader(handle, delimiter="\t")}


def _args(out: Path, **kw) -> argparse.Namespace:
    base = dict(
        discovery_dir=out, run_dir=None, sample=None, minimap2=False,
        host_fasta=None, candidate_fasta=[], panel_fasta=None,
    )  # fmt: skip
    base.update(kw)
    return argparse.Namespace(**base)


def test_states_cli_picks_star_over_kb(tmp_path, monkeypatch):
    out = _discovery(tmp_path, [_candidate(562, "E. coli", 5), _candidate(10239, "Viruses", 1)])
    run = _run(tmp_path)
    (run / "S1" / discovery.STAR_VIRAL_BAM).parent.mkdir()
    (run / "S1" / discovery.STAR_VIRAL_BAM).write_bytes(b"")
    (run / "S1" / discovery.MOLECULE_ASSIGNMENTS).write_bytes(b"")
    calls = []
    monkeypatch.setattr(
        discovery,
        "build_star_evidence",
        lambda d, s, output=None: (
            calls.append("star"),
            _evidence(d, "evidence_star.tsv", 562, "star"),
        )[1],
    )
    monkeypatch.setattr(
        discovery, "build_kb_evidence", lambda d, s, output=None: calls.append("kb")
    )
    record = discovery.run_discover_states(_args(out, run_dir=run))
    assert calls == ["star"] and record["route_i_source"] == "star"
    states = _states(out)
    assert states["562"]["state"] == "corroborated"
    assert states["10239"]["state"] == "observation"
    assert states["10239"]["corroboration"] == "not_assessed"


def test_states_cli_no_evidence_is_not_assessed(tmp_path):
    out = _discovery(tmp_path, [_candidate(562, "E. coli", 5), _candidate(9606, "human", 9, "1")])
    record = discovery.run_discover_states(_args(out))
    assert record["route_i_source"] == "not_assessed" and record["route_ii"] is False
    states = _states(out)
    assert states["562"]["state"] == "observation"  # Kraken2 only, however many groups
    assert states["562"]["reason"] == "no orthogonal evidence supplied"
    assert states["9606"]["state"] == "observation"  # host lineage caps it
    assert not (out / "evidence.tsv").exists()


def test_states_cli_provenance(tmp_path):
    out = _discovery(tmp_path, [_candidate(562, "E. coli", 5)])
    discovery.run_discover_states(_args(out))
    record = json.loads((out / discovery.CANDIDATE_STATES_PROVENANCE_JSON).read_text())
    assert "k" not in record and record["calibrated"] is False
    assert record["discovery_fingerprint"] == "d" * 64
    assert record["constants"]["MIN_TAXON_BEST"] == discovery.MIN_TAXON_BEST


def test_states_cli_refuses_fingerprint_mismatch_and_unpublished(tmp_path):
    out = _discovery(tmp_path, [_candidate(562, "E. coli", 5)])
    with pytest.raises(DiscoveryError, match="fingerprint"):
        discovery.run_discover_states(_args(out, run_dir=_run(tmp_path, fp="b" * 64)))
    (out / discovery.COMPLETE_JSON).unlink()
    with pytest.raises(DiscoveryError, match="not a published discovery"):
        discovery.run_discover_states(_args(out))


@pytest.mark.parametrize("groups", [0, 3, 1648, 10**6])
def test_kraken2_alone_is_always_observation(groups):
    # DISC-12 rule 2 found no k (2026-10-27): no group count lifts a Kraken2-only taxon.
    capped = _candidate(9606, "human", groups, host_lineage="1")
    assert discovery.assign_state(capped, {})["state"] == "observation"
    plain = _candidate(562, "E. coli", groups)
    assert discovery.assign_state(plain, {})["state"] == "observation"
    positive = {
        (562, "corroboration"): dict(
            n_adjudicated=10, n_taxon_best=10, n_shared_cells=1, n_shared_umis=1
        )
    }
    assert discovery.assign_state(plain, positive)["state"] == "corroborated"


def test_negative_evidence_stays_in_reason():
    cand = _candidate(562, "x", 5)
    row = dict.fromkeys(discovery.EVIDENCE_COLUMNS[3:], 0)
    row.update(taxid=562, route="corroboration", n_adjudicated=10, n_taxon_best=0)
    out = discovery.assign_state(cand, {(562, "corroboration"): row})
    assert out["state"] == "observation"
    assert out["reason"] == "corroboration: 0/10 units best for the taxon"

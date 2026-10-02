from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from viralscan import strand
from viralscan.run_safety import (
    RunSafetyError,
    build_run_manifest,
    prepare_output_directory,
    record_strand_inference,
    recorded_strand_inference,
)


def _r(f, r, u):
    return {"forward": f, "reverse": r, "unstranded": u}


@pytest.mark.parametrize(
    "rates",
    [
        _r(6.4, 47.8, 53.6),  # F-020 covid x213 full run
        _r(8.9, 46.2, 54.4),  # F-020 covid x216 full run
        _r(6.8, 53.3, 60.1),  # F-020 HHV-6B 4M
    ],
)
def test_measured_f020_rates_pick_reverse(rates) -> None:
    assert strand.infer_strand(rates) == "reverse"


# The cases below are SYNTHETIC: no 3' forward-library ratio has been measured.
def test_synthetic_forward_library_picks_forward() -> None:
    assert strand.infer_strand(_r(50, 5, 52)) == "forward"


def test_synthetic_unstranded_library_picks_unstranded() -> None:
    assert strand.infer_strand(_r(30, 28, 55)) == "unstranded"


def test_synthetic_both_clear_tau_larger_ratio_wins_tie_goes_unstranded() -> None:
    assert strand.infer_strand(_r(48, 50, 52)) == "reverse"
    assert strand.infer_strand(_r(50, 50, 52)) == "unstranded"


def test_synthetic_zero_unstranded_is_unstranded() -> None:
    assert strand.infer_strand(_r(0, 0, 0)) == "unstranded"


def test_run_pilot_uses_explicit_strands_and_head(tmp_path: Path) -> None:
    fq = tmp_path / "x.fastq"
    fq.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(5)))
    seen = []

    def fake(s, r1, r2, out, *a):
        seen.append((s, r1.read_text().count("@r")))
        return {"forward": 6.4, "reverse": 47.8, "unstranded": 53.6}[s]

    rates = strand.run_pilot(str(fq), str(fq), "i", "g", "10xv3", None, 1, n_reads=2, runner=fake)
    assert seen == [("forward", 2), ("reverse", 2), ("unstranded", 2)]
    assert strand.inference_block(rates)["choice"] == "reverse"


def _args(tmp_path: Path, **changes):
    values = {
        "sample1": str(tmp_path / "R1.fastq"),
        "sample2": str(tmp_path / "R2.fastq"),
        "output": str(tmp_path / "out"),
        "multimap_method": "host-conservative",
        "cell_calling": "emptydrops",
        "called_cells_file": None,
        "resume": False,
        "overwrite": False,
        "yes": False,
        "verbose": False,
        "quiet": False,
        "strand": "auto",
    }
    values.update(changes)
    Path(values["sample1"]).write_text("r1\n")
    Path(values["sample2"]).write_text("r2\n")
    return argparse.Namespace(**values)


def test_resume_reuses_recorded_choice_and_auto_vs_explicit_refused(tmp_path: Path) -> None:
    out = Path(_args(tmp_path).output)
    manifest = build_run_manifest(_args(tmp_path))
    assert manifest["options"]["strand"] == "auto"
    prepare_output_directory(out, manifest, resume=False, overwrite=False, yes=False)
    assert recorded_strand_inference(out, "s") is None
    block = strand.inference_block(_r(6.4, 47.8, 53.6))
    record_strand_inference(out, "s", block)
    saved = json.loads((out / "run_manifest.json").read_text())
    assert saved["run_fingerprint"] == manifest["run_fingerprint"]
    assert saved["strand_inference"]["s"]["choice"] == "reverse"
    # auto vs auto resumes, and the recorded choice is available.
    again = build_run_manifest(_args(tmp_path, resume=True))
    assert prepare_output_directory(out, again, resume=True, overwrite=False, yes=False) == "resume"
    assert recorded_strand_inference(out, "s") == block
    # auto vs explicit is refused with a strand reason.
    explicit = build_run_manifest(_args(tmp_path, resume=True, strand="reverse"))
    with pytest.raises(RunSafetyError, match="--strand differs"):
        prepare_output_directory(out, explicit, resume=True, overwrite=False, yes=False)

"""HOST-05a: ``viralscan host-loss`` writes the D15/D16 truth loss as a lower bound under posthoc/."""

import argparse
import csv
import gzip
import json
from pathlib import Path

import pytest

from viralscan import menu


def _run(tmp_path: Path, rows) -> Path:
    run = tmp_path / "run"
    (run / "host_filtered").mkdir(parents=True)
    with gzip.open(run / "host_filtered" / "fragment_lineage.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["read_id", "filter_decision", "reason"])
        w.writerows(rows)
    return run


def _args(run: Path, truth: Path, **kw) -> argparse.Namespace:
    base = dict(
        run_dir=run,
        truth=str(truth),
        viral_label="viral",
        mixed_manifest=None,
        output=None,
        verbose=False,
        quiet=False,
    )
    base.update(kw)
    return argparse.Namespace(**base)


def test_host_loss_verb_writes_posthoc_json(tmp_path: Path) -> None:
    run = _run(
        tmp_path,
        [("v1", "removed", "host_mapped"), ("v2", "retained", "host_unmapped"),
         ("v3", "removed", "host_mapped"), ("h1", "removed", "host_mapped")],
    )  # fmt: skip
    truth = tmp_path / "truth.tsv"
    truth.write_text(
        "read_id\tlabel\tmolecule_id\nv1\tviral\tm1\nv2\tviral\tm1\nv3\tviral\tm2\nh1\thost\tm9\n"
    )
    menu._run_host_loss(_args(run, truth))
    out = json.loads((run / "posthoc" / "host_loss" / "host_loss.json").read_text())
    assert out["d15_truth_fragments"] == 3 and out["d15_removed_fragments"] == 2
    assert out["d16_truth_molecules"] == 2 and out["d16_lost_molecules"] == 1
    assert out["lower_bound"] is True and out["not_estimable"] is False
    assert out["mixed_manifest"] is None


def test_host_loss_flags_not_estimable_and_dies_without_lineage(tmp_path: Path) -> None:
    run = _run(tmp_path, [("v1", "retained", "host_unmapped")])
    truth = tmp_path / "truth.tsv"
    truth.write_text("read_id\tlabel\nv1\tviral\nv_unseen\tviral\n")
    menu._run_host_loss(_args(run, truth, output=str(tmp_path / "elsewhere")))
    out = json.loads((tmp_path / "elsewhere" / "host_loss.json").read_text())
    assert out["not_estimable"] is True and out["d15_loss_fraction"] is None
    assert out["not_in_lineage"] == 1
    with pytest.raises(SystemExit):
        menu._run_host_loss(_args(tmp_path / "no_such_run", truth))

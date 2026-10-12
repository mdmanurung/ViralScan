"""Generator audit tests (PLAN VAL-01 abundance grid, VAL-02 manifest completeness). Code only.

Everything is written under ``tmp_path``. No acceptance threshold is set here: the grid test checks
the generator's own definition (abundance = distinct planted molecules per virus per sample), and the
manifest test checks that what ``run_manifest.json`` records is what is on disk.
"""

from __future__ import annotations

import hashlib
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import generate_truth_panel as gt  # noqa: E402

from tests.test_generate_truth_panel import read_tsv  # noqa: E402

LEVELS = ["1", "3", "10", "30", "100", "300", "1000"]  # protocol factor viral_abundance


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_abundance_grid_plants_every_target_at_every_level_in_every_chemistry(tmp_path):
    """VAL-01 / protocol limit_of_detection: one probit fit per virus and chemistry is identifiable.

    With contamination and collisions off, a sample at level N holds exactly N distinct (barcode,
    UMI) planted molecules for each of the six targets, so the abundance factor means what D3 says.
    """
    manifest = gt.generate(
        gt.parse_args(
            [
                "--generation-seed",
                "1",
                "--error-rate",
                "0",
                "--datasets",
                "synthetic_factorial",
                "--abundance",
                *LEVELS,
                "--cb-umi-collisions",
                "0",
                "--ambient-index-hopping",
                "0",
                "--out",
                str(tmp_path),
            ]
        )  # fmt: skip
    )
    rows = read_tsv(tmp_path / "truth" / "truth_manifest.tsv")
    molecules: dict[tuple[str, str], set[tuple[str, str]]] = {}
    for r in rows:
        if r["is_planted_molecule"] == "true":
            key = r["biological_sample_id"], r["planted_virus_id"]
            molecules.setdefault(key, set()).add((r["cell_barcode"], r["umi"]))
    levels = {r["biological_sample_id"]: r["viral_abundance_level"] for r in rows}
    grid = Counter()
    for sid, detail in manifest["samples"].items():
        level = levels[sid]
        for target in gt.TARGETS:
            assert len(molecules[(sid, target)]) == int(level), (sid, target)
            grid[(detail["chemistry"], target, level)] += 1
    assert set(grid) == {
        (c, t, level)
        for c in ("10xv2", "10xv3", "drop-seq")
        for t in gt.TARGETS
        for level in LEVELS
    }

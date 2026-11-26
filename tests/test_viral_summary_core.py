"""COLNAMES-01: ``viral_summary_core.tsv`` is a 12-column view of ``viral_summary.tsv``."""

import json
from importlib import resources
from pathlib import Path

import pandas as pd

from viralscan.outputs import V3Outputs
from viralscan.run import Run
from viralscan.scripts.detection import CORE_SUMMARY_COLUMNS, write_tsv_outputs

OLD_NAMES = (
    "background_umis",
    "infected_called",
    "infected_cells",
    "total_cells",
    "pct_infected",
    "molecules_unique_to_virus",
)


def _write(tmp_path: Path) -> V3Outputs:
    stats = {
        "V": {
            "viral_molecules_total_est": 5,
            "infected_barcodes": 1,
            "n_barcodes": 10,
            "pct_infected_barcodes": 10.0,
            "viral_molecules_per_10k_est": 1.0,
        }
    }
    write_tsv_outputs(stats, pd.DataFrame(), str(tmp_path))
    return V3Outputs(tmp_path)


def test_core_view_has_the_same_values_as_the_full_table(tmp_path: Path) -> None:
    outputs = _write(tmp_path)
    full = pd.read_csv(outputs.viral_summary, sep="\t", dtype=str, keep_default_na=False)
    core = pd.read_csv(outputs.viral_summary_core, sep="\t", dtype=str, keep_default_na=False)
    assert list(core.columns) == list(CORE_SUMMARY_COLUMNS) and len(CORE_SUMMARY_COLUMNS) == 12
    assert core.equals(full[list(CORE_SUMMARY_COLUMNS)])


def test_registry_pins_the_core_view_and_the_new_names() -> None:
    tables = json.loads(
        resources.files("viralscan.schemas.v3").joinpath("table_columns.json").read_text("utf-8")
    )["tables"]
    assert tables["viral_summary_core.tsv"] == list(CORE_SUMMARY_COLUMNS)
    for name in ("viral_summary.tsv", "combined_summary.tsv"):
        assert not set(OLD_NAMES) & set(tables[name]), name
    assert {"background_molecules", "infected_called_cells", "infected_barcodes"} <= set(
        tables["viral_summary.tsv"]
    )


def test_run_summary_core_reads_the_core_table(tmp_path: Path) -> None:
    sample = tmp_path / "s1"
    outputs = _write(sample)
    run = Run(tmp_path, sample)
    assert list(run.summary(core=True).columns) == list(CORE_SUMMARY_COLUMNS)
    assert "background_molecules" in run.summary().columns
    assert outputs.viral_summary_core.name == "viral_summary_core.tsv"

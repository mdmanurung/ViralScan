"""API-FACADE-01: ``viralscan.open_run`` reads a finished run and writes nothing."""

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import viralscan
from viralscan.run import Run, open_run


def _make_run(root: Path, sample: str = "s1", matrix: str = "adata_multimap.h5ad") -> Path:
    sample_dir = root / sample
    (sample_dir / "results").mkdir(parents=True)
    (sample_dir / "config.yaml").write_text("output: x\n")
    (root / "run_manifest.json").write_text(json.dumps({"schema_version": "3.0.0"}))
    (root / "run_complete.json").write_text("{}")
    pd.DataFrame({"virus_name": ["EBV", "HHV-6B"], "detected": [True, False]}).to_csv(
        sample_dir / "results" / "viral_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "barcode": ["AAAC", "AAAG", "AAAT"],
            "virus_name": ["EBV"] * 3,
            "is_called_cell": [True, False, True],
        }
    ).to_csv(sample_dir / "results" / "per_cell_viral.tsv", sep="\t", index=False)
    counts = sample_dir / "kb-python" / "counts_unfiltered"
    counts.mkdir(parents=True)
    adata = ad.AnnData(sp.csr_matrix(np.eye(3)))
    adata.obs_names = ["AAAC", "AAAG", "AAAT"]
    adata.var_names = ["EBV_g", "h1", "h2"]
    adata.write_h5ad(counts / matrix)
    return sample_dir


def _digest(root: Path) -> dict[str, str]:
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def test_exported_from_the_package() -> None:
    assert viralscan.open_run is open_run
    assert viralscan.Run is Run


def test_run_root_and_sample_dir_give_the_same_run(tmp_path: Path) -> None:
    sample_dir = _make_run(tmp_path)
    from_root = open_run(tmp_path)
    from_sample = open_run(sample_dir)
    assert from_root.sample == from_sample.sample == "s1"
    assert from_root.root == from_sample.root == tmp_path
    assert from_root.manifest()["schema_version"] == "3.0.0"


def test_summary_cells_and_adata(tmp_path: Path) -> None:
    _make_run(tmp_path)
    run = open_run(tmp_path)
    summary = run.summary()
    assert list(summary["virus_name"]) == ["EBV", "HHV-6B"]
    assert len(run.cells()) == 3
    called = run.cells(called_only=True)
    assert list(called["barcode"]) == ["AAAC", "AAAT"]
    adata = run.adata()
    assert adata.shape == (3, 3)
    assert list(adata.var_names) == ["EBV_g", "h1", "h2"]


def test_adata_falls_back_to_the_uncorrected_matrix(tmp_path: Path) -> None:
    _make_run(tmp_path, matrix="adata.h5ad")
    assert open_run(tmp_path).adata().n_obs == 3


def test_missing_table_and_matrix_raise(tmp_path: Path) -> None:
    sample_dir = _make_run(tmp_path)
    (sample_dir / "results" / "per_cell_viral.tsv").unlink()
    (sample_dir / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad").unlink()
    run = open_run(tmp_path)
    with pytest.raises(FileNotFoundError, match="per_cell_viral.tsv"):
        run.cells()
    with pytest.raises(FileNotFoundError, match="count matrix"):
        run.adata()


def test_several_samples_need_a_name(tmp_path: Path) -> None:
    _make_run(tmp_path, "a")
    _make_run(tmp_path, "b")
    with pytest.raises(ValueError, match="--sample"):
        open_run(tmp_path)
    assert open_run(tmp_path, sample="b").sample == "b"
    with pytest.raises(ValueError, match="not in run"):
        open_run(tmp_path, sample="c")


def test_unfinished_or_missing_run_is_a_value_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="does not exist"):
        open_run(tmp_path / "nope")
    _make_run(tmp_path)
    (tmp_path / "run_complete.json").unlink()
    with pytest.raises(ValueError, match="not a completed"):
        open_run(tmp_path)


def test_reading_never_writes(tmp_path: Path) -> None:
    _make_run(tmp_path)
    before = _digest(tmp_path)
    run = open_run(tmp_path)
    run.manifest(), run.summary(), run.cells(), run.cells(called_only=True), run.adata()
    assert _digest(tmp_path) == before

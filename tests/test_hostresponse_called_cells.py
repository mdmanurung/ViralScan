"""HOSTRESP-CALLED-01: hostresponse uses only the run's called cells."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.errors import ViralScanError
from viralscan.scripts import hostresponse as hr
from viralscan.scripts.hostresponse import run_hostresponse

N_CALLED, N_EMPTY = 80, 60


def _inputs(tmp_path: Path, n_host: int = N_CALLED + N_EMPTY):
    """Host h5ad with cell_0..cell_{n-1}; the virus matrix covers the first N_CALLED only."""
    rng = np.random.default_rng(0)
    host = ad.AnnData(
        X=sp.csr_matrix(rng.negative_binomial(5, 0.3, size=(n_host, 40)).astype(np.float32)),
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_host)]),
        var=pd.DataFrame(index=[f"gene_{j}" for j in range(40)]),
    )
    counts = np.zeros((N_CALLED, 2), dtype=np.float32)
    counts[rng.choice(N_CALLED, 20, replace=False), 0] = 3
    virus = ad.AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(N_CALLED)]),
        var=pd.DataFrame(index=["VIRUS_A", "VIRUS_B"]),
    )
    host.write_h5ad(tmp_path / "host.h5ad")
    virus.write_h5ad(tmp_path / "virus.h5ad")
    (tmp_path / "analysis.txt").write_text("VIRUS_A\nVIRUS_B\n")
    return tmp_path / "virus.h5ad", tmp_path / "host.h5ad", tmp_path / "analysis.txt"


def _called(tmp_path: Path, barcodes) -> Path:
    path = tmp_path / "results" / "called_cells.tsv"
    path.parent.mkdir(exist_ok=True)
    path.write_text("barcode\n" + "".join(f"{b}\n" for b in barcodes))
    return path


def _run(tmp_path: Path, called: Path | None):
    virus, host, acc = _inputs(tmp_path)
    return run_hostresponse(
        str(virus),
        str(host),
        str(acc),
        str(tmp_path / "out"),
        use_hvg=False,
        seeds=[0, 1],
        n_stab_iter=3,
        detection_threshold=1,
        cv_mode="cell",
        called_cells_file=None if called is None else str(called),
    )


def test_only_called_cells_reach_the_models(tmp_path, caplog):
    caplog.set_level("INFO", logger="viralscan")
    called = _called(tmp_path, [f"cell_{i}" for i in range(N_CALLED)])
    seen_rows: list[int] = []
    real = hr._run_l2_regression

    def spy(X, y, *args, **kwargs):
        seen_rows.append(X.shape[0])
        return real(X, y, *args, **kwargs)

    with patch.object(hr, "_run_l2_regression", side_effect=spy):
        prov = _run(tmp_path, called)

    assert seen_rows and set(seen_rows) == {N_CALLED}
    block = prov["called_cells"]
    assert block["n_host_barcodes"] == N_CALLED + N_EMPTY
    assert block["n_called_cells"] == N_CALLED
    assert block["n_host_called"] == N_CALLED
    assert block["n_excluded_not_called"] == N_EMPTY
    assert block["source"] == "called_cells_tsv"
    assert prov["n_cells_aligned"] == N_CALLED
    assert prov["barcode_alignment"]["n_host"] == N_CALLED
    manifest = json.loads((tmp_path / "out" / "hostresponse_manifest.json").read_text())
    assert manifest["called_cells"]["n_excluded_not_called"] == N_EMPTY
    assert any("60 excluded as not called" in r.getMessage() for r in caplog.records)


def test_missing_called_cells_file_is_an_error(tmp_path):
    missing = tmp_path / "results" / "called_cells.tsv"
    with pytest.raises(ViralScanError, match="called_cells.tsv is missing"):
        _run(tmp_path, missing)


def test_no_host_barcode_called_is_an_error(tmp_path):
    called = _called(tmp_path, ["other_1", "other_2"])
    with pytest.raises(ViralScanError, match="none of the"):
        _run(tmp_path, called)


def test_cell_calling_none_keeps_every_barcode(tmp_path):
    """With --cell-calling none detection writes every barcode as called."""
    from viralscan.scripts.cellcalling import write_called_cells

    all_barcodes = [f"cell_{i}" for i in range(N_CALLED + N_EMPTY)]
    write_called_cells(all_barcodes, np.ones(len(all_barcodes), bool), str(tmp_path))
    prov = _run(tmp_path, tmp_path / "results" / "called_cells.tsv")
    assert prov["called_cells"]["n_excluded_not_called"] == 0
    assert prov["called_cells"]["n_host_called"] == N_CALLED + N_EMPTY


def test_rerun_passes_the_run_called_cells_file(tmp_path):
    from tests.test_hostresponse_subcommand import _make_args, _make_sample_dir

    sample_dir = _make_sample_dir(tmp_path)
    fake = sample_dir / "virus.h5ad"
    fake.touch()
    kb = MagicMock()
    kb.current_adata.return_value = fake
    with (
        patch("viralscan.kb_outputs.KbCountOutputs.from_config_output", return_value=kb),
        patch("viralscan.scripts.hostresponse.run_hostresponse") as mock_run,
    ):
        from viralscan.menu import _run_hostresponse_subcommand

        _run_hostresponse_subcommand(_make_args(sample_dir))
    assert Path(mock_run.call_args.kwargs["called_cells_file"]) == (
        sample_dir.resolve() / "results" / "called_cells.tsv"
    )


def test_snakemake_entry_passes_the_run_called_cells_file(tmp_path):
    cfg = MagicMock()
    cfg.output = str(tmp_path) + "/"
    ctx = MagicMock(config=cfg)
    with (
        patch.object(hr, "run_hostresponse") as mock_run,
        patch.object(hr, "resolve_seeds", return_value=[0]),
    ):
        hr.run(ctx, str(tmp_path / "done"))
    assert Path(mock_run.call_args.kwargs["called_cells_file"]) == (
        tmp_path / "results" / "called_cells.tsv"
    )

from __future__ import annotations

import csv
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scripts.summarize_legacy_v2_v3_cells import summarize_cells


def test_summarize_cells_reports_barcode_and_nonzero_overlap(tmp_path: Path) -> None:
    legacy_path = tmp_path / "legacy.h5ad"
    legacy = ad.AnnData(
        X=sparse.csr_matrix((3, 2)),
        obs=pd.DataFrame(index=["a", "b", "c"]),
        var=pd.DataFrame({"is_viral": [True, False]}, index=["virus", "host"]),
    )
    legacy.layers["counts_combined"] = sparse.csr_matrix(
        np.array([[2.0, 0.0], [0.0, 1.0], [1.0, 0.0]])
    )
    legacy.write_h5ad(legacy_path)

    manifest = tmp_path / "manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "run_id",
                "logical_id",
                "sample_class",
                "adata_multimap_path",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "run_id": "run-1",
                "logical_id": "logical-1",
                "sample_class": "skin",
                "adata_multimap_path": legacy_path,
            }
        )

    row_dir = tmp_path / "rows" / "run-1"
    (row_dir / "v3").mkdir(parents=True)
    (row_dir / "status.json").write_text(json.dumps({"status": "success"}), encoding="utf-8")
    pd.DataFrame(
        {
            "barcode": ["a", "b", "d"],
            "v3_unique_viral": [1, 0, 0],
            "v3_equal_viral": [1, 1, 0],
            "v3_host_conservative_viral": [1, 0, 1],
        }
    ).to_csv(row_dir / "v3" / "per_cell.tsv", sep="\t", index=False)

    output = summarize_cells(manifest, tmp_path / "rows", tmp_path / "cells.tsv")
    row = pd.read_csv(output, sep="\t").iloc[0]
    assert row["n_barcode_overlap"] == 2
    assert row["barcode_jaccard"] == 0.5
    assert row["legacy_nonzero_viral_cells"] == 2
    assert row["legacy_v3_unique_nonzero_overlap"] == 1
    assert row["legacy_v3_unique_nonzero_jaccard"] == 0.5
    assert row["legacy_v3_equal_nonzero_overlap"] == 1
    assert row["legacy_v3_equal_nonzero_jaccard"] == 1 / 3
    assert row["legacy_v3_host_conservative_nonzero_overlap"] == 1
    assert row["legacy_v3_host_conservative_nonzero_jaccard"] == 1 / 3

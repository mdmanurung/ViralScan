#!/usr/bin/env python3
"""Exercise packaged assets, doctor and validate-run on a synthetic conserved v3 fixture."""

import argparse
import subprocess
import sys
from importlib.resources import files
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse

from viralscan.run_safety import build_run_manifest, prepare_output_directory, write_run_complete


def smoke(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    for asset in ("Snakefile", "scripts/emptydrops.R", "templates/report.html.j2"):
        if not files("viralscan").joinpath(asset).is_file():
            raise RuntimeError(f"Packaged asset is missing: {asset}")
    r1, r2 = work / "R1.fastq", work / "R2.fastq"
    r1.write_text("@synthetic\nACGT\n+\nIIII\n")
    r2.write_text("@synthetic\nTGCA\n+\nIIII\n")
    args = argparse.Namespace(
        sample1=str(r1),
        sample2=str(r2),
        output=str(work / "run"),
        multimap_method="equal",
        cell_calling="none",
        called_cells_file=None,
        resume=False,
        overwrite=False,
        yes=False,
        verbose=False,
        quiet=False,
    )
    run = Path(args.output)
    prepare_output_directory(
        run, build_run_manifest(args), resume=False, overwrite=False, yes=False
    )
    target = run / "sample" / "kb-python" / "counts_unfiltered"
    target.mkdir(parents=True)
    unique = sparse.csr_matrix([[1.0]])
    adata = ad.AnnData(
        X=unique.copy(),
        obs=pd.DataFrame(index=["BC1"]),
        var=pd.DataFrame(index=["V1"]),
    )
    adata.layers["counts_unique"] = unique
    adata.layers["counts_ambiguous_allocated"] = sparse.csr_matrix((1, 1), dtype=float)
    adata.uns.update(
        {
            "count_schema_version": "3.0.0",
            "quantification_unit": "bustools-resolved-cb-umi-molecule",
            "multimap_method": "equal",
            "molecule_audit": {
                "input_molecules": 1,
                "resolved_molecules": 1,
                "unique_molecules": 1,
                "ambiguous_molecules": 0,
                "unresolved_molecules": 0,
                "ignored_read_multiplicity": 0,
                "allocated_ambiguous_mass": 0.0,
            },
        }
    )
    adata.write_h5ad(target / "adata_multimap.h5ad")
    (run / "sample" / "config.yaml").write_text("output: sample\n")
    write_run_complete(run)
    for command in (("doctor", "--profile", "pip"), ("validate-run", "--run-dir", str(run))):
        subprocess.run([sys.executable, "-m", "viralscan.menu", *command], check=True, cwd=work)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    smoke(parser.parse_args().output.resolve())

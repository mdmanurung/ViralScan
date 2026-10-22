"""Add optional scAutoQC annotations to existing ViralScan count matrices."""

import json
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd
import scanpy as sc

from viralscan.kb_outputs import KbCountOutputs
from viralscan.run_context import RunContext
from viralscan.scautoqc import add_scautoqc
from viralscan.virus_grouping import load_run_identity


def run(ctx: RunContext, done_file: str | Path) -> None:
    """Annotate existing H5AD files in place without changing their matrices."""
    config = ctx.config
    outputs = KbCountOutputs.from_config_output(config.output)
    identity = load_run_identity(config.output)
    paths = [outputs.adata]
    if config.multimapping and outputs.adata_multimap.is_file():
        paths.append(outputs.adata_multimap)
    annotation_path = Path(config.output) / "results" / "scautoqc_gene_annotations.tsv"
    gene_annotations = (
        pd.read_csv(annotation_path, sep="\t", index_col="gene_id")
        if annotation_path.is_file()
        else None
    )

    scored = {}
    for path in paths:
        adata = sc.read_h5ad(str(path))
        shape, obs_names, var_names = adata.shape, adata.obs_names.copy(), adata.var_names.copy()
        add_scautoqc(adata, identity, gene_annotations=gene_annotations)
        if gene_annotations is None:
            gene_annotations = adata.var[
                ["gene_symbol", "scautoqc_mito"]
                + [key for key in ("gene_biotype", "chromosome_name") if key in adata.var]
            ].copy()
            gene_annotations = gene_annotations.rename(columns={"scautoqc_mito": "mito"})
            gene_annotations.index = gene_annotations.index.astype(str).str.replace(
                r"\.\d+$", "", regex=True
            )
            gene_annotations.index.name = "gene_id"
            gene_annotations.to_csv(annotation_path, sep="\t")
        if (
            adata.shape != shape
            or not adata.obs_names.equals(obs_names)
            or not adata.var_names.equals(var_names)
        ):
            raise RuntimeError(
                f"scAutoQC unexpectedly changed the matrix dimensions or indices: {path}"
            )
        scored[path.name] = {"n_cells_scored": int(adata.n_obs), **adata.uns["scautoqc"]}
        fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".h5ad", dir=path.parent)
        os.close(fd)
        try:
            adata.write_h5ad(temporary)
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

    Path(done_file).write_text(json.dumps(scored, indent=2, default=str) + "\n", encoding="utf-8")


snakemake: Any

if "snakemake" in globals():
    run(
        RunContext.from_yaml(snakemake.params.configfile),  # noqa: F821
        snakemake.output[0],  # noqa: F821
    )

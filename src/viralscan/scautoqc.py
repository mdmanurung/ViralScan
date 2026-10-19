"""Optional, in-place scAutoQC annotation for ViralScan AnnData objects."""

from __future__ import annotations

import warnings
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd


def query_gene_annotations(gene_ids: Any) -> pd.DataFrame:
    """Retrieve symbols and mitochondrial membership for Ensembl IDs with Scanpy."""
    import scanpy as sc

    ids = pd.Index(gene_ids.astype(str) if hasattr(gene_ids, "astype") else map(str, gene_ids))
    stable_ids = ids.str.replace(r"\.\d+$", "", regex=True)
    is_mouse = stable_ids.str.startswith("ENSMUSG").any()
    is_human = stable_ids.str.startswith("ENSG").any()
    if is_mouse and is_human:
        raise ValueError("scAutoQC cannot annotate a mixed human/mouse Ensembl feature index.")
    if is_mouse:
        organism = "mmusculus"
    elif stable_ids.str.startswith("ENSG").any():
        organism = "hsapiens"
    else:
        raise ValueError(
            "scAutoQC gene annotation needs gene symbols or human/mouse Ensembl gene IDs."
        )
    try:
        annotations = sc.queries.biomart_annotations(
            organism,
            ["ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name"],
        )
        if "ensembl_gene_id" not in annotations or annotations.empty:
            raise ValueError("empty BioMart response")
        annotations = annotations.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
        mito_ids = set(
            sc.queries.mitochondrial_genes(organism, attrname="ensembl_gene_id")["ensembl_gene_id"]
            .dropna()
            .astype(str)
        )
    except Exception as exc:  # network, BioMart outage, or empty/odd response
        raise RuntimeError(
            "scAutoQC could not retrieve Ensembl BioMart annotations (needs network access to "
            "www.ensembl.org, which serves the latest release). Provide gene symbols in "
            "var['gene_name'] or pre-seed <output>/results/scautoqc_gene_annotations.tsv "
            f"(columns gene_id, gene_symbol, mito). Cause: {exc!r}"
        ) from exc
    mapped = annotations.reindex(stable_ids)
    unmapped = float(mapped["external_gene_name"].isna().mean())
    if unmapped > 0.2:
        warnings.warn(
            f"scAutoQC: {unmapped:.0%} of features have no BioMart symbol (Ensembl release "
            "mismatch with the reference?); mitochondrial/ribosomal sets may be undercounted.",
            stacklevel=2,
        )
    return pd.DataFrame(
        {
            "gene_symbol": mapped["external_gene_name"].fillna("").astype(str).to_numpy(),
            "gene_biotype": mapped["gene_biotype"].fillna("").astype(str).to_numpy(),
            "chromosome_name": mapped["chromosome_name"].fillna("").astype(str).to_numpy(),
            "mito": stable_ids.isin(mito_ids),
        },
        index=stable_ids,
    )


def add_scautoqc(
    adata: Any, identity: Any = None, *, gene_annotations: Any = None, random_state: int = 0
) -> None:
    """Add SCTK cell and cluster QC flags/coordinates to ``adata`` in place.

    QC metrics use host genes only.  Counts, rows, and the input matrix are never
    filtered or rewritten.  ``identity`` is ViralScan's optional gene identity table.
    """
    try:
        import anndata as ad
        import sctk
    except ImportError as exc:
        raise RuntimeError(
            "scAutoQC requires the optional full-workflow dependency sctk; "
            "install it with `pip install 'viralscan[scautoqc]'`."
        ) from exc

    if "gene_role" in adata.var:
        host = adata.var["gene_role"].astype(str).eq("host").to_numpy()
    elif "is_viral" in adata.var:
        host = ~adata.var["is_viral"].astype(bool).to_numpy()
    elif identity is not None:
        by_gene = identity.by_gene()
        host = np.asarray(
            [
                not by_gene.get(str(gene), None).viral if str(gene) in by_gene else False
                for gene in adata.var_names
            ],
            dtype=bool,
        )
    else:
        raise ValueError(
            "scAutoQC needs gene_role/is_viral annotations or ViralScan gene identity."
        )
    if not host.any():
        raise ValueError("scAutoQC found no host genes in the AnnData feature index.")

    symbol_key = next(
        (key for key in ("gene_name", "gene_symbol", "external_gene_name") if key in adata.var),
        None,
    )
    annotation_source = "existing gene symbols"
    if gene_annotations is not None:
        annotation_source = "provided gene annotation table"
        annotation = gene_annotations.reindex(
            adata.var_names.astype(str).str.replace(r"\.\d+$", "", regex=True)
        )
        adata.var["gene_symbol"] = annotation["gene_symbol"].fillna("").astype(str).to_numpy()
        if "mito" in annotation:
            adata.var["scautoqc_mito"] = (
                annotation["mito"].fillna(False).astype(bool).to_numpy() & host
            )
        for key in ("gene_biotype", "chromosome_name"):
            if key in annotation:
                adata.var[key] = annotation[key].fillna("").astype(str).to_numpy()
        symbol_key = "gene_symbol"
    elif symbol_key is None:
        annotation_source = "Scanpy Ensembl BioMart"
        gene_annotations = query_gene_annotations(adata.var_names)
        annotation = gene_annotations.reindex(
            adata.var_names.astype(str).str.replace(r"\.\d+$", "", regex=True)
        )
        adata.var["gene_symbol"] = annotation["gene_symbol"].fillna("").astype(str).to_numpy()
        adata.var["scautoqc_mito"] = annotation["mito"].fillna(False).astype(bool).to_numpy() & host
        for key in ("gene_biotype", "chromosome_name"):
            if key in annotation:
                adata.var[key] = annotation[key].fillna("").astype(str).to_numpy()
        symbol_key = "gene_symbol"
    elif symbol_key is not None:
        values = adata.var[symbol_key].fillna("").astype(str)
        if values[host].str.match(r"^ENS[A-Z]*G\d+(?:\.\d+)?$").mean() > 0.5:
            annotation_source = "Scanpy Ensembl BioMart"
            gene_annotations = query_gene_annotations(adata.var_names)
            annotation = gene_annotations.reindex(
                adata.var_names.astype(str).str.replace(r"\.\d+$", "", regex=True)
            )
            adata.var["gene_symbol"] = annotation["gene_symbol"].fillna("").astype(str).to_numpy()
            adata.var["scautoqc_mito"] = (
                annotation["mito"].fillna(False).astype(bool).to_numpy() & host
            )
            for key in ("gene_biotype", "chromosome_name"):
                if key in annotation:
                    adata.var[key] = annotation[key].fillna("").astype(str).to_numpy()
            symbol_key = "gene_symbol"
    if symbol_key is None:
        raise ValueError("scAutoQC needs gene symbols in var['gene_name'] or var['gene_symbol'].")
    symbols = adata.var[symbol_key].fillna("").astype(str)
    if not symbols[host].ne("").any():
        raise ValueError(f"scAutoQC found no gene symbols in var[{symbol_key!r}].")
    adata.var["gene_symbol"] = symbols.to_numpy()

    mito = (
        adata.var["scautoqc_mito"].astype(bool).to_numpy()
        if "scautoqc_mito" in adata.var
        else symbols.str.upper().str.match(r"^MT-", na=False).to_numpy() & host
    )
    marker_symbols = symbols.str.upper()
    ribo = marker_symbols.str.match(r"^RP[LS]", na=False).to_numpy() & host
    # Curated globins avoid broad ^HB matches such as HBEGF, HBP1, and HBS1L.
    hb = (
        marker_symbols.str.match(
            r"^(?:HBA|HBA[12]|HBB|HBD|HBE1|HBG[12]|HBM|HBQ1|HBZ)(?:P\d+)?$", na=False
        ).to_numpy()
        & host
    )
    for key, values in (("mito", mito), ("ribo", ribo), ("hb", hb)):
        if not values.any():
            raise ValueError(
                f"scAutoQC could not identify any host {key} features from {symbol_key}."
            )
        adata.var[f"scautoqc_{key}"] = values

    # SCTK mutates its AnnData argument. Calculate host-only metrics in a temporary
    # object, then transfer only metadata onto the original complete matrix.
    host_adata = adata[:, host].copy()
    host_adata.var["mito"] = mito[host]
    host_adata.var["ribo"] = ribo[host]
    host_adata.var["hb"] = hb[host]
    sctk.calculate_qc(host_adata, flags={"mito": None, "ribo": None, "hb": None})
    rng_state = np.random.get_state()
    try:
        np.random.seed(random_state)
        try:
            sctk.cellwise_qc(
                host_adata,
                metrics=["n_counts", "n_genes", "percent_mito", "percent_ribo", "percent_hb"],
                cell_qc_key="scautoqc_cell_passed_qc",
            )
        except TypeError as exc:
            # SCTK 0.2.0's all-fail fallback indexes list metrics by name.
            if "list indices must be integers or slices, not str" not in str(exc) or (
                "scautoqc_cell_passed_qc" not in host_adata.obs
                or host_adata.obs["scautoqc_cell_passed_qc"].any()
            ):
                raise
    finally:
        np.random.set_state(rng_state)

    metric_names = [
        "n_counts",
        "n_genes",
        "log1p_n_counts",
        "log1p_n_genes",
        "percent_mito",
        "percent_ribo",
        "percent_hb",
    ]
    for name in metric_names:
        adata.obs[f"scautoqc_{name}"] = host_adata.obs[name].to_numpy()
    adata.obs["scautoqc_cell_passed_qc"] = host_adata.obs["scautoqc_cell_passed_qc"].to_numpy()

    # SCTK's clustering helpers operate on the QC metrics, not expression counts.
    qc_metrics = [
        "scautoqc_log1p_n_counts",
        "scautoqc_log1p_n_genes",
        "scautoqc_percent_mito",
        "scautoqc_percent_ribo",
        "scautoqc_percent_hb",
    ]
    cluster_status = "complete"
    if adata.n_obs >= 6:
        qc_adata = ad.AnnData(X=adata.obs[qc_metrics].to_numpy(), obs=adata.obs.copy())
        sctk.generate_qc_clusters(qc_adata, metrics=qc_metrics)
        sctk.clusterwise_qc(
            qc_adata,
            cell_qc_key="scautoqc_cell_passed_qc",
            key_added="scautoqc_cluster_passed_qc",
        )
        adata.obs["scautoqc_qc_cluster"] = qc_adata.obs["qc_cluster"].to_numpy()
        adata.obs["scautoqc_cluster_passed_qc"] = qc_adata.obs[
            "scautoqc_cluster_passed_qc"
        ].to_numpy()
        adata.obsm["scautoqc_X_umap_qc"] = qc_adata.obsm["X_umap_qc"].copy()
    else:
        cluster_status = "skipped: fewer than six cells"
        adata.obs["scautoqc_qc_cluster"] = "unavailable"
        adata.obs["scautoqc_cluster_passed_qc"] = pd.Series(
            pd.NA, index=adata.obs_names, dtype="boolean"
        )
    adata.uns["scautoqc"] = {
        "method": "sctk.calculate_qc + cellwise_qc + generate_qc_clusters + clusterwise_qc",
        "sctk_version": version("sctk"),
        "gene_symbol_key": symbol_key,
        "gene_annotation_source": annotation_source,
        "host_features": int(host.sum()),
        "qc_metrics": ["n_counts", "n_genes", "percent_mito", "percent_ribo", "percent_hb"],
        "cell_removal_performed": False,
        # Known limitation: scautoqc runs before cell calling, so this is every barcode in the
        # kb H5AD; empty droplets shape the SCTK distributions. See PLAN SCAUTOQC-01.
        "cell_set": "all barcodes in the input H5AD (before cell calling)",
        "n_cells_scored": int(adata.n_obs),
        "clusterwise_status": cluster_status,
        "seed": int(random_state),
    }

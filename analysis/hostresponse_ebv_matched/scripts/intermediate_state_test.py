#!/usr/bin/env python3
"""Depth-aware test for an intermediate EBV host-cell state (todo ideas 1a/4a/4b).

Four disciplinary lenses converged on a hypothesis: the classifier's specificity gap
reflects a real 'primed/intermediate' cell state between EBV-latent and EBV-lytic.
But we now know that gap is largely a sequencing-depth artifact. This test asks the
honest question: is there a distinct intermediate host-transcriptome state that is NOT
merely a depth stratum?

Two passes:
  A. ALL cells — UMAP + Leiden on host HVGs; per-cluster EBV+ rate, EBV UMI, and DEPTH.
     If clusters are organized by depth (a cluster = a depth bin), that is a QC artifact.
  B. DEPTH-MATCHED cells (depth controlled by design) — same clustering; a distinct
     intermediate-EBV cluster here would be genuinely depth-independent.

Writes a UMAP figure (colored by EBV UMI, depth, cluster) and a cluster table.
Run: PYTHONPATH=src <bioenv-python> analysis/hostresponse_ebv_matched/scripts/intermediate_state_test.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, kruskal

BASE = "results/hostresponse_ebv_matched"
FIGDIR = "analysis/hostresponse_ebv_matched/outputs/figures"
THRESH = 10
N_BINS = 20
SEED = 42


def _sum(X):
    return np.asarray(X.sum(1)).ravel() if sp.issparse(X) else np.asarray(X).sum(1).ravel()


def prep():
    host = ad.read_h5ad(f"{BASE}/host_only_matched.h5ad")
    ebv = ad.read_h5ad(f"{BASE}/ebv_burden_matched.h5ad")
    sh = sorted(set(host.obs_names) & set(ebv.obs_names))
    host, ebv = host[sh].copy(), ebv[sh].copy()
    host.obs["ebv_umi"] = _sum(ebv.X)
    host.obs["depth"] = _sum(host.X)
    host.obs["ebv_pos"] = (host.obs["ebv_umi"] >= THRESH).astype(int)
    sc.pp.normalize_total(host, target_sum=1e4)
    sc.pp.log1p(host)
    sc.pp.highly_variable_genes(host, flavor="seurat")
    host = host[:, host.var["highly_variable"]].copy()
    sc.pp.scale(host, max_value=10)
    sc.tl.pca(host, n_comps=30, random_state=SEED)
    return host


def cluster_embed(adata, tag, out):
    from sklearn.cluster import HDBSCAN, KMeans  # noqa: PLC0415
    from sklearn.metrics import silhouette_score  # noqa: PLC0415
    sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    Xp = adata.obsm["X_pca"]
    df = adata.obs
    out.append(f"\n=== {tag}: {adata.n_obs} cells ===")
    out.append(f"  Spearman(depth, ebv_umi) in this set = {spearmanr(df['depth'], df['ebv_umi']).correlation:.3f}")

    # (1) Density-based: does discrete structure exist at all?
    lab = HDBSCAN(min_cluster_size=30).fit_predict(Xp)
    n_hdb = len(set(lab)) - (1 if -1 in lab else 0)
    out.append(f"  HDBSCAN(min_cluster_size=30): {n_hdb} clusters, {(lab == -1).sum()}/{len(lab)} unassigned "
               f"-> {'CONTINUUM (no discrete states)' if n_hdb == 0 else 'discrete structure present'}")

    # (2) Silhouette: how discrete is a forced k-partition (higher = more separated; <0.1 ~ no structure)
    sils = {k: silhouette_score(Xp, KMeans(k, random_state=SEED, n_init=10).fit_predict(Xp)) for k in (2, 3, 4)}
    out.append("  KMeans silhouette (structure strength): " + ", ".join(f"k={k}:{s:.3f}" for k, s in sils.items()))

    # (3) Direct tristable test: force k=3, is there a middle-EBV cluster NOT explained by depth?
    km3 = KMeans(3, random_state=SEED, n_init=10).fit_predict(Xp)
    adata.obs["cluster"] = pd.Categorical([str(c) for c in km3])
    g = df.groupby("cluster", observed=True)
    tab = pd.DataFrame({
        "n": g.size(), "ebv_pos_rate": g["ebv_pos"].mean(),
        "median_ebv_umi": g["ebv_umi"].median(), "median_depth": g["depth"].median(),
    }).sort_values("ebv_pos_rate")
    out.append("  forced KMeans k=3 (EBV vs depth per cluster):")
    out.append("    " + tab.round(2).to_string().replace("\n", "\n    "))
    return tab


def figure(adata, path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    um = adata.obsm["X_umap"]
    for a, key, title in zip(axes, ["ebv_umi", "depth", "cluster"],
                             ["EBV burden (log1p UMI)", "sequencing depth (log1p)", "HDBSCAN clusters"]):
        if key == "cluster":
            for c in adata.obs["cluster"].cat.categories:
                m = (adata.obs["cluster"] == c).to_numpy()
                a.scatter(um[m, 0], um[m, 1], s=3, label=c)
            a.legend(markerscale=3, fontsize=6, ncol=2, loc="best")
        else:
            col = np.log1p(adata.obs[key].to_numpy())
            sctr = a.scatter(um[:, 0], um[:, 1], s=3, c=col, cmap="viridis")
            fig.colorbar(sctr, ax=a, shrink=0.7)
        a.set_title(title)
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def main():
    out = []
    host = prep()
    cluster_embed(host, "ALL CELLS", out)
    figure(host, f"{FIGDIR}/intermediate_umap_allcells.pdf")

    # depth-matched subset (coarsened exact matching on depth quantile bins)
    rng = np.random.default_rng(SEED)
    depth = host.obs["depth"].to_numpy(); y = host.obs["ebv_pos"].to_numpy()
    edges = np.quantile(depth, np.linspace(0, 1, N_BINS + 1)); edges[-1] += 1
    binid = np.digitize(depth, edges[1:-1]); keep = []
    for b in np.unique(binid):
        idx = np.where(binid == b)[0]; pos, neg = idx[y[idx] == 1], idx[y[idx] == 0]
        k = min(len(pos), len(neg))
        if k:
            keep += rng.choice(pos, k, replace=False).tolist() + rng.choice(neg, k, replace=False).tolist()
    hm = host[np.array(sorted(keep))].copy()
    cluster_embed(hm, "DEPTH-MATCHED", out)
    figure(hm, f"{FIGDIR}/intermediate_umap_depthmatched.pdf")

    out.append("\nVerdict rule: if clusters track depth (high depth-Kruskal, EBV+ rate rises with "
               "median_depth) the 'intermediate' is a depth stratum; if a distinct intermediate-EBV "
               "cluster persists in the DEPTH-MATCHED set, it is a genuine depth-independent state.")
    report = "\n".join(out)
    print(report)
    with open(f"{BASE}/intermediate_state_test.txt", "w") as fh:
        fh.write(report + "\n")


if __name__ == "__main__":
    main()

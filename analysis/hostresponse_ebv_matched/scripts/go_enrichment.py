#!/usr/bin/env python3
"""Powered GO enrichment on the FULL depth-robust EBV gene set, with a mito control.

The 5 top-stability depth-robust genes were too few for a statistical GO test. Here we
build the full depth-robust set: for every highly-variable host gene, the depth-adjusted
partial association with EBV status (residualize gene and label on log host-depth, then
correlate), BH-FDR corrected. Genes with FDR<0.05 that also survive an ADDITIONAL control
for percent-mitochondrial content form the depth-and-mito-robust set. We map Ensembl IDs
to symbols (mygene.info) and run GO Biological Process enrichment (Enrichr).

Run: PYTHONPATH=src <bioenv-python> analysis/hostresponse_ebv_matched/scripts/go_enrichment.py
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = "results/hostresponse_ebv_matched"
THRESH = 10
FDR = 0.05


def _sum(X):
    return np.asarray(X.sum(1)).ravel() if sp.issparse(X) else np.asarray(X).sum(1).ravel()


def resid(Y, C):
    """Residualize columns of Y (cells x genes) on covariates C (cells x k, no intercept)."""
    Cc = np.column_stack([np.ones(len(C)), C])
    beta, *_ = np.linalg.lstsq(Cc, Y, rcond=None)
    return Y - Cc @ beta


def partial_assoc(X, y, C):
    """Partial correlation of each gene (cols of X) with y, adjusting for covariates C."""
    Xr = resid(X, C)
    yr = resid(y.reshape(-1, 1).astype(float), C).ravel()
    Xr -= Xr.mean(0); Xr /= (Xr.std(0) + 1e-9)
    yr = (yr - yr.mean()) / (yr.std() + 1e-9)
    r = (Xr * yr[:, None]).mean(0)
    n = len(y)
    t = r * np.sqrt((n - 2) / np.maximum(1 - r**2, 1e-12))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    return r, p


def map_symbols(ensembl_ids):
    ids = [e.split(".")[0] for e in ensembl_ids]
    data = urllib.parse.urlencode({"q": ",".join(ids), "scopes": "ensembl.gene",
                                   "fields": "symbol", "species": "human"}).encode()
    req = urllib.request.Request("https://mygene.info/v3/query", data=data)
    res = json.load(urllib.request.urlopen(req, timeout=60))
    m = {}
    for r in res:
        if "symbol" in r:
            m[r["query"]] = r["symbol"]
    return [m[i] for i in ids if i in m]


def enrichr(genes, library="GO_Biological_Process_2021", topn=8):
    boundary = "----------geneset"
    body = (f"--{boundary}\r\nContent-Disposition: form-data; name=\"list\"\r\n\r\n"
            f"{chr(10).join(genes)}\r\n--{boundary}\r\n"
            f"Content-Disposition: form-data; name=\"description\"\r\n\r\ndepth_robust\r\n--{boundary}--\r\n")
    req = urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=body.encode(),
                                 headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    uid = json.load(urllib.request.urlopen(req, timeout=60))["userListId"]
    time.sleep(1)
    url = f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={library}"
    res = json.load(urllib.request.urlopen(url, timeout=60))[library]
    # columns: rank, term, pval, zscore, combined, genes, adj_pval, ...
    return [(r[1], r[2], r[6], len(r[5])) for r in res[:topn]]


def main():
    host = ad.read_h5ad(f"{BASE}/host_only_matched.h5ad")
    ebv = ad.read_h5ad(f"{BASE}/ebv_burden_matched.h5ad")
    sh = sorted(set(host.obs_names) & set(ebv.obs_names))
    host, ebv = host[sh].copy(), ebv[sh].copy()
    y = (_sum(ebv.X) >= THRESH).astype(int)
    depth = _sum(host.X)
    host.var["mt"] = host.var_names.str.startswith("MT-") | host.var_names.str.startswith("ENSG00000198")
    # pct mito by known mitochondrial Ensembl IDs (MT-* symbols are ENSG000002.. ; use a curated set)
    mt_ensembl = {"ENSG00000198899","ENSG00000198804","ENSG00000198712","ENSG00000228253","ENSG00000198763",
                  "ENSG00000198938","ENSG00000198840","ENSG00000212907","ENSG00000198886","ENSG00000198786",
                  "ENSG00000198695","ENSG00000198727"}
    ids0 = [v.split(".")[0] for v in host.var_names]
    mt_mask = np.array([i in mt_ensembl for i in ids0])
    tot = _sum(host.X)
    mt_counts = _sum(host[:, mt_mask].X) if mt_mask.any() else np.zeros(len(y))
    pct_mt = mt_counts / tot * 100

    sc.pp.normalize_total(host, target_sum=1e4); sc.pp.log1p(host)
    sc.pp.highly_variable_genes(host, flavor="seurat")
    hv = host[:, host.var["highly_variable"]].copy()
    X = hv.X.toarray() if sp.issparse(hv.X) else np.asarray(hv.X)
    genes = np.array(hv.var_names)
    logd = np.log1p(depth)

    out = []
    # depth-adjusted association
    r_d, p_d = partial_assoc(X, y, logd.reshape(-1, 1))
    fdr_d = multipletests(p_d, method="fdr_bh")[1]
    sig_d = fdr_d < FDR
    # + mito-adjusted
    r_m, p_m = partial_assoc(X, y, np.column_stack([logd, pct_mt]))
    fdr_m = multipletests(p_m, method="fdr_bh")[1]
    sig_dm = (fdr_d < FDR) & (fdr_m < FDR)

    out.append(f"HVG candidates: {len(genes)} | pct-mito median {np.median(pct_mt):.1f}%")
    out.append(f"depth-adjusted FDR<{FDR}: {int(sig_d.sum())} genes "
               f"({int((sig_d&(r_d>0)).sum())} up / {int((sig_d&(r_d<0)).sum())} down)")
    out.append(f"depth+mito-adjusted FDR<{FDR}: {int(sig_dm.sum())} genes "
               f"-> {int((sig_d&~sig_dm).sum())} lost when also controlling %mito")
    # is MT-ND4L (ENSG00000212907) affected?
    for probe in ["ENSG00000212907"]:  # MT-ND4L
        j = [k for k, g in enumerate(genes) if g.startswith(probe)]
        if j:
            k = j[0]
            out.append(f"  MT-ND4L: depth-adj FDR={fdr_d[k]:.1e} -> depth+mito-adj FDR={fdr_m[k]:.1e} "
                       f"({'SURVIVES' if sig_dm[k] else 'LOST — likely a mito-QC effect'})")

    up = genes[sig_dm & (r_m > 0)]; dn = genes[sig_dm & (r_m < 0)]
    for label, gs in [("UP in EBV+ (depth+mito robust)", up), ("DOWN in EBV+", dn)]:
        syms = map_symbols(list(gs)) if len(gs) else []
        out.append(f"\n{label}: {len(gs)} genes; e.g. {', '.join(syms[:12])}")
        if len(syms) >= 5:
            out.append("  Top GO Biological Process (Enrichr):")
            for term, p, padj, k in enrichr(syms):
                out.append(f"    {term[:62]:<62} p={p:.1e} FDR={padj:.1e} (n={k})")
        else:
            out.append("  (too few genes for enrichment)")

    report = "\n".join(out)
    print(report)
    with open(f"{BASE}/go_enrichment.txt", "w") as fh:
        fh.write(report + "\n")


if __name__ == "__main__":
    main()

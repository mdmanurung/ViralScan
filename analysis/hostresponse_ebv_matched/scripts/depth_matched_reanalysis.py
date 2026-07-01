#!/usr/bin/env python3
"""Depth-matched re-analysis of the EBV host-response signal (the principled fix).

The raw AUC 0.866 is confounded because the >=10-UMI EBV label tracks sequencing
depth (see depth_confounder_check.py). Here we remove the confound by DESIGN: build a
depth-matched case-control set in which EBV+ and EBV- cells have (near-)identical depth
distributions, so depth alone can no longer discriminate. Then we ask whether the host
transcriptome STILL predicts EBV status. If yes -> a real, depth-independent signal; if
it collapses to chance -> the original result was a depth artifact.

Matching: coarsened exact matching on host-depth quantile bins; within each bin keep
min(n_pos, n_neg) of each class. Reports depth-alone AUC in the matched set (should be
~0.5 if matching worked), then 5-fold CV for the 15 stable genes, the 5 depth-robust
genes, and (as a positive control) depth alone.

Run: PYTHONPATH=src <bioenv-python> analysis/hostresponse_ebv_matched/scripts/depth_matched_reanalysis.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, matthews_corrcoef
from sklearn.preprocessing import StandardScaler

BASE = "results/hostresponse_ebv_matched"
THRESH = 10
N_BINS = 20
SEED = 42
E_ROBUST = 2.0  # E-value cutoff for "depth-robust" genes


def _sum(X):
    return np.asarray(X.sum(1)).ravel() if sp.issparse(X) else np.asarray(X).sum(1).ravel()


def load():
    host = ad.read_h5ad(f"{BASE}/host_only_matched.h5ad")
    ebv = ad.read_h5ad(f"{BASE}/ebv_burden_matched.h5ad")
    shared = sorted(set(host.obs_names) & set(ebv.obs_names))
    host, ebv = host[shared].copy(), ebv[shared].copy()
    depth = _sum(host.X)                        # host-only raw depth
    y = (_sum(ebv.X) >= THRESH).astype(int)
    sc.pp.normalize_total(host, target_sum=1e4)
    sc.pp.log1p(host)
    stab = pd.read_csv(f"{BASE}/Epstein-Barr_virus_stability.csv")
    genes15 = [g for g in stab.loc[stab["stab_prob"] >= 0.6, "gene"] if g in set(host.var_names)]
    # depth-robust genes from the confounder check (E-value >= 2)
    ev = pd.read_csv(f"{BASE}/depth_confounder_gene_evalues.csv")
    genes_rob = [g for g in ev.loc[ev["E_value"] >= E_ROBUST, "gene"] if g in set(host.var_names)]
    Z = host[:, genes15].X
    Z = Z.toarray() if sp.issparse(Z) else np.asarray(Z)
    Zr = host[:, genes_rob].X
    Zr = Zr.toarray() if sp.issparse(Zr) else np.asarray(Zr)
    return y, depth, Z, genes15, Zr, genes_rob


def depth_match(y, depth, seed=SEED):
    """Coarsened exact matching on host-depth quantile bins; equal n per class per bin."""
    rng = np.random.default_rng(seed)
    edges = np.quantile(depth, np.linspace(0, 1, N_BINS + 1))
    edges[-1] += 1
    binid = np.digitize(depth, edges[1:-1])
    keep = []
    for b in np.unique(binid):
        idx = np.where(binid == b)[0]
        pos, neg = idx[y[idx] == 1], idx[y[idx] == 0]
        k = min(len(pos), len(neg))
        if k == 0:
            continue
        keep.extend(rng.choice(pos, k, replace=False).tolist())
        keep.extend(rng.choice(neg, k, replace=False).tolist())
    return np.array(sorted(keep))


def cv(X, y, seed=SEED):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    p = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        s = StandardScaler().fit(X[tr])
        m = LogisticRegression(max_iter=1000, C=1.0).fit(s.transform(X[tr]), y[tr])
        p[te] = m.predict_proba(s.transform(X[te]))[:, 1]
    return roc_auc_score(y, p), matthews_corrcoef(y, (p >= 0.5).astype(int))


def main():
    y, depth, Z, genes15, Zr, genes_rob = load()
    m = depth_match(y, depth)
    ym, dm, Zm, Zrm = y[m], depth[m], Z[m], Zr[m]
    logdm = np.log1p(dm).reshape(-1, 1)

    out = []
    out.append(f"Full set: n={len(y)} ({int(y.sum())}+/{int((1-y).sum())}-)")
    out.append(f"Depth-matched set: n={len(m)} ({int(ym.sum())}+/{int((1-ym).sum())}-)  [{N_BINS} depth bins, equal n/class/bin]")
    # matching quality: class depth medians + Mann-Whitney
    md_pos, md_neg = np.median(dm[ym == 1]), np.median(dm[ym == 0])
    pmw = mannwhitneyu(dm[ym == 1], dm[ym == 0]).pvalue
    out.append(f"  matched class depth medians: EBV+ {md_pos:.0f} vs EBV- {md_neg:.0f}  (Mann-Whitney p={pmw:.2f})")

    auc_d, mcc_d = cv(logdm, ym)
    out.append(f"\n[control] depth-ALONE AUC in matched set = {auc_d:.3f} (MCC {mcc_d:.3f})  <- should be ~0.5 if matching worked")

    auc15, mcc15 = cv(Zm, ym)
    out.append(f"\n[15 stable genes]  matched-set AUC = {auc15:.3f} (MCC {mcc15:.3f})   (confounded all-cell AUC was 0.793; raw-design headline 0.866)")
    aucr, mccr = cv(Zrm, ym)
    out.append(f"[{len(genes_rob)} depth-robust genes] matched-set AUC = {aucr:.3f} (MCC {mccr:.3f})")
    out.append(f"  depth-robust genes: {', '.join(genes_rob)}")

    out.append("\nInterpretation: AUC well above the depth-alone control at matched depth = a real, "
               "depth-independent host-response signal; near the control = the original result was a depth artifact.")

    report = "\n".join(out)
    print(report)
    with open(f"{BASE}/depth_matched_reanalysis.txt", "w") as fh:
        fh.write(report + "\n")


if __name__ == "__main__":
    main()

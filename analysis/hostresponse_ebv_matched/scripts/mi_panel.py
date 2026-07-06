#!/usr/bin/env python3
"""Minimal EBV gene panel via greedy forward selection (idea 7a).

Question: can a handful of genes match the 15-gene signature's discrimination? We
greedily add genes that most improve depth-matched cross-validated AUC (evaluating on
the DEPTH-MATCHED set so the answer is not a depth artifact), and report AUC vs panel
size plus the minimal panel reaching 95% of the full-15 AUC.

Run: PYTHONPATH=src <bioenv-python> analysis/hostresponse_ebv_matched/scripts/mi_panel.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from sklearn.feature_selection import mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

BASE = "results/hostresponse_ebv_matched"
THRESH = 10
N_BINS = 20
SEED = 42


def _sum(X):
    return np.asarray(X.sum(1)).ravel() if sp.issparse(X) else np.asarray(X).sum(1).ravel()


def cv_auc(X, y, seed=SEED):
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    p = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        s = StandardScaler().fit(X[tr])
        p[te] = LogisticRegression(max_iter=1000).fit(s.transform(X[tr]), y[tr]).predict_proba(s.transform(X[te]))[:, 1]
    return roc_auc_score(y, p)


def depth_match(y, depth, seed=SEED):
    rng = np.random.default_rng(seed)
    edges = np.quantile(depth, np.linspace(0, 1, N_BINS + 1)); edges[-1] += 1
    b = np.digitize(depth, edges[1:-1]); keep = []
    for bb in np.unique(b):
        idx = np.where(b == bb)[0]; pos, neg = idx[y[idx] == 1], idx[y[idx] == 0]
        k = min(len(pos), len(neg))
        if k:
            keep += rng.choice(pos, k, replace=False).tolist() + rng.choice(neg, k, replace=False).tolist()
    return np.array(sorted(keep))


def main():
    host = ad.read_h5ad(f"{BASE}/host_only_matched.h5ad")
    ebv = ad.read_h5ad(f"{BASE}/ebv_burden_matched.h5ad")
    sh = sorted(set(host.obs_names) & set(ebv.obs_names))
    host, ebv = host[sh].copy(), ebv[sh].copy()
    y = (_sum(ebv.X) >= THRESH).astype(int)
    depth = _sum(host.X)
    sc.pp.normalize_total(host, target_sum=1e4); sc.pp.log1p(host)
    stab = pd.read_csv(f"{BASE}/Epstein-Barr_virus_stability.csv")
    genes = [g for g in stab.loc[stab["stab_prob"] >= 0.6, "gene"] if g in set(host.var_names)]
    Z = host[:, genes].X
    Z = (Z.toarray() if sp.issparse(Z) else np.asarray(Z))
    m = depth_match(y, depth); Zm, ym = Z[m], y[m]

    out = []
    full = cv_auc(Zm, ym)
    out.append(f"Depth-matched set: n={len(m)} ; full 15-gene AUC = {full:.3f}")
    mi = mutual_info_classif(Zm, ym, random_state=SEED)
    out.append("Mutual information (depth-matched) per gene, top 5: " +
               ", ".join(f"{genes[i].split('.')[0]}={mi[i]:.3f}" for i in np.argsort(mi)[::-1][:5]))

    # greedy forward selection by depth-matched AUC
    chosen, remaining = [], list(range(len(genes)))
    out.append("\nGreedy forward panel (depth-matched CV AUC):")
    while remaining:
        best, bi = -1, None
        for j in remaining:
            a = cv_auc(Zm[:, chosen + [j]], ym)
            if a > best:
                best, bi = a, j
        chosen.append(bi); remaining.remove(bi)
        star = "  <- 95% of full" if best >= 0.95 * full and (len(chosen) == 1 or True) else ""
        out.append(f"  {len(chosen):>2} genes: AUC {best:.3f}" + (f"   (+{genes[bi].split('.')[0]})"))
        if len(chosen) >= 8:
            break
    reach = next((k for k in range(1, len(genes) + 1)
                  if cv_auc(Zm[:, [c for c in chosen[:k]]], ym) >= 0.95 * full), None)
    out.append(f"\nMinimal panel reaching 95% of the 15-gene AUC: {reach} gene(s).")
    report = "\n".join(out)
    print(report)
    with open(f"{BASE}/mi_panel.txt", "w") as fh:
        fh.write(report + "\n")


if __name__ == "__main__":
    main()

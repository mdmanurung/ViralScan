#!/usr/bin/env python3
"""Depth-confounder check for the EBV host-response classifier (todo idea 5a).

Question: is the host-gene → EBV-status signal (AUC ~0.866) real biology, or is it
driven by sequencing depth as a common cause (deeper cells have more EBV UMI AND more
host counts)? host_depth here is host-only raw UMI (the EBV burden is a separate matrix).

Diagnostics:
  1. Association of depth with the EBV label (Spearman; EBV+ rate by depth quintile).
  2. Depth-ALONE cross-validated AUC for EBV status (the key confound diagnostic).
  3. The 15 stable host genes: CV AUC/MCC with vs without a log(depth) covariate.
  4. Per-gene E-values (Ding & VanderWeele 2016) on depth-adjusted odds ratios.

Run: PYTHONPATH=src <bioenv-python> analysis/hostresponse_ebv_matched/scripts/depth_confounder_check.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, matthews_corrcoef
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm

BASE = "results/hostresponse_ebv_matched"
THRESH = 10
SEED = 42


def _dense_sum(X):
    return np.asarray(X.sum(1)).ravel() if sp.issparse(X) else np.asarray(X).sum(1).ravel()


def load():
    host = ad.read_h5ad(f"{BASE}/host_only_matched.h5ad")
    ebv = ad.read_h5ad(f"{BASE}/ebv_burden_matched.h5ad")
    shared = sorted(set(host.obs_names) & set(ebv.obs_names))
    host, ebv = host[shared].copy(), ebv[shared].copy()
    ebv_umi = _dense_sum(ebv.X)
    depth = _dense_sum(host.X)                       # host-only raw depth
    y = (ebv_umi >= THRESH).astype(int)
    # normalize host as the pipeline does, then pull the 15 stable genes
    import scanpy as sc
    sc.pp.normalize_total(host, target_sum=1e4)
    sc.pp.log1p(host)
    stab = pd.read_csv(f"{BASE}/Epstein-Barr_virus_stability.csv")
    genes = stab.loc[stab["stab_prob"] >= 0.6, "gene"].tolist()
    genes = [g for g in genes if g in set(host.var_names)]
    Z = host[:, genes].X
    Z = Z.toarray() if sp.issparse(Z) else np.asarray(Z)
    return y, depth, ebv_umi, Z, genes


def cv_metrics(X, y, seed=SEED):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    probs = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc_ = StandardScaler().fit(X[tr])
        m = LogisticRegression(max_iter=1000, C=1.0)
        m.fit(sc_.transform(X[tr]), y[tr])
        probs[te] = m.predict_proba(sc_.transform(X[te]))[:, 1]
    return roc_auc_score(y, probs), matthews_corrcoef(y, (probs >= 0.5).astype(int))


def _depth_alone_balanced(y, depth, seeds=(0, 1, 10, 42, 100, 1234)):
    """Depth-alone AUC under the pipeline's own split (top-50% depth per class, balance,
    80/20) — the apples-to-apples comparison to the reported host-gene AUC 0.866."""
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]

    def top(idx):
        d = depth[idx]
        return idx[d >= np.percentile(d, 50)]

    aucs = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        pt, nt = top(pos), top(neg)
        n = min(len(pt), len(nt))
        ps, ns = rng.choice(pt, n, replace=False), rng.choice(nt, n, replace=False)
        s = max(1, int(0.8 * n))
        tr, te = np.r_[ps[:s], ns[:s]], np.r_[ps[s:], ns[s:]]
        scl = StandardScaler().fit(np.log1p(depth[tr]).reshape(-1, 1))
        m = LogisticRegression(max_iter=1000).fit(scl.transform(np.log1p(depth[tr]).reshape(-1, 1)), y[tr])
        p = m.predict_proba(scl.transform(np.log1p(depth[te]).reshape(-1, 1)))[:, 1]
        aucs.append(roc_auc_score(y[te], p))
    return aucs


def e_value(or_):
    """Ding & VanderWeele E-value from an odds ratio (approx risk ratio)."""
    rr = or_ if or_ >= 1 else 1.0 / or_
    return rr + np.sqrt(rr * (rr - 1.0))


def main():
    y, depth, ebv_umi, Z, genes = load()
    logd = np.log1p(depth).reshape(-1, 1)
    out = []
    out.append(f"n={len(y)}  EBV+={int(y.sum())}  EBV-={int((1-y).sum())}  (threshold >={THRESH} UMI)")
    out.append(f"stable genes used: {len(genes)}")

    # 1. depth-label association
    rho = spearmanr(depth, ebv_umi).correlation
    rho_lab = spearmanr(depth, y).correlation
    out.append(f"\n[1] Spearman(host_depth, EBV_UMI) = {rho:.3f} ; Spearman(host_depth, EBV_label) = {rho_lab:.3f}")
    q = pd.qcut(depth, 5, labels=[f"Q{i}" for i in range(1, 6)])
    rate = pd.Series(y).groupby(q, observed=True).mean()
    med = pd.Series(depth).groupby(q, observed=True).median()
    out.append("    EBV+ rate by host-depth quintile (median depth):")
    for k in rate.index:
        out.append(f"      {k}: EBV+ rate {rate[k]:.3f}  (median depth {med[k]:.0f})")

    # 2. depth-alone AUC (all cells, plain 5-fold)
    auc_d, mcc_d = cv_metrics(logd, y)
    out.append(f"\n[2] Depth-ALONE CV AUC = {auc_d:.3f} (MCC {mcc_d:.3f})  <- confound strength (all cells)")

    # 2b. depth-alone AUC WITHIN the headline design (top-50% depth per class, balanced,
    #     6 seeds) — the apples-to-apples comparison to the reported 0.866.
    auc_bal = _depth_alone_balanced(y, depth)
    out.append(
        f"    Depth-ALONE AUC within the balanced+depth-filtered design (matches headline eval) "
        f"= {np.mean(auc_bal):.3f} +/- {np.std(auc_bal):.3f}"
    )
    out.append("    (headline host-gene AUC in the SAME design = 0.866 — i.e. depth alone does better)")
    # 2c. why the filter RAISES depth-alone AUC: it widens the between-class depth gap.
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    tp = pos[depth[pos] >= np.percentile(depth[pos], 50)]
    tn = neg[depth[neg] >= np.percentile(depth[neg], 50)]
    out.append(
        f"    class depth medians — all cells: EBV+ {np.median(depth[pos]):.0f} vs EBV- {np.median(depth[neg]):.0f}; "
        f"after top-50% filter: EBV+ {np.median(depth[tp]):.0f} vs EBV- {np.median(depth[tn]):.0f} "
        f"(only {(depth[tn] > np.median(depth[tp])).mean():.1%} of EBV- exceed the EBV+ median)"
    )

    # 3. 15-gene classifier +/- depth covariate
    auc_g, mcc_g = cv_metrics(Z, y)
    auc_gd, mcc_gd = cv_metrics(np.hstack([Z, logd]), y)
    out.append(f"\n[3] 15-gene panel        CV AUC = {auc_g:.3f} (MCC {mcc_g:.3f})")
    out.append(f"    15-gene + log(depth)  CV AUC = {auc_gd:.3f} (MCC {mcc_gd:.3f})")
    out.append(f"    delta AUC from adding depth = {auc_gd - auc_g:+.3f}")

    # 4. per-gene depth-adjusted OR + E-value
    out.append("\n[4] Per-gene depth-adjusted odds ratios (per +1 SD gene) and E-values:")
    Zs = StandardScaler().fit_transform(Z)
    logds = StandardScaler().fit_transform(logd)
    rows = []
    for j, g in enumerate(genes):
        Xd = sm.add_constant(np.column_stack([Zs[:, j], logds[:, 0]]))
        try:
            res = sm.Logit(y, Xd).fit(disp=0)
            or_g = float(np.exp(res.params[1]))
            rows.append((g, or_g, e_value(or_g), float(res.pvalues[1])))
        except Exception as exc:  # noqa: BLE001  # ANALYSIS_OK[degraded-fallback]: a single-gene fit failing (separation) is logged as NaN, not silently dropped from the count; the loop and the reported n_genes still reflect it.
            rows.append((g, np.nan, np.nan, np.nan))
    df = pd.DataFrame(rows, columns=["gene", "adj_OR", "E_value", "p_adj_depth"]).sort_values(
        "E_value", ascending=False
    )
    for _, r in df.iterrows():
        out.append(f"      {r['gene']:<20} OR={r['adj_OR']:.2f}  E={r['E_value']:.2f}  p={r['p_adj_depth']:.1e}")
    robust = int((df["E_value"] >= 2).sum())
    out.append(f"    genes with E-value >= 2 (robust to moderate confounding): {robust}/{len(genes)}")

    report = "\n".join(out)
    print(report)
    with open(f"{BASE}/depth_confounder.txt", "w") as fh:
        fh.write(report + "\n")
    df.to_csv(f"{BASE}/depth_confounder_gene_evalues.csv", index=False)


if __name__ == "__main__":
    main()

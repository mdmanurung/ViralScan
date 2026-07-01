#!/usr/bin/env python3
"""Register the EBV host-response headline values from the summary TSV.

Run with the mycelium register_value helper on PYTHONPATH:
    PYTHONPATH=<mycelium>/skills/core/scripts python analysis/hostresponse_ebv_matched/scripts/register_report_values.py

This writes analysis/hostresponse_ebv_matched/outputs/numbers.json with real
computed_at provenance (this script), not an interactive-session <stdin> path.
"""
from __future__ import annotations

import csv
from pathlib import Path

from register_value import register_value  # provided by mycelium skills/core/scripts

NS = "hostresponse_ebv_matched"
SUMMARY = "results/hostresponse_ebv_matched/hostresponse_summary.tsv"


def main() -> int:
    row = next(csv.DictReader(open(SUMMARY), delimiter="\t"))
    prov = f"{SUMMARY}:col={{col}}"
    reg = lambda key, col, ndigits: register_value(  # noqa: E731
        key, round(float(row[col]), ndigits), namespace=NS, provenance=prov.format(col=col)
    )
    reg("auc_mean", "auc_mean", 3)
    reg("mcc_mean", "mcc_mean", 3)
    reg("mcc_sd", "mcc_sd", 3)
    reg("sensitivity_mean", "sensitivity_mean", 3)
    reg("specificity_mean", "specificity_mean", 3)
    reg("balanced_acc_mean", "balanced_acc_mean", 3)
    register_value("n_cells_matched", int(row["n_cells"]), namespace=NS, provenance=f"{SUMMARY}:col=n_cells")
    register_value("n_ebv_positive", int(row["n_positive"]), namespace=NS, provenance=f"{SUMMARY}:col=n_positive")
    register_value("n_ebv_negative", int(row["n_negative"]), namespace=NS, provenance=f"{SUMMARY}:col=n_negative")
    register_value("n_stable_genes", int(row["n_stable_genes"]), namespace=NS, provenance=f"{SUMMARY}:col=n_stable_genes")
    register_value("detection_threshold_umi", int(row["detection_threshold_umi"]), namespace=NS, provenance=f"{SUMMARY}:col=detection_threshold_umi")
    print("registered host-response report values")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Register headline values from the harmonize_2x2 benchmark comparison.

Writes analysis/reference_strategy_benchmark/outputs/numbers.json in the
mycelium register_value format: a JSON object with "namespace" and "values"
list, each entry having "key", "value", "provenance", and "computed_at".

Run from the repo root:
    python analysis/reference_strategy_benchmark/scripts/register_report_values.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
OUTPUTS = REPO / "analysis" / "reference_strategy_benchmark" / "outputs"
TSV = OUTPUTS / "harmonized_2x2_unique.tsv"
OUT_JSON = OUTPUTS / "numbers.json"
NS = "reference_strategy_benchmark"
PROV_TSV = "analysis/reference_strategy_benchmark/outputs/harmonized_2x2_unique.tsv"
COMPUTED_AT = "analysis/reference_strategy_benchmark/scripts/register_report_values.py"


def _reg(values: list, key: str, value, provenance: str) -> None:
    values.append({
        "key": key,
        "value": value,
        "provenance": provenance,
        "computed_at": COMPUTED_AT,
    })


def main() -> int:
    if not TSV.exists():
        print(f"ERROR: {TSV} not found — run harmonize_2x2.py first", file=sys.stderr)
        return 1

    rows = []
    with TSV.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    if len(rows) != 6:
        print(f"ERROR: expected 6 rows, got {len(rows)}", file=sys.stderr)
        return 1

    values: list = []

    for row in rows:
        ds = row["dataset"]
        rs = row["ref_strat"]
        tag = f"{ds}_{rs}"  # e.g. hhv6b_combined
        prov = f"{PROV_TSV}:dataset={ds},ref_strat={rs}"

        _reg(values, f"shared_anchor_n_{tag}",
             int(row["shared_anchor_n"]), prov)
        _reg(values, f"starsolo_unique_anchor_{tag}",
             round(float(row["starsolo_unique_anchor"]), 1), prov)
        _reg(values, f"viralscan_unique_anchor_{tag}",
             round(float(row["viralscan_unique_anchor"]), 1), prov)
        _reg(values, f"vs_corrected_anchor_{tag}",
             round(float(row["vs_corrected_anchor"]), 1), prov)
        _reg(values, f"multimap_gain_anchor_{tag}",
             round(float(row["multimap_gain_anchor"]), 1), prov)

    # EBV GTF artifact diagnostics (same root cause class as HSV-1)
    _reg(
        values,
        "ebv_starsolo_root_cause_verdict",
        (
            "DUAL GTF ARTIFACT (same class as HSV-1/anellovirus): "
            "(A) 80/94 EBV genes have CDS-only records (no exon) — "
            "structurally zero-counted by STARsolo GeneFull; "
            "(B) the 14 exon-bearing genes have 23 overlapping pairs — "
            "82.8% of exon bases ambiguous, signal concentrates entirely in LMP-1 "
            "(46,343 of 46,419 anchor UMI). "
            "The tools measure near-disjoint gene sets: STAR 99.8% LMP-1, "
            "VS 95% CDS-only genes STAR cannot count. "
            "The single overlap point (LMP-1: 325 vs 46,343 UMI) is unresolved — "
            "sensitive to STAR ambiguity handling and VS unique-layer multimapper exclusion. "
            "The apparent VS/STAR=2x is a reference-completeness artifact "
            "(95% of VS EBV UMI from CDS-only genes). "
            "No clean EBV aligner comparison possible without fixing the GTF."
        ),
        "analysis/reference_strategy_benchmark/outputs/numbers.json",
    )
    _reg(
        values,
        "ebv_starsolo_n_exon_bearing_genes",
        14,
        "analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py:ebv_gtf_artifact",
    )
    _reg(
        values,
        "ebv_starsolo_n_cds_only_genes",
        80,
        "analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py:ebv_gtf_artifact",
    )
    _reg(
        values,
        "ebv_starsolo_exon_ambiguous_fraction",
        0.828,
        "analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py:ebv_gtf_artifact",
    )
    _reg(
        values,
        "ebv_starsolo_lmp1_fraction_of_signal",
        0.998,
        "analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py:ebv_gtf_artifact",
    )
    _reg(
        values,
        "ebv_viralscan_cdsonly_fraction_of_signal",
        0.950,
        "analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py:ebv_gtf_artifact",
    )

    # HSV-1 verdict (text value — string type, registered for completeness)
    _reg(
        values,
        "hsv1_starsolo_root_cause_verdict",
        (
            "DUAL GTF ARTIFACT: (A) 61/79 HHV1 genes have CDS-only records (no exon) "
            "— structurally zero-counted by STARsolo GeneFull; "
            "(B) the 18 exon-bearing genes cluster in terminal repeats with 11 "
            "overlapping pairs — 27.4% of exon bases ambiguous, most reads discarded. "
            "Net: 30 UMI total from 18 genes (19 on anchor). "
            "NOT aligner sensitivity — same mechanism as anellovirus STARsolo=0 (F-005). "
            "Fair HSV-1 comparison requires fixing GTF exon records."
        ),
        "analysis/reference_strategy_benchmark/outputs/hsv1_root_cause.json",
    )
    _reg(
        values,
        "hsv1_starsolo_unique_anchor_combined",
        19,
        f"{PROV_TSV}:dataset=hsv1,ref_strat=combined",
    )
    _reg(
        values,
        "hsv1_starsolo_n_exon_bearing_genes",
        18,
        "analysis/reference_strategy_benchmark/outputs/feature_match_audit.tsv:aligner=starsolo,dataset=hsv1",
    )
    _reg(
        values,
        "hsv1_starsolo_n_cds_only_genes",
        61,
        "analysis/reference_strategy_benchmark/outputs/hsv1_root_cause.json:gtf_diagnosis",
    )
    _reg(
        values,
        "hsv1_starsolo_exon_ambiguous_fraction",
        0.274,
        "analysis/reference_strategy_benchmark/outputs/hsv1_root_cause.json:gtf_diagnosis",
    )

    fragment = {"namespace": NS, "values": sorted(values, key=lambda x: x["key"])}

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSON.open("w", encoding="utf-8") as fh:
        json.dump(fragment, fh, indent=2)

    print(f"Registered {len(values)} values to {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

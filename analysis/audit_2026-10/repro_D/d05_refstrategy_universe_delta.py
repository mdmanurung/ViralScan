"""parsed_benchmark_rows: (a) the fixed barcode universe silently shrinks when one combined row
is incomplete; (b) combined_vs_two_step delta is always '0' on combined rows and has the sign
two_step - combined on two_step rows; (c) incomplete rows report 0 UMI beside the status."""
import tempfile
from pathlib import Path
from unittest.mock import patch

from viralscan import reference_strategy as rs

run = Path(tempfile.mkdtemp())
rs.write_commands(run, rs.default_manifest())
(run / "reference_audit.tsv").write_text("key\n")


def metrics(counts, status="complete", anchor=None):
    return rs.ParsedMetrics(counts, set(anchor if anchor is not None else counts), status,
                            "" if status == "complete" else "boom")


# EBV dataset: STARsolo combined FAILS; ViralScan combined only has the 2 barcodes with viral reads;
# STARsolo two_step saw 3 cells with target UMI (a third barcode ViralScan never listed).
table = {
    "ebv__starsolo__combined": metrics({}, "incomplete"),
    "ebv__starsolo__two_step": metrics({"A": (5.0, 0.0), "B": (5.0, 0.0), "C": (5.0, 0.0)}),
    "ebv__viralscan__combined": metrics({"A": (4.0, 0.0), "B": (4.0, 0.0)}),
    "ebv__viralscan__two_step": metrics({"A": (4.0, 0.0), "B": (4.0, 0.0)}),
}
empty = metrics({}, "incomplete")


def fake(run_dir, row):
    return table.get(row["row_id"], empty)


with patch.object(rs, "_row_metrics", fake):
    rows = rs.parsed_benchmark_rows(run)
for r in rows:
    if r["dataset"] == "ebv":
        print(r["command_id"], "status=", r["status"], "fixed=", r["denominator_fixed_barcodes"],
              "target_umi_fixed=", r["target_umi_fixed"], "delta=", r["combined_vs_two_step_delta_target_umi"])

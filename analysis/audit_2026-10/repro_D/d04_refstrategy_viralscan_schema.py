"""reference_strategy.parse_viralscan_metrics vs the v3 per_cell_viral.tsv schema, and
the DATASETS target/off-target regexes vs the shipped virus catalogue names."""
import re
import tempfile
from pathlib import Path

from viralscan import virus_catalog as vc
from viralscan import virus_identity as vi
from viralscan.reference_strategy import DATASETS, parse_viralscan_metrics

# 1. v3 schema (detection.py:684-693): viral_molecules_total_est, no viral_umi; only infected cells.
tmp = Path(tempfile.mkdtemp())
res = tmp / "results"
res.mkdir()
(res / "per_cell_viral.tsv").write_text(
    "barcode\tvirus_name\tviral_molecules_total_est\tmolecules_total_est\tviral_fraction\tis_called_cell\n"
    "BC1\tEpstein-Barr virus\t3\t1000\t0.003\ttrue\n"
)
(res / "viral_summary.tsv").write_text("virus_name\ttotal_cells\nEpstein-Barr virus\t100\n")
m = parse_viralscan_metrics(tmp, target_regex=DATASETS[1]["target_regex"],
                            off_target_regex=DATASETS[1]["off_target_regex"])
print("v3 per_cell_viral.tsv ->", m.status, "|", m.failure_reason)

# 2. regex cross-matching over catalogue display names
cat = vi._Catalogue(vc.load_catalogue(), {})
names = sorted(set(cat.names.values()))
for d in DATASETS:
    tgt = [n for n in names if re.search(d["target_regex"], n)]
    off = [n for n in names if re.search(d["off_target_regex"], n)]
    print(f"{d['target_virus']}: target matches {len(tgt)}: {tgt[:8]}")
    print(f"{d['target_virus']}: off-target matches {len(off)}: {off[:12]}")
    both = [n for n in tgt if n in off]
    print(f"{d['target_virus']}: matched by BOTH: {both}")

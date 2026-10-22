#!/usr/bin/env python3
"""Generate ``src/viralscan/data/virus_taxonomy.tsv``: the NCBI Taxonomy lineage of every catalogue taxid.

PLAN `TAXO-01`. The catalogue's ``species``/``genus``/``family`` come from each GenBank flatfile, whose
lineage is frozen at submission: EBV's record still says Herpesviridae with no species binomial, while
EBV type 2's says Orthoherpesviridae and "Lymphocryptovirus humangamma4". NCBI Taxonomy puts both under
one species taxid (3050299), so this table, not the catalogue labels, is what counts are summed over.

    python extras/build_virus_taxonomy.py [--email you@example.org] [--check]

``--check`` reports taxids that are missing or whose lineage changed, and writes nothing.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from viralscan.scripts.ncbi_fetch import TAXONOMY_RANKS, fetch_taxonomy  # noqa: E402
from viralscan.virus_catalog import catalogue_path, load_catalogue  # noqa: E402

DEFAULT_OUT = REPO_ROOT / "src" / "viralscan" / "data" / "virus_taxonomy.tsv"
COLUMNS = (
    ["taxid", "name", "rank"]
    + [c for r in TAXONOMY_RANKS for c in (f"{r}_taxid", r)]
    + ["retrieved_utc"]
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--catalogue", default=None, help="virus_catalog.tsv (default: the packaged one)"
    )
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    taxids = sorted(
        {(r.get("taxid") or "").strip() for r in load_catalogue(catalogue_path(a.catalogue))}
        - {""},
        key=int,
    )
    lineage = fetch_taxonomy(taxids, a.email, os.environ.get("NCBI_API_KEY"))
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    rows = [{"taxid": t, **lineage[t], "retrieved_utc": now} for t in taxids if t in lineage]
    missing = [t for t in taxids if t not in lineage]
    if missing:
        print(
            f"{len(missing)} taxid(s) not in NCBI Taxonomy: {', '.join(missing)}", file=sys.stderr
        )
    if a.check:
        old = {}
        if Path(a.out).exists():
            with open(a.out, encoding="utf-8") as fh:
                old = {r["taxid"]: r for r in csv.DictReader(fh, delimiter="\t")}
        changed = [
            r["taxid"]
            for r in rows
            if {k: v for k, v in old.get(r["taxid"], {}).items() if k != "retrieved_utc"}
            != {k: v for k, v in r.items() if k != "retrieved_utc"}
        ]
        print(f"{len(changed)} of {len(rows)} taxid(s) new or changed: {', '.join(changed[:20])}")
        return 1 if changed or missing else 0
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {a.out}: {len(rows)} taxids ({len(missing)} missing)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

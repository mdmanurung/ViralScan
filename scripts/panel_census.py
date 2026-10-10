#!/usr/bin/env python3
"""PANEL-01 (WP1 census): RefSeq viral records whose host is Homo sapiens, against the catalogue (gap check).

    python scripts/panel_census.py [--refresh] [--email you@example.org] [--api-key KEY]
        [--cache analysis/panel_expansion/census_raw.json] [--out analysis/panel_expansion/census.tsv]

Asks NCBI nuccore for viral RefSeq records with a Homo sapiens /host, then marks which ones the catalogue
already holds (and in which `panel`). Records the catalogue lacks are what the user reviews; none is added
automatically. The raw NCBI answer is cached with its date, so a rerun without --refresh is offline and
byte-identical. Pure stdlib, no Biopython (same as ncbi_fetch.py).

Caveat the output cannot remove: NCBI's host qualifier is submitter-provided and many RefSeq viral records
carry none, so this is a lower bound on human-relevant viruses, not the list.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from viralscan.virus_catalog import load_catalogue  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TERM = (
    '"Viruses"[Organism] AND srcdb_refseq[PROP] AND biomol_genomic[PROP] AND "Homo sapiens"[Host]'
)
COLUMNS = ["accession", "title", "organism", "taxid", "length", "in_catalogue_panel"]


def _get(endpoint: str, params: dict[str, str], pause: float) -> dict:
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode({**params, "retmode": "json"})
    with urllib.request.urlopen(url, timeout=60) as response:  # noqa: S310 - fixed https host
        data = json.load(response)
    time.sleep(pause)
    return data


def fetch_raw(email: str | None, api_key: str | None) -> dict:
    """esearch + batched esummary; returns {'query', 'date', 'ids', 'summaries'}."""
    extra = {k: v for k, v in (("email", email), ("api_key", api_key)) if v}
    pause = 0.12 if api_key else 0.4  # NCBI: 10 requests/s with a key, 3 without
    search = _get("esearch.fcgi", {"db": "nuccore", "term": TERM, "retmax": "5000", **extra}, pause)
    ids = search["esearchresult"]["idlist"]
    summaries: dict[str, dict] = {}
    for start in range(0, len(ids), 200):
        batch = ids[start : start + 200]
        result = _get("esummary.fcgi", {"db": "nuccore", "id": ",".join(batch), **extra}, pause)[
            "result"
        ]
        summaries.update({uid: result[uid] for uid in batch if uid in result})
    return {
        "query": TERM,
        "date": datetime.date.today().isoformat(),
        "ids": ids,
        "summaries": summaries,
    }


def census_rows(raw: dict, catalogue: list[dict[str, str]]) -> list[dict[str, str]]:
    panel_of: dict[str, str] = {}
    for row in catalogue:
        panel_of[row["accession"].split(".")[0]] = row["panel"]
    rows = []
    for uid in raw["ids"]:
        s = raw["summaries"].get(uid)
        if not s:
            continue
        accession = s.get("accessionversion") or s.get("caption", "")
        rows.append(
            {
                "accession": accession,
                "title": s.get("title", ""),
                "organism": s.get("organism", ""),
                "taxid": str(s.get("taxid", "")),
                "length": str(s.get("slen", "")),
                "in_catalogue_panel": panel_of.get(accession.split(".")[0], ""),
            }
        )
    rows.sort(key=lambda r: r["accession"])
    return rows


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1] / "analysis" / "panel_expansion"
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cache", type=Path, default=root / "census_raw.json")
    ap.add_argument("--out", type=Path, default=root / "census.tsv")
    ap.add_argument(
        "--refresh", action="store_true", help="query NCBI again instead of using the cache"
    )
    ap.add_argument("--email")
    ap.add_argument("--api-key")
    args = ap.parse_args(argv)

    if args.cache.is_file() and not args.refresh:
        raw = json.loads(args.cache.read_text())
    else:
        raw = fetch_raw(args.email, args.api_key)
        args.cache.parent.mkdir(parents=True, exist_ok=True)
        args.cache.write_text(json.dumps(raw, indent=1, sort_keys=True) + "\n")
    rows = census_rows(raw, load_catalogue())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["in_catalogue_panel"] or "NOT_IN_CATALOGUE"] = (
            counts.get(r["in_catalogue_panel"] or "NOT_IN_CATALOGUE", 0) + 1
        )
    print(
        f"census {raw['date']}: {len(rows)} records: {dict(sorted(counts.items()))}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""CAT-02/CAT-09 human-host RefSeq sweep: Virus-Host DB (host 9606) -> catalogue candidates.

Selects every non-phage Virus-Host DB taxon with a human host that the catalogue
does not yet hold (by RefSeq accession or taxid), fetches each cache-first through
``ncbi_fetch._fetch_one``, re-checks the GenBank /host, collapses to one taxon per
species and writes ``analysis/cat09_sweep/``. Nothing is indexed; the kept
accessions are merged into the catalogue by ``build_virus_catalog.py`` (panel=broad).

Run:  PYTHONPATH=src:extras python extras/cat09_human_host_sweep.py \\
          --vhdb <virushostdb.tsv> --cache <ncbi cache dir> --email you@example.org
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import threading
import xml.etree.ElementTree as ET
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
from build_virus_catalog import catalog_row  # noqa: E402
import requests  # noqa: E402
from viralscan.scripts import ncbi_fetch  # noqa: E402
from viralscan.scripts.ncbi_fetch import NCBIFetchError, genbank_cache_path  # noqa: E402

CATALOGUE = REPO_ROOT / "src" / "viralscan" / "data" / "virus_catalog.tsv"
OUT = REPO_ROOT / "analysis" / "cat09_sweep"
HUMAN_TAXID = "9606"
PHAGE_RE = re.compile(
    r"Microviridae|Microviricetes|Gokushovir|Caudoviricetes|Inoviridae|Fiersviridae|Leviviridae|Leviviricetes|bacterial virus|phage", re.I
)
# same semantics as viral_panel_max_2026-09-28/01_collect.py
DECOY_RE = re.compile(
    r"murine leukemia|xenotropic|squirrel monkey retrovirus|bovine viral diarrh|"
    r"porcine circovirus|bovine polyomavirus|murine type c retrovirus",
    re.I,
)
HUMAN_HOST_RE = re.compile(r"homo sapiens|\bhuman\b", re.I)

_lock = threading.Lock()
_last = [0.0]


def _throttled(real, min_gap=0.34):
    """Global <=3 requests/s across all worker threads."""

    def wrapper(*a, **kw):
        with _lock:
            wait = _last[0] + min_gap - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            _last[0] = time.monotonic()
        return real(*a, **kw)

    return wrapper


def load_candidates(vhdb: Path, catalogue: Path):
    with open(catalogue, newline="") as fh:
        cat = list(csv.DictReader(fh, delimiter="\t"))
    cat_acc = {r["accession"] for r in cat}
    cat_tax = {r["taxid"] for r in cat if r["taxid"]}
    # species, organism and curated common_name: all can become a display name (virus_identity)
    cat_species = {r[c] for r in cat for c in ("species", "organism", "common_name") if r[c]}
    taxa: dict[str, dict] = {}
    n_human = 0
    with open(vhdb, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["host tax id"].strip() != HUMAN_TAXID:
                continue
            t = taxa.setdefault(r["virus tax id"], {"name": r["virus name"], "lineage": r["virus lineage"], "accs": set()})
            t["accs"] |= {a.strip().split(".")[0] for a in r["refseq id"].split(",") if a.strip()}
    n_human = len(taxa)
    nonphage = {k: v for k, v in taxa.items() if not PHAGE_RE.search(v["lineage"] + " " + v["name"])}
    todo = {
        k: v
        for k, v in nonphage.items()
        if v["accs"] and k not in cat_tax and not (v["accs"] & cat_acc)
    }
    return n_human, len(nonphage), todo, cat_species, cat_tax


def species_taxids(taxids: set[str], cache: Path, email: str) -> dict[str, str]:
    """taxid -> species-rank taxid ('' when the lineage has no species rank).

    NCBI Taxonomy efetch in batches of 200, <=3 req/s; results are cached in a TSV so reruns are offline.
    """
    known: dict[str, str] = {}
    if cache.exists():
        with open(cache, newline="") as fh:
            known = {r["taxid"]: r["species_taxid"] for r in csv.DictReader(fh, delimiter="\t")}
    todo = sorted(t for t in taxids if t and t not in known)
    for i in range(0, len(todo), 200):
        batch = todo[i : i + 200]
        for attempt in range(4):
            time.sleep(0.34)
            resp = requests.post(
                ncbi_fetch.EUTILS_BASE,
                data={"db": "taxonomy", "id": ",".join(batch), "retmode": "xml", "email": email, "tool": "ViralScan"},
                timeout=120,
            )
            if resp.status_code == 200:
                break
            time.sleep(2**attempt)
        else:
            continue  # left unresolved; reported by the caller, never silently dropped
        for taxon in ET.fromstring(resp.text).findall("Taxon"):
            tid = taxon.findtext("TaxId")
            sp = tid if taxon.findtext("Rank") == "species" else ""
            for anc in taxon.findall("LineageEx/Taxon"):
                if anc.findtext("Rank") == "species":
                    sp = anc.findtext("TaxId")
            known[tid] = sp or ""
            for alt in taxon.findall("AkaTaxIds/TaxId"):  # merged ids map to the survivor
                known.setdefault(alt.text, sp or "")
    cache.parent.mkdir(parents=True, exist_ok=True)
    with open(cache, "w", newline="") as fh:
        fh.write("taxid\tspecies_taxid\n")
        fh.writelines(f"{t}\t{known[t]}\n" for t in sorted(known, key=int))
    return known


def resolve(acc: str, cache: Path, email: str) -> str:
    hits = sorted(p.name for p in cache.glob(f"{acc}.*") if (p / f"{p.name}.gb").exists())
    if hits:
        return hits[-1]
    text = ncbi_fetch._efetch(acc, "acc", email, None).strip().splitlines()
    return text[0].strip() if text else acc


def fetch(acc: str, cache: Path, email: str):
    try:
        acc_v = resolve(acc, cache, email)
        ncbi_fetch._fetch_one(acc_v, cache, email, None)
        gb = genbank_cache_path(acc_v, cache)
        return acc, catalog_row(acc_v, gb.read_text(), gb), ""
    except (NCBIFetchError, OSError, ValueError) as exc:
        return acc, None, str(exc)[:300]


def host_decision(row: dict) -> tuple[str, str, str]:
    """(keep|review_host, role, reason); every candidate has VHDB human evidence."""
    host = row["host"]
    if DECOY_RE.search(row["species"]):
        return "keep", "decoy", "lab-contaminant decoy"
    if HUMAN_HOST_RE.search(host):
        return "keep", "target", f"host={host}"
    if not host:
        return "keep", "target", "no /host; Virus-Host DB lists Homo sapiens"
    return "review_host", "", f"host={host} organism={row['species']}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vhdb", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--email", required=True)
    ap.add_argument("--catalogue", type=Path, default=CATALOGUE)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args(argv)

    n_human, n_nonphage, todo, cat_species, cat_tax = load_candidates(a.vhdb, a.catalogue)
    acc_taxon = {acc: tx for tx, t in todo.items() for acc in t["accs"]}
    print(f"VHDB human taxa={n_human} non-phage={n_nonphage} not-catalogued={len(todo)} "
          f"accessions={len(acc_taxon)}", flush=True)

    ncbi_fetch._efetch = _throttled(ncbi_fetch._efetch)
    rows, errors = [], []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for i, (acc, row, err) in enumerate(
            pool.map(lambda x: fetch(x, a.cache, a.email), sorted(acc_taxon)), 1
        ):
            if row is None:
                errors.append((acc, acc_taxon[acc], err))
            else:
                row.update(vhdb_taxid=acc_taxon[acc], requested=acc)
                rows.append(row)
            if i % 100 == 0:
                print(f"  {i}/{len(acc_taxon)} {time.time() - t0:.0f}s errors={len(errors)}", flush=True)

    for r in rows:
        r["decision"], r["role"], r["reason"] = host_decision(r)
        # VHDB lineage can omit the phage family (e.g. Gokushovirus, Microviridae): re-check GenBank's.
        if PHAGE_RE.search(r["lineage"] + " " + r["species"]):
            r["decision"], r["role"], r["reason"] = "excluded_phage", "", "GenBank lineage is a phage taxon"
    # Segments of one taxon are all-or-none: a segmented taxon with a kept record keeps every
    # segment, even where one segment's /host names an isolate host (e.g. Snowshoe hare virus L).
    # Non-segmented taxa (several strains under one VHDB taxid) stay per record.
    by_taxon: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_taxon[r["vhdb_taxid"]].append(r)
    for rs in by_taxon.values():
        if any(r["segment"] for r in rs) and any(r["decision"] == "keep" for r in rs):
            for r in rs:
                if r["decision"] == "review_host":
                    r["decision"], r["role"] = "keep", "target"
                    r["reason"] = f"segment of a human-host taxon (segment /host={r['host']})"
    # One representative per true species (NCBI Taxonomy species-rank taxid). A VHDB taxon takes the
    # species of its lowest-accession record; a taxon with no species rank is its own group.
    keepers = [r for r in rows if r["decision"] == "keep"]
    sp_of = species_taxids({r["taxid"] for r in keepers} | cat_tax, a.out / "species_taxids.tsv", a.email)
    cat_sp = {sp_of.get(t, "") for t in cat_tax} - {""}
    by_vhdb: dict[str, list[dict]] = defaultdict(list)
    for r in keepers:
        by_vhdb[r["vhdb_taxid"]].append(r)
    unresolved = sorted({r["taxid"] for r in keepers if r["taxid"] not in sp_of})
    by_species: dict[str, dict[str, list[dict]]] = defaultdict(dict)
    for tx, rs in by_vhdb.items():
        first = min(rs, key=lambda r: r["accession_version"])
        sp = sp_of.get(first["taxid"], "")
        for r in rs:
            r["species_taxid"] = sp_of.get(r["taxid"], "")
        by_species[sp or f"(none:{tx})"][tx] = rs
    kept = []
    for sp, taxa in by_species.items():
        names = {r["organism"] for rs in taxa.values() for r in rs} | {r["species"] for rs in taxa.values() for r in rs}
        why = (
            f"species taxid {sp} already catalogued" if sp in cat_sp
            else "organism/species name already catalogued" if names & cat_species
            else ""
        )
        if why:
            for rs in taxa.values():
                for r in rs:
                    r["decision"], r["reason"] = "skip_species_in_catalogue", why
            continue
        rep = min(taxa, key=lambda tx: min(r["accession_version"] for r in taxa[tx]))
        for tx, rs in taxa.items():
            for r in rs:
                if tx == rep:
                    kept.append(r)
                else:
                    r["decision"], r["reason"] = "collapsed_species", f"species taxid {sp}: representative taxon {rep}"
    kept.sort(key=lambda r: r["accession_version"])

    fixed = ["requested", "vhdb_taxid", "species_taxid", "decision", "role", "reason"]
    cols = fixed + [k for k in (rows[0] if rows else {}) if k not in fixed]
    a.out.mkdir(parents=True, exist_ok=True)

    def write(name, rs):
        with open(a.out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rs)

    allrows = sorted(rows, key=lambda r: (r["decision"], r["accession_version"]))
    for acc, tx, err in errors:
        allrows.append({"requested": acc, "vhdb_taxid": tx, "decision": "fetch_error", "reason": err})
    write("candidates.tsv", allrows)
    write("review_host.tsv", [r for r in allrows if r["decision"] == "review_host"])
    with open(a.out / "fetch_errors.tsv", "w") as fh:
        fh.write("accession\tvhdb_taxid\terror\n")
        fh.writelines(f"{x}\t{y}\t{z}\n" for x, y, z in errors)
    (a.out / "kept_accessions.txt").write_text("".join(r["accession_version"] + "\n" for r in kept))
    print(f"species groups={len(by_species)} unresolved_taxids={len(unresolved)} {unresolved[:10]} "
          f"no_species_rank={sum(k.startswith('(none') for k in by_species)}", flush=True)
    n_dec = sum(r["role"] == "decoy" for r in kept)
    print(f"fetched={len(rows)} errors={len(errors)} host_keep={sum(r['decision'] != 'review_host' for r in rows)} "
          f"review={sum(r['decision'] == 'review_host' for r in rows)} kept_taxa={len({r['vhdb_taxid'] for r in kept})} "
          f"kept_accessions={len(kept)} decoy_kept={n_dec}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""PANEL-01 (WP1b): k-mer sharing between each panel candidate and the built panel (offline, read-only).

    python scripts/panel_kmer_sharing.py [--candidates analysis/panel_expansion/candidates.tsv]
        [--panel-fasta <build>/viral.fa] [--candidate-fasta <cache>/candidates.fasta] [--fetch]
        [--out analysis/panel_expansion/kmer_sharing.tsv] [-k 31]

Run this before curating the candidates (Biomni review finding 18): it says which candidates the index could
tell apart from what is already in the panel, so a candidate that is nearly all shared k-mers is judged before
anyone argues for its relevance. Needs numpy; run it in the `viralscan_test_full` env.

Genomes are the not-excluded rows of `candidates.tsv` plus the panel. The panel is `--panel-fasta`; without it the
catalogue's `shipped` rows are read from the NCBI cache (unmasked, so it differs from a built `viral.fa` wherever
dustmasker wrote N). Candidate sequences come from `--candidate-fasta`; `--fetch` downloads the missing ones from
NCBI (needs NCBI_EMAIL) and appends them to that file. A candidate with no sequence is listed in the log, not scored.

k-mers are canonical (min of forward and reverse complement), counted once per genome, windows with a non-ACGT base
dropped. Origin-wrapping k-mers of circular genomes are not added, as in the index build.

A second file lists, for every candidate with under `--partner-below` of its k-mers outside the panel, its top 3 sharing
partners (panel or candidate), so "same type already shipped" and "mislabelled accession" can be told apart by eye.
`length_ratio` (candidate length / partner length) and `frac_of_partner` (share of the partner's k-mers found in the candidate)
separate a twin from a fragment. `--baseline` also scores every panel genome against the rest of the panel
(`kmer_panel_baseline.tsv`): the sharing the shipped panel already accepts, as a yardstick for the candidate fractions.
Rows `panel_candidates.py` excluded as `kmer_twin_of:<id>` from an earlier run are scored again, so a rerun is stable.

Output columns (fractions are of the candidate's own distinct k-mers; sibling and other can both be non-zero):
  group / group_basis   the unit "sibling" refers to: catalogue `sibling_group`, else (with `--vmr`) the ICTV
                        genus, else `genus`, else the species name. `species` basis means no genus is known,
                        so a same-genus relative counts as `other`.
  frac_unique           in no other genome (panel or candidate)
  frac_not_in_panel     in no panel genome (what the index would gain)
  frac_sibling          in at least one other genome of the same group
  frac_other_group      in at least one genome of a different group
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from viralscan.virus_catalog import load_catalogue  # noqa: E402

COLUMNS = [
    "accession",
    "source",
    "relevance",
    "species",
    "group",
    "group_basis",
    "genome_length",
    "n_kmers",
    "frac_unique",
    "frac_not_in_panel",
    "frac_sibling",
    "frac_other_group",
]
DEFAULT_CANDIDATE_FASTA = (
    Path.home() / ".cache" / "viralscan" / "panel_candidates" / "candidates.fasta"
)
FETCH_BATCH = 50
KMER_TWIN_PREFIX = (
    "kmer_twin_of:"  # exclusion reason written by panel_candidates.py from this script's output
)
PARTNER_COLUMNS = [
    "candidate",
    "source",
    "rank",
    "partner",
    "partner_species",
    "partner_role",
    "shared_kmers",
    "frac_of_candidate",
    "frac_of_partner",
    "length_ratio",
    "same_group",
]
BASELINE_COLUMNS = [
    "accession",
    "species",
    "group",
    "group_basis",
    "n_kmers",
    "frac_unique",
    "frac_sibling",
    "frac_other_group",
]
_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i  # lower case


def base_accession(accession: str) -> str:
    return accession.strip().upper().split(".")[0]


def read_fasta(path: Path):
    """Yield (base accession, upper-case sequence bytes)."""
    name, chunks = None, []
    with path.open("rb") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name, chunks = base_accession(line[1:].split()[0].decode()), []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, b"".join(chunks)


def canonical_kmers(seq: bytes, k: int) -> np.ndarray:
    """Sorted distinct canonical k-mers (2 bits per base, k <= 31) of one sequence."""
    code = _CODE[np.frombuffer(seq, dtype=np.uint8)]
    n = len(code) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64)
    fwd = np.zeros(n, dtype=np.uint64)
    rev = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        window = code[j : j + n]
        fwd = (fwd << np.uint64(2)) | (window & np.uint8(3)).astype(np.uint64)
        rev |= (np.uint64(3) - (window & np.uint8(3)).astype(np.uint64)) << np.uint64(2 * j)
    bad = np.concatenate(([0], np.cumsum(code > 3, dtype=np.int64)))
    ok = (bad[k : k + n] - bad[:n]) == 0
    return np.unique(np.minimum(fwd, rev)[ok])


def sharing(
    kmers: list[np.ndarray],
    group: list[int],
    is_panel: list[bool],
    targets: list[int],
    partner_below: float = 0.0,
    top_n: int = 3,
) -> dict[int, dict]:
    """Per target genome: the four fractions described in the module docstring.

    A target whose `frac_not_in_panel` is under `partner_below` also gets `partners`: up to `top_n`
    (genome index, shared k-mers) pairs, most shared first, ties to the lower index (panel genomes first).
    """
    lens = np.array([len(a) for a in kmers], dtype=np.int64)
    all_k = np.concatenate(kmers)
    gid = np.repeat(np.arange(len(kmers), dtype=np.int32), lens)
    grp = np.asarray(group, dtype=np.int32)[gid]
    order = np.lexsort((grp, all_k))
    k_s, g_s, id_s = all_k[order], grp[order], gid[order]
    del all_k, gid, grp, order
    new_k = np.concatenate(([True], k_s[1:] != k_s[:-1]))
    new_run = new_k | np.concatenate(([True], g_s[1:] != g_s[:-1]))
    kid = np.cumsum(new_k, dtype=np.int32) - 1
    rid = np.cumsum(new_run, dtype=np.int32) - 1
    n_in_kmer = np.bincount(kid)
    n_groups = np.bincount(kid, weights=new_run)
    run_size = np.bincount(rid)
    panel_row = np.asarray(is_panel)[id_s]
    n_panel = np.bincount(kid, weights=panel_row)
    flags = {
        "frac_unique": n_in_kmer[kid] == 1,
        "frac_not_in_panel": (n_panel[kid] - panel_row) == 0,
        "frac_sibling": run_size[rid] > 1,
        "frac_other_group": n_groups[kid] > 1,
    }
    out: dict[int, dict] = {t: {} for t in targets}
    for name, flag in flags.items():
        per_genome = np.bincount(id_s, weights=flag, minlength=len(kmers))
        for t in targets:
            out[t][name] = per_genome[t] / lens[t] if lens[t] else float("nan")
    for t in targets:
        if not out[t]["frac_not_in_panel"] < partner_below:
            continue
        in_target = np.zeros(len(n_in_kmer), dtype=bool)
        in_target[kid[id_s == t]] = True
        shared = np.bincount(id_s[in_target[kid]], minlength=len(kmers))
        shared[t] = 0
        top = np.argsort(-shared, kind="stable")[:top_n]
        out[t]["partners"] = [(int(i), int(shared[i])) for i in top if shared[i]]
    return out


def fetch_missing(accessions: list[str], fasta: Path) -> None:
    from viralscan.scripts.ncbi_fetch import _efetch  # noqa: PLC0415

    email = os.environ.get("NCBI_EMAIL")
    if not email:
        sys.exit("--fetch needs NCBI_EMAIL (NCBI's E-utilities terms)")
    have = {name for name, _ in read_fasta(fasta)} if fasta.is_file() else set()
    todo = [a for a in accessions if base_accession(a) not in have]
    fasta.parent.mkdir(parents=True, exist_ok=True)
    with fasta.open("a") as out:
        for i in range(0, len(todo), FETCH_BATCH):
            batch = todo[i : i + FETCH_BATCH]
            out.write(_efetch(",".join(batch), "fasta", email, os.environ.get("NCBI_API_KEY")))
            out.flush()
            print(f"fetched {i + len(batch)}/{len(todo)}", file=sys.stderr)


_XLSX = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_ACCESSION = re.compile(r"[A-Z]{1,2}_?\d{5,}")
_NAME_TAIL = re.compile(r",? (?:segment|isolate|strain|complete|partial|genomic)\b.*$")


def clean_name(name: str) -> str:
    """Virus name without the isolate, segment and 'complete genome' tail, for matching to ICTV names."""
    return re.sub(r"[^a-z0-9]+", " ", _NAME_TAIL.sub("", name.lower())).strip()


def load_vmr(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """ICTV VMR xlsx -> ({base accession: genus}, {cleaned virus/species name: genus}).

    Reads the sheet with stdlib only (no openpyxl). Rows without a genus are skipped; a name that maps
    to two genera is dropped from the name index rather than guessed.
    """
    with zipfile.ZipFile(path) as z:
        strings = [
            "".join(t.text or "" for t in si.iter(_XLSX + "t"))
            for si in ElementTree.fromstring(z.read("xl/sharedStrings.xml")).findall(_XLSX + "si")
        ]
        sheet = ElementTree.fromstring(z.read("xl/worksheets/sheet2.xml"))  # "VMR MSLnn"
    by_acc: dict[str, str] = {}
    by_name: dict[str, str | None] = {}
    for row in sheet.iter(_XLSX + "row"):
        cells = {}
        for c in row.findall(_XLSX + "c"):
            v = c.find(_XLSX + "v")
            if v is not None:
                col = "".join(ch for ch in c.get("r") if ch.isalpha())
                cells[col] = strings[int(v.text)] if c.get("t") == "s" else v.text
        genus = cells.get("P", "")  # columns: P Genus, R Species, U Virus name(s), X accession
        if not genus or genus == "Genus":
            continue
        for acc in _ACCESSION.findall(cells.get("X", "")):
            by_acc[acc] = genus
        for name in [cells.get("R", ""), *cells.get("U", "").split(";")]:
            key = clean_name(name)
            if key and by_name.setdefault(key, genus) != genus:
                by_name[key] = None
    return by_acc, {k: v for k, v in by_name.items() if v}


def group_of(
    info: dict, vmr: tuple[dict[str, str], dict[str, str]] | None = None, accession: str = ""
) -> tuple[str, str]:
    """The unit "sibling" refers to: catalogue `sibling_group`, else ICTV genus (by accession, then name),
    else catalogue/candidate `genus`, else the species name. ICTV comes before the catalogue genus so that
    panel and candidates are named the same way wherever the VMR knows them."""
    if info.get("sibling_group"):
        return info["sibling_group"], "sibling_group"
    if vmr:
        by_acc, by_name = vmr
        genus = by_acc.get(accession) or by_name.get(clean_name(info.get("species", "")))
        if genus:
            return genus, "ictv_genus"
    if info.get("genus"):
        return info["genus"], "genus"
    return info.get("species") or "unknown", "species"


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1] / "analysis" / "panel_expansion"
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--candidates", type=Path, default=root / "candidates.tsv")
    ap.add_argument("--panel-fasta", type=Path)
    ap.add_argument("--candidate-fasta", type=Path, default=DEFAULT_CANDIDATE_FASTA)
    ap.add_argument(
        "--fetch", action="store_true", help="download missing candidate FASTAs from NCBI"
    )
    ap.add_argument("--out", type=Path, default=root / "kmer_sharing.tsv")
    ap.add_argument("--partners-out", type=Path, default=root / "kmer_partners.tsv")
    ap.add_argument("--partner-below", type=float, default=0.5)
    ap.add_argument(
        "--baseline", action="store_true", help="also score the panel genomes against each other"
    )
    ap.add_argument("--baseline-out", type=Path, default=root / "kmer_panel_baseline.tsv")
    ap.add_argument(
        "--vmr",
        type=Path,
        help="ICTV VMR xlsx (https://ictv.global/vmr): group by ICTV genus before the catalogue genus",
    )
    ap.add_argument("-k", type=int, default=31)
    args = ap.parse_args(argv)

    with args.candidates.open() as handle:
        cands = [
            r
            for r in csv.DictReader(handle, delimiter="\t")
            if not r["exclusion_reason"] or r["exclusion_reason"].startswith(KMER_TWIN_PREFIX)
        ]
    if args.fetch:
        fetch_missing([r["accession"] for r in cands], args.candidate_fasta)

    vmr = load_vmr(args.vmr) if args.vmr else None
    catalogue = load_catalogue()
    info = {base_accession(r["accession"]): r for r in catalogue}
    cand_acc = {base_accession(r["accession"]) for r in cands}

    panel_seqs: list[tuple[str, bytes]] = []
    if args.panel_fasta:
        panel_seqs = list(read_fasta(args.panel_fasta))
    else:
        from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR  # noqa: PLC0415

        for r in catalogue:
            if r["panel"] == "shipped":
                av = r["accession_version"]
                path = DEFAULT_CACHE_DIR / av / f"{av}.fasta"
                if path.is_file():
                    panel_seqs.extend(read_fasta(path))
        print(
            "panel from the NCBI cache (unmasked); pass --panel-fasta for the built viral.fa",
            file=sys.stderr,
        )
    # a panel record that is itself a candidate would match itself
    panel_seqs = [(a, s) for a, s in panel_seqs if a not in cand_acc]

    cand_seqs = dict(read_fasta(args.candidate_fasta)) if args.candidate_fasta.is_file() else {}
    scored = [r for r in cands if base_accession(r["accession"]) in cand_seqs]
    missing = len(cands) - len(scored)
    if missing:
        print(
            f"{missing} candidates have no sequence in {args.candidate_fasta} (use --fetch)",
            file=sys.stderr,
        )

    groups: dict[str, int] = {}
    kmers, group_ids, is_panel, meta, names, lengths = [], [], [], [], [], []
    panel_basis: list[str] = []
    for acc, seq in panel_seqs:
        g, basis = group_of(info.get(acc, {"species": f"unknown:{acc}"}), vmr, acc)
        panel_basis.append(basis)
        lengths.append(len(seq))
        names.append((acc, info.get(acc, {}).get("species", ""), "panel"))
        kmers.append(canonical_kmers(seq, args.k))
        group_ids.append(groups.setdefault(g, len(groups)))
        is_panel.append(True)
        meta.append(None)
    for r in scored:
        base = base_accession(r["accession"])
        g, basis = group_of(
            {
                **info.get(base, {}),
                "genus": r["genus"] or info.get(base, {}).get("genus", ""),
                "species": r["species"],
            },
            vmr,
            base,
        )
        seq = cand_seqs[base]
        names.append((base, r["species"], "candidate"))
        lengths.append(len(seq))
        kmers.append(canonical_kmers(seq, args.k))
        group_ids.append(groups.setdefault(g, len(groups)))
        is_panel.append(False)
        meta.append({**r, "group": g, "group_basis": basis, "genome_length": len(seq)})
    targets = [i for i, m in enumerate(meta) if m]
    stats = sharing(kmers, group_ids, is_panel, targets, args.partner_below)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for t in sorted(targets, key=lambda i: (meta[i]["source"], meta[i]["accession"])):
            row = {c: meta[t].get(c, "") for c in COLUMNS}
            row["n_kmers"] = len(kmers[t])
            row.update({k: f"{v:.4f}" for k, v in stats[t].items() if k != "partners"})
            writer.writerow(row)
    with args.partners_out.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=PARTNER_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for t in sorted(targets, key=lambda i: (meta[i]["source"], meta[i]["accession"])):
            for rank, (p, shared) in enumerate(stats[t].get("partners", []), 1):
                writer.writerow(
                    {
                        "candidate": meta[t]["accession"],
                        "source": meta[t]["source"],
                        "rank": rank,
                        "partner": names[p][0],
                        "partner_species": names[p][1],
                        "partner_role": names[p][2],
                        "shared_kmers": shared,
                        "frac_of_candidate": f"{shared / len(kmers[t]):.4f}",
                        "frac_of_partner": f"{shared / len(kmers[p]):.4f}",
                        "length_ratio": f"{lengths[t] / lengths[p]:.3f}",
                        "same_group": "yes" if group_ids[p] == group_ids[t] else "no",
                    }
                )
    if args.baseline:
        n_panel = len(panel_seqs)
        base = sharing(
            kmers[:n_panel], group_ids[:n_panel], is_panel[:n_panel], list(range(n_panel))
        )
        with args.baseline_out.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=BASELINE_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            group_name = {gid: g for g, gid in groups.items()}
            for i in sorted(range(n_panel), key=lambda i: names[i][0]):
                writer.writerow(
                    {
                        "accession": names[i][0],
                        "species": names[i][1],
                        "group": group_name[group_ids[i]],
                        "group_basis": panel_basis[i],
                        "n_kmers": len(kmers[i]),
                        **{
                            k: f"{base[i][k]:.4f}"
                            for k in ("frac_unique", "frac_sibling", "frac_other_group")
                        },
                    }
                )
        for col in ("frac_unique", "frac_other_group"):
            q = np.nanpercentile([base[i][col] for i in range(n_panel)], [5, 50, 95])
            print(
                f"panel baseline {col} 5/50/95th percentile: {q.round(3).tolist()}", file=sys.stderr
            )
    frac = np.array([stats[t]["frac_not_in_panel"] for t in targets])
    print(
        f"{len(targets)} candidates scored against {len(panel_seqs)} panel genomes (k={args.k}); "
        f"{int((frac < 0.5).sum())} have under half their k-mers outside the panel",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

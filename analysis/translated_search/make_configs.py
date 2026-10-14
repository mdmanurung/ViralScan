#!/usr/bin/env python3
"""Write one arm-P (translated search) config per chemistry group from a DSR manifest.

Inputs are a manifest TSV (dataset, sample, grp, r1, r2, whitelist, ...) and one
chem/<sample>.json per sample (technology under ``chemistry.chemistry``, strand under
``strand_pilot.choice``). The Snakefile takes one technology and at most one whitelist
per config, so samples are grouped by (technology, strand, whitelist). ``strand`` is
recorded for documentation only; the workflow does not pass it to ``kb count``.

ERR13027027's barcode read is its ``_2`` file (28 bp; ``_1`` is 89 bp). The manifest already lists
``_2`` as r1, so no correction is needed; ``--swap`` exists for manifests that get it wrong.
GSM4682311 is 10xv2 in its chem JSON (user decision 2026-10-06) although the manifest hint
says 10xv3; the chem JSON is authoritative, the manifest hint is never read.
Each R1/R2 pair is checked: both files exist and R1's modal read length equals the length
recorded in the chem JSON (``chemistry.r1_length``). Samples failing a check are excluded and the exit status is 1.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

import yaml

TECHNOLOGIES = {"10xv2", "10xv3", "dropseq"}
DEFAULT_SWAP: tuple[str, ...] = ()


def read_lengths(path: Path, n: int) -> Counter:
    opener = gzip.open if str(path).endswith(".gz") else open
    lengths: Counter = Counter()
    with opener(path, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 == 1:
                lengths[len(line.rstrip("\n"))] += 1
                if sum(lengths.values()) >= n:
                    break
    return lengths


def modal(counter: Counter) -> int | None:
    return counter.most_common(1)[0][0] if counter else None


def load_taxa(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with open(path, newline="") as handle:
        return {r["dataset"]: r["taxon_ids"] for r in csv.DictReader(handle, delimiter="\t")}


def expected_taxa(dataset: str, taxa: dict[str, str]) -> str:
    # exact dataset name first, then its first token (kshv_gse190558 -> kshv)
    return taxa.get(dataset) or taxa.get(dataset.split("_")[0], "")


def describe(row: dict, chem: dict, swap: set[str], peek: int) -> dict:
    """Validate one manifest row; returns the sample record or raises ValueError."""
    sample = row["sample"]
    technology = chem["chemistry"]["chemistry"]
    if technology not in TECHNOLOGIES:
        raise ValueError(f"{sample}: technology {technology!r} not one of {sorted(TECHNOLOGIES)}")
    r1, r2 = Path(row["r1"]), Path(row["r2"])
    corrections = []
    if sample in swap:
        r1, r2 = r2, r1
        corrections.append("r1/r2 swapped relative to manifest")
    for p in (r1, r2):
        if not p.is_file():
            raise ValueError(f"{sample}: missing file {p}")
    expected = chem["chemistry"]["r1_length"]
    got1, got2 = modal(read_lengths(r1, peek)), modal(read_lengths(r2, peek))
    if got1 != expected:
        raise ValueError(
            f"{sample}: R1 modal length {got1} != chem r1_length {expected} (R2 {got2})"
        )
    whitelist = row.get("whitelist", "-")
    whitelist = None if whitelist in ("", "-") else whitelist
    if whitelist is not None and not Path(whitelist).is_file():
        raise ValueError(f"{sample}: missing whitelist {whitelist}")
    if chem["chemistry"].get("basis") == "user on-list" and whitelist is None:
        raise ValueError(f"{sample}: chem JSON says user on-list but manifest gives no whitelist")
    strand = (chem.get("strand_pilot") or {}).get("choice", "unknown")
    return {
        "sample": sample,
        "dataset": row["dataset"],
        "technology": technology,
        "strand": strand,
        "whitelist": whitelist,
        "r1": str(r1),
        "r2": str(r2),
        "r1_length": got1,
        "r2_length": got2,
        "corrections": corrections,
    }


def group_name(tech: str, strand: str, whitelist: str | None) -> str:
    name = f"{tech}_{strand}"
    if whitelist:
        name += "_wl-" + re.sub(r"[^A-Za-z0-9.-]+", "-", Path(whitelist).stem)
    return name


def make_configs(
    manifest: Path,
    chem_dir: Path,
    base: Path,
    out_dir: Path,
    taxa_path: Path | None = None,
    swap: tuple[str, ...] = DEFAULT_SWAP,
    peek: int = 1000,
) -> dict:
    base_cfg = yaml.safe_load(base.read_text())
    taxa = load_taxa(taxa_path)
    records, excluded = [], {}
    with open(manifest, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                chem = json.loads((chem_dir / f"{row['sample']}.json").read_text())
                records.append(describe(row, chem, set(swap), peek))
            except (ValueError, OSError, KeyError) as exc:
                excluded[row["sample"]] = str(exc)
    groups: dict[str, list[dict]] = {}
    for rec in records:
        groups.setdefault(
            group_name(rec["technology"], rec["strand"], rec["whitelist"]), []
        ).append(rec)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, recs in groups.items():
        cfg = dict(
            base_cfg,
            technology=recs[0]["technology"],
            strand=recs[0]["strand"],
            samples={r["sample"]: [r["r1"], r["r2"]] for r in recs},
        )
        cfg.pop("whitelist", None)
        if recs[0]["whitelist"]:
            cfg["whitelist"] = recs[0]["whitelist"]
        (out_dir / f"{name}.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    with open(out_dir / "expected_virus.tsv", "w") as handle:
        handle.write("sample\ttaxon_ids\n")
        for rec in records:
            handle.write(f"{rec['sample']}\t{expected_taxa(rec['dataset'], taxa)}\n")
    summary = {
        "groups": {n: [r["sample"] for r in recs] for n, recs in groups.items()},
        "samples": {r["sample"]: r for r in records},
        "excluded": excluded,
    }
    (out_dir / "groups.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--chem-dir", type=Path, required=True)
    parser.add_argument(
        "--base",
        type=Path,
        required=True,
        help="YAML with the shared keys (kb, kallisto, work_dir, ...)",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--taxa", type=Path, help="TSV dataset<TAB>taxon_ids (';'-separated) for expected_virus.tsv"
    )
    parser.add_argument(
        "--swap",
        nargs="*",
        default=list(DEFAULT_SWAP),
        help="samples whose manifest R1/R2 are swapped",
    )
    parser.add_argument(
        "--peek", type=int, default=1000, help="reads peeked per file for length checks"
    )
    args = parser.parse_args()
    summary = make_configs(
        args.manifest,
        args.chem_dir,
        args.base,
        args.out_dir,
        args.taxa,
        tuple(args.swap),
        args.peek,
    )
    for name, samples in summary["groups"].items():
        print(f"{name}: {len(samples)} sample(s)")
    for why in summary["excluded"].values():
        print(f"EXCLUDED {why}", file=sys.stderr)
    return 1 if summary["excluded"] else 0


if __name__ == "__main__":
    raise SystemExit(main())

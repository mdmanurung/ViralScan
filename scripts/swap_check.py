#!/usr/bin/env python3
"""M-03: cross-run barcode-overlap check for index or barcode swapping (tier 1, descriptive).

    python scripts/swap_check.py pool_manifest.tsv --out swap_check.tsv

ViralScan cannot tell which libraries shared a sequencing pool and a barcode whitelist, so the user
declares it in a pool manifest (tab-separated, one row per library, header required):

    library_id   run_dir                  pool_id   whitelist
    donor_lib    runs/donor_lib/out       pool1     10x-3M-february-2018
    recip_lib    runs/recip_lib/out       pool1     10x-3M-february-2018

Two libraries are compared only when `pool_id` AND `whitelist` match and `library_id` differs. A
library is one row: list lanes or technical replicates of one library once (their shared barcodes are
legitimate and would read as 100 % swapping). `run_dir` is a finished ViralScan run directory,
relative paths resolving against the manifest's folder. Not for Drop-seq or any chemistry whose
barcodes are random per library, since it has no shared whitelist.

For each ordered pair (recipient, donor) and each virus with a non-called viral barcode in the
recipient, the test asks whether those barcodes are the donor's virus barcodes more often than chance.
A swapped molecule keeps its donor barcode, so it mostly lands in barcodes the recipient did not call.

    S  = recipient barcodes carrying the virus that are not called cells (per_cell_viral.tsv)
    D  = donor barcodes carrying the virus, called or not (per_cell_viral.tsv)
    U  = every barcode the recipient observed that is not a called cell
         (kb-python/counts_unfiltered/cells_x_genes.barcodes.txt minus results/called_cells.tsv)
    overlap_with_donor = |S & D|, overlap_expected = |S| * |U & D| / |U|,
    overlap_p = upper hypergeometric tail P(X >= |S & D|).

The result is descriptive. Overlap above chance is a signature of a shared origin (index hopping, a
chimeric read, carry-over, a technical replicate declared as two libraries), not a mechanism, and
no overlap with the donors listed says nothing about libraries not listed. It never changes `detected`.
`per_cell_viral.tsv` only lists barcodes at or above `--infected-min-molecules`, so weaker swapping is
invisible. Opt-in --tier2 compares retained corrected BUS keys for the same target;
barcode-label permutations provide a descriptive collision null, never a swapping verdict.
"""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import permutations
from pathlib import Path
from typing import Any

from scipy.stats import hypergeom

from viralscan.scripts.cellcalling import _strip_suffix

#: Conventional one-sided level for the descriptive label. Nothing gates on it.
OVERLAP_ALPHA = 0.05

COLUMNS = [
    "pool_id",
    "whitelist",
    "recipient_library",
    "donor_library",
    "virus_name",
    "n_noncalled_virus_barcodes",
    "n_donor_virus_barcodes",
    "n_noncalled_observed",
    "n_noncalled_in_donor",
    "overlap_with_donor",
    "overlap_expected",
    "overlap_p",
    "swap_assessment",
]
MANIFEST_COLUMNS = ("library_id", "run_dir", "pool_id", "whitelist")


def read_manifest(path: Path) -> list[dict[str, str]]:
    """Validated manifest rows with `run_dir` resolved against the manifest's folder."""
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)
    missing = [c for c in MANIFEST_COLUMNS if c not in header]
    if missing:
        raise ValueError(f"{path}: missing column(s) {missing}; need {list(MANIFEST_COLUMNS)}")
    ids = [r["library_id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{path}: duplicate library_id; list lanes of one library once")
    dirs = [str(Path(path.parent, r["run_dir"]).resolve()) for r in rows]
    if len(set(dirs)) != len(dirs):
        raise ValueError(f"{path}: two libraries share one run_dir")
    for row, run_dir in zip(rows, dirs):
        if not all(row[c].strip() for c in MANIFEST_COLUMNS):
            raise ValueError(f"{path}: library {row['library_id']!r} has an empty field")
        row["run_dir"] = run_dir
    return rows


def _read_barcodes(path: Path, skip_header: bool = False) -> set[str]:
    with path.open(encoding="utf-8") as fh:
        if skip_header:
            fh.readline()
        return {_strip_suffix(line.strip()) for line in fh if line.strip()}


def load_run(run_dir: Path) -> dict[str, Any]:
    """Per-virus viral barcodes (`called`, `noncalled`) and the non-called universe of one run."""
    per_cell = run_dir / "results" / "per_cell_viral.tsv"
    called_file = run_dir / "results" / "called_cells.tsv"
    barcodes_file = run_dir / "kb-python" / "counts_unfiltered" / "cells_x_genes.barcodes.txt"
    for needed in (per_cell, called_file, barcodes_file):
        if not needed.is_file():
            raise FileNotFoundError(
                f"{needed} is missing: {run_dir} is not a finished ViralScan run"
            )
    called = _read_barcodes(called_file, skip_header=True)
    universe = _read_barcodes(barcodes_file) - called
    viral: dict[str, dict[str, set[str]]] = {}
    with per_cell.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            bucket = "called" if row["is_called_cell"] == "True" else "noncalled"
            viral.setdefault(row["virus_name"], {"called": set(), "noncalled": set()})[bucket].add(
                _strip_suffix(row["barcode"])
            )
    for virus, sets in viral.items():
        if not sets["noncalled"] <= universe:
            raise ValueError(
                f"{run_dir}: {virus} has non-called viral barcodes absent from the observed "
                "barcodes list, so the barcode files belong to another run"
            )
    return {"universe": universe, "viral": viral}


def overlap_row(
    n_s: int, n_donor: int, universe_size: int, n_u_in_d: int, overlap: int
) -> dict[str, Any]:
    """The hypergeometric columns for one (recipient, donor, virus)."""
    expected = n_s * n_u_in_d / universe_size if universe_size else 0.0
    p = float(hypergeom.sf(overlap - 1, universe_size, n_u_in_d, n_s)) if universe_size else 1.0
    if n_u_in_d == 0:
        assessment = "not_assessed"  # no donor barcode is observed unclaimed here: no power
    else:
        assessment = "overlap_above_chance" if p < OVERLAP_ALPHA else "not_above_chance"
    return {
        "n_noncalled_virus_barcodes": n_s,
        "n_donor_virus_barcodes": n_donor,
        "n_noncalled_observed": universe_size,
        "n_noncalled_in_donor": n_u_in_d,
        "overlap_with_donor": overlap,
        "overlap_expected": float(f"{expected:.6g}"),
        "overlap_p": float(f"{p:.6g}"),
        "swap_assessment": assessment,
    }


def swap_rows(manifest: list[dict[str, str]]) -> list[dict[str, Any]]:
    runs = {r["library_id"]: load_run(Path(r["run_dir"])) for r in manifest}
    rows: list[dict[str, Any]] = []
    for recipient, donor in permutations(manifest, 2):
        if (recipient["pool_id"], recipient["whitelist"]) != (donor["pool_id"], donor["whitelist"]):
            continue
        rec, don = runs[recipient["library_id"]], runs[donor["library_id"]]
        for virus in sorted(rec["viral"]):
            s = rec["viral"][virus]["noncalled"]
            if not s:
                continue
            d = don["viral"].get(virus, {"called": set(), "noncalled": set()})
            d_all = d["called"] | d["noncalled"]
            rows.append(
                {
                    "pool_id": recipient["pool_id"],
                    "whitelist": recipient["whitelist"],
                    "recipient_library": recipient["library_id"],
                    "donor_library": donor["library_id"],
                    "virus_name": virus,
                    **overlap_row(
                        len(s),
                        len(d_all),
                        len(rec["universe"]),
                        len(rec["universe"] & d_all),
                        len(s & d_all),
                    ),
                }
            )
    return rows


TIER2_COLUMNS = [
    "pool_id",
    "whitelist",
    "recipient_library",
    "donor_library",
    "virus_name",
    "observed_shared_keys",
    "expected_shared_keys",
    "p_value",
    "insert_confirmation",
    "status",
    "reason_codes",
]


def load_molecule_keys(run_dir):
    """Corrected BUS keys and exact taxon compatibility; absent inputs are explicit."""
    import subprocess

    from viralscan.kb_outputs import KbCountOutputs
    from viralscan.multimapping import _resolved_gene_tuple
    from viralscan.run_safety import sha256_file
    from viralscan.runconfig import RunConfig
    from viralscan.scripts.multimap import load_transcripts, read_ec
    from viralscan.virus_grouping import load_run_identity

    kb = KbCountOutputs.from_config_output(run_dir)
    bus = next((p for p in (kb.resolved_bus, kb.corrected_bus) if p.is_file()), None)
    if bus is None:
        return {"status": "not_assessed", "reason_codes": "bus_missing"}
    identity = load_run_identity(run_dir)
    if identity is None:
        return {"status": "not_assessed", "reason_codes": "identity_missing"}
    config = RunConfig.from_yaml(Path(run_dir) / "config.yaml")
    transcripts, t2g = load_transcripts(str(kb.transcripts_txt), config.transcripts)
    genes = kb.genes.read_text().splitlines()
    ec_map = read_ec(str(kb.ec), transcripts, t2g, genes)
    taxa = {g.gene_id: g.virus_name for g in identity.genes if g.viral}
    ecs = {}
    with subprocess.Popen(
        ["bustools", "text", "-p", str(bus)], stdout=subprocess.PIPE, text=True
    ) as process:
        try:
            for line in process.stdout:
                fields = line.split()
                if len(fields) < 4:
                    raise ValueError("Malformed bustools text record")
                cb, umi, ec = fields[:3]
                ecs.setdefault((cb, umi), set()).add(int(ec))
        except BaseException:
            process.kill()
            raise
        if process.wait():
            raise RuntimeError(f"bustools text failed for {bus}")
    targets = {
        key: {taxa[genes[g]] for g in _resolved_gene_tuple(ids, ec_map) if genes[g] in taxa}
        for key, ids in ecs.items()
    }
    return {
        "status": "ok",
        "keys": targets,
        "barcodes": sorted({cb for cb, umi in targets}),
        "bus": str(bus),
        "bus_sha256": sha256_file(bus),
    }


def molecule_overlap(recipient, donor, virus, seed=42, n_permutations=1000):
    """Permutation of recipient barcode labels, preserving within-library loads."""
    import random

    source = sorted(key for key, targets in recipient["keys"].items() if virus in targets)
    dest = {key for key, targets in donor["keys"].items() if virus in targets}
    observed = len(set(source) & dest)
    barcodes = recipient["barcodes"]
    rng = random.Random(seed)
    null = []
    for _ in range(n_permutations):
        shuffled = list(barcodes)
        rng.shuffle(shuffled)
        mapping = dict(zip(barcodes, shuffled))
        null.append(sum((mapping[cb], umi) in dest for cb, umi in source))
    if not source or not dest or min(null) == max(null):
        return {
            "observed_shared_keys": observed,
            "expected_shared_keys": "not_assessed",
            "p_value": "not_assessed",
            "status": "not_assessed",
            "reason_codes": "collision_null_degenerate",
        }
    return {
        "observed_shared_keys": observed,
        "expected_shared_keys": sum(null) / len(null),
        "p_value": (1 + sum(value >= observed for value in null)) / (len(null) + 1),
        "status": "reviewed",
        "reason_codes": "insert_confirmation_not_assessed",
    }


def tier2_rows(manifest, seed=42):
    admitted = [
        (rec, don)
        for rec, don in permutations(manifest, 2)
        if (rec["pool_id"], rec["whitelist"]) == (don["pool_id"], don["whitelist"])
    ]
    ids = {row["library_id"] for pair in admitted for row in pair}
    runs = {
        row["library_id"]: load_molecule_keys(Path(row["run_dir"]))
        for row in manifest
        if row["library_id"] in ids
    }
    rows = []
    for recipient, donor in admitted:
        rec, don = runs[recipient["library_id"]], runs[donor["library_id"]]
        viruses = sorted({v for targets in rec.get("keys", {}).values() for v in targets}) or [""]
        for virus in viruses:
            row = dict.fromkeys(TIER2_COLUMNS, "not_assessed")
            row.update(
                pool_id=recipient["pool_id"],
                whitelist=recipient["whitelist"],
                recipient_library=recipient["library_id"],
                donor_library=donor["library_id"],
                virus_name=virus,
                insert_confirmation="not_assessed",
            )
            if rec["status"] != "ok" or don["status"] != "ok":
                row.update(
                    status="not_assessed",
                    reason_codes=";".join(
                        sorted({r["reason_codes"] for r in (rec, don) if r["status"] != "ok"})
                    ),
                )
            else:
                row.update(molecule_overlap(rec, don, virus, seed))
            rows.append(row)
    receipt = {
        "seed": seed,
        "permutations": 1000,
        "null": "recipient_barcode_label_permutation_v1",
        "p_value": "(1 + null overlaps >= observed) / (1000 + 1)",
        "libraries": {
            name: {k: v for k, v in run.items() if k not in {"keys", "barcodes"}}
            for name, run in runs.items()
        },
    }
    return rows, receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("manifest", type=Path, help="pool manifest TSV (see module docstring)")
    parser.add_argument("--out", type=Path, required=True, help="output TSV")
    parser.add_argument(
        "--tier2", action="store_true", help="Corrected CB–UMI target fingerprints (diagnostic)."
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Tier-2 permutation seed (default 42)."
    )
    args = parser.parse_args(argv)
    try:
        manifest = read_manifest(args.manifest)
        receipt = None
        if args.tier2:
            rows, receipt = tier2_rows(manifest, args.seed)
        else:
            rows = swap_rows(manifest)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"swap_check: {exc}", file=sys.stderr)
        return 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=TIER2_COLUMNS if args.tier2 else COLUMNS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    if receipt is not None:
        import json

        Path(str(args.out) + ".json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"swap_check: {len(rows)} (recipient, donor, virus) row(s) -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

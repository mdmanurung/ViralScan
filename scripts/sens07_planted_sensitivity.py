#!/usr/bin/env python3
"""Per-dataset planted-read sensitivity (PLAN SENS-07; exploratory, not VAL-01).

Plants exact (panel) and held-out (absent from the panel) genomes into one dataset's own
background library, in that library's own geometry, then scores recovery per family.

    plant   background FASTQ pair + FASTA(s) -> planted/ and unplanted/ FASTQ pairs, truth.tsv,
            selection.tsv, plant_manifest.json, and (with --control-genes) the
            --positive-control-* arguments plus a checksum-bound receipt per control gene.
    score   truth.tsv + the count matrix of a ViralScan run on planted/ (and, optionally, on
            unplanted/) -> recovery by condition, per-family capture and LOD.

It reuses ``plant_anello_10x.py`` (read simulation, UMIs, FASTQ writing, background streaming)
and does not write a second planter. What is planted:

* **Geometry** comes from ``--technology`` (barcode and UMI length) and the background's own R2
  length; R1 length is checked against barcode + UMI. An R1 shorter than that is refused. A longer R1
  is refused unless ``--trim-r1-to-geometry`` is given: it then trims R1 to barcode + UMI in BOTH the
  planted and unplanted twin (identical trim), a declared read change recorded in
  ``plant_manifest.json`` (``r1_trim``, ``read_changes``) with the original and trimmed lengths.
* **Cells**: the most frequent barcodes by read count (an outcome-independent anchor, as in the
  truth-panel D1 mode). Each (genome, level) condition gets its own block of ``--cells-per-condition``
  cells, so conditions are told apart by barcode in a single run.
* **Collisions**: a planted (barcode, UMI) never reuses one already in the background.
* **Selection** (``--selection-seed-mode``): ``global`` (default, round 2b) draws every family from
  one RNG in sorted order, so adding a family can change another family's picks; ``per-family``
  (round 2c, Q64) seeds each family's draw from sha256(seed, family), so a family's genomes are the
  same in every dataset whatever other families it plants.
* **Exact** = genomes of ``--exact-fasta`` (the panel's own sequences). **Held-out** = genomes of
  ``--heldout-fasta``: the caller supplies them and their absence from the panel (only the IDs are
  checked here; the partition rule is VAL-RA-DECISION-01's, not this script's). ``--paf`` records
  each held-out genome's nearest-panel identity band.
* Reads are internal fragments, uniform along a linear genome, forward strand, at ``--error-rate``
  substitutions. They are not poly-A 3' ends, and are not mutated copies (that is SENS-CAL-01).

Scoring: recovered molecules = count-matrix molecules in a condition's cells assigned to a gene of
the planted family (``family`` of ``virus_identity.tsv``), minus the same cells' count in
the unplanted twin when given. Capture = recovered / planted. ``lod95_true_molecules`` is
``minimum_molecules_for_detection(0.95, threshold) / capture``. Everything is one library
(n = 1), T1 truth, labelled exploratory, and carries no pass/fail threshold.

    python scripts/sens07_planted_sensitivity.py plant --out-dir D --technology 10xv3 \\
        --background-r1 R1.fastq.gz --background-r2 R2.fastq.gz \\
        --exact-fasta viral.fa --heldout-fasta heldout.fa --families families.tsv
    python scripts/sens07_planted_sensitivity.py score --out-dir S --truth D/truth.tsv \\
        --counts-h5ad run/kb-python/counts_unfiltered/adata_multimap.h5ad --layer X \\
        --identity run/results/virus_identity.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import plant_anello_10x as planter  # noqa: E402

from viralscan import chemistry as chem_mod  # noqa: E402
from viralscan.sensitivity import minimum_molecules_for_detection  # noqa: E402

DEFAULT_CELLS_PER_CONDITION = 20
DEFAULT_ERROR_RATE = 0.0
MODES = ("exact", "heldout")
LABELS = {
    "truth_tier": "T1 exact (planted molecules)",
    "status": "exploratory",
    "n_libraries": 1,
    "thresholds": None,
    "claim_limit": (
        "one library, internal-fragment reads, forward strand; no pass/fail threshold; "
        "not the VAL-01 panel and not a probit LOD (E8)"
    ),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ── inputs ────────────────────────────────────────────────────────────────────────────────────
def read_families(path: Path) -> dict[str, str]:
    """``accession -> family`` from a two-column TSV (a header ``accession`` row is skipped)."""
    out: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#") or row[0] == "accession":
                continue
            if len(row) != 2 or not row[1]:
                raise ValueError(f"{path.name}: expected accession<TAB>family, got {row}")
            out[row[0]] = row[1]
    return out


def geometry(
    technology: str, r1: Path, r2: Path, trim: bool = False
) -> tuple[int, int, int, int]:
    """``(cb_len, umi_len, read_len, r1_len)``: chemistry geometry, the background's own R2 length
    and R1 length. R1 shorter than barcode + UMI is always refused; longer only with *trim*."""
    cb_len, umi_len = chem_mod.cb_umi_geometry(technology)
    need = cb_len + umi_len
    for rec1, rec2 in planter.iter_background(r1, r2, 1):
        r1_len = len(rec1[1].strip())
        if r1_len < need or (r1_len > need and not trim):
            hint = (
                "it is too short to trim"
                if r1_len < need
                else "trim it or pass --trim-r1-to-geometry"
            )
            raise ValueError(
                f"{r1.name}: R1 is {r1_len} nt but {technology} needs {need} (barcode + UMI); {hint}"
            )
        return cb_len, umi_len, len(rec2[1].strip()), r1_len
    raise ValueError(f"{r1.name}: no read pairs")


def scan_background(
    r1: Path, r2: Path, n_pairs: int, cb_len: int, key_len: int, r1_len: int | None = None
):
    """``(barcodes ranked by read count, set of (CB, UMI), pairs read)`` for the first *n_pairs*.

    Every R1 must be *r1_len* nt (default *key_len*); only its first *key_len* bases are used."""
    r1_len = key_len if r1_len is None else r1_len
    counts: dict[str, int] = {}
    taken: set[tuple[str, str]] = set()
    n = 0
    for rec1, _ in planter.iter_background(r1, r2, n_pairs):
        seq = rec1[1].strip()
        if len(seq) != r1_len:
            raise ValueError(f"{r1.name}: an R1 read is {len(seq)} nt, expected {r1_len}")
        cb = seq[:cb_len]
        taken.add((cb, seq[cb_len:key_len]))
        if not set(cb) - set("ACGT"):
            counts[cb] = counts.get(cb, 0) + 1
        n += 1
    ranked = sorted(counts, key=lambda b: (-counts[b], b))
    return ranked, taken, n


def family_rng(seed: int, family: str) -> random.Random:
    """A generator seeded from sha256(seed, family): stable across runs, platforms and datasets."""
    h = hashlib.sha256(f"{seed}\t{family}".encode()).digest()
    return random.Random(int.from_bytes(h[:8], "big"))


def select_genomes(
    fastas: dict[str, dict[str, str]],
    families: dict[str, str],
    per_family: int,
    nearest: dict[str, tuple[float, str]],
    rng: random.Random,
    family_seed: int | None = None,
) -> list[dict[str, Any]]:
    """``per_family`` genomes for each (mode, family), drawn from the sorted accessions.

    With *family_seed* each family draws from its own ``family_rng(family_seed, family)`` instead of
    the shared *rng*, so its picks do not depend on which other families are present."""
    rows: list[dict[str, Any]] = []
    for mode in MODES:
        pools: dict[str, list[str]] = {}
        for acc in sorted(fastas.get(mode, {})):
            if acc in families:
                pools.setdefault(families[acc], []).append(acc)
        for family in sorted(pools):
            draw = rng if family_seed is None else family_rng(family_seed, family)
            for acc in sorted(draw.sample(pools[family], min(per_family, len(pools[family])))):
                seq = fastas[mode][acc]
                ident, target = nearest.get(acc, (None, ""))
                band = (planter.band(ident) if nearest else "unknown") if mode == "heldout" else ""
                rows.append(
                    {
                        "genome": acc,
                        "genus": family,
                        "identity_band": "exact"
                        if mode == "exact"
                        else f"heldout:{band or 'unknown'}",
                        "length": len(seq),
                        "orf1_start": 0,
                        "last_cds_end": len(seq),
                        "mode": mode,
                        "nearest_identity": round(ident, 4) if ident else "",
                        "nearest_panel_genome": target,
                    }
                )
    return rows


def read_gene_spans(path: Path) -> dict[str, tuple[str, int, int]]:
    """``gene_id -> (virus_id, 0-based start, end)`` from the AUD-18 ``--gene-spans`` TSV."""
    out: dict[str, tuple[str, int, int]] = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0] in ("virus_id",) or row[0].startswith("#"):
                continue
            if len(row) != 5:
                raise ValueError(
                    f"{path.name}: expected virus_id, gene_id, start, end, strand: {row}"
                )
            out[row[1]] = (row[0], int(row[2]) - 1, int(row[3]))
    return out


# ── plant ─────────────────────────────────────────────────────────────────────────────────────
def condition_blocks(selection, levels, cells: list[str], per_condition: int):
    """Disjoint cell blocks, one per (genome, "uniform", level), in selection order."""
    keys = [(r["genome"], "uniform", level) for r in selection for level in levels]
    if len(cells) < per_condition * len(keys):
        raise ValueError(
            f"{len(keys)} conditions x {per_condition} cells need {per_condition * len(keys)} "
            f"barcodes, the background has {len(cells)}; lower --cells-per-condition, --per-family "
            "or --levels, or use a bigger --n-background"
        )
    return {k: cells[i * per_condition : (i + 1) * per_condition] for i, k in enumerate(keys)}


def write_receipt(
    out_dir: Path,
    gene: str,
    accession: str,
    sequence: str,
    expected: float,
    index: Path,
    manifest_path: Path,
    count_layer: str,
    count_strategy: str,
) -> Path:
    """A checksum-bound exact-control receipt for *gene* (``sensitivity.verify_exact_control``)."""
    from viralscan.validation import require_schema_valid

    manifest = json.loads(manifest_path.read_text())
    rows = [r for r in manifest["sequences"] if r["accession_version"] == accession]
    if len(rows) != 1:
        raise ValueError(
            f"{manifest_path.name} must list {accession} exactly once, has {len(rows)}"
        )
    if rows[0]["sha256"] != hashlib.sha256(sequence.encode()).hexdigest() or rows[0][
        "length"
    ] != len(sequence):
        raise ValueError(f"planted exact sequence {accession} does not match {manifest_path.name}")
    out_dir.mkdir(parents=True, exist_ok=True)
    evidence = out_dir / f"{gene}.planting_truth.tsv"
    evidence.write_text(f"gene\texpected_molecules\n{gene}\t{expected:g}\n")
    receipt = {
        "schema_version": "3.0.0",
        "gene": gene,
        "accession_version": accession,
        "sequence_sha256": rows[0]["sha256"],
        "index_sha256": sha256_file(index),
        "reference_fasta_sha256": manifest["fasta_sha256"],
        "expected_molecules": expected,
        "count_layer": count_layer,
        "count_strategy": count_strategy,
        "reference_manifest": {
            "path": str(manifest_path.resolve()),
            "sha256": sha256_file(manifest_path),
        },
        "evidence": {"path": evidence.name, "sha256": sha256_file(evidence)},
    }
    path = out_dir / f"{gene}.receipt.json"
    require_schema_valid(receipt, "positive_control_receipt.schema.json", path)
    path.write_text(json.dumps(receipt, indent=2) + "\n")
    return path


def plant(args: argparse.Namespace) -> dict[str, Any]:
    rng = random.Random(args.seed)
    out: Path = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    trim = bool(getattr(args, "trim_r1_to_geometry", False))
    cb_len, umi_len, read_len, r1_len = geometry(
        args.technology, args.background_r1, args.background_r2, trim
    )
    cells, taken, n_bg = scan_background(
        args.background_r1,
        args.background_r2,
        args.n_background or sys.maxsize,
        cb_len,
        cb_len + umi_len,
        r1_len,
    )

    def records(path):
        return {k: v.upper() for k, v in planter.read_fasta(path).items()} if path else {}

    fastas = {"exact": records(args.exact_fasta), "heldout": records(args.heldout_fasta)}
    both = set(fastas["exact"]) & set(fastas["heldout"])
    if both:
        raise ValueError(
            f"held-out genomes must be absent from the panel; in both FASTAs: {sorted(both)[:5]}"
        )
    nearest = planter.nearest_identity(args.paf) if args.paf else {}
    seed_mode = getattr(args, "selection_seed_mode", "global")
    selection = select_genomes(
        fastas,
        read_families(args.families),
        args.per_family,
        nearest,
        rng,
        args.seed if seed_mode == "per-family" else None,
    )
    if not selection:
        raise ValueError("no genome selected: check that --families lists accessions of the FASTAs")
    genomes = {r["genome"]: {"seq": fastas[r["mode"]][r["genome"]]} for r in selection}

    controls: list[dict[str, Any]] = []
    if args.control_genes:
        spans = read_gene_spans(args.gene_spans) if args.gene_spans else {}
        parents = {r["genome"]: r for r in selection if r["mode"] == "exact"}
        for gene in args.control_genes:
            if gene not in spans:
                raise ValueError(f"--control-genes {gene!r} is not in --gene-spans")
            virus, start, end = spans[gene]
            parent = parents.get(virus)
            if parent is None:
                raise ValueError(
                    f"control gene {gene} belongs to {virus}, not a selected exact genome"
                )
            name = f"{virus}|{gene}"
            body = genomes[virus]["seq"][start:end]
            if len(body) < read_len:
                raise ValueError(
                    f"control gene {gene} is {len(body)} bp, shorter than a {read_len} nt read"
                )
            genomes[name] = {"seq": body}
            # The control's gene count must be only the control: plant its parent genome nowhere else.
            selection = [r for r in selection if r["genome"] != virus]
            selection.append(
                {
                    **parent,
                    "genome": name,
                    "identity_band": "exact_control",
                    "length": len(body),
                    "last_cds_end": len(body),
                    "mode": "exact_control",
                }
            )
            controls.append({"gene": gene, "accession": virus, "genome": name})

    levels = tuple(args.levels)
    blocks = condition_blocks(selection, levels, cells, args.cells_per_condition)
    with open(out / "selection.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(selection[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(selection)

    planted, unplanted = out / "planted", out / "unplanted"
    planted.mkdir(exist_ok=True)
    unplanted.mkdir(exist_ok=True)
    with (
        gzip.open(planted / "R1.fastq.gz", "wt") as p1,
        gzip.open(planted / "R2.fastq.gz", "wt") as p2,
        gzip.open(unplanted / "R1.fastq.gz", "wt") as u1,
        gzip.open(unplanted / "R2.fastq.gz", "wt") as u2,
        open(out / "truth.tsv", "w", encoding="utf-8", newline="") as tf,
    ):
        planter.copy_background(
            args.background_r1,
            args.background_r2,
            args.n_background or sys.maxsize,
            [(p1, p2), (u1, u2)],
            r1_len=cb_len + umi_len if trim else None,
        )
        truth = csv.writer(tf, delimiter="\t")
        truth.writerow(
            [
                "read_id",
                "genome",
                "genus",
                "identity_band",
                "window",
                "level",
                "cb",
                "umi",
                "pos",
                "body_len",
                "tail_len",
            ]
        )
        n_planted = planter.plant_reads(
            selection,
            genomes,
            None,
            rng,
            p1,
            p2,
            truth,
            read_len=read_len,
            umi_len=umi_len,
            levels=levels,
            windows=("uniform",),
            error_rate=args.error_rate,
            barcode_blocks=blocks,
            taken=taken,
            circular=False,
        )

    arg_rows = []
    for control in controls:
        expected = float(sum(levels))
        row = {**control, "expected_molecules": expected, "receipt": "", "args": ""}
        words = [
            "--positive-control-gene",
            control["gene"],
            "--positive-control-molecules",
            f"{expected:g}",
        ]
        if args.index and args.reference_manifest and args.identity:
            from viralscan.virus_identity import VirusIdentityTable

            table = VirusIdentityTable.read_tsv(args.identity)
            gene_row = table.by_gene().get(control["gene"])
            if gene_row is None or not gene_row.virus_key:
                raise ValueError(
                    f"{control['gene']} is not an indexed viral gene of {args.identity.name}"
                )
            receipt = write_receipt(
                out / "controls",
                control["gene"],
                control["accession"],
                fastas["exact"][control["accession"]],
                expected,
                args.index,
                args.reference_manifest,
                args.count_layer,
                args.count_strategy,
            )
            row["receipt"] = str(receipt)
            words += [
                "--positive-control-scope",
                "exact_sequence",
                "--positive-control-virus-key",
                gene_row.virus_key,
                "--positive-control-receipt",
                str(receipt),
            ]
        row["args"] = " ".join(words)
        arg_rows.append(row)
    if arg_rows:
        with open(out / "positive_control_args.tsv", "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(arg_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(arg_rows)

    manifest = {
        **LABELS,
        "seed": args.seed,
        "technology": args.technology,
        "cb_len": cb_len,
        "umi_len": umi_len,
        "read_length": read_len,
        "selection_seed_mode": seed_mode,
        "r1_trim": {
            "trim_r1_to_geometry": trim,
            "original_length": r1_len,
            "trimmed_length": cb_len + umi_len,
            "applied_to": ["planted", "unplanted"],
        },
        "read_changes": (
            [
                f"R1 trimmed from {r1_len} to {cb_len + umi_len} nt (barcode + UMI), "
                "identically in planted and unplanted"
            ]
            if trim
            else []
        ),
        "error_rate": args.error_rate,
        "levels": list(levels),
        "cells_per_condition": args.cells_per_condition,
        "n_background_pairs": n_bg,
        "n_planted_reads": n_planted,
        "n_conditions": len(blocks),
        "windows": ["uniform"],
        "genomes_linear": True,
        "inputs": {
            k: {"name": v.name, "sha256": sha256_file(v)}
            for k, v in {
                "background_r1": args.background_r1,
                "background_r2": args.background_r2,
                "exact_fasta": args.exact_fasta,
                "heldout_fasta": args.heldout_fasta,
                "families": args.families,
                "paf": args.paf,
            }.items()
            if v
        },
        "controls": arg_rows,
    }
    (out / "plant_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


# ── score ─────────────────────────────────────────────────────────────────────────────────────
def load_counts(path: Path, layer: str = "X") -> dict[tuple[str, str], float]:
    """``(barcode, gene) -> molecules`` from an h5ad *layer* or a kb ``counts_unfiltered`` dir."""
    counts: dict[tuple[str, str], float] = {}
    if path.is_dir():
        barcodes = (path / "cells_x_genes.barcodes.txt").read_text().split()
        genes = (path / "cells_x_genes.genes.txt").read_text().split()
        with open(path / "cells_x_genes.mtx", encoding="utf-8") as fh:
            lines = (ln for ln in fh if not ln.startswith("%"))
            next(lines)  # the size line
            for line in lines:
                i, j, value = line.split()
                counts[(barcodes[int(i) - 1], genes[int(j) - 1])] = float(value)
        return counts
    import anndata as ad

    adata = ad.read_h5ad(path)
    matrix = adata.X if layer == "X" else adata.layers[layer]
    import scipy.sparse as sp

    coo = matrix.tocoo() if sp.issparse(matrix) else sp.coo_matrix(matrix)
    for i, j, v in zip(coo.row, coo.col, coo.data):
        counts[(str(adata.obs_names[i]), str(adata.var_names[j]))] = float(v)
    return counts


def wilson(k: float, n: float, z: float = 1.959964) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion (planted molecules treated as trials)."""
    if n <= 0:
        return math.nan, math.nan
    p = min(1.0, k / n)
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, centre - half), min(1.0, centre + half)


def recovered(
    counts, cells: set[str], gene_group: dict[str, str], group: str
) -> tuple[float, float]:
    """``(molecules of *group*, molecules of any other viral group)`` in *cells*."""
    mine = other = 0.0
    for (cb, gene), value in counts.items():
        if cb in cells and gene in gene_group:
            if gene_group[gene] == group:
                mine += value
            else:
                other += value
    return mine, other


def score(
    truth_rows: list[dict[str, str]],
    counts,
    gene_group: dict[str, str],
    twin=None,
    threshold: int = 1,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """``(per-condition rows, per-(family, mode) summary rows)``."""
    conditions: dict[tuple[str, str, str, int], dict[str, Any]] = {}
    for row in truth_rows:
        key = (row["genome"], row["genus"], row["identity_band"], int(row["level"]))
        cond = conditions.setdefault(key, {"cells": set(), "planted": 0})
        cond["cells"].add(row["cb"])
        cond["planted"] += 1
    rows = []
    for (genome, group, band, level), cond in sorted(conditions.items()):
        mine, other = recovered(counts, cond["cells"], gene_group, group)
        background = (
            recovered(twin, cond["cells"], gene_group, group)[0] if twin is not None else None
        )
        net = max(0.0, mine - (background or 0.0))
        rows.append(
            {
                "genome": genome,
                "family": group,
                "mode": band.split(":")[0],
                "identity_band": band,
                "level": level,
                "planted": cond["planted"],
                "recovered_raw": mine,
                "background_in_twin": "" if background is None else background,
                "recovered": net,
                "recovery": net / cond["planted"],
                "other_viral": other,
                "detected": net >= threshold,
            }
        )
    summary = []
    for family, mode in sorted({(r["family"], r["mode"]) for r in rows}):
        sel = [r for r in rows if (r["family"], r["mode"]) == (family, mode)]
        planted, got = sum(r["planted"] for r in sel), sum(r["recovered"] for r in sel)
        capture = got / planted
        lo, hi = wilson(got, planted)
        status = "failed" if capture == 0 else "over_recovered" if capture > 1 else "measured"
        lod = (
            minimum_molecules_for_detection(0.95, threshold) / capture
            if status == "measured"
            else ""
        )
        by_level: dict[int, list[bool]] = {}
        for r in sel:
            by_level.setdefault(r["level"], []).append(r["detected"])
        all_found = [lv for lv, found in sorted(by_level.items()) if all(found)]
        summary.append(
            {
                "family": family,
                "mode": mode,
                "n_genomes": len({r["genome"] for r in sel}),
                "planted": planted,
                "recovered": got,
                "capture": capture,
                "capture_ci_low": lo,
                "capture_ci_high": hi,
                "capture_status": status,
                "lod95_true_molecules": lod,
                "lowest_level_all_detected": all_found[0] if all_found else "",
                "detection_threshold": threshold,
            }
        )
    return rows, summary


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else [], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_score(args: argparse.Namespace) -> dict[str, Any]:
    from viralscan.virus_identity import VirusIdentityTable

    if args.group_column != "family":
        raise ValueError("SENS-07 truth records family labels; --group-column must be family")
    table = VirusIdentityTable.read_tsv(args.identity)
    gene_group = {g.gene_id: str(getattr(g, args.group_column)) for g in table.genes if g.viral}
    with open(args.truth, encoding="utf-8") as fh:
        truth_rows = list(csv.DictReader(fh, delimiter="\t"))
    counts = load_counts(args.counts_h5ad, args.layer)
    twin = load_counts(args.twin_h5ad, args.layer) if args.twin_h5ad else None
    rows, summary = score(truth_rows, counts, gene_group, twin, args.threshold)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_dir / "recovery_by_condition.tsv", rows)
    write_tsv(args.out_dir / "family_summary.tsv", summary)
    report = {
        **LABELS,
        "count_layer": args.layer,
        "group_column": args.group_column,
        "twin_subtracted": twin is not None,
        "families": summary,
    }
    (args.out_dir / "sens07_scores.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("plant", help="plant exact and held-out genomes into a background library")
    a.add_argument("--out-dir", type=Path, required=True)
    a.add_argument("--technology", default="10xv3")
    a.add_argument("--background-r1", type=Path, required=True)
    a.add_argument("--background-r2", type=Path, required=True)
    a.add_argument("--exact-fasta", type=Path)
    a.add_argument("--heldout-fasta", type=Path)
    a.add_argument("--families", type=Path, required=True, help="accession<TAB>family")
    a.add_argument("--paf", type=Path, help="minimap2 PAF of held-out vs panel (identity band)")
    a.add_argument("--per-family", type=int, default=1)
    a.add_argument("--levels", type=int, nargs="+", default=list(planter.LEVELS))
    a.add_argument("--cells-per-condition", type=int, default=DEFAULT_CELLS_PER_CONDITION)
    a.add_argument("--n-background", type=int, default=0, help="first N pairs (0 = all)")
    a.add_argument("--error-rate", type=float, default=DEFAULT_ERROR_RATE)
    a.add_argument("--seed", type=int, default=20261012)
    a.add_argument(
        "--selection-seed-mode",
        choices=["global", "per-family"],
        default="global",
        help="global: one RNG over all families (round 2b); per-family: each family's genomes "
        "drawn from sha256(seed, family), identical across datasets (round 2c)",
    )
    a.add_argument(
        "--trim-r1-to-geometry",
        action="store_true",
        help="R1 longer than barcode + UMI: trim it to that length in planted and unplanted "
        "(a declared read change); shorter is always refused",
    )
    a.add_argument(
        "--gene-spans", type=Path, help="virus_id, gene_id, start, end, strand (1-based)"
    )
    a.add_argument(
        "--control-genes", nargs="+", default=[], help="gene IDs planted as exact controls"
    )
    a.add_argument("--index", type=Path)
    a.add_argument("--reference-manifest", type=Path)
    a.add_argument("--identity", type=Path, help="virus_identity.tsv of the reference")
    a.add_argument("--count-layer", default="X")
    a.add_argument("--count-strategy", default="kb_count")
    s = sub.add_parser("score", help="score recovery from a ViralScan run on planted/")
    s.add_argument("--out-dir", type=Path, required=True)
    s.add_argument("--truth", type=Path, required=True)
    s.add_argument(
        "--counts-h5ad", type=Path, required=True, help="h5ad, or a kb counts_unfiltered dir"
    )
    s.add_argument("--twin-h5ad", type=Path, help="the same for the unplanted twin")
    s.add_argument(
        "--layer", default="X", help="count layer; use the one --positive-control-* measures"
    )
    s.add_argument("--identity", type=Path, required=True)
    s.add_argument(
        "--group-column",
        default="family",
        choices=["family"],
        help="match the truth's family labels",
    )
    s.add_argument(
        "--threshold", type=int, default=1, help="molecules for a condition to count as detected"
    )
    return p


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    result = plant(args) if args.command == "plant" else run_score(args)
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()

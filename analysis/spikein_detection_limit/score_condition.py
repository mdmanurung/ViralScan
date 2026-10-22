"""Score one condition's result bundle: one row per arm in {old, flank, star} (PLAN SPIKE-02). Standard library only.

Bundle layout (written by ``run_condition.py``)::

    truth.tsv                         spikein.py truth table (name, cb, class, contigs, ...); absent/empty = unspiked
    star/star_spike_fates.tsv         name, on_truth (0/1: ANY alignment on ANY truth contig), fate
    star/star_contig_reads.tsv        contig, n_any, n_primary for NON-spike reads (the STAR false-positive floor)
    <arm>/kb_molecules.tsv.gz         barcode, gene_id, unique, ambiguous: non-zero viral cells of adata_multimap.h5ad
    <arm>/virus_identity.tsv          gene_id -> genome_accession (BAM contig = VIRUS|<genome_accession>)
    <arm>/multimap_evidence.tsv       the product's per-gene table (evidence tiers; over ALL barcodes)

kb arm count = viral_molecules_unique + viral_molecules_ambiguous_allocated (the multimap step's own columns:
adata layers ``counts_unique_viral`` and ``counts_ambiguous_allocated``) summed over the spike barcodes and the genes
of the truth accessions. Call = kb molecules >= KB_MIN AND STAR reads >= 1 on a truth contig (any of the set).
"unsupported" viral molecules (null-library false-positive floor) = kb molecules, all barcodes, on accessions outside
the truth set that have NO non-spike STAR read.
"""

import argparse
import collections
import csv
import gzip
import sys
from pathlib import Path

KB_MIN = 3
STAR_MIN = 1
TIER_RANK = {
    "candidate_unique": 3,
    "candidate_virus_ambiguous": 2,
    "candidate_host_virus_ambiguous": 1,
    "not_detected": 0,
}
COLUMNS = [
    "cond_id", "arm", "background", "source", "level", "seed", "layout", "cells", "bg_rep",
    "n_spiked", "kb_unique", "kb_ambiguous_allocated", "kb_total", "kb_unique_only", "evidence_tier",
    "evidence_truth_total_all_barcodes", "kb_total_in_baseline",
    "star_n_source_contig_any", "star_n_virus", "star_n_other_virus", "star_n_host", "star_n_absent",
    "call_kb", "call_star", "call_combined",
    "kb_unsupported_molecules", "kb_total_viral_molecules", "n_false_calls",
    "star_nonspike_reads_total", "star_contigs_with_reads",
]  # fmt: skip


def read_tsv(path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_molecules(path):
    with _open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            yield r["barcode"], r["gene_id"], float(r["unique"]), float(r["ambiguous"])


def kb_stats(arm_dir: Path, spike_cbs: set, truth_acc: set, star_any: dict) -> dict:
    """Molecule sums for one kb arm directory (see module doc)."""
    acc_of = {r["gene_id"]: r["genome_accession"] for r in read_tsv(arm_dir / "virus_identity.tsv")}
    per_acc, uniq, amb = collections.Counter(), 0.0, 0.0
    for cb, gene, u, a in read_molecules(arm_dir / "kb_molecules.tsv.gz"):
        acc = acc_of.get(gene)
        if acc is None:
            continue
        per_acc[acc] += u + a
        if acc in truth_acc and cb in spike_cbs:
            uniq, amb = uniq + u, amb + a
    other = {a: m for a, m in per_acc.items() if a not in truth_acc}
    return dict(
        uniq=uniq,
        amb=amb,
        unsupported=sum(m for a, m in other.items() if star_any.get(a, 0) < STAR_MIN),
        total=sum(per_acc.values()),
        false_calls=sum(1 for a, m in other.items() if m >= KB_MIN and star_any.get(a, 0) >= STAR_MIN),
    )  # fmt: skip


def evidence_tier(arm_dir: Path, truth_acc: set):
    """(best evidence tier over the truth genes, their viral_molecules_total_est summed over all barcodes)."""
    acc_of = {r["gene_id"]: r["genome_accession"] for r in read_tsv(arm_dir / "virus_identity.tsv")}
    best, total = "not_detected", 0.0
    for r in read_tsv(arm_dir / "multimap_evidence.tsv"):
        if acc_of.get(r["gene_id"]) in truth_acc:
            total += float(r["viral_molecules_total_est"] or 0)
            if TIER_RANK.get(r["evidence_tier"], 0) > TIER_RANK[best]:
                best = r["evidence_tier"]
    return best, total


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return str(v)
    return f"{v:.6g}" if isinstance(v, float) else str(v)


def score_bundle(bundle: Path, cond: dict, baseline: Path | None = None) -> list[dict]:
    bundle = Path(bundle)
    truth = read_tsv(bundle / "truth.tsv") if (bundle / "truth.tsv").exists() else []
    spiked = bool(truth)
    truth_acc = (
        {c.removeprefix("VIRUS|") for c in truth[0]["contigs"].split(",") if c} if spiked else set()
    )
    spike_cbs = {r["cb"] for r in truth}
    names = {r["name"] for r in truth}
    fates = {r["name"]: r for r in read_tsv(bundle / "star/star_spike_fates.tsv")}
    star_any = {
        r["contig"].removeprefix("VIRUS|"): int(r["n_any"])
        for r in read_tsv(bundle / "star/star_contig_reads.tsv")
    }
    on_truth = sum(1 for n in names if fates.get(n, {}).get("on_truth") == "1")
    cnt = collections.Counter(fates[n]["fate"] if n in fates else "absent" for n in names)
    call_star = (on_truth >= STAR_MIN) if spiked else None
    base = {k: cond.get(k, "") for k in ("cond_id", "background", "source", "level", "seed", "layout", "cells", "bg_rep")}  # fmt: skip
    base["n_spiked"] = len(names)
    rows = []
    for arm in ("old", "flank"):
        s = kb_stats(bundle / arm, spike_cbs, truth_acc, star_any)
        tier, ev_total = evidence_tier(bundle / arm, truth_acc) if spiked else (None, None)
        bg = None
        if spiked and baseline is not None:
            b = kb_stats(Path(baseline) / arm, spike_cbs, truth_acc, star_any)
            bg = b["uniq"] + b["amb"]
        kb_total = s["uniq"] + s["amb"]
        call_kb = (kb_total >= KB_MIN) if spiked else None
        rows.append(
            {
                **base, "arm": arm,
                "kb_unique": s["uniq"] if spiked else None,
                "kb_ambiguous_allocated": s["amb"] if spiked else None,
                "kb_total": kb_total if spiked else None,
                "kb_unique_only": s["uniq"] if spiked else None,
                "evidence_tier": tier, "evidence_truth_total_all_barcodes": ev_total,
                "kb_total_in_baseline": bg,
                "star_n_source_contig_any": on_truth if spiked else None,
                "call_kb": call_kb, "call_star": call_star,
                "call_combined": (call_kb and call_star) if spiked else None,
                "kb_unsupported_molecules": s["unsupported"],
                "kb_total_viral_molecules": s["total"], "n_false_calls": s["false_calls"],
            }
        )  # fmt: skip
    rows.append(
        {
            **base, "arm": "star",
            "star_n_source_contig_any": on_truth if spiked else None,
            "star_n_virus": (cnt["source"] + cnt["other_virus"]) if spiked else None,
            "star_n_other_virus": cnt["other_virus"] if spiked else None,
            "star_n_host": cnt["host"] if spiked else None,
            "star_n_absent": cnt["absent"] if spiked else None,
            "call_star": call_star,
            "star_nonspike_reads_total": sum(star_any.values()),
            "star_contigs_with_reads": sum(1 for v in star_any.values() if v > 0),
        }
    )  # fmt: skip
    return rows


def write_scores(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLUMNS)
        w.writerows([fmt(r.get(c)) for c in COLUMNS] for r in rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--bundle", required=True)
    p.add_argument("--conditions", required=True)
    p.add_argument("--cond-id", required=True)
    p.add_argument(
        "--baseline-bundle", help="bundle of the unspiked replicate on the same background"
    )
    p.add_argument("--out", help="default: <bundle>/scores.tsv")
    a = p.parse_args(argv)
    cond = next((r for r in read_tsv(a.conditions) if r["cond_id"] == a.cond_id), None)
    if cond is None:
        sys.exit(f"score_condition: {a.cond_id} not in {a.conditions}")
    rows = score_bundle(
        Path(a.bundle), cond, Path(a.baseline_bundle) if a.baseline_bundle else None
    )
    write_scores(rows, a.out or Path(a.bundle) / "scores.tsv")


if __name__ == "__main__":
    main()

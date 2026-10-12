#!/usr/bin/env python3
"""VAL-08 golden-run tooling: input preparation, the run command lines and the expected-output checks.

Code only. Nothing here runs ``kb count`` or ViralScan; ``plan`` prints the command lines the owner
runs on the cluster, and ``check`` / ``compare`` judge what came back. Exactness is the definition of
VAL-08 (planted molecules recovered exactly, counts conserved, outputs invariant), not a threshold:
this file adds no acceptance numbers.

Subcommands
-----------
``prepare``   reference FASTA + GTF (+ contig -> virus-id map) -> the ``--viral-fasta`` / ``--gene-spans``
              inputs of ``generate_truth_panel.py`` (validated with the generator's own loader).
``whitelist`` one sample's ``barcode_pool.tsv`` -> a ``-w`` barcode list.
``permute``   one sample's FASTQ pair -> a record-reversed, multi-member-gzip copy (read-order and
              transport-chunk variants carry exactly the same read pairs).
``plan``      the ordered command lines (generator, whitelist, ViralScan on the original and permuted
              FASTQs, ``check``, ``compare``).
``check``     one ViralScan sample directory against the panel truth.
``compare``   two sample directories must carry identical molecule assignments and counts.

``check`` and ``compare`` exit 1 and name every failed check. ``--reference real`` is for a reference
that is not the generator's own templates (the cat42b arm): host molecules then need not be indexed, so
only "no host molecule is called viral" is checked on the host side, not that the host set is exact.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shlex
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generate_truth_panel as gt  # noqa: E402
import v3_score as score  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
GENE_ID = re.compile(r'gene_id "([^"]+)"')


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ── prepare: reference -> generator inputs ───────────────────────────────────


def spans_from_gtf(gtf: Path, contig_to_virus: dict[str, str] | None = None) -> dict[str, Any]:
    """Gene spans of a viral GTF as generator rows; spliced or unstranded genes are skipped, listed.

    The generator plants a contiguous genomic window, so a gene whose exons are not one block would
    not match the spliced transcript the index holds. ``contig_to_virus`` renames contigs (default:
    unchanged); a contig it does not list is not a panel virus and is ignored.
    """
    blocks: dict[tuple[str, str], list[tuple[int, int]]] = {}
    strands: dict[tuple[str, str], set[str]] = {}
    for line in gtf.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("exon", "CDS"):
            continue
        virus = f[0] if contig_to_virus is None else contig_to_virus.get(f[0])
        found = GENE_ID.search(f[8])
        if virus is None or found is None:
            continue
        key = (virus, found.group(1))
        blocks.setdefault(key, []).append((int(f[3]), int(f[4])))
        strands.setdefault(key, set()).add(f[6])
    rows, skipped = [], []
    for (virus, gene), parts in sorted(blocks.items()):
        merged: list[list[int]] = []
        for first, last in sorted(parts):
            if merged and first <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], last)
            else:
                merged.append([first, last])
        sign = strands[(virus, gene)]
        if (
            len(merged) != 1
            or len(sign) != 1
            or not sign <= {"+", "-"}
            or merged[0][1] - merged[0][0] + 1 < gt.MIN_BODY
        ):
            skipped.append(gene)
        else:
            rows.append(
                dict(
                    virus_id=virus,
                    gene_id=gene,
                    start=merged[0][0],
                    end=merged[0][1],
                    strand=next(iter(sign)),
                )
            )
    return dict(rows=rows, skipped=skipped)


def prepare(
    fasta: Path,
    gtf: Path,
    out_fasta: Path,
    out_spans: Path,
    contig_to_virus: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Write the renamed FASTA and gene-span TSV, and validate both with the generator's loader."""
    genomes = gt.read_fasta_records(fasta)
    if contig_to_virus is not None:
        genomes = {contig_to_virus[c]: s for c, s in genomes.items() if c in contig_to_virus}
    found = spans_from_gtf(gtf, contig_to_virus)
    by_virus = Counter(r["virus_id"] for r in found["rows"])
    keep = {
        v for v, n in by_virus.items() if n >= 2 and v in genomes
    }  # the loader needs >= 2 genes
    out_fasta.write_text("".join(f">{v}\n{genomes[v]}\n" for v in sorted(keep)))
    out_spans.write_text(
        "".join(
            "\t".join(str(r[c]) for c in ("virus_id", "gene_id", "start", "end", "strand")) + "\n"
            for r in found["rows"]
            if r["virus_id"] in keep
        )
    )
    gt.load_real_inputs(out_fasta, out_spans, None)  # raises on anything the generator would refuse
    return dict(viruses=sorted(keep), dropped_viruses=sorted(set(genomes) - keep), **found)


# ── permuted FASTQs and whitelist ────────────────────────────────────────────


def permute_fastq(source: Path, dest: Path, chunk: int = 127) -> Path:
    """Record-reversed copy of ``source`` written as ``chunk``-record gzip members (``mtime=0``)."""
    with gzip.open(source, "rt") as handle:
        lines = handle.readlines()
    records = ["".join(lines[i : i + 4]) for i in range(0, len(lines), 4)][::-1]
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(
        b"".join(
            gzip.compress("".join(records[i : i + chunk]).encode(), mtime=0)
            for i in range(0, len(records), chunk)
        )
    )
    return dest


def write_whitelist(panel: Path, sample: str, out: Path) -> int:
    pool = read_tsv(panel / sample / "barcode_pool.tsv")
    out.write_text("".join(row["cell_barcode"] + "\n" for row in pool))
    return len(pool)


# ── check ────────────────────────────────────────────────────────────────────


def gene_virus_map(panel: Path) -> dict[str, str]:
    return {
        r["gene_id"]: r["virus_id"] for r in read_tsv(panel / "templates" / "viral_gene_models.tsv")
    }


def read_assignments(sample_dir: Path) -> list[dict[str, str]]:
    with gzip.open(sample_dir / "results" / "molecule_assignments.tsv.gz", "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def assigned_sets(assignments: list[dict[str, str]], gene_virus: dict[str, str]):
    """Viral ``(cb, umi, virus)``, host ``(cb, umi)`` and unresolved ``(cb, umi)`` keys."""
    viral: set[tuple[str, str, str]] = set()
    host: set[tuple[str, str]] = set()
    unresolved: set[tuple[str, str]] = set()
    for row in assignments:
        key = row["barcode"], row["umi"]
        if row["status"] == "unresolved":
            unresolved.add(key)
            continue
        for gene, weight in zip(row["genes"].split(";"), row["weights"].split(";")):
            if float(weight) > 0:
                if gene in gene_virus:
                    viral.add((*key, gene_virus[gene]))
                else:
                    host.add(key)
    return viral, host, unresolved


def truth_of(panel: Path, sample: str, truth_dir: Path | None = None):
    """``(planted, host, collision_kind)``: planted ``(cb, umi, virus)``, host ``(cb, umi)`` keys."""
    base = truth_dir or panel / "truth"
    rows = [r for r in read_tsv(base / "truth_manifest.tsv") if r["biological_sample_id"] == sample]
    host_rows = [
        r for r in read_tsv(base / "host_only_manifest.tsv") if r["biological_sample_id"] == sample
    ]
    planted = {
        (r["cell_barcode"], r["umi"], r["planted_virus_id"])
        for r in rows
        if r["is_planted_molecule"] == "true"
    }
    host = {
        (r["cell_barcode"], r["umi"]) for r in rows if r["record_granularity"] == "molecule"
    } | {(r["cell_barcode"], r["umi"]) for r in host_rows}
    kinds = {
        (r["cell_barcode"], r["umi"], r["planted_virus_id"]): r["collision_kind"]
        for r in rows
        if r["is_planted_molecule"] == "true"
    }
    if not rows and not host_rows:
        raise ValueError(f"no truth rows for sample {sample!r} under {base}")
    return planted, host, kinds


def count_matrix(adata: Any) -> Any:
    """Counts as a dense array over sorted barcodes and features (row and column order erased)."""
    selected = adata[sorted(adata.obs_names), sorted(adata.var_names)].X
    return selected.toarray() if hasattr(selected, "toarray") else selected


def conservation_failures(adata: Any, assignments: list[dict[str, str]]) -> list[str]:
    import numpy as np

    audit = dict(adata.uns["molecule_audit"])
    status = Counter(row["status"] for row in assignments)
    failures = []

    def expect(label: str, got: float, want: float) -> None:
        if got != want:
            failures.append(f"conservation: {label} is {got}, expected {want}")

    expect("input_molecules vs assignment rows", audit["input_molecules"], len(assignments))
    expect(
        "resolved + unresolved vs input_molecules",
        audit["resolved_molecules"] + audit["unresolved_molecules"],
        audit["input_molecules"],
    )
    expect(
        "unique + ambiguous vs resolved_molecules",
        audit["unique_molecules"] + audit["ambiguous_molecules"],
        audit["resolved_molecules"],
    )
    expect("unresolved_molecules vs rows", audit["unresolved_molecules"], status["unresolved"])
    expect("unique_molecules vs rows", audit["unique_molecules"], status["unique"])
    expect("ambiguous_molecules vs rows", audit["ambiguous_molecules"], status["ambiguous"])
    if abs(float(adata.X.sum()) - audit["resolved_molecules"]) > 1e-6:
        failures.append("conservation: sum of X differs from resolved_molecules")
    layers = (adata.layers["counts_unique"] + adata.layers["counts_ambiguous_allocated"]).toarray()
    if not np.allclose(adata.X.toarray(), layers):
        failures.append("conservation: X differs from counts_unique + counts_ambiguous_allocated")
    return failures


def check_run(
    panel: Path,
    sample: str,
    sample_dir: Path,
    *,
    reference: str = "panel",
    allow_collision_loss: bool = False,
    truth_dir: Path | None = None,
) -> dict[str, Any]:
    """Judge one ViralScan sample directory against the panel truth; ``failures`` is empty if exact.

    ``allow_collision_loss`` is for a panel generated with CB-UMI collisions: a planted molecule may
    be missing only if its (cb, umi) is unresolved and its truth row is a same-key collision (kept in
    recall's denominator, never dropped).
    """
    import anndata as ad
    import numpy as np

    if reference not in ("panel", "real"):
        raise ValueError("reference must be 'panel' or 'real'")
    planted, host_truth, kinds = truth_of(panel, sample, truth_dir)
    assignments = read_assignments(sample_dir)
    viral, host, unresolved = assigned_sets(assignments, gene_virus_map(panel))
    metric = score.precision_recall_f1(planted, viral)
    failures: list[str] = []
    lost = planted - viral
    if metric["fp"]:
        failures.append(f"{metric['fp']} viral calls are not planted molecules")
    if lost:
        explained = allow_collision_loss and all(
            (cb, umi) in unresolved and kinds[(cb, umi, v)] == "same" for cb, umi, v in lost
        )
        if not explained:
            failures.append(f"{len(lost)} planted molecules were not recovered")
    if host_truth & {(cb, umi) for cb, umi, _ in viral}:
        failures.append("host molecules were called viral")
    if host & {(cb, umi) for cb, umi, _ in planted}:
        failures.append("planted viral molecules were called host")
    if reference == "panel":
        if host != host_truth:
            failures.append("the called host molecules differ from the host truth")
        known = host_truth | {(cb, umi) for cb, umi, _ in planted}
        if {(row["barcode"], row["umi"]) for row in assignments} != known:
            failures.append("the molecules that entered counting differ from the truth molecules")
    adata = ad.read_h5ad(sample_dir / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad")
    failures += conservation_failures(adata, assignments)
    return dict(
        sample=sample,
        failures=failures,
        metric=metric,
        n_lost=len(lost),
        n_unresolved=len(unresolved),
        audit={
            k: float(v)
            for k, v in dict(adata.uns["molecule_audit"]).items()
            if isinstance(v, (int, float, np.integer, np.floating))
        },
    )


def compare_runs(first: Path, second: Path) -> list[str]:
    """Differences between two sample directories; none means identical results (determinism,
    read-order and chunk-size invariance are all this comparison on different inputs)."""
    import anndata as ad
    import numpy as np

    def canonical(sample_dir: Path):
        return sorted(
            (
                row["barcode"],
                row["umi"],
                row["status"],
                tuple(sorted(zip(row["genes"].split(";"), row["weights"].split(";")))),
            )
            for row in read_assignments(sample_dir)
        )

    path = Path("kb-python") / "counts_unfiltered" / "adata_multimap.h5ad"
    a, b = ad.read_h5ad(first / path), ad.read_h5ad(second / path)
    failures = []
    if canonical(first) != canonical(second):
        failures.append("molecule assignments differ")
    if sorted(a.obs_names) != sorted(b.obs_names) or sorted(a.var_names) != sorted(b.var_names):
        failures.append("barcode or feature universe differs")
    elif not np.array_equal(count_matrix(a), count_matrix(b)):
        failures.append("count matrices differ")
    if dict(a.uns["molecule_audit"]) != dict(b.uns["molecule_audit"]):
        failures.append("molecule audits differ")
    return failures


# ── plan: the command lines ──────────────────────────────────────────────────


def plan(
    *,
    out: Path,
    fasta: Path,
    spans: Path,
    onlist: Path,
    index: Path,
    t2g: Path,
    gtf: Path,
    sample: str,
    seed: int,
    chemistry: str = "10xv3",
    reference: str = "real",
) -> list[list[str]]:
    """Ordered argv lists: generate (collision-free, error-free golden), whitelist, run the original
    and the permuted FASTQs, check, compare. ``sample`` is read from the panel manifest after step 1."""
    script = str(Path(__file__).resolve())
    generator = str(Path(gt.__file__).resolve())
    panel, work = out / "panel", out / "work"
    r1, r2 = (panel / sample / f"{sample}_R{m}.fastq.gz" for m in (1, 2))
    p1, p2 = (work / "permuted" / f"{sample}_R{m}.fastq.gz" for m in (1, 2))
    whitelist = work / f"{sample}-onlist.txt"

    def viralscan(name: str, a: Path, b: Path) -> list[str]:
        return [
            "viralscan",
            "-o",
            str(work / name),
            "-i",
            str(index),
            "-t",
            str(t2g),
            "-gtf",
            str(gtf),
            "-s1",
            str(a),
            "-s2",
            str(b),
            "-w",
            str(whitelist),
            "-x",
            chemistry,
            "-c",
            "2",
            "--no-visual",
            "--cell-calling",
            "none",
            "--multimap-method",
            "host-conservative",
            "--multimap-molecule-assignments",
            "--yes",
        ]

    def check(name: str) -> list[str]:
        return [
            sys.executable,
            script,
            "check",
            "--panel",
            str(panel),
            "--sample",
            sample,
            "--sample-dir",
            str(work / name / sample),
            "--reference",
            reference,
        ]

    return [
        [
            sys.executable,
            generator,
            "--generation-seed",
            str(seed),
            "--error-rate",
            "0",
            "--out",
            str(panel),
            "--viral-fasta",
            str(fasta),
            "--gene-spans",
            str(spans),
            "--onlist",
            str(onlist),
            "--chemistries",
            chemistry,
            "--datasets",
            "synthetic_factorial",
            "synthetic_mixed_host_virus",
            "--cb-umi-collisions",
            "0",
        ],
        [
            sys.executable,
            script,
            "whitelist",
            "--panel",
            str(panel),
            "--sample",
            sample,
            "--out",
            str(whitelist),
        ],
        viralscan("original", r1, r2),
        viralscan("repeat", r1, r2),
        [
            sys.executable,
            script,
            "permute",
            "--fastq",
            str(r1),
            str(r2),
            "--out-dir",
            str(work / "permuted"),
        ],
        viralscan("permuted", p1, p2),
        check("original"),
        check("repeat"),
        check("permuted"),
        [
            sys.executable,
            script,
            "compare",
            str(work / "original" / sample),
            str(work / "repeat" / sample),
        ],
        [
            sys.executable,
            script,
            "compare",
            str(work / "original" / sample),
            str(work / "permuted" / sample),
        ],
    ]


# ── command line ─────────────────────────────────────────────────────────────


def _report(failures: list[str], payload: dict[str, Any]) -> None:
    print(json.dumps(dict(payload, failures=failures), indent=1, sort_keys=True, default=str))
    if failures:
        sys.exit(1)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("prepare")
    s.add_argument("--fasta", type=Path, required=True)
    s.add_argument("--gtf", type=Path, required=True)
    s.add_argument(
        "--virus-map", type=Path, help="TSV contig<TAB>virus_id; unlisted contigs are dropped"
    )
    s.add_argument("--out-fasta", type=Path, required=True)
    s.add_argument("--out-spans", type=Path, required=True)
    s = sub.add_parser("whitelist")
    s.add_argument("--panel", type=Path, required=True)
    s.add_argument("--sample", required=True)
    s.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("permute")
    s.add_argument("--fastq", type=Path, nargs=2, required=True, metavar=("R1", "R2"))
    s.add_argument("--out-dir", type=Path, required=True)
    s.add_argument("--chunk", type=int, default=127, help="records per gzip member")
    s = sub.add_parser("plan")
    for name in ("out", "fasta", "spans", "onlist", "index", "t2g", "gtf"):
        s.add_argument(f"--{name}", type=Path, required=True)
    s.add_argument("--sample", required=True)
    s.add_argument("--generation-seed", type=int, required=True)
    s.add_argument("--chemistry", default="10xv3")
    s = sub.add_parser("check")
    s.add_argument("--panel", type=Path, required=True)
    s.add_argument("--sample", required=True)
    s.add_argument("--sample-dir", type=Path, required=True)
    s.add_argument("--reference", choices=("panel", "real"), default="panel")
    s.add_argument("--allow-collision-loss", action="store_true")
    s.add_argument("--truth-dir", type=Path, help="default: <panel>/truth")
    s = sub.add_parser("compare")
    s.add_argument("first", type=Path)
    s.add_argument("second", type=Path)
    return p


def main(argv: list[str] | None = None) -> None:
    a = build_parser().parse_args(argv)
    try:
        if a.command == "prepare":
            mapping = {r[0]: r[1] for r in gt.read_rows(a.virus_map)} if a.virus_map else None
            print(json.dumps(prepare(a.fasta, a.gtf, a.out_fasta, a.out_spans, mapping), indent=1))
        elif a.command == "whitelist":
            print(write_whitelist(a.panel, a.sample, a.out))
        elif a.command == "permute":
            for source in a.fastq:
                print(permute_fastq(source, a.out_dir / source.name, a.chunk))
        elif a.command == "plan":
            for line in plan(
                out=a.out,
                fasta=a.fasta,
                spans=a.spans,
                onlist=a.onlist,
                index=a.index,
                t2g=a.t2g,
                gtf=a.gtf,
                sample=a.sample,
                seed=a.generation_seed,
                chemistry=a.chemistry,
            ):
                print(shlex.join(line))
        elif a.command == "check":
            result = check_run(
                a.panel,
                a.sample,
                a.sample_dir,
                reference=a.reference,
                allow_collision_loss=a.allow_collision_loss,
                truth_dir=a.truth_dir,
            )
            _report(result.pop("failures"), result)
        else:
            _report(compare_runs(a.first, a.second), dict(first=str(a.first), second=str(a.second)))
    except (ValueError, OSError, KeyError) as exc:
        sys.exit(f"val08_golden: {exc}")


if __name__ == "__main__":
    main()

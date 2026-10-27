"""kallisto translated search ("arm P"): protein index build and a post-hoc per-run summary.

Ported from the reproducibility repo's ``analysis/translated_search`` scripts (PLAN AA-08).
Method: Luebbert et al., Nature Biotechnology 2025 (doi:10.1038/s41587-025-02614-y); PalmDB:
Edgar et al., Nature 2022 (doi:10.1038/s41586-021-04332-2). The arm is exploratory: it never
changes ``detected``, tiers or any primary output, and writes only under ``posthoc/translated/``.

Host masking is a plain nucleotide D-list (AA-05 decision). kallisto below 0.51.1 does not
translate the D-list, so the mask was inactive (Q155); the index builder refuses such a binary.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from viralscan.run_safety import sha256_file

log = logging.getLogger(__name__)

MIN_KALLISTO = (0, 51, 1)
PALMDB_RAW = "https://raw.githubusercontent.com/pachterlab/LSCHWCP_2023/main/PalmDB/"
#: sha256 of the files the tutorial downloads (fetched 2026-10-14).
PALMDB_SHA256 = {
    "palmdb_rdrp_seqs.fa": "6a774cad43022f801e5dc1015a66344205d0d4c125c487a762706ea9687dcb09",
    "palmdb_clustered_t2g.txt": "3c4d468dcaacdbbd34bfa3ad502f9aac06def0d41f2b5c4a4f603622d68d5248",
    "ID_to_taxonomy_mapping.csv": "a9366dc0a4c271b59cef08643f2158eaa1fb296b250d9c689b6bcb20018e4ee8",
}
RANKS = ("phylum", "class", "order", "family", "genus", "species")
INDEX_MANIFEST = "translated_index.json"
RECEIPT_JSON = "translated_receipt.json"
COMPLETE_JSON = "translated_complete.json"
BARCODE_ABORT_BELOW = 0.9

_SUFFIX = re.compile(r"-\d+$")
_COMP = str.maketrans("ACGTN", "TGCAN")


class TranslatedError(RuntimeError):
    """A user-facing failure of the translated-search verbs."""


# ── barcodes (check_barcodes.py) ─────────────────────────────────────────────


def normalize_barcode(barcode: str) -> str:
    return _SUFFIX.sub("", barcode.strip().upper())


def read_barcodes(path: Path) -> list[str]:
    if str(path).endswith(".h5ad"):
        import anndata as ad

        adata = ad.read_h5ad(path, backed="r")
        try:
            return [str(b) for b in adata.obs_names]
        finally:
            adata.file.close()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return [
            line.split("\t")[0].strip()
            for line in handle
            if line.strip() and not line.startswith("#")
        ]


def _shape(barcodes: list[str]) -> dict[str, Any]:
    return {
        "n": len(barcodes),
        "with_suffix": sum(bool(_SUFFIX.search(b)) for b in barcodes),
        "lowercase": sum(b != b.upper() for b in barcodes),
        "lengths": dict(Counter(len(normalize_barcode(b)) for b in barcodes).most_common(5)),
    }


def compare(observed: list[str], frozen: list[str]) -> dict[str, Any]:
    """Overlap of two barcode lists before/after normalising; ``format_differs`` flags a mismatch."""
    obs, fro = {normalize_barcode(b) for b in observed}, {normalize_barcode(b) for b in frozen}
    raw = len(set(observed) & set(frozen))
    both = len(obs & fro)
    rc = len({b.translate(_COMP)[::-1] for b in obs} & fro)
    return {
        "observed": _shape(observed),
        "frozen": _shape(frozen),
        "raw_overlap": raw,
        "normalized_overlap": both,
        "frozen_found_fraction": both / len(fro) if fro else 0.0,
        "observed_in_frozen_fraction": both / len(obs) if obs else 0.0,
        "revcomp_overlap": rc,
        "format_differs": raw != both,
    }


# ── summary (summarize_translated.py) ────────────────────────────────────────


def taxon_ids(spec: str) -> list[str]:
    return [
        t if t.startswith("taxon_") else f"taxon_{t}"
        for t in spec.replace(",", ";").split(";")
        if t.strip()
    ]


def read_proteins(path: Path) -> dict[str, tuple[str, str, str, str]]:
    with open(path, newline="") as handle:
        return {
            r["target"]: (r["virus_id"], r["species"], r["product"], r["gene"])
            for r in csv.DictReader(handle, delimiter="\t")
        }


def summarize(
    adata: Any,
    transcripts: list[str],
    ec_lines: Iterable[str],
    bus_rows: Iterable[list[str]],
    proteins: dict[str, tuple[str, str, str, str]],
    expected: list[str],
    cells: set[str] | None = None,
    top: int = 10,
) -> dict[str, Any]:
    """``bus_rows``: ``bustools text`` rows (barcode, umi, ec, count); ``cells``: normalised barcodes."""
    if cells is not None:
        adata = adata[[normalize_barcode(b) in cells for b in adata.obs_names]]
    ec = {
        int(a): [transcripts[int(i)] for i in b.split(",")]
        for a, b in (line.rstrip("\n").split("\t")[:2] for line in ec_lines)
    }
    X = adata.X.tocsc() if hasattr(adata.X, "tocsc") else adata.X
    per_virus = {}
    for i, v in enumerate(adata.var_names):
        col = X[:, i]
        total = int(col.sum())
        if total > 0:
            per_virus[v] = (total, int((col > 0).sum()), int((col >= 10).sum()))
    prot_umi: dict[str, Counter[str]] = defaultdict(Counter)
    prot_cells: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    ambiguous: Counter[tuple[str, ...]] = Counter()
    n_records = 0
    for barcode, _umi, e, _count in (r[:4] for r in bus_rows):
        if cells is not None and normalize_barcode(barcode) not in cells:
            continue
        n_records += 1
        targets = ec[int(e)]
        viruses = {proteins[t][0] for t in targets}
        if len(viruses) == 1:
            v = next(iter(viruses))
            for t in targets:
                prot_umi[v][t] += 1
                prot_cells[v][t].add(barcode)
        else:
            ambiguous[tuple(sorted(viruses))] += 1
    label = {v: next((p[1] for p in proteins.values() if p[0] == v), "-") for v in per_virus}

    def top_proteins(v: str, n: int = 4) -> list[list[Any]]:
        return [
            [t.split("__")[-1], proteins[t][2] or proteins[t][3], c, len(prot_cells[v][t])]
            for t, c in prot_umi[v].most_common(n)
        ]

    return {
        "restricted_to_cells": cells is not None,
        "barcodes_in_matrix": int(adata.n_obs),
        "expected": {
            v: {
                "label": label.get(v, "-"),
                "umis": per_virus.get(v, (0, 0, 0))[0],
                "barcodes": per_virus.get(v, (0, 0, 0))[1],
                "barcodes_ge10_umis": per_virus.get(v, (0, 0, 0))[2],
                "distinct_proteins": len(prot_umi[v]),
                "top_proteins": top_proteins(v),
            }
            for v in expected
        },
        "total_virus_umis": int(X.sum()),
        "barcodes_any_virus": int((X.sum(axis=1) > 0).sum()),
        "viruses_nonzero": len(per_virus),
        "bus_records": n_records,
        "ambiguous_multi_virus_records": sum(ambiguous.values()),
        "top_ambiguous": [[list(k), n] for k, n in ambiguous.most_common(5)],
        "top_viruses": [
            [v, label[v], *c, top_proteins(v, 2)]
            for v, c in sorted(per_virus.items(), key=lambda kv: -kv[1][0])[:top]
        ],
    }


def bus_text_rows(bus: Path) -> Iterable[list[str]]:
    """Stream ``bustools text`` rows with the bustools that ``kb`` runs (REL-16)."""
    from viralscan.validation import tool_path

    bustools = tool_path("bustools")
    if bustools is None:
        raise TranslatedError("bustools (via kb) is not available.")
    with subprocess.Popen(
        [bustools, "text", "-p", str(bus)], stdout=subprocess.PIPE, text=True
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            yield line.rstrip("\n").split("\t")
    if proc.returncode:
        raise TranslatedError(f"bustools text failed for {bus} (exit {proc.returncode})")


def summarize_counts_dir(
    counts_dir: Path,
    proteins_tsv: Path,
    expected: list[str],
    cells: set[str] | None = None,
    top: int = 10,
    bus_rows: Iterable[list[str]] | None = None,
) -> dict[str, Any]:
    """``summarize`` over one ``kb count --aa`` output directory."""
    import anndata as ad

    d = Path(counts_dir)
    for rel in (
        "counts_unfiltered/adata.h5ad",
        "transcripts.txt",
        "matrix.ec",
        "output.unfiltered.bus",
    ):
        if not (d / rel).is_file():
            raise TranslatedError(f"{d} is not a kb count --aa output: {rel} missing.")
    with (d / "matrix.ec").open() as ec_lines:
        return summarize(
            ad.read_h5ad(d / "counts_unfiltered" / "adata.h5ad"),
            (d / "transcripts.txt").read_text().split(),
            ec_lines,
            bus_rows if bus_rows is not None else bus_text_rows(d / "output.unfiltered.bus"),
            read_proteins(proteins_tsv),
            expected,
            cells,
            top,
        )


# ── reference (combine_palmdb_reference.py) ──────────────────────────────────


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in seqs:
                    raise TranslatedError(f"Duplicate sequence ID {name} in {path.name}")
                seqs[name] = ""
            elif name is None:
                if line.strip():
                    raise TranslatedError(f"Sequence before header in {path.name}")
            else:
                seqs[name] += line.strip()
    if any(not s for s in seqs.values()):
        raise TranslatedError(f"Empty sequence in {path.name}")
    return seqs


def _read_tsv(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _label(row: dict[str, str], rep: str) -> str:
    """Third t2g column: lowest classified rank ('.' means unclassified), else the gene id."""
    return next((row[k] for k in reversed(RANKS) if row[k] not in ("", ".")), f"palmdb_{rep}")


def fetch_palmdb(dest: Path) -> dict[str, Path]:
    """Download the three pinned PalmDB tutorial files into *dest* (skips files already verified)."""
    from viralscan.scripts.build_reference import _download

    dest.mkdir(parents=True, exist_ok=True)
    out = {}
    for name, digest in PALMDB_SHA256.items():
        path = dest / name
        if not (path.is_file() and sha256_file(path) == digest):
            _download(PALMDB_RAW + name, path)
            if sha256_file(path) != digest:
                raise TranslatedError(f"{name} does not match the pinned tutorial checksum")
        out[name] = path
    return out


def combine_palmdb(annotated_dir: Path, palmdb_dir: Path, output: Path) -> dict[str, Any]:
    """Annotated proteins + PalmDB RdRP in one FASTA/t2g/proteins.tsv under *output*; returns the receipt.

    Gene names: annotated proteins keep ``taxon_<taxid>``; PalmDB sequences map to ``palmdb_<rep_ID>``.
    No sequence is altered. A PalmDB identifier colliding with an annotated one gets the prefix
    ``palmdb_``; any other duplicate is fatal.
    """
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise TranslatedError(f"Output already exists; use a new directory: {output}")
    palm = {n: palmdb_dir / n for n in PALMDB_SHA256}
    ann_seqs = read_fasta(annotated_dir / "proteins.faa")
    ann_t2g = [
        line.split("\t") for line in (annotated_dir / "viruses.t2g").read_text().splitlines()
    ]
    if {t[0] for t in ann_t2g} != set(ann_seqs) or len(ann_t2g) != len(ann_seqs):
        raise TranslatedError("Annotated proteins.faa and viruses.t2g do not list the same targets")
    pal_seqs = read_fasta(palm["palmdb_rdrp_seqs.fa"])
    pal_t2g = [
        line.split("\t") for line in palm["palmdb_clustered_t2g.txt"].read_text().splitlines()
    ]
    pal_rep = {t[0]: t[1] for t in pal_t2g}
    if len(pal_rep) != len(pal_t2g) or set(pal_rep) != set(pal_seqs):
        raise TranslatedError("PalmDB FASTA and clustered t2g do not list the same unique targets")
    tax_rows = _read_tsv(palm["ID_to_taxonomy_mapping.csv"], ",")
    if any(pal_rep.get(r["ID"]) != r["rep_ID"] for r in tax_rows) or len(tax_rows) != len(pal_rep):
        raise TranslatedError("ID_to_taxonomy_mapping.csv disagrees with palmdb_clustered_t2g.txt")
    pal_id = {i: (f"palmdb_{i}" if i in ann_seqs else i) for i in pal_seqs}
    final_ids = [*ann_seqs, *pal_id.values()]
    if len(set(final_ids)) != len(final_ids):
        raise TranslatedError("Sequence IDs are still not unique after prefixing")
    tax: dict[str, dict[str, str]] = {}
    sizes: dict[str, int] = defaultdict(int)
    for r in tax_rows:
        sizes[r["rep_ID"]] += 1
        if r["ID"] == r["rep_ID"]:
            tax[r["rep_ID"]] = r
    if set(tax) != set(sizes):
        raise TranslatedError("A PalmDB cluster representative is missing from its own mapping")
    for r in tax_rows:
        if any(r[k] != tax[r["rep_ID"]][k] for k in RANKS):
            raise TranslatedError(
                f"{r['ID']} taxonomy differs from its representative {r['rep_ID']}"
            )
    gene_rows = [
        {
            "gene": f"palmdb_{rep}",
            "rep_ID": rep,
            "n_sequences": sizes[rep],
            **{k: tax[rep][k] for k in RANKS},
            "strandedness": tax[rep]["strandedness"],
        }
        for rep in sorted(tax, key=lambda x: (len(x), x))
    ]
    overlap = {
        "annotated": {"sequences": len(ann_seqs), "virus_genes": len({t[1] for t in ann_t2g})},
        "palmdb": {
            "sequences": len(pal_seqs),
            "clusters": len(tax),
            "clusters_without_species": sum(1 for r in tax.values() if r["species"] == "."),
            "id_collisions_prefixed": sum(1 for i in pal_seqs if pal_id[i] != i),
        },
        "annotated_proteins_identical_to_a_palmdb_sequence": sum(
            1 for s in ann_seqs.values() if s in set(pal_seqs.values())
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    try:
        with (stage / "combined_proteins.faa").open("w") as handle:
            for i, s in ann_seqs.items():
                handle.write(f">{i}\n{s}\n")
            for i, s in pal_seqs.items():
                handle.write(f">{pal_id[i]}\n{s}\n")
        with (stage / "combined.t2g").open("w") as handle:
            for t in ann_t2g:
                handle.write("\t".join(t) + "\n")
            for i in pal_seqs:
                rep = pal_rep[i]
                handle.write(f"{pal_id[i]}\tpalmdb_{rep}\t{_label(tax[rep], rep)}\n")
        _write_tsv(stage / "palmdb_taxonomy.tsv", gene_rows)
        ann_prot = _read_tsv(annotated_dir / "proteins.tsv")
        ann_cols = list(ann_prot[0])
        _write_tsv(
            stage / "proteins.tsv",
            ann_prot
            + [
                {
                    **dict.fromkeys(ann_cols, ""),
                    "target": pal_id[i],
                    "protein_id": i,
                    "virus_id": f"palmdb_{pal_rep[i]}",
                    "species": _label(tax[pal_rep[i]], pal_rep[i]),
                    "family": tax[pal_rep[i]]["family"],
                    "length_aa": len(pal_seqs[i]),
                }
                for i in pal_seqs
            ],
        )
        receipt: dict[str, Any] = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "palmdb_files": {
                n: {"url": PALMDB_RAW + n, "sha256": sha256_file(p)} for n, p in palm.items()
            },
            "annotated_sha256": {
                n: sha256_file(annotated_dir / n) for n in ("proteins.faa", "viruses.t2g")
            },
            "overlap": overlap,
        }
        receipt["output_sha256"] = {p.name: sha256_file(p) for p in sorted(stage.iterdir())}
        (stage / "reference.json").write_text(json.dumps(receipt, indent=2) + "\n")
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return receipt


# ── index and count (kb --aa) ────────────────────────────────────────────────


def _version_tuple(text: str) -> tuple[int, ...]:
    found = re.search(r"(\d+)\.(\d+)\.(\d+)", text)
    return tuple(int(x) for x in found.groups()) if found else ()


def kallisto_translated_ok() -> tuple[bool, str]:
    """``(ok, version)``: the kallisto that ``kb`` runs translates the D-list (>= 0.51.1)."""
    from viralscan.validation import kb_tools

    tools = kb_tools()
    if "kallisto" not in tools:
        return False, ""
    version = tools["kallisto"][0]
    return _version_tuple(version) >= MIN_KALLISTO, version


def _kb() -> str:
    kb = shutil.which("kb")
    if kb is None:
        raise TranslatedError("kb (kb-python) is not on PATH.")
    return kb


def build_translated_index(
    proteins_faa: Path, host_mask_fa: Path, index: Path, threads: int = 4
) -> list[str]:
    """``kb ref --aa`` with a plain nucleotide D-list; returns the command run. Fails closed on old kallisto."""
    from viralscan.validation import kb_tools

    ok, version = kallisto_translated_ok()
    if not ok:
        raise TranslatedError(
            f"kallisto {version or 'missing'} does not translate the D-list; "
            f">= {'.'.join(map(str, MIN_KALLISTO))} is required (Q155)."
        )
    kallisto = kb_tools()["kallisto"][1]
    cmd = [
        _kb(), "ref", "--kallisto", kallisto, "--aa", "-k", "31", "--workflow", "custom",
        "--d-list", str(host_mask_fa), "-t", str(threads), "-i", str(index), str(proteins_faa),
    ]  # fmt: skip
    subprocess.run(cmd, check=True)
    return cmd


def count_translated(
    r1: str,
    r2: str,
    technology: str,
    whitelist: str | None,
    index: Path,
    t2g: Path,
    out: Path,
    threads: int,
) -> list[str]:
    """``kb count --aa --parity single`` (six-frame search, so no strand flag); returns the command."""
    from viralscan.chemistry import kb_whitelist_arg

    wl = kb_whitelist_arg(technology, whitelist)
    cmd = [
        _kb(), "count", "--aa", "--parity", "single", "-x", technology, "-t", str(threads),
        *(["-w", wl] if wl else []), "--num", "--tmp", str(out.with_name(out.name + "_tmp")),
        "-i", str(index), "-g", str(t2g), "--h5ad", "-o", str(out), r1, r2,
    ]  # fmt: skip
    subprocess.run(cmd, check=True)
    return cmd


# ── verbs ────────────────────────────────────────────────────────────────────


def run_translated_index(args: argparse.Namespace) -> Path:
    """``translated-index``: annotated proteins (+ PalmDB) -> ``kb ref --aa`` index with a manifest."""
    out = Path(args.out_dir).resolve()
    index_dir = out / "index"
    if index_dir.exists() and any(index_dir.iterdir()):
        raise TranslatedError(f"{index_dir} is not empty; use a new --out-dir.")
    proteins_dir = Path(args.proteins_dir).resolve()
    reference = proteins_dir
    if args.palmdb:
        palm = fetch_palmdb(out / "palmdb")
        combine_palmdb(proteins_dir, palm["palmdb_rdrp_seqs.fa"].parent, out / "reference")
        reference = out / "reference"
        fasta, t2g = reference / "combined_proteins.faa", reference / "combined.t2g"
    else:
        fasta, t2g = reference / "proteins.faa", reference / "viruses.t2g"
    for p in (fasta, t2g, reference / "proteins.tsv"):
        if not p.is_file():
            raise TranslatedError(f"{p} is missing.")
    index_dir.mkdir(parents=True, exist_ok=True)
    index = index_dir / "index.idx"
    cmd = build_translated_index(fasta, Path(args.host_mask).resolve(), index, args.threads)
    shutil.copy(t2g, index_dir / "t2g.txt")
    shutil.copy(reference / "proteins.tsv", index_dir / "proteins.tsv")
    manifest = {
        "schema_version": 1,
        "command": cmd,
        "kallisto_version": kallisto_translated_ok()[1],
        "host_mask": {
            "path": str(Path(args.host_mask).resolve()),
            "sha256": sha256_file(args.host_mask),
        },
        "palmdb": bool(args.palmdb),
        "inputs_sha256": {"proteins": sha256_file(fasta), "t2g": sha256_file(t2g)},
        "index_sha256": sha256_file(index),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    (index_dir / INDEX_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")  # last
    return index_dir


def _sample_config(run_dir: Path, sample: str | None) -> tuple[Path, Any]:
    from viralscan.discovery import resolve_sample_dir
    from viralscan.runconfig import RunConfig

    sample_dir = resolve_sample_dir(run_dir, sample)
    return sample_dir, RunConfig.from_yaml(str(sample_dir / "config.yaml"))


def run_translated(args: argparse.Namespace) -> Path:
    """``translated``: count (or reuse ``--counts-dir``) and summarise into ``posthoc/translated``."""
    from viralscan.outputs import posthoc_dir

    run_dir = Path(args.run_dir).resolve()
    index_dir = Path(args.index_dir).resolve()
    if not (index_dir / INDEX_MANIFEST).is_file():
        raise TranslatedError(
            f"{index_dir} is not a finished translated-index ({INDEX_MANIFEST} missing)."
        )
    sample_dir, config = _sample_config(run_dir, getattr(args, "sample", None))
    output = (
        Path(args.output).resolve()
        if getattr(args, "output", None)
        else posthoc_dir(run_dir, "translated", sample_dir.name)
    )
    if output.exists() and any(output.iterdir()):
        raise TranslatedError(f"{output} is not empty.")
    output.mkdir(parents=True, exist_ok=True)
    counts = Path(args.counts_dir).resolve() if getattr(args, "counts_dir", None) else None
    cmd: list[str] | None = None
    if counts is None:
        counts = output / "counts"
        cmd = count_translated(
            config.sample1, config.sample2, config.technology, config.whitelist,
            index_dir / "index.idx", index_dir / "t2g.txt", counts, args.threads,
        )  # fmt: skip
    cells_file = (
        Path(args.cells).resolve()
        if getattr(args, "cells", None)
        else sample_dir / "results" / "called_cells.tsv"
    )
    observed = read_barcodes(counts / "counts_unfiltered" / "adata.h5ad")
    frozen = [b for b in read_barcodes(cells_file) if b.lower() != "barcode"]  # header row
    check = compare(observed, frozen)
    (output / "barcode_check.json").write_text(json.dumps(check, indent=1) + "\n")
    if check["frozen_found_fraction"] < BARCODE_ABORT_BELOW:
        raise TranslatedError(
            f"only {check['frozen_found_fraction']:.1%} of the called cells appear in the translated "
            f"matrix (threshold {BARCODE_ABORT_BELOW:.0%}); barcode formats differ (see barcode_check.json)."
        )
    expected = taxon_ids(args.expected_taxa) if getattr(args, "expected_taxa", None) else []
    proteins_tsv = index_dir / "proteins.tsv"
    summary_all = summarize_counts_dir(counts, proteins_tsv, expected)
    summary_cells = summarize_counts_dir(
        counts, proteins_tsv, expected, {normalize_barcode(b) for b in frozen}
    )
    (output / "summary.all.json").write_text(json.dumps(summary_all, indent=1) + "\n")
    (output / "summary.cells.json").write_text(json.dumps(summary_cells, indent=1) + "\n")
    receipt = {
        "schema_version": 1,
        "sample": sample_dir.name,
        "command": cmd,
        "counts_dir": str(counts),
        "index_manifest": json.loads((index_dir / INDEX_MANIFEST).read_text()),
        "cells_file": str(cells_file),
        "kallisto_translated": dict(zip(("ok", "version"), kallisto_translated_ok())),
        "exploratory": True,
    }
    (output / RECEIPT_JSON).write_text(json.dumps(receipt, indent=2) + "\n")
    (output / COMPLETE_JSON).write_text(
        json.dumps({"outputs": sorted(p.name for p in output.iterdir())}, indent=2) + "\n"
    )  # last
    log.info("Translated search summary written to %s", output)
    return output

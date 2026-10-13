"""VAL-08 (tiny, local): real kallisto/bustools recover the generator's planted molecules exactly.

The golden tiny panel is generated with ``--error-rate 0`` into ``tmp_path`` (nothing is committed),
indexed from the generator's own ``templates/transcripts.fasta``, pseudoaligned with ``kallisto bus``
under ``-x 10xv3`` and read back with ``bustools text``. A (barcode, UMI) is a viral call when every
transcript of its equivalence class belongs to one virus. Those calls must equal the planted
molecules of ``truth_manifest.tsv`` exactly (``v3_score.precision_recall_f1``), no host molecule may
be called viral, and a second generation must reproduce the same bytes.

Planting is keyed at (barcode, UMI, virus): a same-key collision pair is one item there, as the
protocol's collision_rule (one corrected key contributes at most one truth molecule) requires.

The full-workflow arm also runs ``kb count`` and ViralScan with a synthetic on-list and
``--cell-calling none``. It retains default collision-challenge failures in recall's denominator,
then checks exact recovery on a collision-free golden generated with the existing parameter.
Reference/FASTQ order and concatenated gzip-member transport chunks must preserve that result.
These tests establish software bookkeeping only; they do not qualify scientific outcomes.

Run with ``PYTHONPATH=src pytest -m integration tests/integration/test_truth_panel_exact.py``.
"""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import v3_score as score  # noqa: E402
import val08_golden as val08  # noqa: E402

from viralscan.evidence import have_tools  # noqa: E402
from viralscan.run_safety import sha256_file  # noqa: E402
from viralscan.validation import tool_path  # noqa: E402
from viralscan.virus_identity import VirusIdentityTable  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _panel_args():
    return [
        "--generation-seed",
        "1",
        "--error-rate",
        "0",
        "--datasets",
        "synthetic_factorial",
        "synthetic_mixed_host_virus",
        "--chemistries",
        "10xv3",
    ]


def _full_reference(panel: Path, work: Path, *, reverse: bool = False):
    """Use the generator's gene models and genome IDs; host cDNA names its transcript."""
    work.mkdir()
    models = read_tsv(panel / "templates" / "viral_gene_models.tsv")
    gene_virus = {row["gene_id"]: row["virus_id"] for row in models}
    records = gt.read_fasta_records(panel / "templates" / "transcripts.fasta")
    names = sorted(records, reverse=reverse)
    fasta = work / "transcripts.fasta"
    fasta.write_text("".join(f">{name}\n{records[name]}\n" for name in names))
    t2g = work / "t2g.tsv"
    t2g.write_text(
        "".join(
            f"{name}\t{name if name in gene_virus else 'host-' + name}\t{name}\t"
            f"{'viral' if name in gene_virus else 'host'}\t{gene_virus.get(name, name)}\n"
            for name in names
        )
    )
    gtf = work / "viral.gtf"
    gtf.write_text(
        "".join(
            f"{r['virus_id']}\tsynthetic\texon\t{int(r['start']) + 1}\t{r['end']}\t.\t+\t.\t"
            f'gene_id "{r["gene_id"]}"; transcript_id "{r["gene_id"]}";\n'
            for r in models
        )
    )
    index = work / "index.idx"
    subprocess.run(
        [tool_path("kallisto"), "index", "-i", str(index), str(fasta)],
        check=True,
        capture_output=True,
        timeout=600,
    )
    return index, t2g, gtf, gene_virus


def _full_run(panel: Path, sid: str, out: Path, reference, *, reorder: bool = False):
    index, t2g, gtf, gene_virus = reference
    pool = read_tsv(panel / sid / "barcode_pool.tsv")
    whitelist = out.parent / f"{out.name}-onlist.txt"
    whitelist.write_text("".join(row["cell_barcode"] + "\n" for row in pool))
    fastqs = [panel / sid / f"{sid}_R{mate}.fastq.gz" for mate in (1, 2)]
    if reorder:
        directory = out.parent / "reordered FASTQs"
        directory.mkdir()
        changed = []
        for original in fastqs:
            with gzip.open(original, "rt") as handle:
                lines = handle.readlines()
            records = ["".join(lines[i : i + 4]) for i in range(0, len(lines), 4)][::-1]
            path = directory / original.name
            # Multiple gzip members preserve the same paired records across transport chunks.
            path.write_bytes(
                b"".join(
                    gzip.compress("".join(records[i : i + 127]).encode(), mtime=0)
                    for i in range(0, len(records), 127)
                )
            )
            changed.append(path)
        fastqs = changed
    command = [
        sys.executable,
        "-m",
        "viralscan.menu",
        "-o",
        str(out),
        "-i",
        str(index),
        "-t",
        str(t2g),
        "-gtf",
        str(gtf),
        "-s1",
        str(fastqs[0]),
        "-s2",
        str(fastqs[1]),
        "-w",
        str(whitelist),
        "-x",
        "10xv3",
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
    completed = subprocess.run(
        command,
        cwd=Path(__file__).resolve().parents[2],
        capture_output=True,
        text=True,
        timeout=1800,
    )
    (out.parent / f"{out.name}-cli.log").write_text(completed.stdout + completed.stderr)
    assert completed.returncode == 0, completed.stdout[-6000:] + completed.stderr[-6000:]
    sample = out / sid
    identity = VirusIdentityTable.read_tsv(sample / "results" / "virus_identity.tsv")
    actual_identity = {g.gene_id: g.virus_name for g in identity.genes if g.viral}
    assert actual_identity == gene_virus  # genes are pooled by their actual synthetic genome IDs
    assert all(not g.viral for g in identity.genes if g.gene_id.startswith("host-"))
    with gzip.open(sample / "results" / "molecule_assignments.tsv.gz", "rt") as handle:
        assignments = list(csv.DictReader(handle, delimiter="\t"))
    return sample, assignments


def _assigned_sets(assignments, gene_virus):
    viral, host, unresolved = set(), set(), set()
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


def _native_host_anchor(sample, panel, sid, gene_virus, host_truth):
    """The host anchor is declared from inputs, with its independently planted denominator."""
    import numpy as np

    barcodes = [r["cell_barcode"] for r in read_tsv(panel / sid / "barcode_pool.tsv")]
    features = sorted(
        "host-" + transcript
        for transcript in gt.read_fasta_records(panel / "templates" / "transcripts.fasta")
        if transcript not in gene_virus
    )
    manifest = json.loads((sample.parent / "run_manifest.json").read_text())
    counts = sample / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad"
    matrix, provenance = score.read_run_count_layer(
        sample,
        barcodes=barcodes,
        features=features,
        layer="counts_unique",
        denominator=len(host_truth),
        expected_run_fingerprint=manifest["run_fingerprint"],
        expected_counts_sha256=sha256_file(counts),
    )
    expected_per_cell = Counter(cb for cb, _ in host_truth)
    np.testing.assert_array_equal(
        np.asarray(matrix.sum(axis=1)).ravel(), [expected_per_cell[cb] for cb in barcodes]
    )
    assert provenance["source_layer"] == "counts_unique"
    assert provenance["unit"] == "integer-unique-molecule"
    return provenance


@pytest.mark.integration
def test_full_viralscan_truth_and_conservation(tmp_path):
    import anndata as ad
    import numpy as np

    missing = have_tools(["kallisto", "bustools", "kb", "snakemake"])
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")
    work = tmp_path / "golden with spaces"
    work.mkdir()
    receipts = {}
    previous = None
    for name, collision, reorder in [
        ("challenge", None, False),
        ("golden", "0", False),
        ("permuted", "0", True),
    ]:
        panel = work / ("challenge-panel" if collision is None else "golden-panel")
        if not panel.exists():
            args = _panel_args() + (["--cb-umi-collisions", collision] if collision else [])
            manifest = gt.generate(gt.parse_args([*args, "--out", str(panel)]))
        else:
            manifest = json.loads((panel / "run_manifest.json").read_text())
        assert not manifest["holdout_samples"]
        sid = next(
            s
            for s, detail in manifest["samples"].items()
            if detail["dataset"] == "synthetic_factorial"
        )
        reference = _full_reference(panel, work / f"{name}-reference", reverse=reorder)
        sample, assignments = _full_run(
            panel, sid, work / f"{name}-out", reference, reorder=reorder
        )
        rows = [
            r
            for r in read_tsv(panel / "truth" / "truth_manifest.tsv")
            if r["biological_sample_id"] == sid
        ]
        planted = {
            (r["cell_barcode"], r["umi"], r["planted_virus_id"])
            for r in rows
            if r["is_planted_molecule"] == "true"
        }
        host_truth = {
            (r["cell_barcode"], r["umi"]) for r in rows if r["record_granularity"] == "molecule"
        }
        viral, host, unresolved = _assigned_sets(assignments, reference[3])
        metric = score.precision_recall_f1(planted, viral)
        matrix_path = sample / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad"
        adata = ad.read_h5ad(matrix_path)
        audit = dict(adata.uns["molecule_audit"])
        status = Counter(row["status"] for row in assignments)
        known_keys = host_truth | {(cb, umi) for cb, umi, _ in planted}
        assert {(row["barcode"], row["umi"]) for row in assignments} == known_keys
        assert audit["input_molecules"] == len(known_keys)
        assert audit["input_molecules"] == len(assignments)
        assert audit["unresolved_molecules"] == len(unresolved)
        assert audit["unique_molecules"] == status["unique"]
        assert audit["ambiguous_molecules"] == status["ambiguous"]
        assert audit["resolved_molecules"] + len(unresolved) == len(assignments)
        assert float(adata.X.sum()) == pytest.approx(audit["resolved_molecules"])
        np.testing.assert_allclose(
            adata.X.toarray(),
            (adata.layers["counts_unique"] + adata.layers["counts_ambiguous_allocated"]).toarray(),
        )
        np.testing.assert_array_equal(
            adata.layers["counts_ambiguous_allocated"].toarray(),
            adata.layers["counts_multimap_host_conservative"].toarray(),
        )
        assert adata.uns["multimap_method"] == "host-conservative"
        assert adata.uns["quantification_unit"] == "bustools-resolved-cb-umi-molecule"
        assert not host_truth.intersection((cb, umi) for cb, umi, _ in viral)
        assert not host.intersection((cb, umi) for cb, umi, _ in planted)
        assert host == host_truth
        receipts[name] = {
            "metric": metric,
            "lost_planted": sorted(planted - viral),
            "unresolved_keys": sorted(unresolved),
            "n_host_truth": len(host_truth),
            "audit": {
                k: float(v)
                if isinstance(v, np.floating)
                else int(v)
                if isinstance(v, np.integer)
                else v
                for k, v in audit.items()
            },
        }
        (work / "software-receipt.json").write_text(json.dumps(receipts, indent=2))
        # the VAL-08 expected-output checker must agree with the inline assertions above
        checked = val08.check_run(panel, sid, sample, allow_collision_loss=name == "challenge")
        assert checked["failures"] == [], checked["failures"]
        if name == "golden":
            golden_sample = sample
        if name == "permuted":
            assert val08.compare_runs(golden_sample, sample) == []
        if name == "challenge":
            lost = planted - viral
            assert metric["fp"] == 0 and lost
            assert {(cb, umi) for cb, umi, _ in lost} <= unresolved
            # Keep the failures in recall's denominator; disjoint-gene collisions are unresolved.
            assert all(
                any(
                    r["collision_kind"] == "same"
                    and (r["cell_barcode"], r["umi"], r["planted_virus_id"]) == key
                    for r in rows
                )
                for key in lost
            )
        else:
            assert viral == planted, metric
            assert not unresolved
            receipts[name]["native_host_anchor"] = _native_host_anchor(
                sample, panel, sid, reference[3], host_truth
            )
            (work / "software-receipt.json").write_text(json.dumps(receipts, indent=2))
            selected = adata[sorted(adata.obs_names), sorted(adata.var_names)].X.toarray()
            canonical = sorted(
                (
                    row["barcode"],
                    row["umi"],
                    row["status"],
                    tuple(sorted(zip(row["genes"].split(";"), row["weights"].split(";")))),
                )
                for row in assignments
            )
            if previous is not None:
                np.testing.assert_array_equal(selected, previous[0])
                assert canonical == previous[1]
            previous = selected, canonical


@pytest.mark.integration
def test_kallisto_recovers_exactly_the_planted_molecules(tmp_path: Path) -> None:
    missing = have_tools(["kallisto", "bustools"])
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")
    args = [
        "--generation-seed",
        "1",
        "--error-rate",
        "0",
        "--datasets",
        "synthetic_factorial",
        "synthetic_mixed_host_virus",
        "--chemistries",
        "10xv3",
    ]
    panel = tmp_path / "panel"
    manifest = gt.generate(gt.parse_args([*args, "--out", str(panel)]))
    again = gt.generate(gt.parse_args([*args, "--out", str(tmp_path / "again")]))
    assert manifest["samples"] == again["samples"]  # determinism, including the FASTQ hashes
    sid = next(s for s, d in manifest["samples"].items() if d["dataset"] == "synthetic_factorial")

    index = tmp_path / "index.idx"
    subprocess.run(
        [
            tool_path("kallisto"),
            "index",
            "-i",
            str(index),
            str(panel / "templates" / "transcripts.fasta"),
        ],
        check=True,
        capture_output=True,
    )
    bus = tmp_path / "bus"
    subprocess.run(
        [
            tool_path("kallisto"),
            "bus",
            "-i",
            str(index),
            "-o",
            str(bus),
            "-x",
            "10xv3",
            "-t",
            "1",
            str(panel / sid / f"{sid}_R1.fastq.gz"),
            str(panel / sid / f"{sid}_R2.fastq.gz"),
        ],
        check=True,
        capture_output=True,
    )
    text = subprocess.run(
        [tool_path("bustools"), "text", "-p", str(bus / "output.bus")],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    names = (bus / "transcripts.txt").read_text().split()
    classes = {
        int(ec): [names[int(i)] for i in members.split(",")]
        for ec, members in (line.split() for line in (bus / "matrix.ec").read_text().splitlines())
    }
    viral: set[tuple[str, str, str]] = set()
    host: set[tuple[str, str]] = set()
    for line in text.splitlines():
        cb, umi, ec = line.split()[:3]
        members = classes[int(ec)]
        viruses = {m.split("_G")[0] for m in members if m.split("_G")[0] in gt.TARGETS}
        if len(viruses) == 1 and all(m.split("_G")[0] in gt.TARGETS for m in members):
            viral.add((cb, umi, next(iter(viruses))))
        elif not viruses:
            host.add((cb, umi))

    rows = [
        r
        for r in read_tsv(panel / "truth" / "truth_manifest.tsv")
        if r["biological_sample_id"] == sid
    ]
    truth = [r for r in rows if r["is_planted_molecule"] == "true"]
    assert any(
        r["origin"] == "hopped" for r in truth
    )  # the golden panel includes a hopped molecule
    planted = {(r["cell_barcode"], r["umi"], r["planted_virus_id"]) for r in truth}
    assert planted
    result = score.precision_recall_f1(planted, viral)
    assert (result["tp"], result["fp"], result["fn"]) == (len(planted), 0, 0), result
    host_truth = {
        (r["cell_barcode"], r["umi"]) for r in rows if r["record_granularity"] == "molecule"
    }
    allocation = score.host_virus_allocation_error(
        host_truth, {(cb, umi) for cb, umi, _ in planted}, {(cb, umi) for cb, umi, _ in viral}, host
    )
    assert allocation["host_to_viral"] == 0.0

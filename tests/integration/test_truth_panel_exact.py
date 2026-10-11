"""VAL-08 (tiny, local): real kallisto/bustools recover the generator's planted molecules exactly.

The golden tiny panel is generated with ``--error-rate 0`` into ``tmp_path`` (nothing is committed),
indexed from the generator's own ``templates/transcripts.fasta``, pseudoaligned with ``kallisto bus``
under ``-x 10xv3`` and read back with ``bustools text``. A (barcode, UMI) is a viral call when every
transcript of its equivalence class belongs to one virus. Those calls must equal the planted
molecules of ``truth_manifest.tsv`` exactly (``v3_score.precision_recall_f1``), no host molecule may
be called viral, and a second generation must reproduce the same bytes.

Planting is keyed at (barcode, UMI, virus): a same-key collision pair is one item there, as the
protocol's collision_rule (one corrected key contributes at most one truth molecule) requires.

Scope: this covers the generator and ``kallisto bus``. It does not run ``kb count`` or ViralScan's
cell calling, so ``VAL-08`` stays open for the full-workflow arm.

Run with ``PYTHONPATH=src pytest -m integration tests/integration/test_truth_panel_exact.py``.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import v3_score as score  # noqa: E402

from viralscan.evidence import have_tools  # noqa: E402
from viralscan.validation import tool_path  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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

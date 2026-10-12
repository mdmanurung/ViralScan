"""VAL-08 golden-run tooling (scripts/val08_golden.py): command lines and expected-output checks, offline.

No ``kb count`` or ViralScan runs here. A tiny panel is generated into ``tmp_path`` and a stand-in
ViralScan sample directory is written from its truth (``fake_run``); each defect a cluster run could
show (a lost molecule, a host molecule called viral, a broken audit, an output that depends on read
order) is injected and the checker must name it. The real-tool arm is
``tests/integration/test_truth_panel_exact.py``.
"""

from __future__ import annotations

import csv
import gzip
import json
import random
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import val08_golden as val  # noqa: E402

from tests.test_generate_truth_panel import fastq_records, read_tsv  # noqa: E402

ARGS = [
    "--generation-seed",
    "1",
    "--error-rate",
    "0",
    "--datasets",
    "synthetic_factorial",
    "--chemistries",
    "10xv3",
]


@pytest.fixture(scope="module")
def golden(tmp_path_factory):
    out = tmp_path_factory.mktemp("golden")
    manifest = gt.generate(gt.parse_args([*ARGS, "--cb-umi-collisions", "0", "--out", str(out)]))
    return out, next(iter(manifest["samples"]))


@pytest.fixture(scope="module")
def challenge(tmp_path_factory):
    out = tmp_path_factory.mktemp("challenge")
    manifest = gt.generate(gt.parse_args([*ARGS, "--out", str(out)]))
    return out, next(iter(manifest["samples"]))


def fake_run(
    panel: Path,
    sample: str,
    dest: Path,
    *,
    drop_planted: int = 0,
    unresolve_same_collisions: bool = False,
    shuffle: int | None = None,
    include_host: bool = True,
    extra_viral: bool = False,
    host_as_viral: bool = False,
    audit_shift: int = 0,
    swap_gene: bool = False,
) -> Path:
    """Write ``<dest>/<sample>/`` as ViralScan would for a perfect (or deliberately broken) run."""
    import anndata as ad
    import numpy as np
    from scipy import sparse

    rows = [
        r
        for r in read_tsv(panel / "truth" / "truth_manifest.tsv")
        if r["biological_sample_id"] == sample
    ]
    pool = [r["cell_barcode"] for r in read_tsv(panel / sample / "barcode_pool.tsv")]
    planted: dict[tuple[str, str, str], tuple[str, str]] = {}
    for r in rows:
        if r["is_planted_molecule"] == "true":
            planted[(r["cell_barcode"], r["umi"], r["planted_virus_id"])] = (
                r["planted_gene_id"],
                r["collision_kind"],
            )
    molecules: list[tuple[str, str, str, str]] = []  # barcode, umi, status, gene ("" = unresolved)
    seen: set[tuple[str, str]] = set()
    for n, ((cb, umi, _), (gene, kind)) in enumerate(sorted(planted.items())):
        if n < drop_planted:
            continue
        if unresolve_same_collisions and kind == "same":
            if (cb, umi) not in seen:
                molecules.append((cb, umi, "unresolved", ""))
        elif (cb, umi) not in seen:
            molecules.append((cb, umi, "unique", gene))
        seen.add((cb, umi))
    viral_gene = sorted(val.gene_virus_map(panel))[0]
    if include_host:
        host = {
            (r["cell_barcode"], r["umi"]): "host-" + r["host_gene_id"]
            for r in rows
            if r["record_granularity"] == "molecule"
        }
        for n, ((cb, umi), gene) in enumerate(sorted(host.items())):
            molecules.append((cb, umi, "unique", viral_gene if host_as_viral and n == 0 else gene))
    if extra_viral:
        molecules.append((pool[0], "A" * 12, "unique", viral_gene))
    if swap_gene:
        cb, umi, status, gene = molecules[0]
        molecules[0] = (cb, umi, status, "host-swapped")
    if shuffle is not None:
        random.Random(shuffle).shuffle(molecules)
    sample_dir = dest / sample
    (sample_dir / "results").mkdir(parents=True)
    with gzip.open(
        sample_dir / "results" / "molecule_assignments.tsv.gz", "wt", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["barcode", "umi", "ec_ids", "status", "genes", "weights"])
        for cb, umi, status, gene in molecules:
            writer.writerow([cb, umi, "0", status, gene, "1" if gene else ""])
    genes = sorted({g for _, _, _, g in molecules if g})
    cells = pool if shuffle is None else random.Random(shuffle).sample(pool, len(pool))
    unique = np.zeros((len(cells), len(genes)))
    for cb, _, _status, gene in molecules:
        if gene:
            unique[cells.index(cb), genes.index(gene)] += 1
    adata = ad.AnnData(sparse.csr_matrix(unique))
    adata.obs_names, adata.var_names = cells, genes
    adata.layers["counts_unique"] = sparse.csr_matrix(unique)
    adata.layers["counts_ambiguous_allocated"] = sparse.csr_matrix(unique.shape)
    resolved = sum(1 for m in molecules if m[3])
    adata.uns["molecule_audit"] = dict(
        source="viralscan",  # a real audit carries non-numeric fields too
        input_molecules=len(molecules) + audit_shift,
        resolved_molecules=resolved,
        unique_molecules=resolved,
        ambiguous_molecules=0,
        unresolved_molecules=len(molecules) - resolved,
    )
    target = sample_dir / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad"
    target.parent.mkdir(parents=True)
    adata.write_h5ad(target)
    return sample_dir


def failures(panel, sample, tmp_path, name="run", **kwargs):
    check = {k: kwargs.pop(k) for k in ("reference", "allow_collision_loss") if k in kwargs}
    sample_dir = fake_run(panel, sample, tmp_path / name, **kwargs)
    return val.check_run(panel, sample, sample_dir, **check)["failures"]


# ── check ────────────────────────────────────────────────────────────────────


def test_a_perfect_run_has_no_failures_and_exact_metrics(golden, tmp_path):
    panel, sample = golden
    sample_dir = fake_run(panel, sample, tmp_path / "ok")
    result = val.check_run(panel, sample, sample_dir)
    assert result["failures"] == []
    assert (result["metric"]["precision"], result["metric"]["recall"]) == (1.0, 1.0)
    assert result["n_lost"] == 0 and result["metric"]["tp"] > 0


def test_a_lost_planted_molecule_is_a_failure(golden, tmp_path):
    panel, sample = golden
    found = failures(panel, sample, tmp_path, drop_planted=1)
    assert any("not recovered" in f for f in found)


def test_a_false_viral_call_is_a_failure(golden, tmp_path):
    panel, sample = golden
    assert any("not planted" in f for f in failures(panel, sample, tmp_path, extra_viral=True))


def test_a_host_molecule_called_viral_is_a_failure(golden, tmp_path):
    panel, sample = golden
    assert any("called viral" in f for f in failures(panel, sample, tmp_path, host_as_viral=True))


@pytest.mark.parametrize("shift", [1, -1])
def test_a_broken_molecule_audit_is_a_failure(golden, tmp_path, shift):
    panel, sample = golden
    found = failures(panel, sample, tmp_path, audit_shift=shift)
    assert any("conservation" in f for f in found)


def test_the_count_matrix_must_equal_the_layers(golden, tmp_path):
    import anndata as ad

    panel, sample = golden
    sample_dir = fake_run(panel, sample, tmp_path / "x")
    path = sample_dir / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad"
    adata = ad.read_h5ad(path)
    adata.X = adata.X.copy()
    adata.X[0, 0] += 1
    adata.write_h5ad(path)
    assert any(
        "X differs" in f or "sum of X" in f
        for f in val.check_run(panel, sample, sample_dir)["failures"]
    )


def test_real_reference_does_not_require_the_host_set_to_be_indexed(golden, tmp_path):
    panel, sample = golden
    assert failures(panel, sample, tmp_path, "a", include_host=False)  # panel mode needs the host
    assert failures(panel, sample, tmp_path, "b", include_host=False, reference="real") == []
    # but a host molecule called viral still fails against a real reference
    assert failures(panel, sample, tmp_path, "c", reference="real", host_as_viral=True)


def test_collision_losses_are_allowed_only_when_unresolved_same_key(challenge, tmp_path):
    panel, sample = challenge
    assert failures(panel, sample, tmp_path, "a", unresolve_same_collisions=True)
    assert (
        failures(
            panel, sample, tmp_path, "b", unresolve_same_collisions=True, allow_collision_loss=True
        )
        == []
    )
    # a missing non-collision molecule is not excused by the option
    dropped = failures(panel, sample, tmp_path, "c", drop_planted=1, allow_collision_loss=True)
    assert any("not recovered" in f for f in dropped)


def test_missing_truth_is_an_error_not_a_pass(golden, tmp_path):
    panel, _ = golden
    with pytest.raises(ValueError, match="no truth rows"):
        val.truth_of(panel, "s000000000000")


# ── compare (determinism, read-order and chunk-size invariance) ─────────────


def test_identical_results_compare_equal_whatever_the_row_and_barcode_order(golden, tmp_path):
    panel, sample = golden
    first = fake_run(panel, sample, tmp_path / "a")
    assert val.compare_runs(first, fake_run(panel, sample, tmp_path / "b", shuffle=3)) == []
    assert val.compare_runs(first, fake_run(panel, sample, tmp_path / "c", shuffle=4)) == []


def test_results_that_differ_are_named(golden, tmp_path):
    panel, sample = golden
    first = fake_run(panel, sample, tmp_path / "a")
    assert "molecule assignments differ" in val.compare_runs(
        first, fake_run(panel, sample, tmp_path / "b", drop_planted=1)
    )
    assert "molecule assignments differ" in val.compare_runs(
        first, fake_run(panel, sample, tmp_path / "c", swap_gene=True)
    )
    audit = val.compare_runs(first, fake_run(panel, sample, tmp_path / "d", audit_shift=1))
    assert audit == ["molecule audits differ"]


# ── permuted FASTQs and whitelist ────────────────────────────────────────────


@pytest.mark.parametrize("chunk", [1, 3, 127, 100_000])
def test_permuted_fastq_keeps_the_same_pairs_in_reverse_order(golden, tmp_path, chunk):
    panel, sample = golden
    pairs = {}
    for mate in (1, 2):
        source = panel / sample / f"{sample}_R{mate}.fastq.gz"
        moved = val.permute_fastq(source, tmp_path / source.name, chunk)
        pairs[mate] = (fastq_records(source), fastq_records(moved))
    for original, moved in pairs.values():
        assert moved == original[::-1] and len(original) > 3
    # mates stay aligned: the same read name at every position
    assert [n for n, _ in pairs[1][1]] == [n for n, _ in pairs[2][1]]
    members = (tmp_path / f"{sample}_R1.fastq.gz").read_bytes().count(b"\x1f\x8b\x08")
    assert (
        members >= len(pairs[1][0]) // chunk
    )  # chunk-sized members, not one stream (for small chunks)


def test_whitelist_is_the_barcode_pool(golden, tmp_path):
    panel, sample = golden
    n = val.write_whitelist(panel, sample, tmp_path / "w.txt")
    lines = (tmp_path / "w.txt").read_text().split()
    assert (
        n
        == len(lines)
        == len({r["cell_barcode"] for r in read_tsv(panel / sample / "barcode_pool.tsv")})
    )


# ── prepare: reference FASTA + GTF -> generator inputs ──────────────────────


def _reference(tmp_path: Path):
    rng = random.Random(5)
    genome = lambda: "".join(rng.choice("ACGT") for _ in range(1800))  # noqa: E731
    fasta = tmp_path / "ref.fa"
    fasta.write_text(
        "".join(f">{c}\n{genome()}\n" for c in ("NC_1.1", "NC_2.1", "NC_3.1", "NC_9.1"))
    )

    def exon(contig, start, end, strand, gene):
        return f'{contig}\tref\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}"; transcript_id "{gene}.1";\n'

    gtf = tmp_path / "ref.gtf"
    gtf.write_text(
        exon("NC_1.1", 100, 500, "+", "A1")
        + exon("NC_1.1", 600, 1000, "-", "A2")
        + exon("NC_1.1", 1100, 1200, "+", "ASPLICED")  # two blocks: skipped
        + exon("NC_1.1", 1300, 1400, "+", "ASPLICED")
        + exon("NC_1.1", 1500, 1510, "+", "ASHORT")  # < 40 bp: skipped
        + exon("NC_2.1", 100, 500, "+", "B1")
        + exon("NC_2.1", 600, 1000, "+", "B2a")  # adjacent exons merge into one block
        + exon("NC_2.1", 1001, 1200, "+", "B2a")
        + exon("NC_3.1", 100, 500, "+", "C1")  # one usable gene: virus dropped
        + exon("NC_9.1", 100, 500, "+", "Z1")  # contig not in the map
    )
    names = tmp_path / "map.tsv"
    names.write_text("NC_1.1\tHSV-1\nNC_2.1\tHSV-2\nNC_3.1\tEBV\n")
    return fasta, gtf, names


def test_prepare_builds_inputs_the_generator_accepts(tmp_path):
    fasta, gtf, names = _reference(tmp_path)
    mapping = {r[0]: r[1] for r in gt.read_rows(names)}
    result = val.prepare(fasta, gtf, tmp_path / "v.fa", tmp_path / "spans.tsv", mapping)
    assert result["viruses"] == ["HSV-1", "HSV-2"] and result["dropped_viruses"] == ["EBV"]
    assert sorted(result["skipped"]) == ["ASHORT", "ASPLICED"]
    spans = {
        tuple(line.split("\t")[:2]): line.split("\t")[2:]
        for line in (tmp_path / "spans.tsv").read_text().splitlines()
    }
    assert spans[("HSV-1", "A2")] == ["600", "1000", "-"]
    assert spans[("HSV-2", "B2a")] == ["600", "1200", "+"]  # merged
    assert ("EBV", "C1") not in spans and not any(k[1] == "Z1" for k in spans)
    manifest = gt.generate(
        gt.parse_args(
            [
                *ARGS,
                "--viral-fasta",
                str(tmp_path / "v.fa"),
                "--gene-spans",
                str(tmp_path / "spans.tsv"),
                "--sibling-pairs",
                str(_pairs(tmp_path)),
                "--out",
                str(tmp_path / "panel"),
            ]
        )
    )
    assert manifest["samples"]
    assert {
        r["virus_id"] for r in read_tsv(tmp_path / "panel" / "templates" / "viral_gene_models.tsv")
    } == {"HSV-1", "HSV-2"}


def _pairs(tmp_path: Path) -> Path:
    path = tmp_path / "pairs.tsv"
    path.write_text("HSV-1\tHSV-2\n")
    return path


def test_prepare_refuses_what_the_generator_would_refuse(tmp_path):
    fasta, gtf, names = _reference(tmp_path)
    with pytest.raises(ValueError, match="no records|needs at least 2"):
        val.prepare(fasta, gtf, tmp_path / "v.fa", tmp_path / "s.tsv", {"NC_3.1": "EBV"})


# ── plan: the command lines parse with the real parsers ─────────────────────


def _plan(tmp_path: Path, sample: str = "s0123456789ab"):
    return val.plan(
        out=tmp_path / "out",
        fasta=tmp_path / "v.fa",
        spans=tmp_path / "spans.tsv",
        onlist=tmp_path / "onlist.txt",
        index=tmp_path / "index.idx",
        t2g=tmp_path / "t2g.tsv",
        gtf=tmp_path / "ref.gtf",
        sample=sample,
        seed=7,
    )


def test_plan_commands_parse_with_the_real_parsers(tmp_path):
    from viralscan.menu import create_help

    commands = _plan(tmp_path)
    assert len(commands) == 11
    cfg = gt.parse_args(commands[0][2:])
    assert cfg["generation_seed"] == 7 and cfg["chemistries"] == ["10xv3"]
    assert cfg["error_rate"] == 0 and cfg["covariates"]["cb_umi_collisions"] == "0"
    assert cfg["viral_fasta"] and cfg["gene_spans"] and cfg["onlist"] and not cfg["production"]
    runs = [c for c in commands if c[0] == "viralscan"]
    assert len(runs) == 3
    for argv in runs:
        with patch("sys.argv", argv):
            args = create_help()
        assert args.multimap_molecule_assignments and args.cell_calling == "none"
    parser = val.build_parser()
    for argv in commands:
        if argv[0] != "viralscan" and Path(argv[1]).name == "val08_golden.py":
            parser.parse_args(argv[2:])
    # the permuted run reads the permuted FASTQs; the others the generator's own
    assert any("permuted" in part for part in runs[2]) and not any(
        "permuted" in p for p in runs[0][:12]
    )


def test_plan_prints_one_shell_quoted_command_per_line(tmp_path, capsys):
    argv = [
        "plan",
        "--out", str(tmp_path / "an out"),
        "--fasta", "v.fa", "--spans", "s.tsv", "--onlist", "o.txt", "--index", "i.idx",
        "--t2g", "t.tsv", "--gtf", "r.gtf", "--sample", "s1", "--generation-seed", "3",
    ]  # fmt: skip
    val.main(argv)
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 11 and "'" in lines[0]  # the path with a space is quoted


# ── command line ─────────────────────────────────────────────────────────────


def test_cli_check_and_compare_exit_codes(golden, tmp_path, capsys):
    panel, sample = golden
    good = fake_run(panel, sample, tmp_path / "good")
    bad = fake_run(panel, sample, tmp_path / "bad", drop_planted=2)
    val.main(["check", "--panel", str(panel), "--sample", sample, "--sample-dir", str(good)])
    assert json.loads(capsys.readouterr().out)["failures"] == []
    with pytest.raises(SystemExit) as raised:
        val.main(["check", "--panel", str(panel), "--sample", sample, "--sample-dir", str(bad)])
    assert raised.value.code == 1
    assert "not recovered" in capsys.readouterr().out
    val.main(["compare", str(good), str(good)])
    with pytest.raises(SystemExit) as raised:
        val.main(["compare", str(good), str(bad)])
    assert raised.value.code == 1


def test_cli_reports_a_missing_run_as_an_error(golden, tmp_path):
    panel, sample = golden
    with pytest.raises(SystemExit) as raised:
        val.main(
            [
                "check",
                "--panel",
                str(panel),
                "--sample",
                sample,
                "--sample-dir",
                str(tmp_path / "nope"),
            ]
        )
    assert "val08_golden" in str(raised.value.code)

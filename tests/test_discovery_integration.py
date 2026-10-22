"""DISC-10: `viralscan discover` end to end with real Kraken2 on a hand-built taxonomy database.

The database, the reads and the barcodes (real 10x v3 on-list) are generated at test time, seeded and
offline; every per-UMI outcome is designed, so the published outputs are compared to truth.
"""

import gzip
import json
import random
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from unittest.mock import patch

import pytest

from viralscan import discovery

TOOLS = ("kraken2-build", "kraken2", "k2", "bustools")
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        any(shutil.which(t) is None for t in TOOLS), reason=f"needs {', '.join(TOOLS)} on PATH"
    ),
]

SEED = 20261021
GENOME_LEN = 30_000
READ_LEN = 90
SHARED = int(0.4 * GENOME_LEN)  # E. coli / Shigella share the first 40 %
# taxid, parent, rank, name
TAXONOMY = [
    (1, 1, "no rank", "root"),
    (131567, 1, "no rank", "cellular organisms"),
    (2, 131567, "superkingdom", "Bacteria"),
    (543, 2, "family", "Enterobacteriaceae"),
    (561, 543, "genus", "Escherichia"),
    (562, 561, "species", "Escherichia coli"),
    (620, 543, "genus", "Shigella"),
    (623, 620, "species", "Shigella flexneri"),
    (1279, 2, "genus", "Staphylococcus"),
    (1280, 1279, "species", "Staphylococcus aureus"),
    (2759, 131567, "superkingdom", "Eukaryota"),
    (9605, 2759, "genus", "Homo"),
    (9606, 9605, "species", "Homo sapiens"),
    (10239, 1, "superkingdom", "Viruses"),
    (10292, 10239, "family", "Herpesviridae"),
    (10376, 10292, "species", "Human gammaherpesvirus 4"),
]
SOURCE_TAXID = {"human": 9606, "ecoli": 562, "shigella": 623, "saureus": 1280, "ebv": 10376}
CONFIDENCE = 0.5
UNIQUE = (SHARED, GENOME_LEN - READ_LEN)  # E. coli bases Shigella lacks


def _genomes(rng: random.Random) -> dict[str, str]:
    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    ecoli = rand(GENOME_LEN)
    return {
        "human": rand(GENOME_LEN),
        "ecoli": ecoli,
        "shigella": ecoli[:SHARED] + rand(GENOME_LEN - SHARED),
        "saureus": rand(GENOME_LEN),
        "ebv": rand(GENOME_LEN),
        "novel": rand(GENOME_LEN),  # not in the database
    }


def _build_db(work: Path, seqs: dict[str, str]) -> Path:
    db = work / "db"
    tax = db / "taxonomy"
    tax.mkdir(parents=True)
    with (tax / "nodes.dmp").open("w") as nodes, (tax / "names.dmp").open("w") as names:
        for taxid, parent, rank, name in TAXONOMY:
            nodes.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
            names.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    library = work / "library.fa"
    library.write_text(
        "".join(f">{k}|kraken:taxid|{t}\n{seqs[k]}\n" for k, t in SOURCE_TAXID.items())
    )
    for cmd in (
        ["kraken2-build", "--add-to-library", str(library), "--db", str(db), "--no-masking"],
        ["kraken2-build", "--build", "--db", str(db), "--threads", "2"],
    ):
        subprocess.run(cmd, check=True, capture_output=True)
    return db


def _read(seqs, source, rng, region=(0, GENOME_LEN - READ_LEN)):
    start = rng.randrange(*region)
    return seqs[source][start : start + READ_LEN]


#: 45 clean E. coli bases then 45 N: 11 clean 35-mers and 45 ambiguous ones, `562:11 A:45`.
def _partly_ambiguous(seqs) -> str:
    return seqs["ecoli"][SHARED + 1000 : SHARED + 1045] + "N" * 45


def _fixture(work: Path):
    """Run directory, designed truth {(cb, umi): taxid or None} and called/off-cell barcodes."""
    from viralscan.chemistry import onlist_path

    rng = random.Random(SEED)
    seqs = _genomes(rng)
    db = _build_db(work, seqs)
    with gzip.open(onlist_path("10xv3"), "rt") as handle:
        onlist = [next(handle).strip() for _ in range(12)]
    called, off = onlist[:10], onlist[10:]
    # kind -> (reads per UMI as source/region tuples, designed resolved taxid)
    cases = {
        "ebv": ([("ebv", None)] * 3, 10376),
        "ecoli": ([("ecoli", UNIQUE)] * 2, 562),
        "shigella": ([("shigella", (SHARED, GENOME_LEN - READ_LEN))], 623),
        "saureus": ([("saureus", None)], 1280),
        "human": ([("human", None)] * 2, 9606),
        "family": ([("ecoli", (0, SHARED - READ_LEN))], 543),  # shared E. coli / Shigella
        "novel": ([("novel", None)] * 2, None),
        "deepest": ([("ecoli", UNIQUE), ("ecoli", (0, SHARED - READ_LEN))], 562),  # 562 + 543
        "lca": ([("ecoli", UNIQUE), ("saureus", None)], 2),  # conflicting: LCA Bacteria
        "ambiguous": ([("ebv", None), ("human", None)], None),  # LCA root: never a taxon
        "denominator": ([("partial", None)], None),  # 562:11 A:45 at confidence 0.5
    }
    fragments, truth = [], {}
    for cb in called + off:
        for kind, (sources, taxid) in cases.items():
            if cb in off and kind not in ("ebv", "ecoli"):
                continue
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            for source, region in sources:
                read = (
                    _partly_ambiguous(seqs)
                    if source == "partial"
                    else _read(seqs, source, rng, region or (0, GENOME_LEN - READ_LEN))
                )
                fragments.append((cb, umi, read))
            truth[cb, umi] = taxid
    rng.shuffle(fragments)

    run = work / "run"
    sample = run / "S1"
    hf = sample / "host_filtered"
    (sample / "results").mkdir(parents=True)
    hf.mkdir()
    (sample / "config.yaml").write_text("output: x\ntechnology: 10xv3\n")
    (run / "run_manifest.json").write_text("{}\n")
    (run / "run_complete.json").write_text("{}\n")
    (sample / "results" / "called_cells.tsv").write_text(
        "barcode\n" + "".join(f"{c}-1\n" for c in called)
    )
    with gzip.open(hf / "R1.fastq.gz", "wt") as r1, gzip.open(hf / "R2.fastq.gz", "wt") as r2:
        for i, (cb, umi, read) in enumerate(fragments):
            r1.write(f"@f{i} 1:N\n{cb}{umi}\n+\n{'I' * 28}\n")
            r2.write(f"@f{i} 2:N\n{read}\n+\n{'I' * len(read)}\n")
    n = len(fragments)
    (hf / "host_filter_audit.tsv").write_text(
        "category\tfragments\tinterpretation\n"
        f"input\t{n}\tx\nretained_host_unmapped\t{n}\tx\nstar_param_set\tpinned\tx\n"
    )
    with gzip.open(hf / "fragment_lineage.tsv.gz", "wt") as handle:
        handle.write("read_id\tfilter_decision\treason\n")
        handle.writelines(f"f{i}\tretained\thost_unmapped\n" for i in range(n))
    return run, db, truth, called, off, seqs


@pytest.fixture(scope="module")
def discovered(tmp_path_factory):
    work = tmp_path_factory.mktemp("disc_integration")
    run, db, truth, called, off, seqs = _fixture(work)
    out = work / "out"
    from viralscan.menu import main

    argv = ["viralscan", "discover", "--run-dir", str(run), "-o", str(out), "--kraken2-db", str(db)]
    argv += ["--confidence", str(CONFIDENCE), "--minimum-hit-groups", "2", "--threads", "2"]
    with patch("sys.argv", argv):
        try:
            main()
            code = 0
        except SystemExit as exc:
            code = exc.code or 0
    assert code == 0
    return out, db, truth, called, off, seqs


def _tsv(path: Path):
    import csv

    with gzip.open(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_published_outputs_equal_designed_truth(discovered):
    import anndata as ad

    out, _db, truth, called, off, _seqs = discovered
    names = {p.name for p in out.iterdir()}
    assert {discovery.CELL_TAXON_H5AD, discovery.COMPLETE_JSON, discovery.TAXONOMY_TSV} <= names
    assert discovery.PARTIAL not in names
    groups = {(r["corrected_cb"], r["umi"]): r for r in _tsv(out / discovery.UMI_GROUPS_TSV)}
    assert set(groups) == set(truth)
    for key, taxid in truth.items():
        got = groups[key]["resolved_taxid"]
        assert got == ("" if taxid is None else str(taxid)), (key, got, taxid)
    classes = Counter(g["group_class"] for g in groups.values())
    assert classes["ambiguous"] == len(called)  # one EBV + human UMI per called cell
    assert classes["unclassified"] == 2 * len(called)  # the novel and the 562:11 A:45 UMI
    # Matrix: one unit per resolved UMI group of a called cell, at its resolved taxid.
    adata = ad.read_h5ad(out / discovery.CELL_TAXON_H5AD)
    assert list(adata.obs_names) == sorted(called)
    expected = Counter(
        (cb, str(t)) for (cb, _umi), t in truth.items() if t is not None and cb in set(called)
    )
    observed = Counter()
    dense = adata.X.toarray()
    for i, cb in enumerate(adata.obs_names):
        for j, taxid in enumerate(adata.var_names):
            if dense[i, j]:
                observed[cb, taxid] = int(dense[i, j])
    assert observed == expected
    offcell = Counter()
    for row in _tsv(out / discovery.OFFCELL_TSV):
        offcell[row["corrected_cb"], row["taxid"]] = int(row["n_groups"])
    assert offcell == Counter(
        (cb, str(t)) for (cb, _umi), t in truth.items() if t is not None and cb in set(off)
    )
    qc = json.loads((out / discovery.QC_JSON).read_text())
    assert qc["cells"]["n_called_cells"] == len(called)
    assert qc["called_cells_overlap"] == 1.0


def test_confidence_denominator_counts_ambiguous_kmers(discovered, tmp_path):
    """`562:11 A:45` is unclassified at 0.5 (11 / 56 < 0.5) and classified at 0.05."""
    out, db, _truth, _called, _off, seqs = discovered
    # Published: the partly ambiguous read of every called cell is its own unclassified group.
    with gzip.open(out / discovery.KRAKEN_OUT, "rt") as handle:
        raw = handle.read()
    assert "562:11 A:45" in raw
    assert all(line.split("\t")[0] == "U" for line in raw.splitlines() if "562:11 A:45" in line), (
        raw
    )
    fasta = tmp_path / "partial.fa"
    fasta.write_text(f">0\n{_partly_ambiguous(seqs)}\n")
    for confidence, status in ((0.5, "U"), (0.05, "C")):
        result = tmp_path / f"out_{confidence}.txt"
        subprocess.run(
            [
                "kraken2",
                "--db",
                str(db),
                "--confidence",
                str(confidence),
                "--minimum-hit-groups",
                "2",
                "--output",
                str(result),
                "--report",
                str(tmp_path / "report.txt"),
                str(fasta),
            ],
            check=True,
            capture_output=True,
        )
        fields = result.read_text().split("\t")
        assert fields[0] == status and "562:11 A:45" in fields[4], fields

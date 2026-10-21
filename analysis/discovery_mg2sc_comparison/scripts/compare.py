"""Synthetic truth benchmark: ViralScan discovery versus mg2sc/Cellsnake parameters and voting.

Usage (needs kraken2 2.17.2 + k2 + kraken2-build, kb's bustools, and the ViralScan source):

    python compare.py WORK_DIR OUTPUT_TSV

WORK_DIR receives the test database, the fake completed run and one discovery output per
Kraken2 setting; OUTPUT_TSV gets one row per (setting, aggregation). Everything is seeded.

What is compared, with everything else held equal:
- Kraken2 settings: ViralScan default (confidence 0.0, minimum hit groups 2), Cellsnake's
  documented default (0.05, 4; scrna-workflow config.yaml:50-51) and the Snakefile fallback
  (0.01, 3; scrna-workflow workflow/Snakefile:167-168).
- Aggregation per (CB, UMI): ViralScan's LCA of classified reads (`viralscan discover`, run end
  to end), its audit-only `deepest_consistent_taxid` column (not the approved default), versus mg2sc's majority vote `max(set(labels), key=labels.count)` over all reads,
  unclassified ("0") included (mg2sc k2sc.py:119-121, Cellsnake fork e744c2c), applied to the
  same classifications and corrected barcodes. The vote runs under several PYTHONHASHSEED values
  because its tie-break follows set iteration order.
Not emulated: mg2sc's BAM input (Cell Ranger CB/UB tags, samtools -f 4) and the fork's optional
fastp/bowtie2/second-Kraken2 host filters.
"""

import csv
import gzip
import json
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

SEED = 20261021
GENOME_LEN = 30_000
READ_LEN = 90
ERROR_RATE = 0.01
N_CELLS = 60
UMIS_PER_CELL = 20
SETTINGS = {
    "viralscan_0.0_2": (0.0, 2),
    "cellsnake_0.05_4": (0.05, 4),
    "cellsnake_fallback_0.01_3": (0.01, 3),
}
HASH_SEEDS = range(10)

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
PARENT = {t: p for t, p, _r, _n in TAXONOMY}


def lineage(taxid: int) -> list[int]:
    out = [taxid]
    while PARENT[out[-1]] != out[-1]:
        out.append(PARENT[out[-1]])
    return out


def genomes(rng: random.Random) -> dict[str, str]:
    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    human = rand(GENOME_LEN)
    ecoli = rand(GENOME_LEN)
    shared = int(0.4 * GENOME_LEN)  # E. coli / Shigella share 40 % -> family-level reads
    shigella = ecoli[:shared] + rand(GENOME_LEN - shared)
    ebv = rand(GENOME_LEN)
    ebv = ebv[:10_000] + human[5_000:7_000] + ebv[12_000:]  # 2 kb human-identical (EVE-like)
    return {
        "human": human,
        "ecoli": ecoli,
        "shigella": shigella,
        "saureus": rand(GENOME_LEN),
        "ebv": ebv,
        "novel": rand(GENOME_LEN),  # not in the database
    }


SOURCE_TAXID = {"human": 9606, "ecoli": 562, "shigella": 623, "saureus": 1280, "ebv": 10376}


def build_db(work: Path, seqs: dict[str, str]) -> Path:
    db = work / "db"
    tax = db / "taxonomy"
    tax.mkdir(parents=True)
    with (tax / "nodes.dmp").open("w") as nodes, (tax / "names.dmp").open("w") as names:
        for taxid, parent, rank, name in TAXONOMY:
            nodes.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
            names.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    lib = work / "library.fa"
    with lib.open("w") as handle:
        for key, taxid in SOURCE_TAXID.items():
            handle.write(f">{key}|kraken:taxid|{taxid}\n{seqs[key]}\n")
    for cmd in (
        ["kraken2-build", "--add-to-library", str(lib), "--db", str(db), "--no-masking"],
        ["kraken2-build", "--build", "--db", str(db), "--threads", "4"],
    ):
        subprocess.run(cmd, check=True, capture_output=True)
    return db


def mutate(read: str, rng: random.Random) -> str:
    bases = [rng.choice("ACGT".replace(b, "")) if rng.random() < ERROR_RATE else b for b in read]
    if rng.random() < 0.05:  # an N run, as in low-quality tails
        start = rng.randrange(READ_LEN - 15)
        bases[start : start + 15] = "N" * 15
    return "".join(bases)


def draw(seqs, source, rng, region=None):
    lo, hi = region or (0, GENOME_LEN - READ_LEN)
    start = rng.randrange(lo, hi)
    return mutate(seqs[source][start : start + READ_LEN], rng)


def simulate(seqs, rng, onlist):
    """Fragments (cb, umi, r2) and truth per (cb, umi)."""
    kinds = [
        ("ebv", 0.25),
        ("ebv_homolog", 0.05),
        ("ecoli", 0.12),
        ("ecoli_shared", 0.08),
        ("shigella", 0.08),
        ("saureus", 0.1),
        ("human", 0.12),
        ("novel", 0.1),
        ("collision", 0.1),
    ]
    fragments, truth = [], {}
    for cell in onlist[:N_CELLS]:
        umis = set()
        while len(umis) < UMIS_PER_CELL:
            umis.add("".join(rng.choice("ACGT") for _ in range(12)))
        for umi in sorted(umis):
            kind = rng.choices([k for k, _w in kinds], weights=[w for _k, w in kinds])[0]
            n = rng.randint(1, 5)
            if kind == "collision":
                a, b = rng.sample(["ebv", "ecoli", "saureus", "human"], 2)
                reads = [draw(seqs, a if i % 2 == 0 else b, rng) for i in range(max(n, 2))]
                truth[cell, umi] = ("collision", (SOURCE_TAXID[a], SOURCE_TAXID[b]))
            elif kind == "ebv_homolog":
                reads = [draw(seqs, "ebv", rng, (10_000, 12_000 - READ_LEN)) for _ in range(n)]
                truth[cell, umi] = ("ebv_homolog", 10376)
            elif kind == "ecoli_shared":
                reads = [
                    draw(seqs, "ecoli", rng, (0, int(0.4 * GENOME_LEN) - READ_LEN))
                    for _ in range(n)
                ]
                truth[cell, umi] = ("ecoli_shared", 562)
            elif kind == "novel":
                reads = [draw(seqs, "novel", rng) for _ in range(n)]
                truth[cell, umi] = ("novel", None)
            else:
                reads = [draw(seqs, kind, rng) for _ in range(n)]
                truth[cell, umi] = (kind, SOURCE_TAXID[kind])
            fragments += [(cell, umi, r) for r in reads]
    rng.shuffle(fragments)
    return fragments, truth


def write_run(work: Path, fragments, called) -> Path:
    run = work / "run"
    hf = run / "S1" / "host_filtered"
    hf.mkdir(parents=True)
    (run / "S1" / "results").mkdir()
    (run / "S1" / "config.yaml").write_text("output: x\ntechnology: 10xv3\n")
    (run / "run_manifest.json").write_text("{}\n")
    (run / "run_complete.json").write_text("{}\n")
    (run / "S1" / "results" / "called_cells.tsv").write_text(
        "barcode\n" + "".join(f"{c}-1\n" for c in called)
    )
    with gzip.open(hf / "R1.fastq.gz", "wt") as r1, gzip.open(hf / "R2.fastq.gz", "wt") as r2:
        for i, (cb, umi, read) in enumerate(fragments):
            r1.write(f"@f{i} 1:N\n{cb}{umi}\n+\n{'I' * 28}\n")
            r2.write(f"@f{i} 2:N\n{read}\n+\n{'I' * len(read)}\n")
    (hf / "host_filter_audit.tsv").write_text(
        "category\tfragments\tinterpretation\n"
        f"input\t{len(fragments)}\tx\nretained_host_unmapped\t{len(fragments)}\tx\nstar_param_set\tpinned\tx\n"
    )
    with gzip.open(hf / "fragment_lineage.tsv.gz", "wt") as handle:
        handle.write("read_id\tfilter_decision\treason\n")
        handle.writelines(f"f{i}\tretained\thost_unmapped\n" for i in range(len(fragments)))
    return run


def staged(out: Path, name: str) -> Path:
    """Discovery outputs: published, or staged in .partial/ before DISC-09 publication."""
    return out / name if (out / name).exists() else out / ".partial" / name


def rows(path: Path):
    with gzip.open(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


VOTE = """
import json, sys
groups = json.load(sys.stdin)
print(json.dumps({k: max(set(v), key=v.count) for k, v in groups.items()}))
"""


def mg2sc_votes(out: Path) -> list[dict[str, str]]:
    """mg2sc's vote per (CB, UMI) on discovery's own classifications, one dict per hash seed."""
    labels = defaultdict(list)
    calls = rows(staged(out, "classifications.tsv.gz"))
    for barcode in rows(staged(out, "barcodes.tsv.gz")):
        call = next(calls)
        if barcode["cb_status"] in ("onlist", "corrected"):
            labels[f"{barcode['corrected_cb']}|{barcode['umi']}"].append(call["taxid"])
    payload = json.dumps(labels)
    votes = []
    for seed in HASH_SEEDS:
        env = {**os.environ, "PYTHONHASHSEED": str(seed)}
        done = subprocess.run(
            [sys.executable, "-c", VOTE],
            input=payload,
            text=True,
            capture_output=True,
            env=env,
            check=True,
        )
        votes.append(json.loads(done.stdout))
    return votes


def score(assignments: dict[tuple[str, str], int | None], truth) -> Counter:
    """Per-UMI outcome against truth; None = no taxon (unclassified or ambiguous)."""
    out = Counter()
    for key, (kind, true) in truth.items():
        got = assignments.get(key)
        if kind == "novel":
            out["novel_false_positive" if got else "novel_correctly_unassigned"] += 1
            continue
        if kind == "collision":
            if got is None:
                out["collision_unassigned"] += 1
            elif got in true:
                out["collision_species_overcall"] += 1  # one transcript's species claimed for both
            else:
                out["collision_ancestor"] += 1
            continue
        group = "ebv_homolog" if kind == "ebv_homolog" else "single"
        out[f"{group}_{outcome(got, true)}"] += 1
    return out


def outcome(got: int | None, true: int) -> str:
    if got is None:
        return "unassigned"
    if got == true:
        return "exact"
    if got in (1, 131567):
        return "root_label"  # a root-level taxon reported as a feature
    if got in lineage(true):
        return "correct_coarser"
    return "wrong_lineage"  # e.g. EBV reads called human


def write_detail(path: Path, truth, per_setting) -> None:
    """One row per UMI: truth and each (setting, aggregation) assignment ('' = no taxon)."""
    keys = sorted(per_setting)
    with gzip.GzipFile(path, "wb", mtime=0) as raw:
        lines = ["\t".join(["cb", "umi", "kind", "true_taxid", *[f"{s}:{a}" for s, a in keys]])]
        for (cb, umi), (kind, true) in sorted(truth.items()):
            true_s = "+".join(map(str, true)) if isinstance(true, tuple) else str(true or "")
            got = [str(per_setting[k].get((cb, umi)) or "") for k in keys]
            lines.append("\t".join([cb, umi, kind, true_s, *got]))
        raw.write(("\n".join(lines) + "\n").encode())


def main(work: Path, output: Path) -> None:
    rng = random.Random(SEED)
    work.mkdir(parents=True, exist_ok=True)
    seqs = genomes(rng)
    db = build_db(work, seqs)
    onlist_path = subprocess.run(
        [
            sys.executable,
            "-c",
            "from viralscan.chemistry import onlist_path; print(onlist_path('10xv3'))",
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    with gzip.open(onlist_path, "rt") as handle:
        onlist = [next(handle).strip() for _ in range(N_CELLS)]
    fragments, truth = simulate(seqs, rng, onlist)
    run = write_run(work, fragments, onlist[:N_CELLS])
    results, per_setting = [], {}
    for name, (confidence, hit_groups) in SETTINGS.items():
        out = work / f"discover_{name}"
        cmd = [
            sys.executable,
            "-m",
            "viralscan.menu",
            "discover",
            "--run-dir",
            str(run),
            "-o",
            str(out),
            "--kraken2-db",
            str(db),
            "--confidence",
            str(confidence),
            "--minimum-hit-groups",
            str(hit_groups),
            "--threads",
            "4",
            "--quiet",
        ]
        subprocess.run(
            cmd, capture_output=True
        )  # exits 1 at the DISC-09 stop until publication lands
        ours, deepest = {}, {}
        for row in rows(staged(out, "umi_groups.tsv.gz")):
            key = row["corrected_cb"], row["umi"]
            ours[key] = int(row["resolved_taxid"]) if row["resolved_taxid"] else None
            # Audit-only variant (not the approved default): the deepest taxon when every
            # classified read lies on one root-to-leaf path, else the LCA.
            deep = row["deepest_consistent_taxid"]
            consistent = deep and row["lineage_consistent"] == "1" and int(deep) not in (1, 131567)
            deepest[key] = int(deep) if consistent else ours[key]
        votes = mg2sc_votes(out)
        first = {tuple(k.split("|")): (int(v) or None) for k, v in votes[0].items()}
        unstable = sum(len({v[k] for v in votes}) > 1 for k in votes[0])
        for aggregation, assignments in (
            ("viralscan_lca", ours),
            ("viralscan_deepest_consistent_audit", deepest),
            ("mg2sc_vote", first),
        ):
            per_setting[name, aggregation] = assignments
            row = {
                "setting": name,
                "aggregation": aggregation,
                "n_umis": len(truth),
                **score(assignments, truth),
            }
            row["vote_unstable_across_hash_seeds"] = unstable if aggregation == "mg2sc_vote" else 0
            results.append(row)
    columns = sorted({k for r in results for k in r} - {"setting", "aggregation", "n_umis"})
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, ["setting", "aggregation", "n_umis", *columns], delimiter="\t", restval=0
        )
        writer.writeheader()
        writer.writerows(results)
    write_detail(output.with_name("umi_detail.tsv.gz"), truth, per_setting)
    kinds = Counter(kind for kind, _t in truth.values())
    print(
        json.dumps(
            {"n_fragments": len(fragments), "n_umis": len(truth), "kinds": kinds}, sort_keys=True
        )
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))

#!/usr/bin/env python3
"""Read-level confirmation and artefact checks for arm-P (translated search) hits.

Arm P (``kb count --aa``) reports protein/virus UMIs but keeps no read identity. This tool works
from a BUS file written by ``kallisto bus --aa --num`` (the BUS flag column then holds the 0-based
read number), pulls the reads that support chosen proteins/viruses, and checks them:

* nucleotide alignment (minimap2, short-read preset) to the host genome, the host transcriptome
  and the viral genomes;
* translated best hit (blastx) against the viral protein set and against the human proteome;
* low complexity (DUST) and poly-A/T/G tails (tails are trimmed before alignment);
* a per-read verdict, a per-group verdict and barcode-level concentration.

Subcommands: ``build-viral-genomes`` (FASTA from the hash-verified GenBank cache), ``run`` and
``recurrence`` (proteins hit in several libraries).
Group specs for ``run --group``: a protein id (``AAD38970.1``), a taxon (``taxon_10376``) or a
gene/product within a taxon (``taxon_10376:BHLF1``); ``--top-viruses N`` adds the N viruses with
most UMIs. Host capture (what ``bustools capture --complement`` would remove) is measured as the share
of a group's reads that also align the host; nothing is filtered.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import random
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path

VERDICTS = (
    "viral_supported",
    "host_derived",
    "divergent_homolog",
    "protein_only_viral",
    "misassigned_virus",
    "adapter_artefact",
    "low_complexity",
    "unresolved",
)
DUST_MAX = 2.0  # raw DUST score (sum c(c-1)/2 / (l-1) over triplets); random 90 nt is ~0.7
MIN_TRIMMED = 30  # shorter reads after tail trimming carry no usable evidence
_TAILS = (  # (regex, side): poly-A / poly-G tails, and their reverse complements at the start
    (re.compile(r"A{10,}$"), "3"),
    (re.compile(r"G{8,}$"), "3"),
    (re.compile(r"^T{10,}"), "5"),
    (re.compile(r"^C{8,}"), "5"),
)


# --------------------------------------------------------------------------- sequence checks
# 10x 5' template-switch oligo / Drop-seq SMART primer core and its reverse complement: reads that
# contain it and little else are library artefacts (TSO-poly(T) products), not transcripts.
_TSO = "GGTATCAACGCAGAGTAC"
_TSO_RC = _TSO[::-1].translate(str.maketrans("ACGT", "TGCA"))


def strip_tso(seq: str) -> tuple[str, bool]:
    """Remove the TSO core and everything on its adapter side; return (rest, found)."""
    seq = seq.upper()
    i = seq.find(_TSO)
    if i >= 0:
        return seq[i + len(_TSO) :], True
    j = seq.find(_TSO_RC)
    if j >= 0:
        return seq[:j], True
    return seq, False


def trim_tails(seq: str) -> tuple[str, bool]:
    """Remove poly-A/poly-G tails (and reverse-complement heads); return (trimmed, was_trimmed)."""
    seq = seq.upper()
    out = seq
    for rx, _ in _TAILS:
        out = rx.sub("", out)
    return out, out != seq


def dust_score(seq: str) -> float:
    """Raw symmetric-DUST-style score over ACGT triplets (0 for sequences under 4 triplets)."""
    counts = Counter(
        seq[i : i + 3] for i in range(len(seq) - 2) if set(seq[i : i + 3]) <= set("ACGT")
    )
    n = sum(counts.values())
    if n < 4:
        return 0.0
    return sum(c * (c - 1) / 2 for c in counts.values()) / (n - 1)


def longest_run(seq: str) -> int:
    return max((len(m.group()) for m in re.finditer(r"(.)\1*", seq)), default=0)


def complexity(seq: str) -> dict:
    stripped, tso = strip_tso(seq)
    trimmed, was_trimmed = trim_tails(stripped)
    dust = dust_score(trimmed)
    return {
        "trimmed_seq": trimmed,
        "tso_adapter": tso,
        "tail_trimmed": was_trimmed,
        "dust": round(dust, 3),
        "low_complexity": len(trimmed) < MIN_TRIMMED
        or dust >= DUST_MAX
        or longest_run(trimmed) >= 15,
    }


# --------------------------------------------------------------------------- BUS / reads
def read_proteins(path: Path) -> dict[str, dict]:
    with open(path, newline="") as handle:
        return {r["target"]: r for r in csv.DictReader(handle, delimiter="\t")}


def bus_rows(bus: Path) -> Iterator[tuple[str, str, int, int]]:
    """(barcode, umi, ec, read number) from ``bustools text -p -f`` (needs ``kallisto bus --num``)."""
    cmd = ["bustools", "text", "-p", "-f", str(bus)]
    with subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            bc, umi, ec, _count, flag = line.rstrip("\n").split("\t")[:5]
            yield bc, umi, int(ec), int(flag)
    if proc.returncode:
        raise RuntimeError(f"bustools text failed for {bus}")


def read_ec(ec_lines: Iterable[str], transcripts: list[str]) -> dict[int, list[str]]:
    out = {}
    for line in ec_lines:
        a, b = line.rstrip("\n").split("\t")[:2]
        out[int(a)] = [transcripts[int(i)] for i in b.split(",")]
    return out


def parse_group_specs(specs: list[str], proteins: dict[str, dict]) -> dict[str, set[str]]:
    """group name -> set of protein targets. Raises on an unknown spec (no silent empty group)."""
    groups: dict[str, set[str]] = {}
    for spec in specs:
        if spec.startswith("taxon_") and ":" not in spec:
            hit = {t for t, r in proteins.items() if r["virus_id"] == spec}
        elif spec.startswith("taxon_"):
            taxon, key = spec.split(":", 1)
            k = key.lower()
            hit = {
                t
                for t, r in proteins.items()
                if r["virus_id"] == taxon
                and k
                in (
                    r["gene"].lower().split("/")
                    + [r["protein_id"].lower()]
                    + r["product"].lower().split()
                )
            }
        else:
            hit = {t for t, r in proteins.items() if r["protein_id"] == spec or t == spec}
        if not hit:
            raise SystemExit(f"group spec {spec!r} matches no protein in proteins.tsv")
        groups[spec] = hit
    return groups


def collect_hits(
    rows: Iterable[tuple[str, str, int, int]], ec: dict[int, list[str]], proteins: dict
) -> dict:
    """Per protein target and per virus: list of (barcode, umi, read number).

    As in ``summarize_translated.py`` only records whose equivalence class is within one virus
    count; every protein of that class receives the record (protein multimapping)."""
    prot: dict[str, list] = defaultdict(list)
    virus: dict[str, list] = defaultdict(list)
    n_ambiguous = 0
    for bc, umi, e, flag in rows:
        targets = ec[e]
        viruses = {proteins[t]["virus_id"] for t in targets}
        if len(viruses) != 1:
            n_ambiguous += 1
            continue
        virus[next(iter(viruses))].append((bc, umi, flag))
        for t in targets:
            prot[t].append((bc, umi, flag))
    return {"protein": prot, "virus": virus, "ambiguous_records": n_ambiguous}


def umi_stats(records: list[tuple[str, str, int]]) -> dict:
    """Barcode-level view of one group: ambient/systematic (many barcodes, 1-2 UMIs each) or cells."""
    umis = {(bc, umi) for bc, umi, _ in records}
    per_bc = Counter(bc for bc, _ in umis)
    n_umis, n_bc = len(umis), len(per_bc)
    if not n_umis:
        return {"reads": 0, "umis": 0, "barcodes": 0, "pattern": "none"}
    top10 = sum(c for _, c in per_bc.most_common(10)) / n_umis
    le2 = sum(1 for c in per_bc.values() if c <= 2) / n_bc
    if n_umis < 20:
        pattern = "sparse"
    elif top10 >= 0.5:
        pattern = "concentrated"
    elif n_bc >= 50 and le2 >= 0.9:
        pattern = "diffuse"
    else:
        pattern = "intermediate"
    return {
        "reads": len(records),
        "umis": n_umis,
        "barcodes": n_bc,
        "reads_per_umi": round(len(records) / n_umis, 3),
        "frac_barcodes_1_umi": round(sum(1 for c in per_bc.values() if c == 1) / n_bc, 4),
        "frac_barcodes_le2_umi": round(le2, 4),
        "median_umis_per_barcode": statistics.median(per_bc.values()),
        "max_umis_per_barcode": max(per_bc.values()),
        "top10_barcode_umi_share": round(top10, 4),
        "pattern": pattern,
    }


def sample_reads(flags: list[int], cap: int, seed: int) -> tuple[list[int], str]:
    flags = sorted(set(flags))
    if len(flags) <= cap:
        return flags, "all reads"
    return sorted(
        random.Random(seed).sample(flags, cap)
    ), f"seeded random sample of {cap} reads (seed {seed})"


def fetch_reads(fastq: Path, wanted: set[int]) -> dict[int, tuple[str, str]]:
    """read number (0-based) -> (name, sequence) from a (gzipped) FASTQ, stopping when all are found."""
    out: dict[int, tuple[str, str]] = {}
    if not wanted:
        return out
    opener = gzip.open if str(fastq).endswith(".gz") else open
    last = max(wanted)
    with opener(fastq, "rt") as handle:
        for i, name in enumerate(handle):
            seq = handle.readline().strip()
            handle.readline()
            handle.readline()
            if i in wanted:
                out[i] = (name.split()[0][1:], seq)
            if i >= last:
                break
    return out


# --------------------------------------------------------------------------- external tools
def need(tool: str) -> str:
    path = shutil.which(tool)
    if not path:
        raise SystemExit(f"{tool} not on PATH")
    return path


def parse_paf(lines: Iterable[str], per_target: bool = False) -> dict[str, dict]:
    """Best alignment per read (highest matches): identity = matches / block length.

    ``per_target``: read -> {target: best alignment}, so a read's best hit within one virus's genomes is
    not hidden by a slightly better hit to another virus."""
    best: dict[str, dict] = {}
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        qlen, qs, qe, matches, block = int(f[1]), int(f[2]), int(f[3]), int(f[9]), int(f[10])
        hit = {
            "target": f[5],
            "identity": matches / block if block else 0.0,
            "qcov": (qe - qs) / qlen if qlen else 0.0,
            "matches": matches,
        }
        slot = best.setdefault(f[0], {}) if per_target else best
        key = f[5] if per_target else f[0]
        if key not in slot or matches > slot[key]["matches"]:
            slot[key] = hit
    return best


def minimap2(
    ref: Path, query: Path, threads: int, batch: str = "300M", per_target: bool = False
) -> dict[str, dict]:
    """``-I`` bounds index memory for a human genome (batches are scored independently; we keep the best)."""
    cmd = [
        need("minimap2"),
        "-x",
        "sr",
        "-c",
        "--secondary=no",
        "-I",
        batch,
        "-t",
        str(threads),
        str(ref),
        str(query),
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return parse_paf(res.stdout.splitlines(), per_target)


def blastx(db: Path, query: Path, threads: int, evalue: float) -> dict[str, dict]:
    fmt = "6 qseqid sseqid pident length qstart qend bitscore evalue qlen"
    cmd = [need("blastx"), "-db", str(db), "-query", str(query), "-outfmt", fmt, "-evalue", str(evalue),
           "-max_target_seqs", "5", "-num_threads", str(threads)]  # fmt: skip
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    best: dict[str, dict] = {}
    for line in res.stdout.splitlines():
        q, s, pid, _length, qs, qe, bits, ev, qlen = line.split("\t")
        hit = {
            "target": s,
            "pident": float(pid),
            "qcov": (abs(int(qe) - int(qs)) + 1) / int(qlen),
            "bitscore": float(bits),
            "evalue": float(ev),
        }
        if q not in best or hit["bitscore"] > best[q]["bitscore"]:
            best[q] = hit
    return best


def make_blast_db(fasta: Path, db: Path) -> Path:
    """makeblastdb for a protein FASTA (plain or .gz via stdin); skipped if the db already exists."""
    if Path(f"{db}.pdb").exists() or Path(f"{db}.phr").exists():
        return db
    db.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        need("makeblastdb"),
        "-in",
        "-",
        "-dbtype",
        "prot",
        "-out",
        str(db),
        "-parse_seqids",
        "-title",
        db.name,
    ]
    opener = gzip.open if str(fasta).endswith(".gz") else open
    with opener(fasta, "rb") as src:  # gzip objects have no usable fileno: pipe the decoded bytes
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
        assert proc.stdin is not None
        shutil.copyfileobj(src, proc.stdin)
        proc.stdin.close()
        if proc.wait():
            raise RuntimeError(f"makeblastdb failed for {fasta}")
    return db


# --------------------------------------------------------------------------- verdicts
def read_verdict(r: dict, min_id: float, min_qcov: float) -> tuple[str, str]:
    """(verdict, reason) from one read's evidence dict."""

    def ok(hit):
        return bool(hit) and hit["identity"] >= min_id and hit["qcov"] >= min_qcov

    host = ok(r.get("host_genome")) or ok(r.get("host_tx"))
    same = ok(r.get("viral_same"))
    other = ok(r.get("viral_other"))
    if r["low_complexity"]:
        if host and not same:  # a low-complexity read that aligns the host is still host-derived
            return "host_derived", "nt_host;low_complexity"
        if r.get("tso_adapter"):
            return "adapter_artefact", "tso_adapter;low_complexity_remainder"
        return "low_complexity", "dust/tail/short"
    if same and host:
        return "unresolved", "nt_viral_and_host"
    if same:
        return "viral_supported", "nt_viral"
    if host:
        return "host_derived", "nt_host"
    if other:
        return "misassigned_virus", "nt_other_virus"
    if r.get("prot_human"):
        return "divergent_homolog", "no_nt_match;human_protein_hit"
    if r.get("prot_viral"):
        return "protein_only_viral", "no_nt_match;viral_protein_hit_only"
    if r.get("tso_adapter"):
        return "adapter_artefact", "tso_adapter;remainder_unmatched"
    return "unresolved", "no_nt_or_protein_hit"


def group_verdict(counts: Counter) -> str:
    n = sum(counts.values())
    if not n:
        return "no_reads"
    verdict, c = counts.most_common(1)[0]
    return verdict if c / n >= 0.5 else "mixed"


# --------------------------------------------------------------------------- build viral genomes
def genbank_sequence(path: Path) -> str:
    lines, on = [], False
    for line in path.read_text().splitlines():
        if line.startswith("ORIGIN"):
            on = True
        elif line.startswith("//"):
            break
        elif on:
            lines.append(re.sub(r"[^A-Za-z]", "", line))
    return "".join(lines).upper()


def build_viral_genomes(proteins: dict[str, dict], cache_dir: Path, out: Path) -> dict:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
    from viralscan.scripts.ncbi_fetch import genbank_cache_path  # noqa: PLC0415

    seen = {}
    for r in proteins.values():
        seen[r["accession"]] = r["genbank_sha256"]
    n = 0
    with open(out, "w") as fh:
        for acc, sha in sorted(seen.items()):
            gb = genbank_cache_path(acc, cache_dir)
            if not gb.is_file() or hashlib.sha256(gb.read_bytes()).hexdigest() != sha:
                raise SystemExit(f"GenBank cache for {acc} missing or hash mismatch")
            seq = genbank_sequence(gb)
            fh.write(f">{acc}\n{seq}\n")
            n += len(seq)
    return {"accessions": len(seen), "bases": n}


# --------------------------------------------------------------------------- main run
def run_confirm(
    *,
    hits: dict,
    proteins: dict,
    groups: dict[str, set[str]],
    reads_fastq: Path,
    work: Path,
    refs: dict,
    max_reads: int = 2000,
    seed: int = 1,
    threads: int = 2,
    min_id: float = 0.9,
    min_qcov: float = 0.5,
    evalue: float = 1e-3,
) -> dict:
    """Align + classify the sampled reads of every group. ``refs``: viral_genomes, host_genome,
    host_tx, viral_db, human_db (the last two are BLAST protein db prefixes)."""
    work.mkdir(parents=True, exist_ok=True)
    plan = {}
    for name, targets in groups.items():
        if name.startswith("taxon_") and ":" not in name:
            records = hits["virus"].get(name, [])
        else:
            records = sorted({x for t in targets for x in hits["protein"].get(t, [])})
        flags, how = sample_reads([f for _, _, f in records], max_reads, seed)
        plan[name] = {"records": records, "flags": flags, "sampling": how, "targets": targets}
    wanted = {f for p in plan.values() for f in p["flags"]}
    reads = fetch_reads(reads_fastq, wanted)
    cx = {i: complexity(seq) for i, (_, seq) in reads.items()}
    fa = work / "reads.fa"
    with open(fa, "w") as fh:
        for i, c in cx.items():
            if len(c["trimmed_seq"]) >= 20:
                fh.write(f">r{i}\n{c['trimmed_seq']}\n")
    nt = {k: minimap2(Path(refs[k]), fa, threads) for k in ("host_genome", "host_tx")}
    nt["viral_genomes"] = minimap2(Path(refs["viral_genomes"]), fa, threads, per_target=True)
    pv = blastx(Path(refs["viral_db"]), fa, threads, evalue) if fa.stat().st_size else {}
    ph = blastx(Path(refs["human_db"]), fa, threads, evalue) if fa.stat().st_size else {}
    acc_of_virus = defaultdict(set)
    for t, r in proteins.items():
        acc_of_virus[r["virus_id"]].add(r["accession"])
    out_groups, out_reads = {}, []
    for name, p in plan.items():
        viruses = {proteins[t]["virus_id"] for t in p["targets"]}
        accs = set().union(*(acc_of_virus[v] for v in viruses))  # every genome of the group's virus
        umi_of = defaultdict(set)
        for bc, umi, f in p["records"]:
            umi_of[f].add((bc, umi))
        counts: Counter = Counter()
        agg: Counter = Counter()
        for f in p["flags"]:
            key = f"r{f}"
            c = cx[f]
            vh = nt["viral_genomes"].get(key, {})
            same = max(
                (h for a, h in vh.items() if a in accs), key=lambda h: h["matches"], default=None
            )
            other = max(
                (h for a, h in vh.items() if a not in accs),
                key=lambda h: h["matches"],
                default=None,
            )
            row = {
                "group": name,
                "read": f,
                "barcode": sorted(umi_of[f])[0][0],
                "umi": sorted(umi_of[f])[0][1],
                "read_len": len(reads[f][1]),
                "trimmed_len": len(c["trimmed_seq"]),
                "dust": c["dust"],
                "tail_trimmed": c["tail_trimmed"],
                "low_complexity": c["low_complexity"],
                "tso_adapter": c["tso_adapter"],
                "host_genome": nt["host_genome"].get(key),
                "host_tx": nt["host_tx"].get(key),
                "viral_same": same,
                "viral_other": other,
                "prot_viral": pv.get(key),
                "prot_human": ph.get(key),
            }
            # protein-level class: best hit by bitscore (ties go to human: the conservative call)
            v, h = row["prot_viral"], row["prot_human"]
            row["prot_best"] = (
                "none"
                if not v and not h
                else "viral"
                if v and (not h or v["bitscore"] > h["bitscore"])
                else "human"
            )
            row["verdict"], row["reason"] = read_verdict(row, min_id, min_qcov)
            counts[row["verdict"]] += 1

            def good(hit):
                return bool(hit) and hit["identity"] >= min_id and hit["qcov"] >= min_qcov

            agg["host_genome"] += good(row["host_genome"])
            agg["host_tx"] += good(row["host_tx"])
            agg["host_any"] += good(row["host_genome"]) or good(row["host_tx"])
            agg["viral_same"] += good(row["viral_same"])
            agg["prot_human"] += row["prot_best"] == "human"
            agg["prot_viral"] += row["prot_best"] == "viral"
            agg["prot_none"] += row["prot_best"] == "none"
            agg["low_complexity"] += bool(row["low_complexity"])
            agg["tso_adapter"] += bool(row["tso_adapter"])
            agg["tail_trimmed"] += bool(row["tail_trimmed"])
            out_reads.append(row)
        n = len(p["flags"]) or 1
        out_groups[name] = {
            "group": name,
            "viruses": ";".join(sorted(viruses)),
            "proteins": len(p["targets"]),
            "reads_total": len({f for _, _, f in p["records"]}),
            "reads_checked": len(p["flags"]),
            "sampling": p["sampling"],
            **{
                f"pct_{k}": round(100 * agg[k] / n, 1)
                for k in (
                    "host_genome",
                    "host_tx",
                    "host_any",
                    "viral_same",
                    "prot_viral",
                    "prot_human",
                    "prot_none",
                    "low_complexity",
                    "tso_adapter",
                    "tail_trimmed",
                )
            },  # fmt: skip
            **{f"n_{v}": counts[v] for v in VERDICTS},
            "verdict": group_verdict(counts),
            **{f"bc_{k}": v for k, v in umi_stats(p["records"]).items()},
        }
    return {"groups": out_groups, "reads": out_reads}


def _cell(x):
    if isinstance(x, dict):
        return f"{x['target']}|id={x.get('identity', x.get('pident', 0)):.3f}|qcov={x['qcov']:.2f}"
    return "" if x is None else x


def write_outputs(result: dict, out_dir: Path, sample: str, meta: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    groups = list(result["groups"].values())
    with open(out_dir / f"{sample}.confirm_groups.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(groups[0]), delimiter="\t")
        w.writeheader()
        w.writerows(groups)
    keys = ["group", "read", "barcode", "umi", "read_len", "trimmed_len", "dust", "tail_trimmed",
            "low_complexity", "tso_adapter", "host_genome", "host_tx", "viral_same", "viral_other", "prot_viral",
            "prot_human", "prot_best", "verdict", "reason"]  # fmt: skip
    with open(out_dir / f"{sample}.confirm_reads.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(keys)
        for r in result["reads"]:
            w.writerow([_cell(r[k]) for k in keys])
    lines = [f"# Hit confirmation: {sample}", "", f"Parameters: `{json.dumps(meta)}`", "",
             "| group | reads (checked/total) | %host genome | %host tx | %viral nt | %prot viral/human/none | barcode pattern | verdict |",
             "|---|---|---|---|---|---|---|---|"]  # fmt: skip
    for g in groups:
        lines.append(
            f"| {g['group']} | {g['reads_checked']}/{g['reads_total']} | {g['pct_host_genome']} | {g['pct_host_tx']} "
            f"| {g['pct_viral_same']} | {g['pct_prot_viral']}/{g['pct_prot_human']}/{g['pct_prot_none']} "
            f"| {g.get('bc_pattern')} ({g.get('bc_umis')} UMIs, {g.get('bc_barcodes')} barcodes) | {g['verdict']} |"
        )
    (out_dir / f"{sample}.confirm.md").write_text("\n".join(lines) + "\n")


def top_umis(hits_by_key: dict, n: int) -> list[str]:
    umi = {k: len({(b, u) for b, u, _ in recs}) for k, recs in hits_by_key.items()}
    return [k for k, _ in sorted(umi.items(), key=lambda kv: (-kv[1], kv[0]))[:n]]


def top_group_specs(hits: dict, proteins: dict, n_viruses: int, n_proteins: int) -> list[str]:
    specs = []
    for v in top_umis(hits["virus"], n_viruses):
        specs.append(v)
        mine = {t: x for t, x in hits["protein"].items() if proteins[t]["virus_id"] == v}
        specs += top_umis(mine, n_proteins)
    return specs


def recurrence(samples: dict[str, tuple[dict, dict]], proteins: dict) -> list[dict]:
    """Per protein: UMIs/barcodes in each sample, the number of samples with hits, and each sample's share.

    ``samples``: name -> (hits from ``collect_hits``, unused). A protein seen in every sample with
    similar per-sample UMI share points to reference/host homology or a reagent, not biology."""
    names = list(samples)
    totals = {
        n: sum(len({(b, u) for b, u, _ in r}) for r in samples[n][0]["virus"].values())
        for n in names
    }
    rows = []
    for t in sorted({t for n in names for t in samples[n][0]["protein"]}):
        row = {"target": t, "virus_id": proteins[t]["virus_id"], "species": proteins[t]["species"],
               "gene": proteins[t]["gene"] or proteins[t]["product"]}  # fmt: skip
        present = 0
        for n in names:
            recs = samples[n][0]["protein"].get(t, [])
            u = {(b, x) for b, x, _ in recs}
            row[f"{n}_umis"] = len(u)
            row[f"{n}_barcodes"] = len({b for b, _ in u})
            row[f"{n}_share_pct"] = round(100 * len(u) / totals[n], 3) if totals[n] else 0.0
            present += bool(u)
        row["samples_present"] = present
        rows.append(row)
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build-viral-genomes")
    b.add_argument("--proteins-tsv", type=Path, required=True)
    b.add_argument("--cache-dir", type=Path, required=True)
    b.add_argument("--out", type=Path, required=True)
    c = sub.add_parser(
        "recurrence", help="proteins hit in several libraries (reference/reagent artefact check)"
    )
    c.add_argument("--sample", action="append", required=True, metavar="NAME=COUNTS_DIR")
    c.add_argument("--reference-dir", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    r = sub.add_parser("run")
    r.add_argument(
        "--counts-dir",
        type=Path,
        required=True,
        help="kallisto bus --aa --num output (matrix.ec, transcripts.txt)",
    )
    r.add_argument(
        "--bus", type=Path, help="BUS file (default: <counts-dir>/output.unfiltered.bus)"
    )
    r.add_argument(
        "--reads", type=Path, required=True, help="cDNA FASTQ (R2) the BUS file was made from"
    )
    r.add_argument(
        "--reference-dir", type=Path, required=True, help="holds proteins.tsv and proteins.faa"
    )
    r.add_argument("--viral-genomes", type=Path, required=True)
    r.add_argument("--host-genome", type=Path, required=True)
    r.add_argument("--host-tx", type=Path, required=True)
    r.add_argument("--human-proteins", type=Path, required=True)
    r.add_argument("--work-dir", type=Path, required=True)
    r.add_argument("--out-dir", type=Path, required=True)
    r.add_argument("--sample", required=True)
    r.add_argument("--group", action="append", default=[])
    r.add_argument("--top-viruses", type=int, default=0)
    r.add_argument("--top-proteins-per-virus", type=int, default=0)
    r.add_argument("--max-reads", type=int, default=2000)
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--threads", type=int, default=2)
    r.add_argument("--min-identity", type=float, default=0.9)
    r.add_argument("--min-qcov", type=float, default=0.5)
    a = ap.parse_args(argv)
    if a.cmd == "build-viral-genomes":
        proteins = read_proteins(a.proteins_tsv)
        print(json.dumps(build_viral_genomes(proteins, a.cache_dir, a.out)))
        return 0
    proteins = read_proteins(a.reference_dir / "proteins.tsv")
    if a.cmd == "recurrence":
        samples = {}
        for item in a.sample:
            name, _, path = item.partition("=")
            cd = Path(path)
            ec = read_ec((cd / "matrix.ec").open(), (cd / "transcripts.txt").read_text().split())
            samples[name] = (
                collect_hits(bus_rows(cd / "output.unfiltered.bus"), ec, proteins),
                None,
            )
        rows = recurrence(samples, proteins)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        return 0
    d = a.counts_dir
    ec = read_ec((d / "matrix.ec").open(), (d / "transcripts.txt").read_text().split())
    hits = collect_hits(bus_rows(a.bus or d / "output.unfiltered.bus"), ec, proteins)
    specs = list(a.group) + top_group_specs(hits, proteins, a.top_viruses, a.top_proteins_per_virus)
    groups = parse_group_specs(list(dict.fromkeys(specs)), proteins)
    viral_db = make_blast_db(a.reference_dir / "proteins.faa", a.work_dir / "db" / "viral_proteins")
    human_db = make_blast_db(a.human_proteins, a.work_dir / "db" / "human_proteins")
    refs = {"viral_genomes": a.viral_genomes, "host_genome": a.host_genome, "host_tx": a.host_tx,
            "viral_db": viral_db, "human_db": human_db}  # fmt: skip
    result = run_confirm(hits=hits, proteins=proteins, groups=groups, reads_fastq=a.reads, work=a.work_dir / a.sample,
                         refs=refs, max_reads=a.max_reads, seed=a.seed, threads=a.threads,
                         min_id=a.min_identity, min_qcov=a.min_qcov)  # fmt: skip
    meta = {"max_reads": a.max_reads, "seed": a.seed, "min_identity": a.min_identity, "min_qcov": a.min_qcov,
            "ambiguous_multi_virus_records": hits["ambiguous_records"]}  # fmt: skip
    write_outputs(result, a.out_dir, a.sample, meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

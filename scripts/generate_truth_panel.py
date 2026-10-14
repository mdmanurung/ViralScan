#!/usr/bin/env python3
"""Seeded generator for the v3 truth panel, tiny golden scale (PLAN VAL-01/02/03/05).

Code only: nothing here runs an outcome and the protocol is read, never edited. The levels, column
contracts, sibling pairs, holdout fraction and split seed all come from
``analysis/v3_validation/protocol.yaml``; the only value that is not in it is ``--generation-seed``
(``seeds.generation`` is still pending, G3 is the user's), so it is a required argument.

What it writes, per sample (``<out>/<sample_id>/``): paired FASTQs (``{id}_R1/_R2.fastq.gz``, gzip
``mtime=0`` so the bytes are reproducible) and, for training, ``barcode_pool.tsv``. Panel-wide:
``<out>/truth/{truth_manifest,host_only_manifest,host_homology_manifest,mixed_host_virus_manifest}.tsv``
(training rows only; columns are ``v3_score.required_columns(name)`` plus appended extras),
``<out>/templates/`` and ``<out>/run_manifest.json``.

Holdout (protocol D6, ``partitions.outcome_blinding``): holdout samples use opaque read names
``r<position>`` and a seeded shuffled read order. Their truth manifests, barcode pools and run
manifest go to ``--holdout-dir`` (owner-only, outside the repo); holdout generation also mixes in the
owner-held secret from ``--holdout-salt-file``. Next to the public holdout FASTQs sit only their
sha256 values.

Design defaults a reviewer must confirm (the task list in the PR note repeats them):

* the background is SYNTHETIC with known truth, so there is no ``unknown`` origin and no unplanted twin
  (``twin_sample_id`` is empty, ``background_library_id`` is the sample id: the unit is the sample);
  ``--background-r1/-r2`` (D1) is the exception: a real library, truth ``unknown``, an unplanted twin
  per planted sample and the library as the unit (design doc section 8c). Q121: each library carries
  the full abundance x host-virus-homology grid, a ``repeat_homology`` sample adding REF-07 challenge reads;
* host molecules (every background host read, including empty-droplet ambient host) are truth rows of
  ``host_only_manifest`` (one row per molecule, ``read_id`` = the molecule's first read in file order,
  extra ``host_molecule_id`` / ``n_reads``), not of ``truth_manifest``, which holds every non-host read;
* a one-sample stratum goes to training (protocol "unsplittable") BEFORE the shared split runs, so
  ``v3_partition.apportion`` is untouched; the strata are reported in ``run_manifest.json``;
* mixed host-virus molecules are an extra population (10 % of planted molecules, ceil), not part of the
  planted abundance; their basis is always ``chimeric-template``;
* synthetic barcodes are random, not on a 10x on-list, and partition membership of barcodes and UMIs is
  a hash of the sequence (a sequencing error can flip it; the truth keeps the source sequence).

ponytail: reads of the whole panel are built in memory (a hopped molecule needs its recipient's
barcode pool, and holdout needs a global shuffle). Fine at tiny scale; the production 50 M-pair
samples need a streaming rewrite, to be sized after the pilot (D5).
"""

from __future__ import annotations

import argparse
import functools
import gzip
import hashlib
import json
import math
import os
import random
import stat
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plant_anello_10x as planter  # noqa: E402
import v3_partition as part  # noqa: E402
import v3_score as score  # noqa: E402
from validate_v3_protocol import require_execution_allowed  # noqa: E402

from viralscan import __version__ as VIRALSCAN_VERSION  # noqa: E402
from viralscan import chemistry as chem_mod  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
GENERATOR_VERSION = "viralscan-truth-panel/1"
#: protocol spelling -> ``viralscan.chemistry`` name.
PACKAGE_CHEMISTRY = {"10xv2": "10xv2", "10xv3": "10xv3", "drop-seq": "dropseq"}
TARGETS = ("EBV", "HHV-6A", "HHV-6B", "HSV-1", "HSV-2", "KSHV")
VIRAL_DATASETS = ("synthetic_factorial", "synthetic_mixed_host_virus")
DATASETS = (*VIRAL_DATASETS, "synthetic_host_only", "synthetic_host_homology")
ARTEFACT_CLASSES = ("polyG", "polyA", "CAG", "TSOpolyA")
ARTEFACT_KIND = dict(zip(ARTEFACT_CLASSES, ("polyG", "polyA", "CAG", "TSO")))
MIN_BODY = 40  # a planted 3' fragment keeps >= 40 templated bases, so any aligner can see it
GENE_SPANS = ((0, 700), (600, 1300), (1200, 1800))  # overlapping, so some fragments are ambiguous
GENOME_LEN, SEGMENT_LEN, FLANK = 1800, 120, 150
BASES = planter.BASES
MANIFESTS = (
    "truth_manifest",
    "host_only_manifest",
    "host_homology_manifest",
    "mixed_host_virus_manifest",
)
ROWKEY = dict(zip(MANIFESTS, ("truth", "host_only", "homology", "mixed")))
EXTRA = {
    "truth_manifest": [
        "truth_class",
        "collision_group",
        "collision_kind",
        "templated_len",
        "host_gene_id",
        "host_molecule_id",
        "n_reads",
    ],
    "host_only_manifest": ["host_molecule_id", "n_reads"],
    "host_homology_manifest": [],
    "mixed_host_virus_manifest": [],
}
#: Appended to ``truth_manifest`` only with ``--hop-preserve-cb`` (M04-GEN): where a hopped molecule
#: came from and where its barcode landed. Absent by default, so the default file is unchanged.
HOP_EXTRA = ["hop_donor_sample_id", "hop_donor_molecule_id", "hop_cb_preserved", "hop_landing"]
COVARIATES = (
    ("infected_cell_fraction", "infected_cell_fraction_level", "0.15"),
    ("sibling_virus_similarity", "sibling_virus_similarity_level", "co_infected"),
    ("low_complexity", "low_complexity_level", "0.01"),
    ("pcr_duplication", "pcr_duplication_level", "2"),
    ("cb_umi_collisions", "cb_umi_collision_level", "0.02"),
    ("ambient_index_hopping", "ambient_index_hopping_level", "0.05"),
)
VIRAL_ONLY = {
    "infected_cell_fraction",
    "sibling_virus_similarity",
    "cb_umi_collisions",
    "ambient_index_hopping",
}


def digest(*parts: Any) -> str:
    return hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()


def side(sequence: str) -> str:
    """Partition a barcode or UMI sequence belongs to (hash of the sequence, ~30 % holdout)."""
    return part.HOLDOUT if int(digest("seq", sequence)[:8], 16) % 10 < 3 else part.TRAINING


def rand_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def err(seq: str, rate: float, rng: random.Random) -> str:
    return planter.mutate(seq, rate, rng) if rate > 0 else seq


def dup_count(rng: random.Random, mean: float) -> int:
    """Reads per source molecule: 1 + geometric (failures before success), mean ``mean``."""
    n = 1
    while rng.random() > 1 / mean:
        n += 1
    return n


def low_complexity_fraction(seq: str) -> float:
    """Share of bases in homopolymer runs of >= 6 (a coarse stand-in for the REF-07 measure)."""
    covered, i = 0, 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        covered += (j - i) if j - i >= 6 else 0
        i = j
    return covered / len(seq)


def tsv_flag(value: bool) -> str:
    return "true" if value else "false"


@dataclass
class Read:
    cb: str
    umi: str
    r1: str
    r2: str
    rows: dict[str, dict[str, Any]]  # manifest key -> columns specific to this read
    name: str = ""
    q1: str = ""  # real FASTQ qualities of a D1 background read; empty means the generator's "I"
    q2: str = ""


@dataclass
class Mol:
    virus: str
    gene: str
    cb: str
    umi: str
    seq: str
    body: int
    cls: str
    mid: str
    origin: str = "cell"
    group: str = ""
    kind: str = ""
    hop: dict[str, str] = field(default_factory=dict)  # hop_* truth columns (--hop-preserve-cb)


@dataclass
class Sample:
    dataset: str
    chemistry: str
    stratum: str
    replicate: int
    abundance: str
    homology: str
    sid: str = ""
    partition: str = part.TRAINING
    levels: dict[str, str] = field(default_factory=dict)
    background: str = (
        ""  # D1: the real background library id; empty means the sample is its own unit
    )
    twin: str = ""  # D1: the unplanted twin of this planted sample
    twin_of: str = ""  # D1: set on the twin; the planted sample it mirrors


@dataclass
class Background:
    """A real 10x library (D1): reads kept verbatim, cells anchored on read counts alone."""

    library_id: str
    chemistry: str
    reads: list[tuple[str, str, str, str, str, str]]  # cb, umi, r1, q1, r2, q2
    cells: list[str]
    empties: list[str]
    used: set[tuple[str, str]]
    entry: dict[str, Any]


# ── protocol ─────────────────────────────────────────────────────────────────


@functools.lru_cache(maxsize=1)
def load_protocol() -> dict[str, Any]:
    document: dict[str, Any] = yaml.safe_load(score.PROTOCOL.read_text())
    return document


def factor_levels(protocol: dict[str, Any]) -> dict[str, list[str]]:
    return {f["id"]: [str(level) for level in f["levels"]] for f in protocol["factors"]}


def assign_partitions(samples: list[Sample], fraction: float) -> list[list[str]]:
    """Set ``partition`` on every sample; return the unsplittable (one-sample) strata.

    Protocol "unsplittable" rule: a stratum of exactly one sample is training. ``apportion`` would
    give some such strata a seat, so they are removed before the shared split runs.
    """
    unsplittable = []
    for dataset in sorted({s.dataset for s in samples}):
        mine = [s for s in samples if s.dataset == dataset]
        strata: dict[str, list[str]] = {}
        for s in mine:
            strata.setdefault(s.stratum, []).append(s.sid)
        splittable = {k: v for k, v in strata.items() if len(v) > 1}
        unsplittable += [[dataset, k] for k, v in sorted(strata.items()) if len(v) == 1]
        labels = part.split(splittable, fraction) if splittable else {}
        for s in mine:
            s.partition = labels.get(s.sid, part.TRAINING)
    return unsplittable


# ── templates ────────────────────────────────────────────────────────────────


def build_templates(
    seed: int,
    n_host: int,
    per_virus: int,
    divergence: float,
    identity: float,
    fraction: float,
    real: RealInputs | None = None,
) -> dict[str, Any]:
    """Viral genomes (stand-ins, or ``real``), host transcripts and host loci homologous to them."""
    rng = random.Random(int(digest("templates", seed)[:16], 16))
    pair_identity: dict[str, float] = {}
    if real is None:
        genomes = {v: rand_seq(rng, GENOME_LEN) for v in TARGETS}
        for first, second in json_pairs():
            genomes[second] = planter.mutate(genomes[first], divergence, rng)
            pair_identity[f"{first}/{second}"] = round(
                sum(a == b for a, b in zip(genomes[first], genomes[second])) / GENOME_LEN, 4
            )
        genes = {v: {f"{v}_G{i + 1}": span for i, span in enumerate(GENE_SPANS)} for v in TARGETS}
        strand = {g: "+" for gs in genes.values() for g in gs}
        viruses: tuple[str, ...] = TARGETS
    else:  # real genomes are not mutated: pair identity is not computed (different lengths, indels)
        genomes, genes, strand = real.genomes, real.genes, real.strand
        viruses = tuple(sorted(genomes))
    host = {f"ENST{i:011d}": rand_seq(rng, rng.randrange(500, 900)) for i in range(n_host)}
    loci = {}
    for index, virus in enumerate(v for v in viruses for _ in range(per_virus)):
        start = rng.randrange(0, len(genomes[virus]) - SEGMENT_LEN)
        segment = planter.mutate(genomes[virus][start : start + SEGMENT_LEN], 1 - identity, rng)
        matched = sum(a == b for a, b in zip(segment, genomes[virus][start : start + SEGMENT_LEN]))
        loci[f"hc{index:03d}"] = dict(
            virus=virus,
            locus=f"synthetic_locus_{index:03d}",
            seq=rand_seq(rng, FLANK) + segment + rand_seq(rng, FLANK),
            identity=round(matched / SEGMENT_LEN, 4),
        )
    # D4: host and challenge templates are assigned whole to one partition; viral genomes are shared.
    host_side = part.split({"host": sorted(host)}, fraction)
    cluster_side = part.split({"homology": sorted(loci)}, fraction)
    return dict(
        genomes=genomes,
        genes=genes,
        strand=strand,
        host=host,
        loci=loci,
        pair_identity=pair_identity,
        host_partition=host_side,
        cluster_partition=cluster_side,
    )


def json_pairs() -> list[tuple[str, str]]:
    return [tuple(p) for p in load_protocol()["partitions"]["sibling_virus_pairs"]]  # type: ignore[misc]


# ── real reference inputs (AUD-18) ───────────────────────────────────────────

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_rows(path: Path) -> list[list[str]]:
    """Non-blank, non-``#`` lines of a small TSV, split on tabs."""
    return [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]


@dataclass
class RealInputs:
    """Real genomes and gene models; ``genes`` spans are 0-based half-open, ``strand`` is per gene."""

    genomes: dict[str, str]
    genes: dict[str, dict[str, tuple[int, int]]]
    strand: dict[str, str]
    pairs: list[tuple[str, str]]
    spans: list[dict[str, Any]]
    inputs: dict[str, Any]


def read_fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0] if line[1:].split() else ""
            if not name:
                raise ValueError(f"{path.name}: a record has an empty id")
            if name in records:
                raise ValueError(f"{path.name}: duplicate record id {name!r}")
            records[name] = ""
        elif line.strip():
            if not name:
                raise ValueError(f"{path.name}: sequence before the first '>' header")
            records[name] += line.strip().upper()
    for name, seq in records.items():
        bad = sorted(set(seq) - set("ACGTN"))
        if bad:
            raise ValueError(f"{path.name}: {name} has bases outside ACGTN: {bad}")
        if len(seq) <= SEGMENT_LEN:
            raise ValueError(f"{path.name}: {name} is {len(seq)} bp, needs > {SEGMENT_LEN}")
    if not records:
        raise ValueError(f"{path.name}: no records")
    return records


def load_real_inputs(fasta: Path, spans: Path, sibling_pairs: Path | None) -> RealInputs:
    """Validate and parse ``--viral-fasta`` / ``--gene-spans`` / ``--sibling-pairs``.

    ponytail: a window that happens to contain ``N`` is planted as is (``N`` never matches), so exact
    recovery needs gene bodies without ``N``; screen the FASTA rather than the generator.
    """
    genomes = read_fasta_records(fasta)
    rows = read_rows(spans)
    if rows and rows[0][0] == "virus_id":
        rows = rows[1:]
    genes: dict[str, dict[str, tuple[int, int]]] = {v: {} for v in genomes}
    strand: dict[str, str] = {}
    listed: list[dict[str, Any]] = []
    for row in rows:
        if len(row) != 5:
            raise ValueError(f"{spans.name}: expected virus_id, gene_id, start, end, strand: {row}")
        virus, gene, first, last, sign = row
        if virus not in genomes:
            raise ValueError(f"{spans.name}: virus {virus!r} is not in {fasta.name}")
        if any(gene in g for g in genes.values()):
            raise ValueError(f"{spans.name}: duplicate gene_id {gene!r}")
        if sign not in ("+", "-"):
            raise ValueError(f"{spans.name}: strand of {gene} must be + or -, got {sign!r}")
        try:
            start, end = int(first), int(last)
        except ValueError:
            raise ValueError(f"{spans.name}: non-integer start/end for {gene}") from None
        if not 1 <= start <= end <= len(genomes[virus]):
            raise ValueError(
                f"{spans.name}: {gene} {start}-{end} is outside {virus} (1-{len(genomes[virus])})"
            )
        if end - start + 1 < MIN_BODY:
            raise ValueError(f"{spans.name}: {gene} is shorter than {MIN_BODY} bp")
        genes[virus][gene] = (start - 1, end)
        strand[gene] = sign
        listed.append(dict(virus_id=virus, gene_id=gene, start=start, end=end, strand=sign))
    for virus, found in genes.items():
        if len(found) < 2:  # a collision moves a molecule to a different gene of the same virus
            raise ValueError(f"{spans.name}: {virus} needs at least 2 gene spans, has {len(found)}")
    if sibling_pairs is not None:
        pairs = [tuple(r) for r in read_rows(sibling_pairs)]
        for pair in pairs:
            if len(pair) != 2 or any(v not in genomes for v in pair):
                raise ValueError(f"{sibling_pairs.name}: pair {pair} is not two ids from the FASTA")
        if len({v for p in pairs for v in p}) != 2 * len(pairs):
            raise ValueError(f"{sibling_pairs.name}: a virus appears in more than one pair")
    else:  # the protocol's pairs that the FASTA carries; the rest are absent (negative) viruses
        pairs = [p for p in json_pairs() if all(v in genomes for v in p)]
    if not pairs:
        raise ValueError(
            "no sibling pair: the FASTA holds no complete protocol pair; pass --sibling-pairs"
        )
    inputs: dict[str, Any] = {
        "viral_fasta": dict(name=fasta.name, sha256=sha256_file(fasta), n_viruses=len(genomes)),
        "gene_spans": dict(name=spans.name, sha256=sha256_file(spans), n_genes=len(listed)),
    }
    if sibling_pairs is not None:
        inputs["sibling_pairs"] = dict(name=sibling_pairs.name, sha256=sha256_file(sibling_pairs))
    return RealInputs(genomes, genes, strand, pairs, listed, inputs)  # type: ignore[arg-type]


def load_onlist(path: Path) -> tuple[list[str], dict[str, Any]]:
    """Sorted, validated barcodes of ``--onlist`` (plain or ``.gz``) and its manifest entry."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:  # type: ignore[call-overload]
        barcodes = [line.strip().upper() for line in handle if line.strip()]
    if not barcodes:
        raise ValueError(f"{path.name}: empty on-list")
    if len({len(b) for b in barcodes}) != 1:
        raise ValueError(f"{path.name}: barcodes differ in length")
    if set("".join(barcodes)) - set("ACGT"):
        raise ValueError(f"{path.name}: barcodes must be ACGT only")
    if len(set(barcodes)) != len(barcodes):
        raise ValueError(f"{path.name}: duplicate barcodes")
    entry = dict(
        name=path.name, sha256=sha256_file(path), n_barcodes=len(barcodes), length=len(barcodes[0])
    )
    return sorted(barcodes), entry


def load_background(
    r1: Path, r2: Path, chemistry: str, n_pairs: int, n_cells: int, n_empty: int
) -> Background:
    """Read the first ``n_pairs`` (0 = all) of a real FASTQ pair and anchor its cells (D1).

    Cells are the ``n_cells`` most frequent barcodes by read count and the empty droplets the next
    ``n_empty``: a count rank, so no ViralScan output decides where molecules are planted. It stands
    in for the protocol's ``cell_universe`` anchor, which replaces it on a production panel.

    ponytail: the pairs are held in memory, and every sample (planted and twin) copies them. Fine at
    golden scale; cap with ``--n-background``, and stream for a production library.
    """
    geometry = chem_mod.get(PACKAGE_CHEMISTRY[chemistry])
    key = geometry.cb_len + geometry.umi_len
    reads = []
    for rec1, rec2 in planter.iter_background(r1, r2, n_pairs or sys.maxsize):
        seq1, seq2 = rec1[1].strip(), rec2[1].strip()
        if (
            len(seq1) != key
        ):  # a shorter or longer R1 would sit beside planted reads of length `key`
            raise ValueError(
                f"{r1.name}: R1 is {len(seq1)} nt but {chemistry} needs {key} (barcode + UMI); trim it"
            )
        reads.append(
            (
                seq1[: geometry.cb_len],
                seq1[geometry.cb_len :],
                seq1,
                rec1[3].strip(),
                seq2,
                rec2[3].strip(),
            )
        )
    if not reads:
        raise ValueError(f"{r1.name}: no read pairs")
    counts: dict[str, int] = {}
    for cb, *_ in reads:
        if not set(cb) - set("ACGT"):
            counts[cb] = counts.get(cb, 0) + 1
    ranked = sorted(counts, key=lambda b: (-counts[b], b))
    if len(ranked) < n_cells + n_empty:
        raise ValueError(
            f"{r1.name}: {len(ranked)} distinct ACGT barcodes, {n_cells + n_empty} "
            f"(= n_cells + n_empty) needed"
        )
    sums = (sha256_file(r1), sha256_file(r2))
    entry = dict(
        library_id="bg" + digest(*sums)[:12],
        r1=dict(name=r1.name, sha256=sums[0]),
        r2=dict(name=r2.name, sha256=sums[1]),
        chemistry=chemistry,
        n_pairs=len(reads),
        n_cells=n_cells,
        n_empty=n_empty,
        cell_anchor="barcodes ranked by read count",
        min_cell_reads=counts[ranked[n_cells - 1]],
    )
    return Background(
        entry["library_id"],
        chemistry,
        reads,
        ranked[:n_cells],
        ranked[n_cells : n_cells + n_empty],
        {(cb, umi) for cb, umi, *_ in reads},
        entry,
    )


def fasta(records: dict[str, str]) -> str:
    return "".join(f">{name}\n{seq}\n" for name, seq in sorted(records.items()))


# ── panel ────────────────────────────────────────────────────────────────────


class Panel:
    def __init__(self, cfg: dict[str, Any]):
        self.cfg = cfg
        self.protocol = load_protocol()
        self.levels = factor_levels(self.protocol)
        self.fraction = float(self.protocol["partitions"]["holdout_fraction"])
        self.split_seed = int(cfg["split_seed"] or self.protocol["seeds"]["split"]["value"])
        self.salt = cfg["salt"]
        real: RealInputs | None = cfg.get("real")
        self.pairs = real.pairs if real else json_pairs()
        self.tpl = build_templates(
            cfg["generation_seed"],
            cfg["n_host_transcripts"],
            cfg["clusters_per_virus"],
            cfg["sibling_divergence"],
            cfg["homology_identity"],
            self.fraction,
            real,
        )
        self.viruses: tuple[str, ...] = tuple(self.tpl["genomes"])
        self.onlist: dict[str, list[str]] | None = cfg.get("onlist_pools")
        self.bg: Background | None = cfg.get("background")
        self.samples: dict[str, Sample] = {}
        self.unsplittable: list[list[str]] = []
        self._layout: dict[str, dict[str, Any]] = {}
        self.incoming: dict[str, list[Read]] = {}
        self.realised: dict[str, dict[str, int]] = {}
        self.enumerate_samples()

    # seeds ---------------------------------------------------------------
    def rng(self, s: Sample, stage: str) -> random.Random:
        """Stage seed per ``seeds.derivation``; holdout samples also mix in the secret salt."""
        text = [
            self.protocol["protocol_version"],
            s.dataset,
            s.sid,
            f"rep{s.replicate}",
            s.stratum,
            "generator",
            stage,
            self.cfg["generation_seed"],
        ]
        if s.partition == part.HOLDOUT:
            text.append(self.salt.hex())
        return random.Random(int(digest(*text)[:16], 16))

    # enumeration ---------------------------------------------------------
    def enumerate_samples(self) -> None:
        cfg, levels = self.cfg, self.levels
        for name, factor in (("abundance", "viral_abundance"), ("chemistries", "chemistry")):
            for level in cfg[name]:
                if str(level) not in levels[factor]:
                    raise ValueError(
                        f"{factor} level {level!r} is not in the protocol {levels[factor]}"
                    )
        for factor, _, _ in COVARIATES:
            wanted = cfg["covariates"][factor]
            if wanted != "cycle" and wanted not in levels[factor]:
                raise ValueError(
                    f"{factor} level {wanted!r} is not in the protocol {levels[factor]}"
                )
        for dataset in cfg["datasets"]:
            if dataset not in DATASETS:
                raise ValueError(f"unknown dataset {dataset!r}")
        if self.bg is not None:
            if sorted(cfg["datasets"]) != ["synthetic_factorial"]:
                raise ValueError("a real background (D1) plants synthetic_factorial only")
            if sorted(cfg["chemistries"]) != [self.bg.chemistry]:
                raise ValueError(
                    f"a real background is one chemistry: --chemistries {self.bg.chemistry}"
                )
        out = []
        for dataset in sorted(cfg["datasets"]):
            abundances = (
                [str(a) for a in cfg["abundance"]] if dataset in VIRAL_DATASETS else ["none"]
            )
            for chem in sorted(cfg["chemistries"]):
                for abundance in sorted(abundances):
                    for homology in self.homology_levels(dataset):
                        keys = [abundance] if dataset in VIRAL_DATASETS else []
                        if dataset == "synthetic_host_homology" or self.crosses_homology(dataset):
                            keys.append(homology)
                        stratum = "|".join([*keys, chem])
                        for rep in range(1, cfg["replicates"] + 1):
                            sid = part.sample_id(self.split_seed, dataset, stratum, rep)
                            out.append(
                                Sample(dataset, chem, stratum, rep, abundance, homology, sid)
                            )
        if (
            self.bg is not None
        ):  # D1: every planted sample gets an unplanted twin of the same library
            for planted in list(out):
                planted.background = self.bg.library_id
                twin_sid = part.sample_id(
                    self.split_seed, planted.dataset, planted.stratum + "|twin", planted.replicate
                )
                planted.twin = twin_sid
                out.append(
                    Sample(
                        planted.dataset,
                        planted.chemistry,
                        planted.stratum,
                        planted.replicate,
                        "none",
                        "none",
                        twin_sid,
                        background=self.bg.library_id,
                        twin_of=planted.sid,
                    )
                )
        self.samples = {s.sid: s for s in sorted(out, key=lambda s: s.sid)}
        if len(self.samples) != len(out):
            raise ValueError("sample id collision")
        if self.bg is None:
            self.unsplittable = assign_partitions(list(self.samples.values()), self.fraction)
        else:
            # D1: the unit is the donor library, never the sample (SCI-05 r13 M5). The label is the
            # one split_libraries draws over every donor library of the chemistry, shared by every
            # panel run; without --donor-libraries one library cannot fill a holdout seat, so it
            # is training.
            label = part.TRAINING
            if donors := cfg.get("donor_libraries"):
                if self.bg.library_id not in donors:
                    raise ValueError(
                        f"--donor-libraries must list this library ({self.bg.library_id})"
                    )
                label = part.split_libraries(
                    {self.bg.chemistry: donors}, self.split_seed, self.fraction
                )[self.bg.library_id]
            for s in self.samples.values():
                s.partition = label
            self.unsplittable = []
        for s in self.samples.values():
            s.levels = self.covariate_levels(s)
        if self.onlist is not None:
            self.check_onlist()

    def crosses_homology(self, dataset: str) -> bool:
        """Q121: on a real background every donor library carries the abundance x homology grid."""
        return self.bg is not None and dataset == "synthetic_factorial"

    def homology_levels(self, dataset: str) -> list[str]:
        if self.crosses_homology(dataset):
            return ["none", "repeat_homology"]
        return ["repeat_homology" if dataset == "synthetic_host_homology" else "none"]

    def check_onlist(self) -> None:
        """The on-list must cover each chemistry's barcode length and each partition's pool."""
        assert self.onlist is not None
        length = len(next(b for pool in self.onlist.values() for b in pool))
        need = self.cfg["n_cells"] + self.cfg["n_empty"]
        for chem in sorted({s.chemistry for s in self.samples.values()}):
            cb_len = chem_mod.get(PACKAGE_CHEMISTRY[chem]).cb_len
            if cb_len != length:
                raise ValueError(
                    f"--onlist barcodes are {length} nt but {chem} uses {cb_len} nt; "
                    f"restrict --chemistries"
                )
        for side_ in sorted({s.partition for s in self.samples.values()}):
            if len(self.onlist[side_]) < need:
                raise ValueError(
                    f"--onlist holds {len(self.onlist[side_])} {side_} barcodes (the partition is a "
                    f"hash of the sequence); {need} (= n_cells + n_empty) are needed"
                )

    def covariate_levels(self, s: Sample) -> dict[str, str]:
        out = {}
        for factor, column, _ in COVARIATES:
            wanted = self.cfg["covariates"][factor]
            if (s.dataset not in VIRAL_DATASETS or s.twin_of) and factor in VIRAL_ONLY:
                out[column] = "none"
            elif self.bg is not None and factor == "low_complexity":
                out[column] = "0"  # nothing is added to a real background except planted molecules
            elif wanted == "cycle":
                order = list(range(len(self.levels[factor])))
                random.Random(
                    int(
                        digest("cycle", self.cfg["generation_seed"], s.dataset, s.stratum, factor)[
                            :16
                        ],
                        16,
                    )
                ).shuffle(order)
                out[column] = self.levels[factor][order[(s.replicate - 1) % len(order)]]
            else:
                out[column] = wanted
        return out

    # layout --------------------------------------------------------------
    def layout(self, s: Sample) -> dict[str, Any]:
        """Barcode pool and designated infected cells; shared with hop donors, so cached."""
        if s.sid in self._layout:
            return self._layout[s.sid]
        cfg = self.cfg
        cb_len = chem_mod.get(PACKAGE_CHEMISTRY[s.chemistry]).cb_len
        rng = self.rng(s, "barcodes")
        pool: list[str] = []
        if self.bg is not None:  # D1: the real library's own cells and empty droplets
            pool = self.bg.cells + self.bg.empties
        elif self.onlist is not None:  # seeded draw without replacement from the real on-list
            pool = rng.sample(self.onlist[s.partition], cfg["n_cells"] + cfg["n_empty"])
        while len(pool) < cfg["n_cells"] + cfg["n_empty"]:
            cb = rand_seq(rng, cb_len)
            if cb not in pool and side(cb) == s.partition:
                pool.append(cb)
        cells, empties = pool[: cfg["n_cells"]], pool[cfg["n_cells"] :]
        designated: dict[str, list[str]] = {}
        absent: list[str] = []
        if s.dataset in VIRAL_DATASETS and not s.twin_of:
            rng = self.rng(s, "designation")
            k = max(1, math.ceil(float(s.levels["infected_cell_fraction_level"]) * len(cells)))
            mode = s.levels["sibling_virus_similarity_level"]
            if mode == "co_infected":
                units = [list(p) for p in self.pairs]
            elif mode == "disjoint_cells":
                units = [[v] for p in self.pairs for v in p]
            else:  # one_member_only: a seeded member of each pair
                units = [[rng.choice(p)] for p in self.pairs]
            if len(units) * k > len(cells):
                raise ValueError(
                    f"{len(units)} disjoint sets of {k} cells do not fit {len(cells)} cells"
                )
            order = cells[:]
            rng.shuffle(order)
            for i, unit in enumerate(units):
                for virus in unit:
                    designated[virus] = order[i * k : (i + 1) * k]
            absent = [v for v in self.viruses if v not in designated]
        infected = {c for cs in designated.values() for c in cs}
        out = dict(
            cells=cells,
            empties=empties,
            designated=designated,
            absent=absent,
            uninfected=[c for c in cells if c not in infected] + empties,
        )
        self._layout[s.sid] = out
        return out

    # reads ---------------------------------------------------------------
    def fresh_umi(self, rng: random.Random, s: Sample, cb: str, used: set[tuple[str, str]]) -> str:
        n = chem_mod.get(PACKAGE_CHEMISTRY[s.chemistry]).umi_len
        while True:
            umi = rand_seq(rng, n)
            if (
                len(set(umi)) > 1
                and "AAAAAA" not in umi
                and (cb, umi) not in used
                and side(umi) == s.partition
            ):
                used.add((cb, umi))
                return umi

    def emit(
        self,
        rng: random.Random,
        cb: str,
        umi: str,
        seq: str,
        rows: dict[str, dict[str, Any]],
        mean: float,
    ) -> list[Read]:
        rate = self.cfg["error_rate"]
        return [
            Read(cb, umi, err(cb + umi, rate, rng), err(seq, rate, rng), dict(rows))
            for _ in range(dup_count(rng, mean))
        ]

    def fragment(
        self, rng: random.Random, virus: str, gene: str | None = None
    ) -> tuple[str, str, str, int, str]:
        """3' fragment of a viral gene: ``(gene, read, truth_class, templated_len, ...)``."""
        genes, strand = self.tpl["genes"][virus], self.tpl["strand"]
        gene = gene or rng.choice(sorted(genes))
        first, end = genes[gene]
        if strand[gene] == "+":
            start = rng.randrange(max(first, end - 300), end - MIN_BODY + 1)
            read, body = planter.polya_read(self.tpl["genomes"][virus], start, end, rng)
            lo, hi = start, min(start + planter.READ_LEN, end)
        else:  # the transcript is the reverse complement; its 3' end is the gene's left edge
            transcript = revcomp(self.tpl["genomes"][virus][first:end])
            length = end - first
            start = rng.randrange(max(0, length - 300), length - MIN_BODY + 1)
            read, body = planter.polya_read(transcript, start, length, rng)
            lo, hi = end - min(start + planter.READ_LEN, length), end - start
        # ambiguous: the read's genomic footprint fits more than one same-strand gene model
        compatible = [
            g for g, (a, b) in genes.items() if strand[g] == strand[gene] and a <= lo and hi <= b
        ]
        return gene, read, "unique" if len(compatible) == 1 else "ambiguous", body, virus

    def molecule(
        self,
        rng: random.Random,
        s: Sample,
        virus: str,
        cb: str,
        used: set[tuple[str, str]],
        counter: list[int],
        gene: str | None = None,
    ) -> Mol:
        gene, read, cls, body, _ = self.fragment(rng, virus, gene)
        counter[0] += 1
        return Mol(
            virus,
            gene,
            cb,
            self.fresh_umi(rng, s, cb, used),
            read,
            body,
            cls,
            f"{s.sid}:m{counter[0]}",
        )

    @staticmethod
    def viral_rows(m: Mol) -> dict[str, dict[str, Any]]:
        return {
            "truth": dict(
                origin=m.origin,
                planted_virus_id=m.virus,
                planted_gene_id=m.gene,
                planted_molecule_id=m.mid,
                is_infected_cell=tsv_flag(m.origin == "cell"),
                is_planted_molecule="true",
                truth_class=m.cls,
                collision_group=m.group,
                collision_kind=m.kind,
                templated_len=m.body,
                **m.hop,
            )
        }

    def host_reads(
        self,
        rng: random.Random,
        s: Sample,
        cb: str,
        n_mol: int,
        used: set[tuple[str, str]],
        pool: list[str],
        weights: list[float],
        counter: list[int],
    ) -> list[Read]:
        out = []
        for _ in range(n_mol):
            tx = rng.choices(pool, weights)[0]
            seq = self.tpl["host"][tx]
            start = rng.randrange(max(0, len(seq) - 300), len(seq) - MIN_BODY)
            read, _ = planter.polya_read(seq, start, len(seq), rng)
            counter[0] += 1
            mid = f"h:{tx}:{counter[0]}"
            out += self.emit(
                rng,
                cb,
                self.fresh_umi(rng, s, cb, used),
                read,
                {"host_only": dict(host_gene_id=tx, host_molecule_id=mid)},
                self.pcr(s),
            )
        return out

    @staticmethod
    def pcr(s: Sample) -> float:
        return float(s.levels["pcr_duplication_level"])

    def build_reads(self, s: Sample) -> list[Read]:
        """Every read of one sample; hopped copies for other samples go to ``self.incoming``."""
        if self.bg is not None:
            return self.background_reads(s)
        cfg, lay = self.cfg, self.layout(s)
        reads: list[Read] = []
        used: set[tuple[str, str]] = set()
        counter = [0]
        cells, empties = lay["cells"], lay["empties"]
        # host background: expression-weighted synthetic transcripts of this sample's partition
        rng = self.rng(s, "host")
        pool = sorted(t for t, p in self.tpl["host_partition"].items() if p == s.partition)
        if not pool:
            raise ValueError(
                f"no host transcript in partition {s.partition}; raise --n-host-transcripts"
            )
        weights = [rng.lognormvariate(0, 1) for _ in pool]
        per_cell = max(1, round(cfg["reads_per_cell"] / self.pcr(s)))
        for cb in cells:
            reads += self.host_reads(rng, s, cb, per_cell, used, pool, weights, counter)
        # ambient-only content: log-normal UMI total, median partitions.generated_sample_structure
        # .empty_median_umi (SCI-05 r13: was a hard-coded 50)
        median = self.protocol["partitions"]["generated_sample_structure"]["empty_median_umi"]
        for cb in empties:
            n = max(1, round(rng.lognormvariate(math.log(median), 0.5)))
            reads += self.host_reads(rng, s, cb, n, used, pool, weights, counter)
        if s.dataset in VIRAL_DATASETS:
            reads += self.viral_reads(s, lay, used, counter)
        if s.dataset == "synthetic_host_homology":
            reads += self.challenge_reads(s, lay, used)
        reads += self.artefact_reads(s, lay, used)
        return reads

    def background_reads(self, s: Sample) -> list[Read]:
        """D1: the real library's reads verbatim (truth ``unknown``), then the planted molecules.

        The twin is the same reads with nothing added, so a call in it is background signal.
        """
        assert self.bg is not None
        truth = dict(origin="unknown", is_infected_cell="false", is_planted_molecule="false")
        reads = [
            Read(cb, umi, r1, r2, {"truth": dict(truth)}, q1=q1, q2=q2)
            for cb, umi, r1, q1, r2, q2 in self.bg.reads
        ]
        if s.twin_of:
            return reads
        used = set(self.bg.used)
        reads += self.viral_reads(s, self.layout(s), used, [0])
        if (
            s.homology == "repeat_homology"
        ):  # Q121: REF-07 challenge reads beside the planted molecules
            reads += self.challenge_reads(s, self.layout(s), used)
        return reads

    def viral_reads(
        self, s: Sample, lay: dict[str, Any], used: set[tuple[str, str]], counter: list[int]
    ) -> list[Read]:
        rng = self.rng(s, "viral")
        abundance, lv = int(s.abundance), s.levels
        mols: dict[str, list[Mol]] = {
            v: [self.molecule(rng, s, v, rng.choice(cs), used, counter) for _ in range(abundance)]
            for v, cs in sorted(lay["designated"].items())
        }
        flat = [m for v in sorted(mols) for m in mols[v]]
        claimed: set[str] = set()
        done = dict(
            collision_same=0,
            collision_disjoint=0,
            ambient=0,
            hopped=0,
            mixed=0,
            requested_collision=0,
            requested_ambient=0,
            requested_hopped=0,
        )
        # CB/UMI collisions: n = max(2, ceil(level * molecules)), half "same" key and half 1 Hamming away
        level = float(lv["cb_umi_collision_level"])
        n_coll = max(2, math.ceil(level * len(flat))) if level > 0 else 0
        done["requested_collision"] = n_coll
        for k in range(n_coll):
            virus = sorted(mols)[k % len(mols)]
            free = [m for m in mols[virus] if m.mid not in claimed]
            if len(free) < 2:
                continue
            keep, moved = rng.sample(free, 2)
            same = k % 2 == 0
            others = [g for g in self.tpl["genes"][virus] if g != keep.gene]
            gene, read, cls, body, _ = self.fragment(rng, virus, rng.choice(others))
            umi = keep.umi
            if not same:
                options = [
                    keep.umi[:i] + b + keep.umi[i + 1 :]
                    for i in range(len(keep.umi))
                    for b in BASES
                    if b != keep.umi[i]
                ]
                rng.shuffle(options)
                options = [
                    u for u in options if side(u) == s.partition and (keep.cb, u) not in used
                ]
                if not options:
                    continue
                umi = options[0]
            used.add((keep.cb, umi))
            moved.cb, moved.umi, moved.gene, moved.seq, moved.cls, moved.body = (
                keep.cb,
                umi,
                gene,
                read,
                cls,
                body,
            )
            moved.group = keep.group = f"{s.sid}:c{k}"
            moved.kind = keep.kind = "same" if same else "disjoint"
            claimed |= {keep.mid, moved.mid}
            done["collision_same" if same else "collision_disjoint"] += 1
        # ambient (relocated within the sample) and hopped (copied to a same-partition sample)
        level = float(lv["ambient_index_hopping_level"])
        for virus in sorted(mols):
            n_total = max(2, math.ceil(level * abundance)) if level > 0 else 0
            done["requested_ambient"] += n_total // 2 + n_total % 2
            done["requested_hopped"] += n_total // 2
            free = [m for m in mols[virus] if m.mid not in claimed]
            rng.shuffle(free)
            for m in free[: n_total // 2 + n_total % 2]:
                if not lay["uninfected"]:
                    break
                m.cb, m.origin = rng.choice(lay["uninfected"]), "ambient"
                m.umi = self.fresh_umi(rng, s, m.cb, used)
                claimed.add(m.mid)
                done["ambient"] += 1
            for m in free[n_total // 2 + n_total % 2 : n_total]:
                if self.hop(rng, s, m):
                    done["hopped"] += 1
        reads: list[Read] = []
        for m in flat:
            reads += self.emit(rng, m.cb, m.umi, m.seq, self.viral_rows(m), self.pcr(s))
        if s.dataset == "synthetic_mixed_host_virus":
            reads += self.mixed_reads(s, lay, rng, used, counter, flat, done)
        self.realised[s.sid] = done
        return reads

    def hop(self, rng: random.Random, s: Sample, m: Mol) -> bool:
        """Copy ``m`` into another viral-dataset sample of the same chemistry and partition."""
        options = [
            o
            for o in self.samples.values()
            if o.dataset in VIRAL_DATASETS
            and (o.chemistry, o.partition) == (s.chemistry, s.partition)
            and (o.background or o.sid) != (s.background or s.sid)
        ]
        if not options:
            return False
        target = rng.choice(options)  # `samples` is sorted by id, so the choice is stable
        if self.cfg.get("hop_preserve_cb"):
            return self.hop_keeping_barcode(rng, s, m, target)
        copy = Mol(
            m.virus,
            m.gene,
            rng.choice(self.layout(target)["uninfected"]),
            m.umi,
            m.seq,
            m.body,
            m.cls,
            f"{m.mid}>{target.sid}",
            "hopped",
        )
        self.incoming.setdefault(target.sid, []).extend(
            self.emit(rng, copy.cb, copy.umi, copy.seq, self.viral_rows(copy), self.pcr(s))
        )
        return True

    def hop_keeping_barcode(self, rng: random.Random, s: Sample, m: Mol, target: Sample) -> bool:
        """Index hopping proper: the copy keeps the donor's ``(CB, UMI)`` and sequence.

        Where it lands depends only on the recipient's barcode space: a called cell, an empty droplet,
        or another barcode that the recipient's matrix holds as non-called. The generator does not
        steer it into cells; with random barcodes it is almost always ``other``, and with ``--onlist``
        it is a cell by on-list coincidence, as it is in a real pooled run.
        """
        lay = self.layout(target)
        landing = "cell" if m.cb in lay["cells"] else "empty" if m.cb in lay["empties"] else "other"
        copy = Mol(
            m.virus,
            m.gene,
            m.cb,
            m.umi,
            m.seq,
            m.body,
            m.cls,
            f"{m.mid}>{target.sid}",
            "hopped",
            hop=dict(
                hop_donor_sample_id=s.sid,
                hop_donor_molecule_id=m.mid,
                hop_cb_preserved="true",
                hop_landing=landing,
            ),
        )
        self.incoming.setdefault(target.sid, []).extend(
            self.emit(rng, copy.cb, copy.umi, copy.seq, self.viral_rows(copy), self.pcr(s))
        )
        return True

    def mixed_reads(
        self,
        s: Sample,
        lay: dict[str, Any],
        rng: random.Random,
        used: set[tuple[str, str]],
        counter: list[int],
        flat: list[Mol],
        done: dict[str, int],
    ) -> list[Read]:
        pool = sorted(t for t, p in self.tpl["host_partition"].items() if p == s.partition)
        reads: list[Read] = []
        n_mixed = max(1, math.ceil(0.1 * len(flat)))
        first = rng.randrange(2)
        for k in range(n_mixed):
            virus = rng.choice(sorted(lay["designated"]))
            m = self.molecule(rng, s, virus, rng.choice(lay["designated"][virus]), used, counter)
            truth_class = ("host", "virus")[(k + first) % 2]
            tx = rng.choice(pool)
            seq = self.tpl["host"][tx]
            hread, _ = planter.polya_read(
                seq, rng.randrange(max(0, len(seq) - 300), len(seq) - MIN_BODY), len(seq), rng
            )
            mixed = dict(
                truth_class=truth_class,
                host_gene_id=tx,
                viral_virus_id=virus,
                viral_gene_id=m.gene,
                mixture_basis="chimeric-template",
                corrected_cell_barcode=m.cb,
                corrected_umi=m.umi,
            )
            vrows = (
                self.viral_rows(m)
                if truth_class == "virus"
                else {
                    "truth": dict(
                        origin="artefact:chimera",
                        is_infected_cell="false",
                        is_planted_molecule="false",
                        collision_kind="",
                        templated_len=m.body,
                    )
                }
            )
            reads += self.emit(rng, m.cb, m.umi, m.seq, {**vrows, "mixed": mixed}, self.pcr(s))
            reads += self.emit(rng, m.cb, m.umi, hread, {"mixed": mixed}, self.pcr(s))
            done["mixed"] += 1
        return reads

    def challenge_reads(
        self, s: Sample, lay: dict[str, Any], used: set[tuple[str, str]]
    ) -> list[Read]:
        rng = self.rng(s, "challenge")
        reads: list[Read] = []
        for cluster in sorted(
            c for c, p in self.tpl["cluster_partition"].items() if p == s.partition
        ):
            locus = self.tpl["loci"][cluster]
            for _ in range(self.cfg["challenge_molecules"]):
                cb = rng.choice(lay["cells"])
                start = FLANK + rng.randrange(0, SEGMENT_LEN - planter.READ_LEN + 1)
                seq = locus["seq"][start : start + planter.READ_LEN]
                rows = {
                    "truth": dict(origin="challenge", templated_len=planter.READ_LEN),
                    "homology": dict(
                        host_source_locus=locus["locus"],
                        homology_cluster_id=cluster,
                        homologous_virus_id=locus["virus"],
                        homology_identity=locus["identity"],
                        low_complexity_fraction=round(low_complexity_fraction(seq), 4),
                        planted_as_host="true",
                    ),
                }
                reads += self.emit(
                    rng, cb, self.fresh_umi(rng, s, cb, used), seq, rows, self.pcr(s)
                )
        return reads

    def artefact_reads(
        self, s: Sample, lay: dict[str, Any], used: set[tuple[str, str]]
    ) -> list[Read]:
        level = float(s.levels["low_complexity_level"])
        if level <= 0:
            return []
        rng = self.rng(s, "artefacts")
        total = max(4, math.ceil(level * self.cfg["n_cells"] * self.cfg["reads_per_cell"]))
        reads = []
        for i in range(total):
            kind = ARTEFACT_CLASSES[i % 4]
            cb = rng.choice(lay["cells"] + lay["empties"])
            umi = self.fresh_umi(rng, s, cb, used)
            seq = planter.artefact_read(ARTEFACT_KIND[kind], rng)
            reads.append(
                Read(
                    cb,
                    umi,
                    cb + umi,
                    seq,
                    {
                        "truth": dict(
                            origin=f"artefact:{kind}",
                            is_infected_cell="false",
                            is_planted_molecule="false",
                            templated_len=0,
                        )
                    },
                )
            )
        return reads


# ── output ───────────────────────────────────────────────────────────────────


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fastq_gz(path: Path, records: list[tuple[str, str, str]]) -> str:
    """Write gzip with ``mtime=0`` and no filename; return the sha256 of the uncompressed stream."""
    raw = hashlib.sha256()
    with (
        open(path, "wb") as handle,
        gzip.GzipFile(filename="", fileobj=handle, mode="wb", mtime=0, compresslevel=6) as gz,
    ):
        for name, seq, qual in records:
            data = f"@{name}\n{seq}\n+\n{qual or 'I' * len(seq)}\n".encode()
            raw.update(data)
            gz.write(data)
    return raw.hexdigest()


def write_text(path: Path, text: str, private: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700 if private else 0o777)
    flags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
    with os.fdopen(os.open(path, flags, 0o600 if private else 0o666), "w", newline="") as handle:
        handle.write(text)
    if private:
        os.chmod(path, stat.S_IRUSR | stat.S_IWUSR)


def tsv(columns: list[str], rows: list[dict[str, Any]]) -> str:
    lines = ["\t".join(columns)]
    lines += ["\t".join(str(row.get(c, "")) for c in columns) for row in rows]
    return "\n".join(lines) + "\n"


def assemble(panel: Panel, s: Sample, reads: list[Read]) -> list[Read]:
    """Final read order and names. Holdout: seeded shuffle and opaque ``r<n>``; training keeps the
    generation order and names host reads by their molecule (``h:<ENST>:<n>:<k>``)."""
    reads = reads + panel.incoming.get(s.sid, [])
    if s.partition == part.HOLDOUT:
        panel.rng(s, "shuffle").shuffle(reads)
        for i, read in enumerate(reads):
            read.name = f"r{i}"
        return reads
    seen: dict[str, int] = {}
    for i, read in enumerate(reads):
        host = read.rows.get("host_only")
        if host:
            seen[host["host_molecule_id"]] = seen.get(host["host_molecule_id"], 0) + 1
            read.name = f"{host['host_molecule_id']}:{seen[host['host_molecule_id']]}"
        else:
            read.name = f"p{i}"
    return reads


def sample_rows(panel: Panel, s: Sample, reads: list[Read]) -> dict[str, list[dict[str, Any]]]:
    common = dict(
        dataset_id=s.dataset,
        condition_id=s.stratum,
        biological_sample_id=s.sid,
        replicate_id=f"rep{s.replicate}",
        background_library_id=s.background or s.sid,
        twin_sample_id=s.twin,
        chemistry=s.chemistry,
        viral_abundance_level=s.abundance,
        host_virus_homology_level=s.homology,
        generator_version=GENERATOR_VERSION,
        generator_seed=panel.cfg["generation_seed"],
        partition=s.partition,
        **s.levels,
    )
    counts: dict[str, int] = {}
    for read in reads:
        host = read.rows.get("host_only")
        if host:
            counts[host["host_molecule_id"]] = counts.get(host["host_molecule_id"], 0) + 1
    # is_infected_cell is a fact about the cell: it received >= 1 planted molecule of origin "cell"
    infected = {
        r.cb
        for r in reads
        if r.rows.get("truth", {}).get("origin") == "cell"
        and r.rows["truth"].get("is_planted_molecule") == "true"
    }
    empties = set(panel.layout(s)["empties"])
    out: dict[str, list[dict[str, Any]]] = {m: [] for m in MANIFESTS}
    emitted: set[str] = set()
    for read in reads:
        base = common | dict(
            cell_barcode=read.cb,
            umi=read.umi,
            read_id=read.name,
            record_granularity="read",
            is_infected_cell=tsv_flag(read.cb in infected),
        )
        if "truth" in read.rows:
            out["truth_manifest"].append(
                base | read.rows["truth"] | dict(is_infected_cell=base["is_infected_cell"])
            )
        if "homology" in read.rows:
            out["host_homology_manifest"].append(base | read.rows["homology"])
        if "mixed" in read.rows:
            out["mixed_host_virus_manifest"].append(base | read.rows["mixed"])
        host = read.rows.get("host_only")
        if host and host["host_molecule_id"] not in emitted:
            emitted.add(host["host_molecule_id"])
            molecule = (
                base
                | host
                | dict(record_granularity="molecule", n_reads=counts[host["host_molecule_id"]])
            )
            if s.dataset == "synthetic_host_only":
                out["host_only_manifest"].append(molecule)
            else:  # host molecules of a planted sample are molecule rows of the truth manifest
                out["truth_manifest"].append(
                    molecule
                    | dict(
                        origin="ambient" if read.cb in empties else "cell",
                        is_planted_molecule="false",
                    )
                )
    return out


def production_structure_errors(cfg: dict[str, Any], protocol: dict[str, Any]) -> list[str]:
    """Where a ``--production`` config departs from the frozen panel structure (SCI-05 r12 M4).

    ``partitions.generated_sample_structure`` gives cells and reads per cell; ``seeds.split`` is
    the split seed. The empty-barcode count of a synthetic background (B15: 3,000) is prose in the
    protocol, so it is read only from a structured ``empty_barcodes_per_sample`` field and a
    synthetic run is refused while that field is absent. A real background brings its own uncalled
    barcodes (not a count), so ``n_empty`` is not checked there.
    """
    structure = protocol["partitions"]["generated_sample_structure"]
    errors = []
    for arg, key in (("n_cells", "cells_per_sample"), ("reads_per_cell", "reads_per_cell")):
        if cfg[arg] != structure[key]:
            errors.append(f"{arg}={cfg[arg]} differs from the frozen {key}={structure[key]}")
    if not cfg.get("background_r1"):
        if "empty_barcodes_per_sample" not in structure:
            errors.append(
                "the protocol has no structured empty_barcodes_per_sample (B15 is prose only); "
                "a synthetic background cannot be checked"
            )
        elif cfg["n_empty"] != structure["empty_barcodes_per_sample"]:
            errors.append(
                f"n_empty={cfg['n_empty']} differs from the frozen "
                f"empty_barcodes_per_sample={structure['empty_barcodes_per_sample']}"
            )
    frozen_split = int(protocol["seeds"]["split"]["value"])
    if cfg["split_seed"] and cfg["split_seed"] != frozen_split:
        errors.append(f"split_seed={cfg['split_seed']} differs from the frozen seeds.split")
    return errors


def production_design_errors(cfg: dict[str, Any], protocol: dict[str, Any]) -> list[str]:
    """Where a ``--production`` config does not cover the frozen design (SCI-05 r13 M4, M5).

    ``partitions.minimum_samples_rationale`` makes VAL-01 recompute the floors from the frozen level
    counts, and the LOD estimator needs every target planted at every abundance level within each
    chemistry, so a production panel must carry every abundance and chemistry level, every
    ``training-and-holdout`` dataset, and every covariate level (``cycle`` with at least as many
    replicates as levels). ``synthetic_factorial`` without a real background is refused (M5, r14), and a real background is refused outright (blocker ``real_background_anchor``): the protocol anchors its cells on
    the shared host-only EmptyDrops universe and draws its empties from barcodes outside that
    anchor with a host-only UMI total in [10, ``host_only_emptydrops.lower``) (B15), and the
    generator has neither input; ranking barcodes by read count is only a golden-scale stand-in.
    """
    levels = factor_levels(protocol)
    errors = []
    if cfg.get("background_r1"):
        errors.append(
            "a real background (--background-r1) is refused under --production: its cells must be "
            "the frozen host-only EmptyDrops anchor (harmonization.cell_universe) and its empties "
            "the B15 barcodes outside it, not barcodes ranked by read count"
        )
    if "synthetic_factorial" in (cfg.get("datasets") or []) and not cfg.get("background_r1"):
        # SCI-05 r14 M5: its planting_rule needs real donor libraries (no abundance x homology
        # crossing exists on a synthetic background), so clearing the ``simulator`` blocker alone
        # must not open a production path around the real-background refusal above.
        errors.append(
            "synthetic_factorial needs a real healthy-donor PBMC background under --production "
            "(synthetic_factorial.planting_rule); none was supplied"
        )
    for arg, factor in (("abundance", "viral_abundance"), ("chemistries", "chemistry")):
        have, want = sorted(map(str, cfg.get(arg) or [])), sorted(levels[factor])
        if have != want:
            errors.append(f"{arg}={have} does not cover every frozen {factor} level {want}")
    roles = protocol["partitions"]["applies_to_dataset_roles"]
    want_sets = sorted(d["id"] for d in protocol["datasets"] if d["analysis_role"] in roles)
    if sorted(cfg.get("datasets") or []) != want_sets:
        errors.append(f"datasets={sorted(cfg.get('datasets') or [])} differ from {want_sets}")
    for factor, wanted in (cfg.get("covariates") or {}).items():
        if wanted != "cycle":
            errors.append(f"covariate {factor}={wanted!r} pins one level; production needs 'cycle'")
        elif cfg.get("replicates", 0) < len(levels[factor]):
            errors.append(
                f"replicates={cfg.get('replicates')} cannot realise all {len(levels[factor])} "
                f"{factor} levels"
            )
    return errors


def generate(cfg: dict[str, Any]) -> dict[str, Any]:
    """Generate the panel described by ``cfg`` (see :func:`main` for the keys) and return the manifest."""
    out = Path(cfg["out"]).resolve()
    holdout_dir = Path(cfg["holdout_dir"]).resolve() if cfg.get("holdout_dir") else None
    for label, path in (("--out", out), ("--holdout-dir", holdout_dir)):
        if path is not None and (path == REPO or REPO in path.parents):
            raise ValueError(
                f"{label} must be outside the repository tree (governance: no generated panels in git)"
            )
    if cfg["production"]:
        # SCI-05 r12 M4: the full gate (status, freeze record, blockers, factors, seeds) first.
        try:
            require_execution_allowed("training")
        except SystemExit as exc:
            raise ValueError(f"production run refused: {exc}") from exc
        protocol = load_protocol()
        pending = [f["id"] for f in protocol["factors"] if f["status"] != "frozen"]
        if pending:
            raise ValueError(f"production run refused: factors not frozen (G3): {pending}")
        if cfg["replicates"] < protocol["partitions"]["minimum_samples_per_stratum"]:
            raise ValueError("production run refused: replicates below minimum_samples_per_stratum")
        if errors := production_structure_errors(cfg, protocol) + production_design_errors(
            cfg, protocol
        ):
            raise ValueError("production run refused: " + "; ".join(errors))
    cfg = dict(cfg, salt=b"")
    inputs: dict[str, Any] = {}
    if bool(cfg.get("viral_fasta")) != bool(cfg.get("gene_spans")):
        raise ValueError("--viral-fasta and --gene-spans must be given together")
    if cfg.get("sibling_pairs") and not cfg.get("viral_fasta"):
        raise ValueError("--sibling-pairs needs --viral-fasta")
    if cfg.get("viral_fasta"):
        pairs_file = Path(cfg["sibling_pairs"]) if cfg.get("sibling_pairs") else None
        cfg["real"] = load_real_inputs(
            Path(cfg["viral_fasta"]), Path(cfg["gene_spans"]), pairs_file
        )
        inputs.update(cfg["real"].inputs)
    if bool(cfg.get("background_r1")) != bool(cfg.get("background_r2")):
        raise ValueError("--background-r1 and --background-r2 must be given together")
    if cfg.get("background_r1"):
        if cfg.get("onlist"):
            raise ValueError("--onlist is for synthetic barcodes; a real background brings its own")
        cfg["background"] = load_background(
            Path(cfg["background_r1"]),
            Path(cfg["background_r2"]),
            cfg["background_chemistry"],
            cfg["n_background"],
            cfg["n_cells"],
            cfg["n_empty"],
        )
        inputs["background"] = cfg["background"].entry
    if cfg.get("onlist"):
        barcodes, inputs["onlist"] = load_onlist(Path(cfg["onlist"]))
        cfg["onlist_pools"] = {
            p: [b for b in barcodes if side(b) == p] for p in (part.TRAINING, part.HOLDOUT)
        }
    # salt is needed iff some sample lands in holdout, which needs the partition first
    probe = Panel(cfg)
    has_holdout = any(s.partition == part.HOLDOUT for s in probe.samples.values())
    if has_holdout:
        if holdout_dir is None or not cfg.get("holdout_salt_file"):
            raise ValueError(
                "this panel has holdout samples: --holdout-dir and --holdout-salt-file are required"
            )
        salt_path = Path(cfg["holdout_salt_file"])
        if salt_path.stat().st_mode & 0o077:
            raise ValueError(f"{salt_path} must be owner-only (chmod 600)")
        cfg["salt"] = salt_path.read_bytes().strip()
        if len(cfg["salt"]) < 16:
            raise ValueError("holdout salt must be at least 16 bytes")
    if has_holdout:
        assert holdout_dir is not None
        holdout_dir.mkdir(parents=True, exist_ok=True)
        holdout_dir.chmod(0o700)
    panel = Panel(cfg)
    order = list(panel.samples.values())  # sorted by sample id
    own = {s.sid: panel.build_reads(s) for s in order}

    public: dict[str, Any] = {"samples": {}, "holdout_samples": []}
    private: dict[str, Any] = {"samples": {}}
    rows: dict[str, dict[str, list[dict[str, Any]]]] = {
        side_: {m: [] for m in MANIFESTS} for side_ in (part.TRAINING, part.HOLDOUT)
    }
    for s in order:
        reads = assemble(panel, s, own[s.sid])
        sdir = out / s.sid
        sdir.mkdir(parents=True, exist_ok=True)
        r1, r2 = sdir / f"{s.sid}_R1.fastq.gz", sdir / f"{s.sid}_R2.fastq.gz"
        raw1 = write_fastq_gz(r1, [(r.name, r.r1, r.q1) for r in reads])
        raw2 = write_fastq_gz(r2, [(r.name, r.r2, r.q2) for r in reads])
        files = dict(
            r1_sha256=sha256_file(r1),
            r2_sha256=sha256_file(r2),
            r1_raw_sha256=raw1,
            r2_raw_sha256=raw2,
            n_read_pairs=len(reads),
        )
        lay = panel.layout(s)
        pool = tsv(
            ["cell_barcode", "role"],
            [{"cell_barcode": c, "role": "cell"} for c in lay["cells"]]
            + [{"cell_barcode": c, "role": "empty"} for c in lay["empties"]],
        )
        for manifest, found in sample_rows(panel, s, reads).items():
            rows[s.partition][manifest] += found
        detail = dict(
            dataset=s.dataset,
            chemistry=s.chemistry,
            stratum=s.stratum,
            replicate=s.replicate,
            levels=s.levels,
            absent_viruses=lay["absent"],
            realised=panel.realised.get(s.sid, {}),
            **files,
        )
        if s.background:  # D1: pair each sample with its twin, and name the unit that is resampled
            detail.update(
                background_library_id=s.background, twin_sample_id=s.twin, twin_of=s.twin_of
            )
        if s.partition == part.HOLDOUT:
            assert holdout_dir is not None
            write_text(holdout_dir / s.sid / "barcode_pool.tsv", pool, True)
            write_text(
                sdir / f"{s.sid}.sha256",
                f"{files['r1_sha256']}  {r1.name}\n{files['r2_sha256']}  {r2.name}\n",
                False,
            )
            public["holdout_samples"].append(
                dict(sample_id=s.sid, r1_sha256=files["r1_sha256"], r2_sha256=files["r2_sha256"])
            )
            private["samples"][s.sid] = detail
        else:
            write_text(sdir / "barcode_pool.tsv", pool, False)
            public["samples"][s.sid] = dict(detail, partition=s.partition)
    for side_, base, private_ in (
        (part.TRAINING, out / "truth", False),
        (part.HOLDOUT, holdout_dir and holdout_dir / "truth", True),
    ):
        if base is None or (side_ == part.HOLDOUT and not has_holdout):
            continue
        for manifest in MANIFESTS:
            columns = score.required_columns(manifest) + EXTRA[manifest]
            if cfg.get("hop_preserve_cb") and manifest == "truth_manifest":
                columns += HOP_EXTRA
            write_text(base / f"{manifest}.tsv", tsv(columns, rows[side_][manifest]), private_)
            if manifest == "truth_manifest":
                score.read_truth_manifest(
                    base / f"{manifest}.tsv"
                )  # the generator's own output must validate
    templates = out / "templates"
    write_text(templates / "viral_genomes.fasta", fasta(panel.tpl["genomes"]), False)
    write_text(templates / "host_transcripts.fasta", fasta(panel.tpl["host"]), False)
    write_text(
        templates / "host_homology_loci.fasta",
        fasta({k: v["seq"] for k, v in panel.tpl["loci"].items()}),
        False,
    )
    gene_seqs = {
        g: (revcomp if panel.tpl["strand"][g] == "-" else str)(panel.tpl["genomes"][v][a:b])
        for v in panel.viruses
        for g, (a, b) in panel.tpl["genes"][v].items()
    }
    write_text(templates / "viral_genes.fasta", fasta(gene_seqs), False)
    models = (  # real input: the user's 1-based inclusive spans and strand, echoed back
        tsv(["virus_id", "gene_id", "start", "end", "strand"], cfg["real"].spans)
        if cfg.get("real")
        else tsv(
            ["virus_id", "gene_id", "start", "end"],
            [
                dict(virus_id=v, gene_id=g, start=a, end=b)
                for v in TARGETS
                for g, (a, b) in panel.tpl["genes"][v].items()
            ],
        )
    )
    write_text(templates / "viral_gene_models.tsv", models, False)
    write_text(
        templates / "transcripts.fasta",
        fasta(gene_seqs | panel.tpl["host"] | {k: v["seq"] for k, v in panel.tpl["loci"].items()}),
        False,
    )
    manifest_doc = dict(
        generator_version=GENERATOR_VERSION,
        viralscan_version=VIRALSCAN_VERSION,
        generation_seed=cfg["generation_seed"],
        split_seed=panel.split_seed,
        protocol_version=panel.protocol["protocol_version"],
        protocol_sha256=sha256_file(score.PROTOCOL),
        panel={
            k: sorted(cfg[k]) if k in ("datasets", "chemistries", "abundance") else cfg[k]
            for k in (
                "datasets",
                "chemistries",
                "abundance",
                "replicates",
                "n_cells",
                "n_empty",
                "reads_per_cell",
                "error_rate",
                "covariates",
                "n_host_transcripts",
                "clusters_per_virus",
                "challenge_molecules",
                "sibling_divergence",
                "homology_identity",
            )
        },
        sibling_pair_identity=panel.tpl["pair_identity"],
        unsplittable_strata=panel.unsplittable,
        templates={p.name: sha256_file(p) for p in sorted(templates.iterdir())},
        truth={m: sha256_file(out / "truth" / f"{m}.tsv") for m in MANIFESTS},
        training_host_transcripts=sorted(
            t for t, p in panel.tpl["host_partition"].items() if p == part.TRAINING
        ),
        training_homology_clusters=sorted(
            t for t, p in panel.tpl["cluster_partition"].items() if p == part.TRAINING
        ),
        **public,
    )
    if inputs:
        manifest_doc["inputs"] = inputs
    if cfg.get("hop_preserve_cb"):
        manifest_doc["hop_mode"] = "cb-preserving"
    if has_holdout:
        assert holdout_dir is not None
        manifest_doc["holdout_salt_sha256"] = hashlib.sha256(cfg["salt"]).hexdigest()
        private.update(
            truth={m: sha256_file(holdout_dir / "truth" / f"{m}.tsv") for m in MANIFESTS},
            holdout_host_transcripts=sorted(
                t for t, p in panel.tpl["host_partition"].items() if p == part.HOLDOUT
            ),
            holdout_homology_clusters=sorted(
                t for t, p in panel.tpl["cluster_partition"].items() if p == part.HOLDOUT
            ),
        )
        write_text(
            holdout_dir / "run_manifest.json",
            json.dumps(private, indent=1, sort_keys=True) + "\n",
            True,
        )
        for directory in (d for d in holdout_dir.rglob("*") if d.is_dir()):
            directory.chmod(0o700)  # mkdir(parents=True) leaves the intermediate ones at the umask
    write_text(
        out / "run_manifest.json", json.dumps(manifest_doc, indent=1, sort_keys=True) + "\n", False
    )
    return manifest_doc


def parse_args(argv: list[str] | None = None) -> dict[str, Any]:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--generation-seed",
        type=int,
        required=True,
        help="public generator seed (seeds.generation is pending)",
    )
    p.add_argument("--out", required=True, help="output directory, outside the repository")
    p.add_argument(
        "--holdout-dir",
        help="owner-only directory for holdout truth (required if any sample is holdout)",
    )
    p.add_argument("--holdout-salt-file", help="owner-only file holding the secret holdout salt")
    p.add_argument("--split-seed", type=int, help="default: seeds.split of the protocol")
    p.add_argument(
        "--datasets",
        nargs="+",
        default=None,
        choices=DATASETS,
        help="default: all four; with --background-r1, synthetic_factorial",
    )
    p.add_argument(
        "--chemistries",
        nargs="+",
        default=None,
        help="default: 10xv2 10xv3 drop-seq; with --background-r1, --background-chemistry",
    )
    p.add_argument(
        "--abundance", nargs="+", default=["10"], help="planted molecules per virus per sample"
    )
    p.add_argument("--replicates", type=int, default=1, help="samples per stratum")
    p.add_argument("--n-cells", type=int, default=50)
    p.add_argument("--n-empty", type=int, default=20)
    p.add_argument("--reads-per-cell", type=int, default=200)
    p.add_argument(
        "--error-rate", type=float, default=0.001, help="per-base substitution rate on every read"
    )
    p.add_argument("--n-host-transcripts", type=int, default=40)
    p.add_argument("--clusters-per-virus", type=int, default=3)
    p.add_argument(
        "--challenge-molecules", type=int, default=3, help="per homology cluster per sample"
    )
    p.add_argument("--sibling-divergence", type=float, default=0.10)
    p.add_argument("--homology-identity", type=float, default=0.90)
    p.add_argument(
        "--viral-fasta",
        help="real viral genomes (record id = virus id); with --gene-spans replaces the stand-in "
        "random genomes. Viruses in no sibling pair are planted nowhere (absent controls)",
    )
    p.add_argument(
        "--gene-spans",
        help="TSV virus_id, gene_id, start, end, strand (1-based inclusive, strand + or -); at "
        "least 2 genes per virus. truth_class unique/ambiguous follows these models",
    )
    p.add_argument(
        "--sibling-pairs",
        help="TSV of two virus ids per line (disjoint pairs); default with real genomes: the "
        "protocol pairs whose members are both in --viral-fasta",
    )
    p.add_argument(
        "--onlist",
        help="barcode on-list, one ACGT barcode per line (.gz ok); cell, empty-droplet, hopped and "
        "collision barcodes are drawn from it (seeded). Its length must match every chemistry in "
        "--chemistries (16 for 10xv2/10xv3)",
    )
    p.add_argument(
        "--background-r1",
        help="D1: R1 of a real 10x library (gzip, FASTQ). With --background-r2 it replaces the "
        "synthetic host background: every planted sample is the library's reads verbatim plus planted "
        "molecules, and has an unplanted twin (the same reads, nothing added). Cells are the "
        "--n-cells most frequent barcodes, empty droplets the next --n-empty",
    )
    p.add_argument("--background-r2", help="D1: R2 of the real library")
    p.add_argument(
        "--donor-libraries",
        nargs="+",
        help="D1: ids of every donor library of this chemistry in the full design (this library "
        "included). The library's partition is the one v3_partition.split_libraries assigns; "
        "default: training",
    )
    p.add_argument(
        "--background-chemistry",
        default="10xv3",
        choices=sorted(PACKAGE_CHEMISTRY),
        help="D1: chemistry of the real library; its R1 must be exactly barcode + UMI long",
    )
    p.add_argument(
        "--n-background",
        type=int,
        default=0,
        help="D1: use the first N read pairs of the library (0 = all), the same block in both arms",
    )
    p.add_argument(
        "--hop-preserve-cb",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="hopped molecules keep the donor (CB, UMI), as index hopping does, instead of taking a "
        "random uninfected barcode of the recipient; truth_manifest gains hop_donor_sample_id, "
        "hop_donor_molecule_id, hop_cb_preserved and hop_landing (cell, empty or other). Default "
        "off: the output is then byte-identical to a run without the option",
    )
    p.add_argument(
        "--production", action="store_true", help="enforce the frozen-level and sample-floor gates"
    )
    for factor, _, default in COVARIATES:
        p.add_argument(
            f"--{factor.replace('_', '-')}", default=default, help="a protocol level, or 'cycle'"
        )
    ns = vars(p.parse_args(argv))
    ns["covariates"] = {f: ns.pop(f) for f, _, _ in COVARIATES}
    real = bool(ns["background_r1"])
    ns["datasets"] = ns["datasets"] or (["synthetic_factorial"] if real else list(DATASETS))
    ns["chemistries"] = ns["chemistries"] or (
        [ns["background_chemistry"]] if real else ["10xv2", "10xv3", "drop-seq"]
    )
    return ns


def main(argv: list[str] | None = None) -> None:
    try:
        manifest = generate(parse_args(argv))
    except (ValueError, OSError) as exc:
        sys.exit(f"generate_truth_panel: {exc}")
    n = len(manifest["samples"]) + len(manifest["holdout_samples"])
    print(
        f"wrote {n} samples ({len(manifest['holdout_samples'])} holdout) under the output directory"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Optional decoy layer for false-discovery calibration (PLAN CAL-13; software only, nothing is run here).

    make   panel FASTA(s) [+ unrelated-host viral FASTA(s)] -> decoys.fa, decoys.gtf, decoys_manifest.tsv
    fdr    per-library target/decoy call table -> stratified empirical p, target-decoy q, expected
           false calls per library at a score threshold, with an interval.

``make``
    Each panel genome is replaced by a **k-let-preserving shuffle** (default k = 6, the uShuffle
    6-let of the discovery protocol): a uniformly random sequence with exactly the same k-mer counts
    and the same first and last k-1 bases (Altschul and Erickson 1985, Eulerian path on the
    (k-1)-mer de Bruijn graph; the random last-exit tree is drawn with Wilson's algorithm, so every
    admissible sequence is equally likely). stdlib only; ``ushuffle`` is not needed. Each decoy gets
    its own RNG stream from sha256(seed, k, decoy id), so output does not depend on record order.
    ``--unrelated-host-fasta`` registers genomes of viruses of unrelated hosts, taken as they are
    (no download code), as the second decoy class. Every decoy id is ``DECOY_`` + source id (+ ``_<rep>``
    for shuffles), so it cannot collide with a real accession or virus key. The GTF has one
    gene/transcript/exon per decoy over the whole sequence (``gene_biotype "whole_genome"``, the same
    writer as ``build-ref``). Build the decoy index separately, with the same host and D-list.

``fdr``  (input: ``--calls`` TSV with library_id, target_id, is_decoy, stratum, score)
    One row per searched entry, **zero scores included**, one library at a time (a decoy that is
    missing from a library would shrink m; ``--manifest`` checks that all manifest decoys are present).
    Everything is per library, so the null is depth-matched; nothing is pooled across libraries.
    * ``p_empirical = (b + 1) / (m + 1)``: m decoys in the same library and stratum, b of them with
      score >= the target's.
    * ``q_td`` (stratified target-decoy, the protocol's ``target_decoy_stratified``): per library and
      stratum, FDR(s) = min(1, (D(s) + 1) / (r * T(s))), D and T the decoys and targets with score >= s,
      r = decoys / targets of the stratum; q = lowest FDR over all thresholds <= s. The +1 is the
      conservative estimator, so one call among 4 decoys cannot reach q below 1 / T.
      A stratum with targets but no decoys gets q = 1.
    * ``q_bh``: Benjamini-Hochberg on ``p_empirical`` over all targets of the library (strata pooled).
    * Expected false calls at ``--threshold`` (score >= threshold is a call): sum over strata of
      n_targets * D / m. Interval: Monte Carlo over Beta(D + 1/2, m - D + 1/2) (Jeffreys), so it is not
      zero-width when D = 0. Strata with targets but no decoys are excluded and counted.
    No threshold, target or interval is built in: the threshold is an argument; the pre-registered
    target and interval come from ``--protocol`` (``false_discovery.expected_false_L2_calls_per_library``
    and ``false_discovery.interval`` = [lo, hi]; values still starting with ``TBD`` are ignored).

    python scripts/cal13_decoys.py make --panel-fasta panel.fa --out-dir D --seed 1 [--k 6] [--copies 1]
    python scripts/cal13_decoys.py fdr --calls calls.tsv --out-dir S --threshold 3 --seed 1
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import random
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from viralscan.enrichment import _bh_adjust  # noqa: E402
from viralscan.run_safety import sha256_file  # noqa: E402
from viralscan.scripts.build_reference import (  # noqa: E402
    _fasta_records,
    _genome_as_transcript_gtf,
    _write_fasta,
)

DEFAULT_K = 6
PREFIX = "DECOY_"
MANIFEST_COLUMNS = (
    "decoy_id",
    "class",
    "source_id",
    "seed",
    "k",
    "rep",
    "length",
    "sha256",
    "source_sha256",
)
CALL_COLUMNS = ("library_id", "target_id", "is_decoy", "stratum", "score")
TRUE, FALSE = {"true", "1", "yes", "t", "y"}, {"false", "0", "no", "f", "n"}


def _die(message: str) -> None:
    sys.exit(f"cal13_decoys: {message}")


# ---------------------------------------------------------------- k-let shuffle


def klet_shuffle(seq: str, k: int, rng: random.Random) -> str:
    """Uniform random sequence with the same k-let counts and the same first/last k-1 symbols."""
    n = len(seq)
    if k < 1:
        raise ValueError("k must be >= 1")
    if n < k:
        raise ValueError(f"sequence length {n} < k = {k}")
    out: dict[str, list[str]] = defaultdict(list)  # (k-1)-mer -> last symbols of its out-edges
    for i in range(n - k + 1):
        out[seq[i : i + k - 1]].append(seq[i + k - 1])
    end = seq[n - k + 1 :]
    # Wilson: a uniform spanning in-tree towards `end`; each vertex's tree edge is its last exit.
    # ponytail: pure-Python dict walk, ~1e6 bases/s (200 kb in 0.2 s); numpy codes only if Gb-scale.
    last: dict[str, str] = {}
    in_tree = {end}
    for root in sorted(out):
        u = root
        while u not in in_tree:
            last[u] = rng.choice(out[u])
            u = (u + last[u])[1:]
        u = root
        while u not in in_tree:
            in_tree.add(u)
            u = (u + last[u])[1:]
    for vertex, symbols in out.items():
        if vertex in last:
            symbols.remove(last[vertex])  # one copy: parallel edges are distinct
        symbols.sort()  # fixed order before the shuffle, so the stream is reproducible
        rng.shuffle(symbols)
        if vertex in last:
            symbols.append(last[vertex])
    walk = [seq[: k - 1]]
    vertex = seq[: k - 1]
    for _ in range(n - k + 1):
        symbol = out[vertex].pop(0)
        walk.append(symbol)
        vertex = (vertex + symbol)[1:]
    return "".join(walk)


def _rng(seed: int, k: int, decoy_id: str) -> random.Random:
    digest = hashlib.sha256(f"{seed}:{k}:{decoy_id}".encode()).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def _sha(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def make_decoys(
    panel: list[tuple[str, str]],
    unrelated: list[tuple[str, str]],
    *,
    k: int,
    seed: int,
    copies: int,
) -> list[dict[str, Any]]:
    """Decoy records (manifest columns plus ``sequence``). Ids are checked for collisions."""
    rows: list[dict[str, Any]] = []
    for cls, records in (("kmer_shuffle", panel), ("unrelated_host", unrelated)):
        for source_id, seq in records:
            if source_id.startswith(PREFIX):
                raise ValueError(f"source id {source_id!r} already starts with {PREFIX}")
            if cls == "kmer_shuffle":
                for rep in range(1, copies + 1):
                    did = f"{PREFIX}{source_id}_{rep}"
                    shuffled = klet_shuffle(seq, k, _rng(seed, k, did))
                    rows.append(_row(did, cls, source_id, seed, k, rep, shuffled, seq))
            else:
                rows.append(_row(f"{PREFIX}{source_id}", cls, source_id, "", "", "", seq, seq))
    ids = [r["decoy_id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate decoy ids (a source id repeats across the inputs)")
    return rows


def _row(did: str, cls: str, src: str, seed: Any, k: Any, rep: Any, seq: str, source: str):
    return {
        "decoy_id": did,
        "class": cls,
        "source_id": src,
        "seed": seed,
        "k": k,
        "rep": rep,
        "length": len(seq),
        "sha256": _sha(seq),
        "source_sha256": _sha(source),
        "sequence": seq,
    }


def cmd_make(args: argparse.Namespace) -> None:
    def load(paths: list[Path]) -> list[tuple[str, str]]:
        return [rec for p in paths for rec in _fasta_records(p)]

    panel, unrelated = load(args.panel_fasta), load(args.unrelated_host_fasta)
    if not panel and not unrelated:
        _die("no input sequences")
    try:
        rows = make_decoys(panel, unrelated, k=args.k, seed=args.seed, copies=args.copies)
    except ValueError as exc:
        _die(str(exc))
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    _write_fasta(out / "decoys.fa", [(r["decoy_id"], r["sequence"]) for r in rows])
    gtf = [
        _genome_as_transcript_gtf(f">{r['decoy_id']}\n{r['sequence']}\n", r["decoy_id"])
        for r in rows
    ]
    (out / "decoys.gtf").write_text("\n".join(gtf) + "\n")
    with open(out / "decoys_manifest.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, MANIFEST_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} decoys to {out} (fasta sha256 {sha256_file(out / 'decoys.fa')[:12]})")


# ---------------------------------------------------------------- scoring


def empirical_p(score: float, decoy_scores_sorted: list[float]) -> tuple[int, float]:
    """(b, (b + 1) / (m + 1)) with b = decoys scoring >= ``score``."""
    m = len(decoy_scores_sorted)
    b = m - bisect_left(decoy_scores_sorted, score)
    return b, (b + 1) / (m + 1)


def td_qvalues(targets: list[float], decoys: list[float]) -> list[float]:
    """Stratified target-decoy q (one stratum): see module docstring. Order follows ``targets``."""
    if not targets:
        return []
    if not decoys:
        return [1.0] * len(targets)
    ts, ds = sorted(targets), sorted(decoys)
    ratio = len(ds) / len(ts)
    q_at: dict[float, float] = {}
    running = 1.0
    for s in sorted(set(ts)):  # ascending: FDR falls as the threshold rises
        d, t = len(ds) - bisect_left(ds, s), len(ts) - bisect_left(ts, s)
        running = min(running, 1.0, (d + 1) / (ratio * t))
        q_at[s] = running
    return [q_at[s] for s in targets]


def read_calls(path: Path) -> dict[str, list[dict[str, Any]]]:
    """{library_id: rows}; validates columns, booleans, finite scores and unique (library, id)."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(CALL_COLUMNS) - set(reader.fieldnames or ())
        if missing:
            _die(f"{path}: missing columns {sorted(missing)}")
        libs: dict[str, list[dict[str, Any]]] = defaultdict(list)
        seen: set[tuple[str, str]] = set()
        for n, rec in enumerate(reader, start=2):
            flag = rec["is_decoy"].strip().lower()
            if flag not in TRUE | FALSE:
                _die(f"{path}:{n}: is_decoy {rec['is_decoy']!r} is not a boolean")
            try:
                score = float(rec["score"])
            except ValueError:
                score = math.nan
            if not math.isfinite(score):
                _die(f"{path}:{n}: score {rec['score']!r} is not a finite number")
            if not rec["stratum"].strip():
                _die(f"{path}:{n}: empty stratum (use one constant value for no stratification)")
            key = (rec["library_id"], rec["target_id"])
            if key in seen:
                _die(f"{path}:{n}: duplicate (library_id, target_id) {key}")
            seen.add(key)
            libs[rec["library_id"]].append(
                {
                    "target_id": rec["target_id"],
                    "is_decoy": flag in TRUE,
                    "stratum": rec["stratum"].strip(),
                    "score": score,
                }
            )
    return libs


def score_library(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Target rows of one library with n_decoys, b, p_empirical, q_td, q_bh."""
    decoys: dict[str, list[float]] = defaultdict(list)
    targets: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        if r["is_decoy"]:
            decoys[r["stratum"]].append(r["score"])
        else:
            targets[r["stratum"]].append(r)
    out: list[dict[str, Any]] = []
    for stratum, trows in targets.items():
        dsorted = sorted(decoys.get(stratum, []))
        qs = td_qvalues([t["score"] for t in trows], dsorted)
        for t, q in zip(trows, qs, strict=True):
            b, p = empirical_p(t["score"], dsorted)
            out.append({**t, "n_decoys": len(dsorted), "b": b, "p_empirical": p, "q_td": q})
    if out:
        q_bh = _bh_adjust([r["p_empirical"] for r in out])
        for r, q in zip(out, q_bh, strict=True):
            r["q_bh"] = float(q)
    return out


def expected_false_calls(
    rows: list[dict[str, Any]],
    threshold: float,
    *,
    rng: np.random.Generator,
    n_draws: int,
    level: float,
) -> dict[str, Any]:
    """Expected false calls at ``score >= threshold`` and a Jeffreys-posterior interval."""
    strata: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])  # n_targets, n_decoys, decoy hits
    called = 0
    for r in rows:
        hit = r["score"] >= threshold
        s = strata[r["stratum"]]
        if r["is_decoy"]:
            s[1] += 1
            s[2] += hit
        else:
            s[0] += 1
            called += hit
    expected, draws, excluded = 0.0, np.zeros(n_draws), 0
    for n_t, n_d, d in strata.values():
        if not n_t:
            continue
        if not n_d:
            excluded += n_t
            continue
        expected += n_t * d / n_d
        draws += n_t * rng.beta(d + 0.5, n_d - d + 0.5, size=n_draws)
    lo, hi = np.quantile(draws, [(1 - level) / 2, 1 - (1 - level) / 2])
    return {
        "n_targets": sum(s[0] for s in strata.values()),
        "n_decoys": sum(s[1] for s in strata.values()),
        "n_called": called,
        "expected_false_calls": expected,
        "ci_lo": float(lo),
        "ci_hi": float(hi),
        "n_targets_in_strata_without_decoys": excluded,
    }


def load_prereg(path: Path | None) -> tuple[float | None, tuple[float, float] | None]:
    """Pre-registered (target, [lo, hi]) from a protocol YAML; ``TBD*`` placeholders are unset."""
    if path is None:
        return None, None
    import yaml

    block = (yaml.safe_load(path.read_text()) or {}).get("false_discovery") or {}

    def num(v: Any) -> float | None:
        return None if v is None or isinstance(v, str) else float(v)

    target = num(block.get("expected_false_L2_calls_per_library"))
    iv = block.get("interval")
    interval = None
    if isinstance(iv, (list, tuple)):
        if len(iv) != 2 or num(iv[0]) is None or num(iv[1]) is None:
            _die(f"{path}: false_discovery.interval must be [lo, hi]")
        interval = (float(iv[0]), float(iv[1]))
    return target, interval


def cmd_fdr(args: argparse.Namespace) -> None:
    if not 0 < args.level < 1:
        _die("--level must be in (0, 1)")
    libs = read_calls(args.calls)
    if args.manifest:
        with open(args.manifest, newline="") as fh:
            want = {r["decoy_id"] for r in csv.DictReader(fh, delimiter="\t")}
        for lib, rows in libs.items():
            absent = want - {r["target_id"] for r in rows if r["is_decoy"]}
            if absent:
                _die(f"library {lib}: {len(absent)} manifest decoys missing (e.g. {min(absent)})")
    target, interval = load_prereg(args.protocol)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    t_cols = (
        "library_id",
        "target_id",
        "stratum",
        "score",
        "n_decoys",
        "b",
        "p_empirical",
        "q_td",
        "q_bh",
    )
    l_cols = (
        "library_id", "threshold", "n_targets", "n_decoys", "n_called", "expected_false_calls",
        "ci_lo", "ci_hi", "level", "n_targets_in_strata_without_decoys", "fdr_at_threshold",
        "prereg_target", "prereg_lo", "prereg_hi", "expected_within_prereg_interval",
        "observed_within_prereg_interval",
    )  # fmt: skip
    with (
        open(args.out_dir / "target_pvalues.tsv", "w", newline="") as tf,
        open(args.out_dir / "library_summary.tsv", "w", newline="") as lf,
    ):
        tw = csv.DictWriter(tf, t_cols, delimiter="\t", extrasaction="ignore")
        lw = csv.DictWriter(lf, l_cols, delimiter="\t")
        tw.writeheader()
        lw.writeheader()
        for lib in sorted(libs):
            rows = libs[lib]
            for r in score_library(rows):
                tw.writerow({"library_id": lib, **r})
            seed = int.from_bytes(hashlib.sha256(f"{args.seed}:{lib}".encode()).digest()[:8], "big")
            s = expected_false_calls(
                rows,
                args.threshold,
                rng=np.random.default_rng(seed),
                n_draws=args.n_draws,
                level=args.level,
            )
            s["fdr_at_threshold"] = (
                min(1.0, s["expected_false_calls"] / s["n_called"]) if s["n_called"] else "NA"
            )
            lo, hi = interval if interval else (None, None)
            s.update(
                library_id=lib, threshold=args.threshold, level=args.level,
                prereg_target="NA" if target is None else target,
                prereg_lo="NA" if lo is None else lo, prereg_hi="NA" if hi is None else hi,
                expected_within_prereg_interval=(
                    "NA" if lo is None else lo <= s["expected_false_calls"] <= hi
                ),
                observed_within_prereg_interval="NA" if lo is None else lo <= s["n_called"] <= hi,
            )  # fmt: skip
            lw.writerow(s)
    print(f"scored {len(libs)} libraries -> {args.out_dir}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    mk = sub.add_parser("make", help="write decoy FASTA, GTF and manifest")
    mk.add_argument("--panel-fasta", type=Path, nargs="+", default=[], help="panel genome FASTA(s)")
    mk.add_argument("--unrelated-host-fasta", type=Path, nargs="+", default=[],
                    help="FASTA(s) of unrelated-host viral genomes, registered unshuffled")  # fmt: skip
    mk.add_argument("--out-dir", type=Path, required=True)
    mk.add_argument("--seed", type=int, required=True)
    mk.add_argument("--k", type=int, default=DEFAULT_K, help="k-let size (default 6)")
    mk.add_argument("--copies", type=int, default=1, help="shuffles per panel genome")
    mk.set_defaults(func=cmd_make)
    fd = sub.add_parser("fdr", help="empirical p, target-decoy q, expected false calls")
    fd.add_argument("--calls", type=Path, required=True, help="library_id target_id is_decoy stratum score")  # fmt: skip
    fd.add_argument("--out-dir", type=Path, required=True)
    fd.add_argument("--threshold", type=float, required=True, help="score >= threshold is a call")
    fd.add_argument("--seed", type=int, required=True, help="interval Monte Carlo seed")
    fd.add_argument("--n-draws", type=int, default=10000)
    fd.add_argument("--level", type=float, default=0.95, help="interval level")
    fd.add_argument("--manifest", type=Path, help="decoys_manifest.tsv: require every decoy per library")  # fmt: skip
    fd.add_argument("--protocol", type=Path, help="protocol YAML with the pre-registered false_discovery target and interval")  # fmt: skip
    fd.set_defaults(func=cmd_fdr)
    args = parser.parse_args(argv)
    if args.command == "make" and (args.k < 1 or args.copies < 1):
        parser.error("--k and --copies must be >= 1")
    args.func(args)


if __name__ == "__main__":
    main()

"""Canonical k-mer sharing between genomes: the curation-time twin machinery, importable from the package.

Moved out of `scripts/panel_kmer_sharing.py` (AUD3-05) so the build can reuse it; the script re-exports these names.
Needs numpy. No `__future__` import here (kept safe for Snakemake-adjacent callers).

`KMER_TWIN_BELOW`: a genome with under this fraction of its distinct canonical k-mers outside the other genomes
is a twin (`kmer_twin_of:` in `panel_candidates.py`). `KMER_TWIN_K` is the k the screen uses. Neither is a
homology threshold: calibrating those is PLAN REF-08.
"""

from typing import Any

import numpy as np

KMER_TWIN_BELOW = 0.05
KMER_TWIN_K = 31

_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i  # lower case


def canonical_kmers(seq: bytes, k: int) -> np.ndarray:
    """Sorted distinct canonical k-mers (2 bits per base, k <= 31) of one sequence."""
    code = _CODE[np.frombuffer(seq, dtype=np.uint8)]
    n = len(code) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64)
    fwd = np.zeros(n, dtype=np.uint64)
    rev = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        window = code[j : j + n]
        fwd = (fwd << np.uint64(2)) | (window & np.uint8(3)).astype(np.uint64)
        rev |= (np.uint64(3) - (window & np.uint8(3)).astype(np.uint64)) << np.uint64(2 * j)
    bad = np.concatenate(([0], np.cumsum(code > 3, dtype=np.int64)))
    ok = (bad[k : k + n] - bad[:n]) == 0
    return np.unique(np.minimum(fwd, rev)[ok])


def sharing(
    kmers: list[np.ndarray],
    group: list[int],
    is_panel: list[bool],
    targets: list[int],
    partner_below: float = 0.0,
    top_n: int = 3,
) -> dict[int, dict[str, Any]]:
    """Per target genome: the four fractions described in the module docstring.

    A target whose `frac_not_in_panel` is under `partner_below` also gets `partners`: up to `top_n`
    (genome index, shared k-mers) pairs, most shared first, ties to the lower index (panel genomes first).
    """
    lens = np.array([len(a) for a in kmers], dtype=np.int64)
    all_k = np.concatenate(kmers)
    gid = np.repeat(np.arange(len(kmers), dtype=np.int32), lens)
    grp = np.asarray(group, dtype=np.int32)[gid]
    order = np.lexsort((grp, all_k))
    k_s, g_s, id_s = all_k[order], grp[order], gid[order]
    del all_k, gid, grp, order
    new_k = np.concatenate(([True], k_s[1:] != k_s[:-1]))
    new_run = new_k | np.concatenate(([True], g_s[1:] != g_s[:-1]))
    kid = np.cumsum(new_k, dtype=np.int32) - 1
    rid = np.cumsum(new_run, dtype=np.int32) - 1
    n_in_kmer = np.bincount(kid)
    n_groups = np.bincount(kid, weights=new_run)
    run_size = np.bincount(rid)
    panel_row = np.asarray(is_panel)[id_s]
    n_panel = np.bincount(kid, weights=panel_row)
    flags = {
        "frac_unique": n_in_kmer[kid] == 1,
        "frac_not_in_panel": (n_panel[kid] - panel_row) == 0,
        "frac_sibling": run_size[rid] > 1,
        "frac_other_group": n_groups[kid] > 1,
    }
    out: dict[int, dict[str, Any]] = {t: {} for t in targets}
    for name, flag in flags.items():
        per_genome = np.bincount(id_s, weights=flag, minlength=len(kmers))
        for t in targets:
            out[t][name] = per_genome[t] / lens[t] if lens[t] else float("nan")
    for t in targets:
        if not out[t]["frac_not_in_panel"] < partner_below:
            continue
        in_target = np.zeros(len(n_in_kmer), dtype=bool)
        in_target[kid[id_s == t]] = True
        shared = np.bincount(id_s[in_target[kid]], minlength=len(kmers))
        shared[t] = 0
        top = np.argsort(-shared, kind="stable")[:top_n]
        out[t]["partners"] = [(int(i), int(shared[i])) for i in top if shared[i]]
    return out


def _sibling_group(name: str, group_of: dict[str, str]) -> str:
    return group_of.get(name) or group_of.get(name.split(".")[0], "")


def near_duplicate_pairs(
    records: list[tuple[str, str]],
    *,
    viral_ids: set[str] | None = None,
    group_of: dict[str, str] | None = None,
    k: int = KMER_TWIN_K,
) -> list[tuple[str, str, float]]:
    """Records that are twins of another record, on both strands, as ``(record, partner, shared_fraction)``.

    A record is a twin when under `KMER_TWIN_BELOW` of its distinct canonical k-mers lie outside every other record
    (the curation-time `kmer_twin_of` rule); the partner is the record sharing most with it. With *viral_ids* only
    viral records are tested (against host and viral partners), so host-host paralogs and isoforms are ignored and
    host records sharing no k-mer with any viral record are never kept. A partner in the same *group_of* group
    (declared siblings, `identifier or bare accession -> group`) is not reported.
    """
    # ponytail: only the top-ranked partner is judged, so a record covered jointly by a sibling and a third
    # unrelated record is not reported; look at ranks 2-3 if that case shows up.
    group_of = group_of or {}
    kmers = {
        name: canonical_kmers(seq.encode(), k)
        for name, seq in records
        if viral_ids is None or name in viral_ids
    }
    targets = list(kmers)
    if viral_ids is not None and kmers:
        viral_union = np.unique(np.concatenate(list(kmers.values())))
        for name, seq in records:
            if name not in viral_ids:
                host = canonical_kmers(seq.encode(), k)
                if np.isin(host, viral_union, assume_unique=True).any():
                    kmers[name] = host
    if not sum(len(a) for a in kmers.values()):
        return []
    names = list(kmers)
    index = {n: i for i, n in enumerate(names)}
    stats = sharing(
        [kmers[n] for n in names],
        list(range(len(names))),
        [True] * len(names),
        [index[t] for t in targets],
        partner_below=KMER_TWIN_BELOW,
    )
    pairs: dict[frozenset[str], tuple[str, str, float]] = {}
    for t in targets:
        top = stats[index[t]].get("partners")
        if not top:
            continue
        partner, shared = names[top[0][0]], top[0][1]
        group = _sibling_group(t, group_of)
        if group and group == _sibling_group(partner, group_of):
            continue
        pairs.setdefault(frozenset((t, partner)), (t, partner, shared / len(kmers[t])))
    return list(pairs.values())

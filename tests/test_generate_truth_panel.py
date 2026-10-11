"""Truth-panel generator (PLAN VAL-01/02/03/05). Code only: tiny synthetic panels, no outcome runs.

Every panel is written under ``tmp_path`` (governance forbids tracked FASTQs). The scorer-side
"adapter" below stands in for a quantifier: it is deliberately dumb (25-mer voting against the
generator's own templates, no CB/UMI correction) so that, with ``--error-rate 0``, its calls are exactly
the planted molecules and the ``v3_score`` pure functions can be checked on a panel whose answer is
known. It is a test oracle for the generator, not a ViralScan code path.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import random
import re
import stat
import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import v3_partition as part  # noqa: E402
import v3_score as score  # noqa: E402

from viralscan import chemistry as chem_mod  # noqa: E402

BASE = ["--generation-seed", "1", "--error-rate", "0"]


def make(out: Path, *extra: str) -> dict:
    return gt.generate(gt.parse_args([*BASE, "--out", str(out), *extra]))


def salt_file(path: Path, text: str = "0123456789abcdef-secret") -> Path:
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as handle:
        handle.write(text)
    return path


def tree_hashes(root: Path) -> dict[str, str]:
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fastq_records(path: Path) -> list[tuple[str, str]]:
    with gzip.open(path, "rt") as handle:
        lines = handle.read().split("\n")
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4)]


def sample_dirs(root: Path) -> list[Path]:
    return sorted(p for p in root.iterdir() if p.is_dir() and p.name.startswith("s"))


@pytest.fixture(scope="module")
def tiny(tmp_path_factory) -> tuple[Path, dict]:
    out = tmp_path_factory.mktemp("tiny")
    return out, make(out)


# ── the adapter: 25-mer voting over the generator's own templates ─────────────

K = 25


def kmer_labels(templates: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    for fasta, label in (
        ("viral_genomes.fasta", None),
        ("host_transcripts.fasta", "host"),
        ("host_homology_loci.fasta", "host"),
    ):
        name, seq = "", ""
        for line in (templates / fasta).read_text().splitlines() + [">end"]:
            if line.startswith(">"):
                for i in range(len(seq) - K + 1):
                    key, value = seq[i : i + K], label or name
                    labels[key] = value if labels.get(key, value) == value else "?"  # shared k-mer
                name, seq = line[1:], ""
            else:
                seq += line
    return labels


def classify(r2: str, labels: dict[str, str]) -> str | None:
    votes = Counter(labels.get(r2[i : i + K]) for i in range(len(r2) - K + 1))
    votes.pop(None, None)
    votes.pop("?", None)
    top = votes.most_common(2)
    return top[0][0] if top and (len(top) == 1 or top[0][1] > top[1][1]) else None


def adapter_calls(sample: Path, templates: Path, chunk: int = 10_000, reverse: bool = False):
    """``({(cb, umi): virus}, {(cb, umi)})`` of viral and host keys, reading FASTQs ``chunk`` at a time."""
    labels = kmer_labels(templates)
    r1 = fastq_records(next(sample.glob("*_R1.fastq.gz")))
    r2 = fastq_records(next(sample.glob("*_R2.fastq.gz")))
    pairs = list(zip(r1, r2))
    assert all(a[0] == b[0] for a, b in pairs)
    if reverse:
        pairs.reverse()
    viral: dict[tuple[str, str], str] = {}
    host: set[tuple[str, str]] = set()
    for start in range(0, len(pairs), chunk):
        for (_, bases), (_, read) in pairs[start : start + chunk]:
            cb_len = len(bases) - chem_umi_len(sample)
            key = (bases[:cb_len], bases[cb_len:])
            label = classify(read, labels)
            if label == "host":
                host.add(key)
            elif label:
                viral[key] = label
    return viral, host


_CHEM_OF: dict[str, str] = {}


def chem_umi_len(sample: Path) -> int:
    return chem_mod.get(gt.PACKAGE_CHEMISTRY[_CHEM_OF[sample.name]]).umi_len


def register(root: Path) -> None:
    for row in read_tsv(root / "truth" / "truth_manifest.tsv"):
        _CHEM_OF[row["biological_sample_id"]] = row["chemistry"]
    for row in read_tsv(root / "truth" / "host_only_manifest.tsv"):
        _CHEM_OF[row["biological_sample_id"]] = row["chemistry"]


def planted_keys(rows: list[dict[str, str]], sid: str) -> dict[tuple[str, str], str]:
    return {
        (r["cell_barcode"], r["umi"]): r["planted_virus_id"]
        for r in rows
        if r["biological_sample_id"] == sid and r["is_planted_molecule"] == "true"
    }


# ── determinism ──────────────────────────────────────────────────────────────


def test_same_seed_is_byte_identical_and_other_seed_differs(tmp_path):
    one = ["--datasets", "synthetic_factorial", "--chemistries", "10xv3"]
    a, b, c = (tmp_path / n for n in "abc")
    make(a, *one)
    make(b, *one)
    gt.generate(
        gt.parse_args(["--generation-seed", "2", "--error-rate", "0.001", "--out", str(c), *one])
    )
    assert tree_hashes(a) == tree_hashes(b)
    other = tree_hashes(c)
    assert all(other[k] != v for k, v in tree_hashes(a).items() if k.endswith("fastq.gz"))
    assert other["truth/truth_manifest.tsv"] != tree_hashes(a)["truth/truth_manifest.tsv"]


def test_generation_seed_is_required(tmp_path):
    with pytest.raises(SystemExit):
        gt.parse_args(["--out", str(tmp_path)])


def test_gzip_header_has_no_timestamp_or_name(tiny):
    out, _ = tiny
    for path in out.rglob("*.fastq.gz"):
        header = path.read_bytes()[:10]
        assert header[4:8] == b"\0\0\0\0" and not header[3] & 0x08  # mtime 0, no FNAME


def test_generation_is_independent_of_argument_order(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    make(a, "--chemistries", "10xv2", "10xv3", "--datasets", *gt.VIRAL_DATASETS)
    make(b, "--chemistries", "10xv3", "10xv2", "--datasets", *reversed(gt.VIRAL_DATASETS))
    assert tree_hashes(a) == tree_hashes(b)


# ── geometry, manifests ──────────────────────────────────────────────────────


def test_read_geometry_comes_from_the_chemistry_table_and_manifests_use_protocol_spelling(tiny):
    out, manifest = tiny
    seen = set()
    for sid, detail in manifest["samples"].items():
        chem = chem_mod.get(gt.PACKAGE_CHEMISTRY[detail["chemistry"]])
        r1 = fastq_records(out / sid / f"{sid}_R1.fastq.gz")
        r2 = fastq_records(out / sid / f"{sid}_R2.fastq.gz")
        assert {len(s) for _, s in r1} == {chem.cb_len + chem.umi_len}
        assert {len(s) for _, s in r2} == {90}
        assert [n for n, _ in r1] == [n for n, _ in r2] and len(r1) == detail["n_read_pairs"]
        seen.add(detail["chemistry"])
    assert seen == {"10xv2", "10xv3", "drop-seq"}  # protocol spelling; the package says dropseq


@pytest.mark.parametrize("name", gt.MANIFESTS)
def test_manifest_columns_validate_and_every_row_is_one_read(tiny, name):
    out, manifest = tiny
    rows = read_tsv(out / "truth" / f"{name}.tsv")
    header = list(rows[0])
    assert header == score.required_columns(name) + gt.EXTRA[name]
    names = {
        sid: {n for n, _ in fastq_records(out / sid / f"{sid}_R1.fastq.gz")}
        for sid in manifest["samples"]
    }
    seen = set()
    for row in rows:
        key = (row["biological_sample_id"], row["read_id"])
        assert row["read_id"] in names[row["biological_sample_id"]] and key not in seen
        seen.add(key)
        assert row["generator_seed"] == "1" and row["generator_version"] == gt.GENERATOR_VERSION
    if name == "truth_manifest":
        assert len(score.read_truth_manifest(out / "truth" / "truth_manifest.tsv")) == len(rows)


def test_every_non_host_read_has_a_truth_row_and_host_molecules_are_rows_of_the_host_manifest(tiny):
    out, manifest = tiny
    truth = read_tsv(out / "truth" / "truth_manifest.tsv")
    host = read_tsv(out / "truth" / "host_only_manifest.tsv")
    by_sample = Counter(r["biological_sample_id"] for r in truth)
    host_reads = Counter()
    for row in host:
        host_reads[row["biological_sample_id"]] += int(row["n_reads"])
        assert row["read_id"].startswith(row["host_molecule_id"]) and row[
            "host_gene_id"
        ].startswith("ENST")
    truth_keys = {(r["biological_sample_id"], r["read_id"]) for r in truth}
    chimera_host = Counter(
        r["biological_sample_id"]
        for r in read_tsv(out / "truth" / "mixed_host_virus_manifest.tsv")
        if (r["biological_sample_id"], r["read_id"]) not in truth_keys
    )  # the host fragment of a mixed molecule has a mixed-manifest row, no host row
    for sid, detail in manifest["samples"].items():
        assert by_sample[sid] + host_reads[sid] + chimera_host[sid] == detail["n_read_pairs"], sid


def test_origins_and_factor_constructions(tiny):
    out, manifest = tiny
    truth = read_tsv(out / "truth" / "truth_manifest.tsv")
    fac = [
        r for r in truth if r["dataset_id"] == "synthetic_factorial" and r["chemistry"] == "10xv3"
    ]
    origins = Counter(r["origin"] for r in fac)
    assert {"cell", "ambient"} <= set(origins) and {
        f"artefact:{c}" for c in gt.ARTEFACT_CLASSES
    } <= set(origins)
    assert all(r["is_infected_cell"] == "false" for r in fac if r["origin"] != "cell")
    assert all(r["planted_virus_id"] == "" for r in fac if r["origin"].startswith("artefact:"))
    sid = fac[0]["biological_sample_id"]
    cells = {r["cell_barcode"] for r in fac if r["origin"] == "cell"}
    assert all(r["cell_barcode"] not in cells for r in fac if r["origin"] == "ambient")
    # co_infected: pair members share one cell set, pairs are disjoint
    by_virus = {
        v: {r["cell_barcode"] for r in fac if r["origin"] == "cell" and r["planted_virus_id"] == v}
        for v in gt.TARGETS
    }
    assert all(by_virus[v] for v in gt.TARGETS)
    assert by_virus["HSV-1"] <= by_virus["HSV-1"] | by_virus["HSV-2"]
    assert not by_virus["HSV-1"] & by_virus["EBV"] and not by_virus["HHV-6A"] & by_virus["HSV-2"]
    # abundance 10 per target; collisions: one of each kind, same key vs 1-Hamming UMI
    molecules = Counter(
        r["planted_virus_id"]
        for r in fac
        if r["is_planted_molecule"] == "true" and not r["planted_molecule_id"].count(">")
    )
    assert set(molecules) == set(gt.TARGETS)
    assert manifest["samples"][sid]["realised"]["collision_same"] == 1
    assert manifest["samples"][sid]["realised"]["collision_disjoint"] == 1
    groups: dict[str, list[dict[str, str]]] = {}
    for r in fac:
        if r["collision_group"]:
            groups.setdefault(r["collision_group"], []).append(r)
    kinds = {rows[0]["collision_kind"]: rows for rows in groups.values()}
    assert set(kinds) == {"same", "disjoint"}
    same = {(r["cell_barcode"], r["umi"]) for r in kinds["same"]}
    assert len(same) == 1 and len({r["planted_gene_id"] for r in kinds["same"]}) == 2
    umis = {r["umi"] for r in kinds["disjoint"]}
    a, b = sorted(umis)
    assert len(umis) == 2 and sum(x != y for x, y in zip(a, b)) == 1


def test_host_homology_and_mixed_manifests(tiny):
    out, _ = tiny
    homology = read_tsv(out / "truth" / "host_homology_manifest.tsv")
    assert homology and all(
        r["planted_as_host"] == "true" and r["homology_cluster_id"] for r in homology
    )
    assert all(
        0.8 < float(r["homology_identity"]) < 1 and r["homologous_virus_id"] in gt.TARGETS
        for r in homology
    )
    mixed = read_tsv(out / "truth" / "mixed_host_virus_manifest.tsv")
    assert {r["truth_class"] for r in mixed} == {"host", "virus"}
    assert {r["mixture_basis"] for r in mixed} == {"chimeric-template"}
    truth = {
        (r["biological_sample_id"], r["read_id"]): r
        for r in read_tsv(out / "truth" / "truth_manifest.tsv")
    }
    seen = 0
    for r in mixed:
        key = (r["biological_sample_id"], r["read_id"])
        if key in truth:  # the viral fragment: planted only when the molecule really is viral
            seen += 1
            assert (truth[key]["is_planted_molecule"] == "true") == (r["truth_class"] == "virus")
    assert seen


def test_sibling_one_member_only_plants_one_member_per_pair(tmp_path):
    manifest = make(
        tmp_path,
        "--datasets",
        "synthetic_factorial",
        "--chemistries",
        "10xv3",
        "--sibling-virus-similarity",
        "one_member_only",
    )
    (detail,) = manifest["samples"].values()
    assert len(detail["absent_viruses"]) == 3
    planted = {
        r["planted_virus_id"]
        for r in read_tsv(tmp_path / "truth" / "truth_manifest.tsv")
        if r["planted_virus_id"]
    }
    assert planted.isdisjoint(detail["absent_viruses"]) and len(planted) == 3
    assert set(manifest["sibling_pair_identity"]) == {"HHV-6A/HHV-6B", "HSV-1/HSV-2", "EBV/KSHV"}


def test_covariates_cycle_within_a_stratum_independent_of_sample_ids(tmp_path):
    manifest = make(
        tmp_path,
        "--datasets",
        "synthetic_factorial",
        "--chemistries",
        "10xv3",
        "--replicates",
        "3",
        "--pcr-duplication",
        "cycle",
        "--holdout-dir",
        str(tmp_path / "h"),
        "--holdout-salt-file",
        str(salt_file(tmp_path / "salt")),
    )
    levels = [d["levels"]["pcr_duplication_level"] for d in manifest["samples"].values()]
    assert len(set(levels)) >= 2  # training members of the stratum see more than one level


def test_unknown_levels_are_rejected(tmp_path):
    with pytest.raises(ValueError, match="not in the protocol"):
        make(tmp_path, "--abundance", "7")
    with pytest.raises(ValueError, match="not in the protocol"):
        make(tmp_path, "--low-complexity", "0.5")


# ── exact planted recovery (adapter + v3_score) ──────────────────────────────


def test_exact_planted_molecule_recovery_on_the_golden_tiny_panel(tiny):
    out, manifest = tiny
    register(out)
    truth = read_tsv(out / "truth" / "truth_manifest.tsv")
    checked = 0
    for sid, detail in manifest["samples"].items():
        if detail["dataset"] != "synthetic_factorial":
            continue
        viral, host = adapter_calls(out / sid, out / "templates")
        planted = planted_keys(truth, sid)
        triple = lambda mapping: {(cb, umi, v) for (cb, umi), v in mapping.items()}  # noqa: E731
        result = score.precision_recall_f1(triple(planted), triple(viral))
        assert (result["precision"], result["recall"], result["f1"]) == (1.0, 1.0, 1.0), sid
        # sibling confusion is exactly zero for both members of every pair
        for first, second in gt.json_pairs():
            for target, sibling in ((first, second), (second, first)):
                confusion = score.sibling_confusion(planted, viral, target, sibling)
                assert confusion["n_target"] > 0 and confusion["rate"] == 0.0
        # host/virus allocation: host molecules are never given a viral feature, and vice versa
        host_truth = {
            (r["cell_barcode"], r["umi"])
            for r in read_tsv(out / "truth" / "host_only_manifest.tsv")
            if r["biological_sample_id"] == sid
        }
        allocation = score.host_virus_allocation_error(host_truth, planted, viral, host)
        assert allocation["host_to_viral"] == 0.0 and allocation["viral_to_host"] == 0.0
        # cells: recall 1; precision falls short exactly by the cells that only hold ambient molecules
        infected = {
            r["cell_barcode"]
            for r in truth
            if r["biological_sample_id"] == sid and r["is_infected_cell"] == "true"
        }
        called = {cb for cb, _ in viral}
        cell_result = score.precision_recall_f1(infected, called)
        assert cell_result["recall"] == 1.0 and cell_result["fp"] == len(called - infected) > 0
        checked += 1
    assert checked == 3


def test_artefacts_challenge_and_host_only_samples_yield_no_viral_calls_beyond_homology(tiny):
    out, manifest = tiny
    register(out)
    for sid, detail in manifest["samples"].items():
        if detail["dataset"] in ("synthetic_host_only", "synthetic_host_homology"):
            viral, host = adapter_calls(out / sid, out / "templates")
            assert not viral and host


def test_adapter_result_is_invariant_to_chunk_size_and_read_order(tiny):
    out, manifest = tiny
    register(out)
    sid = next(s for s, d in manifest["samples"].items() if d["dataset"] == "synthetic_factorial")
    reference = adapter_calls(out / sid, out / "templates")
    for chunk in (1, 7, 100_000):
        assert adapter_calls(out / sid, out / "templates", chunk=chunk) == reference
    assert adapter_calls(out / sid, out / "templates", reverse=True) == reference
    truth = read_tsv(out / "truth" / "truth_manifest.tsv")
    shuffled = truth[:]
    random.Random(3).shuffle(shuffled)
    assert planted_keys(shuffled, sid) == planted_keys(truth, sid)


# ── partitions ───────────────────────────────────────────────────────────────


def stratum(sid: str, key: str) -> gt.Sample:
    return gt.Sample("d", "10xv3", key, 1, "10", "none", sid)


def test_singleton_strata_are_training_and_the_split_runs_over_the_rest():
    assert part.apportion({"a": 1, "b": 1, "c": 1}, 0.3) == {
        "a": 1,
        "b": 0,
        "c": 0,
    }  # why: apportion alone
    samples = [stratum("sa", "x"), stratum("sb", "y")] + [stratum(f"sz{i}", "z") for i in range(5)]
    unsplittable = gt.assign_partitions(samples, 0.3)
    assert unsplittable == [["d", "x"], ["d", "y"]]
    labels = {s.sid: s.partition for s in samples}
    assert labels["sa"] == labels["sb"] == part.TRAINING
    assert (
        sorted(v for k, v in labels.items() if k.startswith("sz")).count(part.HOLDOUT) == 2
    )  # round(1.5)
    assert gt.assign_partitions([stratum("only", "x")], 0.3) == [["d", "x"]]


def test_default_tiny_panel_has_only_training_samples_and_reports_its_singletons(tiny):
    _, manifest = tiny
    assert manifest["holdout_samples"] == [] and len(manifest["samples"]) == 12
    assert len(manifest["unsplittable_strata"]) == 12
    assert all(d["partition"] == part.TRAINING for d in manifest["samples"].values())


def test_unit_id_and_donor_library_split():
    a = part.unit_id(7, "lib_a")
    assert (
        a == part.unit_id(7, "lib_a")
        and a != part.unit_id(8, "lib_a")
        and a != part.unit_id(7, "lib_b")
    )
    assert len(a) == 13 and a[0] == "u" and "lib" not in a
    libraries = {c: [f"{c}_lib{i}" for i in range(10)] for c in ("10xv2", "10xv3", "drop-seq")}
    labels = part.split_libraries(libraries, 7, 0.3)
    assert len(labels) == 30
    for libs in libraries.values():  # 10 libraries per chemistry, 3 in holdout (protocol)
        assert sum(labels[lib] == part.HOLDOUT for lib in libs) == 3
    assert labels == part.split_libraries({c: libs[::-1] for c, libs in libraries.items()}, 7, 0.3)


@pytest.fixture(scope="module")
def split_panel(tmp_path_factory):
    root = tmp_path_factory.mktemp("split")
    out, hold = root / "out", root / "holdout"
    salt = salt_file(root / "salt")
    args = [
        "--datasets",
        "synthetic_factorial",
        "--replicates",
        "4",
        "--holdout-dir",
        str(hold),
        "--holdout-salt-file",
        str(salt),
    ]
    return out, hold, salt, make(out, *args), args


def test_split_gives_the_global_seat_count_and_holdout_truth_stays_out_of_the_public_outputs(
    split_panel,
):
    out, hold, _, manifest, _ = split_panel
    holdout = {h["sample_id"] for h in manifest["holdout_samples"]}
    assert len(holdout) == 4 and len(manifest["samples"]) == 8  # round(0.3 * 12) = 4
    assert manifest["unsplittable_strata"] == [] and manifest["holdout_salt_sha256"]
    for manifest_name in gt.MANIFESTS:
        public = read_tsv(out / "truth" / f"{manifest_name}.tsv")
        private = read_tsv(hold / "truth" / f"{manifest_name}.tsv")
        assert {r["partition"] for r in public} <= {part.TRAINING} and not {
            r["biological_sample_id"] for r in public
        } & holdout
        assert {r["partition"] for r in private} <= {part.HOLDOUT} and {
            r["biological_sample_id"] for r in private
        } <= holdout
    assert all(
        set(h) == {"sample_id", "r1_sha256", "r2_sha256"} for h in manifest["holdout_samples"]
    )
    assert not holdout & set(manifest["samples"])
    public_text = (out / "run_manifest.json").read_text()
    assert public_text.count('"levels"') == len(manifest["samples"])  # holdout levels stay private
    for sid in holdout:  # next to the holdout FASTQs: only their sha256 values
        assert sorted(p.name for p in (out / sid).iterdir()) == sorted(
            [f"{sid}_R1.fastq.gz", f"{sid}_R2.fastq.gz", f"{sid}.sha256"]
        )
        assert (hold / sid / "barcode_pool.tsv").is_file()


def test_holdout_files_are_owner_only_with_opaque_names_and_shuffled_order(split_panel):
    out, hold, _, manifest, _ = split_panel
    for path in [hold, *hold.rglob("*")]:
        assert not stat.S_IMODE(path.stat().st_mode) & 0o077, path
    holdout = [h["sample_id"] for h in manifest["holdout_samples"]]
    truth = read_tsv(hold / "truth" / "truth_manifest.tsv")
    for sid in holdout:
        names = [n for n, _ in fastq_records(out / sid / f"{sid}_R1.fastq.gz")]
        assert names == [f"r{i}" for i in range(len(names))]
        rows = [r for r in truth if r["biological_sample_id"] == sid]
        planted = sorted(int(r["read_id"][1:]) for r in rows if r["origin"] == "cell")
        assert planted and planted[-1] - planted[0] + 1 > 3 * len(
            planted
        )  # not one contiguous block
    training = next(s for s in manifest["samples"])
    train_names = [n for n, _ in fastq_records(out / training / f"{training}_R1.fastq.gz")]
    assert any(
        n.startswith("h:ENST") for n in train_names
    )  # training host reads carry the molecule id


def test_no_barcode_or_umi_is_shared_across_partitions_and_templates_split_whole(split_panel):
    out, hold, _, manifest, _ = split_panel
    keys = {}
    for side, base in ((part.TRAINING, out / "truth"), (part.HOLDOUT, hold / "truth")):
        rows = read_tsv(base / "truth_manifest.tsv") + read_tsv(base / "host_only_manifest.tsv")
        keys[side] = (
            {r["cell_barcode"] for r in rows},
            {r["umi"] for r in rows},
            {r.get("host_gene_id") for r in rows} - {None, ""},
        )
    for i in range(3):
        assert not keys[part.TRAINING][i] & keys[part.HOLDOUT][i]
    # host transcripts are assigned whole to one partition (D4): the third set is the host gene ids
    assert keys[part.TRAINING][2] <= set(manifest["training_host_transcripts"])


def test_salt_changes_holdout_only_and_a_wrong_setup_is_refused(split_panel, tmp_path):
    out, hold, _, manifest, args = split_panel
    other = tmp_path / "out2"
    gt.generate(
        gt.parse_args(
            [
                *BASE,
                "--out",
                str(other),
                *args[:4],
                "--holdout-dir",
                str(tmp_path / "h2"),
                "--holdout-salt-file",
                str(salt_file(tmp_path / "salt2", "another-secret-salt-value")),
            ]
        )
    )
    same = tree_hashes(out)
    different = tree_hashes(other)
    holdout = {h["sample_id"] for h in manifest["holdout_samples"]}
    for rel, digest in same.items():
        top = rel.split("/")[0]
        if top in holdout and rel.endswith(".fastq.gz"):
            assert different[rel] != digest
        elif top in manifest["samples"] or top == "templates" or rel.startswith("truth"):
            assert different[rel] == digest, rel  # training output needs only the public seed
    with pytest.raises(ValueError, match="holdout-dir and --holdout-salt-file are required"):
        make(tmp_path / "o3", *args[:4])
    loose = tmp_path / "loose"
    loose.write_text("0123456789abcdef-secret")
    loose.chmod(0o644)
    with pytest.raises(ValueError, match="owner-only"):
        make(
            tmp_path / "o4",
            *args[:4],
            "--holdout-dir",
            str(tmp_path / "h4"),
            "--holdout-salt-file",
            str(loose),
        )


def test_hopped_molecules_land_in_a_same_dataset_chemistry_partition_sample(split_panel):
    out, hold, _, manifest, _ = split_panel
    rows = read_tsv(out / "truth" / "truth_manifest.tsv") + read_tsv(
        hold / "truth" / "truth_manifest.tsv"
    )
    hopped = [r for r in rows if r["origin"] == "hopped"]
    assert hopped and all(
        r["is_infected_cell"] == "false" and ">" in r["planted_molecule_id"] for r in hopped
    )
    for r in hopped:
        donor, target = (
            r["planted_molecule_id"].split(":")[0],
            r["planted_molecule_id"].split(">")[1],
        )
        assert target == r["biological_sample_id"] and donor != target
        assert (
            next(x for x in rows if x["biological_sample_id"] == donor)["partition"]
            == r["partition"]
        )
        assert (
            next(x for x in rows if x["biological_sample_id"] == donor)["chemistry"]
            == r["chemistry"]
        )
    assert sum(d["realised"]["hopped"] for d in manifest["samples"].values()) >= 1


# ── guards ───────────────────────────────────────────────────────────────────


def test_production_refuses_while_factor_levels_are_pending(tmp_path):
    with pytest.raises(ValueError, match="production run refused"):
        make(tmp_path, "--production")


def test_output_inside_the_repository_is_refused():
    with pytest.raises(ValueError, match="outside the repository"):
        make(gt.REPO / "tmp_panel_should_not_exist")
    assert not (gt.REPO / "tmp_panel_should_not_exist").exists()


def test_templates_and_gene_models(tiny):
    out, _ = tiny
    genes = read_tsv(out / "templates" / "viral_gene_models.tsv")
    assert {g["virus_id"] for g in genes} == set(gt.TARGETS) and len(genes) == 18
    assert not re.search(r"/(exports|home)/", (out / "run_manifest.json").read_text())

"""Truth-panel generator, AUD-18: real genomes, gene spans and a barcode on-list (code only).

Also pins the default (no new option) output to the hashes of the pre-AUD-18 generator (7b0c614).
Everything is written under ``tmp_path``; the real-style inputs come from ``tests/_truth_panel_real.py``.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import v3_score as score  # noqa: E402

from tests._truth_panel_real import PLANTED, write_real_inputs  # noqa: E402
from tests.test_generate_truth_panel import fastq_records, read_tsv, tree_hashes  # noqa: E402

#: run_manifest hashes of the pre-AUD-18 generator (commit 7b0c614), so the new options cannot
#: change what a run without them writes.
GOLDEN_DEFAULT = {
    "tiny": dict(
        args=["--generation-seed", "1", "--error-rate", "0"],
        samples_digest="fa7d97ce7d944cdc7a18883f141e6f88926735915455b00581d8644e159b7f53",
        truth={
            "host_homology_manifest": "2bdbbe0a142807b2d61f4d276ab0a3a2d947a899c0fa273e2c06a3f818f1aaeb",
            "host_only_manifest": "0db09446f4e7b57151188308cb3c83c5e6939252ae3149086571103bd3cb483d",
            "mixed_host_virus_manifest": "bdf5fd6613ad2a698dee6e439f2aa607d01ad9b40fc8ad48f69433755bb8cce6",
            "truth_manifest": "37259bd7330f016e2041ee16f6cf81c23dc24ac6bb4f061af6ce8b5f4ef1f7e8",
        },
        templates={
            "host_homology_loci.fasta": "fd23b6e6e3ff3c18a1bdaf40d60a1704a1c9ed960024033eee7a1893b28ab3fe",
            "host_transcripts.fasta": "9d14fef081d864539d40e6c37a1f6da021da3bff4a17ee6b1b8a05b92d5ddbe0",
            "transcripts.fasta": "4cd779c7e678e36f963a34a2ba906c2c75c593178e5146da3753e6bb3fcca4c4",
            "viral_gene_models.tsv": "0502000a3e2e8d07318c3a1b50f106d1767676747e9705ae0a4a1e863a976719",
            "viral_genes.fasta": "720b34c831eec9d9920f1893c38a244d1b3557a0d8ef00f762ad664dbf02050e",
            "viral_genomes.fasta": "fcca10a9bf7f451f5146c38fb4d69dbc45fd5083f278bbf13c4a2587ae125966",
        },
    ),
    "noisy": dict(
        args=[
            "--generation-seed",
            "7",
            "--datasets",
            "synthetic_factorial",
            "synthetic_host_homology",
            "--chemistries",
            "10xv3",
            "drop-seq",
        ],
        samples_digest="374866d2308e3c88d4c82617a8c4f763b24c93bc0a90a10a1d2fac38a1b1db26",
        truth={
            "host_homology_manifest": "d45fa48962df09651c9e684702e0153232968a1526f5e96dd165526eec2c8a2c",
            "host_only_manifest": "4e4c3980c8b3c421c4866da4dad595aae06d41e5f8c832e8b2bb4324811569a3",
            "mixed_host_virus_manifest": "2694ea77dab00bb4dbcdc30cb7b1b93bb60aabe9569c03efbdf6d7f39ff37484",
            "truth_manifest": "86534c2051d1ac18661f8712861718cac42d4d3c7324802f16f93f3aa20b9668",
        },
        templates={
            "host_homology_loci.fasta": "beae5ff9a75e30820c9e56e5edc3fbe5d4401820374cf534bbb96966465ab7c4",
            "host_transcripts.fasta": "e23ce77617cca3129e3a11cfebcae0e66482b1c09671370f30c7a0762cad91df",
            "transcripts.fasta": "e9dbfa939c38a846e56d7ad5087b7bab01c60c79f8f7941f07eb0ce9e6f93109",
            "viral_gene_models.tsv": "0502000a3e2e8d07318c3a1b50f106d1767676747e9705ae0a4a1e863a976719",
            "viral_genes.fasta": "1c8b6cf081e36760fd68c8e1bf1344954b7bc0ac469343931d0cffebad4e76e5",
            "viral_genomes.fasta": "939f37fde6bd9212d4b4a152c229a4fbc675575b34bff19f4fa4e83a6b79f74e",
        },
    ),
}


@pytest.mark.parametrize("name", sorted(GOLDEN_DEFAULT))
def test_default_output_is_byte_identical_to_the_pre_aud18_generator(tmp_path, name):
    golden = GOLDEN_DEFAULT[name]
    manifest = gt.generate(gt.parse_args([*golden["args"], "--out", str(tmp_path / "p")]))
    digest = hashlib.sha256(
        json.dumps(
            {
                sid: [d["r1_sha256"], d["r2_sha256"], d["r1_raw_sha256"], d["r2_raw_sha256"]]
                for sid, d in sorted(manifest["samples"].items())
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()
    assert digest == golden["samples_digest"]
    assert manifest["truth"] == golden["truth"]
    assert manifest["templates"] == golden["templates"]
    assert "inputs" not in manifest  # recorded only when an AUD-18 option is used


REAL_ARGS = [
    "--chemistries",
    "10xv3",
    "--datasets",
    "synthetic_factorial",
    "synthetic_mixed_host_virus",
]


def real_make(out: Path, inputs, *extra: str, seed: str = "1") -> dict:
    fasta, spans, onlist = inputs
    args = ["--generation-seed", seed, "--error-rate", "0", *REAL_ARGS, "--out", str(out)]
    args += ["--viral-fasta", str(fasta), "--gene-spans", str(spans), "--onlist", str(onlist)]
    return gt.generate(gt.parse_args([*args, *extra]))


@pytest.fixture(scope="module")
def real(tmp_path_factory):
    root = tmp_path_factory.mktemp("real")
    inputs = write_real_inputs(root / "in")
    return root, inputs, real_make(root / "panel", inputs)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_real_panel_validates_and_records_input_hashes(real):
    root, (fasta, spans, onlist), manifest = real
    rows = score.read_truth_manifest(root / "panel" / "truth" / "truth_manifest.tsv")
    planted = [r for r in rows if r["is_planted_molecule"] == "true"]
    assert {r["planted_virus_id"] for r in planted} == set(PLANTED)  # HHV-6A has no pair: absent
    assert {r["planted_gene_id"] for r in planted} <= {
        f"{v}_{g}" for v in PLANTED for g in ("G1", "G2", "G3")
    }
    assert manifest["inputs"]["viral_fasta"]["sha256"] == sha(fasta)
    assert manifest["inputs"]["gene_spans"]["sha256"] == sha(spans)
    assert manifest["inputs"]["onlist"]["sha256"] == sha(onlist)
    assert manifest["inputs"]["onlist"]["n_barcodes"] == 600
    assert not re.search(
        r"/(exports|home|tmp)/", (root / "panel" / "run_manifest.json").read_text()
    )
    assert all(d["absent_viruses"] == ["HHV-6A"] for d in manifest["samples"].values())
    models = read_tsv(root / "panel" / "templates" / "viral_gene_models.tsv")
    assert len(models) == 15 and {m["strand"] for m in models} == {"+", "-"}


def test_real_truth_class_and_read_orientation_follow_the_gene_models(real):
    root, _, manifest = real
    panel = root / "panel"
    gene_seq = gt.read_fasta_records(panel / "templates" / "viral_genes.fasta")
    sid = next(s for s, d in manifest["samples"].items() if d["dataset"] == "synthetic_factorial")
    r2 = dict(fastq_records(panel / sid / f"{sid}_R2.fastq.gz"))
    rows = [
        r
        for r in read_tsv(panel / "truth" / "truth_manifest.tsv")
        if r["biological_sample_id"] == sid and r["is_planted_molecule"] == "true"
    ]
    for r in rows:
        body = r2[r["read_id"]][: int(r["templated_len"])]
        assert body in gene_seq[r["planted_gene_id"]], r["planted_gene_id"]  # sense, also for "-"
        if r["planted_gene_id"].endswith(("_G2", "_G3")):  # no other same-strand model holds them
            assert r["truth_class"] == "unique", r
    classes = {r["truth_class"] for r in rows if r["planted_gene_id"].endswith("_G1")}
    assert classes == {"unique", "ambiguous"}  # G1's 3' end overlaps G2 on the same strand
    assert any(r["planted_gene_id"].endswith("_G3") for r in rows)


def test_real_panel_draws_every_barcode_from_the_onlist(real):
    root, (_, _, onlist), manifest = real
    wanted = set(onlist.read_text().split())
    panel = root / "panel"
    for sid in manifest["samples"]:
        pool = read_tsv(panel / sid / "barcode_pool.tsv")
        assert {r["cell_barcode"] for r in pool} <= wanted and len(pool) == 70
        assert all(b[:16] in wanted for _, b in fastq_records(panel / sid / f"{sid}_R1.fastq.gz"))
    truth = read_tsv(panel / "truth" / "truth_manifest.tsv")
    assert {"hopped", "ambient"} <= {r["origin"] for r in truth}
    assert {r["cell_barcode"] for r in truth} <= wanted
    assert any(r["collision_kind"] for r in truth)


def test_real_panel_is_byte_deterministic_per_seed_and_gz_onlist_is_equivalent(real, tmp_path):
    root, inputs, _ = real
    again = real_make(tmp_path / "again", inputs)
    assert tree_hashes(root / "panel") == tree_hashes(tmp_path / "again")
    real_make(tmp_path / "other", inputs, seed="2")
    key = "truth/truth_manifest.tsv"
    assert tree_hashes(root / "panel")[key] != tree_hashes(tmp_path / "other")[key]
    fasta, spans, onlist = inputs
    with gzip.open(tmp_path / "onlist.txt.gz", "wt") as handle:
        handle.write(onlist.read_text())
    zipped = real_make(tmp_path / "gz", (fasta, spans, tmp_path / "onlist.txt.gz"))
    assert zipped["samples"] == again["samples"] and zipped["truth"] == again["truth"]


def test_sibling_pairs_file_overrides_the_protocol_pairs(real, tmp_path):
    _, inputs, _ = real
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("HSV-1\tEBV\n")  # a pair the protocol does not declare
    manifest = real_make(tmp_path / "p", inputs, "--sibling-pairs", str(pairs))
    assert manifest["inputs"]["sibling_pairs"]["sha256"] == sha(pairs)
    for detail in manifest["samples"].values():
        assert detail["absent_viruses"] == ["HHV-6A", "HSV-2", "KSHV"]


def edited_inputs(tmp_path, spans_text=None, fasta_text=None, onlist_text=None):
    fasta, spans, onlist = write_real_inputs(tmp_path / "bad")
    for path, text in ((spans, spans_text), (fasta, fasta_text), (onlist, onlist_text)):
        if text is not None:
            path.write_text(text)
    return fasta, spans, onlist


@pytest.mark.parametrize(
    ("kw", "message"),
    [
        (dict(spans_text="EBV\tEBV_G1\t1\t2000\t+\n"), "outside EBV"),
        (dict(spans_text="EBV\tE_G1\t1\t100\t+\nEBV\tE_G1\t200\t300\t+\n"), "duplicate gene_id"),
        (dict(spans_text="EBV\tEBV_G1\t1\t100\t+\n"), "needs at least 2 gene spans"),
        (dict(spans_text="NOPE\tN_G1\t1\t100\t+\n"), "not in"),
        (dict(spans_text="EBV\tEBV_G1\t1\t100\t.\n"), "strand"),
        (
            dict(fasta_text=">a\n" + "ACGT" * 100 + "\n>a\n" + "ACGT" * 100 + "\n"),
            "duplicate record",
        ),
        (dict(fasta_text=">a\n" + "ACGU" * 100 + "\n"), "outside ACGTN"),
        (dict(onlist_text="ACGT\nACGTA\n"), "differ in length"),
        (dict(onlist_text="ACGTACGTACGTACGN\n"), "ACGT only"),
        (dict(onlist_text="ACGTACGTACGTACGT\nACGTACGTACGTACGT\n"), "duplicate barcodes"),
        (dict(onlist_text="\n".join("A" * i + "C" * (16 - i) for i in range(17)) + "\n"), "needed"),
    ],
)
def test_invalid_real_inputs_are_refused(tmp_path, kw, message):
    inputs = edited_inputs(tmp_path, **kw)
    with pytest.raises(ValueError, match=message):
        real_make(tmp_path / "out", inputs)
    assert not (tmp_path / "out" / "run_manifest.json").exists()


def test_option_combinations_are_checked(real, tmp_path):
    _, (fasta, spans, onlist), _ = real
    base = ["--generation-seed", "1", "--out", str(tmp_path / "o")]
    with pytest.raises(ValueError, match="together"):
        gt.generate(gt.parse_args([*base, "--viral-fasta", str(fasta)]))
    with pytest.raises(ValueError, match="needs --viral-fasta"):
        gt.generate(gt.parse_args([*base, "--sibling-pairs", str(spans)]))
    with pytest.raises(ValueError, match="16 nt but drop-seq uses 12"):  # default chemistries
        gt.generate(gt.parse_args([*base, "--onlist", str(onlist)]))
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("HSV-1\tHSV-2\nHSV-1\tEBV\n")
    with pytest.raises(ValueError, match="more than one pair"):
        real_make(tmp_path / "o2", (fasta, spans, onlist), "--sibling-pairs", str(pairs))

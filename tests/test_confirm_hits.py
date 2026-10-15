import gzip
import random
import shutil
import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis" / "translated_search"))
import confirm_hits as ch  # noqa: E402

rng = random.Random(7)


def rand_dna(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


PROTEINS = {
    "V1.1__P1": {"protein_id": "P1", "virus_id": "taxon_1", "gene": "G1", "product": "tegument protein UL82", "accession": "V1.1", "species": "sp"},
    "V1.1__P2": {"protein_id": "P2", "virus_id": "taxon_1", "gene": "G2/G3", "product": "capsid", "accession": "V1.1", "species": "sp"},
    "V2.1__P3": {"protein_id": "P3", "virus_id": "taxon_2", "gene": "", "product": "replication protein", "accession": "V2.1", "species": "sp"},
}  # fmt: skip


def test_trim_and_complexity():
    core = rand_dna(60)
    trimmed, was = ch.trim_tails(core + "A" * 25)
    assert was and trimmed == core
    assert ch.trim_tails("G" * 12 + core)[0].startswith("G" * 12)  # leading G is not a tail
    assert ch.complexity(core)["low_complexity"] is False
    assert ch.complexity("ACACACACACACACACACACACACACACACACACACACACACACAC")["low_complexity"]
    assert ch.complexity(core[:20] + "G" * 12)["low_complexity"]  # short after trimming
    assert ch.dust_score("A" * 60) > ch.DUST_MAX > ch.dust_score(core)


def test_tso_is_stripped_and_flagged():
    core = rand_dna(60)
    read = "AAGCAGTGGTATCAACGCAGAGTAC" + "T" * 28 + rand_dna(12)
    c = ch.complexity(read)
    assert c["tso_adapter"] and c["low_complexity"]  # little left after the adapter and poly(T)
    c2 = ch.complexity(rc(read))  # antisense orientation
    assert c2["tso_adapter"] and c2["low_complexity"]
    c3 = ch.complexity(
        "AAGCAGTGGTATCAACGCAGAGTAC" + core
    )  # real insert behind the adapter survives
    assert c3["tso_adapter"] and c3["trimmed_seq"] == core and not c3["low_complexity"]
    assert not ch.complexity(core)["tso_adapter"]


def test_group_specs_and_unknown_spec():
    g = ch.parse_group_specs(["P3", "taxon_1", "taxon_1:G3", "taxon_1:UL82"], PROTEINS)
    assert g["P3"] == {"V2.1__P3"}
    assert g["taxon_1"] == {"V1.1__P1", "V1.1__P2"}
    assert g["taxon_1:G3"] == {"V1.1__P2"}  # slash-separated gene aliases
    assert g["taxon_1:UL82"] == {"V1.1__P1"}  # token of the product name
    with pytest.raises(SystemExit):
        ch.parse_group_specs(["nope"], PROTEINS)


def test_collect_hits_drops_multi_virus_ec():
    ec = {0: ["V1.1__P1"], 1: ["V1.1__P1", "V1.1__P2"], 2: ["V1.1__P1", "V2.1__P3"]}
    rows = [("AAA", "U1", 0, 0), ("AAA", "U2", 1, 1), ("CCC", "U3", 2, 2)]
    h = ch.collect_hits(rows, ec, PROTEINS)
    assert h["ambiguous_records"] == 1
    assert [x[2] for x in h["protein"]["V1.1__P1"]] == [0, 1]
    assert [x[2] for x in h["protein"]["V1.1__P2"]] == [1]
    assert len(h["virus"]["taxon_1"]) == 2 and "taxon_2" not in h["virus"]
    assert ch.top_group_specs(h, PROTEINS, 1, 1) == ["taxon_1", "V1.1__P1"]


def test_umi_stats_patterns():
    diffuse = [(f"B{i}", "U", i) for i in range(200)]
    s = ch.umi_stats(diffuse)
    assert s["pattern"] == "diffuse" and s["frac_barcodes_1_umi"] == 1.0
    conc = [(f"B{i % 3}", f"U{i}", i) for i in range(60)]
    assert ch.umi_stats(conc)["pattern"] == "concentrated"
    assert ch.umi_stats([("B", "U", 0)])["pattern"] == "sparse"
    assert ch.umi_stats([])["pattern"] == "none"


def test_sample_reads_is_seeded_and_reports_how():
    flags = list(range(5000))
    a, how = ch.sample_reads(flags, 100, seed=3)
    assert a == ch.sample_reads(flags, 100, seed=3)[0] and len(a) == 100 and "seed 3" in how
    assert ch.sample_reads([5, 1, 1], 100, 1) == ([1, 5], "all reads")


def test_fetch_reads_zero_based(tmp_path):
    fq = tmp_path / "r.fastq.gz"
    with gzip.open(fq, "wt") as fh:
        for i in range(10):
            fh.write(f"@read{i} 2:N\n{'ACGT' * (i + 1)}\n+\n{'I' * 4 * (i + 1)}\n")
    got = ch.fetch_reads(fq, {0, 3})
    assert got[0] == ("read0", "ACGT") and got[3][0] == "read3" and len(got[3][1]) == 16


def test_parse_paf_keeps_best_and_identity():
    paf = [
        "r1\t100\t0\t50\t+\tT1\t999\t0\t50\t45\t50\t60",
        "r1\t100\t0\t90\t+\tT2\t999\t0\t90\t88\t90\t60",
        "short",
    ]
    best = ch.parse_paf(paf)
    assert best["r1"]["target"] == "T2" and best["r1"]["identity"] == pytest.approx(88 / 90)
    assert best["r1"]["qcov"] == pytest.approx(0.9)
    per = ch.parse_paf(paf, per_target=True)
    assert set(per["r1"]) == {"T1", "T2"} and per["r1"]["T1"]["matches"] == 45


def hit(identity=0.99, qcov=0.9):
    return {"target": "X", "identity": identity, "qcov": qcov, "bitscore": 50.0}


def base(**kw):
    r = {"low_complexity": False, "host_genome": None, "host_tx": None, "viral_same": None,
         "viral_other": None, "prot_viral": None, "prot_human": None}  # fmt: skip
    return r | kw


@pytest.mark.parametrize(
    "kw,expected",
    [
        ({"viral_same": hit()}, "viral_supported"),
        ({"viral_same": hit(), "host_genome": hit()}, "unresolved"),
        ({"host_tx": hit()}, "host_derived"),
        ({"host_genome": hit(0.8)}, "unresolved"),  # below the identity gate
        ({"viral_other": hit()}, "misassigned_virus"),
        ({"prot_viral": hit(), "prot_human": hit()}, "divergent_homolog"),
        ({"prot_viral": hit()}, "protein_only_viral"),
        ({"viral_same": hit(), "low_complexity": True}, "low_complexity"),
        ({"host_genome": hit(), "low_complexity": True}, "host_derived"),
        ({"tso_adapter": True, "low_complexity": True}, "adapter_artefact"),
        ({"tso_adapter": True}, "adapter_artefact"),  # adapter, remainder matches nothing
        ({"tso_adapter": True, "prot_viral": hit()}, "protein_only_viral"),
        ({}, "unresolved"),
    ],
)
def test_read_verdict(kw, expected):
    assert ch.read_verdict(base(**kw), 0.9, 0.5)[0] == expected


def test_group_verdict_needs_majority():
    assert ch.group_verdict(Counter(viral_supported=6, host_derived=4)) == "viral_supported"
    assert ch.group_verdict(Counter(viral_supported=4, host_derived=3, unresolved=3)) == "mixed"
    assert ch.group_verdict(Counter()) == "no_reads"


def test_genbank_sequence(tmp_path):
    gb = tmp_path / "x.gb"
    gb.write_text("LOCUS x\nORIGIN\n        1 acgtacgtac gg\n       61 tt\n//\nTAIL 1 aaaa\n")
    assert ch.genbank_sequence(gb) == "ACGTACGTACGGTT"


@pytest.mark.skipif(
    not all(shutil.which(t) for t in ("minimap2", "blastx", "makeblastdb")),
    reason="minimap2/BLAST missing",
)
def test_end_to_end_synthetic(tmp_path):
    """Viral reads align to the viral genome, host reads to the host genome, a polyA-tailed viral read survives."""
    viral = rand_dna(3000)
    host = rand_dna(6000)
    codon = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG"}  # fmt: skip
    prot = "".join(rng.choice("ACDEFGHIKL") for _ in range(120))
    (tmp_path / "v.fa").write_text(f">V1.1\n{viral}\n")
    (tmp_path / "h.fa").write_text(f">chr1\n{host}\n")
    (tmp_path / "tx.fa").write_text(f">tx1\n{host[1000:2500]}\n")
    (tmp_path / "viral.faa").write_text(f">V1.1__P1\n{prot}\n")
    with gzip.open(tmp_path / "human.faa.gz", "wt") as fh:  # gzipped, like the Ensembl peptide file
        fh.write(">H1\n" + "".join(rng.choice("ACDEFGHIKL") for _ in range(120)) + "\n")
    reads = [viral[100 + 5 * i : 190 + 5 * i] for i in range(8)]  # 8 viral reads
    reads.append(rc(viral[900:960]) + "A" * 30)  # viral, antisense, polyA tail
    reads += [host[1100 + 7 * i : 1190 + 7 * i] for i in range(6)]  # 6 host reads
    reads.append(
        "".join(codon[a] for a in prot[:30])
    )  # protein-only read (synthetic codons: no nt match)
    fq = tmp_path / "r.fastq"
    fq.write_text("".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(reads)))
    ec = {0: ["V1.1__P1"]}
    rows = [(f"B{i % 4}", f"U{i}", 0, i) for i in range(len(reads))]
    proteins = {
        "V1.1__P1": {
            "protein_id": "P1",
            "virus_id": "taxon_1",
            "gene": "G",
            "product": "p",
            "accession": "V1.1",
            "species": "sp",
        }
    }
    hits = ch.collect_hits(rows, ec, proteins)
    refs = {
        "viral_genomes": tmp_path / "v.fa",
        "host_genome": tmp_path / "h.fa",
        "host_tx": tmp_path / "tx.fa",
        "viral_db": ch.make_blast_db(tmp_path / "viral.faa", tmp_path / "db" / "v"),
        "human_db": ch.make_blast_db(tmp_path / "human.faa.gz", tmp_path / "db" / "h"),
    }
    res = ch.run_confirm(
        hits=hits,
        proteins=proteins,
        groups={"taxon_1": {"V1.1__P1"}},
        reads_fastq=fq,
        work=tmp_path / "w",
        refs=refs,
        max_reads=100,
        threads=1,
    )
    g = res["groups"]["taxon_1"]
    assert g["reads_checked"] == len(reads) and g["reads_total"] == len(reads)
    assert g["n_viral_supported"] == 9  # includes the polyA-tailed antisense read
    assert g["n_host_derived"] == 6 and g["pct_host_any"] == pytest.approx(37.5)
    assert g["pct_tail_trimmed"] == pytest.approx(100 / 16, abs=0.1)
    ch.write_outputs(res, tmp_path / "out", "S", {})
    assert (tmp_path / "out" / "S.confirm.md").read_text().count("taxon_1") >= 1


def test_recurrence_counts_samples_and_share():
    ec = {0: ["V1.1__P1"], 1: ["V2.1__P3"]}
    a = ch.collect_hits([("A", "U1", 0, 0), ("B", "U2", 0, 1), ("B", "U3", 1, 2)], ec, PROTEINS)
    b = ch.collect_hits([("C", "U1", 0, 0)], ec, PROTEINS)
    rows = {r["target"]: r for r in ch.recurrence({"a": (a, None), "b": (b, None)}, PROTEINS)}
    assert rows["V1.1__P1"]["samples_present"] == 2 and rows["V1.1__P1"]["a_umis"] == 2
    assert rows["V1.1__P1"]["b_share_pct"] == 100.0 and rows["V2.1__P3"]["samples_present"] == 1

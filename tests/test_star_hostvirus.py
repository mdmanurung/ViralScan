"""STARsolo host-unmapped alignment against a combined HOST|/VIRUS| index (STAR-HOSTVIRUS-01)."""

import gzip
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from viralscan import star_hostvirus, star_run
from viralscan.star_hostvirus import (
    STAR_SUMMARY_COLUMNS,
    align_cmd,
    align_host_unmapped,
    build_star_hostvirus_index,
    prune_to_viral,
    resolve_whitelist,
    summarise_by_virus,
    write_reference,
    write_summary,
)

needs_samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools")


def _seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_write_reference_prefixes(tmp_path):
    (tmp_path / "h.fa").write_text(">chr1 x\nACGT\n")
    (tmp_path / "v.fa").write_text(">AB1.1\nACGT\n")
    (tmp_path / "h.gtf").write_text('chr1\ts\texon\t1\t4\t.\t+\t.\tgene_id "g";\n')
    (tmp_path / "v.gtf").write_text('AB1.1\ts\texon\t1\t4\t.\t+\t.\tgene_id "v";\n')
    fa, gtf = write_reference(
        tmp_path / "h.fa", tmp_path / "v.fa", tmp_path / "h.gtf", tmp_path / "v.gtf", tmp_path / "o"
    )
    assert [x for x in fa.read_text().splitlines() if x[0] == ">"] == [
        ">HOST|chr1 x",
        ">VIRUS|AB1.1",
    ]
    assert [x.split("\t")[0] for x in gtf.read_text().splitlines()] == ["HOST|chr1", "VIRUS|AB1.1"]


def test_align_cmd_geometry_and_multimappers():
    cmd = align_cmd("STAR", "idx", "r1.fq.gz", "r2.fq.gz", "10xv2", "out/", 4)
    assert cmd[cmd.index("--readFilesIn") + 1 : cmd.index("--readFilesIn") + 3] == [
        "r2.fq.gz",
        "r1.fq.gz",
    ]
    assert cmd[cmd.index("--soloCBlen") + 1] == "16" and cmd[cmd.index("--soloUMIlen") + 1] == "10"
    assert cmd[cmd.index("--outFilterMultimapNmax") + 1] == "100"
    assert cmd[cmd.index("--soloCBwhitelist") + 1] == "None"
    assert cmd[cmd.index("--soloStrand") + 1] == "Unstranded"


@pytest.mark.parametrize(
    "strand,solo", [("forward", "Forward"), ("reverse", "Reverse"), ("unstranded", "Unstranded")]
)
def test_align_cmd_strand(strand, solo):
    cmd = align_cmd("STAR", "idx", "r1", "r2", "10xv2", "o/", 1, strand=strand)
    assert cmd.count("--soloStrand") == 1 and cmd[cmd.index("--soloStrand") + 1] == solo


def test_align_cmd_rejects_unknown_strand():
    with pytest.raises(ValueError, match="strand"):
        align_cmd("STAR", "idx", "r1", "r2", "10xv2", "o/", 1, strand="Forward")


def test_resolve_whitelist_default_explicit_and_none(tmp_path, monkeypatch):
    gz = tmp_path / "list.txt.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("ACGT\nTTTT\n")
    monkeypatch.setattr(star_hostvirus, "onlist_path", lambda tech: str(gz))
    plain = resolve_whitelist("10xv2", None, tmp_path)
    assert plain == str(tmp_path / "whitelist.txt")  # gunzipped: STAR cannot read .gz
    assert (tmp_path / "whitelist.txt").read_text() == "ACGT\nTTTT\n"
    assert resolve_whitelist("10xv2", "/x/mine.txt", tmp_path) == "/x/mine.txt"  # explicit wins


def test_resolve_whitelist_none_without_list(tmp_path):
    # Drop-seq has no on-list; a geometry string is not a named chemistry (SW-21).
    assert resolve_whitelist("dropseq", None, tmp_path) is None
    assert resolve_whitelist("0,0,12:0,12,20:1,0,0", None, tmp_path) is None


FINAL_ZERO = "Uniquely mapped reads number |\t0\nNumber of reads mapped to multiple loci |\t0\n"
FINAL_MAPPED = "Uniquely mapped reads number |\t3\nNumber of reads mapped to multiple loci |\t0\n"


def _fake_run(monkeypatch, work, codes, final):
    calls = []

    def run(cmd, check=True):
        calls.append(list(cmd))
        code = codes[len(calls) - 1]
        if code:
            raise subprocess.CalledProcessError(code, cmd)
        (work / "Log.final.out").write_text(final)

    monkeypatch.setattr(star_run.subprocess, "run", run)
    return calls


def test_star139_retry_uses_barcode_length_dummy(tmp_path, monkeypatch):
    calls = _fake_run(monkeypatch, tmp_path, [139, 0], FINAL_ZERO)
    star_run.run_star(["STAR", "--soloCBwhitelist", "None"], tmp_path, cb_len=16)
    assert len(calls) == 2
    dummy = tmp_path / "dummy_whitelist.txt"
    assert calls[1][calls[1].index("--soloCBwhitelist") + 1] == str(dummy)
    assert dummy.read_text() == "A" * 16 + "\n"
    assert (tmp_path / "star139_fallback.txt").is_file()


def test_star139_retry_that_maps_reads_is_refused(tmp_path, monkeypatch):
    _fake_run(monkeypatch, tmp_path, [139, 0], FINAL_MAPPED)
    with pytest.raises(RuntimeError, match="refusing to drop barcodes"):
        star_run.run_star(["STAR", "--soloCBwhitelist", "None"], tmp_path)


@pytest.mark.parametrize("code,whitelist", [(139, "wl.txt"), (1, "None")])
def test_star139_other_failures_reraised(tmp_path, monkeypatch, code, whitelist):
    calls = _fake_run(monkeypatch, tmp_path, [code], FINAL_ZERO)
    with pytest.raises(subprocess.CalledProcessError):
        star_run.run_star(["STAR", "--soloCBwhitelist", whitelist], tmp_path)
    assert len(calls) == 1


def _capture_align(tmp_path, monkeypatch, technology, whitelist=None):
    """Run align_host_unmapped with STAR/samtools stubbed; return (cmd, run_star kwargs)."""
    seen = {}

    def fake_run_star(cmd, work, cb_len=12):
        seen["cmd"], seen["cb_len"] = cmd, cb_len

    monkeypatch.setattr(star_hostvirus, "run_star", fake_run_star)
    monkeypatch.setattr(star_hostvirus.subprocess, "run", lambda *a, **k: None)
    align_host_unmapped(
        "r1", "r2", "idx", tmp_path / "o", technology, whitelist=whitelist, threads=1
    )
    return seen["cmd"], seen["cb_len"]


def test_align_host_unmapped_default_onlist_reaches_star(tmp_path, monkeypatch):
    gz = tmp_path / "list.txt.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("ACGT\n")
    monkeypatch.setattr(star_hostvirus, "onlist_path", lambda tech: str(gz))
    cmd, _ = _capture_align(tmp_path, monkeypatch, "10xv2")
    wl = cmd[cmd.index("--soloCBwhitelist") + 1]
    assert wl == str(tmp_path / "o" / "whitelist.txt")  # not STAR's raw "None"


@pytest.mark.parametrize(
    "technology,cb_len", [("10xv2", 16), ("dropseq", 12), ("0,0,14:0,14,24:1,0,0", 14)]
)
def test_align_host_unmapped_passes_barcode_length_to_run_star(
    tmp_path, monkeypatch, technology, cb_len
):
    monkeypatch.setattr(star_hostvirus, "onlist_path", lambda tech: "None")
    cmd, got = _capture_align(tmp_path, monkeypatch, technology, whitelist="wl.txt")
    assert got == cb_len
    assert cmd[cmd.index("--soloCBlen") + 1] == str(cb_len)


def _rec(q, flag, ref, pos, nh, hi, cb, ub):
    return (
        f"{q}\t{flag}\t{ref}\t{pos}\t255\t10M\t*\t0\t0\t{'A' * 10}\t{'I' * 10}"
        f"\tNH:i:{nh}\tHI:i:{hi}\tCB:Z:{cb}\tUB:Z:{ub}\n"
    )


def _synthetic_bam(tmp_path):
    """Hand-computable STAR-style BAM; see test_summarise_by_virus for the expected counts."""
    head = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(
        f"@SQ\tSN:{c}\tLN:1000\n" for c in ("HOST|chr1", "VIRUS|V1", "VIRUS|V2")
    )
    recs = [
        _rec("u1", 0, "VIRUS|V1", 10, 1, 1, "AAAA", "UUU1"),
        _rec("u2", 16, "VIRUS|V1", 20, 1, 1, "AAAA", "UUU1"),  # same molecule as u1
        _rec("u3", 0, "VIRUS|V1", 30, 1, 1, "BBBB", "UUU2"),
        _rec("mh", 0, "VIRUS|V1", 40, 2, 1, "CCCC", "UUU3"),  # V1 primary, host partner
        _rec("mh", 256, "HOST|chr1", 50, 2, 2, "CCCC", "UUU3"),
        _rec("mv", 0, "VIRUS|V1", 60, 2, 1, "DDDD", "UUU4"),  # V1 primary, V2 secondary
        _rec("mv", 256, "VIRUS|V2", 60, 2, 2, "DDDD", "UUU4"),
        _rec("nocb", 0, "VIRUS|V2", 70, 1, 1, "-", "-"),
        _rec("sh", 0, "HOST|chr1", 80, 2, 1, "EEEE", "UUU5"),  # host primary, V2 secondary
        _rec("sh", 256, "VIRUS|V2", 90, 2, 2, "EEEE", "UUU5"),
        _rec("hostonly", 0, "HOST|chr1", 100, 1, 1, "FFFF", "UUU6"),
    ]
    sam = tmp_path / "x.sam"
    sam.write_text(head + "".join(recs))
    bam = tmp_path / "Aligned.sortedByCoord.out.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    return bam


@needs_samtools
def test_summarise_by_virus(tmp_path):
    rows = summarise_by_virus(_synthetic_bam(tmp_path))
    # V1: u1 u2 u3 mh mv (all primary on V1); molecules (A,1) (B,2) (C,3) (D,4);
    #     host_competing = mh; other_virus = mv.
    # V2: mv (secondary), nocb (primary, CB "-"), sh (secondary); molecules (D,4) (E,5);
    #     host_competing = sh; other_virus = mv.
    assert rows == [
        {"accession": "V1", "star_reads_primary": 5, "star_reads_any": 5, "star_molecules": 4,
         "star_host_competing": 1, "star_other_virus": 1},
        {"accession": "V2", "star_reads_primary": 1, "star_reads_any": 3, "star_molecules": 2,
         "star_host_competing": 1, "star_other_virus": 1},
    ]  # fmt: skip
    tsv = write_summary(rows, tmp_path / "s.tsv").read_text().splitlines()
    assert tsv[0].split("\t") == list(STAR_SUMMARY_COLUMNS)
    assert tsv[1].split("\t") == ["V1", "5", "5", "4", "1", "1"]


@needs_samtools
def test_summarise_by_virus_without_viral_reads(tmp_path):
    sam = tmp_path / "h.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:HOST|chr1\tLN:1000\n@SQ\tSN:VIRUS|V1\tLN:1000\n"
        + _rec("h", 0, "HOST|chr1", 5, 1, 1, "AAAA", "UUU1")
    )
    bam = tmp_path / "h.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    assert summarise_by_virus(bam) == []
    viral = prune_to_viral(bam)
    out = subprocess.run(["samtools", "view", str(viral)], capture_output=True, text=True)
    assert out.stdout == "" and not bam.exists()


@needs_samtools
def test_prune_keeps_only_viral_records(tmp_path):
    bam = _synthetic_bam(tmp_path)
    (tmp_path / "_STARtmp").mkdir()
    (tmp_path / "_STARtmp" / "junk").write_text("x")
    viral = prune_to_viral(bam)
    assert viral.is_file() and viral.with_name(viral.name + ".bai").is_file()
    assert not bam.exists() and not bam.with_name(bam.name + ".bai").exists()
    assert not (tmp_path / "_STARtmp").exists()
    view = subprocess.run(["samtools", "view", str(viral)], capture_output=True, text=True)
    refs = [ln.split("\t")[2] for ln in view.stdout.splitlines()]
    assert len(refs) == 8 and all(r.startswith("VIRUS|") for r in refs)


@pytest.mark.skipif(
    shutil.which("STAR") is None or shutil.which("samtools") is None, reason="needs STAR, samtools"
)
def test_end_to_end_viral_and_host_reads(tmp_path):
    host, virus = _seq(20000, 1), _seq(3000, 2)
    (tmp_path / "h.fa").write_text(f">chr1\n{host}\n")
    (tmp_path / "v.fa").write_text(f">V1.1\n{virus}\n")
    (tmp_path / "h.gtf").write_text(
        'chr1\ts\texon\t100\t900\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    )
    (tmp_path / "v.gtf").write_text(
        'V1.1\ts\texon\t100\t900\t.\t+\t.\tgene_id "v"; transcript_id "tv";\n'
    )
    fa, gtf = write_reference(
        tmp_path / "h.fa",
        tmp_path / "v.fa",
        tmp_path / "h.gtf",
        tmp_path / "v.gtf",
        tmp_path / "ref",
    )
    index = build_star_hostvirus_index(fa, gtf, tmp_path / "idx", threads=1, sa_index_nbases=6)
    cb, umi = "ACGTACGTACGTACGT", "TTGCAATTGC"
    r1, r2 = [], []
    for i, src in enumerate([virus] * 3 + [host] * 2):
        pos = 200 + 10 * i
        r2.append(f"@r{i}\n{src[pos : pos + 90]}\n+\n{'I' * 90}\n")
        r1.append(f"@r{i}\n{cb}{umi}\n+\n{'I' * 26}\n")
    (tmp_path / "R1.fq").write_text("".join(r1))
    (tmp_path / "R2.fq").write_text("".join(r2))
    (tmp_path / "wl.txt").write_text(f"{cb}\n")
    bam = align_host_unmapped(
        tmp_path / "R1.fq",
        tmp_path / "R2.fq",
        index,
        tmp_path / "aln",
        "10xv2",
        whitelist=tmp_path / "wl.txt",
        threads=1,
        strand="forward",
    )
    assert bam.with_name(bam.name + ".bai").exists()
    sam = subprocess.run(["samtools", "view", str(bam)], capture_output=True, text=True, check=True)
    rows = [line.split("\t") for line in sam.stdout.splitlines()]
    assert sorted(r[2] for r in rows) == ["HOST|chr1"] * 2 + ["VIRUS|V1.1"] * 3
    assert all(f"CB:Z:{cb}" in r and f"UB:Z:{umi}" in r for r in rows)
    # The same BAM through the summary and the prune path.
    assert summarise_by_virus(bam) == [{
        "accession": "V1.1", "star_reads_primary": 3, "star_reads_any": 3, "star_molecules": 1,
        "star_host_competing": 0, "star_other_virus": 0,
    }]  # fmt: skip
    pruned = prune_to_viral(bam)
    view = subprocess.run(["samtools", "view", str(pruned)], capture_output=True, text=True)
    assert [ln.split("\t")[2] for ln in view.stdout.splitlines()] == ["VIRUS|V1.1"] * 3


# ── star-index verb + manifest (STAR-03) ────────────────────────────────────────────
def _inputs(tmp_path):
    files = {
        "host_fasta": (">chr1\nACGTACGTAC\n>chr2\nTTGGCCAATT\n", "h.fa"),
        "host_gtf": ('chr1\tx\tgene\t1\t9\t.\t+\t.\tgene_id "g";\n', "h.gtf"),
        "viral_fasta": (">NC_1.1\nACGTACGTAC\n", "v.fa"),
        "viral_gtf": ('NC_1.1\tx\tgene\t1\t9\t.\t+\t.\tgene_id "v";\n', "v.gtf"),
    }
    out = {}
    for key, (text, name) in files.items():
        (tmp_path / name).write_text(text)
        out[key] = str(tmp_path / name)
    return out


def _fake_build(monkeypatch):
    def build(fasta, gtf, index_dir, threads=8, ram_bytes=0, star="STAR", sa_index_nbases=14):
        index_dir.mkdir(parents=True, exist_ok=True)
        for f in star_hostvirus._INDEX_FILES:
            (index_dir / f).write_text("x")
        return index_dir

    monkeypatch.setattr(star_hostvirus, "build_star_hostvirus_index", build)
    monkeypatch.setattr(star_hostvirus, "star_version", lambda star="STAR": "2.7.11b")


def _run_verb(monkeypatch, tmp_path, *extra, which="/usr/bin/STAR"):
    import shutil as _shutil

    from viralscan import menu

    ins = _inputs(tmp_path)
    argv = ["viralscan", "star-index", "--out-dir", str(tmp_path / "out"), *extra]
    for key, path in ins.items():
        argv += ["--" + key.replace("_", "-"), path]
    monkeypatch.setattr("sys.argv", argv)
    monkeypatch.setattr(_shutil, "which", lambda name: which)
    menu.main()
    return ins


def test_star_index_parser_options():
    from viralscan.cli.parsers import build_parser

    ns = build_parser().parse_args(
        [
            "star-index",
            "--host-fasta",
            "a",
            "--host-gtf",
            "b",
            "--viral-fasta",
            "c",
            "--viral-gtf",
            "d",
            "--out-dir",
            "o",
        ]  # fmt: skip
    )
    assert ns._subcommand == "star-index"
    assert (ns.threads, ns.force) == (8, False)
    assert ns.ram_bytes == star_hostvirus.DEFAULT_RAM_BYTES
    with pytest.raises(SystemExit):
        build_parser().parse_args(["star-index", "--out-dir", "o"])


def test_star_index_verb_writes_manifest(monkeypatch, tmp_path):
    _fake_build(monkeypatch)
    ins = _run_verb(monkeypatch, tmp_path, "--threads", "3")
    index = tmp_path / "out" / "index"
    m = star_hostvirus.read_index_manifest(index)
    assert star_hostvirus.validate_index_manifest(index, verify_inputs=True) == []
    assert m["contigs"] == {"host": 2, "virus": 1}
    assert m["star_version"] == "2.7.11b" and m["threads"] == 3
    import hashlib

    digest = hashlib.sha256(Path(ins["host_fasta"]).read_bytes()).hexdigest()
    assert m["inputs"]["host_fasta"]["sha256"] == digest
    cmd = m["genome_generate_args"]
    assert cmd[cmd.index("--runMode") + 1] == "genomeGenerate"
    assert cmd[cmd.index("--runThreadN") + 1] == "3"
    assert cmd[cmd.index("--genomeDir") + 1] == str(index)


def test_star_index_verb_requires_star(monkeypatch, tmp_path, capsys):
    _fake_build(monkeypatch)
    with pytest.raises(SystemExit) as exc:
        _run_verb(monkeypatch, tmp_path, which=None)
    assert exc.value.code == 1
    assert "STAR" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


def test_star_index_refuses_overwrite_then_force(monkeypatch, tmp_path):
    _fake_build(monkeypatch)
    _run_verb(monkeypatch, tmp_path)
    with pytest.raises(SystemExit):
        _run_verb(monkeypatch, tmp_path)
    _run_verb(monkeypatch, tmp_path, "--force")
    assert star_hostvirus.validate_index_manifest(tmp_path / "out" / "index") == []


def test_star_index_force_spares_foreign_dir(monkeypatch, tmp_path):
    _fake_build(monkeypatch)
    foreign = tmp_path / "out" / "index"
    foreign.mkdir(parents=True)
    (foreign / "precious.txt").write_text("keep")
    with pytest.raises(SystemExit):
        _run_verb(monkeypatch, tmp_path, "--force")
    assert (foreign / "precious.txt").exists()


def test_validate_index_manifest_reports_problems(monkeypatch, tmp_path):
    assert "missing" in star_hostvirus.validate_index_manifest(tmp_path)[0]
    _fake_build(monkeypatch)
    ins = _run_verb(monkeypatch, tmp_path)
    index = tmp_path / "out" / "index"
    (index / "SA").unlink()
    with Path(ins["viral_gtf"]).open("a") as handle:
        handle.write("# changed\n")
    problems = star_hostvirus.validate_index_manifest(index, verify_inputs=True)
    assert any("SA missing" in p for p in problems)
    assert any("viral_gtf changed" in p for p in problems)


@pytest.mark.parametrize(
    "tamper,expect",
    [
        (lambda m, p, ins: p.write_text("{"), "not valid JSON"),
        (lambda m, p, ins: m.update(schema="other/9"), "schema 'other/9'"),
        (lambda m, p, ins: m.pop("threads"), "manifest field 'threads' missing"),
        (lambda m, p, ins: m.update(contigs={"host": 0, "virus": 3}), "no HOST| contigs"),
        (lambda m, p, ins: m.update(contigs={"host": 3, "virus": 0}), "no VIRUS| contigs"),
        (lambda m, p, ins: Path(ins["viral_gtf"]).unlink(), "viral_gtf not readable"),
    ],
)
def test_validate_index_manifest_problem_branches(monkeypatch, tmp_path, tamper, expect):
    import json

    _fake_build(monkeypatch)
    ins = _run_verb(monkeypatch, tmp_path)
    index = tmp_path / "out" / "index"
    assert star_hostvirus.validate_index_manifest(index, verify_inputs=True) == []
    path = index / star_hostvirus.INDEX_MANIFEST
    m = json.loads(path.read_text())
    tamper(m, path, ins)
    if path.read_text() != "{":
        path.write_text(json.dumps(m))
    problems = star_hostvirus.validate_index_manifest(index, verify_inputs=True)
    assert any(expect in p for p in problems), problems


def test_write_index_manifest_rejects_wrong_input_keys(tmp_path):
    with pytest.raises(ValueError, match="exactly the keys"):
        star_hostvirus.write_index_manifest(
            tmp_path, {"host_fasta": "x"}, {"host": 1, "virus": 1}, ["STAR"], "2.7", 1, 1
        )


def _bam_from(tmp_path, recs, contigs):
    head = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{c}\tLN:1000\n" for c in contigs)
    sam = tmp_path / "y.sam"
    sam.write_text(head + "".join(recs))
    bam = tmp_path / "Aligned.sortedByCoord.out.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    return bam


def _strain_bam(tmp_path):
    """V1a and V1b are two accessions of virus V1; V2 is another virus; V3 is not in the map."""
    recs = [
        _rec("d1", 0, "VIRUS|V1a", 10, 2, 1, "AAAA", "U1"),  # same-virus pair: not "other virus"
        _rec("d1", 256, "VIRUS|V1b", 5, 2, 2, "AAAA", "U1"),
        _rec("d2", 0, "VIRUS|V1a", 20, 2, 1, "BBBB", "U2"),  # V1 primary, V2 secondary
        _rec("d2", 256, "VIRUS|V2", 20, 2, 2, "BBBB", "U2"),
        _rec("m1", 0, "VIRUS|V1a", 30, 1, 1, "-", "U9"),  # UB but no CB: not a molecule
        _rec("m2", 0, "VIRUS|V1a", 40, 1, 1, "GGGG", "-"),  # CB but no UB: not a molecule
        _rec("m3", 0, "VIRUS|V1a", 50, 1, 1, "HHHH", "U8"),
        _rec("lo", 256, "VIRUS|V1a", 6, 2, 2, "CCCC", "U3"),  # secondary sorts before its primary
        _rec("lo", 0, "VIRUS|V2", 500, 2, 1, "CCCC", "U3"),
        _rec("x", 0, "VIRUS|V3", 60, 1, 1, "DDDD", "U4"),
    ]
    return _bam_from(
        tmp_path, recs, ("HOST|chr1", "VIRUS|V1a", "VIRUS|V1b", "VIRUS|V2", "VIRUS|V3")
    )


@needs_samtools
def test_summarise_by_virus_counts_each_read_once_per_virus(tmp_path):
    rows = summarise_by_virus(_strain_bam(tmp_path), virus_of={"V1a": "V1", "V1b": "V1"})
    got = {r["accession"]: [r[c] for c in STAR_SUMMARY_COLUMNS[1:]] for r in rows}
    # columns: primary, any, molecules, host_competing, other_virus
    assert got == {
        "V1": [5, 6, 4, 0, 2],  # d1 counted once; other_virus = d2, lo (not d1)
        "V2": [1, 2, 2, 0, 2],
        "V3": [1, 1, 1, 0, 0],  # accession missing from the map keeps its own name
    }
    assert [r["accession"] for r in rows] == ["V1", "V2", "V3"]


@needs_samtools
def test_summarise_by_virus_per_accession_default(tmp_path):
    got = {r["accession"]: r for r in summarise_by_virus(_strain_bam(tmp_path))}
    assert got["V1a"]["star_other_virus"] == 3 and got["V1a"]["star_molecules"] == 4
    assert got["V1b"]["star_reads_any"] == 1


@needs_samtools
def test_align_host_unmapped_prune_writes_summary_before_pruning(tmp_path, monkeypatch):
    out = tmp_path / "o"
    src = _strain_bam(tmp_path)

    def fake_run_star(cmd, work, cb_len=12):
        shutil.copy(src, work / "Aligned.sortedByCoord.out.bam")

    monkeypatch.setattr(star_hostvirus, "run_star", fake_run_star)
    viral = align_host_unmapped(
        "r1", "r2", "idx", out, "10xv3", whitelist="wl.txt", prune=True,
        virus_of={"V1a": "V1", "V1b": "V1"},
    )  # fmt: skip
    assert viral.name == "Aligned.viral.bam" and viral.is_file()
    assert not (out / "Aligned.sortedByCoord.out.bam").exists()
    tsv = (out / "star_support_by_virus.tsv").read_text().splitlines()
    assert tsv[0].split("\t") == list(STAR_SUMMARY_COLUMNS)
    assert [ln.split("\t")[0] for ln in tsv[1:]] == ["V1", "V2", "V3"]
    assert tsv[1].split("\t") == ["V1", "5", "6", "4", "0", "2"]


FINAL_MULTI = "Uniquely mapped reads number |\t0\nNumber of reads mapped to multiple loci |\t5\n"


def test_star139_retry_that_maps_only_multimappers_is_refused(tmp_path, monkeypatch):
    _fake_run(monkeypatch, tmp_path, [139, 0], FINAL_MULTI)
    with pytest.raises(RuntimeError, match="mapped 5 reads"):
        star_run.run_star(["STAR", "--soloCBwhitelist", "None"], tmp_path)


def test_star139_retry_clears_the_crashed_runs_leftovers(tmp_path, monkeypatch):
    for leftover in ("_STARtmp", "_STARgenome", "_STARpass1", "Solo.out"):
        (tmp_path / leftover).mkdir()
    seen = []
    codes = iter([139, 0])

    def run(cmd, check=True):
        seen.append([(tmp_path / n).exists() for n in ("_STARtmp", "Solo.out")])
        if next(codes):
            raise subprocess.CalledProcessError(139, cmd)
        (tmp_path / "Log.final.out").write_text(FINAL_ZERO)

    monkeypatch.setattr(star_run.subprocess, "run", run)
    star_run.run_star(["STAR", "--soloCBwhitelist", "None"], tmp_path)
    assert seen == [[True, True], [False, False]]  # gone before the retry started

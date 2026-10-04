"""ANDET-09: anellovirus alignment branch — pure helpers on synthetic SAM."""

from __future__ import annotations

from pathlib import Path

import pytest

from viralscan import anello_align as aa
from viralscan.defaults import DEFAULTS
from viralscan.runconfig import RunConfig


def sam(qname, rname, pos, cigar, seq, nh=1, nm=0, flag=0, cb="AAAC", ub="GGGG"):
    tags = [f"NH:i:{nh}", f"NM:i:{nm}", f"CB:Z:{cb}", f"UB:Z:{ub}"]
    return "\t".join([qname, str(flag), rname, str(pos), "255", cigar, "*", "0", "0", seq, "*", *tags])


def parse(lines):
    return [a for a in (aa.parse_sam_line(x) for x in lines) if a is not None]


SEQ10 = "ACGTACGTAC"


def test_parse_skips_header_and_unmapped():
    assert aa.parse_sam_line("@HD\tVN:1.4") is None
    assert aa.parse_sam_line(sam("r", "*", 0, "*", SEQ10, flag=4)) is None
    a = aa.parse_sam_line(sam("r", "A1", 5, "10M", SEQ10, nh=3, nm=1))
    assert (a.rname, a.pos, a.nh, a.tags["CB"]) == ("A1", 5, 3, "AAAC")


def test_accession_metrics_ambiguity_breadth_identity():
    rows = parse(
        [
            sam("u1", "A1", 1, "10M", SEQ10, nm=1),  # unique, identity 0.9
            sam("m1", "A1", 11, "10M", SEQ10, nh=2),  # placed on A1 and A2
            sam("m1", "A2", 1, "10M", SEQ10, nh=2, flag=256),
            sam("s1", "A1", 1, "5M10N5M", "A" * 10, flag=16),  # spliced, reverse
            sam("h1", "A1", 1, "20M", "G" * 20),  # poly-G artefact read
        ]
    )
    got = {r["accession"]: r for r in aa.accession_metrics(rows, {"A1": 40, "A2": 40})}
    a1, a2 = got["A1"], got["A2"]
    assert a1["reads"] == 4 and a1["unique_reads"] == 3
    assert a2["reads"] == 1 and a2["unique_reads"] == 0
    assert a2["weighted_reads"] == 0.5
    # A2 is covered only by an NH=2 read, so unique-only breadth stays 0.
    assert a2["breadth"] == 0.25 and a2["breadth_unique"] == 0.0
    # A1: positions 1-20 (u1, m1, h1, s1 first block) plus 16-20 from the spliced read.
    assert a1["breadth"] == 0.5
    assert a1["splice_reads"] == 1
    assert a1["homopolymer_fraction"] == 0.25
    assert a1["sense_fraction"] == 0.75
    assert a1["median_identity"] == 1.0  # identities 0.9, 1, 1, 1


def test_virus_molecules_genus_unique_and_umi_dedup():
    acc_to_virus = {"A1": "Alpha", "A2": "Alpha", "B1": "Beta"}
    rows = parse(
        [
            sam("r1", "A1", 1, "10M", SEQ10, nh=2, ub="U1"),
            sam("r1", "A2", 1, "10M", SEQ10, nh=2, ub="U1", flag=256),  # same genus
            sam("r2", "A1", 1, "10M", SEQ10, ub="U1"),  # PCR duplicate of r1's molecule
            sam("r3", "A1", 1, "10M", SEQ10, nh=2, ub="U3"),
            sam("r3", "B1", 1, "10M", SEQ10, nh=2, ub="U3", flag=256),  # two genera
            sam("r4", "B1", 1, "10M", SEQ10, cb="-", ub="U4"),  # no valid barcode
            sam("r5", "B1", 1, "10M", SEQ10, cb="CCCC", ub="U5"),
        ]
    )
    got = aa.virus_molecules(rows, acc_to_virus)
    assert got["Alpha"] == {"molecules": 1, "cells": 1}
    assert got["Beta"] == {"molecules": 1, "cells": 1}


def test_virus_summary_aggregates_accessions():
    acc_rows = [
        {"accession": "A1", "reads": 3, "unique_reads": 2, "median_identity": 0.98,
         "start_sites": 2, "homopolymer_fraction": 0.0, "splice_reads": 1},
        {"accession": "A2", "reads": 1, "unique_reads": 0, "median_identity": 0.9,
         "start_sites": 1, "homopolymer_fraction": 1.0, "splice_reads": 0},
    ]
    got = aa.virus_summary(acc_rows, {"Alpha": {"molecules": 2, "cells": 1}},
                           {"A1": "Alpha", "A2": "Alpha"})["Alpha"]
    assert got["alignment_reads"] == 4
    assert got["alignment_accessions"] == 2
    assert got["alignment_molecules_unique"] == 2
    assert got["alignment_homopolymer_fraction"] == 0.25


BASE = {"virus_name": "", "viral_molecules_total_est": 5, "n_called_cells": 100,
        "total_cells": 900, "n_comparable_cells": 80}


def _row(name):
    return dict(BASE, virus_name=name)


def test_merge_adds_columns_and_alignment_only_rows():
    rows = [_row("EBV"), _row("Alphatorquevirus")]
    evidence = {
        "Alphatorquevirus": {"alignment_molecules_unique": "3", "alignment_reads": "9"},
        "Betatorquevirus": {"alignment_molecules_unique": "1", "alignment_reads": "2"},
        "Gammatorquevirus": {"alignment_molecules_unique": "0", "alignment_reads": "4"},
    }
    names = {"Alphatorquevirus", "Betatorquevirus", "Gammatorquevirus"}
    template = dict(BASE, viral_molecules_total_est=0)
    out = {r["virus_name"]: r for r in aa.merge_summary_rows(rows, evidence, names, aa.STATUS_OK, template)}
    assert out["EBV"]["detection_source"] == "kallisto"
    assert out["EBV"]["alignment_status"] == "" and out["EBV"]["alignment_reads"] == ""
    assert out["Alphatorquevirus"]["detection_source"] == "kallisto+alignment"
    assert out["Alphatorquevirus"]["alignment_reads"] == "9"
    # >= 1 unique molecule and no kallisto row -> a new row; 0 molecules -> none.
    assert out["Betatorquevirus"]["detection_source"] == "alignment_only"
    assert out["Betatorquevirus"]["viral_molecules_total_est"] == 0
    assert out["Betatorquevirus"]["n_called_cells"] == 100
    assert "Gammatorquevirus" not in out


def test_merge_records_skip_status_without_new_rows():
    rows = [_row("Alphatorquevirus")]
    evidence = {"Betatorquevirus": {"alignment_molecules_unique": "5"}}
    out = aa.merge_summary_rows(rows, evidence, {"Alphatorquevirus", "Betatorquevirus"},
                                aa.STATUS_NO_HOST_FILTER, {})
    assert len(out) == 1
    assert out[0]["alignment_status"] == "skipped_no_host_filter"
    assert out[0]["detection_source"] == "kallisto"
    assert out[0]["alignment_reads"] == ""


def test_merge_ran_but_no_evidence_reports_zero():
    out = aa.merge_summary_rows([_row("Alphatorquevirus")], {}, {"Alphatorquevirus"},
                                aa.STATUS_OK, {})
    assert out[0]["alignment_reads"] == 0
    assert out[0]["alignment_status"] == "ok"


@pytest.mark.parametrize(
    "enabled,host,index,want",
    [
        (False, "h", "i", "disabled"),
        (True, None, "i", "skipped_no_host_filter"),
        (True, "h", None, "skipped_no_anello_index"),
        (True, "h", "i", "ok"),
    ],
)
def test_status_for(enabled, host, index, want):
    assert aa.status_for(enabled, host, index) == want


def test_write_star_reference_selects_anellovirus_records(tmp_path: Path):
    fa = tmp_path / "viral.fa"
    fa.write_text(">AB1.1 TTV\nACGT\nACGT\n>NC_9.1 EBV\nGGGG\n>KP2.1 TTMV\nTTTT\n")
    fasta, gtf, n = aa.write_star_reference(fa, {"AB1.1", "KP2"}, tmp_path / "idx")
    assert n == 2
    assert [name for name, _ in aa.iter_fasta(fasta)] == ["AB1.1", "KP2.1"]
    exons = [ln.split("\t") for ln in gtf.read_text().splitlines() if "\texon\t" in ln]
    assert [(e[0], e[4]) for e in exons] == [("AB1.1", "8"), ("KP2.1", "4")]
    with pytest.raises(ValueError):
        aa.write_star_reference(fa, {"NOPE"}, tmp_path / "idx2")


def test_align_cmd_reads_cdna_first_and_caps_match():
    cmd = aa.align_cmd("STAR", "/idx", "R2.fastq.gz", "R1.fastq.gz", ["--soloType", "X"], "p/", 2)
    i = cmd.index("--readFilesIn")
    assert cmd[i + 1 : i + 3] == ["R2.fastq.gz", "R1.fastq.gz"]
    assert cmd[cmd.index("--readFilesCommand") + 1] == "zcat"
    filt = cmd[cmd.index("--outFilterMultimapNmax") + 1]
    assert filt == cmd[cmd.index("--outSAMmultNmax") + 1] == "100"
    assert cmd[cmd.index("--soloStrand") + 1] == "Unstranded"


def test_resolve_index_and_runconfig(tmp_path: Path):
    kb_index = tmp_path / "panel.idx"
    kb_index.write_text("")
    cfg = {"output": str(tmp_path / "o") + "/", "index": str(kb_index), "transcripts": "t",
           "sample1": "s1", "sample2": "s2", "technology": "10xv3", "visual": "True",
           "multimapping": "True", "gtf": "None", "fasta": "None", "f1": "None",
           "reference": "False", "umap": "False", "whitelist": "None", "emptydrops_seed": 100}
    # No anello_star/ beside the index: nothing to resolve, whatever the default.
    assert aa.resolve_index(str(kb_index)) is None
    rc = RunConfig.from_snakemake_config({**cfg, "anello_align": "true"})
    assert rc.anello_align is True and rc.anello_index is None

    (tmp_path / "anello_star").mkdir()
    (tmp_path / "anello_star" / "SA").write_text("")
    rc = RunConfig.from_snakemake_config({**cfg, "anello_align": "true"})
    assert rc.anello_index == str((tmp_path / "anello_star").resolve())
    # Disabled, so the index is not resolved even though it is present.
    assert RunConfig.from_snakemake_config({**cfg, "anello_align": "false"}).anello_index is None
    # The shipped default decides what an unset config does (ANDET-09e: off).
    assert RunConfig.from_snakemake_config(cfg).anello_align is DEFAULTS["anello_align"]


# ── Read-side artefact measures (F-019 update, 2026-10-04) ───────────────────
# The point of these: the covid TTV reads and a genuine 3' end read are
# indistinguishable by pileup, poly-A content or NM. They must be separable by
# body complexity and TSO content, and genuine reads must never be flagged.
_TSO_RC = "GTACTCTGCGTTGATACCACTGCTT"


def _viral_body(n=40, seed=3):
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_a_genuine_3p_read_has_a_complex_body_despite_its_polya_tail():
    """[complex viral 3'UTR][untemplated poly-A] — the shape a real read has."""
    read = _viral_body(40) + "A" * 50
    assert aa.read_body(read) == _viral_body(40)
    assert aa.is_complex_body(read) is True
    assert aa.has_tso(read) is False


def test_the_covid_chimera_shape_has_no_complex_body():
    """[poly-A][TSO-rc] — what 20/30 of the covid reads actually look like."""
    read = "A" * 65 + _TSO_RC
    assert aa.is_complex_body(read) is False
    assert aa.has_tso(read) is True


def test_a_complex_body_is_not_rescued_by_carrying_the_tso():
    """Both measures are reported; neither overrides the other."""
    read = _viral_body(40) + "A" * 25 + _TSO_RC
    assert aa.is_complex_body(read) is True   # the body is genuinely complex
    assert aa.has_tso(read) is True           # but it cannot be a real molecule


def test_the_artefact_classes_all_fail_the_complexity_gate():
    for seq in ("A" * 90, "G" * 90, ("CAG" * 30)[:90], ("AC" * 45)):
        assert aa.is_complex_body(seq) is False, seq[:12]


@pytest.mark.parametrize("body_len", [25, 30, 40, 60])
def test_no_genuine_body_is_discarded_at_any_length(body_len):
    """The permissive setting: a real 3' end read must never read as artefact.

    These bodies are random ACGT, which is the hardest honest case — a real
    3'UTR is not more repetitive than chance. Pinning this at several lengths
    catches a future threshold raise that would silently cost sensitivity.
    """
    import random

    rng = random.Random(11)
    bodies = ["".join(rng.choice("ACGT") for _ in range(body_len)) for _ in range(400)]
    reads = [b + "A" * (90 - body_len) for b in bodies]
    kept = sum(aa.is_complex_body(r) for r in reads)
    assert kept == len(reads), f"{len(reads) - kept} genuine {body_len} nt bodies discarded"


def test_the_threshold_sits_in_the_measured_gap():
    """Artefact classes top out at 1.70 bits; genuine 20 nt bodies floor at 2.21."""
    worst_artefact = max(
        aa.dinucleotide_entropy(s)
        for s in ("A" * 40, "G" * 40, ("CAG" * 14)[:40], "AC" * 20, "AAAAC" * 8)
    )
    assert worst_artefact < aa.MIN_BODY_ENTROPY <= 2.21


def test_the_tso_is_found_in_both_orientations_and_tolerates_two_mismatches():
    assert aa.has_tso("CCCC" + aa.TSO + "CCCC") is True
    assert aa.has_tso("CCCC" + _TSO_RC + "CCCC") is True
    mutated = list(aa.TSO)
    mutated[3] = "T" if mutated[3] != "T" else "A"
    mutated[17] = "T" if mutated[17] != "T" else "A"
    assert aa.has_tso("".join(mutated)) is True
    assert aa.has_tso(_viral_body(60, seed=9)) is False


def test_query_coverage_exposes_a_soft_clipped_perfect_match():
    """identity 1.0 on 34 of 90 nt is what made the covid NM=0 misleading."""
    line = "\t".join(
        ["r1", "0", "MZ286238.1", "2830", "3", "34M56S",
         "*", "0", "0", "A" * 34 + "C" * 56, "I" * 90, "NM:i:0", "NH:i:1"]
    )
    a = aa.parse_sam_line(line)
    assert a.identity() == 1.0
    assert a.query_coverage() == pytest.approx(34 / 90, abs=1e-4)

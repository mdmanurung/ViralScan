"""ANDET-09: anellovirus alignment branch — pure helpers on synthetic SAM."""

from __future__ import annotations

from pathlib import Path

import pytest

from viralscan import anello_align as aa
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
    assert aa.resolve_index(str(kb_index)) is None
    rc = RunConfig.from_snakemake_config(cfg)
    assert rc.anello_align is True and rc.anello_index is None

    (tmp_path / "anello_star").mkdir()
    (tmp_path / "anello_star" / "SA").write_text("")
    rc = RunConfig.from_snakemake_config(cfg)
    assert rc.anello_index == str((tmp_path / "anello_star").resolve())
    assert RunConfig.from_snakemake_config({**cfg, "anello_align": "false"}).anello_index is None

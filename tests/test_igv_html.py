"""EVID-04: ``viralscan evidence --html`` writes one self-contained igv.js page.

The fixtures are the synthetic evidence directory of ``test_splice``. There is no browser in CI: the "render
with networking off" test is static (the page, moved to another directory, holds no network target outside the
vendored library and the licence text). A headless Firefox run was done by hand and is recorded in PLAN.
"""

import argparse
import hashlib
import html
import json
import re
import shutil
from pathlib import Path

import pytest

from tests.test_splice import GENOME, _evidence_dir, _fasta, _gtf
from viralscan import igv_html, splice

needs_tools = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and samtools are required",
)

#: Pinned vendor files (npm igv@2.15.11 dist/igv.min.js and LICENSE). A changed file fails here.
IGV_MIN_JS_SHA256 = "7ea7c3a666f02694900e3eccc6d662f7881cb90f9e3c889a5ab60a2729b54434"
IGV_LICENSE_SHA256 = "194adb0bf01eedc74fb6def807e2b0047bd80521f74f87fba068ed6f5c24a87e"


def _sha(name: str) -> str:
    return hashlib.sha256(igv_html.vendor_text(name).encode()).hexdigest()


def _receipt(page: str) -> dict:
    match = re.search(
        r'<script id="viralscan-receipt" type="application/json">(.*?)</script>', page, re.S
    )
    assert match
    return json.loads(match.group(1))


def _data(page: str) -> dict:
    match = re.search(r'<script id="igv-data" type="application/json">(.*?)</script>', page, re.S)
    assert match
    return json.loads(match.group(1))


def _own_markup(page: str) -> str:
    """The page without the vendored library and the licence text."""
    page = re.sub(r'<script id="igv-lib">.*?</script>', "", page, flags=re.S)
    return re.sub(r"<pre>.*?</pre>", "", page, flags=re.S)


@pytest.fixture
def built(tmp_path: Path) -> tuple[Path, Path, Path]:
    """(evidence dir, viral fasta, gtf) with a finished --splice run."""
    ev = _evidence_dir(tmp_path)
    fasta, gtf = _fasta(tmp_path), _gtf(tmp_path, "exon")
    splice.run_splice(ev, viral_fasta=fasta, annotation_gtfs=[gtf])
    return ev, Path(fasta), gtf


# --- vendored library ------------------------------------------------------------------------------


def test_vendored_igv_is_the_pinned_release() -> None:
    assert _sha("igv.min.js") == IGV_MIN_JS_SHA256
    assert _sha("LICENSE") == IGV_LICENSE_SHA256
    assert igv_html.IGV_VERSION == "2.15.11"
    assert "The MIT License" in igv_html.vendor_text("LICENSE")
    assert "2.15.11" in igv_html.vendor_text("igv.min.js")


def test_vendor_files_are_package_data() -> None:
    import tomllib

    root = Path(__file__).resolve().parents[1]
    data = tomllib.loads((root / "pyproject.toml").read_text())["tool"]["setuptools"][
        "package-data"
    ]
    assert {"vendor/igv/igv.min.js", "vendor/igv/LICENSE", "templates/*.html"} <= set(
        data["viralscan"]
    )


# --- the page --------------------------------------------------------------------------------------


@needs_tools
def test_page_is_self_contained_after_moving_it(
    built: tuple[Path, Path, Path], tmp_path: Path
) -> None:
    ev, fasta, gtf = built
    path = igv_html.write_igv_page(ev, viral_fasta=fasta, annotation_gtfs=[gtf])
    assert path == ev / "igv_browser.html"
    moved = tmp_path / "elsewhere" / "copy.html"
    moved.parent.mkdir()
    shutil.move(path, moved)
    page = moved.read_text()
    own = _own_markup(page)
    # no network target in anything but the library and the licence text
    assert not re.search(r"https?://|//cdn|src=|href=|@import|url\(\s*[^d)]", own), own[:2000]
    assert "loadDefaultGenomes: false" in own
    assert "connect-src data: blob:" in html.unescape(
        own
    ) and "default-src 'none'" in html.unescape(own)
    data = _data(page)
    for key in ("fasta", "fai", "bam", "bai", "splice_bam", "splice_bai", "bed"):
        assert data[key].startswith("data:application/octet-stream;base64,")
    assert data["contigs"] == ["chrV"]
    assert "showCoverage: false" in own  # igv.js would draw a depth track computed from the sample


@needs_tools
def test_page_states_sampling_rule_omissions_and_shows_no_metric(
    built: tuple[Path, Path, Path],
) -> None:
    ev, fasta, gtf = built
    page = igv_html.write_igv_page(ev, viral_fasta=fasta, annotation_gtfs=[gtf]).read_text()
    receipt = _receipt(page)
    assert f"fixed seed {igv_html.SEED}" in page
    assert f"at most {igv_html.MAX_READS} of those reads" in page
    assert "longer than 20000 bp left out" in page
    assert "count nothing here" in page
    assert "splice/splice.bam" in page
    # population (3 virus-best molecules, 5 member reads) is stated from the full tables
    assert receipt["category_counts"]["counted_virus_best"] == 3
    assert receipt["member_reads"] == 5
    assert receipt["member_reads_primary_on_host"] == 0
    assert receipt["contigs"] == [
        {
            "name": "chrV",
            "length": len(GENOME),
            "reads_available": 5,
            "reads_drawn": 5,
            "splice_available": 3,
            "splice_drawn": 3,
        }
    ]
    assert receipt["igv_js"]["sha256"] == IGV_MIN_JS_SHA256
    assert receipt["inputs"]["competitive_reads.raw.bam"]
    # only counts of the population and of the draw: no rate, fraction or percentage in the page text
    text = re.sub(r"<[^>]+>", " ", _own_markup(page))
    assert not re.search(r"\d\s*%|fraction of|per cent|mean |median ", text)


@needs_tools
def test_same_inputs_give_the_same_bytes(built: tuple[Path, Path, Path], tmp_path: Path) -> None:
    ev, fasta, gtf = built
    a = igv_html.write_igv_page(
        ev, viral_fasta=fasta, annotation_gtfs=[gtf], out_dir=tmp_path / "a"
    )
    b = igv_html.write_igv_page(
        ev, viral_fasta=fasta, annotation_gtfs=[gtf], out_dir=tmp_path / "b"
    )
    assert a.read_bytes() == b.read_bytes()


@needs_tools
def test_genome_falls_back_to_the_splice_target_genome_and_needs_one(
    built: tuple[Path, Path, Path], tmp_path: Path
) -> None:
    ev, _, _ = built
    page = igv_html.write_igv_page(ev).read_text()  # splice/target_genome.fa
    assert "No gene track" in page  # no GTF given
    assert "bed" not in _data(page)
    (ev / "splice" / "target_genome.fa").unlink()
    with pytest.raises(igv_html.IgvError, match="viral genome"):
        igv_html.write_igv_page(ev)


@needs_tools
def test_no_splice_run_gives_no_splice_track(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path)
    page = igv_html.write_igv_page(ev, viral_fasta=_fasta(tmp_path)).read_text()
    assert "splice_bam" not in _data(page)
    assert "No splice track" in page


# --- sampling and the size cap ---------------------------------------------------------------------


def _many(n: int) -> list:
    """*n* single-read counted virus-best molecules, every read 60 bp of the genome."""
    return [
        ((f"CB{i}", f"U{i}"), True, "virus_best", [(i + 1, 0, GENOME[(i * 7) % 400 :][:60])])
        for i in range(n)
    ]


@needs_tools
def test_read_sample_is_capped_seeded_and_reported(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path, _many(40))
    fasta = _fasta(tmp_path)
    page_a = igv_html.write_igv_page(ev, viral_fasta=fasta, max_reads=10, out_dir=tmp_path / "a")
    page_b = igv_html.write_igv_page(ev, viral_fasta=fasta, max_reads=10, out_dir=tmp_path / "b")
    contig = _receipt(page_a.read_text())["contigs"][0]
    assert (contig["reads_available"], contig["reads_drawn"]) == (40, 10)
    assert page_a.read_bytes() == page_b.read_bytes()
    assert "at most 10 of those reads" in page_a.read_text()
    # a different cap is a different, still uniform, draw
    other = igv_html.write_igv_page(ev, viral_fasta=fasta, max_reads=20, out_dir=tmp_path / "c")
    assert _receipt(other.read_text())["contigs"][0]["reads_drawn"] == 20


@needs_tools
def test_size_cap_halves_the_sample_and_says_so(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path, _many(400))
    fasta = _fasta(tmp_path)
    full = igv_html.write_igv_page(ev, viral_fasta=fasta, max_reads=400, out_dir=tmp_path / "full")
    size = full.stat().st_size
    capped = igv_html.write_igv_page(
        ev,
        viral_fasta=fasta,
        max_reads=400,
        max_html_bytes=size - 3000,
        out_dir=tmp_path / "capped",
    )
    page = capped.read_text()
    receipt = _receipt(page)
    assert capped.stat().st_size <= size - 3000
    assert receipt["reduced_for_size"] is True
    assert receipt["read_cap_used"] < 400
    contig = receipt["contigs"][0]
    assert contig["reads_available"] == 400 and contig["reads_drawn"] == receipt["read_cap_used"]
    assert "reduced from 400 to" in page
    assert _receipt(full.read_text())["reduced_for_size"] is False


@needs_tools
def test_page_that_cannot_fit_is_refused_not_truncated(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path)
    with pytest.raises(igv_html.IgvError, match="above the"):
        igv_html.write_igv_page(ev, viral_fasta=_fasta(tmp_path), max_html_bytes=1000)


@needs_tools
def test_no_unit_is_not_evaluable(tmp_path: Path) -> None:
    ev = _evidence_dir(tmp_path)
    (ev / "molecule_verdicts.tsv.gz").unlink()
    with pytest.raises(igv_html.IgvNotEvaluable):
        igv_html.write_igv_page(ev, viral_fasta=_fasta(tmp_path))
    (tmp_path / "h").mkdir()
    only_host = _evidence_dir(
        tmp_path / "h", [(("CB", "U"), True, "host_best", [(1, 0, GENOME[:60])])]
    )
    with pytest.raises(igv_html.IgvNotEvaluable, match="no counted virus-best"):
        igv_html.write_igv_page(only_host, viral_fasta=_fasta(tmp_path))


# --- gene track ------------------------------------------------------------------------------------


def test_bed12_rules(tmp_path: Path) -> None:
    def row(seq: str, kind: str, a: int, b: int, attrs: str) -> str:
        return f"{seq}\tt\t{kind}\t{a}\t{b}\t.\t+\t.\t{attrs}"

    gtf = tmp_path / "g.gtf"
    gtf.write_text(
        "\n".join(
            [
                "# comment",
                # two-exon transcript: gene row skipped, one row per transcript, blocks from exons
                row("VIRUS|v", "gene", 1, 300, 'gene_id "A"; gene_name "geneA";'),
                row(
                    "v", "transcript", 1, 300, 'gene_id "A"; transcript_id "TA"; gene_name "geneA";'
                ),
                row("v", "exon", 1, 100, 'gene_id "A"; transcript_id "TA"; gene_name "geneA";'),
                row("v", "exon", 201, 300, 'gene_id "A"; transcript_id "TA"; gene_name "geneA";'),
                # CDS-only gene in both repeats, 50 kb apart: two copies, each reaches its gene-row UTR ends
                row("v", "gene", 1001, 1500, 'gene_id "B"; gene_name "geneB";'),
                row("v", "CDS", 1101, 1400, 'gene_id "B"; gene_name "geneB";'),
                row("v", "gene", 51001, 51500, 'gene_id "B"; gene_name "geneB";'),
                row("v", "CDS", 51101, 51400, 'gene_id "B"; gene_name "geneB";'),
                # a feature over MAX_FEATURE_BP: left out and counted
                row("v", "gene", 100, 60000, 'gene_id "L"; gene_name "long";'),
                # another contig: never read
                row("other", "gene", 1, 50, 'gene_id "X"; gene_name "x";'),
            ]
        )
        + "\n"
    )
    text, left_out, n_rows = igv_html.bed12([gtf], ["v"])
    lines = [line.split("\t") for line in text.splitlines()]
    assert [(f[0], int(f[1]), int(f[2]), f[3]) for f in lines] == [
        ("v", 0, 300, "geneA"),
        ("v", 1000, 1500, "geneB"),
        ("v", 51000, 51500, "geneB"),
    ]
    a = lines[0]
    assert (a[9], a[10], a[11]) == ("2", "100,100", "0,200")  # BED12 blocks from the two exons
    b = lines[1]
    assert (b[6], b[7]) == ("1100", "1400")  # thick part is the CDS; the block spans the UTRs
    assert left_out == [("v", "long", 59901)]
    assert n_rows == 9  # every gene/transcript/exon/CDS row of "v" is read


def test_fasta_and_index_are_consistent() -> None:
    fasta, fai = igv_html.fasta_and_index({"a": "ACGT" * 40, "b": "TT"})
    lines = fai.decode().splitlines()
    name, length, offset, bases, width = lines[0].split("\t")
    assert (name, int(length), int(bases), int(width)) == ("a", 160, 60, 61)
    assert fasta[int(offset) : int(offset) + 4] == b"ACGT"
    off_b = int(lines[1].split("\t")[2])
    assert fasta[off_b : off_b + 2] == b"TT"


# --- wiring ----------------------------------------------------------------------------------------


def test_html_flag_is_off_by_default_on_both_verbs() -> None:
    from viralscan.cli.parsers import parse_args

    base = ["evidence", "--run-dir", "r", "--virus", "EBV"]
    assert parse_args(base).html is False
    assert parse_args([*base, "--html"]).html is True
    assert parse_args(["evidence-plot", "--evidence-dir", "e"]).html is False
    ns = parse_args(["evidence-plot", "--evidence-dir", "e", "--html", "--viral-fasta", "v.fa"])
    assert ns.html is True and ns.viral_fasta == "v.fa"


def test_html_without_viral_fasta_is_refused(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    from viralscan.scripts.evidence_run import run_evidence

    args = argparse.Namespace(
        run_dir=str(tmp_path), virus="EBV", html=True, viral_fasta=None, review_candidates=False
    )
    with pytest.raises(SystemExit):
        run_evidence(args)
    assert "--html requires --viral-fasta" in caplog.text


@needs_tools
def test_evidence_plot_html_on_an_existing_directory(
    built: tuple[Path, Path, Path], tmp_path: Path
) -> None:
    from viralscan.competitive_plot import run_evidence_plot

    ev, fasta, gtf = built
    args = argparse.Namespace(
        evidence_dir=[str(ev)],
        output=str(tmp_path / "plots"),
        html=True,
        viral_fasta=str(fasta),
        annotation=str(gtf),
        cross_call=None,
        coverage_map=False,
        target=None,
    )
    run_evidence_plot(args)
    assert (tmp_path / "plots" / ev.name / "igv_browser.html").is_file()
    assert not list(ev.glob("competitive_alignment.*"))  # --html replaces the figures
    args.target = "chrV"
    with pytest.raises(ValueError, match="--html cannot be combined"):
        run_evidence_plot(args)


@needs_tools
def test_evidence_manifest_lists_the_page(built: tuple[Path, Path, Path]) -> None:
    from viralscan.scripts.evidence_run import _write_evidence_manifest

    ev, fasta, gtf = built
    igv_html.write_igv_page(ev, viral_fasta=fasta, annotation_gtfs=[gtf])
    replay = {"strand": None, "technology": "10xv3", "kallisto_technology": "10XV3"}
    path = _write_evidence_manifest(
        ev,
        target_label="T",
        target_genes=["g"],
        method="m",
        run_fingerprint="f",
        references={},
        replay=replay,
    )
    assert "igv_browser.html" in json.loads(path.read_text())["outputs"]

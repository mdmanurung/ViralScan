"""PROG-13: merge 31-mer-identical repeat copies to one gene in the t2g."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

from viralscan.repeat_merges import (
    apply_repeat_merges,
    find_t2g,
    load_repeat_merges,
)
from viralscan.sensitivity import DEFAULT_K

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA = REPO_ROOT / "src" / "viralscan" / "data"
NCBI_CACHE = Path.home() / ".cache" / "viralscan" / "ncbi"

#: (gtf, accession, gene-id prefix, [(copy, survivor)]) -- the pairs PROG-13 merges.
PAIRS = [
    (
        "Human_herpesvirus_1_NC_001806.gtf",
        "NC_001806.2",
        "HUM_HERP1_",
        [
            ("HHV1gp00s01", "HHV1gp00s02"),  # LAT
            ("HHV1gp00p76", "HHV1gp00p17"),  # RL2 / ICP0
            ("HHV1gp00p01", "HHV1gp00p15"),  # RS1 / ICP4
        ],
    ),
    (
        "Varicella_zoster_virus_NC_001348.gtf",
        "NC_001348.1",
        "VARICELLA_HHV3_",
        [("gp72", "gp63"), ("gp71", "gp64")],  # ORF62, ORF63
    ),
]


def _t2g_text() -> str:
    return (
        "ENST00000632585.1\tENSG00000282172.1\t\t\tENST00000632585.1\t1\t408\t+\n"
        "NC_001806.2-unassigned_transcript_1\tHHV1gp00s01\t\t\tNC_001806.2\t1\t7569\t-\n"
        "NC_001806.2-unassigned_transcript_62\tHHV1gp00s02\t\t\tNC_001806.2\t118805\t127151\t+\n"
        "HHV3_gp72\tHHV3_gp72\t\t\tNC_001348.1\t120764\t124756\t+\n"
        "HHV3_gp63\tHHV3_gp63\t\t\tNC_001348.1\t105141\t109133\t-\n"
        "HHV3_gp05\tHHV3_gp05\t\t\tNC_001348.1\t1\t2\t+\n"
    )


class TestTable:
    def test_loads_and_is_idempotent_by_construction(self) -> None:
        merges = load_repeat_merges()
        assert merges, "the packaged table is empty"
        assert not set(merges) & set(merges.values()), "a survivor is also a copy"

    def test_lists_every_prog13_pair_in_both_spellings(self) -> None:
        merges = load_repeat_merges()
        for _, _, prefix, pairs in PAIRS:
            for copy, survivor in pairs:
                assert merges[prefix + copy] == prefix + survivor
                bare_copy = copy
                bare_survivor = survivor
                if prefix == "VARICELLA_HHV3_":
                    bare_copy, bare_survivor = "HHV3_" + copy, "HHV3_" + survivor
                assert merges[bare_copy] == bare_survivor

    def test_rejects_a_chain(self, tmp_path: Path) -> None:
        bad = tmp_path / "m.tsv"
        bad.write_text("gene_id\tmerged_into\tvirus\trationale\na\tb\tv\tr\nb\tc\tv\tr\n")
        with pytest.raises(ValueError, match="both a copy and a survivor"):
            load_repeat_merges(bad)

    def test_rejects_wrong_columns(self, tmp_path: Path) -> None:
        bad = tmp_path / "m.tsv"
        bad.write_text("gene\tinto\na\tb\n")
        with pytest.raises(ValueError, match="expected columns"):
            load_repeat_merges(bad)


class TestRewrite:
    MERGES = {"HHV1gp00s01": "HHV1gp00s02", "HHV3_gp72": "HHV3_gp63"}

    def test_rewrites_only_the_gene_column_of_listed_copies(self, tmp_path: Path) -> None:
        path = tmp_path / "t2g.txt"
        path.write_text(_t2g_text())
        assert apply_repeat_merges(path, self.MERGES) == 2
        lines = path.read_text().splitlines()
        original = _t2g_text().splitlines()
        changed = [i for i, (a, b) in enumerate(zip(original, lines)) if a != b]
        assert changed == [1, 3]
        assert lines[1] == original[1].replace("\tHHV1gp00s01\t", "\tHHV1gp00s02\t")
        # transcript (col 1) keeps its name; host rows and the other virus are untouched
        assert lines[1].split("\t")[0] == "NC_001806.2-unassigned_transcript_1"
        assert lines[3].split("\t")[:2] == ["HHV3_gp72", "HHV3_gp63"]

    def test_second_application_changes_nothing(self, tmp_path: Path) -> None:
        path = tmp_path / "t2g.txt"
        path.write_text(_t2g_text())
        apply_repeat_merges(path, self.MERGES)
        once = path.read_bytes()
        assert apply_repeat_merges(path, self.MERGES) == 0
        assert path.read_bytes() == once

    def test_unlisted_index_is_left_byte_identical(self, tmp_path: Path) -> None:
        path = tmp_path / "t2g.txt"
        path.write_bytes(b"ENST1\tENSG1\t\t\tENST1\r\nHHV3_gp05\tHHV3_gp05")  # CRLF, no final \n
        assert apply_repeat_merges(path, self.MERGES) == 0
        assert path.read_bytes() == b"ENST1\tENSG1\t\t\tENST1\r\nHHV3_gp05\tHHV3_gp05"

    def test_preserves_line_endings_of_rewritten_rows(self, tmp_path: Path) -> None:
        path = tmp_path / "t2g.txt"
        path.write_bytes(b"tx1\tHHV1gp00s01\t\t\tNC_001806.2\r\nENST1\tENSG1\t\t\tENST1\r\n")
        assert apply_repeat_merges(path, self.MERGES) == 1
        assert (
            path.read_bytes() == b"tx1\tHHV1gp00s02\t\t\tNC_001806.2\r\nENST1\tENSG1\t\t\tENST1\r\n"
        )

    def test_packaged_table_merges_both_spellings(self, tmp_path: Path) -> None:
        path = tmp_path / "t2g.txt"
        path.write_text(
            "t1\tHUM_HERP1_HHV1gp00s01\t\t\tNC_001806.2\nt2\tHHV3_gp71\t\t\tNC_001348.1\n"
        )
        assert apply_repeat_merges(path) == 2
        assert path.read_text().splitlines()[0].split("\t")[1] == "HUM_HERP1_HHV1gp00s02"
        assert path.read_text().splitlines()[1].split("\t")[1] == "HHV3_gp64"

    def test_find_t2g_prefers_t2g_txt_then_panel_t2g(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            find_t2g(tmp_path)
        (tmp_path / "panel.t2g").write_text("")
        assert find_t2g(tmp_path).name == "panel.t2g"
        (tmp_path / "t2g.txt").write_text("")
        assert find_t2g(tmp_path).name == "t2g.txt"


class TestCli:
    def _script(self):
        sys.path.insert(0, str(REPO_ROOT / "scripts"))
        import apply_repeat_gene_merges

        return apply_repeat_gene_merges

    def test_standalone_cli_rewrites_an_index_dir(self, tmp_path: Path, capsys) -> None:
        (tmp_path / "panel.t2g").write_text(_t2g_text())
        assert self._script().main(["--index-dir", str(tmp_path)]) == 0
        assert "2 row(s) rewritten" in capsys.readouterr().out
        assert self._script().main(["--index-dir", str(tmp_path)]) == 0
        assert "nothing to do" in capsys.readouterr().out

    def test_standalone_cli_fails_without_a_t2g(self, tmp_path: Path) -> None:
        with pytest.raises(SystemExit) as exc:
            self._script().main(["--index-dir", str(tmp_path)])
        assert "no t2g.txt or panel.t2g" in str(exc.value)

    def test_builder_hooks_the_rewrite_after_kb_ref(self) -> None:
        source = (REPO_ROOT / "scripts" / "build_bundled_panel_ref.py").read_text()
        assert source.index("subprocess.run(cmd, check=True)") < source.index(
            "apply_repeat_merges("
        )


# ---------------------------------------------------------------- 31-mer evidence

_COMP = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _canonical_kmers(seq: str, k: int) -> set[str]:
    out = set()
    for i in range(len(seq) - k + 1):
        word = seq[i : i + k]
        if set(word) <= set("ACGT"):
            out.add(min(word, _revcomp(word)))
    return out


def _gene_sequence(gtf: Path, genome: str, gene_id: str) -> str:
    """Spliced, strand-aware sequence; exon rows, else CDS, else the gene span."""
    rows: dict[str, list[tuple[int, int, str]]] = {"exon": [], "CDS": [], "gene": []}
    for line in gtf.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 9 or cols[2] not in rows:
            continue
        match = re.search(r'gene_id "([^"]*)"', cols[8])
        if match and match.group(1) == gene_id:
            rows[cols[2]].append((int(cols[3]), int(cols[4]), cols[6]))
    blocks = rows["exon"] or rows["CDS"] or rows["gene"]
    assert blocks, f"{gene_id} not in {gtf.name}"
    seq = "".join(genome[a - 1 : b] for a, b, _ in sorted(blocks))
    return _revcomp(seq) if blocks[0][2] == "-" else seq


def _read_fasta(path: Path) -> str:
    return "".join(ln.strip().upper() for ln in path.read_text().splitlines() if ln[:1] != ">")


@pytest.mark.parametrize("gtf_name, accession, prefix, pairs", PAIRS)
def test_repeat_copies_share_their_31mers(gtf_name, accession, prefix, pairs) -> None:
    """The evidence for merging: the smaller copy's 31-mers all occur in the larger copy.

    Canonical k-mers (kallisto's), of the strand-aware sequence, since the copies sit on
    opposite strands. VZV's copies are identical; HSV-1's differ by the tail the longer copy
    adds (LAT +407, ICP0 +27, ICP4 +1), so containment, not equality, is the invariant.
    """
    gtf = DATA / gtf_name
    fasta = NCBI_CACHE / accession / f"{accession}.fasta"
    if not gtf.is_file() or not fasta.is_file():
        pytest.skip("bundled viral GTFs are gitignored / NCBI cache absent in this checkout")
    genome = _read_fasta(fasta)
    merges = load_repeat_merges()
    for copy, survivor in pairs:
        kept = _canonical_kmers(_gene_sequence(gtf, genome, prefix + survivor), DEFAULT_K)
        gone = _canonical_kmers(_gene_sequence(gtf, genome, prefix + copy), DEFAULT_K)
        assert gone <= kept, f"{copy}: {len(gone - kept)} 31-mers are not in {survivor}"
        # the survivor is the copy with at least as many k-mers, so no unique sequence is lost
        assert len(kept) >= len(gone)
        assert merges[prefix + copy] == prefix + survivor

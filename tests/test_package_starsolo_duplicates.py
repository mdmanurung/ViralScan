"""REF-13 residue: duplicate FASTA IDs must fail before samtools/STAR see them."""

import importlib.util
from pathlib import Path

import pytest

from viralscan.evidence import write_competitive_fasta

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "package_starsolo_viral_references.py"
_spec = importlib.util.spec_from_file_location("package_starsolo_viral_references", SCRIPT)
pkg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pkg)

GTF = 'NC_002076.2\tx\texon\t1\t4\t.\t+\t.\tgene_id "g";\n'


def test_validate_pair_rejects_duplicate_fasta_ids(tmp_path):
    fa = tmp_path / "v.fa"
    fa.write_text(">NC_002076.2 a\nACGT\n>NC_002076.2 b\nACGT\n")
    gtf = tmp_path / "v.gtf"
    gtf.write_text(GTF)
    with pytest.raises(SystemExit, match="NC_002076.2"):
        pkg.validate_pair(fa, gtf)


def test_validate_pair_accepts_unique_ids(tmp_path):
    fa = tmp_path / "v.fa"
    fa.write_text(">NC_002076.2\nACGT\n")
    gtf = tmp_path / "v.gtf"
    gtf.write_text(GTF)
    pkg.validate_pair(fa, gtf)


def test_competitive_fasta_rejects_duplicate_ids(tmp_path):
    host = tmp_path / "h.fa"
    host.write_text(">chr1\nACGT\n")
    virus = tmp_path / "v.fa"
    virus.write_text(">NC_002076.2\nACGT\n>NC_002076.2 dup\nACGT\n")
    with pytest.raises(ValueError, match="NC_002076.2"):
        write_competitive_fasta(str(host), str(virus), str(tmp_path / "o.fa"))


def test_competitive_fasta_same_id_across_host_and_virus_is_fine(tmp_path):
    host = tmp_path / "h.fa"
    host.write_text(">x\nACGT\n")
    virus = tmp_path / "v.fa"
    virus.write_text(">x\nACGT\n")
    out = write_competitive_fasta(str(host), str(virus), str(tmp_path / "o.fa"))
    assert Path(out).read_text().count(">") == 2

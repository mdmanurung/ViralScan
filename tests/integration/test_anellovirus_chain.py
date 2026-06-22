"""Integration test: anellovirus GTF-generation → gene_id labeling chain.

Validates that the anellovirus reference pipeline produces GTFs whose
gene_ids are correctly resolved to genus labels by the detection layer —
without requiring kb, bustools, or any network access.

The test exercises the chain that matters for end-to-end correctness:

  synthetic FASTA
       ↓  _gtf_from_merged_fasta
  anellovirus.gtf (gene_id "{acc}_geneN")
       ↓  re.findall (analysis.py pattern)
  gene_id list
       ↓  virus_name_for_gene(gid, merged_name_map())
  genus labels

A full kb-python pipeline (index → kb count → detection) is an operational
step that requires installed binaries and is not a missing code path.

Run with::

    PYTHONPATH=src python -m pytest -m integration \\
        tests/integration/test_anellovirus_chain.py -v
"""

from __future__ import annotations

import re
import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest

from viralscan.anellovirus import merged_name_map
from viralscan.scripts.build_reference import (
    _gtf_from_merged_fasta,
    build_anellovirus_reference,
)
from viralscan.virus_grouping import group_genes_by_virus, virus_name_for_gene

# Two real accessions from the packaged TSV covering different genera.
_ACC_ALPHA = "NC_002076.2"   # Alphatorquevirus  (viralscan-refseq source)
_ACC_BETA = "AB026929.1"    # Betatorquevirus   (clareaulab source)

_SYNTHETIC_FASTA = textwrap.dedent(f"""\
    >{_ACC_ALPHA} Torque teno virus 1 genome
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    >{_ACC_BETA} Torque teno virus clone SANBAN genome
    TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
""")

_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def _extract_gene_ids(gtf_text: str) -> set[str]:
    """Mirror the pattern used by analysis.py to find gene_ids in a GTF."""
    return set(_GENE_ID_RE.findall(gtf_text))


# ── Tests ──────────────────────────────────────────────────────────────────────


@pytest.mark.integration
class TestAnellovirusLabelingChain:
    """Validate GTF generation → gene_id extraction → genus resolution."""

    def test_gtf_from_merged_fasta_produces_expected_gene_ids(self, tmp_path):
        fasta_path = tmp_path / "anellovirus.fa"
        fasta_path.write_text(_SYNTHETIC_FASTA)
        gtf_path = tmp_path / "anellovirus.gtf"

        _gtf_from_merged_fasta(fasta_path, gtf_path)

        gtf_text = gtf_path.read_text()
        gene_ids = _extract_gene_ids(gtf_text)

        assert f"{_ACC_ALPHA}_gene1" in gene_ids, gene_ids
        assert f"{_ACC_BETA}_gene1" in gene_ids, gene_ids

    def test_gene_ids_resolve_to_correct_genus(self, tmp_path):
        fasta_path = tmp_path / "anellovirus.fa"
        fasta_path.write_text(_SYNTHETIC_FASTA)
        gtf_path = tmp_path / "anellovirus.gtf"
        _gtf_from_merged_fasta(fasta_path, gtf_path)

        gene_ids = _extract_gene_ids(gtf_path.read_text())
        name_map = merged_name_map()

        resolved = {gid: virus_name_for_gene(gid, name_map) for gid in gene_ids}

        assert resolved[f"{_ACC_ALPHA}_gene1"] == "Alphatorquevirus", resolved
        assert resolved[f"{_ACC_BETA}_gene1"] == "Betatorquevirus", resolved

    def test_group_genes_by_virus_separates_genera(self, tmp_path):
        fasta_path = tmp_path / "anellovirus.fa"
        fasta_path.write_text(_SYNTHETIC_FASTA)
        gtf_path = tmp_path / "anellovirus.gtf"
        _gtf_from_merged_fasta(fasta_path, gtf_path)

        gene_ids = list(_extract_gene_ids(gtf_path.read_text()))
        name_map = merged_name_map()

        groups, detected = group_genes_by_virus(gene_ids, name_map)

        assert "Alphatorquevirus" in detected, detected
        assert "Betatorquevirus" in detected, detected
        assert f"{_ACC_ALPHA}_gene1" in groups["Alphatorquevirus"]
        assert f"{_ACC_BETA}_gene1" in groups["Betatorquevirus"]

    def test_build_anellovirus_reference_produces_labelable_gtf(self, tmp_path):
        """Full builder run (mocked NCBI) → GTF gene_ids resolve to genus."""
        fake_fasta = tmp_path / "ncbi" / "merged.fasta"
        fake_fasta.parent.mkdir(parents=True)
        fake_fasta.write_text(_SYNTHETIC_FASTA)
        fake_gtf = tmp_path / "ncbi" / "merged.gtf"
        fake_gtf.write_text("")

        with patch(
            "viralscan.scripts.ncbi_fetch.fetch_reference",
            return_value=(fake_fasta, fake_gtf),
        ):
            result = build_anellovirus_reference(
                out_dir=tmp_path / "out",
                accessions=[_ACC_ALPHA, _ACC_BETA],
                mask=False,
                cluster=False,
                run_kb_ref=False,
            )

        gtf_text = result["gtf"].read_text()
        gene_ids = _extract_gene_ids(gtf_text)
        name_map = merged_name_map()

        assert virus_name_for_gene(f"{_ACC_ALPHA}_gene1", name_map) == "Alphatorquevirus"
        assert virus_name_for_gene(f"{_ACC_BETA}_gene1", name_map) == "Betatorquevirus"
        assert f"{_ACC_ALPHA}_gene1" in gene_ids
        assert f"{_ACC_BETA}_gene1" in gene_ids

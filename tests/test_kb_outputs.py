"""Unit tests for the kb-python output layout (viralscan.kb_outputs).

This is the layout that used to be hand-built as f-strings in multimap.py,
detection.py, umap.py and the Snakefile. These tests pin the convention so a
change to where kb-python writes is a one-place, test-covered edit.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from viralscan.kb_outputs import KbCountOutputs


class TestPaths:
    def test_named_paths_match_kb_python_convention(self) -> None:
        kb = KbCountOutputs.from_config_output("/out/sample/")
        assert kb.root == Path("/out/sample/kb-python")
        assert kb.counts_dir == Path("/out/sample/kb-python/counts_unfiltered")
        assert kb.bus == Path("/out/sample/kb-python/output.bus")
        assert kb.bus_txt == Path("/out/sample/kb-python/output.bus.txt")
        assert kb.ec == Path("/out/sample/kb-python/matrix.ec")
        assert kb.transcripts_txt == Path("/out/sample/kb-python/transcripts.txt")
        assert kb.adata == Path("/out/sample/kb-python/counts_unfiltered/adata.h5ad")
        assert kb.adata_multimap == Path(
            "/out/sample/kb-python/counts_unfiltered/adata_multimap.h5ad"
        )
        assert kb.barcodes == Path(
            "/out/sample/kb-python/counts_unfiltered/cells_x_genes.barcodes.txt"
        )
        assert kb.genes == Path("/out/sample/kb-python/counts_unfiltered/cells_x_genes.genes.txt")
        assert kb.gene_names == Path(
            "/out/sample/kb-python/counts_unfiltered/cells_x_genes.genes.names.txt"
        )

    @pytest.mark.parametrize("output", ["/out/sample/", "/out/sample"])
    def test_trailing_separator_is_normalised(self, output: str) -> None:
        """A trailing os.sep on config['output'] must not change the targets.

        The old f-strings produced a harmless double slash; pathlib normalises it.
        Both forms must resolve to the same file.
        """
        kb = KbCountOutputs.from_config_output(output)
        assert kb.adata == Path("/out/sample/kb-python/counts_unfiltered/adata.h5ad")

    def test_parity_with_legacy_fstring_target(self) -> None:
        """str(kb.adata) must point at the same file the old f-string did."""
        output = "/out/sample/"
        legacy = f"{output}/kb-python/counts_unfiltered/adata.h5ad"
        assert Path(str(KbCountOutputs.from_config_output(output).adata)) == Path(legacy)


class TestCurrentAdata:
    def test_resolves_to_multimap_when_flag_set(self) -> None:
        kb = KbCountOutputs.from_config_output("/out/sample/")
        assert kb.current_adata(multimapping=True) == kb.adata_multimap

    def test_resolves_to_plain_when_flag_unset(self) -> None:
        kb = KbCountOutputs.from_config_output("/out/sample/")
        assert kb.current_adata(multimapping=False) == kb.adata

    def test_resolution_is_by_flag_not_file_existence(self, tmp_path: Path) -> None:
        """Even with no files on disk, the flag alone decides — no I/O."""
        kb = KbCountOutputs.from_config_output(str(tmp_path))
        assert kb.current_adata(multimapping=True).name == "adata_multimap.h5ad"
        assert kb.current_adata(multimapping=False).name == "adata.h5ad"

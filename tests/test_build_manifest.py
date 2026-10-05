"""Index build manifest (PLAN ``DEF-03``): written by both build paths, read by the Run."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from tests.test_virus_identity import ANELLO_GENUS, CATALOGUE, _combined_t2g
from viralscan import reference_strategy, virus_identity
from viralscan.scripts import analysis
from viralscan.scripts.build_reference import _write_index_manifest
from viralscan.virus_identity import (
    HOST,
    MANIFEST_SCHEMA_VERSION,
    BuildManifestContradiction,
    build_identity_table,
    deversion,
    gtf_gene_ids,
    load_build_manifest,
    manifest_path_for_index,
    write_build_manifest,
    write_build_manifest_from_t2g,
)

HOST_IDS = ["ENSG1.1", "ENSG2.1", "MYSTERY_g1"]
VIRAL_IDS = ["NC_007605.1_gene1", "NC_007605.1_gene2", "NC_009334_gene1", "MY_VIRUS_g1"]


def _gtf(path: Path, gene_ids) -> Path:
    path.write_text(
        "".join(
            f'{gid}_seq\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id "{gid}"; transcript_id "{gid}_t";\n'
            for gid in gene_ids
        )
    )
    return path


def _manifest(tmp_path: Path, host=HOST_IDS, viral=VIRAL_IDS) -> Path:
    return write_build_manifest(tmp_path / "index.idx", host, viral, {"builder": "test"})


def _table(tmp_path: Path, gtf_genes, manifest: Path | None = None):
    return build_identity_table(
        _combined_t2g(tmp_path),
        gtf_genes,
        catalogue_rows=CATALOGUE,
        anello_genus=ANELLO_GENUS,
        build_manifest=manifest,
    )


class TestManifestFile:
    def test_path_sits_next_to_the_index(self, tmp_path):
        assert manifest_path_for_index(tmp_path / "x" / "index.idx") == (
            tmp_path / "x" / "index.idx.build_manifest.json"
        )

    def test_round_trip_and_schema(self, tmp_path):
        path = _manifest(tmp_path)
        data = json.loads(path.read_text())
        assert data["schema_version"] == MANIFEST_SCHEMA_VERSION
        assert data["viralscan_version"]
        assert data["provenance"] == {"builder": "test"}
        manifest = load_build_manifest(path)
        assert manifest.viral_gene_ids == {deversion(g) for g in VIRAL_IDS}
        assert manifest.host_gene_ids == {deversion(g) for g in HOST_IDS}

    def test_gene_in_both_sets_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="both host and viral"):
            write_build_manifest(tmp_path / "i.idx", ["G.1"], ["G.2"])

    def test_unknown_schema_major_is_refused(self, tmp_path):
        path = _manifest(tmp_path)
        data = json.loads(path.read_text())
        data["schema_version"] = "9.0"
        path.write_text(json.dumps(data))
        with pytest.raises(ValueError, match="unsupported build manifest schema_version"):
            load_build_manifest(path)

    def test_deversion_strips_only_a_trailing_number(self):
        assert deversion("ENSG00000123.4") == "ENSG00000123"
        assert deversion("NC_007605.1_gene1") == "NC_007605.1_gene1"


class TestBuildersWriteTheManifest:
    def test_build_reference_helper_splits_host_and_viral_gtfs(self, tmp_path):
        host_gtf = _gtf(tmp_path / "host.gtf", ["ENSG1.1", "ENSG2.1"])
        viral_gtf = _gtf(tmp_path / "viral.gtf", ["V_gene1"])
        fasta = tmp_path / "combined.fa"
        fasta.write_text(">x\nACGT\n")
        combined = tmp_path / "combined.gtf"
        combined.write_text(host_gtf.read_text() + viral_gtf.read_text())
        index = tmp_path / "index.idx"
        out = _write_index_manifest(
            index,
            host_gtf=host_gtf,
            viral_gtf=viral_gtf,
            builder="viralscan build-ref",
            extra={"host_species": "homo_sapiens"},
            fasta=fasta,
            gtf=combined,
        )
        assert out == manifest_path_for_index(index)
        manifest = load_build_manifest(out)
        assert manifest.host_gene_ids == {"ENSG1", "ENSG2"}
        assert manifest.viral_gene_ids == {"V_gene1"}
        prov = manifest.provenance
        assert prov["builder"] == "viralscan build-ref"
        assert prov["host_species"] == "homo_sapiens"
        assert len(prov["gtf"]["sha256"]) == 64

    def test_reference_path_splits_by_the_t2g_structural_guard(self, tmp_path):
        from viralscan.menu import _write_reference_manifest

        t2g = _combined_t2g(tmp_path)
        # A combined GTF passed to --reference: host genes must not become viral.
        gtf = _gtf(tmp_path / "all.gtf", ["ENSG1.1", "ENSG2.1", *VIRAL_IDS])
        fasta = tmp_path / "ref.fa"
        fasta.write_text(">x\nACGT\n")
        index = tmp_path / "index.idx"
        _write_reference_manifest(str(index), str(t2g), str(fasta), str(gtf))
        manifest = load_build_manifest(manifest_path_for_index(index))
        assert manifest.viral_gene_ids == {deversion(g) for g in VIRAL_IDS}
        assert manifest.host_gene_ids == {deversion(g) for g in HOST_IDS}

    def test_from_t2g_direct(self, tmp_path):
        path = write_build_manifest_from_t2g(
            tmp_path / "index.idx", _combined_t2g(tmp_path), ["MY_VIRUS_g1"]
        )
        assert "MY_VIRUS_g1" in load_build_manifest(path).viral_gene_ids

    def test_gtf_gene_ids_reads_plain_and_gz(self, tmp_path):
        import gzip

        plain = _gtf(tmp_path / "a.gtf", ["A", "B"])
        gz = tmp_path / "a.gtf.gz"
        gz.write_bytes(gzip.compress(plain.read_bytes()))
        assert gtf_gene_ids(plain) == gtf_gene_ids(gz) == {"A", "B"}


class TestManifestOverridesGtf:
    def test_matching_gtf_gives_the_same_table_as_no_manifest(self, tmp_path):
        manifest = _manifest(tmp_path)
        with_m = _table(tmp_path, ["MY_VIRUS_g1"], manifest)
        without = _table(tmp_path, ["MY_VIRUS_g1"])
        assert with_m == without

    def test_manifest_supplies_the_viral_set_when_gtf_is_empty_of_it(self, tmp_path):
        # MY_VIRUS_g1 is uncatalogued: without the manifest and with a GTF that
        # lacks it, it is host. With a manifest and the same GTF the Run is refused.
        manifest = _manifest(tmp_path)
        assert _table(tmp_path, []).by_gene()["MY_VIRUS_g1"].status == HOST
        with pytest.raises(BuildManifestContradiction, match="not viral under --gtf"):
            _table(tmp_path, [], manifest)

    def test_host_gene_in_gtf_contradicts_the_manifest(self, tmp_path):
        # MYSTERY_g1 is a host gene in the manifest; a --gtf naming it makes it
        # viral (its column 5 is not an index transcript), so the Run is refused.
        manifest = _manifest(tmp_path)
        with pytest.raises(BuildManifestContradiction) as exc:
            _table(tmp_path, ["MY_VIRUS_g1", "MYSTERY_g1"], manifest)
        assert "MYSTERY_g1" in str(exc.value)
        assert "--gtf" in str(exc.value)

    def test_version_drift_is_accepted(self, tmp_path):
        drifted_host = ["ENSG1.9", "ENSG2.3", "MYSTERY_g1"]
        drifted_viral = ["NC_007605.1_gene1", "NC_007605.1_gene2", "NC_009334_gene1", "MY_VIRUS_g1"]
        manifest = _manifest(tmp_path, host=drifted_host, viral=drifted_viral)
        table = _table(tmp_path, ["MY_VIRUS_g1.7"], manifest)
        assert table.by_gene()["MY_VIRUS_g1"].viral

    def test_viral_gene_version_drift_between_gtf_and_manifest_is_accepted(self, tmp_path):
        # Gene IDs are compared de-versioned, so a manifest written with a versioned
        # viral ID still matches the unversioned t2g/GTF ID.
        manifest = _manifest(
            tmp_path,
            viral=["NC_007605.1_gene1", "NC_007605.1_gene2", "NC_009334_gene1", "MY_VIRUS_g1.2"],
        )
        assert _table(tmp_path, ["MY_VIRUS_g1"], manifest).by_gene()["MY_VIRUS_g1"].viral

    def test_manifest_not_describing_the_index_is_refused(self, tmp_path):
        manifest = write_build_manifest(tmp_path / "index.idx", ["ENSG1.1"], ["MY_VIRUS_g1"])
        with pytest.raises(BuildManifestContradiction, match="neither host nor viral"):
            _table(tmp_path, ["MY_VIRUS_g1"], manifest)

    def test_manifest_genes_absent_from_the_index_are_ignored(self, tmp_path):
        manifest = _manifest(tmp_path, viral=[*VIRAL_IDS, "NOT_INDEXED_g1"])
        assert _table(tmp_path, ["MY_VIRUS_g1"], manifest).by_gene()["MY_VIRUS_g1"].viral

    def test_none_is_unchanged(self, tmp_path):
        # The pre-manifest behaviour: the --gtf set alone decides rule 2.
        genes = _table(tmp_path, ["MY_VIRUS_g1", "MYSTERY_g1"], None).by_gene()
        assert genes["MYSTERY_g1"].viral
        assert genes["ENSG1.1"].status == HOST


class TestAnalysisDiscovery:
    def _config(self, tmp_path, t2g: Path):
        return SimpleNamespace(
            transcripts=str(t2g), index=str(tmp_path / "index.idx"), output=str(tmp_path)
        )

    @pytest.fixture(autouse=True)
    def _packaged_catalogue(self, monkeypatch):
        # Keep the test independent of the packaged catalogue.
        real = build_identity_table

        def patched(t2g, gtf, **kw):
            kw.setdefault("catalogue_rows", CATALOGUE)
            kw.setdefault("anello_genus", ANELLO_GENUS)
            return real(t2g, gtf, **kw)

        monkeypatch.setattr(virus_identity, "build_identity_table", patched)

    def test_manifest_next_to_the_index_is_used(self, tmp_path):
        t2g = _combined_t2g(tmp_path)
        _manifest(tmp_path)
        with pytest.raises(BuildManifestContradiction):
            analysis.write_identity_table(self._config(tmp_path, t2g), {"MYSTERY_g1"})

    def test_absent_manifest_falls_back_to_gtf_and_warns(self, tmp_path, caplog):
        t2g = _combined_t2g(tmp_path)
        with caplog.at_level(logging.WARNING):
            out = analysis.write_identity_table(self._config(tmp_path, t2g), {"MY_VIRUS_g1"})
        assert out is not None and out.exists()
        assert "No index build manifest" in caplog.text


class TestReferenceStrategyPassesViralGtf:
    def test_combined_strategy_passes_the_viral_only_gtf(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            reference_strategy,
            "fastq_paths_for_dataset",
            lambda manifest, row: (tmp_path / "r1.fq", tmp_path / "r2.fq"),
        )
        manifest = {
            "references": {
                "viralscan": {
                    "all_virus": {
                        "kallisto_index": "v.idx",
                        "t2g": "v.t2g",
                        "gtf": "viral_only.gtf",
                    },
                    "combined": {
                        "kallisto_index": "c.idx",
                        "t2g": "c.t2g",
                        "gtf": "combined_with_host.gtf",
                    },
                }
            }
        }
        row = {"row_id": "r", "method": "viralscan", "technology": "10xv3"}
        for strategy, index in (("combined", "c.idx"), ("all_virus", "v.idx")):
            (cmd,) = reference_strategy.commands_for_row(
                {**row, "reference_strategy": strategy}, manifest, tmp_path
            )
            assert cmd[cmd.index("-i") + 1] == index
            assert cmd[cmd.index("-gtf") + 1] == "viral_only.gtf"

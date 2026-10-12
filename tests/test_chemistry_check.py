"""Library diagnostics (CHEM-01): bulk vs single-cell, 3' vs 5', post-kb sanity gate."""

import random

from viralscan import chemistry_check as cc

random.seed(7)


def _rand(n):
    return "".join(random.choice("ACGT") for _ in range(n))


_SHARED_TAIL = _rand(2)  # one tail per import: every call shares it, as the old default did


def _barcoded(n_cells=50, reads_per_cell=40, tail=_SHARED_TAIL):
    cells = [_rand(16) for _ in range(n_cells)]
    return [random.choice(cells) + _rand(10) + tail for _ in range(n_cells * reads_per_cell)]


class TestLibraryKind:
    def test_barcoded_reads_are_single_cell(self):
        assert cc.library_kind(_barcoded())[0] == "single-cell"

    def test_on_list_match_alone_is_single_cell(self):
        assert (
            cc.library_kind([_rand(100) for _ in range(500)], on_list_rate=0.9)[0] == "single-cell"
        )

    def test_genomic_reads_are_bulk(self):
        assert cc.library_kind([_rand(100) for _ in range(2000)])[0] == "bulk"

    def test_cdna_mate_as_r1_is_bulk(self):
        transcripts = [_rand(98) for _ in range(30)]
        reads = [random.choice(transcripts) for _ in range(3000)]
        kind, reason = cc.library_kind(reads)
        assert kind == "bulk" and "swapped" in reason

    def test_dropseq_barcode_is_found_at_12bp(self):
        cells = [_rand(12) for _ in range(40)]
        reads = [random.choice(cells) + _rand(8) for _ in range(2000)]
        assert cc.library_kind(reads)[0] == "single-cell"

    def test_short_unique_reads_are_unclear(self):
        assert cc.library_kind([_rand(20) for _ in range(2000)])[0] == "unclear"


class TestInferEnd:
    F5 = {"forward": 12.0, "reverse": 50.0, "unstranded": 58.0}  # F-020 5' shape
    F3 = {"forward": 70.0, "reverse": 8.0, "unstranded": 72.0}

    def test_r1_adapter_is_decisive(self):
        assert cc.infer_end("5p", None)[0] == "5p"

    def test_trimmed_r1_uses_the_pilot(self):
        assert cc.infer_end(None, self.F5) == ("5p", "strand pilot (R1 is trimmed)")
        assert cc.infer_end(None, self.F3)[0] == "3p"

    def test_conflict_is_reported_not_resolved(self):
        end, basis = cc.infer_end("3p", self.F5)
        assert end is None and basis.startswith("conflict")

    def test_nothing_to_go_on(self):
        assert cc.infer_end(None, None)[0] is None


class TestSanityGate:
    def test_low_rate_is_an_error_that_points_at_check_chemistry(self):  # F-020: 8.2 %
        (f,) = cc.sanity_gate({"p_pseudoaligned": 8.2})
        assert f["level"] == "error" and "check-chemistry" in f["message"]

    def test_no_strand_is_not_blamed_on_a_forward_default(self):
        # kallisto 0.52.0 (pinned): no flag is unstranded, so the "forward-strand default"
        # explanation of F-020 no longer applies (M-01).
        for strand in (None, "unstranded", "reverse"):
            (f,) = cc.sanity_gate({"p_pseudoaligned": 8.2}, strand=strand)
            assert "forward-strand default" not in f["message"]
            assert "unstranded under kallisto 0.52.0" in f["message"]

    def test_explicit_forward_strand_is_named_as_the_cause(self):
        (f,) = cc.sanity_gate({"p_pseudoaligned": 8.2}, strand="forward")
        assert "--strand forward was used" in f["message"]
        assert "unstranded under kallisto" not in f["message"]

    def test_sanity_report_passes_the_strand_through(self, tmp_path):
        (tmp_path / "run_info.json").write_text('{"p_pseudoaligned": 8.2}')
        report = cc.sanity_report(tmp_path, strand="forward")
        assert "--strand forward was used" in report["findings"][0]["message"]

    def test_moderate_is_a_warning(self):
        assert cc.sanity_gate({"p_pseudoaligned": 25.0})[0]["level"] == "warning"

    def test_healthy_run_is_silent(self):
        assert cc.sanity_gate({"p_pseudoaligned": 73.4}) == []

    def test_host_filtered_index_is_not_judged(self):
        assert cc.sanity_gate({"p_pseudoaligned": 3.0}, host_in_index=False) == []

    def test_missing_run_info_is_silent(self, tmp_path):
        assert cc.sanity_gate_from_dir(tmp_path) == []

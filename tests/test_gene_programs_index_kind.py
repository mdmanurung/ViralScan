"""PROG-INDEXKIND-01: gene_program_summary.tsv records the index kind and host filter."""

import inspect

from viralscan.scripts import gene_programs


def test_index_state_records_index_kind_and_host_filter() -> None:
    assert gene_programs.index_state({"index_kind": "combined"}, None) == ("combined", "none")
    assert gene_programs.index_state({"index_kind": "virus_only"}, "starsolo") == (
        "virus_only",
        "starsolo",
    )
    # an H5AD without the multimap stamp is labelled, not guessed
    assert gene_programs.index_state({}, "") == ("unknown", "none")


def test_main_appends_the_columns_after_existing_ones() -> None:
    source = inspect.getsource(gene_programs.main)
    assert (
        source.index('summary["n_called_cells"]')
        < source.index('summary["index_kind"], summary["host_filter"]')
        < source.index("write_program_outputs")
    )

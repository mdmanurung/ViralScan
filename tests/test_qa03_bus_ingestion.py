"""QA-03: BUS ingestion is chunk-size and order invariant (property tests on tiny BUS text).

Ingestion is already one pass over ``bustools text`` rows (SW-03: ``_iter_bus_text_records`` feeding
``iter_cb_umi_molecules``), so there is no chunk boundary that can split a CB-UMI molecule; a separate
chunked mode is not added. These tests prove the properties a chunked reader would have to keep, over
seeded random fixtures (``hypothesis`` is not a dependency): every read-buffer size, including ones that
cut a molecule's records mid-row, gives the same layers and audit, and unsorted text is refused rather
than silently depending on row order.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.test_multimapping_properties import build, random_fixture
from viralscan.defaults import MULTIMAP_METHODS

BUFFER_SIZES = (1, 2, 3, 5, 7, 64, 8192)


def write_bus_text(fixture, path, *, ordered: bool = True, seed: int = 0):
    frame = fixture[0]
    frame = (
        frame.sort_values(["barcode", "umi", "ec"])
        if ordered
        else frame.sample(frac=1, random_state=seed)
    )
    frame.to_csv(path, sep="\t", header=False, index=False)
    return path


@pytest.mark.parametrize("method", MULTIMAP_METHODS)
@pytest.mark.parametrize("seed", range(5))
def test_every_buffer_size_gives_the_same_layers_and_audit(method, seed, tmp_path):
    fixture = random_fixture(seed)
    path = write_bus_text(fixture, tmp_path / "bus.txt")
    reference = build(fixture, method=method, frame=path)
    in_memory = build(fixture, method=method)
    assert reference.audit == in_memory.audit
    for size in BUFFER_SIZES:
        streamed = build(fixture, method=method, frame=path, bus_buffer_size=size)
        assert streamed.audit == reference.audit, size
        for layer in ("unique", "corrected", "equal"):
            np.testing.assert_array_equal(
                getattr(streamed, layer).toarray(), getattr(reference, layer).toarray()
            )


@pytest.mark.parametrize("seed", range(8))
def test_unsorted_text_is_refused_not_silently_order_dependent(seed, tmp_path):
    fixture = random_fixture(seed)
    path = write_bus_text(fixture, tmp_path / "bus.txt", ordered=False, seed=seed)
    # out-of-order on-list rows must raise; a shuffle that leaves them sorted is accepted and must
    # then equal the sorted result
    try:
        shuffled = build(fixture, method="host-conservative", frame=path)
    except ValueError as exc:
        assert "sorted" in str(exc)
    else:
        sorted_path = write_bus_text(fixture, tmp_path / "sorted2.txt")
        expected = build(fixture, method="host-conservative", frame=sorted_path)
        np.testing.assert_array_equal(shuffled.unique.toarray(), expected.unique.toarray())
        assert shuffled.audit == expected.audit


@pytest.mark.parametrize("seed", range(8))
def test_row_order_of_an_in_memory_frame_does_not_matter(seed):
    fixture = random_fixture(seed)
    reference = build(fixture, method="host-conservative")
    for shuffle in range(3):
        frame = fixture[0].sample(frac=1, random_state=100 * seed + shuffle)
        other = build(fixture, method="host-conservative", frame=frame)
        assert other.audit == reference.audit
        np.testing.assert_array_equal(other.corrected.toarray(), reference.corrected.toarray())

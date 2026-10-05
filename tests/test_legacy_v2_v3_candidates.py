from pathlib import Path

import pandas as pd

from scripts.summarize_legacy_v2_v3_candidates import summarize_candidates


def test_candidate_summary_classifies_transitions_and_lane_repeat(
    tmp_path: Path,
) -> None:
    runs = pd.DataFrame(
        {
            "run_id": ["r1", "r2"],
            "logical_id": ["sample_S1_L001", "sample_S1_L002"],
            "sample_class": ["skin", "skin"],
        }
    )
    legacy = pd.DataFrame(
        {
            "run_id": ["r1", "r2"],
            "record_type": ["virus", "virus"],
            "identifier": ["virus-a", "virus-b"],
            "reconstructed_value": [2.0, 3.0],
        }
    )
    virus = pd.DataFrame(
        {
            "run_id": ["r1", "r1", "r2", "r2"],
            "virus_name": ["virus-a", "virus-c", "virus-a", "virus-b"],
            "v3_unique": [0, 1, 0, 1],
            "v3_equal": [1, 1, 1, 1],
            "v3_host_conservative": [1, 1, 1, 0],
        }
    )
    paths = {
        "legacy": tmp_path / "legacy.tsv",
        "virus": tmp_path / "virus.tsv",
        "runs": tmp_path / "runs.tsv",
    }
    for name, frame in [
        ("legacy", legacy),
        ("virus", virus),
        ("runs", runs),
    ]:
        frame.to_csv(paths[name], sep="\t", index=False)

    output = summarize_candidates(
        paths["legacy"], paths["virus"], paths["runs"], tmp_path / "out.tsv"
    )
    result = pd.read_csv(output, sep="\t")
    first = result[(result["run_id"] == "r1") & (result["virus_name"] == "virus-a")]
    assert first.iloc[0]["v3_unique_transition"] == "lost"
    assert first.iloc[0]["v3_host_conservative_transition"] == "persisted"
    assert first.iloc[0]["host_repeated_across_lanes"]
    new = result[(result["run_id"] == "r1") & (result["virus_name"] == "virus-c")]
    assert new.iloc[0]["v3_host_conservative_transition"] == "new-nonzero"
    assert new.iloc[0]["v3_unique_transition"] == "new-nonzero"
    absent = result[(result["run_id"] == "r2") & (result["virus_name"] == "virus-a")]
    assert absent.iloc[0]["v3_unique_transition"] == "absent"

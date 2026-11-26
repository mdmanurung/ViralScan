"""CONFIG-LOAD-01: one loader (`RunConfig.load`) and a validated `with_overrides`."""

import pytest

from viralscan.runconfig import RunConfig


@pytest.fixture
def config_yaml(tmp_path):
    # RunConfig() holds the defaults; from_snakemake_config(to_dict()) is the validated twin.
    cfg = RunConfig.from_snakemake_config(RunConfig(output=str(tmp_path / "out") + "/").to_dict())
    path = tmp_path / "config.yaml"
    cfg.to_yaml(path)
    return cfg, path


def test_load_reads_what_to_yaml_wrote(config_yaml):
    cfg, path = config_yaml
    assert RunConfig.load(path) == cfg
    assert RunConfig.from_yaml(path) == cfg


def test_with_overrides_round_trips_without_changes(config_yaml):
    cfg, _ = config_yaml
    assert cfg.with_overrides() == cfg


def test_with_overrides_sets_values_and_leaves_the_source_alone(config_yaml):
    cfg, _ = config_yaml
    new = cfg.with_overrides(multimap_method="em-global", cores=3)
    assert (new.multimap_method, new.cores) == ("em-global", 3)
    assert (cfg.multimap_method, cfg.cores) != ("em-global", 3)


def test_with_overrides_coerces_legacy_string_booleans(config_yaml):
    cfg, _ = config_yaml
    assert cfg.with_overrides(gene_programs="True").gene_programs is True
    assert cfg.with_overrides(gene_programs="False").gene_programs is False


def test_with_overrides_runs_the_constructor_checks(config_yaml):
    cfg, _ = config_yaml
    with pytest.raises(ValueError, match="cores must be >= 1"):
        cfg.with_overrides(cores=0)
    with pytest.raises(ValueError, match="multimap_pseudocount must be > 0"):
        cfg.with_overrides(multimap_pseudocount=0)


def test_with_overrides_rejects_an_unknown_field(config_yaml):
    cfg, _ = config_yaml
    with pytest.raises(TypeError):
        cfg.with_overrides(not_a_field=1)


def test_load_keeps_refusing_virus_key(config_yaml, tmp_path):
    _, path = config_yaml
    path.write_text(path.read_text() + "positive_control_scope: virus_key\n")
    with pytest.raises(ValueError, match="virus_key|VAL-RA-CAL"):
        RunConfig.load(path)

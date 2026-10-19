"""total_budget_v1: --cores bounds every native worker request (legacy request when unset)."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from viralscan import menu
from viralscan.resources import bounded_workers
from viralscan.scripts import build_reference as br

LEGACY_CORES = [None, 1, 2, 3, 6, 8, 32]


@pytest.mark.parametrize("legacy", [8, 3])
@pytest.mark.parametrize("cores", [1, 2, 3, 6, 8, 32])
def test_bounded_workers_table(legacy, cores):
    assert bounded_workers(legacy, cores) == min(legacy, cores)
    assert bounded_workers(legacy, None) == legacy


@pytest.mark.parametrize("cores", LEGACY_CORES)
@pytest.mark.parametrize("entrypoint", ["menu", "finish", "finish_dlist"])
def test_reference_entrypoints_bound_workers(tmp_path, monkeypatch, entrypoint, cores):
    fasta, gtf = tmp_path / "reference.fa", tmp_path / "reference.gtf"
    fasta.write_text(">V1\nACGT\n")
    gtf.write_text('V1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "V1";\n')
    commands, receipts = [], []

    def run(command, out_dir):
        commands.append(command)
        return {"elapsed_seconds": 1.0}

    monkeypatch.setattr(br, "_run_kb_ref", run)
    monkeypatch.setattr(br, "record_reference_build", lambda *a, **kw: receipts.append(kw))
    monkeypatch.setattr(
        br, "prepare_reference_inputs", lambda *a, **kw: (fasta, tmp_path / "manifest.json")
    )
    monkeypatch.setattr(br, "_write_index_manifest", lambda *a, **kw: None)
    monkeypatch.setattr(br, "record_index_inspect", lambda *a, **kw: None)
    monkeypatch.setattr(br.shutil, "which", lambda tool: "/native/kb")
    monkeypatch.setattr(menu, "_prepare_kb_ref_inputs", lambda *a: (str(fasta), str(gtf)))
    monkeypatch.setattr(menu, "_write_reference_manifest", lambda *a: None)
    if entrypoint == "menu":
        menu._build_kb_ref(tmp_path / "out", str(fasta), str(gtf), cores)
    else:
        br._finish_reference_build(
            tmp_path,
            fasta,
            gtf,
            tmp_path / "manifest.json",
            run_kb_ref=True,
            genome_dlist=tmp_path / "host.fa" if entrypoint == "finish_dlist" else None,
            host_gtf=None,
            viral_gtf=gtf,
            builder="test",
            extra={},
            cores=cores,
        )
    command = commands[0]
    assert command[1:4] == ["ref", "-t", str(8 if cores is None else min(8, cores))]
    assert command.count("-t") == 1
    assert command[-2:] == [str(fasta), str(gtf)]
    assert ("--d-list" in command) == (entrypoint == "finish_dlist")
    assert receipts[0]["command"] == command


@pytest.mark.parametrize("cores", LEGACY_CORES)
@pytest.mark.parametrize("fails", [False, True])
def test_homology_bounds_workers_and_preserves_search_and_errors(
    tmp_path, monkeypatch, fails, cores
):
    fasta = tmp_path / "viral.fa"
    fasta.write_text(">V1\n" + "ACGT" * 25 + "\n")
    host, output = tmp_path / "host.fa", tmp_path / "out.tsv"
    commands = []
    failure = subprocess.CalledProcessError(17, ["native"], stderr="complete stderr")

    def run(command, **kwargs):
        commands.append(command)
        assert kwargs == {"check": True, "capture_output": True, "text": True}
        if fails:
            raise failure
        return subprocess.CompletedProcess(command, 0, stdout="")

    monkeypatch.setattr(br.shutil, "which", lambda tool: "/native/minimap2")
    monkeypatch.setattr(br.subprocess, "run", run)
    if fails:
        with pytest.raises(subprocess.CalledProcessError) as caught:
            br.measure_host_homology(fasta, host, output, cores)
        assert caught.value is failure
        assert not output.exists()
    else:
        br.measure_host_homology(fasta, host, output, cores)
        assert output.is_file()
    assert commands == [
        [
            "/native/minimap2",
            "-t",
            str(3 if cores is None else min(3, cores)),
            "-c",
            "-k",
            "15",
            "-w",
            "10",
            "-N",
            "200",
            "-p",
            "0.01",
            "-s",
            "40",
            "-m",
            "20",
            str(fasta),
            str(host),
        ]
    ]


@pytest.mark.parametrize("cores", LEGACY_CORES)
def test_clustering_bounds_cpu_and_keeps_memory_and_identity(tmp_path, monkeypatch, cores):
    fasta, output = tmp_path / "in.fa", tmp_path / "out.fa"
    commands = []

    def run(command, **kwargs):
        commands.append(command)
        output.write_text(">V1\nACGT\n")

    monkeypatch.setattr(br.shutil, "which", lambda tool: "/native/cd-hit-est")
    monkeypatch.setattr(br.subprocess, "run", run)
    assert br._run_cdhit_est(fasta, output, cores=cores)
    assert commands == [
        [
            "/native/cd-hit-est",
            "-i",
            str(fasta),
            "-o",
            str(output),
            "-c",
            "0.95",
            "-n",
            "8",
            "-M",
            "8000",
            "-T",
            "0" if cores is None else str(cores),
            "-d",
            "0",
        ]
    ]


@pytest.mark.integration
@pytest.mark.parametrize("cores", [1, 6, 8, 32])
@pytest.mark.parametrize("whitelist", [False, True])
@pytest.mark.parametrize("override", [False, True])
def test_count_legacy_request_and_library_environment(tmp_path, cores, whitelist, override):
    """Execute the actual rule with an argv/env recorder; native bytes have a separate gate."""
    if not shutil.which("snakemake"):
        pytest.skip("snakemake not installed")
    out, binary = tmp_path / "out", tmp_path / "bin"
    (out / "log").mkdir(parents=True)
    binary.mkdir()
    recorder = binary / "kb"
    recorder.write_text(
        f"#!{sys.executable}\n"
        "import json,os,pathlib,sys\n"
        "out=pathlib.Path(sys.argv[sys.argv.index('-o')+1])\n"
        "(out/'counts_unfiltered').mkdir(parents=True)\n"
        "(out/'counts_unfiltered/adata.h5ad').write_bytes(b'fixture')\n"
        "(out/'request.json').write_text(json.dumps(dict(argv=sys.argv[1:],env={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']})))\n"
    )
    recorder.chmod(0o755)
    onlist = tmp_path / "onlist.txt"
    onlist.write_text("ACGT\n")
    config = dict(
        output=str(out) + "/",
        index="/fake/index.idx",
        transcripts="/fake/t2g.txt",
        technology="10xv3",
        whitelist=str(onlist) if whitelist else None,
        host_index=None,
        host_filter_aligner=None,
        cores=cores,
        kb_r1="/fake/R1.fastq",
        kb_r2="/fake/R2.fastq",
    )
    configfile = out / "config.yaml"
    configfile.write_text(yaml.safe_dump(config))
    command = [
        "snakemake",
        "--snakefile",
        str(Path(__file__).parents[1] / "src/viralscan/Snakefile"),
        str(out / "log/kb.done"),
        "--cores",
        str(cores),
        "--configfile",
        str(configfile),
    ]
    if override:
        command += ["--set-threads", "kb_count=2"]
    expected = min(8, cores, 2) if override else min(8, cores)
    env = {
        **os.environ,
        "PATH": str(binary) + os.pathsep + os.environ["PATH"],
        "XDG_CACHE_HOME": str(tmp_path / "cache"),
        "PYTHONPATH": str(Path(__file__).parents[1] / "src"),
    }
    done = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr
    request = json.loads((out / "kb-python/request.json").read_text())
    assert request["argv"][:3] == ["count", "-t", str(expected)]
    assert request["argv"].count("-t") == 1
    assert ("-w" in request["argv"]) is whitelist
    assert set(request["env"].values()) == {"1"}  # numerical libraries stay pinned at one worker

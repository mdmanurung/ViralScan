"""total_budget_v1: --cores bounds every native worker request (legacy request when unset)."""

import argparse
import gzip
import json
import os
import random
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


def test_main_passes_cli_cores_to_reference_build(tmp_path, monkeypatch):
    fasta, gtf = tmp_path / "r.fa", tmp_path / "r.gtf"
    fasta.write_text(">V1\nACGTACGTAC\n")
    gtf.write_text('V1\tt\texon\t1\t4\t.\t+\t.\tgene_id "V1";\n')
    for name in ("R1", "R2"):
        (tmp_path / f"{name}.fastq.gz").touch()
    seen = []

    def stop(*args):
        seen.append(args)
        raise SystemExit(0)

    monkeypatch.setattr(menu, "_build_kb_ref", stop)
    monkeypatch.setattr(menu, "_check_required_tools", lambda: None)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        ["viralscan", "--reference", "--fasta", str(fasta), "--gtf", str(gtf), "-o", "out"]
        + ["-s1", "R1.fastq.gz", "-s2", "R2.fastq.gz", "-x", "10xv3", "--force-technology"]
        + ["--yes", "--cores", "2", "--cell-calling", "knee"],
    )
    with pytest.raises(SystemExit) as done:
        menu.main()
    assert done.value.code == 0
    assert seen[0][3] == 2


def _recording_native_run(commands):
    """Stand-in for every native tool: record argv, create the file each would write."""

    def run(command, **kwargs):
        commands.append(list(command))
        tool = Path(command[0]).name
        if tool == "cd-hit-est":
            Path(command[command.index("-o") + 1]).write_text(
                Path(command[command.index("-i") + 1]).read_text()
            )
        elif tool == "kb":
            Path(command[command.index("-i") + 1]).write_text("index")
            Path(command[command.index("-g") + 1]).write_text("t2g")
        return subprocess.CompletedProcess(command, 0, stdout="")

    return run


@pytest.mark.parametrize("panel", ["anellovirus", "combined"])
def test_build_ref_main_threads_cores_to_every_native_tool(tmp_path, monkeypatch, panel):
    """viralscan build-ref --cores 2: kb, minimap2 and cd-hit-est all get 2 (not 8/3/all)."""
    from viralscan import data_fetch
    from viralscan.scripts import ncbi_fetch

    def dna(n, seed):
        rng = random.Random(seed)
        return "".join(rng.choice("ACGT") for _ in range(n))

    host = tmp_path / "host.fa"
    host.write_text(">chr1\n" + dna(200, 1) + "\n")
    commands: list[list[str]] = []
    monkeypatch.setattr(br.subprocess, "run", _recording_native_run(commands))
    monkeypatch.setattr(br.shutil, "which", lambda tool: f"/native/{tool}")
    monkeypatch.setattr(br, "record_reference_build", lambda *a, **kw: None)
    monkeypatch.setattr(br, "_write_index_manifest", lambda *a, **kw: None)
    monkeypatch.setattr(br, "record_index_inspect", lambda *a, **kw: None)
    args = argparse.Namespace(
        output=tmp_path / "out", genome_dlist=str(host), no_mask=True, cores=2, cluster=True
    )
    viral = tmp_path / "viral.fa"
    viral.write_text(f">NC_000001.1 virus\n{dna(400, 2)}\n")
    if panel == "anellovirus":
        monkeypatch.setattr(data_fetch, "bundled_anellovirus_fasta", lambda cache: viral)
        args.reference_panel = "anellovirus"
    else:
        cdna, host_gtf = tmp_path / "h.cdna.fa.gz", tmp_path / "h.gtf.gz"
        with gzip.open(cdna, "wt") as fh:
            fh.write(">ENST1.1 cdna gene:ENSG1.1\n" + dna(80, 3) + "\n")
        with gzip.open(host_gtf, "wt") as fh:
            fh.write('chr1\tEnsembl\texon\t1\t8\t.\t+\t.\tgene_id "HOST1";\n')
        viral_gtf = tmp_path / "viral.gtf"
        viral_gtf.write_text('NC_000001.1\tNCBI\texon\t1\t8\t.\t+\t0\tgene_id "V";\n')
        monkeypatch.setattr(br, "fetch_host_cdna", lambda *a, **kw: (cdna, host_gtf))
        monkeypatch.setattr(ncbi_fetch, "fetch_reference", lambda *a, **kw: (viral, viral_gtf))
        args.reference_panel, args.host = None, "human"
        args.virus_accessions = ["NC_000001.1"]
    br.build_ref_main(args)
    by_tool = {Path(c[0]).name: c for c in commands}
    assert by_tool["kb"][by_tool["kb"].index("-t") + 1] == "2"
    assert by_tool["minimap2"][by_tool["minimap2"].index("-t") + 1] == "2"
    if panel == "anellovirus":
        assert by_tool["cd-hit-est"][by_tool["cd-hit-est"].index("-T") + 1] == "2"


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

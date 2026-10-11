"""evidence_run.run_evidence: invalid flag combos are rejected only AFTER the replay,
and are silently ignored when the replay finds 0 reads. Stubs copy
tests/test_evidence_subcommand.py::_run_stubbed_evidence."""
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from viralscan.scripts import evidence_run


def run(viral_reads, **extra):
    tmp = Path(tempfile.mkdtemp())
    run_dir = tmp / "run"
    run_dir.mkdir()
    (run_dir / "config.yaml").write_text("technology: 10xv3\n")
    kb_dir = run_dir / "kb-python"
    kb_dir.mkdir()
    for n in ("resolved.bus.txt", "genes.txt", "transcripts.txt", "kb.ec"):
        (kb_dir / n).write_text("")
    r1 = tmp / "R1.fq.gz"
    r2 = tmp / "R2.fq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    idx = tmp / "i.idx"
    idx.write_bytes(b"")
    config = SimpleNamespace(sample1=str(r1), sample2=str(r2), kb_r1=str(r1), kb_r2=str(r2),
                             technology="10xv3", transcripts=str(kb_dir / "transcripts.txt"),
                             index=str(idx), whitelist=None, strand=None, multimap_method="equal")
    kb = SimpleNamespace(resolved_bus_txt=kb_dir / "resolved.bus.txt", genes=kb_dir / "genes.txt",
                         transcripts_txt=kb_dir / "transcripts.txt", ec=kb_dir / "kb.ec", root=kb_dir)
    flagged = tmp / "f.txt"
    flagged.write_text("")
    identity = SimpleNamespace(viral_gene_ids=lambda: {"g1"})
    stats = SimpleNamespace(viral_reads=viral_reads, total_reads=10)
    args = SimpleNamespace(run_dir=str(run_dir), output=str(tmp / "ev"), virus="EBV", cores=1,
                           verbose=False, quiet=True, viral_fasta=None, **extra)
    with (
        patch.object(evidence_run.RunConfig, "from_yaml", return_value=config),
        patch.object(evidence_run.KbCountOutputs, "from_config_output", return_value=kb),
        patch.object(evidence_run, "load_transcripts", return_value=(["tx1"], {"tx1": "g1"})),
        patch.object(evidence_run, "_replay_ec_map", return_value={}),
        patch.object(evidence_run, "load_run_identity", return_value=identity),
        patch.object(evidence_run, "resolve_viral_target", return_value=("EBV", {"g1"})),
        patch.object(evidence_run, "replay_exact_target_bus", return_value=flagged) as replay,
        patch.object(evidence_run, "parse_flagged_target_bus", return_value={}),
        patch.object(evidence_run, "extract_exact_reads_by_number", return_value=stats),
    ):
        try:
            evidence_run.run_evidence(args)
            outcome = "returned normally (flag silently ignored)"
        except SystemExit as exc:
            outcome = f"SystemExit({exc.code})"
    return outcome, replay.called


print("--blast without --viral-fasta, replay finds 0 reads:", run(0, blast=True))
print("--blast without --viral-fasta, replay finds 5 reads:", run(5, blast=True))

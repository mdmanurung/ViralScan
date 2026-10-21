"""Real-library acceptance of the star_confirm rule on EBV GSM4796271 (PLAN STAR-ACCEPT-01).

Runs the REAL Snakefile rule and script on the Track A host-filtered reads and compares the pruned BAM with
the earlier manual STAR run (same index, no whitelist, unstranded). Whitelist and strand change barcode tags
and Solo counts only, so per-contig alignment records must be identical; anything else is a regression.
Inputs under viralscan_work are read-only: the index is adopted through a symlink view plus a manifest,
and the run directory is new.
"""

import subprocess
import sys
from pathlib import Path

import yaml

from viralscan import star_hostvirus as sh

M = Path("/exports/archive/hg-funcgenom-research/mdmanurung")
W = M / "viralscan_work"
SRC_RUN = W / "trackA_round1/runs/ebv/twostep/GSM4796271/GSM4796271"
OLD_IDX = W / "star_hostvirus/index"
OLD_BAM = W / "star_hostvirus/GSM4796271/Aligned.sortedByCoord.out.bam"
HG = Path(
    "/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A"
)
ACC = Path(sys.argv[1])
SNAKEFILE = Path(sh.__file__).parent / "Snakefile"


def idxstats(bam: Path) -> dict[str, int]:
    out = subprocess.run(
        ["samtools", "idxstats", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    return {
        f[0]: int(f[2])
        for f in (ln.split("\t") for ln in out.splitlines())
        if f[0].startswith("VIRUS|")
    }


def main() -> int:
    idx = ACC / "index"
    idx.mkdir(parents=True, exist_ok=True)
    for f in OLD_IDX.iterdir():
        if not (idx / f.name).exists():
            (idx / f.name).symlink_to(f)
    fa = W / "star_hostvirus/ref/hostvirus.fa"
    inputs = {
        "host_fasta": HG / "fasta/genome.fa",
        "host_gtf": HG / "genes/genes.gtf",
        "viral_fasta": M / "viral_ref_cat42d/build/viral.fa",
        "viral_gtf": M / "viral_ref_cat42d/build/viral_panel.gtf",
    }
    if not (idx / sh.INDEX_MANIFEST).exists():
        sh.write_index_manifest(
            idx,
            inputs,
            sh.count_contigs(fa),
            sh.genome_generate_cmd(
                "STAR", idx, fa, W / "star_hostvirus/ref/hostvirus.gtf", 8, 85_000_000_000
            ),
            sh.star_version(),
            8,
            85_000_000_000,
        )
    problems = sh.validate_index_manifest(idx)
    assert not problems, problems

    run = ACC / "run"
    (run / "results").mkdir(parents=True, exist_ok=True)
    (run / "log").mkdir(exist_ok=True)
    (run / "results/virus_identity.tsv").write_bytes(
        (SRC_RUN / "results/virus_identity.tsv").read_bytes()
    )
    cfg = yaml.safe_load((SRC_RUN / "config.yaml").read_text())
    cfg.update(output=str(run) + "/", star_confirm=True, star_hostvirus_index=str(idx), cores=8)
    (run / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=True))
    done = f"{run}/log/star_confirm.done"
    subprocess.run(
        [
            "snakemake",
            done,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(run / "config.yaml"),
            "--cores",
            "8",
            "--nolock",
            # only the new rule: a fresh run dir would otherwise re-run host_filter/kb_count/analysis (hours)
            "--allowed-rules",
            "star_confirm",
            "--rerun-triggers",
            "mtime",
        ],
        check=True,
        cwd=run,
    )

    old, new = idxstats(OLD_BAM), idxstats(run / "star_hostvirus/viral.bam")
    bad = {
        c: (old.get(c, 0), new.get(c, 0))
        for c in set(old) | set(new)
        if old.get(c, 0) != new.get(c, 0)
    }
    print("contigs compared:", len(set(old) | set(new)), "mismatching:", len(bad))
    for c, (o, n) in sorted(bad.items())[:20]:
        print("  MISMATCH", c, "manual", o, "star_confirm", n)
    print((run / "results/star_support_by_virus.tsv").read_text()[:2000])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

"""Run ONE grid condition end to end with the real product, score it, delete the big intermediates (PLAN SPIKE-02).

``python run_condition.py --row N`` (1-based row of conditions.tsv; ``grid.sbatch`` passes SLURM_ARRAY_TASK_ID) or
``--cond-id ID``; ``--dry-run`` prints every command and touches nothing. Steps:

 1. background subsample for (background, bg_rep): seeded ``spikein.py bg --frac`` built once under a file lock and
    cached in ROOT/bg; also the frozen called-cell list (top N_CALLED on-list barcodes of that subsample)
 2. ``spikein.py spike`` for this level/seed/cells (spiked conditions only)
 3. input = background followed by the spikes, concatenated gzip (valid gzip), named ``<cond_id>_R{1,2}.fastq.gz``
 4. ``viralscan run`` per kb arm (old, then flank), --host-filter starsolo; --star-confirm on the FIRST arm only (the
    STAR arm reads the host-filtered reads, which do not depend on the kb index); extract a small bundle; delete the run
 5. score_condition.py -> bundle/scores.tsv; delete the inputs; write ``DONE`` LAST.
Resumable: a finished condition is skipped, a finished arm (``stage_<arm>.done``) is skipped, a half-finished arm is
wiped (``viralscan run`` refuses a non-empty -o) and redone. Every failure raises; nothing is swallowed.
Deletion goes through ``rmtree_inside`` (refuses anything outside ROOT/cond), never an unquoted shell ``rm``.

Product behaviour that would break or mislead on spiked data, and what this runner does about it:
 * chemistry check (first 100k R1 reads, on-list match >= 0.5, UMI length): spikes are APPENDED after the background
   and have the background R1 length, so they are not in the first 100k; -x and --strand are passed explicitly.
 * read ids: paired ids must match and be unique; spikes are ``SPIKE:<class>:<i>`` (unique, never in the background).
 * barcodes must be on the chemistry on-list or kb/STAR drop them: spike barcodes are drawn from the on-list (--whitelist).
 * cell calling: the default emptyDrops needs star_tmp/Solo.out and varies with the Monte-Carlo seed; the runner freezes a
   HEADERLESS called list per background replicate (top N_CALLED barcodes) and passes ``--cell-calling external``; it
   asserts every spike barcode is in it (else the called-cell gate would silently be 0). Scoring does NOT use the
   ``detected`` column anyway: it reads per-barcode molecules from adata_multimap.h5ad.
 * viral_summary / multimap_evidence totals are over ALL barcodes (background viral reads count): kb scoring is
   restricted to the spike barcodes and the truth accessions. The unspiked replicate may not exist yet when a
   spiked condition finishes (array order), so the baseline is subtracted in ``aggregate.apply_baseline``, not here.
 * ``--gtf`` is mandatory for these manifest-less indexes (no build manifest: the bundled panel fallback is absent).
 * old and flank runs cannot share an output dir or --resume (index hash in the manifest): separate dirs, host filter
   repeated (the dominant cost); the STAR columns of the flank run are not requested.
 * spikes that map to the host under the pinned host filter are removed before kb: the limit includes that loss
   (counted in host_filtered/fragment_lineage.tsv.gz; not extracted here).
 * one cwd per run (Snakemake writes .snakemake to the cwd; concurrent array tasks must not share it).
 * UNVERIFIED: that spike read names survive the product's host-filter step into viral.bam (the smoke test checked this
   only for ``star_hostvirus align``); the scorer counts a spike absent from the BAM as absent, so a rename would show
   as zero STAR recovery, which the done-gate below turns into a loud failure for spiked conditions at level >= 100.
"""

import argparse
import fcntl
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grid  # noqa: E402
import score as spike_score  # noqa: E402
import score_condition  # noqa: E402

PY = str(grid.ENV_BIN / "python")
SPIKEIN = str(HERE / "spikein.py")
KEEP_FROM_RUN = (
    "results/viral_summary.tsv",
    "results/multimap_evidence.tsv",
    "results/virus_identity.tsv",
    "results/star_support_by_virus.tsv",
    "results/cell_calling_input.json",
    "results/chemistry_sanity.json",
    "count_audit.tsv",
    "config.yaml",
    "kb-python/run_info.json",
    "host_filtered/host_filter_audit.tsv",
    "star_hostvirus/Log.final.out",
)


def rmtree_inside(path, root=None):
    """rmtree that refuses anything not strictly inside ROOT/cond."""
    base = Path(root or grid.ROOT / "cond").resolve()
    p = Path(path).resolve()
    if base not in p.parents:
        sys.exit(f"refusing to delete {p}: not inside {base}")
    if p.exists():
        shutil.rmtree(p)


class Ctx:
    def __init__(self, cond: dict, dry: bool, root: Path | None = None):
        self.c, self.dry = cond, dry
        self.root = Path(root or grid.ROOT)
        self.bg = grid.BACKGROUNDS[cond["background"]]
        self.src = grid.SOURCES.get(cond["source"]) if cond["source"] else None
        self.id = cond["cond_id"]
        self.d = self.root / "cond" / self.id
        self.bgd = self.root / "bg" / f"{cond['background']}-r{cond['bg_rep']}"
        self.level, self.seed, self.cells = (
            int(cond["level"]),
            int(cond["seed"]),
            int(cond["cells"]),
        )
        self.env = {
            **os.environ,
            "PATH": f"{grid.ENV_BIN}:{os.environ.get('PATH', '')}",
            "PYTHONPATH": str(grid.SRC),
        }

    def sh(self, cmd, log=None, cwd=None):
        print("+ " + shlex.join(str(x) for x in cmd), flush=True)
        if self.dry:
            return
        out = open(log, "ab") if log else None  # noqa: SIM115
        try:
            subprocess.run([str(x) for x in cmd], check=True, env=self.env, cwd=cwd, stdout=out,
                           stderr=subprocess.STDOUT if out else None)  # fmt: skip
        finally:
            if out:
                out.close()

    def note(self, msg):
        print("# " + msg, flush=True)

    def onlist(self) -> str:
        if self.dry:
            return f"<onlist:{self.bg['technology']}>"
        sys.path.insert(0, str(grid.SRC))
        from viralscan.chemistry import onlist_path

        p = onlist_path(self.bg["technology"])
        if not p:
            sys.exit(f"no bundled on-list for {self.bg['technology']}")
        return p


def count_pairs(ctx: Ctx, r1: Path) -> int:
    """Pair count of R1, cached per background (zcat | wc -l, list form, no shell)."""
    cache = ctx.root / "bg" / f"{ctx.c['background']}.total_pairs"
    cache.parent.mkdir(parents=True, exist_ok=True)
    # one lock per BACKGROUND (stage_bg locks per replicate, so replicates race here); atomic cache write
    with open(cache.with_name(cache.name + ".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if cache.exists():
            return int(cache.read_text())
        z = subprocess.Popen(["zcat", str(r1)], stdout=subprocess.PIPE)
        w = subprocess.run(["wc", "-l"], stdin=z.stdout, capture_output=True, text=True, check=True)
        if z.wait() != 0:
            sys.exit(f"zcat {r1} failed")
        n = int(w.stdout.split()[0]) // 4
        tmp = cache.with_name(cache.name + ".part")
        tmp.write_text(str(n))
        os.replace(tmp, cache)
        return n


def stage_bg(ctx: Ctx):
    """Cached seeded background subsample + frozen called list (built once; array tasks serialize on a lock)."""
    b, d = ctx.bg, ctx.bgd
    cmd = [PY, "-I", SPIKEIN, "bg", "--r1", b["r1"], "--r2", b["r2"], "--seed", ctx.c["bg_rep"],
           "--out-prefix", d / "tmp"]  # fmt: skip
    if ctx.dry:
        ctx.note(f"background {d.name}: built once under a lock if {d / 'READY'} is absent")
        ctx.sh([*cmd, "--frac", f"<{grid.BG_TARGET_PAIRS}/total pairs of R1>"])
        ctx.note(f"called list = top {grid.N_CALLED} on-list barcodes -> {d / 'called_cells.txt'}")
        return
    d.mkdir(parents=True, exist_ok=True)
    with open(d.parent / f"{d.name}.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (d / "READY").exists():
            return
        total = count_pairs(ctx, b["r1"])
        frac = min(1.0, grid.BG_TARGET_PAIRS / total)
        ctx.sh([*cmd, "--frac", f"{frac:.8f}"], log=d / "bg.log")
        for m in ("R1", "R2"):
            os.replace(d / f"tmp_{m}.fastq.gz", d / f"bg_{m}.fastq.gz")
        import spikein

        onl = spikein.load_whitelist(ctx.onlist())
        cbs = spikein.top_barcodes(
            str(d / "bg_R1.fastq.gz"), b["cb_len"], grid.N_CALLED, onlist=onl
        )
        (d / "called_cells.txt").write_text("".join(c + "\n" for c in cbs))  # HEADERLESS
        (d / "READY").write_text(
            json.dumps({"total_pairs": total, "frac": frac, "seed": ctx.c["bg_rep"]})
        )


def stage_spike(ctx: Ctx):
    if ctx.src is None:
        return
    s, b = ctx.src, ctx.bg
    spikes = ctx.d / "spikes"
    cmd = [PY, "-I", SPIKEIN, "spike", "--bam", s["bam"], "--gtf", s["gtf"],
           "--bg-r1", ctx.bgd / "bg_R1.fastq.gz", "--out-dir", spikes, "--contig", s["contigs"],
           "--max-nh", s["max_nh"], "--levels", ctx.level, "--seeds", ctx.seed, "--cells", ctx.cells,
           "--whitelist", ctx.onlist(), "--cb-len", b["cb_len"], "--umi-len", b["umi_len"],
           "--r1-len", b["r1_len"], "--r2-len", b["r2_len"], "--source-orient", s["orient"],
           "--bg-orient", b["orient"]]  # fmt: skip
    truth = spikes / f"spike_L{ctx.level}_s{ctx.seed}.truth.tsv"
    if not ctx.dry and truth.exists() and (spikes / "SPIKES_OK").exists():
        return
    if not ctx.dry:
        rmtree_inside(spikes, ctx.root / "cond")
        (ctx.d / "logs").mkdir(parents=True, exist_ok=True)
    ctx.sh(cmd, log=ctx.d / "logs" / "spike.log")
    if ctx.dry:
        return
    called = set(ctx.bgd.joinpath("called_cells.txt").read_text().split())
    cbs = {ln.split("\t")[1] for ln in truth.read_text().splitlines()[1:]}
    if not cbs <= called:
        sys.exit(
            f"{ctx.id}: spike barcodes {sorted(cbs - called)[:3]}... are not in the frozen called list"
        )
    (spikes / "SPIKES_OK").write_text("ok\n")


def inputs(ctx: Ctx) -> tuple[Path, Path]:
    return tuple(ctx.d / "in" / f"{ctx.id}_R{m}.fastq.gz" for m in (1, 2))


def stage_inputs(ctx: Ctx):
    """Background then spikes, concatenated (concatenated gzip members are one valid gzip stream)."""
    for m, dest in zip((1, 2), inputs(ctx)):
        parts = [ctx.bgd / f"bg_R{m}.fastq.gz"]
        if ctx.src:
            parts.append(ctx.d / "spikes" / f"spike_L{ctx.level}_s{ctx.seed}_R{m}.fastq.gz")
        ctx.note(f"cat {' '.join(map(str, parts))} > {dest}")
        if ctx.dry or dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".part")
        with open(tmp, "wb") as out:
            for p in parts:
                with open(p, "rb") as fh:
                    shutil.copyfileobj(fh, out, 1 << 24)
        os.replace(tmp, dest)


def viralscan_cmd(ctx: Ctx, arm: str, out: Path) -> list:
    a, b = grid.ARMS[arm], ctx.bg
    r1, r2 = inputs(ctx)
    cmd = [PY, "-m", "viralscan.menu", "run", "--r1", r1, "--r2", r2, "-o", out, "-c", grid.CORES,
           "-i", a["idx"], "-t", a["t2g"], "--gtf", a["gtf"], "-x", b["technology"],
           "--strand", b["strand"], "--host-filter", "starsolo", "--host-index", grid.HOST_INDEX,
           "--cell-calling", "external", "--called-cells-file", ctx.bgd / "called_cells.txt",
           "--no-visual", "--yes"]  # fmt: skip
    if a["star_confirm"]:
        cmd += ["--star-confirm", "--star-hostvirus-index", grid.STAR_HV_INDEX]
    return cmd


def extract_kb(sample: Path, dest: Path, spike_cbs=None, min_level=0):
    """Non-zero viral entries of adata_multimap.h5ad -> kb_molecules.tsv.gz (barcode, gene_id, unique, ambiguous).

    Sanity gates (fail loudly, a silent empty table would score as "not reached"): at least one viral column must match
    virus_identity.tsv (matched on ``var['gene_id']``: the h5ad index is the kb matrix name, gene_id is a separate
    column); when ``spike_cbs`` is given and ``min_level`` >= 100, at least one spike barcode must be in obs_names.
    A trailing ``-<digit>`` on a barcode is stripped."""
    import anndata
    import pandas as pd

    ad = anndata.read_h5ad(sample / "kb-python/counts_unfiltered/adata_multimap.h5ad")
    viral = set(pd.read_csv(sample / "results/virus_identity.tsv", sep="\t")["gene_id"])
    gene = list(ad.var["gene_id"]) if "gene_id" in ad.var else list(ad.var_names)
    cols = [i for i, g in enumerate(gene) if g in viral]
    if not cols:
        sys.exit(
            f"extract_kb: no viral gene of virus_identity.tsv found in {sample} adata_multimap.h5ad"
        )
    obs = [re.sub(r"-\d+$", "", str(o)) for o in ad.obs_names]
    if spike_cbs and min_level >= 100 and not set(obs) & set(spike_cbs):
        sys.exit(
            f"extract_kb: none of {len(spike_cbs)} spike barcodes is in obs_names of {sample}: barcode mismatch"
        )
    u = ad.layers["counts_unique_viral"][:, cols].tocoo()
    a = ad.layers["counts_ambiguous_allocated"][:, cols].tocoo()
    cell: dict[tuple[int, int], list[float]] = {}
    for m, k in ((u, 0), (a, 1)):
        for r, c, v in zip(m.row, m.col, m.data):
            cell.setdefault((int(r), int(c)), [0.0, 0.0])[k] += float(v)
    with gzip.open(dest / "kb_molecules.tsv.gz", "wt") as fh:
        fh.write("barcode\tgene_id\tunique\tambiguous\n")
        for (r, c), (uu, aa) in sorted(cell.items()):
            fh.write(f"{obs[r]}\t{gene[cols[c]]}\t{uu:g}\t{aa:g}\n")


def extract_star_lines(lines, contigs: str, dest: Path):
    """Spike fates + non-spike per-contig read counts from ``samtools view viral.bam`` lines."""
    lines = list(lines)
    fate, on_contig = spike_score.star_fates(lines, contigs)
    any_n, prim_n = {}, {}
    for ln in lines:
        if ln.startswith("SPIKE"):
            continue
        f = ln.split("\t", 4)
        any_n[f[2]] = any_n.get(f[2], 0) + 1
        if not int(f[1]) & 0x100:
            prim_n[f[2]] = prim_n.get(f[2], 0) + 1
    dest.mkdir(parents=True, exist_ok=True)
    with open(dest / "star_spike_fates.tsv", "w") as fh:
        fh.write("name\ton_truth\tfate\n")
        names = set(fate) | on_contig
        fh.writelines(
            f"{n}\t{int(n in on_contig)}\t{fate.get(n, 'absent')}\n" for n in sorted(names)
        )
    with open(dest / "star_contig_reads.tsv", "w") as fh:
        fh.write("contig\tn_any\tn_primary\n")
        fh.writelines(f"{c}\t{n}\t{prim_n.get(c, 0)}\n" for c, n in sorted(any_n.items()))


def extract_star(ctx: Ctx, sample: Path, dest: Path):
    contigs = ctx.src["contigs"] if ctx.src else "VIRUS|NONE"
    p = subprocess.Popen(["samtools", "view", str(sample / "star_hostvirus/viral.bam")],
                         stdout=subprocess.PIPE, text=True, env=ctx.env)  # fmt: skip
    extract_star_lines(p.stdout, contigs, dest)
    if p.wait() != 0:
        sys.exit("samtools view viral.bam failed")


def run_arm(ctx: Ctx, arm: str):
    marker = ctx.d / f"stage_{arm}.done"
    if marker.exists():
        ctx.note(f"arm {arm}: done, skipping")
        return
    work = ctx.d / f"work_{arm}"
    out, sample = work / "out", work / "out" / ctx.id
    ctx.note(f"arm {arm}: wiping any half-finished {work}")
    if not ctx.dry:
        rmtree_inside(work, ctx.root / "cond")
        work.mkdir(parents=True)
        (ctx.d / "logs").mkdir(exist_ok=True)
    ctx.sh(viralscan_cmd(ctx, arm, out), log=ctx.d / "logs" / f"viralscan_{arm}.log", cwd=work)
    if ctx.dry:
        ctx.note(f"extract {sample} -> {ctx.d / 'bundle' / arm}; delete {work}; write {marker}")
        return
    dest = ctx.d / "bundle" / arm
    dest.mkdir(parents=True, exist_ok=True)
    for rel in KEEP_FROM_RUN:
        if (sample / rel).exists():
            shutil.copy2(sample / rel, dest / Path(rel).name)
    for rel in ("run_manifest.json", "run_complete.json"):
        if (out / rel).exists():
            shutil.copy2(out / rel, dest / rel)
    truth = ctx.d / "spikes" / f"spike_L{ctx.level}_s{ctx.seed}.truth.tsv"
    cbs = {ln.split("\t")[1] for ln in truth.read_text().splitlines()[1:]} if ctx.src else set()
    extract_kb(sample, dest, cbs, ctx.level)
    if grid.ARMS[arm]["star_confirm"]:
        extract_star(ctx, sample, ctx.d / "bundle" / "star")
    rmtree_inside(work, ctx.root / "cond")
    marker.write_text(time.strftime("%Y-%m-%dT%H:%M:%S\n"))


def finish(ctx: Ctx):
    cmd = [PY, "-I", HERE / "score_condition.py", "--bundle", ctx.d / "bundle", "--conditions",
           ctx.root / "conditions.tsv", "--cond-id", ctx.id]  # fmt: skip
    stem = ctx.d / "spikes" / f"spike_L{ctx.level}_s{ctx.seed}"
    if ctx.dry or not ctx.src:
        ctx.note(f"copy {stem}.truth.tsv and .json into the bundle (spiked conditions)")
    elif stem.with_suffix(
        ".truth.tsv"
    ).exists():  # absent only when a crashed finish already moved it
        shutil.copy2(stem.with_suffix(".truth.tsv"), ctx.d / "bundle" / "truth.tsv")
        shutil.copy2(stem.with_suffix(".json"), ctx.d / "bundle" / "spike_meta.json")
    ctx.sh(cmd, log=ctx.d / "logs" / "score.log")
    if ctx.dry:
        ctx.note(f"delete {ctx.d / 'in'} and spikes; write {ctx.d / 'DONE'} last")
        return
    rows = score_condition.read_tsv(ctx.d / "bundle" / "scores.tsv")
    star = next(r for r in rows if r["arm"] == "star")
    if ctx.src and ctx.level >= 100 and star["star_n_source_contig_any"] == "0":
        sys.exit(
            f"{ctx.id}: 0 STAR reads on the truth contig at level {ctx.level}: spike names lost or run broken"
        )
    rmtree_inside(ctx.d / "in", ctx.root / "cond")
    rmtree_inside(ctx.d / "spikes", ctx.root / "cond")  # truth + meta are in the bundle
    done = {"finished": time.strftime("%Y-%m-%dT%H:%M:%S"), "cond_id": ctx.id}
    (ctx.d / "DONE").write_text(json.dumps(done))


def run(cond: dict, dry=False, root=None):
    ctx = Ctx(cond, dry, root)
    if not dry and (ctx.d / "DONE").exists():
        print(f"{ctx.id}: DONE, skipping")
        return ctx
    if not dry:
        ctx.d.mkdir(parents=True, exist_ok=True)
    todo = [a for a in grid.ARMS if dry or not (ctx.d / f"stage_{a}.done").exists()]
    if todo:
        stage_bg(ctx)
        stage_spike(ctx)
        stage_inputs(ctx)
        for arm in grid.ARMS:
            run_arm(ctx, arm)
    finish(ctx)
    return ctx


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--row", type=int, help="1-based row of conditions.tsv (SLURM_ARRAY_TASK_ID)")
    g.add_argument("--cond-id")
    p.add_argument("--conditions", default=str(grid.ROOT / "conditions.tsv"))
    p.add_argument("--dry-run", action="store_true", help="print the commands, execute nothing")
    a = p.parse_args(argv)
    rows = (
        grid.load_conditions(a.conditions)
        if Path(a.conditions).exists()
        else grid.build_conditions()
    )
    if a.row:
        if not 1 <= a.row <= len(rows):
            sys.exit(f"--row {a.row} outside 1..{len(rows)}")
        cond = rows[a.row - 1]
    else:
        cond = next((r for r in rows if r["cond_id"] == a.cond_id), None)
        if cond is None:
            sys.exit(f"unknown condition {a.cond_id}")
    if not a.dry_run and not Path(a.conditions).exists():
        sys.exit(f"{a.conditions} missing: run grid.py first")
    run({k: str(v) for k, v in cond.items()}, dry=a.dry_run)


if __name__ == "__main__":
    main()

"""Host-subtraction loss of spiked reads (PLAN DISC-12-HH; protocol docs/plans/2026-10-22-discovery-benchmark.md §1 Q1, §4).

Question: what fraction of truth viral fragments does the product's STAR host filter remove? Two pools:
``outside_truth``  the host-homologous / EVE-like reads (real HHV-6B / EBV reads that also align to HOST| or another
                   virus; the grid's ``alignment_outside_truth_set`` drop, user decision Q16), ``--max-nh 20``;
``default``        the ordinary pool of the SPIKE-02 grid (same sources and NH policy), for the baseline number.
Each pool needs its own 30 conditions (own ROOT); nothing here submits anything.

``python hosthomolog.py grid --pool outside_truth``  write ROOT/conditions.tsv (5 levels x 3 seeds x c100 x 2 backgrounds)
``python hosthomolog.py run --row N | --cond-id ID [--dry-run]``   one condition (``hosthomolog.sbatch`` passes the array id)

A condition = the CACHED background of the SPIKE-02 grid (``grid.ROOT/bg/<bg>-r<seed>``, read-only, must be READY) + the
spikes written by ``spikein.py spike --pool ...`` (100 cells), then ONLY ``host_filter._starsolo_filter`` (the product's
code path, ``pinned`` params, then ``star-default`` as the protocol's second arm), then ``lost_truth_counts`` with
label = viral against the filter's ``fragment_lineage.tsv.gz``. Output per condition: ``ROOT/cond/<id>/host_loss.tsv``
(one row per STAR param set: pool size, spiked fragments, removed fraction D15, lost molecules D16, ``not_in_lineage``).
A zero or too-small pool fails loudly (spikein exits); ``not_in_lineage`` > 0 fails the row (not estimable).
Resumable: ``host_loss.tsv`` is written LAST; half-finished work is wiped.
"""

import argparse
import csv
import json
import os
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grid  # noqa: E402
import run_condition  # noqa: E402

ROOT = Path(os.environ.get("HOSTLOSS_ROOT", grid.W / "disc12_hosthomolog_2026-10-23"))
PARAM_SETS = (
    "pinned",
    "star-default",
)  # protocol §3: pinned primary, star-default second arm for Q1
MAX_NH = 20  # outside_truth pool: STAR hostvirus keeps up to 100 loci; 20 = the product host filter's multimap limit
TAGS = {"outside_truth": "HH", "default": "DP"}
CELLS, LAYOUT = 100, "c100"
COLUMNS = ["cond_id", "background", "source", "pool", "level", "seed", "layout", "cells", "bg_rep"]
OUT_COLUMNS = [
    "cond_id", "background", "pool", "level", "seed", "star_params", "pool_size", "spiked_fragments",
    "removed_fragments", "d15_loss_fraction", "truth_molecules", "lost_molecules", "d16_loss_fraction",
    "not_in_lineage", "input_pairs", "retained_pairs",
]  # fmt: skip


def build_conditions(pool: str = "outside_truth") -> list[dict]:
    if pool not in TAGS:
        sys.exit(f"unknown pool {pool!r}")
    return [
        dict(
            cond_id=f"{bg}-{TAGS[pool]}-L{level:03d}-s{seed}-{LAYOUT}",
            background=bg,
            source=b["source"],
            pool=pool,
            level=level,
            seed=seed,
            layout=LAYOUT,
            cells=CELLS,
            bg_rep=seed,
        )
        for bg, b in grid.BACKGROUNDS.items()
        for level in grid.LEVELS
        for seed in grid.SEEDS
    ]


def host_loss(truth_tsv, lineage, work: Path) -> dict:
    """``lost_truth_counts`` of a spike truth table (``name``/``cb``/``umi``) against a fragment lineage."""
    sys.path.insert(0, str(grid.SRC))
    from viralscan.scripts.host_filter import lost_truth_counts

    work.mkdir(parents=True, exist_ok=True)
    conv = work / "truth_for_loss.tsv"
    with open(truth_tsv, newline="") as fh, open(conv, "w") as out:
        out.write("read_id\tlabel\tmolecule_id\n")
        for r in csv.DictReader(fh, delimiter="\t"):
            out.write(f"{r['name']}\tviral\t{r['cb']}-{r['umi']}\n")
    return lost_truth_counts(str(conv), str(lineage), "viral")


def run_star_filter(ctx, param: str, out_dir: Path):
    """The product's host filter on the condition input (STAR only; stubbed in tests)."""
    sys.path.insert(0, str(grid.SRC))
    os.environ["PATH"] = f"{grid.ENV_BIN}:{os.environ.get('PATH', '')}"
    from viralscan.scripts import host_filter

    r1, r2 = run_condition.inputs(ctx)
    out_dir.mkdir(parents=True)
    host_filter._starsolo_filter(
        str(r1), str(r2), str(grid.HOST_INDEX), ctx.bg["technology"], ctx.onlist(), out_dir,
        str(out_dir / "R1.fastq.gz"), str(out_dir / "R2.fastq.gz"), grid.CORES, param,
    )  # fmt: skip


def stage_spike(ctx, pool: str):
    spikes = ctx.d / "spikes"
    s, b = ctx.src, ctx.bg
    if (spikes / "SPIKES_OK").exists() and not ctx.dry:
        return
    if not ctx.dry:
        run_condition.rmtree_inside(spikes, ctx.root / "cond")
        (ctx.d / "logs").mkdir(parents=True, exist_ok=True)
    max_nh = MAX_NH if pool == "outside_truth" else s["max_nh"]
    ctx.sh([run_condition.PY, "-I", run_condition.SPIKEIN, "spike", "--bam", s["bam"], "--gtf", s["gtf"],
            "--bg-r1", ctx.bgd / "bg_R1.fastq.gz", "--out-dir", spikes, "--contig", s["contigs"],
            "--max-nh", max_nh, "--pool", pool, "--levels", ctx.level, "--seeds", ctx.seed,
            "--cells", ctx.cells, "--whitelist", ctx.onlist(), "--cb-len", b["cb_len"],
            "--umi-len", b["umi_len"], "--r1-len", b["r1_len"], "--r2-len", b["r2_len"],
            "--source-orient", s["orient"], "--bg-orient", b["orient"]],
           log=ctx.d / "logs" / "spike.log")  # fmt: skip
    if not ctx.dry:
        (spikes / "SPIKES_OK").write_text("ok\n")


def run(cond: dict, dry=False, root=None, params=PARAM_SETS):
    root = Path(root or ROOT)
    ctx = run_condition.Ctx(cond, dry, root)
    ctx.bgd = (
        grid.ROOT / "bg" / f"{cond['background']}-r{cond['bg_rep']}"
    )  # the grid's cached background, read-only
    result = ctx.d / "host_loss.tsv"
    if not dry:
        if result.exists():
            print(f"{ctx.id}: done, skipping")
            return ctx
        if not (ctx.bgd / "READY").exists():
            sys.exit(
                f"{ctx.id}: cached background {ctx.bgd} is not READY (run the SPIKE-02 grid first)"
            )
        ctx.d.mkdir(parents=True, exist_ok=True)
    stage_spike(ctx, cond["pool"])
    run_condition.stage_inputs(ctx)
    if dry:
        for p in params:
            ctx.note(f"STAR host filter ({p}) -> lost_truth_counts -> {result}")
        return ctx
    meta = json.loads((ctx.d / "spikes" / f"spike_L{ctx.level}_s{ctx.seed}.json").read_text())
    pool_size = sum(meta["pool_class_mix"].values())
    truth = ctx.d / "spikes" / f"spike_L{ctx.level}_s{ctx.seed}.truth.tsv"
    rows = []
    for p in params:
        work = ctx.d / f"work_{p}"
        run_condition.rmtree_inside(work, root / "cond")
        run_star_filter(ctx, p, work)
        lc = host_loss(truth, work / "fragment_lineage.tsv.gz", work)
        audit = {
            r[0]: r[1] for r in csv.reader(open(work / "host_filter_audit.tsv"), delimiter="\t")
        }
        if lc["not_in_lineage"]:
            sys.exit(
                f"{ctx.id} [{p}]: {lc['not_in_lineage']} spiked reads absent from the lineage: not estimable"
            )
        rows.append(
            dict(
                cond_id=ctx.id, background=cond["background"], pool=cond["pool"], level=ctx.level,
                seed=ctx.seed, star_params=p, pool_size=pool_size, spiked_fragments=lc["d15_truth_fragments"],
                removed_fragments=lc["d15_removed_fragments"], d15_loss_fraction=lc["d15_loss_fraction"],
                truth_molecules=lc["d16_truth_molecules"], lost_molecules=lc["d16_lost_molecules"],
                d16_loss_fraction=lc["d16_loss_fraction"], not_in_lineage=lc["not_in_lineage"],
                input_pairs=audit.get("input"), retained_pairs=audit.get("retained_host_unmapped"),
            )
        )  # fmt: skip
        run_condition.rmtree_inside(work, root / "cond")
    run_condition.rmtree_inside(ctx.d / "in", root / "cond")
    tmp = result.with_name("host_loss.tsv.part")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, result)  # LAST
    shutil.copy2(truth, ctx.d / "truth.tsv")
    print(f"{ctx.id}: done {time.strftime('%Y-%m-%dT%H:%M:%S')}")
    return ctx


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("grid", help="write conditions.tsv and the slurm/ log dir")
    g.add_argument("--pool", choices=sorted(TAGS), default="outside_truth")
    g.add_argument("--out", default=None)
    r = sub.add_parser("run")
    m = r.add_mutually_exclusive_group(required=True)
    m.add_argument("--row", type=int, help="1-based row of conditions.tsv")
    m.add_argument("--cond-id")
    r.add_argument("--conditions", default=None)
    r.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)
    if a.cmd == "grid":
        rows = build_conditions(a.pool)
        out = Path(a.out or ROOT / "conditions.tsv")
        out.parent.mkdir(parents=True, exist_ok=True)
        (out.parent / "slurm").mkdir(exist_ok=True)
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(
            f"{len(rows)} conditions -> {out}; sbatch --array=1-{len(rows)}%6 hosthomolog.sbatch",
            file=sys.stderr,
        )
        return
    path = Path(a.conditions or ROOT / "conditions.tsv")
    if not path.exists():
        sys.exit(f"{path} missing: run `hosthomolog.py grid` first")
    rows = grid.load_conditions(path)
    if a.row:
        if not 1 <= a.row <= len(rows):
            sys.exit(f"--row {a.row} outside 1..{len(rows)}")
        cond = rows[a.row - 1]
    else:
        cond = next((x for x in rows if x["cond_id"] == a.cond_id), None)
        if cond is None:
            sys.exit(f"unknown condition {a.cond_id}")
    run({k: str(v) for k, v in cond.items()}, dry=a.dry_run)


if __name__ == "__main__":
    main()

"""Bounded exploratory read-cleanup pilot; package defaults remain unchanged."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import plant_anello_10x as planter

from viralscan.anello_align import _edge_len, _revcomp
from viralscan.evidence import _open_maybe_gzip, parse_flagged_target_bus
from viralscan.run_safety import sha256_file
from viralscan.scripts import read_filter as rf
from viralscan.scripts.host_filter import canonical_read_id


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def fresh_output(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite pilot output: {output}")
    return output


def records(path):
    with _open_maybe_gzip(str(path)) as handle:
        while True:
            rec = [handle.readline() for _ in range(4)]
            if not rec[0]:
                return
            if not all(rec) or not rec[0].startswith("@") or not rec[2].startswith("+"):
                raise ValueError("Malformed or truncated FASTQ")
            if len(rec[1].rstrip()) != len(rec[3].rstrip()):
                raise ValueError("FASTQ sequence/quality length mismatch")
            yield rec


def pairs(r1, r2, limit=None):
    seen = set()
    stream = itertools.zip_longest(records(r1), records(r2))
    if limit is not None:
        stream = itertools.islice(stream, limit)
    for a, b in stream:
        if a is None or b is None or canonical_read_id(a[0]) != canonical_read_id(b[0]):
            raise ValueError("FASTQ mate mismatch")
        name = canonical_read_id(a[0])
        if name in seen:
            raise ValueError(f"Duplicate fragment ID: {name}")
        seen.add(name)
        yield a, b


def sample(r1, r2, output, limit, geometry):
    output = fresh_output(output)
    n = 0
    with (
        gzip.open(output / "R1.fastq.gz", "wt", compresslevel=1) as a,
        gzip.open(output / "R2.fastq.gz", "wt", compresslevel=1) as b,
    ):
        for rec1, rec2 in pairs(r1, r2, limit):
            if len(rec1[1].rstrip()) != geometry:
                raise ValueError("Unsupported R1 geometry; barcode read must be exact CB+UMI")
            a.writelines(rec1)
            b.writelines(rec2)
            n += 1
    if n != limit:
        raise ValueError(f"Requested {limit} pairs but found {n}")
    write_json(
        output / "sample.json",
        {
            "pairs": n,
            "selection": "first N, not random",
            "r1_sha256": sha256_file(output / "R1.fastq.gz"),
            "r2_sha256": sha256_file(output / "R2.fastq.gz"),
        },
    )


def transform(r1, r2, output, mode, orientation="sense", minimum=31):
    output = Path(output)
    if mode not in ("off", "tso", "cleanup") or orientation not in ("sense", "antisense"):
        raise ValueError("Unsupported cleanup recipe")
    output = fresh_output(output)
    counts = Counter()
    if mode == "cleanup":
        with gzip.open(output / "tso_input.fastq.gz", "wt", compresslevel=1) as pre:
            for _, rec2 in pairs(r1, r2):
                seq, quality = rec2[1].rstrip(), rec2[3].rstrip()
                n = _edge_len(seq.upper(), rf.TSO_FULL)
                seq, quality = seq[n:], quality[n:]
                pre.writelines([rec2[0], seq + "\n", rec2[2], quality + "\n"])
        tail = (
            ["-a", "polyA=A{1000}X;min_overlap=15;max_error_rate=0"]
            if orientation == "sense"
            else ["-g", "polyT=XT{1000};min_overlap=15;max_error_rate=0"]
        )
        cmd = [
            sys.executable,
            "-m",
            "cutadapt",
            "-j",
            "1",
            "-e",
            "0.1",
            "-n",
            "2",
            *tail,
            "-a",
            "adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAX;min_overlap=12",
            "--info-file",
            str(output / "cutadapt_info.tsv"),
            "--json",
            str(output / "cutadapt.json"),
            "-o",
            str(output / "cleaned_R2.fastq.gz"),
            str(output / "tso_input.fastq.gz"),
        ]
        write_json(output / "command.json", cmd)
        with open(output / "cutadapt.log", "w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        cleaned = records(output / "cleaned_R2.fastq.gz")
        info_handle = open(output / "cutadapt_info.tsv")
        info = itertools.groupby(
            csv.reader(info_handle, delimiter="\t"), key=lambda row: canonical_read_id(row[0])
        )
    with (
        gzip.open(output / "R1.fastq.gz", "wt", compresslevel=1) as a,
        gzip.open(output / "R2.fastq.gz", "wt", compresslevel=1) as b,
        gzip.open(output / "fragment_lineage.tsv.gz", "wt") as lf,
    ):
        writer = csv.writer(lf, delimiter="\t")
        writer.writerow(
            ["read_id", "filter_decision", "reason", "original_length", "trim_left", "trim_right"]
        )
        for rec1, rec2 in pairs(r1, r2):
            name, seq = canonical_read_id(rec2[0]), rec2[1].rstrip()
            out = rec2
            left = 0
            if mode == "tso":
                n = _edge_len(seq.upper(), rf.TSO_FULL)
                left = n
                out = [rec2[0], seq[n:] + "\n", rec2[2], rec2[3].rstrip()[n:] + "\n"]
            elif mode == "cleanup":
                out = next(cleaned, None)
                if out is None or canonical_read_id(out[0]) != name:
                    raise ValueError("Trimming changed read identity/order")
                info_name, matches = next(info, (None, []))
                if info_name != name:
                    raise ValueError("Cutadapt audit differs from original read order")
                left = _edge_len(seq.upper(), rf.TSO_FULL) + sum(
                    int(row[3]) for row in matches if row[1] != "-1" and row[7] == "polyT"
                )
            result = out[1].rstrip()
            if (
                result != seq[left : left + len(result)]
                or out[3].rstrip() != rec2[3].rstrip()[left : left + len(result)]
            ):
                raise ValueError("Trimmed sequence/qualities are not a contiguous original slice")
            right = len(seq) - left - len(result)
            keep = mode == "off" or len(result) >= minimum
            reason = "short_after_trim" if not keep else "trimmed" if left or right else "untouched"
            counts["input"] += 1
            counts[reason] += 1
            writer.writerow(
                [name, "retained" if keep else "removed", reason, len(seq), left, right]
            )
            if keep:
                counts["retained"] += 1
                a.writelines(rec1)
                b.writelines(out)
        if mode == "cleanup" and next(cleaned, None) is not None:
            raise ValueError("Trimmer emitted extra records")
        if mode == "cleanup":
            if next(info, None) is not None:
                raise ValueError("Trimmer audit has extra records")
            info_handle.close()
    counts["removed"] = counts["input"] - counts["retained"]
    counts["r1_unchanged"] = True
    write_json(output / "audit.json", dict(counts))


def assess(config, sample_id, arm):
    cfg = json.loads(Path(config).read_text())
    root = Path(cfg["output"])
    sample_cfg = cfg["samples"][sample_id]
    out = root / sample_id / arm
    q = out / "quant"
    transcript_genes = {}
    for line in Path(cfg["t2g"]).read_text().splitlines():
        f = line.split("\t")
        transcript_genes[f[0]] = f[1]
    viral_rows = [line.split("\t") for line in Path(cfg["viral_t2g"]).read_text().splitlines()]
    viral = {row[0] for row in viral_rows}
    target = {row[0] for row in viral_rows if row[4] in sample_cfg["accessions"]}
    transcripts = (q / "transcripts.txt").read_text().splitlines()
    target_ids = {i for i, t in enumerate(transcripts) if t in target}
    if not target_ids:
        raise ValueError("Selected accessions have no indexed transcripts")
    viral_ids = {i for i, t in enumerate(transcripts) if t in viral}
    ec_map = {}
    for line in (q / "matrix.ec").read_text().splitlines():
        key, values = line.split("\t")
        ec_map[int(key)] = [int(x) for x in values.split(",")]
    with open(q / "bus.tsv") as lines:
        selected = parse_flagged_target_bus(
            lines, ec_map, target_ids, viral_ids, "host-conservative"
        )
    by_class = Counter(row["ambiguity_class"] for row in selected.values())
    keys = {
        (r["cb"], r["ub"]) for r in selected.values() if r["ambiguity_class"] == "unique_target"
    }
    with (
        open(out / "candidate_reads.fastq", "w") as reads,
        open(out / "candidates.tsv", "w") as table,
    ):
        writer = csv.writer(table, delimiter="\t")
        writer.writerow(["sample", "read_id", "barcode", "umi", "candidate_class"])
        for i, (a, b) in enumerate(pairs(out / "R1.fastq.gz", out / "R2.fastq.gz")):
            if i not in selected:
                continue
            row = selected[i]
            writer.writerow(
                [sample_id, canonical_read_id(a[0]), row["cb"], row["ub"], row["ambiguity_class"]]
            )
            if row["ambiguity_class"] == "unique_target":
                reads.writelines([f"@{sample_id}|{row['cb']}_{row['ub']}_{i}\n", b[1], "+\n", b[3]])
    write_json(
        out / "candidate_summary.json",
        {
            "sample": sample_id,
            "arm": arm,
            "read_classes": dict(by_class),
            "unique_target_raw_cb_umi_keys": len(keys),
            "barcode_correction": False,
            "cell_calling": False,
            "claim": "candidate read mapping; not production molecule counts or sensitivity",
        },
    )


def truth(output):
    import random

    output = fresh_output(output)
    rng = random.Random(20261014)
    genome = "".join(rng.choice("ACGT") for _ in range(5000))
    rows = []
    with gzip.open(output / "R1.fastq.gz", "wt") as a, gzip.open(output / "R2.fastq.gz", "wt") as b:
        for i in range(120):
            body = "CG" + genome[i * 20 : i * 20 + 58] + "CG"
            kind = ["clean", "polyA", "TSO_polyA", "polyT", "internal_A", "templated_terminal_A"][
                i % 6
            ]
            seq = body
            if kind == "polyA":
                seq, _ = planter.polya_read(body, 0, len(body), rng, len(body) + 25)
            elif kind == "TSO_polyA":
                seq = rf.TSO_FULL + body + "A" * 25
            elif kind == "polyT":
                seq = _revcomp(body + "A" * 25)
            elif kind == "internal_A":
                seq = body[:32] + "A" * 20 + body[32:]
            elif kind == "templated_terminal_A":
                seq = body + "A" * 25
            name = f"truth_{i}"
            planter.fastq(a, name, "ACGTACGTACGTACGT" + planter.umi(rng, 12))
            planter.fastq(b, name, seq)
            rows.append([name, kind, len(body), kind == "templated_terminal_A"])
        for i, seq in enumerate(
            ["A" * 90, rf.TSO_FULL + "A" * 60, "T" * 90, genome[:20] + "A" * 50]
        ):
            name = f"short_{i}"
            planter.fastq(a, name, "ACGTACGTACGTACGT" + planter.umi(rng, 12))
            planter.fastq(b, name, seq)
            rows.append([name, "too_short_or_reagent", 0 if i < 3 else 20, False])
    with open(output / "truth.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["read_id", "kind", "templated_body_length", "terminal_A_is_templated"])
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("sample")
    p.add_argument("r1")
    p.add_argument("r2")
    p.add_argument("output")
    p.add_argument("--pairs", dest="limit", type=int, required=True)
    p.add_argument("--geometry", type=int, required=True)
    p = sub.add_parser("transform")
    p.add_argument("r1")
    p.add_argument("r2")
    p.add_argument("output")
    p.add_argument("--mode", choices=["off", "tso", "cleanup"], required=True)
    p.add_argument("--orientation", choices=["sense", "antisense"], default="sense")
    p.add_argument("--minimum", type=int, default=31)
    p = sub.add_parser("assess")
    p.add_argument("config")
    p.add_argument("sample_id")
    p.add_argument("arm")
    p = sub.add_parser("truth")
    p.add_argument("output")
    args = vars(parser.parse_args())
    command = args.pop("command")
    globals()[command](**args)


if __name__ == "__main__":
    main()

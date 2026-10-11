"""Append findings to findings_D_remaining.tsv from a JSON file: add.py batch.json
(batch.json = list of 8-element lists: id, file:line, function, class, severity, evidence, fix, confidence)."""
import csv
import json
import sys

rows = json.load(open(sys.argv[1]))
with open("analysis/audit_2026-10/findings_D_remaining.tsv", "a", newline="") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    for r in rows:
        assert len(r) == 8, r[0]
        w.writerow([str(c).replace("\t", " ").replace("\n", " ") for c in r])
print("appended", len(rows))

"""reference_strategy._presence samples only head+tail (plain) or head (gz); and
_validate_provenance accepts bool counts and fetched+missing < expected."""
import gzip
import tempfile
from pathlib import Path

from viralscan import reference_strategy as rs

tmp = Path(tempfile.mkdtemp())
filler = "chr1\tx\tgene\t1\t2\t.\t+\t.\tgene_id \"HUMAN\";\n" * 100_000  # ~4.3 MB
target = "NC_007605\tx\tgene\t1\t2\t.\t+\t.\tgene_id \"EPSTEIN_HHV4_BZLF1\";\n"
for name, text, opener in (
    ("mid.gtf", filler + target + filler, open),
    ("tail_gz.gtf.gz", filler + target, lambda p, m: gzip.open(p, "wt")),
):
    p = tmp / name
    with opener(p, "w") as h:
        h.write(text)
    print(name, "contains EBV gene ->", rs._presence(p))

base = {"build_commands": {k: "cmd" for k in rs.REQUIRED_BUILD_COMMANDS}}
for label, vp in (
    ("bool counts", {"anellovirus_expected_count": True, "anellovirus_fetched_count": True,
                     "anellovirus_missing_count": False}),
    ("fetched+missing=10 of expected=2021", {"anellovirus_expected_count": 2021,
                                            "anellovirus_fetched_count": 5,
                                            "anellovirus_missing_count": 5}),
):
    try:
        rs._validate_provenance({**base, "viral_panel": vp})
        print(label, "-> accepted")
    except rs.BenchmarkContractError as exc:
        print(label, "-> rejected:", exc)

"""Mocked-HTTP probes of viralscan.scripts.ncbi_fetch (no network). Run with PYTHONPATH=src."""
import os, sys, tempfile, time, requests
from pathlib import Path
from unittest import mock
from viralscan.scripts import ncbi_fetch as nf

class R:
    def __init__(self, code, text): self.status_code, self.text = code, text
nf.time.sleep = lambda s: None  # no real sleeping
nf._RATE_INTERVAL_S = nf._RATE_INTERVAL_WITH_KEY_S = 0

GB = """LOCUS       NC_000001               100 bp    DNA     linear   VRL 01-JAN-2000
VERSION     NC_000001.2
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="Foo virus"
     CDS             1..30
                     /gene="a"
                     /product="a long product name that
                     continues here"
ORIGIN
        1 acgt
//
"""
def p(tag, *a): print(f"[{tag}]", *a)

# 1. API key / email leak in exception text on connection failure (real requests error text)
nf.EUTILS_BASE = "http://127.0.0.1:9/entrez/eutils/efetch.fcgi"   # closed port, loopback only
try:
    nf._efetch("NC_000001.1", "fasta", "me@lumc.nl", "SECRETKEY123")
except nf.NCBIFetchError as e:
    p("leak", "SECRETKEY123" in str(e), "me%40lumc.nl" in str(e) or "me@lumc.nl" in str(e)); p("msg", str(e)[:300])

# 2. truncated 200 FASTA is accepted and cached
def fake_get(url, params=None, timeout=None):
    if params["rettype"] == "fasta": return R(200, ">NC_000001.2 Foo\nACGTACGT")   # no trailing newline, truncated vs LOCUS 100 bp
    return R(200, GB)
with tempfile.TemporaryDirectory() as d, mock.patch("requests.get", fake_get):
    f, g = nf._fetch_one("NC_000001", Path(d), "a@b.c", None)
    p("trunc", "accepted+cached; FASTA len", sum(len(l) for l in f.read_text().splitlines()[1:]), "vs LOCUS 100")
    # 3. FASTA records glued when chunk lacks trailing newline; duplicate accession accepted
    out = Path(d) / "out"
    fa, gt = nf.fetch_reference(["NC_000001", "NC_000001"], out, email="a@b.c", cache_dir=d)
    p("dup-accession", fa.read_text().count(">"), "FASTA records;", gt.read_text().count("exon"), "exon rows")
    fa2, _ = nf.fetch_reference(["NC_000001", "NC_000001.2"], out, email="a@b.c", cache_dir=d)
    p("glue", repr(fa2.read_text()[:60]))
    # 4. cached accession still requires an e-mail
    try: nf.fetch_reference(["NC_000001"], out, email=None, cache_dir=d)
    except nf.NCBIFetchError as e: p("cached-needs-email", str(e)[:60])
    # 5. multi-line /product truncated to the first physical line
    p("product", [l.split("product")[1][:40] for l in g.read_text().splitlines()[:1]])

# 6. 429 forever -> number of sleeps / calls; 404 -> immediate
calls = []
def f429(url, params=None, timeout=None): calls.append(1); return R(429, "slow down")
with mock.patch("requests.get", f429):
    try: nf._efetch("NC_000001.1", "fasta", "a@b.c", None)
    except nf.NCBIFetchError as e: p("429", len(calls), "calls;", str(e)[:80])
def f200empty(url, params=None, timeout=None): return R(200, "")
with mock.patch("requests.get", f200empty):
    try: nf._efetch("NC_000001.1", "fasta", "a@b.c", None)
    except nf.NCBIFetchError as e: p("200-empty", str(e)[:90])
# 7. unreachable cache dir -> traceback type
try: nf._fetch_one("NC_000001", Path("/proc/nope/cache"), "a@b.c", None)
except Exception as e: p("unwritable-cache", type(e).__name__)
# 8. accession validation
for a in ["NC_001664", "NC_001664.1", "nc_001664.1", "NZ_CP012345.1", "KX020937.1", "NC_0016 64", "NC_1.1.1", "../etc"]:
    try: nf._validate_accession(a); ok=True
    except nf.NCBIFetchError: ok=False
    p("acc", a, ok)

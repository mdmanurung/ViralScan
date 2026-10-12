"""build-ref calls `kb ref` WITHOUT --overwrite. kb_python 0.30.2 (ref.py:642-674) then skips cDNA/t2g generation and
`kallisto index` when cdna.fa / index.idx already exist, so a rebuild into a non-empty --output keeps the OLD index
while combined.fa/gtf and reference_manifest.json are rewritten. No network, no kb run: _run_kb_ref and the fetchers
are stubbed. Also shows Ensembl cDNA cache reuse across releases. Run with PYTHONPATH=src."""
import gzip, tempfile
from pathlib import Path
from unittest import mock
from viralscan.scripts import build_reference as br, ncbi_fetch as nf

with tempfile.TemporaryDirectory() as d:
    d = Path(d); out = d / "out"
    (d / "h.fa.gz").write_bytes(gzip.compress(b">ENST1.1 cdna gene:ENSG1.1\n" + b"ACGT" * 50 + b"\n"))
    vf = d / "v.fasta"; vf.write_text(">NC_1.1 v\n" + "ACGTTGCA" * 40 + "\n")
    vg = d / "v.gtf"; vg.write_text('NC_1.1\tNCBI\texon\t1\t90\t.\t+\t0\tgene_id "NC_1.1_a"; transcript_id "NC_1.1_a_t1";\n')
    seen = {}
    def fake_kb(cmd, out_dir):
        seen["cmd"] = cmd; return {}
    with mock.patch.object(br, "fetch_host_cdna", lambda *a, **k: (d / "h.fa.gz", d / "h.fa.gz")), \
         mock.patch.object(nf, "fetch_reference", lambda *a, **k: (vf, vg)), \
         mock.patch.object(br, "_run_kb_ref", fake_kb), \
         mock.patch.object(br.shutil, "which", lambda n: "/bin/kb"):
        try:
            br.build_combined_reference("human", ["NC_1.1"], out, email="a@b.c", run_kb_ref=True, mask=False)
        except Exception as e:
            print("build raised (index file absent in stub):", type(e).__name__, str(e)[:80])
    print("kb ref cmd:", " ".join(seen.get("cmd", [])))
    print("'--overwrite' in cmd:", "--overwrite" in seen.get("cmd", []))

# Ensembl cDNA cache: filename is release-agnostic, GTF filename carries the release
with tempfile.TemporaryDirectory() as d:
    d = Path(d); dl = []
    def fake_dl(url, dest, **k): dl.append(url); dest.write_bytes(b"x"); return dest
    def listing(sp, base, retries=3):
        return ["Homo_sapiens.GRCh38.cdna.all.fa.gz"] if "/fasta/" in base else ["Homo_sapiens.GRCh38.%s.gtf.gz" % base.split("release-")[1].split("/")[0]]
    for rel in (116, 117):
        with mock.patch.object(br, "_ensembl_release", lambda r=rel: r), mock.patch.object(br, "_list_ensembl_files", listing), mock.patch.object(br, "_download", fake_dl):
            br.fetch_host_cdna("human", d / f"o{rel}", cache_dir=d / "c")
    print("downloads:", [u.split("/pub/")[1] for u in dl])

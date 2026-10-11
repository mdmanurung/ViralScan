"""read_t2g: one short row flips the whole t2g to legacy; manifest catalogue drift is warn-only; normal manifest round trip."""
import sys, logging, json
from pathlib import Path
from viralscan.virus_identity import read_t2g, write_build_manifest, build_identity_table, manifest_path_for_index
logging.basicConfig(level=logging.WARNING)
w=Path("/tmp/claude-149488/-exports-archive-hg-funcgenom-research-mdmanurung-ViralScan/a0c3f363-e67e-4594-9357-4e376c79d18a/scratchpad/auditA/t2g"); w.mkdir(exist_ok=True)
(w/"t2g.txt").write_text("T1\tG1\t\t\tNC_007605.1\nT2\tG2\t\t\tNC_007605.1\nSHORT\tG3\n")
r=read_t2g(w/"t2g.txt"); print("legacy flag with one 2-col row:",r.legacy,"accessions kept:",dict(r.accession))
# CAT-06 drift: wrong sha only warns
(w/"t2g2.txt").write_text("T1\tG1\t\t\tNC_007605.1\n")
idx=w/"index.idx"; write_build_manifest(idx,["H1"],["G1"])
mp=manifest_path_for_index(idx); d=json.loads(mp.read_text()); d["catalogue_sha256"]="0"*64; mp.write_text(json.dumps(d))
t=build_identity_table(w/"t2g2.txt",["G1"],build_manifest=mp)
print("drifted manifest still builds table:",len(t.genes),"genes; virus:",t.genes[0].virus_name)

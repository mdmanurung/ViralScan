"""Real-data reconciliation: does the in-package hostresponse reproduce the
external scripts' AUC pattern on the showcase EBV run? (advisor's discriminating check)"""
import tempfile, os, warnings
warnings.filterwarnings("ignore")
import pandas as pd
from viralscan.scripts.hostresponse import run_hostresponse

BASE = "results/hostresponse_ebv_matched"
HOST = f"{BASE}/host_only_matched.h5ad"
VIRUS = f"{BASE}/ebv_burden_matched.h5ad"

def run(label, depth_match, tag):
    out = tempfile.mkdtemp(prefix=f"rec_{tag}_")
    atxt = os.path.join(out, "analysis.txt")
    open(atxt, "w").write("Epstein-Barr virus\n")
    run_hostresponse(
        virus_h5ad=VIRUS, host_h5ad=HOST, viral_accessions_file=atxt, out_dir=out,
        use_hvg=True, n_stab_iter=20, stab_min_prob=0.6, detection_threshold=10,
        label=label, depth_match=depth_match, control_mito=True,
    )
    m = pd.read_csv(os.path.join(out, "hostresponse_metrics.csv")).iloc[0]
    return m

print(f"{'config':<22}{'model_auc':>10}{'depth_alone':>13}{'n_pos':>8}")
for label, dm, tag in [("raw", False, "raw"), ("cpm", False, "cpm"), ("raw", True, "match")]:
    try:
        m = run(label, dm, tag)
        print(f"{f'{label},depth_match={dm}':<22}{m.get('auc_mean',float('nan')):>10.3f}"
              f"{m.get('depth_alone_auc_mean',float('nan')):>13.3f}{int(m['n_positive']):>8}")
    except Exception as e:
        print(f"{label},dm={dm}: FAILED {e}")

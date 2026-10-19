"""Paired original/new preprocessing fits; expected scientific differences, immutable input."""
import hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parent;W=Path('/tmp/viralscan-shifted-log-fit');C=Path('/tmp/viralscan-shifted-log-candidate');B=Path('/exports/archive/hg-funcgenom-research/mdmanurung/vs_work/main')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
if len(sys.argv)>1:
 sys.path.insert(0,sys.argv[1]);from viralscan.scripts import hostresponse as hr
 assert Path(hr.__file__).resolve().is_relative_to(Path(sys.argv[1]).resolve())
 hr.run_hostresponse(str(W/'inputs/virus.h5ad'),str(W/'inputs/host.h5ad'),str(W/'inputs/ids.txt'),str(W/'out'),use_hvg=False,seeds=[0,1],n_stab_iter=3,detection_threshold=1,cv_mode='cell',do_enrichment=False)
else:
 assert not W.exists();(W/'inputs').mkdir(parents=True)
 sys.path[:0]=[str(C/'src'),str(C)]
 from tests.test_hostresponse import _make_host_adata,_make_virus_adata
 _make_host_adata().write_h5ad(W/'inputs/host.h5ad');_make_virus_adata().write_h5ad(W/'inputs/virus.h5ad');(W/'inputs/ids.txt').write_text('VIRUS_A\nVIRUS_B\n')
 inputs={str(f):sha(f) for f in (W/'inputs').iterdir()}
 env=dict(os.environ,PATH='/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_test_full/bin:'+os.environ['PATH'],PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',NUMBA_CACHE_DIR='/tmp/viralscan-shifted-log-cache/numba',MPLCONFIGDIR='/tmp/viralscan-shifted-log-cache/mpl',XDG_CACHE_HOME='/tmp/viralscan-shifted-log-cache/xdg')
 reports={}
 for arm,root in [('original_a',B),('original_b',B),('candidate_a',C),('candidate_b',C)]:
  if (W/'out').exists():shutil.rmtree(W/'out')
  dest=P/'paired_fit'/arm;assert not dest.exists();dest.parent.mkdir(exist_ok=True)
  with (dest.parent/(arm+'.stdout')).open('wb') as out,(dest.parent/(arm+'.stderr')).open('wb') as err:
   subprocess.run([sys.executable,__file__,str(root/'src')],env=env,stdout=out,stderr=err,check=True)
  shutil.copytree(W/'out',dest);reports[arm]={str(f.relative_to(dest)):sha(f) for f in dest.rglob('*') if f.is_file()}
  assert inputs=={f:sha(Path(f)) for f in inputs}
 assert reports['original_a']==reports['original_b'],'original A/A not deterministic'
 assert reports['candidate_a']==reports['candidate_b'],'candidate A/A not deterministic'
 assert reports['original_a'].keys()==reports['candidate_a'].keys()
 import pandas as pd
 summaries={}
 for arm in ['original_a','candidate_a']:
  d=P/'paired_fit'/arm;m=pd.read_csv(d/'hostresponse_metrics.csv');w=pd.read_csv(d/'VIRUS_A_gene_weights.csv');stab=pd.read_csv(d/'VIRUS_A_stability.csv')
  summaries[arm]=dict(auc_mean=m['auc_mean'].tolist(),weight_absolute_sum=float(w.weight_mean.abs().sum()),stability_sum=float(stab.stab_prob.sum()),stable_genes=int(stab.stable.sum()))
 result=dict(input_sha256=inputs,raw_manifests=reports,original_AA_raw_equal=True,candidate_AA_raw_equal=True,different_files=[f for f in reports['original_a'] if reports['original_a'][f]!=reports['candidate_a'][f]],summaries=summaries,interpretation='intentional new preprocessing; bounded synthetic technical evidence, no scientific improvement claim')
 (P/'paired_fit_comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summaries,indent=2))

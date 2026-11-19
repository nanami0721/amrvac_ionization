#!/usr/bin/env python3
"""Paired legacy/midpoint transport runs from already isolated source snapshots."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--plumbing-only',action='store_true',help='Check default forward and opt-in reverse against retained full-pair outputs');p.add_argument('--face-probe',action='store_true',help='HLL+MC two-rank pairs only; use analytic coarse MPI checks');p.add_argument('workdir',type=Path);p.add_argument('--legacy-source',type=Path,required=True);p.add_argument('--candidate-source',type=Path,required=True);p.add_argument('--reference-results',type=Path,required=True)
a=p.parse_args();repo=Path(__file__).resolve().parents[3];r=a.workdir.resolve()
assert r!=repo and repo not in r.parents
r.mkdir(parents=True,exist_ok=False);results={};status={}
guard=repo/'tests/yinyang/interface_transport/convergence/memory_guard.py'
# Compare the same complete physical configuration in both directions.
par=(a.reference_results/'yy-p0-s2-d1-hllmc/amrvac.par').read_text()
prior=json.loads((a.reference_results/'yy-p0-s2-d1-hllmc-memory/summary.json').read_text())
estimate=1.25*max(prior['peak_rss_mib'],prior['peak_footprint_mib'])
assert estimate<950
(r/'memory-prediction.json').write_text(json.dumps(dict(estimate_mib=estimate,ceiling_mib=950)))
(r/'driver.py').write_bytes(Path(__file__).read_bytes())
def job(name,cmd,cwd,env,limit=1024):
 with (r/(name+'.log')).open('w') as log:
  proc=subprocess.run([sys.executable,str(guard),'--out',str(r/(name+'-memory')),'--limit-mib',str(limit),'--',*map(str,cmd)],cwd=cwd,env=env,stdout=log,stderr=log)
 status[name]=proc.returncode;(r/'status.json').write_text(json.dumps(status,indent=2));print(name,proc.returncode,flush=True)
 if proc.returncode:raise SystemExit(proc.returncode)
for label,source in [('legacy',a.legacy_source),('midpoint',a.candidate_source)]:
 source=source.resolve();assert source!=repo and repo not in source.parents
 env=dict(os.environ,AMRVAC_DIR=str(source),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PATH='/usr/local/openmpi/bin:'+os.environ['PATH'])
 (r/(label+'-manifest.json')).write_text(json.dumps({str(f.relative_to(source)):hashlib.sha256(f.read_bytes()).hexdigest() for f in source.joinpath('src').rglob('*.t')},indent=2))
 d=r/(label+'-build');d.mkdir();shutil.copy2(repo/'tests/yinyang/ct_transport/mod_usr.t',d/'mod_usr.t');(d/'amrvac.par').write_text(par)
 job(label+'-setup',[source/'setup.pl','-d=3','-arch=debug'],d,env)
 job(label+'-build',['make','F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096','-j1'],d,env,4096)
 cases=[(method,direction,2) for method in (['hllmc'] if a.face_probe else ['base','hllmc']) for direction in [1,-1]]
 if a.plumbing_only:
  assert a.face_probe
  cases=[('hllmc',1 if label=='legacy' else -1,2)]
 if label=='midpoint' and not a.face_probe:cases.append(('hllmc',1,4))
 for method,direction,ranks in cases:
  name=f'{label}-{method}-d{direction}-np{ranks}';out=r/name;out.mkdir();shutil.copy2(d/'amrvac',out/'amrvac');config=par.replace('direction=1',f'direction={direction}')
  if a.face_probe and label=='midpoint':config=config.replace('&methodlist','&methodlist\n yy_ct_face_moments=T')
  if method=='base':config=config.replace("20*'hll'","20*'tvdlf'").replace("20*'woodward'","20*'vanleer'")
  (out/'amrvac.par').write_text(config);job(name,['mpirun','--oversubscribe','-np',ranks,'./amrvac'],out,env)
  assert (r/(name+'-memory/command.err')).stat().st_size==0
  h=[[float(x) for x in line.split()] for line in (out/'transport-history.txt').read_text().splitlines()];assert len(h)==301 and h[-1][0]==300 and abs(h[-1][1]-.6)<1e-12
  hashes={}
  for prefix in ['transport','faces']:
   lines=sorted(line for f in out.glob(prefix+'-rank-*.txt') for line in f.read_text().splitlines());hashes[prefix]=hashlib.sha256('\n'.join(lines).encode()).hexdigest()
  results[name]=dict(initial=h[0],final=h[-1],l2=(h[-1][15]/h[-1][16])**.5,max_error=h[-1][6]/.001,max_active_divergence=max(row[2] for row in h),hashes=hashes)
  (r/'results.json').write_text(json.dumps(results,indent=2))
  if ranks==4:assert hashes==results[f'{label}-{method}-d{direction}-np2']['hashes']
print('Paired representation transport tests complete',flush=True)

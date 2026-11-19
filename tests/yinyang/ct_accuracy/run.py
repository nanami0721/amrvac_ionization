#!/usr/bin/env python3
"""Build/run analytic Yin-Yang CT transfer diagnostics in an external source copy."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('workdir',type=Path);p.add_argument('--pilot',action='store_true');p.add_argument('--cache-source',type=Path,required=True);a=p.parse_args()
s=Path(__file__).resolve().parents[3];r=a.workdir.resolve();b=r/'source'
if r==s or s in r.parents:p.error('external directory required')
r.mkdir(parents=True,exist_ok=False)
shutil.copytree(s,b,ignore=shutil.ignore_patterns('.git','lib','__pycache__'))
shutil.copytree(a.cache_source/'lib',b/'lib')
changed=[]
for f in b.joinpath('src').rglob('*.t'):
 prior=a.cache_source/f.relative_to(b)
 if not prior.exists() or f.read_bytes()!=prior.read_bytes():
  changed.append(str(f.relative_to(b)));f.touch()
  for lib in b.joinpath('lib').iterdir():
   if lib.is_dir():
    for suffix in ('.f','.o'):(lib/(f.stem+suffix)).unlink(missing_ok=True)
(r/'cache-invalidation.json').write_text(json.dumps(changed))
files=[*s.joinpath('src').rglob('*.t'),s/'src/makefile',*Path(__file__).parent.glob('*')]
(r/'manifest.json').write_text(json.dumps({str(f.relative_to(s)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files if f.is_file()},indent=2))
(r/'source.diff').write_bytes(subprocess.check_output(['git','diff','HEAD'],cwd=s));(r/'commit.txt').write_bytes(subprocess.check_output(['git','rev-parse','HEAD'],cwd=s))
guard=b/'tests/yinyang/interface_transport/convergence/memory_guard.py'
env=dict(os.environ,AMRVAC_DIR=str(b),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PATH='/usr/local/openmpi/bin:'+os.environ['PATH'])
status={};results={}
def job(name,cmd,cwd,limit=1024):
 with (r/(name+'.log')).open('w') as log:
  proc=subprocess.run([sys.executable,str(guard),'--out',str(r/(name+'-memory')),'--limit-mib',str(limit),'--',*map(str,cmd)],cwd=cwd,env=env,stdout=log,stderr=log)
 status[name]=proc.returncode;(r/'status.json').write_text(json.dumps(status,indent=2));print(name,proc.returncode,flush=True)
 if proc.returncode:raise SystemExit(proc.returncode)
compiler='F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096'
par=(s/'tests/yinyang/ct_preflight/amrvac.par').read_text().replace('it_max=5','it_max=0')
par=par.replace('&methodlist','&methodlist\n ct_stage_split=T\n yy_ct_experimental=T')
par=par.replace('&meshlist','&meshlist\n max_blocks=128').replace('xprobmin1=1.d0','xprobmin1=.75d0').replace('xprobmax1=2.d0','xprobmax1=2.25d0')
for arch in (['debug'] if a.pilot else ['debug','default']):
 d=r/arch;d.mkdir();shutil.copy2(b/'tests/yinyang/ct_accuracy/mod_usr.t',d/'mod_usr.t');(d/'amrvac.par').write_text(par)
 job('setup-'+arch,[b/'setup.pl','-d=3','-arch='+arch],d)
 job('build-'+arch,['make',compiler,'-j1'],d,4096)
 cases=[(8,2,8)] if a.pilot else [(8,2,8),(12,2,8),(14,2,8),(14,2,16),(8,4,8)]
 if arch=='default':cases=[(8,2,8),(12,2,8),(14,2,8)]
 for nr,ranks,quad in cases:
  name=f'{arch}-n{nr}-np{ranks}-q{quad}';out=r/name;out.mkdir();shutil.copy2(d/'amrvac',out/'amrvac')
  config=par.replace('domain_nx1=16',f'domain_nx1={nr}').replace('domain_nx2=16',f'domain_nx2={2*nr}').replace('domain_nx3=48',f'domain_nx3={6*nr}')
  for c in [1,2,3]:config=config.replace(f'block_nx{c}=8',f'block_nx{c}={nr}')
  config+=f'\n&accuracy_test\n quadrature_order={quad}\n/\n';(out/'amrvac.par').write_text(config)
  if nr>8:
   coarse=max(n for n in [8,12,14] if n<nr and (r/f'{arch}-n{n}-np2-q8-memory/summary.json').exists())
   prior=json.loads((r/f'{arch}-n{coarse}-np2-q8-memory/summary.json').read_text());estimate=1.25*max(prior['peak_rss_mib'],prior['peak_footprint_mib'])*((nr+4)/(coarse+4))**3
   (out/'memory-prediction.json').write_text(json.dumps(dict(estimate_mib=estimate,ceiling_mib=950)))
   if estimate>950:raise RuntimeError('analytic refinement memory forecast exceeds ceiling')
  job(name,['mpirun','--oversubscribe','-np',ranks,'./amrvac'],out)
  mem=r/(name+'-memory');assert (mem/'command.err').stat().st_size==0
  rows=[line.split() for line in (mem/'command.log').read_text().splitlines() if line.startswith(('EDGE_ACCURACY ','GHOST_ACCURACY '))]
  assert len(rows)==28
  results[name]=[dict(kind=row[0],field=int(row[1]),mode=int(row[2]),values=list(map(float,row[3:]))) for row in rows]
  (r/'results.json').write_text(json.dumps(results,indent=2))
print('Analytic CT attribution cases completed',flush=True)

#!/usr/bin/env python3
"""Build/run the experimental Yin-Yang CT weak-loop transport in an external source copy."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('workdir',type=Path);p.add_argument('--radius',type=float,default=.35);p.add_argument('--radial-min',type=float,default=1.);p.add_argument('--radial-max',type=float,default=2.);p.add_argument('--resume',action='store_true');p.add_argument('--seam-pair',action='store_true');p.add_argument('--pilot',action='store_true');p.add_argument('--cache-source',type=Path,required=True);a=p.parse_args()
s=Path(__file__).resolve().parents[3];r=a.workdir.resolve();b=r/'source'
if r==s or s in r.parents:p.error('external directory required')
geometry=dict(radius=a.radius,radial_min=a.radial_min,radial_max=a.radial_max)
assert 0<a.radial_min<1.5-a.radius and (1.5**2+.3**2)**.5+a.radius<a.radial_max,'analytic support must stay inside shell'
if not a.resume:
 r.mkdir(parents=True,exist_ok=False)
 (r/'geometry.json').write_text(json.dumps(geometry,indent=2))
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
else:
 assert json.loads((r/'geometry.json').read_text())==geometry
 for rel,digest in json.loads((r/'manifest.json').read_text()).items():
  if rel.endswith('/run.py'):continue
  assert hashlib.sha256((s/rel).read_bytes()).hexdigest()==digest,rel
 shutil.copy2(__file__,r/'resumed-run.py')
guard=b/'tests/yinyang/interface_transport/convergence/memory_guard.py'
env=dict(os.environ,AMRVAC_DIR=str(b),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PATH='/usr/local/openmpi/bin:'+os.environ['PATH'])
status=json.loads((r/'status.json').read_text()) if a.resume else {};results=json.loads((r/'results.json').read_text()) if a.resume else {}
def job(name,cmd,cwd,limit=1024):
 with (r/(name+'.log')).open('w') as log:
  proc=subprocess.run([sys.executable,str(guard),'--out',str(r/(name+'-memory')),'--limit-mib',str(limit),'--',*map(str,cmd)],cwd=cwd,env=env,stdout=log,stderr=log)
 status[name]=proc.returncode;(r/'status.json').write_text(json.dumps(status,indent=2));print(name,proc.returncode,flush=True)
 if proc.returncode:raise SystemExit(proc.returncode)
compiler='F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096'
par=(s/'tests/yinyang/ct_preflight/amrvac.par').read_text()
par=par.replace('xprobmin1=1.d0',f'xprobmin1={a.radial_min}').replace('xprobmax1=2.d0',f'xprobmax1={a.radial_max}')
par=par.replace('it_max=5','it_max=300').replace('ditsave_dat=5','ditsave_dat=1000')
par=par.replace('&methodlist','&methodlist\n ct_stage_split=T\n yy_ct_experimental=T').replace("20*'minmod'","20*'vanleer'")
par=par.replace('&paramlist','&paramlist\n dtpar=0.002d0')
par=par.replace('&meshlist','&meshlist\n max_blocks=128')
par=par.replace('domain_nx2=16','domain_nx2=24').replace('domain_nx3=48','domain_nx3=72')
par=par.replace('block_nx2=8','block_nx2=12').replace('block_nx3=8','block_nx3=12')
for arch in (['debug'] if a.pilot or a.seam_pair else ['debug','default']):
 d=r/arch;d.mkdir(exist_ok=a.resume);shutil.copy2(b/'tests/yinyang/ct_transport/mod_usr.t',d/'mod_usr.t');(d/'amrvac.par').write_text(par)
 if status.get('setup-'+arch)!=0:job('setup-'+arch,[b/'setup.pl','-d=3','-arch='+arch],d)
 if status.get('build-'+arch)!=0:job('build-'+arch,['make',compiler,'-j1'],d,4096)
 cases=[(1,1,1,2,1)] if a.pilot else [(1,1,1,2,1),(-1,1,1,2,1),(1,0,1,2,1),(-1,0,1,2,1),(1,1,2,2,1),(-1,1,2,2,1),(1,0,2,2,1),(-1,0,2,2,1),(1,1,1,2,0.5),(1,1,1,4,1)]
 if a.seam_pair:cases=[(1,1,1,2,1),(1,1,2,2,1)]
 if arch=='default':cases=[(1,1,1,2,1)]
 for direction,placement,scale,ranks,amp in cases:
  name=f'{arch}-d{direction}-p{placement}-s{scale}-np{ranks}-a{amp}';out=r/name
  if name in results:continue
  out.mkdir(exist_ok=a.resume);shutil.copy2(d/'amrvac',out/'amrvac')
  config=par
  if scale==2:
   config=config.replace('block_nx1=8','block_nx1=16').replace('domain_nx2=24','domain_nx2=32').replace('domain_nx3=72','domain_nx3=96').replace('block_nx2=12','block_nx2=16').replace('block_nx3=12','block_nx3=16')
   prior=json.loads((r/f'{arch}-d1-p1-s1-np2-a1-memory/summary.json').read_text())
   estimate=1.25*max(prior['peak_rss_mib'],prior['peak_footprint_mib'])*(20/16)**2*(20/24)
   (out/'memory-prediction.json').write_text(json.dumps(dict(estimate_mib=estimate,ceiling_mib=950)))
   if estimate>950:raise RuntimeError('refinement memory forecast exceeds 950 MiB')
  config+=f'\n&transport_test\n radius={a.radius}\n direction={direction}\n placement={placement}\n amplitude={amp*1e-3}\n/\n'
  (out/'amrvac.par').write_text(config)
  job(name,['mpirun','--oversubscribe','-np',ranks,'./amrvac'],out)
  mem=r/(name+'-memory');assert (mem/'command.err').stat().st_size==0
  history=[[float(v) for v in line.split()] for line in (out/'transport-history.txt').read_text().splitlines()]
  assert len(history)==301 and history[-1][0]==300 and abs(history[-1][1]-.6)<1e-12
  hashes={}
  for prefix in ['transport','faces']:
   lines=sorted(line for f in out.glob(prefix+'-rank-*.txt') for line in f.read_text().splitlines())
   hashes[prefix]=hashlib.sha256('\n'.join(lines).encode()).hexdigest()
  results[name]=dict(initial=history[0],final=history[-1],hashes=hashes)
  (r/'results.json').write_text(json.dumps(results,indent=2))
  if ranks==4:assert hashes==results[f'{arch}-d{direction}-p{placement}-s{scale}-np2-a{amp}']['hashes']
print('CT transport suite completed',flush=True)

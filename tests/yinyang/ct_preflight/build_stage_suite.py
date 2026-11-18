#!/usr/bin/env python3
"""External, restartable serial build/validation of the staged CT API.

Use the project Python with NumPy. --cache-source is an optional prior external
source+lib tree. Changed templates invalidate generated files in every shape.
Resume only when source/template hashes match the saved manifest.
"""
import argparse,hashlib,importlib.util,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('workdir',type=Path);p.add_argument('--cache-source',type=Path)
p.add_argument('--resume',action='store_true');a=p.parse_args()
s=Path(__file__).resolve().parents[3];r=a.workdir.resolve();b=r/'source'
if r==s or s in r.parents:p.error('external workdir required')
if importlib.util.find_spec('numpy') is None:p.error('NumPy is required')
files=[*s.joinpath('src').rglob('*.t'),*Path(__file__).parent.glob('*')]
manifest={str(f.relative_to(s)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files if f.is_file()}
if a.resume:
 if json.loads((r/'manifest.json').read_text())!=manifest:p.error('source changed; choose a fresh workdir')
else:
 r.mkdir(parents=True,exist_ok=False)
 shutil.copytree(s,b,ignore=shutil.ignore_patterns('.git','lib','__pycache__'))
 (r/'manifest.json').write_text(json.dumps(manifest,indent=2))
 (r/'source.diff').write_bytes(subprocess.check_output(['git','diff','HEAD'],cwd=s))
 (r/'commit.txt').write_bytes(subprocess.check_output(['git','rev-parse','HEAD'],cwd=s))
 if a.cache_source:
  old=a.cache_source.resolve()
  if old==s or s in old.parents:p.error('cache must be external')
  shutil.copytree(old/'lib',b/'lib')
  changed=[]
  for f in b.joinpath('src').rglob('*.t'):
   prior=old/f.relative_to(b)
   if not prior.exists() or f.read_bytes()!=prior.read_bytes():
    changed.append(str(f.relative_to(b)));f.touch()
    for lib in b.joinpath('lib').iterdir():
     if lib.is_dir():
      for suffix in ('.f','.o'):(lib/(f.stem+suffix)).unlink(missing_ok=True)
  (r/'cache-invalidation.json').write_text(json.dumps(dict(cache=str(old),templates=changed),indent=2))
env=dict(os.environ,AMRVAC_DIR=str(b),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PATH='/usr/local/openmpi/bin:'+os.environ['PATH'])
guard=b/'tests/yinyang/interface_transport/convergence/memory_guard.py'
compiler='F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096'
status=json.loads((r/'status.json').read_text()) if (r/'status.json').exists() else {}
def job(name,cmd,cwd,limit=1024):
 if status.get(name,{}).get('passed'):return
 attempt=1
 while (r/f'{name}-attempt{attempt}').exists():attempt+=1
 out=r/f'{name}-attempt{attempt}'
 with (r/f'{name}-attempt{attempt}.log').open('w') as log:
  proc=subprocess.run([sys.executable,str(guard),'--limit-mib',str(limit),'--out',str(out),'--',*map(str,cmd)],cwd=cwd,env=env,stdout=log,stderr=log)
 status[name]=dict(passed=proc.returncode==0,returncode=proc.returncode,evidence=str(out))
 (r/'status.json').write_text(json.dumps(status,indent=2));print(name,proc.returncode,flush=True)
 if proc.returncode:raise SystemExit(proc.returncode)
for arch in ('debug','default'):
 d=r/f'ct-{arch}';d.mkdir(exist_ok=True)
 for name in ('mod_usr.t','amrvac.par'):shutil.copy2(b/'tests/yinyang/ct_preflight'/name,d/name)
 job('setup-'+arch,[b/'setup.pl','-d=3','-arch='+arch],d)
 job('build-'+arch,['make',compiler,'-j1'],d,4096)
 suite=r/('stage-suite-'+arch)
 if suite.exists() and not status.get('suite-'+arch,{}).get('passed'):
  suffix=1
  while suite.with_name(suite.name+f'-failed{suffix}').exists():suffix+=1
  suite.rename(suite.with_name(suite.name+f'-failed{suffix}'))
 job('suite-'+arch,[sys.executable,b/'tests/yinyang/ct_preflight/stage_suite.py',d/'amrvac',suite],r)
for nd,nv in ((1,1),(1,2),(1,3),(2,2),(2,3),(3,3)):
 d=b/'lib'/f'{nd}d{nv}_debug';d.mkdir(parents=True,exist_ok=True)
 if not (d/'amrvac.h').exists():(d/'amrvac.h').touch()
 previous=env.copy();env.update(NDIM=str(nd),NDIR=str(nv),ARCH='debug')
 job(f'matrix-{nd}d{nv}',['make',compiler,'-j1','-f',b/'arch/lib.make','libamrvac'],d,4096)
 env.clear();env.update(previous)
print('Staged CT tests and six-configuration build passed',flush=True)

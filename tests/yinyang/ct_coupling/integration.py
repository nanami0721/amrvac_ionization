#!/usr/bin/env python3
"""Build/run the CT edge and ghost constraint adapters in an external source copy."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('workdir',type=Path);p.add_argument('--cache-source',type=Path,required=True);a=p.parse_args()
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
status={}
def job(name,cmd,cwd,limit=1024):
 with (r/(name+'.log')).open('w') as log:
  proc=subprocess.run([sys.executable,str(guard),'--out',str(r/(name+'-memory')),'--limit-mib',str(limit),'--',*map(str,cmd)],cwd=cwd,env=env,stdout=log,stderr=log)
 status[name]=proc.returncode;(r/'status.json').write_text(json.dumps(status,indent=2));print(name,proc.returncode,flush=True)
 if proc.returncode:raise SystemExit(proc.returncode)
compiler='F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096'
par=(s/'tests/yinyang/rho/amrvac.par').read_text().replace('it_max=5','it_max=0').replace('&meshlist','&meshlist\n max_blocks=128')
for arch in ['debug','default']:
 d=r/arch;d.mkdir();shutil.copy2(b/'tests/yinyang/ct_coupling/mod_usr.t',d/'mod_usr.t');(d/'amrvac.par').write_text(par)
 job('setup-'+arch,[b/'setup.pl','-d=3','-arch='+arch],d)
 job('build-'+arch,['make',compiler,'-j1'],d,4096)
 for stretch in [False,True]:
  for ranks in [2,4,6]:
   name=f'{arch}-stretch{int(stretch)}-np{ranks}';out=r/name;out.mkdir();shutil.copy2(d/'amrvac',out/'amrvac')
   config=par.replace('&meshlist',"&meshlist\n stretch_dim(1)='uni'") if stretch else par
   (out/'amrvac.par').write_text(config)
   job(name,['mpirun','--oversubscribe','-np',ranks,'./amrvac'],out)
   mem=r/(name+'-memory');assert (mem/'command.err').stat().st_size==0
   assert (mem/'command.log').read_text().count('CT_COUPLING_MESH ')==1
for nd,nv in [(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]:
 d=b/'lib'/f'{nd}d{nv}_debug'
 env.update(NDIM=str(nd),NDIR=str(nv),ARCH='debug')
 job(f'matrix-{nd}d{nv}',['make',compiler,'-j1','-f',b/'arch/lib.make','libamrvac'],d,4096)
print('Mesh integration and full-library dimension checks passed',flush=True)

#!/usr/bin/env python3
"""Build/run the same-patch CT edge and ordered magnetic ghost synchronization in an external source copy."""
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
status={};results={}
def job(name,cmd,cwd,limit=1024):
 with (r/(name+'.log')).open('w') as log:
  proc=subprocess.run([sys.executable,str(guard),'--out',str(r/(name+'-memory')),'--limit-mib',str(limit),'--',*map(str,cmd)],cwd=cwd,env=env,stdout=log,stderr=log)
 status[name]=proc.returncode;(r/'status.json').write_text(json.dumps(status,indent=2));print(name,proc.returncode,flush=True)
 if proc.returncode:raise SystemExit(proc.returncode)
compiler='F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096'
par=(s/'tests/yinyang/rho/amrvac.par').read_text().replace('it_max=5','it_max=0').replace('&meshlist','&meshlist\n max_blocks=128\n iprob=0')
for arch in ['debug','default']:
 d=r/arch;d.mkdir();shutil.copy2(b/'tests/yinyang/ct_patch_sync/mod_usr.t',d/'mod_usr.t');(d/'amrvac.par').write_text(par)
 job('setup-'+arch,[b/'setup.pl','-d=3','-arch='+arch],d)
 job('build-'+arch,['make',compiler,'-j1'],d,4096)
 for stretch in [False,True]:
  for layout,b2,b3 in [('split',8,8),('wide',16,16)]:
   for ranks in [2,4,6]:
    name=f'{arch}-stretch{int(stretch)}-{layout}-np{ranks}';out=r/name;out.mkdir();shutil.copy2(d/'amrvac',out/'amrvac')
    config=par.replace('domain_nx1=16','domain_nx1=8')
    config=config.replace('&meshlist','&meshlist\n max_blocks=128').replace('block_nx1=8','block_nx1=4').replace('block_nx2=8',f'block_nx2={b2}').replace('block_nx3=8',f'block_nx3={b3}')
    if stretch:config=config.replace('&meshlist',"&meshlist\n stretch_dim(1)='uni'")
    (out/'amrvac.par').write_text(config)
    job(name,['mpirun','--oversubscribe','-np',ranks,'./amrvac'],out)
    mem=r/(name+'-memory');assert (mem/'command.err').stat().st_size==0
    rows=[x.split() for x in (mem/'command.log').read_text().splitlines() if x.startswith('CT_PATCH_SYNC ')]
    assert len(rows)==1,rows
    x=rows[0];assert int(x[1])==0 and all(int(v)>0 for v in x[2:5]) and max(map(float,x[5:]))<2e-11,rows
    lines=sorted(line for f in out.glob('sync-rank-*.txt') for line in f.read_text().splitlines())
    results[name]=dict(bad=int(x[1]),corners=int(x[2]),holes=int(x[3]),mixed=int(x[4]),edge_error=float(x[5]),divb=float(x[6]),sha256=hashlib.sha256('\n'.join(lines).encode()).hexdigest())
    (r/'results.json').write_text(json.dumps(results,indent=2))
   keys=[f'{arch}-stretch{int(stretch)}-{layout}-np{x}' for x in [2,4,6]]
   assert len({results[k]['sha256'] for k in keys})==1
 # Invalid receiver must fail before any collective exchange, without hanging.
 name='reject-nonowner-'+arch;out=r/name;out.mkdir();shutil.copy2(d/'amrvac',out/'amrvac')
 (out/'amrvac.par').write_text(config.replace('iprob=0','iprob=1'))
 with (r/(name+'.log')).open('w') as log:
  proc=subprocess.run([sys.executable,str(guard),'--out',str(r/(name+'-memory')),'--limit-mib','1024','--timeout-seconds','60','--','mpirun','--oversubscribe','-np','2','./amrvac'],cwd=out,env=env,stdout=log,stderr=log)
 status[name]=proc.returncode;(r/'status.json').write_text(json.dumps(status,indent=2))
 message=(r/(name+'-memory/command.log')).read_text()+(r/(name+'-memory/command.err')).read_text()
 assert proc.returncode!=0 and 'CT overset receiver must be the geometric owner' in message
 assert json.loads((r/(name+'-memory/summary.json')).read_text())['stop_reason'] is None
 print(name,'expected rejection',flush=True)
for nd,nv in [(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]:
 d=b/'lib'/f'{nd}d{nv}_debug'
 env.update(NDIM=str(nd),NDIR=str(nv),ARCH='debug')
 job(f'matrix-{nd}d{nv}',['make',compiler,'-j1','-f',b/'arch/lib.make','libamrvac'],d,4096)
print('Mesh integration and full-library dimension checks passed',flush=True)

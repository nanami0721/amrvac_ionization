#!/usr/bin/env python3
"""Build/run the magnetic ghost fill and divergence constraints in an external source copy."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--identity-only',action='store_true',help='One mesh per layout for exact constant/affine reproduction');p.add_argument('workdir',type=Path);p.add_argument('--cache-source',type=Path,required=True);a=p.parse_args()
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
par=(s/'tests/yinyang/rho/amrvac.par').read_text().replace('it_max=5','it_max=0')
for arch in ['debug','default']:
 d=r/arch;d.mkdir();shutil.copy2(b/'tests/yinyang/ct_ghost_fill/mod_usr.t',d/'mod_usr.t');(d/'amrvac.par').write_text(par)
 job('setup-'+arch,[b/'setup.pl','-d=3','-arch='+arch],d)
 job('build-'+arch,['make',compiler,'-j1'],d,4096)
 for stretch in [False,True]:
  for scale in ([1] if a.identity_only else [1,2,3]):
   for ranks in [2,4,6]:
    name=f'{arch}-stretch{int(stretch)}-s{scale}-np{ranks}';out=r/name;out.mkdir();shutil.copy2(d/'amrvac',out/'amrvac')
    config=par.replace('domain_nx1=16','domain_nx1=8').replace('domain_nx2=16',f'domain_nx2={16*scale}').replace('domain_nx3=48',f'domain_nx3={48*scale}')
    config=config.replace('&meshlist','&meshlist\n max_blocks=128').replace('block_nx1=8','block_nx1=4').replace('block_nx2=8','block_nx2=16').replace('block_nx3=8','block_nx3=16')
    if stretch:config=config.replace('&meshlist',"&meshlist\n stretch_dim(1)='uni'")
    (out/'amrvac.par').write_text(config)
    if scale>1:
     prior=r/f'{arch}-stretch{int(stretch)}-s{1 if scale==2 else 2}-np{ranks}-memory/summary.json'
     measured=json.loads(prior.read_text())['peak_rss_mib']
     estimate=measured*4*1.2
     if scale==3:
      coarse=json.loads((r/f'{arch}-stretch{int(stretch)}-s1-np{ranks}-memory/summary.json').read_text())['peak_rss_mib']
      estimate=1.25*(coarse+max(0,measured-coarse)*8/3)
     if estimate>950:raise RuntimeError(f'{name}: predicted process-tree RSS {estimate:.1f} MiB exceeds safety gate')
    job(name,['mpirun','--oversubscribe','-np',ranks,'./amrvac'],out)
    mem=r/(name+'-memory');assert (mem/'command.err').stat().st_size==0
    rows=[x.split() for x in (mem/'command.log').read_text().splitlines() if x.startswith('MAGNETIC_GHOST ')]
    assert len(rows)==4
    data={f'{x[1]}-{x[2]}':dict(max_error=float(x[-3]),rms_error=float(x[-2]),divb=float(x[-1]),receivers=int(x[3])) for x in rows}
    lines=sorted(line for f in out.glob('ghost-rank-*.txt') for line in f.read_text().splitlines())
    data['sha256']=hashlib.sha256('\n'.join(lines).encode()).hexdigest()
    results[name]=data;(r/'results.json').write_text(json.dumps(results,indent=2))
   keys=[f'{arch}-stretch{int(stretch)}-s{scale}-np{x}' for x in [2,4,6]]
   assert len({results[k]['sha256'] for k in keys})==1
  # Constant and affine solenoidal fluxes are reproduction identities now;
  # convergence slopes of roundoff-level errors would be meaningless.
  for key,data in results.items():
   if not key.startswith(arch+'-'):continue
   for field in [0,1]:
    for closure in [0,1]:assert data[f'{field}-{closure}']['max_error']<1e-10
for nd,nv in [(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]:
 d=b/'lib'/f'{nd}d{nv}_debug'
 env.update(NDIM=str(nd),NDIR=str(nv),ARCH='debug')
 job(f'matrix-{nd}d{nv}',['make',compiler,'-j1','-f',b/'arch/lib.make','libamrvac'],d,4096)
print('Mesh integration and full-library dimension checks passed',flush=True)

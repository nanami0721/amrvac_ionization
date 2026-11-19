#!/usr/bin/env python3
"""Build/run the experimental Yin-Yang CT static-field evolution in an external source copy."""
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
par=(s/'tests/yinyang/ct_preflight/amrvac.par').read_text()
par=par.replace('it_max=5','it_max=50').replace('ditsave_dat=5','ditsave_dat=1000')
par=par.replace('&methodlist','&methodlist\n ct_stage_split=T\n yy_ct_experimental=T')
par=par.replace('&paramlist','&paramlist\n dtpar=0.001d0')
par=par.replace('&meshlist','&meshlist\n max_blocks=128').replace('block_nx1=8','block_nx1=4').replace('domain_nx1=16','domain_nx1=8')
for arch in ['debug','default']:
 d=r/arch;d.mkdir();shutil.copy2(b/'tests/yinyang/ct_static/mod_usr.t',d/'mod_usr.t');(d/'amrvac.par').write_text(par)
 job('setup-'+arch,[b/'setup.pl','-d=3','-arch='+arch],d)
 job('build-'+arch,['make',compiler,'-j1'],d,4096)
 cases=[(-1,1,0,2),(0,1,0,2),(1,1,0,2),(1,1,0,4),(1,2,0,2),(1,2,0,4),(1,1,1,2)]
 if arch=='default':cases.append((2,1,0,2))
 for field,scale,stretch,ranks in cases:
  name=f'{arch}-f{field}-s{scale}-r{stretch}-np{ranks}';out=r/name;out.mkdir();shutil.copy2(d/'amrvac',out/'amrvac')
  config=par.replace('&meshlist',f'&meshlist\n iprob={field}').replace('domain_nx2=16',f'domain_nx2={16*scale}').replace('domain_nx3=48',f'domain_nx3={48*scale}')
  if scale==2:
   config=config.replace('block_nx2=8','block_nx2=16').replace('block_nx3=8','block_nx3=16')
   prior=json.loads((r/f'{arch}-f1-s1-r0-np{ranks}-memory/summary.json').read_text())
   # Same block count/radial extent; include two ghosts on each angular side.
   growth=((16+4)/(8+4))**2
   estimate=1.25*max(prior['peak_rss_mib'],prior['peak_footprint_mib'])*growth
   if ranks==4:
    fine=json.loads((r/f'{arch}-f1-s2-r0-np2-memory/summary.json').read_text())
    coarse=json.loads((r/f'{arch}-f1-s1-r0-np2-memory/summary.json').read_text())
    peak=lambda v:max(v['peak_rss_mib'],v['peak_footprint_mib'])
    # Fixed global block count: take measured fine allocation, plus twice the
    # observed coarse 2-to-4-rank overhead and 25% margin.
    estimate=1.25*(peak(fine)+2*max(0,peak(prior)-peak(coarse)))
   (out/'memory-prediction.json').write_text(json.dumps(dict(estimate_mib=estimate,ceiling_mib=950)))
   if estimate>950:raise RuntimeError('refinement memory forecast exceeds 950 MiB')
  steps=500 if field==2 else 50
  if field==2:config=config.replace('it_max=50','it_max=500').replace('iprob=2','iprob=1').replace('&meshlist','&meshlist\n lb_diagnose=T')
  if field==-1:config=config.replace('iprob=-1','iprob=0').replace('yy_ct_experimental=T','yy_ct_experimental=F')+'\n&static_test\n spherical_control=T\n/\n'
  if stretch:config=config.replace('&meshlist',"&meshlist\n stretch_dim(1)='uni'")
  (out/'amrvac.par').write_text(config)
  job(name,['mpirun','--oversubscribe','-np',ranks,'./amrvac'],out)
  mem=r/(name+'-memory');assert (mem/'command.err').stat().st_size==0
  history=[[float(v) for v in line.split()] for line in (out/'static-history.txt').read_text().splitlines()]
  assert len(history)==steps+1 and history[-1][0]==steps and abs(history[-1][1]-steps*.001)<1e-12
  lines=sorted(line for f in out.glob('static-rank-*.txt') for line in f.read_text().splitlines())
  results[name]=dict(initial=history[0],final=history[-1],sha256=hashlib.sha256('\n'.join(lines).encode()).hexdigest())
  (r/'results.json').write_text(json.dumps(results,indent=2))
  if ranks==4:assert results[name]['sha256']==results[f'{arch}-f{field}-s{scale}-r{stretch}-np2']['sha256']
 for label,replacement,expected in [('disabled',('yy_ct_experimental=T','yy_ct_experimental=F'),'requires explicit yy_ct_experimental'),('unstaged',('ct_stage_split=T','ct_stage_split=F'),'YY CT needs staged'),('contact',("type_ct='uct_hll'","type_ct='uct_contact'"),'YY CT needs staged')]:
  name=f'reject-{arch}-{label}';out=r/name;out.mkdir();shutil.copy2(d/'amrvac',out/'amrvac')
  (out/'amrvac.par').write_text(par.replace(*replacement))
  with (r/(name+'.log')).open('w') as log:
   proc=subprocess.run([sys.executable,str(guard),'--out',str(r/(name+'-memory')),'--limit-mib','1024','--timeout-seconds','60','--','mpirun','--oversubscribe','-np','2','./amrvac'],cwd=out,env=env,stdout=log,stderr=log)
  message=(r/(name+'-memory/command.log')).read_text()+(r/(name+'-memory/command.err')).read_text()
  status[name]=proc.returncode;(r/'status.json').write_text(json.dumps(status,indent=2))
  assert proc.returncode!=0 and expected in message
  assert json.loads((r/(name+'-memory/summary.json')).read_text())['stop_reason'] is None
  print(name,'expected rejection',flush=True)
for nd,nv in [(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]:
 d=b/'lib'/f'{nd}d{nv}_debug'
 env.update(NDIM=str(nd),NDIR=str(nv),ARCH='debug')
 job(f'matrix-{nd}d{nv}',['make',compiler,'-j1','-f',b/'arch/lib.make','libamrvac'],d,4096)
print('Static MHD CT suite and full-library dimension checks passed',flush=True)

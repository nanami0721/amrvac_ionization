#!/usr/bin/env python3
"""Build/run the experimental Yin-Yang CT weak-loop transport in an external source copy."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description='Same-orientation native controls and one-factor CT accuracy probes')
p.add_argument('workdir',type=Path)
p.add_argument('--resume',action='store_true')
p.add_argument('--pilot',action='store_true')
p.add_argument('--cache-source',type=Path,required=True,help='validated YY wide-transport source; parent must retain the fine two-rank memory record')
a=p.parse_args()
s=Path(__file__).resolve().parents[3];r=a.workdir.resolve();b=r/'source'
if r==s or s in r.parents:p.error('external directory required')
geometry=dict(radius=.65,radial_min=.75,radial_max=2.25)
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
  if rel==str(Path(__file__).resolve().relative_to(s)):continue
  assert hashlib.sha256((s/rel).read_bytes()).hexdigest()==digest,rel
 shutil.copy2(__file__,r/'resumed-run_controls.py')
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
par=par.replace('xprobmin1=1.d0','xprobmin1=.75d0').replace('xprobmax1=2.d0','xprobmax1=2.25d0')
par=par.replace('it_max=5','it_max=300').replace('ditsave_dat=5','ditsave_dat=1000')
par=par.replace('&methodlist','&methodlist\n ct_stage_split=T\n yy_ct_experimental=T').replace("20*'minmod'","20*'vanleer'")
par=par.replace('&paramlist','&paramlist\n dtpar=.002d0').replace('&meshlist','&meshlist\n max_blocks=128')
par=par.replace('domain_nx2=16','domain_nx2=24').replace('domain_nx3=48','domain_nx3=72').replace('block_nx2=8','block_nx2=12').replace('block_nx3=8','block_nx3=12')
d=r/'debug';d.mkdir(exist_ok=a.resume);shutil.copy2(b/'tests/yinyang/ct_transport/mod_usr.t',d/'mod_usr.t');(d/'amrvac.par').write_text(par)
if status.get('setup-debug')!=0:job('setup-debug',[b/'setup.pl','-d=3','-arch=debug'],d)
if status.get('build-debug')!=0:job('build-debug',['make',compiler,'-j1'],d,4096)
# native orientations, baseline YY samples, and one-factor accuracy probes.
cases=[(True,0,1,1,'base'),(True,1,1,1,'base'),(True,0,2,1,'base'),(True,1,2,1,'base'),(True,0,2,-1,'base'),(True,1,2,-1,'base'),(False,0,2,1,'base'),(False,0,2,-1,'base'),(True,1,2,1,'halfdt'),(True,1,2,1,'mc'),(False,0,2,1,'mc'),(True,1,2,1,'radial'),(True,1,2,1,'hll'),(False,0,2,1,'hll'),(True,1,2,1,'hllmc'),(False,0,2,1,'hllmc'),(False,0,2,1,'halfdt')]
if a.pilot:cases=cases[:1]
for native,patch,scale,direction,variant in cases:
 name=f'{"native" if native else "yy"}-p{patch}-s{scale}-d{direction}-{variant}';out=r/name
 if name in results:continue
 out.mkdir(exist_ok=a.resume);shutil.copy2(d/'amrvac',out/'amrvac');config=par
 if scale==2:config=config.replace('block_nx1=8','block_nx1=16').replace('domain_nx2=24','domain_nx2=32').replace('domain_nx3=72','domain_nx3=96').replace('block_nx2=12','block_nx2=16').replace('block_nx3=12','block_nx3=16')
 if native:
  for side in ['min2','max2','min3','max3']:config=config.replace(f"typeboundary_{side}=8*'cont'",f"typeboundary_{side}=5*'special',3*'cont'")
  config=config.replace('yy_ct_experimental=T','yy_ct_experimental=F')
  if patch==0:
   config=config.replace('xprobmin3=-0.390625d0','xprobmin3=0.130208333333333333d0').replace('xprobmax3=0.390625d0','xprobmax3=0.651041666666666667d0').replace(f'domain_nx3={72 if scale==1 else 96}',f'domain_nx3={48 if scale==1 else 64}')
  else:
   config=config.replace('xprobmin2=0.109375d0','xprobmin2=0.0390625d0').replace('xprobmax2=0.390625d0','xprobmax2=0.3203125d0').replace('xprobmin3=-0.390625d0','xprobmin3=-0.130208333333333333d0').replace('xprobmax3=0.390625d0','xprobmax3=0.130208333333333333d0').replace(f'domain_nx3={72 if scale==1 else 96}',f'domain_nx3={24 if scale==1 else 32}')
 steps=300
 if variant=='halfdt':config=config.replace('dtpar=.002d0','dtpar=.001d0').replace('it_max=300','it_max=600');steps=600
 if variant in ['mc','hllmc']:config=config.replace("20*'vanleer'","20*'woodward'")
 if variant=='radial':config=config.replace('domain_nx1=16','domain_nx1=32')
 if variant in ['hll','hllmc']:config=config.replace("20*'tvdlf'","20*'hll'")
 config+=f'\n&transport_test\n radius=.65\n direction={direction}\n placement=1\n amplitude=.001\n native_control={"T" if native else "F"}\n reference_patch={patch}\n save_samples=T\n/\n'
 (out/'amrvac.par').write_text(config)
 # All native allocations fit below the previous measured YY fine layout.
 if not native:
  prior=json.loads((a.cache_source.parent/'debug-d1-p1-s2-np2-a1-memory/summary.json').read_text())
  estimate=1.25*max(prior['peak_rss_mib'],prior['peak_footprint_mib'])
  (out/'memory-prediction.json').write_text(json.dumps(dict(estimate_mib=estimate,ceiling_mib=950)))
  if estimate>950:raise RuntimeError('YY memory prediction exceeds ceiling')
 job(name,['mpirun','--oversubscribe','-np','2','./amrvac'],out)
 assert (r/(name+'-memory/command.err')).stat().st_size==0
 h=[[float(v) for v in line.split()] for line in (out/'transport-history.txt').read_text().splitlines()]
 assert len(h)==steps+1 and h[-1][0]==steps and abs(h[-1][1]-.6)<1e-12
 results[name]=dict(initial=h[0],final=h[-1]);(r/'results.json').write_text(json.dumps(results,indent=2))
print('Native control and accuracy probes completed',flush=True)

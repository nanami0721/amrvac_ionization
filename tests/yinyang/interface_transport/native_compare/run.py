#!/usr/bin/env python3
"""External-only convergence experiments with explicit immutable parameters."""
import argparse, hashlib, json, shutil, subprocess, time, os, sys
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('root',type=Path);p.add_argument('exe',type=Path)
p.add_argument('--grid',choices=['spherical'],default='spherical')
p.add_argument('--scale',type=int,required=True)
p.add_argument('--courant',type=float,default=.3)
p.add_argument('--width',type=float,default=.18)
p.add_argument('--limiter',choices=['woodward','cada3'],default='woodward')
p.add_argument('--extended',action='store_true')
p.add_argument('--no-divbfix',action='store_true')
p.add_argument('--debug',action='store_true')
p.add_argument('--no-native',action='store_true')
p.add_argument('--ranks',type=int,choices=[1,2,4],default=1)
p.add_argument('--direction',type=int,choices=[-1,1],default=1)
a=p.parse_args();r=a.root.resolve();source=Path(__file__).resolve().parents[4]
if r==source or source in r.parents:p.error('external directory required')
if a.scale<1 or not 0<a.courant<=.3 or a.width<=0:p.error('invalid numerical parameter')
if a.grid=='cartesian' and a.extended:p.error('extended is spherical only')
# Direct invocation also uses the watchdog; a wrapped child inherits this marker.
if os.environ.get('AMRVAC_MEMORY_GUARD_ACTIVE')!='1':
    guard=Path(__file__).resolve().parent.parent/'convergence/memory_guard.py'
    result=subprocess.run([sys.executable,str(guard),'--out',str(r.with_name(r.name+'-memory')),'--',sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]])
    raise SystemExit(result.returncode)
r.mkdir(parents=True,exist_ok=False);shutil.copy2(a.exe,r/'amrvac')
nd=3 if a.grid=='spherical' else 1;s=a.scale
bc='\n'.join(f" typeboundary_{side}{axis}=8*'special'" for axis in range(1,nd+1) for side in ['min','max'])
mesh=(f''' block_nx1=8
 block_nx2=4
 block_nx3=4
 domain_nx1={(24 if a.extended else 16)*s}
 domain_nx2={20*s}
 domain_nx3={52*s}
 xprobmin1={.75 if a.extended else 1.}
 xprobmax1={2.25 if a.extended else 2.}
 xprobmin2=.09375d0
 xprobmax2=.40625d0
 xprobmin3=0.d0
 xprobmax3=.8125d0''' if nd==3 else f''' block_nx1=16
 domain_nx1={32*s}
 xprobmin1=-2.356194490192345d0
 xprobmax1=2.356194490192345d0''')
fix='2*.true.,4*.false.' if nd==3 and not a.no_divbfix else f'{2*nd}*.false.'
par=f'''&filelist
 base_filename='adv'
 autoconvert=.false.
/
&savelist
 itsave(1,1)=0
 itsave(1,2)=0
 ditsave_log=20
 dtsave_dat=.15d0
/
&stoplist
 it_max={3 if a.debug else 1000000}
 time_max=.6d0
/
&methodlist
 time_stepper='threestep'
 flux_scheme=20*'tvdlf'
 limiter=20*'{a.limiter}'
/
&boundlist
{bc}
/
&meshlist
 refine_max_level=1
{mesh}
/
&paramlist
 courantpar={a.courant}
/
&mhd_list
 typedivbfix='glm'
 mhd_glm_alpha=.5d0
 boundary_divbfix={fix}
/
&transport_test
 direction={a.direction}
 pulse_width={a.width}
 native_diagnostics={'.false.' if a.no_native else '.true.'}
 native_every_step={'.true.' if a.debug else '.false.'}
{(' test_case=1' + chr(10) + ' placement=1') if nd==3 else ''}
/
'''
(r/'amrvac.par').write_text(par)
meta={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
meta.update(memory_limit_mib=float(os.environ['AMRVAC_MEMORY_LIMIT_MIB']),executable_sha256=hashlib.sha256((r/'amrvac').read_bytes()).hexdigest(),ranks=a.ranks if a.ranks is not None else (4 if nd==3 else 1))
(r/'parameters.json').write_text(json.dumps(meta,indent=2))
t=time.monotonic()
with (r/'run.log').open('w') as out,(r/'run.err').open('w') as err:
    try:
        proc=subprocess.run(['/usr/local/openmpi/bin/mpirun','--oversubscribe','-np',str(meta['ranks']),'./amrvac'],cwd=r,stdout=out,stderr=err,timeout=3600)
        code=proc.returncode
    except subprocess.TimeoutExpired:code=124
meta.update(returncode=code,wall=time.monotonic()-t,stderr=(r/'run.err').stat().st_size)
(r/'run.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
if code or meta['stderr']:raise SystemExit('run failed; inspect external logs')

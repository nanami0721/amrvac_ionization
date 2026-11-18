#!/usr/bin/env python3
"""Run exact entropy / force-free magnetic transport through and away from seams."""
from pathlib import Path
import argparse,subprocess,shutil,json,time,hashlib,os,sys
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('exe',type=Path);p.add_argument('--scale',type=int,choices=[1,2,3],default=1);p.add_argument('--control',action='store_true');p.add_argument('--debug',action='store_true');p.add_argument('--radial-divbfix',action='store_true');p.add_argument('--radial-bc',type=int,choices=[0,1],default=0);p.add_argument('--no-native',action='store_true');p.add_argument('--courant',type=float,default=0.3)
p.add_argument('--case',type=int,choices=[0,1],default=None)
p.add_argument('--limiter',choices=['minmod','woodward'],default='minmod')
p.add_argument('--direction',type=int,choices=[-1,1],default=None)
p.add_argument('--glm-alpha',type=float,default=0.5)
p.add_argument('--extended-shell',action='store_true',help='r=0.75..2.25 at unchanged radial spacing')
a=p.parse_args()
if a.radial_divbfix and a.radial_bc!=0:p.error('radial divB fix and characteristic Bn prescription are alternatives')
if not 0<a.courant<=0.3:p.error('diagnostic courant must be in (0,0.3]')
if not 0<=a.glm_alpha<=1:p.error('diagnostic GLM alpha must be in [0,1]')
a.root=a.root.resolve();src=Path(__file__).resolve().parents[3]
if a.root==src or src in a.root.parents:p.error('external run root required')
if os.environ.get('AMRVAC_MEMORY_GUARD_ACTIVE')!='1':
 guard=Path(__file__).resolve().parent/'convergence/memory_guard.py'
 result=subprocess.run([sys.executable,str(guard),'--out',str(a.root.with_name(a.root.name+'-memory')),'--',sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]])
 raise SystemExit(result.returncode)
a.root.mkdir(parents=True,exist_ok=True)
# Core changes are not made by this runner. Parameters mirror the earlier benchmark.
# Keep the actual parameter source explicit instead of evaluating another runner.
for case in ([a.case] if a.case is not None else [0,1]):
 for direction in ([a.direction] if a.direction is not None else [1,-1]):
  place=0 if a.control else 1
  name=f'c{case}-p{place}-d{direction}-s{a.scale}';d=a.root/name;d.mkdir()
  shutil.copy2(a.exe,d/'amrvac');s=a.scale
  par=f'''&filelist
 base_filename='adv'
 autoconvert=.false.
/
&savelist
 itsave(1,1)=0
 itsave(1,2)=0
 ditsave_log=20
 dtsave_dat=0.3d0
/
&stoplist
 it_max={3 if a.debug else 100000}
 time_max=0.6d0
/
&methodlist
 time_stepper='threestep'
 flux_scheme=20*'tvdlf'
 limiter=20*'{a.limiter}'
/
&boundlist
 typeboundary_min1=8*'special'
 typeboundary_max1=8*'special'
 typeboundary_min2=8*'cont'
 typeboundary_max2=8*'cont'
 typeboundary_min3=8*'cont'
 typeboundary_max3=8*'cont'
/
&meshlist
 refine_max_level=1
 block_nx1=8
 block_nx2=4
 block_nx3=4
 domain_nx1={(24 if a.extended_shell else 16)*s}
 domain_nx2={20*s}
 domain_nx3={52*s}
 xprobmin1={0.75 if a.extended_shell else 1.0}
 xprobmax1={2.25 if a.extended_shell else 2.0}
 xprobmin2=0.09375d0
 xprobmax2=0.40625d0
 xprobmin3=-0.40625d0
 xprobmax3=0.40625d0
/
&paramlist
 courantpar={a.courant}
/
&mhd_list
 typedivbfix='glm'
 mhd_glm_alpha={a.glm_alpha}
 boundary_divbfix={"2*.true.,4*.false." if a.radial_divbfix else "6*.false."}
/
&transport_test
 test_case={case}
 direction={direction}
 placement={place}
 radial_bc={a.radial_bc}
 native_diagnostics={'.false.' if a.no_native else '.true.'}
 native_every_step={'.true.' if a.debug else '.false.'}
/
'''
  (d/'amrvac.par').write_text(par);t=time.monotonic()
  with (d/'run.log').open('w') as out,(d/'run.err').open('w') as err:
   pr=subprocess.run(['/usr/local/openmpi/bin/mpirun','--oversubscribe','-np','4','./amrvac'],cwd=d,stdout=out,stderr=err,timeout=1800)
  result=dict(executable_sha256=hashlib.sha256((d/'amrvac').read_bytes()).hexdigest(),case=case,direction=direction,placement=place,scale=s,debug=a.debug,radial_divbfix=a.radial_divbfix,radial_bc=a.radial_bc,native_diagnostics=not a.no_native,courant=a.courant,limiter=a.limiter,glm_alpha=a.glm_alpha,extended_shell=a.extended_shell,returncode=pr.returncode,wall=time.monotonic()-t,stderr=(d/'run.err').stat().st_size)
  (d/'run.json').write_text(json.dumps(result,indent=2));print(name,json.dumps(result),flush=True)
  if pr.returncode or result['stderr']:raise SystemExit('run failed; inspect stored logs')

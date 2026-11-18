#!/usr/bin/env python3
"""Guarded hydrostatic/dipole equilibrium runs; external directories only."""
from pathlib import Path
import argparse,hashlib,json,math,os,shutil,subprocess,sys,time
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('exe',type=Path)
p.add_argument('--case',type=int,choices=[0,1,2,3],default=2)
p.add_argument('--grid',choices=['yinyang','spherical'],default='yinyang')
p.add_argument('--ranks',type=int,choices=[1,2,4],default=None)
p.add_argument('--crossings',type=float,default=1)
p.add_argument('--cfl',type=float,default=.3);p.add_argument('--beta',type=float,default=4.);p.add_argument('--output-interval',type=float,default=None)
p.add_argument('--refine-direction',choices=['theta','phi']);p.add_argument('--radial-refine',action='store_true');p.add_argument('--debug',action='store_true')
p.add_argument('--no-diagnostics',action='store_true');p.add_argument('--restart',type=Path)
a=p.parse_args();root=a.root.resolve();source=Path(__file__).resolve().parents[3]
if root==source or source in root.parents:p.error('external run root required')
if not 0<a.cfl<=.3:p.error('0<cfl<=0.3 required')
if not math.isfinite(a.beta) or a.beta<=0 or not 0<a.crossings<=10:p.error('finite positive beta and 0<crossings<=10 required')
interval=a.output_interval if a.output_interval is not None else min(1.,a.crossings)
if not 0<interval<=a.crossings or abs(a.crossings/interval-round(a.crossings/interval))>1e-8:p.error('output interval must divide duration')
a.output_interval=interval
ranks=a.ranks if a.ranks is not None else (2 if a.grid=='yinyang' else 1)
if a.grid=='yinyang' and ranks<2:p.error('Yin-Yang needs two patch groups')
if os.environ.get('AMRVAC_MEMORY_GUARD_ACTIVE')!='1':
 guard=Path(__file__).resolve().parents[1]/'interface_transport/convergence/memory_guard.py'
 proc=subprocess.run([sys.executable,str(guard),'--out',str(root.with_name(root.name+'-memory')),'--',sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]])
 raise SystemExit(proc.returncode)
root.mkdir(parents=True,exist_ok=False);shutil.copy2(a.exe,root/'amrvac');tau=1/math.sqrt(5/3*.5)
bc='\n'.join(f" typeboundary_{side}{axis}=8*'{('special' if axis==1 or a.grid=='spherical' else 'cont')}'" for axis in [1,2,3] for side in ['min','max'])
par=f'''&filelist
 base_filename='eq'
 autoconvert=.false.
/
&savelist
 itsave(1,1)=0
 itsave(1,2)=0
 ditsave_log=20
 dtsave_dat={interval*tau:.17g}
/
&stoplist
 it_max={3 if a.debug else 1000000}
 time_max={a.crossings*tau:.17g}
/
&methodlist
 time_stepper='threestep'
 time_integrator='ssprk3'
 ssprk_order=3
 flux_scheme=20*'tvdlf'
 limiter=20*'woodward'
 small_values_method='error'
 check_small_values=.true.
 fix_small_values=.false.
/
&boundlist
{bc}
/
&meshlist
 refine_max_level=1
 block_nx1=8
 block_nx2=4
 block_nx3=4
 domain_nx1={32 if a.radial_refine else 16}
 domain_nx2={40 if a.refine_direction=='theta' else 20}
 domain_nx3={104 if a.refine_direction=='phi' else 52}
 xprobmin1=1.d0
 xprobmax1=2.d0
 xprobmin2=.09375d0
 xprobmax2=.40625d0
 xprobmin3={'-.40625d0' if a.grid=='yinyang' else '0.d0'}
 xprobmax3={'.40625d0' if a.grid=='yinyang' else '.8125d0'}
/
&paramlist
 courantpar={a.cfl:.17g}
/
&mhd_list
 typedivbfix='glm'
 mhd_glm_alpha=.5d0
 mhd_gravity={'.true.' if a.case in [0,2] else '.false.'}
 boundary_divbfix=6*.false.
/
&eos_list
 gamma=1.6666666666666667d0
/
&grav_list
 grav_split=.false.
/
&equilibrium_test
 equilibrium_case={a.case}
 beta_reference={a.beta:.17g}
 output_interval={interval:.17g}
 spherical_control={'.true.' if a.grid=='spherical' else '.false.'}
 diagnostics={'.false.' if a.no_diagnostics else '.true.'}
 diagnostic_every_step={'.true.' if a.debug else '.false.'}
/
'''
(root/'amrvac.par').write_text(par)
cmd=['/usr/local/openmpi/bin/mpirun','--oversubscribe','-np',str(ranks),'./amrvac']
if a.restart:
 prefix='eq_yin' if a.grid=='yinyang' else 'eq'
 snapshots=sorted(a.restart.glob(prefix+'[0-9][0-9][0-9][0-9].dat'))
 if not snapshots:raise SystemExit('restart directory has no snapshots')
 number=snapshots[-1].stem[-4:]
 for f in a.restart.glob(f'eq*{number}.dat*'):shutil.copy2(f,root/f.name)
 cmd+=['-if',f'eq_{{patch}}{number}.dat' if a.grid=='yinyang' else f'eq{number}.dat']
meta={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()};meta.update(ranks=ranks,tau=tau,executable_sha256=hashlib.sha256((root/'amrvac').read_bytes()).hexdigest())
(root/'parameters.json').write_text(json.dumps(meta,indent=2));start=time.monotonic()
with (root/'run.log').open('w') as out,(root/'run.err').open('w') as err:
 try:code=subprocess.run(cmd,cwd=root,stdout=out,stderr=err,timeout=3600).returncode
 except subprocess.TimeoutExpired:code=124
meta.update(returncode=code,wall=time.monotonic()-start,stderr=(root/'run.err').stat().st_size)
(root/'run.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
if code or meta['stderr']:raise SystemExit('run failed; inspect bounded logs')

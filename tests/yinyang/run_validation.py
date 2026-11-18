#!/usr/bin/env python3
"""Build and run outside the source repository; requires MPI, make and NumPy.

python run_validation.py --workdir /private/tmp/amrvac-tests/yy-UNIQUE
Use --extended for convergence and perturbation tests, --compile-matrix for
full dependency compilation of all six dimension/component configurations.
Logs and copied sources remain in workdir. No commit/push is performed.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import hashlib
import importlib.util

sys.dont_write_bytecode = True

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--workdir',type=Path,required=True)
ap.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[2])
ap.add_argument('--mpiexec',default='mpirun')
ap.add_argument('--build-cache',type=Path,help='external lib cache copied by value; all copied objects are rebuilt')
ap.add_argument('--extended',action='store_true')
ap.add_argument('--compile-matrix',action='store_true')
ap.add_argument('--arch',default='debug')
a=ap.parse_args();src=a.source.resolve();root=a.workdir.resolve()
if root==src or src in root.parents:ap.error('workdir must be outside source')
if importlib.util.find_spec('numpy') is None:ap.error('NumPy required; use the project scientific Python environment')
root.mkdir(parents=True,exist_ok=False)
build=root/'build-source'
shutil.copytree(src,build,ignore=shutil.ignore_patterns('.git','lib','__pycache__'))
env=dict(os.environ,AMRVAC_DIR=str(build),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');results={}
guard=src/'tests/yinyang/interface_transport/convergence/memory_guard.py'
worker=src/'tests/yinyang/validation_worker.py'
if a.build_cache:
    cache=a.build_cache.resolve()
    if cache==src or src in cache.parents:ap.error('cache must be external')
    shutil.copytree(cache,build/'lib')
    # Regenerate objects against this exact source, never trust cache identity.
    for item in (build/'lib').rglob('*.o'):item.unlink()
    for item in (build/'lib').rglob('*.a'):item.unlink()
    for item in (build/'lib').rglob('*.f'):item.unlink()
manifest={str(p.relative_to(build)):hashlib.sha256(p.read_bytes()).hexdigest()
          for p in build.rglob('*') if p.is_file() and 'lib' not in p.relative_to(build).parts}
(root/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(root/'source-commit.txt').write_text(subprocess.check_output(['git','rev-parse','HEAD'],cwd=src,text=True))
(root/'source.diff').write_text(subprocess.check_output(['git','diff','HEAD'],cwd=src,text=True))
# The copied source also preserves untracked implementation/test files.

def command(cmd,cwd,stem,environment=env,limit=1024):
    outdir=cwd/(stem+'-memory')
    p=subprocess.run([sys.executable,str(guard),'--limit-mib',str(limit),
                      '--out',str(outdir),'--',*map(str,cmd)],cwd=cwd,env=environment)
    if p.returncode:raise SystemExit(p.returncode)
    return (outdir/'command.err').stat().st_size

def analysis(name,operation,**kwargs):
    request=root/(name+'.request.json');output=root/(name+'.json')
    request.write_text(json.dumps(dict(operation=operation,source=str(build),output=str(output),**kwargs)))
    command([sys.executable,str(worker),str(request)],root,name)
    return json.loads(output.read_text())

compiler='F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096'

if a.compile_matrix:
    for nd,nv in [(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]:
        d=build/'lib'/f'{nd}d{nv}_{a.arch}';d.mkdir(parents=True,exist_ok=True);(d/'amrvac.h').touch()
        command(['make',compiler,'-j1','-f',str(build/'arch/lib.make'),'libamrvac'],d,'compile',dict(env,NDIM=str(nd),NDIR=str(nv),ARCH=a.arch),limit=4096)
        print(f'compiled {nd}D{nv}V',flush=True)

case=root/'build-case';case.mkdir()
for name in ['mod_usr.t','amrvac.par']:shutil.copy2(src/'tests/yinyang/mhd'/name,case/name)
command([str(build/'setup.pl'),'-d=3','-arch='+a.arch],case,'setup')
command(['make',compiler,'-j1'],case,'build',limit=4096)
par=(case/'amrvac.par').read_text()

def run(name,np_,parameters=par,restart=None):
    d=root/name;d.mkdir();shutil.copy2(case/'amrvac',d/'amrvac');(d/'amrvac.par').write_text(parameters)
    cmd=[a.mpiexec,'--oversubscribe','-np',str(np_),'./amrvac']
    if restart:
        for f in restart.glob('*.dat*'):shutil.copy2(f,d/f.name)
        cmd+=['-if','yy_{patch}0001.dat']
    size=command(cmd,d,'run')
    if size:raise AssertionError(f'nonempty runtime stderr: {d}/run.err ({size} bytes)')
    print(name+' passed',flush=True)
    return d

runs={n:run(f'np{n}',n) for n in (2,4)}
for n in (4,):
    results[f'np{n}_difference']=analysis(f'compare-np{n}','compare',left=str(runs[2]),right=str(runs[n]),snapshot=1)['max_difference']
continuous=run('continuous',2,par.replace('it_max=5','it_max=10'))
restart=run('restart',2,par.replace('it_max=5','it_max=10'),runs[2])
results['restart_difference']=analysis('compare-restart','compare',left=str(continuous),right=str(restart),snapshot=2)['max_difference']
results['uniform']=analysis('uniform-analysis','diagnose',run=str(runs[2]),snapshot=1)
results['overlap']=analysis('overlap-analysis','overlap',run=str(runs[2]),snapshot=1)
if a.extended:
    convergence=[]
    for n in (16,24,32):
        p=par.replace('it_max=5','it_max=10000\n time_max=0.026').replace('ditsave_dat=5','ditsave_dat=10000')
        for axis,old,new in [(1,16,n),(2,16,n),(3,48,3*n)]:p=p.replace(f'domain_nx{axis}={old}',f'domain_nx{axis}={new}')
        # Reject the next grid before allocation when the measured estimate is too large.
        baseline=json.loads((runs[4]/'run-memory/summary.json').read_text())
        peak=max(baseline['peak_rss_mib'],baseline['peak_footprint_mib'])
        estimate=peak*(n/16)**3
        if estimate>0.75*1024:
            raise RuntimeError(f'convergence grid {n} exceeds memory preflight: {estimate:.1f} MiB')
        d=run(f'conv{n}',4,p);convergence.append(analysis(f'conv{n}-analysis','diagnose',run=str(d),snapshot=1))
    assert all(y['uniform_error_l1']<x['uniform_error_l1']*.75 for x,y in zip(convergence,convergence[1:]))
    results['convergence']=convergence
    for mode in (1,2):
        d=run(f'case{mode}',4,par.replace('it_max=5','it_max=50').replace('ditsave_dat=5','ditsave_dat=50')+f'\n&yy_test\n initial_case={mode}\n/\n')
        analysis(f'case{mode}-analysis','positive',run=str(d),snapshot=1)
results['arch']=a.arch
(root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(root/'results.json')

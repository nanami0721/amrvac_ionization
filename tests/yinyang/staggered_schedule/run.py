#!/usr/bin/env python3
"""Build actual geometry-only module and test MPI tables outside the repository."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('workdir',type=Path);a=p.parse_args()
s=Path(__file__).resolve().parents[3];r=a.workdir.resolve()
if r==s or s in r.parents:p.error('external directory required')
r.mkdir(parents=True,exist_ok=False);b=r/'source';b.mkdir()
files=[s/'src/mod_yinyang_staggered.t',s/'src/vacpp.pl',Path(__file__).with_name('check.f90'),Path(__file__).resolve()]
for f in files:shutil.copy2(f,b/f.name)
(r/'manifest.json').write_text(json.dumps({str(f.relative_to(s)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},indent=2))
(r/'source.diff').write_bytes(subprocess.check_output(['git','diff','HEAD'],cwd=s))
(r/'commit.txt').write_bytes(subprocess.check_output(['git','rev-parse','HEAD'],cwd=s))
guard=s/'tests/yinyang/interface_transport/convergence/memory_guard.py'
shutil.copy2(guard,b/guard.name);guard=b/guard.name
env=dict(os.environ,AMRVAC_DIR=str(b),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PATH='/usr/local/openmpi/bin:'+os.environ['PATH'])
def job(name,cmd,cwd,limit):
 with (r/(name+'.log')).open('w') as log:
  proc=subprocess.run([sys.executable,str(guard),'--out',str(r/(name+'-memory')),'--limit-mib',str(limit),'--',*map(str,cmd)],cwd=cwd,env=env,stdout=log,stderr=log)
 if proc.returncode:raise RuntimeError(f'{name} failed: {proc.returncode}; see {r/(name+"-memory")}')
 print(name,'passed',flush=True)
results={}
for arch in ['debug','optimized']:
 d=r/arch;d.mkdir()
 # Expand the exact template through VACPP, no flattened-index substitute.
 with (d/'mod_yinyang_staggered.f90').open('w') as out:
  subprocess.run(['perl',str(b/'vacpp.pl'),'-d=3','-v=3',str(b/'mod_yinyang_staggered.t')],stdout=out,check=True)
 flags=['-O2','-g','-fcheck=all','-fbacktrace','-ffpe-trap=invalid,zero,overflow','-finit-real=snan'] if arch=='debug' else ['-O3']
 job('build-'+arch,['mpif90','-ffree-form','-ffree-line-length-none',*flags,d/'mod_yinyang_staggered.f90',b/'check.f90','-o',d/'check'],d,4096)
 for family in [1,2]:
  for mode in [0,1]:
   errors=[]
   for scale in [1,2,4]:
    hashes=[]
    for ranks in [2,4,6]:
     name=f'{arch}-f{family}-m{mode}-s{scale}-np{ranks}';out=r/name;out.mkdir()
     job(name,['mpirun','--oversubscribe','-np',ranks,d/'check',scale,family,mode],out,1024)
     mem=r/(name+'-memory');assert (mem/'command.err').stat().st_size==0,name
     text=(mem/'command.log').read_text();rows=[x.split() for x in text.splitlines() if x.startswith('RESULT ')]
     assert len(rows)==2 and 'CHECK        0' in text
     # Sorting restores canonical physical target order, independent of ranks.
     lines=sorted(line for f in out.glob('rank-*.txt') for line in f.read_text().splitlines())
     digest=hashlib.sha256('\n'.join(lines).encode()).hexdigest();hashes.append(digest)
     results[name]={'sha256':digest,'maximum_error':float(rows[1][-2]),'rms_error':float(rows[1][-1]),'uniform_maximum':float(rows[0][-2]),'receivers':int(rows[1][-3])}
    assert len(set(hashes))==1,(arch,family,mode,scale)
    errors.append(results[name]['rms_error'])
   import math
   orders=[math.log(errors[i]/errors[i+1],2) for i in range(2)]
   results[f'{arch}-f{family}-m{mode}-orders']=orders
   assert min(orders)>1.7,orders
 (r/'results.json').write_text(json.dumps(results,indent=2))
print('All MPI, coverage, uniform preservation and refinement checks passed',flush=True)

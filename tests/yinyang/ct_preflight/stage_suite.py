#!/usr/bin/env python3
"""Compare staged/native CT including a nontrivial edge modification.

Provide a compiled ct_preflight executable and a fresh external output directory.
All MPI runs and comparisons are bounded by the existing memory watchdog.
"""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('executable',type=Path);p.add_argument('workdir',type=Path)
a=p.parse_args();source=Path(__file__).resolve().parents[3];root=a.workdir.resolve()
if root==source or source in root.parents:p.error('external workdir required')
if os.environ.get('AMRVAC_MEMORY_GUARD_ACTIVE')!='1':
 guard=source/'tests/yinyang/interface_transport/convergence/memory_guard.py'
 result=subprocess.run([sys.executable,str(guard),'--out',str(root.with_name(root.name+'-memory')),'--',sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]])
 raise SystemExit(result.returncode)
root.mkdir(parents=True,exist_ok=False)
par=Path(__file__).with_name('amrvac.par').read_text().replace('domain_nx3=48','domain_nx3=40')
# Twenty blocks give unequal block counts on three ranks: the stage collective
# must still execute once per rank, not once per local block.
exe=a.executable.resolve();results=dict(executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),runs={},comparisons={})
def run(name,ranks,split,case,scale=1.,stage_scale=False,ct='uct_hll',expect=None,flux='tvdlf'):
 d=root/name;d.mkdir();shutil.copy2(exe,d/'amrvac')
 config=par.replace('&methodlist',f'&methodlist\n ct_stage_split={".true." if split else ".false."}')
 config=config.replace("type_ct='uct_hll'",f"type_ct='{ct}'").replace("flux_scheme=20*'tvdlf'",f"flux_scheme=20*'{flux}'")
 config+=f'\n&ct_test\n ct_test_case={case}\n emf_scale={scale}\n scale_in_stage={".true." if stage_scale else ".false."}\n/\n'
 (d/'amrvac.par').write_text(config)
 with (d/'run.log').open('w') as out,(d/'run.err').open('w') as err:
  proc=subprocess.Popen(['/usr/local/openmpi/bin/mpirun','--oversubscribe','-np',str(ranks),'./amrvac'],cwd=d,stdout=out,stderr=err,env=dict(os.environ,OMP_NUM_THREADS='1'))
  try:code=proc.wait(timeout=120)
  except subprocess.TimeoutExpired:
   # Let the outer guard kill verified descendants, never orphan MPI children.
   raise RuntimeError(f'{name}: MPI timeout')
 errsize=(d/'run.err').stat().st_size
 if (d/'run.log').stat().st_size>262144 or errsize>262144:raise AssertionError('unbounded logging')
 log=(d/'run.log').read_text()
 if expect is not None:
  assert code!=0 and expect in log+(d/'run.err').read_text()
 else:
  assert code==0 and errsize==0,(name,code,errsize)
  rows=[x.split() for x in log.splitlines() if x.startswith('CT_PREFLIGHT')]
  assert len(rows)==11
  assert max(float(x[-1]) for x in rows)<1e-10
  calls=[int(x.split()[-1]) for x in log.splitlines() if x.startswith('CT_STAGE_READY')]
  assert calls==(list(range(1,16)) if split else [])
 results['runs'][name]=dict(ranks=ranks,returncode=code,stderr_bytes=errsize,split=split,case=case,scale=scale,scale_in_stage=stage_scale)
 (root/'progress.json').write_text(json.dumps(results,indent=2)+'\n')
 return d

def identical(label,left,right):
 hashes=[]
 for number in (0,1):
  f=f'ct{number:04d}.dat'
  lh=hashlib.sha256((left/f).read_bytes()).hexdigest();rh=hashlib.sha256((right/f).read_bytes()).hexdigest()
  assert lh==rh,(label,number)
  hashes.append(lh)
 results['comparisons'][label]=dict(binary_identical=True,sha256=hashes)

for case in (0,1):
 for ranks in (1,2,3):
  left=run(f'c{case}-np{ranks}-native',ranks,False,case)
  right=run(f'c{case}-np{ranks}-split',ranks,True,case)
  identical(f'c{case}-np{ranks}',left,right)
  if ranks>1:identical(f'c{case}-np1-vs-np{ranks}',root/f'c{case}-np1-split',right)
for ranks in (1,3):
 left=run(f'scale-np{ranks}-block',ranks,False,1,.5,False)
 right=run(f'scale-np{ranks}-stage',ranks,True,1,.5,True)
 identical(f'scale-np{ranks}',left,right)
 # The nontrivial modification must actually change the evolved solution.
 assert (left/'ct0001.dat').read_bytes()!=(root/f'c1-np{ranks}-native/ct0001.dat').read_bytes()
left=run('hll-native',2,False,1,flux='hll')
right=run('hll-split',2,True,1,flux='hll')
identical('hll',left,right)
run('reject-contact',2,True,0,ct='uct_contact',expect='split CT requires UCT-HLL')
(root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(dict(runs=len(results['runs']),comparisons=len(results['comparisons']),result='passed')))

#!/usr/bin/env python3
"""Verify that staged CT caches can be rebuilt after a native CT restart."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('executable',type=Path);p.add_argument('suite',type=Path);p.add_argument('workdir',type=Path)
a=p.parse_args();source=Path(__file__).resolve().parents[3];root=a.workdir.resolve()
if root==source or source in root.parents:p.error('external workdir required')
if not (a.suite/'results.json').exists():p.error('completed stage suite required')
if os.environ.get('AMRVAC_MEMORY_GUARD_ACTIVE')!='1':
 guard=source/'tests/yinyang/interface_transport/convergence/memory_guard.py'
 q=subprocess.run([sys.executable,str(guard),'--out',str(root.with_name(root.name+'-memory')),'--',sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]])
 raise SystemExit(q.returncode)
root.mkdir(parents=True,exist_ok=False)
prior=a.suite/'c1-np3-split';par=(prior/'amrvac.par').read_text().replace('it_max=5','it_max=10')
for label in ('continuous','restart'):
 d=root/label;d.mkdir();shutil.copy2(a.executable,d/'amrvac');(d/'amrvac.par').write_text(par)
 cmd=['/usr/local/openmpi/bin/mpirun','--oversubscribe','-np','3','./amrvac']
 if label=='restart':
  shutil.copy2(prior/'ct0001.dat',d/'ct0001.dat');cmd+=['-if','ct0001.dat']
 with (d/'run.log').open('w') as out,(d/'run.err').open('w') as err:
  q=subprocess.run(cmd,cwd=d,stdout=out,stderr=err,timeout=120,env=dict(os.environ,OMP_NUM_THREADS='1'))
 assert q.returncode==0 and (d/'run.err').stat().st_size==0
 expected=30 if label=='continuous' else 15
 calls=[int(x.split()[-1]) for x in (d/'run.log').read_text().splitlines() if x.startswith('CT_STAGE_READY')]
 assert calls==list(range(1,expected+1))
hashes=[]
for label in ('continuous','restart'):
 with (root/label/'ct0002.dat').open('rb') as f:hashes.append(hashlib.file_digest(f,'sha256').hexdigest())
assert hashes[0]==hashes[1]
(root/'results.json').write_text(json.dumps(dict(binary_identical=True,sha256=hashes[0],ranks=3),indent=2)+'\n')
print('staged CT restart binary-identical')

#!/usr/bin/env python3
"""Exercise v2 checkpoint rejection with bounded logs, runtime and memory."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('validated_run',type=Path,help='completed run_validation directory')
p.add_argument('workdir',type=Path)
a=p.parse_args();source=Path(__file__).resolve().parents[2];root=a.workdir.resolve()
if root==source or source in root.parents:p.error('external workdir required')
r=a.validated_run.resolve()
if not (r/'results.json').exists():p.error('positive restart tests must complete first')
root.mkdir(parents=True,exist_ok=False)
guard=source/'tests/yinyang/interface_transport/convergence/memory_guard.py'
results={}
for name,expected in [('missing_completion','missing yinyang checkpoint completion'),
                      ('stale_completion','stale yinyang checkpoint completion'),
                      ('wrong_size','yinyang checkpoint data size mismatch'),
                      ('wrong_mode','yinyang checkpoint mesh, patch, MPI or field mode mismatch'),
                      ('unknown_version','unsupported yinyang checkpoint version'),
                      ('missing_metadata','missing yinyang checkpoint metadata')]:
    d=root/name;d.mkdir()
    shutil.copy2(r/'build-case/amrvac',d/'amrvac')
    shutil.copy2(r/'continuous/amrvac.par',d/'amrvac.par')
    for f in (r/'np2').glob('*0001.dat*'):shutil.copy2(f,d/f.name)
    metadata=d/'yy_yang0001.dat.yy';completion=d/'yy_yang0001.dat.yy.complete'
    if name=='missing_completion':completion.unlink()
    elif name=='stale_completion':completion.write_text('0 0.0 1\n')
    elif name=='missing_metadata':metadata.unlink()
    else:
        lines=metadata.read_text().splitlines();header=lines[0].split()
        if name=='unknown_version':header[0]='99'
        elif name=='wrong_mode':header[6]='1'
        elif name=='wrong_size':
            row=lines[1].split();row[2]=str(int(row[2])+1);lines[1]=' '.join(row)
        lines[0]=' '.join(header);metadata.write_text('\n'.join(lines)+'\n')
    cmd=[sys.executable,str(guard),'--timeout-seconds','30','--out',str(d/'guard'),'--',
         '/usr/local/openmpi/bin/mpirun','--oversubscribe','-np','2','./amrvac','-if','yy_{patch}0001.dat']
    result=subprocess.run(cmd,cwd=d,env=dict(os.environ,OMP_NUM_THREADS='1'))
    summary=json.loads((d/'guard/summary.json').read_text())
    text=''
    for f in ('command.log','command.err'):
        path=d/'guard'/f
        if path.stat().st_size>262144:raise AssertionError('unbounded failure logging')
        text+=path.read_text(errors='replace')
    assert result.returncode!=0 and summary['stop_reason'] is None,(name,summary)
    assert expected in text,(name,text[-2000:])
    results[name]=dict(returncode=result.returncode,diagnostic=expected,wall=summary['wall'])
    (root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))

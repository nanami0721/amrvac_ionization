#!/usr/bin/env python3
from pathlib import Path
import argparse,json,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('workdir',type=Path);p.add_argument('--time-step-control',action='store_true');a=p.parse_args();r=a.workdir.resolve();s=Path(__file__).resolve().parents[3];defs=s/'tests/yinyang/hydrostatic';guard=s/'tests/yinyang/interface_transport/convergence/memory_guard.py'
results={}
def guarded(cmd,out):
 with out.with_suffix('.driver.log').open('w') as f:
  q=subprocess.run([sys.executable,str(guard),'--out',str(out),'--',*map(str,cmd)],stdout=f,stderr=f)
 if q.returncode:raise SystemExit('Stopped: '+str(out))
for direction in (['halfdt'] if a.time_step_control else ['radial','theta','phi']):
 for grid in ['yinyang','spherical']:
  name=f'refine-{direction}-{grid}';d=r/name
  flags=['--radial-refine'] if direction=='radial' else ['--refine-direction',direction]
  if direction=='halfdt':flags=['--cfl',.15]
  guarded([sys.executable,defs/'run.py',d,r/'case-default/amrvac','--case',3,'--grid',grid,'--beta',.01,'--crossings',.1,'--output-interval',.02,*flags],r/(name+'-guard'))
  guarded([sys.executable,defs/'analyze.py',s,d],r/(name+'-analysis-guard'))
  q=json.loads((d/'analysis.json').read_text());results[name]={'final':[v for v in q['rows'] if v['region']=='common'][-1],'engineering_pass':q['engineering_pass']};(r/('time-control.json' if a.time_step_control else 'refinement.json')).write_text(json.dumps(results,indent=2));print(name,q['engineering_pass'],flush=True)

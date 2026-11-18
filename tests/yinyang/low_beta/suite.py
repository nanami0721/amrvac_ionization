#!/usr/bin/env python3
from pathlib import Path
import argparse,json,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('workdir',type=Path);a=p.parse_args()
r=a.workdir.resolve();source=Path(__file__).resolve().parents[3];defs=source/'tests/yinyang/hydrostatic';guard=source/'tests/yinyang/interface_transport/convergence/memory_guard.py'
assert (r/'build-provenance.json').exists()
results={};blocked=set()
def save(): (r/'suite.json').write_text(json.dumps(results,indent=2))
def guarded(cmd,out):
 with out.with_suffix('.driver.log').open('w') as log:
  q=subprocess.run([sys.executable,str(guard),'--out',str(out),'--',*map(str,cmd)],stdout=log,stderr=log)
 summary=json.loads((out/'summary.json').read_text())
 if summary['stop_reason']:raise SystemExit('Memory/monitor stop: '+str(out))
 return q.returncode

def run(name,case,grid,beta,duration,debug=False):
 d=r/name;cmd=[sys.executable,defs/'run.py',d,r/('case-debug' if debug else 'case-default')/'amrvac','--case',case,'--grid',grid,'--beta',beta,'--crossings',duration,'--output-interval',.02 if duration==.1 else 1]
 if duration==1 and beta!=4:cmd[-1]=.1
 if debug:cmd+=['--debug']
 code=guarded(cmd,r/(name+'-guard'))
 if code:
  results[name]={'run_failed':True};save();return None
 if debug:
  meta=json.loads((d/'run.json').read_text());assert meta['stderr']==0 and len(list(d.glob('equilibrium_*.bin')))==4*meta['ranks'];return d
 if guarded([sys.executable,defs/'analyze.py',source,d],r/(name+'-analysis-guard')):raise SystemExit('Analysis failure: '+name)
 q=json.loads((d/'analysis.json').read_text());common=[x for x in q['rows'] if x['region']=='common']
 assert abs(common[-1]['crossings']-duration)<1e-10
 passed=all(x['mach_rms']<.01 and x['mach_max']<.05 and x['rho_L1']<.02 and x['pressure_L1']<.02 and x['B_L1']<.02 for x in common)
 results[name]={'run_failed':False,'all_saved_times_pass':passed,'final':common[-1]};save();print(name,'pass',passed,flush=True);return d

for case in [3,2]:
 for grid in ['yinyang','spherical']:
  if run(f'debug-c{case}-{grid}',case,grid,.004,.1,True) is None:raise SystemExit('Debug run failed')
results['debug']='4 runs; stderr empty; repeated native dumps verified';save()
d=run('baseline-beta4',2,'yinyang',4,1)
assert d
old=r.parent/'yinyang-gravity-20261117'
out=r/'baseline-comparison'
assert guarded([sys.executable,old/'compare_snapshots.py',source,d,old/'short-yinyang-c2','last','eq_yin','eq_yang'],out)==0
results['baseline']=json.loads((out/'command.log').read_text());save()
for beta in [1,.1,.01,.004]:
 for case in [3,2]:
  for grid in ['yinyang','spherical']:
   name=f'short-c{case}-{grid}-b{beta:g}'
   if (case,grid) in blocked:results[name]={'skipped':'lower beta blocked after numerical failure'};save();continue
   if run(name,case,grid,beta,.1) is None:blocked.add((case,grid))
for beta in [1,.1,.01,.004]:
 for case in [3,2]:
  pair=[results[f'short-c{case}-{grid}-b{beta:g}'] for grid in ['yinyang','spherical']]
  if not all(x.get('all_saved_times_pass',False) for x in pair):continue
  for grid in ['yinyang','spherical']:run(f'long-c{case}-{grid}-b{beta:g}',case,grid,beta,1)
results['completed']=True;save()

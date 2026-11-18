#!/usr/bin/env python3
"""Explicit checks; does not silently certify formal second-order accuracy."""
from pathlib import Path
import json,math,sys
root=Path(sys.argv[1]);rows=json.loads((root/'results.json').read_text());final=[r for r in rows if r['snapshot']==2];by={(r['case'],r['placement'],r['direction'],r['scale']):r for r in final}
checks={};orders={}
checks['expected_complete_runs']=len(final)==16
debug=[json.loads(p.read_text()) for p in (root/'debug-runs').glob('*/run.json')]
checks['four_debug_runs']=len(debug)==4
checks['debug_success_and_empty_stderr']=all(r['returncode']==0 and r['stderr']==0 for r in debug)
checks['finite_diagnostics']=all(math.isfinite(v) for r in rows for v in r.values() if isinstance(v,float))
checks['initial_exact_solution']=max(r['relative_L1'] for r in rows if r['snapshot']==0)<1e-12
checks['positive_pressure_density']=all(r['pressure_min']>0 and r['rho_min']>0 for r in rows)
checks['empty_runtime_stderr']=all(r['stderr']==0 for r in final)
for case in [0,1]:
 for direction in [1,-1]:
  rr=[by[(case,1,direction,s)] for s in [1,2,3]]
  for metric in ['relative_L1','relative_L2','roi_L1','seam_L1_contribution','div_defect_rms']:
   ee=[r[metric] for r in rr];checks[f'c{case}_d{direction}_{metric}_decreases']=ee[0]>ee[1]>ee[2]
   orders[f'c{case}_d{direction}_{metric}']=[math.log(ee[0]/ee[1])/math.log(2),math.log(ee[1]/ee[2])/math.log(1.5)]
result=dict(checks=checks,observed_orders=orders,all_pass=all(checks.values()))
(root/'assessment.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
if not result['all_pass']:raise SystemExit(1)

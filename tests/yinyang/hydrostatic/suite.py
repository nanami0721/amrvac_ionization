#!/usr/bin/env python3
"""Serial guarded validation. Requires the external build phase to have succeeded."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys
sys.dont_write_bytecode=True
p=argparse.ArgumentParser();p.add_argument('workdir',type=Path);a=p.parse_args();r=a.workdir.resolve();defs=Path(__file__).resolve().parent;source=defs.parents[2];guard=defs.parent/'interface_transport/convergence/memory_guard.py'
for label in ['equilibrium-default','equilibrium-debug','uniform-default','transport-default']:
 if not (r/label/'amrvac').exists():raise SystemExit('Complete guarded prepare_build.py before suite.py')
if not (r/'build-provenance.json').exists():raise SystemExit('Missing completed build provenance')
results={}
def guarded(cmd,out,cwd=r):
 q=subprocess.run([sys.executable,str(guard),'--out',str(out),'--',*map(str,cmd)],cwd=cwd)
 if q.returncode:raise SystemExit(q.returncode)
 return json.loads((out/'summary.json').read_text())
def equilibrium(name,case,grid='yinyang',arch='default',extra=()):
 d=r/name;guarded([sys.executable,defs/'run.py',d,r/f'equilibrium-{arch}/amrvac','--case',case,'--grid',grid,*extra],r/(name+'-guard'))
 return d
def analyze(d):
 guarded([sys.executable,defs/'analyze.py',source,d],r/(d.name+'-analysis-guard'))
 return json.loads((d/'analysis.json').read_text())
# The analysis/comparison subprocess uses its own memory guard as well.
def compare(left,right,prefixes,time=None):
 script=r/'compare_snapshots.py'
 if not script.exists():
  script.write_text('''from pathlib import Path
import json,sys,numpy as np
sys.path.insert(0,str(Path(sys.argv[1])/'tests/yinyang'));import check_mhd
u=check_mhd.reader(sys.argv[1]);left,right=map(Path,sys.argv[2:4]);time=None if sys.argv[4]=='last' else float(sys.argv[4]);maximum=0.
def select(d,prefix):
 candidates=[]
 for f in sorted(d.glob(prefix+'[0-9][0-9][0-9][0-9].dat')):
  with f.open('rb') as stream:h=u.get_header(stream)
  if time is None or abs(h['time']-time)<1e-12:candidates.append((h['time'],f))
 assert candidates
 return max(candidates)[1]
for prefix in sys.argv[5:]:
 ha,a=check_mhd.load(u,select(left,prefix));hb,b=check_mhd.load(u,select(right,prefix));assert ha['time']==hb['time'];maximum=max(maximum,float(abs(a-b).max()));assert np.allclose(a,b,atol=1e-12,rtol=1e-12)
print(json.dumps(dict(max_difference=maximum)))
''')
 out=r/('compare-'+left.name+'-'+right.name)
 guarded([sys.executable,script,source,left,right,'last' if time is None else time,*prefixes],out)
 return json.loads((out/'command.log').read_text())
# Debug callbacks repeat every step; the source self-check uses nonzero velocity.
for case in [0,1,2]:
 d=equilibrium(f'debug-c{case}',case,arch='debug',extra=['--debug']);meta=json.loads((d/'run.json').read_text());assert meta['stderr']==0
 assert len(list(d.glob('equilibrium_*.bin')))==4*meta['ranks']
# Uniform and forward/reverse transport exercise the unchanged no-gravity path.
for name,exe,par,np_ in [('uniform-2','uniform-default',(source/'tests/yinyang/mhd/amrvac.par').read_text(),2),('uniform-4','uniform-default',(source/'tests/yinyang/mhd/amrvac.par').read_text(),4)]:
 d=r/name;d.mkdir();shutil.copy2(r/exe/'amrvac',d/'amrvac');(d/'amrvac.par').write_text(par);guarded(['/usr/local/openmpi/bin/mpirun','--oversubscribe','-np',np_,'./amrvac'],d/'guard',d);assert not (d/'guard/command.err').stat().st_size
results['uniform_partition']=compare(r/'uniform-2',r/'uniform-4',['yy_yin','yy_yang'])
baseline=r.parent/'yinyang-glm-20261116/verified-fix/runs/c1-p1-d1-s1'
for direction in [1,-1]:
 d=r/f'transport-{direction}';d.mkdir();shutil.copy2(r/'transport-default/amrvac',d/'amrvac');par=(baseline/'amrvac.par').read_text().replace('direction=1',f'direction={direction}').replace('native_diagnostics=.true.','native_diagnostics=.false.');(d/'amrvac.par').write_text(par)
 guarded(['/usr/local/openmpi/bin/mpirun','--oversubscribe','-np','2','./amrvac'],d/'guard',d);assert not (d/'guard/command.err').stat().st_size
results['transport_baseline']=compare(r/'transport-1',baseline,['adv_yin','adv_yang'])
# First collect all separated short controls before deciding on long runs.
passes=[]
for case in [0,1,2]:
 for grid in ['yinyang','spherical']:
  d=equilibrium(f'short-{grid}-c{case}',case,grid);v=analyze(d);passes.append(v['engineering_pass']);results[d.name]=v['engineering_pass']
# Read-only diagnostics and MPI partition checks at identical final time.
off=equilibrium('short-no-diagnostics',2,extra=['--no-diagnostics'])
results['diagnostics_invariant']=compare(r/'short-yinyang-c2',off,['eq_yin','eq_yang'])
four=equilibrium('short-four-ranks',2,extra=['--ranks','4'])
results['equilibrium_partition']=compare(r/'short-yinyang-c2',four,['eq_yin','eq_yang'])
(r/'suite.json').write_text(json.dumps(results,indent=2))
if not all(passes):raise SystemExit('Short equilibrium criteria failed; diagnose separated cases before long evolution')
for case in [0,1,2]:
 for grid in ['yinyang','spherical']:
  d=equilibrium(f'long-{grid}-c{case}',case,grid,extra=['--crossings','10']);v=analyze(d);results[d.name]=v['engineering_pass'];(r/'suite.json').write_text(json.dumps(results,indent=2))
restart=equilibrium('restart-c2',2,extra=['--crossings','3','--restart',str(r/'short-yinyang-c2')]);analyze(restart)
tau=json.loads((restart/'run.json').read_text())['tau'];results['restart']=compare(r/'long-yinyang-c2',restart,['eq_yin','eq_yang'],3*tau)
peak=json.loads((r/'short-yinyang-c2-guard/summary.json').read_text())['peak_rss_mib']
if 2*peak<.7*1024:
 d=equilibrium('radial-refined-c2',2,extra=['--radial-refine']);results['radial_refined']=analyze(d)['engineering_pass']
else:results['radial_refined']='skipped: insufficient measured memory margin'
(r/'suite.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))

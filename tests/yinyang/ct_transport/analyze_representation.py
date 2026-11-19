#!/usr/bin/env python3
"""Audit paired evolved EMF-representation comparisons without assuming improvement."""
import argparse,json,math
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plumbing-reference',type=Path);p.add_argument('--face-probe',action='store_true');p.add_argument('--legacy-paired-results',type=Path);p.add_argument('workdir',type=Path);p.add_argument('--reference-results',type=Path,required=True);a=p.parse_args();r=a.workdir
results=json.loads((r/'results.json').read_text())
if a.plumbing_reference:
 assert a.face_probe and len(results)==2
 prior=json.loads((a.plumbing_reference/'results.json').read_text())
 for key,row in results.items():
  assert row['hashes']==prior[key]['hashes'],key
  assert row['initial']==prior[key]['initial'] and row['final']==prior[key]['final'],key
 audit=dict(default_and_opt_in_bitwise=True,matched_cases=list(results))
 (r/'comparison.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2));raise SystemExit(0)
assert len(results)==(4 if a.face_probe else 9)
comparisons={}
for method in (['hllmc'] if a.face_probe else ['base','hllmc']):
 for direction in [1,-1]:
  old=f'legacy-{method}-d{direction}-np2';new=f'midpoint-{method}-d{direction}-np2';b,c=results[old],results[new]
  # Initial active states must be unchanged; compare exact printed cell samples.
  def initial(name):return sorted(line for f in (r/name).glob('sample-0-rank-*.txt') for line in f.read_text().splitlines())
  before,after=initial(old),initial(new);assert before and before==after
  assert all(math.isfinite(v) for row in [b['initial'],b['final'],c['initial'],c['final']] for v in row)
  assert c['max_active_divergence']<1e-10 and min(c['final'][13:15])>0
  if a.legacy_paired_results:
   prior=json.loads((a.legacy_paired_results/'results.json').read_text())[f'midpoint-{method}-d{direction}-np2']
   assert max(abs(x-y) for x,y in zip(prior['final'],b['final']))<1e-12
  elif method=='base' or direction==1:
   golden=a.reference_results/f'yy-p0-s2-d{direction}-{method}'/'transport-history.txt'
   h=[[float(x) for x in line.split()] for line in golden.read_text().splitlines()]
   assert max(abs(x-y) for x,y in zip(h[-1],b['final']))<1e-12
  comparisons[f'{method}-d{direction}']=dict(old_l2=b['l2'],new_l2=c['l2'],relative_l2_reduction=1-c['l2']/b['l2'],old_max_error=b['max_error'],new_max_error=c['max_error'],max_active_divergence=c['max_active_divergence'])
if not a.face_probe:assert results['midpoint-hllmc-d1-np2']['hashes']==results['midpoint-hllmc-d1-np4']['hashes']
mem={p.parent.name:json.loads(p.read_text()) for p in r.glob('*-memory/summary.json')}
runtime=[v for k,v in mem.items() if '-np' in k]
audit=dict(comparisons=comparisons,mpi_bitwise=(None if a.face_probe else True),initial_samples_bitwise=True,peak_runtime_mib=max(max(x['peak_rss_mib'],x['peak_footprint_mib']) for x in runtime),peak_system_percent=max(x['peak_system_percent'] for x in mem.values()))
(r/'comparison.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))

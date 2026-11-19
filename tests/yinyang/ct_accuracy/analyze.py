#!/usr/bin/env python3
"""Summarize transfer errors; sample RMS is not a volume-weighted conservation norm."""
import argparse,json,math
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('workdir',type=Path);a=p.parse_args();root=a.workdir
raw=json.loads((root/'results.json').read_text());out={}
for name,rows in raw.items():
 out[name]={}
 for row in rows:
  v=row['values'];field=row['field'];amp=1 if field<2 else .001
  def rms(s,n):return math.sqrt(s/n)/amp
  if row['kind']=='EDGE_ACCURACY':
   d=dict(edge_rms=rms(v[0],v[3]),midpoint_rms=rms(v[2],v[3]),edge_max=v[14],seam_curl_rms=rms(v[8],v[11]),midpoint_seam_curl_rms=rms(v[10],v[11]),curl_change_max=v[15],bulk_curl_change_max=v[16])
  else:
   d=dict(face_rms=rms(v[0],v[2]),face_max=v[14],divergence_max=v[15],component_rms=[rms(v[3*c],v[3*c+2]) for c in [1,2,3]])
  out[name][f"{row['kind']}:{field}:{row['mode']}"]=d
# These contracts test actual identities, not desired transport performance.
for name,metrics in out.items():
 for field in range(4):
  assert metrics[f'EDGE_ACCURACY:{field}:0']['bulk_curl_change_max']==0
  assert metrics[f'GHOST_ACCURACY:{field}:0']['face_max']<1e-11
  for stage in [0,2]:assert metrics[f'GHOST_ACCURACY:{field}:{stage}']['divergence_max']<1e-10
  if field<2:
   for stage in [1,2]:assert metrics[f'GHOST_ACCURACY:{field}:{stage}']['face_max']<1e-10
 native=metrics['EDGE_ACCURACY:0:3']
 assert abs(native['seam_curl_rms']/native['midpoint_seam_curl_rms']-1)<1e-8
 for field in range(4):assert metrics[f'EDGE_ACCURACY:{field}:3']['bulk_curl_change_max']==0
 constant=metrics['EDGE_ACCURACY:0:2']
 assert abs(constant['seam_curl_rms']/constant['midpoint_seam_curl_rms']-1)<1e-3
orders={}
for n0,n1 in [(8,12),(12,14),(8,14)]:
 lo=out[f'debug-n{n0}-np2-q8'];hi=out[f'debug-n{n1}-np2-q8']
 orders[f'{n0}-{n1}']={key:{m:math.log(v/hi[key][m])/math.log(n1/n0) for m,v in val.items() if isinstance(v,float) and v>1e-12 and hi[key][m]>1e-12} for key,val in lo.items()}
def records(name):
 # Unique global geometry keys remove same-patch block-boundary duplication.
 result={}
 for f in (root/name).glob('accuracy-rank-*.txt'):
  for line in f.open():
   row=line.split();key=tuple(row[:8]);value=tuple(map(float,row[8:]));prev=result.get(key)
   if prev is not None:assert max(abs(x-y) for x,y in zip(prev,value))<1e-13
   result[key]=value
 return result
def compare(n1,n2):
 r1,r2=records(n1),records(n2);assert r1.keys()==r2.keys()
 delta=max(abs(x-y) for key in r1 for x,y in zip(r1[key],r2[key]));return dict(records=len(r1),max_absolute_difference=delta)
checks=dict(mpi=compare('debug-n8-np2-q8','debug-n8-np4-q8'),quadrature=compare('debug-n14-np2-q8','debug-n14-np2-q16'),optimized={str(n):compare(f'debug-n{n}-np2-q8',f'default-n{n}-np2-q8') for n in [8,12,14]})
mem=[json.loads(f.read_text()) for f in root.glob('*-memory/summary.json')]
checks['memory']=dict(peak_process_mib=max(max(x['peak_rss_mib'],x['peak_footprint_mib']) for x in mem),peak_system_percent=max(x['peak_system_percent'] for x in mem))
assert checks['mpi']['max_absolute_difference']<1e-13
assert checks['quadrature']['max_absolute_difference']<1e-9
assert all(v['max_absolute_difference']<1e-12 for v in checks['optimized'].values())
(root/'analysis.json').write_text(json.dumps(dict(metrics=out,orders=orders,checks=checks),indent=2))
print(json.dumps(checks,indent=2))
for field in [1,2,3]:
 print('field',field)
 for n in [8,12,14]:
  m=out[f'debug-n{n}-np2-q8'];print(n,'E RMS/curl',*[m[f'EDGE_ACCURACY:{field}:0'][k] for k in ['edge_rms','seam_curl_rms']], 'B raw/closed',*[m[f'GHOST_ACCURACY:{field}:{stage}']['face_rms'] for stage in [1,2]])

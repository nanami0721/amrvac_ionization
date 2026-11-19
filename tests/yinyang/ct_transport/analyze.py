#!/usr/bin/env python3
"""Small-memory summaries of native CT transport histories (no snapshot loading)."""
import argparse,json,math,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('workdir',type=Path);p.add_argument('--plot',action='store_true');p.add_argument('--complete',action='store_true');a=p.parse_args();r=a.workdir
results={}
for f in sorted(r.glob('*-d*-p*-s*-np*-a*/transport-history.txt')):
 name=f.parent.name
 m=re.fullmatch(r'(debug|default)-d(-?1)-p([01])-s([12])-np([24])-a([\d.]+)',name)
 if not m:continue
 arch,d,placement,scale,ranks,amp=m.groups();d=int(d);placement=int(placement)
 h=[[float(v) for v in line.split()] for line in f.read_text().splitlines()]
 if h[-1][0]!=300:continue
 def metrics(v):
  phi=math.radians(146.25*placement);velocity=[-d*math.sin(phi),d*math.cos(phi),0]
  center=[1.5*math.cos(phi)+(v[1]-.3)*velocity[0],1.5*math.sin(phi)+(v[1]-.3)*velocity[1],0]
  centroid=[v[i]/v[22] for i in range(19,22)]
  return dict(time=v[1],relative_l2=math.sqrt(v[15]/v[16]),relative_l1=v[17]/v[18],peak_fraction=v[11]/(float(amp)*.001),centroid=centroid,centroid_phi_deg=math.degrees(math.atan2(centroid[1],centroid[0])),centroid_error=math.dist(center,centroid),rho_min=v[13],p_min=v[14],velocity_error_max=v[7],yin_reference_energy_fraction=v[25]/v[26])
 results[name]=dict(initial=metrics(h[0]),final=metrics(h[-1]),max_active_divb=max(v[2] for v in h),max_ghost_divb=max(v[4] for v in h),max_seam_divb=max(v[3] for v in h))
 phi=math.radians(146.25*placement);velocity=[-d*math.sin(phi),d*math.cos(phi),0]
 results[name]['centroid_displacement_along_flow']=sum((results[name]['final']['centroid'][i]-results[name]['initial']['centroid'][i])*velocity[i] for i in range(3))
 if placement:
  crossing=.3+1.5*math.tan(math.radians(140.625-146.25))/d
  results[name]['center_crossing_time']=crossing
  results[name]['at_crossing']=metrics(min(h,key=lambda v:abs(v[1]-crossing)))
checks={}
for d in [1,-1]:
 for placement in [0,1]:
  keys=[f'debug-d{d}-p{placement}-s{s}-np2-a1' for s in [1,2]]
  if all(k in results for k in keys):checks[f'd{d}-p{placement}-fine_over_coarse_l2']=results[keys[1]]['final']['relative_l2']/results[keys[0]]['final']['relative_l2']
full='debug-d1-p1-s1-np2-a1';half='debug-d1-p1-s1-np2-a0.5'
if full in results and half in results:checks['amplitude_halving_relative_l2_change']=results[half]['final']['relative_l2']/results[full]['final']['relative_l2']-1
memory=[json.loads(f.read_text()) for f in r.glob('*-memory/summary.json')]
geometry=json.loads((r/'geometry.json').read_text()) if (r/'geometry.json').exists() else dict(radius=.35,radial_min=1.,radial_max=2.)
summary=dict(geometry=geometry,cases=results,comparisons=checks,memory=dict(peak_system_percent=max(v['peak_system_percent'] for v in memory),peak_build_rss_mib=max(v['peak_rss_mib'] for v in memory if v['limit_mib']==4096),peak_run_rss_mib=max(v['peak_rss_mib'] for v in memory if 'mpirun' in v['command'])))
(r/'analysis.json').write_text(json.dumps(summary,indent=2))
for k,v in results.items():print(k,'L2',round(v['initial']['relative_l2'],4),'->',round(v['final']['relative_l2'],4),'peak',round(v['final']['peak_fraction'],4),'div',format(v['max_active_divb'],'.3e'))
print(json.dumps(checks));print(json.dumps(summary['memory']))

if a.plot:
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(1,3,figsize=(13,3.7),layout='constrained')
 for name in results:
  if not name.startswith('debug-') or '-np2-a1' not in name:continue
  h=[[float(v) for v in line.split()] for line in (r/name/'transport-history.txt').read_text().splitlines()]
  fine='-s2-' in name;seam='-p1-' in name;reverse='-d-1-' in name
  label=('seam' if seam else 'interior')+(' reverse' if reverse else ' forward')+(' fine' if fine else ' coarse')
  style=dict(color=('tab:blue' if seam else 'tab:orange'),linestyle=('--' if reverse else '-'),alpha=(1 if fine else .45),linewidth=(1.8 if fine else 1))
  t=[v[1] for v in h]
  axes[0].plot(t,[math.sqrt(v[15]/v[16]) for v in h],label=label,**style)
  axes[1].plot(t,[v[11]/.001 for v in h],**style)
  axes[2].semilogy(t,[max(v[2],1e-20) for v in h],**style)
 for ax,title in zip(axes,['Relative L2 magnetic error','Peak |B| / analytic peak','Native max active |div B|']):
  ax.set(xlabel='Time',title=title);ax.grid(alpha=.2)
 axes[0].legend(fontsize=6.5)
 fig.suptitle(f"Loop radius {geometry['radius']}; shell [{geometry['radial_min']}, {geometry['radial_max']}]")
 fig.savefig(r/'transport-summary.png',dpi=160)

if a.complete:
 assert len(results)==11, f'expected 11 completed cases, got {len(results)}'
 assert all(v['max_active_divb']<1e-10 and v['max_ghost_divb']<1e-10 for v in results.values())
 assert all(v['final']['rho_min']>0 and v['final']['p_min']>0 for v in results.values())
 assert all(checks[f'd{d}-p{p}-fine_over_coarse_l2']<1 for d in [1,-1] for p in [0,1])
 assert abs(checks['amplitude_halving_relative_l2_change'])<1e-3
 raw=json.loads((r/'results.json').read_text())
 assert raw['debug-d1-p1-s1-np4-a1']['hashes']==raw[full]['hashes']
 assert abs(results['default-d1-p1-s1-np2-a1']['final']['relative_l2']-results[full]['final']['relative_l2'])<1e-8
 assert all(v['stop_reason'] is None and v['returncode']==0 for v in memory)
 assert all((r/(name+'-memory/command.err')).stat().st_size==0 for name in results)
 print('Completed matrix: stability, refinement trend, weak-field feedback, MPI and build-mode checks passed.')

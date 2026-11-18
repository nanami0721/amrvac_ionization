#!/usr/bin/env python3
"""Process one completed run; use common definitions without retaining all snapshots."""
import sys,json
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import check_mhd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from diagnostics import basis,rotation
s,d=map(Path,sys.argv[1:3]);u=check_mhd.reader(s)
meta=json.loads((d/'run.json').read_text());assert meta['returncode']==0 and meta['stderr']==0
width=meta['width'];direction=meta.get('direction',1);rows=[]
for f in sorted(d.glob('adv[0-9][0-9][0-9][0-9].dat')):
 h,w=check_mhd.load(u,f);nd=len(h['domain_nx']);dx=(h['xmax']-h['xmin'])/h['domain_nx']
 axes=[a+(np.arange(n)+.5)*dd for a,n,dd in zip(h['xmin'],h['domain_nx'],dx)]
 if nd==1:
  xi=axes[0]+direction*(.3-h['time']);phase=np.pi/2*np.exp(-(xi/width)**2)
  ex=np.stack([np.zeros_like(phase),.3*np.cos(phase),.3*np.sin(phase)],-1)
  B=w[...,5:8];vol=np.full(xi.shape,dx[0]);roi=abs(xi)<.54
  signal=np.linalg.norm(ex-[0,.3,0],axis=-1);masks=dict(roi=roi,all=np.ones_like(roi))
 else:
  r,t,p=np.meshgrid(*axes,indexing='ij');b=basis(t,p,0);pos=r[...,None]*b[..., :,0];rot=rotation(1)
  center=rot@np.array([direction*(h['time']-.3),0,1.5]);offset=pos-center
  phase=np.pi/2*np.exp(-(np.sum(offset*rot[:,0],axis=-1)/width)**2)
  ex=np.stack([np.zeros_like(phase),.3*np.cos(phase),.3*np.sin(phase)],-1)@rot.T
  B=np.einsum('...ij,...j->...i',b,w[...,5:8]);vol=((r+dx[0]/2)**3-(r-dx[0]/2)**3)/3*(np.cos(t-dx[1]/2)-np.cos(t+dx[1]/2))*dx[2]
  roi=np.linalg.norm(offset,axis=-1)<.54;core=(r>1.2)&(r<1.8)
  signal=np.linalg.norm(ex-.3*rot[:,1],axis=-1)
  masks=dict(roi=roi,roi_core=roi&core,fixed_core=core,radial_edge=~core,all=np.ones_like(roi))
 error=np.linalg.norm(B-ex,axis=-1)
 pressure=(h['gamma']-1)*(w[...,4]-.5*np.sum(w[...,1:4]**2,axis=-1)/w[...,0]-.5*np.sum(w[...,5:8]**2,axis=-1))
 assert pressure.min()>0 and w[...,0].min()>0
 for name,mask in masks.items():
  v=vol[mask];e=error[mask];sig=signal[mask]
  rows.append(dict(time=h['time'],steps=h['it'],region=name,L1=float(np.sum(e*v)/np.sum(sig*v)),L2=float(np.sqrt(np.sum(e*e*v)/np.sum(sig*sig*v))),max=float(e.max()),volume=float(v.sum()),signal_integral=float(np.sum(sig*v)),rho_min=float(w[...,0].min()),pressure_min=float(pressure.min()),B_strength_error_max=float(abs(np.linalg.norm(B,axis=-1)-.3).max()),psi_max=float(abs(w[...,8]).max()) if w.shape[-1]>8 else 0.))
(d/'errors.json').write_text(json.dumps(rows,indent=2)+'\n')
assert max(x['L1'] for x in rows if x['time']==0)<1e-12
print(d.name,[(x['region'],x['L1']) for x in rows if abs(x['time']-.6)<1e-12],flush=True)

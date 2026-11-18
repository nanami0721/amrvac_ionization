#!/usr/bin/env python3
import sys,json,re
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import check_mhd
src,root=map(Path,sys.argv[1:3]);util=check_mhd.reader(src)

def rotation(tr):
 if tr==0:return np.eye(3)
 if tr==1:return np.array([[0,0,1],[0,1,0],[-1,0,0.]])
 p=np.deg2rad(146.25);return np.array([[-np.sin(p),0,np.cos(p)],[np.cos(p),0,np.sin(p)],[0,1,0.]])

def geometry(h,patch,yy):
 dx=(h['xmax']-h['xmin'])/h['domain_nx'];axes=[a+(np.arange(n)+.5)*d for a,n,d in zip(h['xmin'],h['domain_nx'],dx)]
 r,t,p=np.meshgrid(*axes,indexing='ij');st,ct,sp,cp=np.sin(t),np.cos(t),np.sin(p),np.cos(p)
 basis=np.stack([np.stack([st*cp,st*sp,ct],-1),np.stack([ct*cp,ct*sp,-st],-1),np.stack([-sp,cp,np.zeros_like(p)],-1)],-1)
 if patch:basis=basis[...,[0,2,1],:];basis[...,0,:]*=-1
 solid=(np.cos(t[0]-dx[1]/2)-np.cos(t[0]+dx[1]/2))*dx[2]
 if yy:
  solid=np.zeros_like(t[0]);ns=16
  for jt in range(ns):
   tt=axes[1][:,None]+((jt+.5)/ns-.5)*dx[1];ds=np.cos(tt-dx[1]/(2*ns))-np.cos(tt+dx[1]/(2*ns))
   for jp in range(ns):
    pp=axes[2][None,:]+((jp+.5)/ns-.5)*dx[2]
    ot=np.arctan2(np.hypot(np.sin(tt)*np.cos(pp),np.cos(tt)),np.sin(tt)*np.sin(pp));op=np.arctan2(np.cos(tt),-np.sin(tt)*np.cos(pp))
    overlap=(ot>=h['xmin'][1])&(ot<=h['xmax'][1])&(op>=h['xmin'][2])&(op<=h['xmax'][2])
    solid+=ds*dx[2]/ns*np.where(overlap,.5,1.)
 vol=((r+dx[0]/2)**3-(r-dx[0]/2)**3)/3*solid[None]
 return r[...,None]*basis[...,:,0],basis,vol

results=[]
for d in sorted(root.glob('*-s*-t*')):
 if not (d/'run.json').exists():continue
 run=json.loads((d/'run.json').read_text())
 if run['returncode']:continue
 yy=d.name.startswith('yinyang');tr=int(d.name[-1]);rot=rotation(tr);vel=rot@np.array([1,0,0.]);bg=rot@np.array([.2,.1,.3]);origin=rot@np.array([-.3,0,1.5])
 for snap in [0,1,2]:
  l1=l2=den1=den2=volsum=roierr=roiden=mass=0.;peak=exactpeak=velerr=berr=0.;centroid=np.zeros(3);pmin=1e99;rhomin=1e99;linf=0;divnorm=divvol=0.;divmax=0.
  for patch,suffix in enumerate(['_yin','_yang'] if yy else ['']):
   h,w=check_mhd.load(util,d/f'adv{suffix}{snap:04d}.dat');pos,basis,vol=geometry(h,patch,yy);center=origin+vel*h['time'];r2=np.sum((pos-center)**2,-1);pulse=.5*np.exp(-r2/.18**2);err=abs(w[...,0]-1-pulse);roi=r2<(.54)**2
   v=np.einsum('...ij,...j->...i',basis,w[...,1:4])/w[...,0,None];b=np.einsum('...ij,...j->...i',basis,w[...,5:8]);pr=(h['gamma']-1)*(w[...,4]-.5*w[...,0]*np.sum(v*v,-1)-.5*np.sum(b*b,-1))
   l1+=np.sum(err*vol);l2+=np.sum(err**2*vol);den1+=np.sum(pulse*vol);den2+=np.sum(pulse**2*vol);volsum+=vol.sum();roierr+=np.sum(err[roi]*vol[roi]);roiden+=np.sum(pulse[roi]*vol[roi]);linf=max(linf,float(err.max()))
   peak=max(peak,float(w[...,0].max()-1));exactpeak=max(exactpeak,float(pulse.max()));velerr+=np.sum(np.sum((v-vel)**2,-1)*vol);berr+=np.sum(np.sum((b-bg)**2,-1)*vol)
   dx=(h['xmax']-h['xmin'])/h['domain_nx'];rad=np.linalg.norm(pos,axis=-1);theta=h['xmin'][1]+(np.arange(h['domain_nx'][1])+.5)*dx[1];st=np.sin(theta)[None,:,None]
   div=np.gradient(rad*rad*w[...,5],dx[0],axis=0)/(rad*rad)+np.gradient(st*w[...,6],dx[1],axis=1)/(rad*st)+np.gradient(w[...,7],dx[2],axis=2)/(rad*st)
   core=np.s_[1:-1,1:-1,1:-1];divnorm+=np.sum(div[core]**2*vol[core]);divvol+=np.sum(vol[core]);divmax=max(divmax,float(np.max(abs(div[core]))))
   excess=np.maximum(w[...,0]-1,0)*vol*roi;mass+=excess.sum();centroid+=np.sum(pos*excess[...,None],axis=(0,1,2));pmin=min(pmin,float(pr.min()));rhomin=min(rhomin,float(w[...,0].min()))
  log=(d/'run.log').read_text();times=re.findall(r'Total timeloop took\s*:\s*([\d.]+)',log)
  table=np.loadtxt(d/('adv_yin.log' if yy else 'adv.log'),skiprows=1,usecols=(0,1,2))
  res=dict(case=d.name,snapshot=snap,time=h['time'],steps=h['it'],cells=int(np.prod(h['domain_nx']))*(2 if yy else 1),dt_initial=float(table[0,2]),loop_seconds=max(map(float,times)),wall_seconds=run['wall'],normalized_L1=l1/den1,normalized_L2=np.sqrt(l2/den2),roi_normalized_L1=roierr/roiden,density_Linf=linf,peak_fraction=peak/.5,sampled_exact_peak_fraction=exactpeak/.5,centroid_error=float(np.linalg.norm(centroid/mass-center)),velocity_rms=float(np.sqrt(velerr/volsum)),B_rms=float(np.sqrt(berr/volsum)),pressure_min=pmin,rho_min=rhomin,volume_bias=volsum/(28*np.pi/3)-1,divB_fd_rms=float(np.sqrt(divnorm/divvol)),divB_fd_max=divmax)
  results.append(res)
(root.parent/'results.json').write_text(json.dumps(results,indent=2)+'\n')
for r in results:
 if r['snapshot']==2:print(f"{r['case']:20s} steps={r['steps']:4d} time={r['loop_seconds']:7.2f}s L1={r['normalized_L1']:.4f} ROI={r['roi_normalized_L1']:.4f} peak={r['peak_fraction']:.3f} B={r['B_rms']:.2e}")

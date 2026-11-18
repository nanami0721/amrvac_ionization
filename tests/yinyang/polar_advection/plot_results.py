#!/usr/bin/env python3
"""Figures from completed benchmarks. Run only after timing jobs finish."""
from pathlib import Path
import sys,json,importlib.util
sys.dont_write_bytecode=True
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));import check_mhd
source,root,out=map(Path,sys.argv[1:4]);out.mkdir(exist_ok=True,parents=True);u=check_mhd.reader(source)
results=json.loads((root.parent/'results.json').read_text());rows=[r for r in results if r['snapshot']==2];by={r['case']:r for r in rows}
colors={'spherical':'#b04b39','yinyang':'#197c9b'}
fig,axes=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
for geom in colors:
 for tr,label in [(0,'polar'),(1,'equatorial')]:
  rr=[by[f'{geom}-s{s}-t{tr}'] for s in (1,2)];ls='-' if tr==0 else '--';marker='o' if tr==0 else 's'
  axes[0].plot([32,64],[r['steps'] for r in rr],ls,marker=marker,color=colors[geom],label=f'{geom} / {label}')
  axes[1].plot([32,64],[r['loop_seconds'] for r in rr],ls,marker=marker,color=colors[geom])
  axes[2].plot([r['loop_seconds'] for r in rr],[r['roi_normalized_L1'] for r in rr],ls,marker=marker,color=colors[geom])
if 'yinyang-s3-t0' in by:
 r=by['yinyang-s3-t0'];axes[2].scatter(r['loop_seconds'],r['roi_normalized_L1'],marker='*',s=130,color=colors['yinyang'],zorder=5);axes[2].annotate('YY follow-up',(r['loop_seconds'],r['roi_normalized_L1']),xytext=(8,7),textcoords='offset points',fontsize=8)
 axes[2].axhline(by['spherical-s2-t0']['roi_normalized_L1'],ls=':',color='#777777',lw=.8)
for ax in axes:ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
axes[0].set(xlabel='Equivalent full-sphere Nθ',ylabel='Time steps to t = 0.6',yscale='log',xticks=[32,64],title='CFL cost')
axes[1].set(xlabel='Equivalent full-sphere Nθ',ylabel='Evolution wall time [s]',yscale='log',xticks=[32,64],title='Four MPI ranks; matched cells')
axes[2].set(xlabel='Evolution wall time [s]',ylabel='Pulse-region relative L1 error',xscale='log',title='Measured cost–accuracy trade-off')
axes[0].legend(fontsize=8);fig.savefig(out/'cost_accuracy.png',dpi=180);fig.savefig(out/'cost_accuracy.pdf');plt.close(fig)

def rot(tr):
 return np.eye(3) if tr==0 else np.array([[0,0,1],[0,1,0],[-1,0,0.]])

def load_interp(case):
 yy=case.startswith('yinyang');data=[];bounds=[]
 for patch,suffix in enumerate(['_yin','_yang'] if yy else ['']):
  h,w=check_mhd.load(u,root/case/f'adv{suffix}0002.dat');dx=(h['xmax']-h['xmin'])/h['domain_nx'];axes=[a+(np.arange(n)+.5)*d for a,n,d in zip(h['xmin'],h['domain_nx'],dx)];rho=w[...,0]
  if not yy:
   # Across each true pole: reflect theta and shift longitude by pi.
   rho=np.concatenate([np.roll(rho[:,:1,:],len(axes[2])//2,axis=2),rho,np.roll(rho[:,-1:,:],len(axes[2])//2,axis=2)],axis=1)
   axes[1]=np.r_[-axes[1][0],axes[1],2*np.pi-axes[1][-1]]
   rho=np.concatenate([rho[:,:,-1:],rho,rho[:,:,:1]],axis=2);axes[2]=np.r_[axes[2][-1]-2*np.pi,axes[2],axes[2][0]+2*np.pi]
  data.append(RegularGridInterpolator(axes,rho,bounds_error=False,fill_value=None));bounds.append((axes[1][0],axes[1][-1],axes[2][0],axes[2][-1]))
 def sample(points):
  xyz=np.atleast_2d(points);coords=[];score=[]
  for patch in range(len(data)):
   local=xyz.copy()
   if patch:local=local[:,[0,2,1]];local[:,0]*=-1
   r=np.linalg.norm(local,axis=1);t=np.arctan2(np.linalg.norm(local[:,:2],axis=1),local[:,2]);p=np.arctan2(local[:,1],local[:,0])
   if not yy:p=p%(2*np.pi)
   coords.append(np.c_[r,t,p]);a,b,c,d=bounds[patch];score.append(np.minimum.reduce([t-a,b-t,p-c,d-p]))
  choice=np.argmax(score,axis=0);v=np.empty(len(xyz))
  for patch in range(len(data)):
   mask=choice==patch
   if mask.any():v[mask]=data[patch](coords[patch][mask])
  return v
 return sample

xx=np.linspace(-.15,.75,260);zz=np.linspace(1.05,1.95,260);X,Z=np.meshgrid(xx,zz);canonical=np.c_[X.ravel(),np.zeros(X.size),Z.ravel()];exact=.5*np.exp(-((X-.3)**2+(Z-1.5)**2)/.18**2)
fig,axes=plt.subplots(2,3,figsize=(10.2,6.8),layout='constrained',sharex=True,sharey=True)
for tr in [0,1]:
 for col,geom in enumerate(['spherical','yinyang','analytic']):
  value=exact if col==2 else load_interp(f'{geom}-s2-t{tr}')(canonical@rot(tr).T).reshape(X.shape)-1
  valid=(np.linalg.norm(canonical,axis=1)<=2).reshape(X.shape);value=np.where(valid,value,np.nan)
  ax=axes[tr,col];im=ax.pcolormesh(X,Z,value,vmin=0,vmax=.5,cmap='viridis',shading='auto',rasterized=True);ax.contour(X,Z,exact,levels=[.1,.3],colors='white',linewidths=.6,linestyles='--');ax.set_aspect('equal');ax.set_title(geom)
  if col==0:ax.set_ylabel(('Polar path' if tr==0 else 'Equatorial path')+'\ntransverse coordinate')
  if tr==1:ax.set_xlabel('Along-trajectory coordinate')
fig.colorbar(im,ax=axes,label='Density excess ρ − 1',shrink=.8);fig.suptitle('Same physical plane after undoing rotation; fine grid, t = 0.6\nWhite contours: exact solution')
fig.savefig(out/'density_comparison.png',dpi=180);fig.savefig(out/'density_comparison.pdf');plt.close(fig)

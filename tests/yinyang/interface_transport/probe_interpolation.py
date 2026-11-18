#!/usr/bin/env python3
"""Static analytic-field probe only: no MHD update or positivity guarantee."""
from pathlib import Path
import sys,json,argparse
sys.dont_write_bytecode=True
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));import check_mhd
from diagnostics import geometry,basis,exact,local_coords
from analyze_native import fv_from_padded
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('runs',type=Path);p.add_argument('output',type=Path);a=p.parse_args();u=check_mhd.reader(a.source);rows=[]
for scale in [1,2,3]:
 d=a.runs/f'c1-p1-d1-s{scale}';geoms=[]
 for patch,name in enumerate(['yin','yang']):
  h,_=check_mhd.load(u,d/f'adv_{name}0000.dat');g=geometry(h,patch);g['patch']=patch;geoms.append(g)
 for direction in [1,-1]:
  for time in [0.,.3,.6]:
   for patch,g in enumerate(geoms):
    other=geoms[1-patch];cart=exact(other['pos'],time,1,direction,1)[1]
    axes=[np.r_[x[0]-dx,x,x[-1]+dx] for x,dx in zip(g['axes'],g['dx'])]
    r,t,phi=np.meshgrid(*axes,indexing='ij');bb=basis(t,phi,patch);pos=r[...,None]*bb[..., :,0];ex=exact(pos,time,1,direction,1)[1];pad=np.einsum('...ji,...j->...i',bb,ex)
    ideal,_=fv_from_padded(g,pad)
    mask=np.zeros(r.shape,bool);mask[1:-1,[0,-1],1:-1]=True;mask[1:-1,1:-1,[0,-1]]=True
    coords=local_coords(pos[mask],1-patch);q=(coords-np.array([x[0] for x in other['axes']]))/other['dx'];ir=np.rint(q[:,0]).astype(int)
    assert np.max(abs(ir-q[:,0]))<1e-10
    core=np.zeros(g['r'].shape,bool);core[:,[0,-1],:]=True;core[:,:,[0,-1]]=True;core[0]=False;core[-1]=False
    for degree in [1,2]:
     ids=[];weights=[]
     for axis in [1,2]:
      if degree==1:
       left=np.clip(np.floor(q[:,axis]).astype(int),0,len(other['axes'][axis])-2);ii=left[:,None]+np.arange(2)
      else:
       center=np.clip(np.rint(q[:,axis]).astype(int),1,len(other['axes'][axis])-2);ii=center[:,None]+np.arange(-1,2)
      ww=np.ones(ii.shape)
      for j in range(degree+1):
       for k in range(degree+1):
        if j!=k:ww[:,j]*=(q[:,axis]-ii[:,k])/(ii[:,j]-ii[:,k])
      ids.append(ii);weights.append(ww)
     val=np.zeros((len(ir),3))
     for j in range(degree+1):
      for k in range(degree+1):val+=weights[0][:,j,None]*weights[1][:,k,None]*cart[ir,ids[0][:,j],ids[1][:,k]]
     trial=pad.copy();trial[mask]=np.einsum('...ji,...j->...i',bb[mask],val);div,_=fv_from_padded(g,trial);err=np.linalg.norm(val-ex[mask],axis=-1);v=g['vol'][core];delta=(div-ideal)[core]
     rows.append(dict(scale=scale,direction=direction,time=time,patch=patch,degree=degree,ghost_B_max=float(err.max()),ghost_B_rms=float(np.sqrt(np.mean(err**2))),induced_div_max=float(abs(delta).max()),induced_div_rms=float(np.sqrt(np.sum(delta*delta*v)/v.sum()))))
a.output.write_text(json.dumps(rows,indent=2)+'\n')
print('Static probe records:',len(rows))

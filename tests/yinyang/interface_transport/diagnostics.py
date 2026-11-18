"""Geometry, exact fields and seam-inclusive diagnostic divergence (not solver divB)."""
import numpy as np
from scipy.interpolate import RegularGridInterpolator

def rotation(placement):
 p=np.deg2rad(146.25 if placement else 0.)
 return np.array([[-np.sin(p),0,np.cos(p)],[np.cos(p),0,np.sin(p)],[0,1,0.]])

def basis(t,p,patch):
 st,ct,sp,cp=np.sin(t),np.cos(t),np.sin(p),np.cos(p)
 b=np.stack([np.stack([st*cp,st*sp,ct],-1),np.stack([ct*cp,ct*sp,-st],-1),np.stack([-sp,cp,np.zeros_like(p)],-1)],-1)
 if patch:b=b[...,[0,2,1],:];b[...,0,:]*=-1
 return b

def geometry(h,patch):
 dx=(h['xmax']-h['xmin'])/h['domain_nx'];axes=[a+(np.arange(n)+.5)*d for a,n,d in zip(h['xmin'],h['domain_nx'],dx)]
 r,t,p=np.meshgrid(*axes,indexing='ij');b=basis(t,p,patch);pos=r[...,None]*b[...,:,0]
 solid=np.zeros(t.shape[1:]);ns=16
 for jt in range(ns):
  tt=axes[1][:,None]+((jt+.5)/ns-.5)*dx[1];ds=np.cos(tt-dx[1]/(2*ns))-np.cos(tt+dx[1]/(2*ns))
  for jp in range(ns):
   pp=axes[2][None,:]+((jp+.5)/ns-.5)*dx[2]
   ot=np.arctan2(np.hypot(np.sin(tt)*np.cos(pp),np.cos(tt)),np.sin(tt)*np.sin(pp));op=np.arctan2(np.cos(tt),-np.sin(tt)*np.cos(pp))
   overlap=(ot>=h['xmin'][1])&(ot<=h['xmax'][1])&(op>=h['xmin'][2])&(op<=h['xmax'][2])
   solid+=ds*dx[2]/ns*np.where(overlap,.5,1.)
 vol=((r+dx[0]/2)**3-(r-dx[0]/2)**3)/3*solid[None]
 # Fixed angular-width physical bands on BOTH sides of BOTH patch edges.
 # Selecting only native exterior rows would bias inward vs outward transport.
 seam=np.zeros(r.shape,bool);width=np.deg2rad(10)
 for frame in [0,1]:
  coord=local_coords(pos,frame);theta,phi=coord[...,1],coord[...,2];lo,hi=h['xmin'],h['xmax']
  tb=(np.minimum(abs(theta-lo[1]),abs(theta-hi[1]))<width)&(phi>=lo[2]-width)&(phi<=hi[2]+width)
  pb=(np.minimum(abs(phi-lo[2]),abs(phi-hi[2]))<width)&(theta>=lo[1]-width)&(theta<=hi[1]+width)
  seam|=tb|pb
 return dict(axes=axes,dx=dx,r=r,t=t,p=p,basis=b,pos=pos,vol=vol,seam=seam)

def exact(pos,time,case,direction,placement):
 rot=rotation(placement);vel=direction*rot[:,0];center=rot@np.array([-.3*direction,0,1.5])+vel*time;off=pos-center
 if case==0:
  rho=1+.5*np.exp(-np.sum(off*off,-1)/.18**2);field=np.broadcast_to(rot@np.array([.2,.1,.3]),pos.shape)
 else:
  rho=np.ones(pos.shape[:-1]);phase=np.pi/2*np.exp(-(np.sum(off*rot[:,0],-1)/.18)**2)
  local=np.stack([np.zeros_like(phase),.3*np.cos(phase),.3*np.sin(phase)],-1);field=local@rot.T
 return rho,field,vel,center

def local_coords(pos,patch):
 q=pos.copy()
 if patch:q=q[...,[0,2,1]];q[...,0]*=-1
 return np.stack([np.linalg.norm(q,axis=-1),np.arctan2(np.linalg.norm(q[...,:2],axis=-1),q[...,2]),np.arctan2(q[...,1],q[...,0])],-1)

def divergence(g,localB,other_interp,patch):
 # Reconstruct one angular ghost layer from the other patch, in Cartesian B.
 # Radial derivatives omit the first/last radial layer, so no radial extrapolation.
 axes=g['axes'];dx=g['dx'];nt,np_=len(axes[1]),len(axes[2]);t=np.r_[axes[1][0]-dx[1],axes[1],axes[1][-1]+dx[1]];p=np.r_[axes[2][0]-dx[2],axes[2],axes[2][-1]+dx[2]]
 r,tt,pp=np.meshgrid(axes[0],t,p,indexing='ij');pad=np.empty((*r.shape,3));pad[:,1:-1,1:-1]=localB
 edge=np.ones(r.shape,bool);edge[:,1:-1,1:-1]=False
 b=basis(tt[edge],pp[edge],patch);pos=r[edge,None]*b[:,:,0];coords=local_coords(pos,1-patch)
 # Donors must stay inside angular interpolation range; radial values coincide.
 for axis in [1,2]:
  assert np.all(coords[:,axis]>=other_interp.grid[axis][0]-1e-12)
  assert np.all(coords[:,axis]<=other_interp.grid[axis][-1]+1e-12)
 cart=other_interp(coords);pad[edge]=np.einsum('...ji,...j->...i',b,cart)
 divr=np.gradient(g['r']**2*localB[...,0],dx[0],axis=0)/g['r']**2
 st=np.sin(tt)
 divt=(st[:,2:,1:-1]*pad[:,2:,1:-1,1]-st[:,:-2,1:-1]*pad[:,:-2,1:-1,1])/(2*dx[1]*g['r']*np.sin(g['t']))
 divp=(pad[:,1:-1,2:,2]-pad[:,1:-1,:-2,2])/(2*dx[2]*g['r']*np.sin(g['t']))
 return divr+divt+divp

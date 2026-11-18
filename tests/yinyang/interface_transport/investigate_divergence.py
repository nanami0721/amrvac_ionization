from pathlib import Path
import sys,json,argparse
sys.dont_write_bytecode=True
p=argparse.ArgumentParser(description="Localize divergence and isolate analytic ghost interpolation error")
p.add_argument("source",type=Path);p.add_argument("runs",type=Path);p.add_argument("output",type=Path);p.add_argument("--scales",type=int,nargs="+",default=[1,2,3]);a=p.parse_args()
import numpy as np
from scipy.interpolate import RegularGridInterpolator
src=a.source.resolve()
if a.output.resolve().is_relative_to(src):p.error('output must be external to source')
sys.path.insert(0,str(src/'tests/yinyang'));import check_mhd
sys.path.insert(0,str(src/'tests/yinyang/interface_transport'))
from diagnostics import geometry,exact,divergence,basis,local_coords
u=check_mhd.reader(src); rows=[]
for scale in a.scales:
 for snap in [0,2]:
  root=a.runs/f'c1-p1-d1-s{scale}'
  data=[]
  for patch,name in enumerate(['yin','yang']):
   h,w=check_mhd.load(u,root/f'adv_{name}{snap:04d}.dat');g=geometry(h,patch);e=exact(g['pos'],h['time'],1,1,1)[1];b=np.einsum('...ij,...j->...i',g['basis'],w[...,5:8]);data.append((h,w,g,e,b))
  for patch,(h,w,g,e,b) in enumerate(data):
   other=data[1-patch]
   num=divergence(g,w[...,5:8],RegularGridInterpolator(other[2]['axes'],other[4],bounds_error=False,fill_value=None),patch)
   el=np.einsum('...ji,...j->...i',g['basis'],e)
   dex=divergence(g,el,RegularGridInterpolator(other[2]['axes'],other[3],bounds_error=False,fill_value=None),patch)
   class Oracle:
    grid=tuple(np.array([-np.inf,np.inf]) for _ in range(3))
    def __call__(self,coords):
     bb=basis(coords[:,1],coords[:,2],1-patch);pos=coords[:,0,None]*bb[:,:,0]
     return exact(pos,h['time'],1,1,1)[1]
   oracle=divergence(g,el,Oracle(),patch)
   # Ghost interpolation alone: omit unused double-angular corners.
   axes=g['axes'];dx=g['dx']
   aa=[axes[0],np.r_[axes[1][0]-dx[1],axes[1],axes[1][-1]+dx[1]],np.r_[axes[2][0]-dx[2],axes[2],axes[2][-1]+dx[2]]]
   rr,tt,pp=np.meshgrid(*aa,indexing='ij');face=np.zeros(rr.shape,bool)
   face[:,[0,-1],1:-1]=True;face[:,1:-1,[0,-1]]=True
   bb=basis(tt[face],pp[face],patch);pos=rr[face,None]*bb[:,:,0];coords=local_coords(pos,1-patch)
   sampled=RegularGridInterpolator(other[2]['axes'],other[3],bounds_error=False,fill_value=None)(coords)
   ghost_error=np.linalg.norm(sampled-exact(pos,h['time'],1,1,1)[1],axis=-1)
   core=np.zeros(num.shape,bool);core[1:-1]=True
   edge=np.zeros(num.shape,bool);edge[:,[0,-1],:]=True;edge[:,:,[0,-1]]=True
   for region,mask in [('edge',core&edge),('interior',core&~edge),('seam',core&g['seam']),('original_shell_seam',core&g['seam']&(g['r']>1+g['dx'][0])&(g['r']<2-g['dx'][0])),('bulk',core&(g['r']>1.2)&(g['r']<1.8)),('bulk_seam',core&g['seam']&(g['r']>1.2)&(g['r']<1.8))]:
    vol=g['vol'][mask];idx=np.unravel_index(np.argmax(np.where(mask,abs(num),-1)),num.shape)
    rec=dict(scale=scale,time=h['time'],patch=patch,region=region,index=list(map(int,idx)),position=g['pos'][idx].tolist(),r=float(g['r'][idx]),ghost_B_error_max=float(ghost_error.max()),ghost_B_error_rms=float(np.sqrt(np.mean(ghost_error**2))))
    for key,v in [('num',num),('sample',dex),('oracle',oracle),('interp',dex-oracle),('defect',num-dex)]:rec[key+'_max']=float(abs(v[mask]).max());rec[key+'_rms']=float(np.sqrt(np.sum(v[mask]**2*vol)/sum(vol)))
    rows.append(rec)
a.output.write_text(json.dumps(rows,indent=2))
for r in rows:
 if r['region']=='edge':print(r['scale'],r['time'],r['patch'],'num/sample/oracle/interp max',*[round(r[k+'_max'],5) for k in ['num','sample','oracle','interp']],r['index'],round(r['r'],3))

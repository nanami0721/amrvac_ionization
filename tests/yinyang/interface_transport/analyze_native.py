#!/usr/bin/env python3
"""Read solver get_divb dumps; independently check its exact-field FV reference."""
from pathlib import Path
import argparse,json,sys
sys.dont_write_bytecode=True
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));import check_mhd
from diagnostics import geometry,exact,basis

def exact_fv(g,time,info):
 axes=[np.r_[a[0]-d,a,a[-1]+d] for a,d in zip(g['axes'],g['dx'])]
 r,t,p=np.meshgrid(*axes,indexing='ij');bb=basis(t,p,g['patch']);pos=r[...,None]*bb[..., :,0]
 field=exact(pos,time,info['case'],info['direction'],info['placement'])[1]
 b=np.einsum('...ji,...j->...i',bb,field)
 return fv_from_padded(g,b)

def fv_from_padded(g,b):
 c=b[1:-1,1:-1,1:-1];r=g['r'];t=g['t'];dr,dt,dp=g['dx'];solid=2*np.sin(t)*np.sin(dt/2)*dp
 vol=(r*r+dr*dr/12)*dr*solid
 div=((r+dr/2)**2*solid*(c[...,0]+b[2:,1:-1,1:-1,0])-(r-dr/2)**2*solid*(c[...,0]+b[:-2,1:-1,1:-1,0]))/2
 div+=r*dr*dp*(np.sin(t+dt/2)*(c[...,1]+b[1:-1,2:,1:-1,1])-np.sin(t-dt/2)*(c[...,1]+b[1:-1,:-2,1:-1,1]))/2
 div+=r*dr*dt*(b[1:-1,1:-1,2:,2]-b[1:-1,1:-1,:-2,2])/2
 return div/vol,vol

def main():
 p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('runs',type=Path);a=p.parse_args();u=check_mhd.reader(a.source);rows=[]
 for d in sorted(a.runs.glob('c*-p*-d*-s*')):
  info=json.loads((d/'run.json').read_text())
  if info['debug']:continue
  for snap in [0,1,2]:
   patches=[];times=[];ch=[];roundtrip=0.
   for patch,name in enumerate(['yin','yang']):
    h,w=check_mhd.load(u,d/f'adv_{name}{snap:04d}.dat');g=geometry(h,patch);g['patch']=patch
    data=np.full((*tuple(h['domain_nx']),10),np.nan);coverage=np.zeros(h['domain_nx'],dtype=np.int32)
    for f in d.glob(f'native_p{patch}_r*_i{h["it"]:07d}.bin'):
     with f.open('rb') as stream:
      head=np.fromfile(stream,dtype=np.int32,count=5);stamp=np.fromfile(stream,dtype=np.float64,count=2);values=np.fromfile(stream,dtype=np.float64)
     assert tuple(head[:2])==(1,patch) and head[3]==h['it']
     assert len(values)==10*head[4] and abs(stamp[0]-h['time'])<1e-12
     v=values.reshape(-1,10);index=np.rint((v[:,:3]-h['xmin'])/g['dx']-.5).astype(int)
     assert np.all(index>=0) and np.all(index<h['domain_nx'])
     assert np.max(abs((index+.5)*g['dx']+h['xmin']-v[:,:3]))<1e-12
     ix=tuple(index.T);data[ix]=v;np.add.at(coverage,ix,1);times.append(stamp[0]);ch.append(stamp[1])
    assert np.all(coverage==1) and np.all(np.isfinite(data)),str(d)
    ref,vol=exact_fv(g,h['time'],info);roundtrip=max(roundtrip,float(np.max(abs(ref-data[...,4]))))
    assert np.allclose(vol,data[...,5],rtol=2e-13,atol=1e-16)
    patches.append((g,data))
   assert roundtrip<2e-12,roundtrip
   assert max(ch)==min(ch)
   for region in ['all','seam','core','core_seam','bulk_seam','bulk','bulk_patch_edge','bulk_off_edge','inner','outer']:
    sq=exsq=dsq=volume=psi2=incoming2=outgoing2=0.;mx=-1.;peak=None
    for patch,(g,v) in enumerate(patches):
     rad=g['r'];mask=np.ones(rad.shape,bool)
     if 'core' in region:mask[0]=False;mask[-1]=False
     if 'seam' in region:mask &= g['seam']
     if region.startswith('bulk'):mask &= (rad>1.2)&(rad<1.8)
     if region in ['bulk_patch_edge','bulk_off_edge']:
      edge=np.zeros(rad.shape,bool);edge[:,[0,-1],:]=True;edge[:,:,[0,-1]]=True
      mask &= edge if region=='bulk_patch_edge' else ~edge
     if region=='inner':mask &= rad<g['axes'][0][0]-g['dx'][0]/2+.1
     if region=='outer':mask &= rad>g['axes'][0][-1]+g['dx'][0]/2-.1
     vol=g['vol'][mask];num=v[...,3][mask];ex=v[...,4][mask];psi=v[...,6][mask];db=v[...,7][mask]
     sign=-1 if region=='inner' else 1
     sq+=float(np.sum(num*num*vol));exsq+=float(np.sum(ex*ex*vol));dsq+=float(np.sum((num-ex)**2*vol));volume+=float(vol.sum());psi2+=float(np.sum(psi*psi*vol))
     incoming2+=float(np.sum((psi-sign*ch[0]*db)**2*vol));outgoing2+=float(np.sum((psi+sign*ch[0]*db)**2*vol))
     idx=np.unravel_index(np.argmax(np.where(mask,abs(v[...,3]),-1)),rad.shape)
     if abs(v[...,3][idx])>mx:mx=float(abs(v[...,3][idx]));peak=dict(patch=patch,r=float(rad[idx]),position=g['pos'][idx].tolist())
    rows.append(dict(name=d.name,**info,time=h['time'],snapshot=snap,region=region,div_rms=np.sqrt(sq/volume),div_max=mx,exact_rms=np.sqrt(exsq/volume),defect_rms=np.sqrt(dsq/volume),psi_rms=np.sqrt(psi2/volume),incoming_rms=np.sqrt(incoming2/volume) if region in ['inner','outer'] else None,outgoing_rms=np.sqrt(outgoing2/volume) if region in ['inner','outer'] else None,peak=peak,reference_operator_max_difference=roundtrip))
 (a.runs.parent/'native_results.json').write_text(json.dumps(rows,indent=2)+'\n')
 for x in rows:
  if x['snapshot']==2 and x['region'] in ['all','core_seam','outer']:print(x['name'],x['region'],'RMS',x['div_rms'],'max',x['div_max'])
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Native/YY controls, matching actual grid points and separating gain/shape error."""
import argparse,json,re,math
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('workdir',type=Path);p.add_argument('--complete',action='store_true');p.add_argument('--plot',action='store_true');a=p.parse_args();root=a.workdir
pi=np.pi

def config(d):
 text=(d/'amrvac.par').read_text()
 def val(key):
  x=re.search(r'(?m)^\s*'+key+r'\s*=\s*([^\s,\n/]+)',text).group(1)
  return float(x.lower().replace('d','e'))
 nn=np.array([int(val('domain_nx'+str(i))) for i in [1,2,3]])
 lower=np.array([val('xprobmin'+str(i)) for i in [1,2,3]])*[1,2*pi,2*pi]
 upper=np.array([val('xprobmax'+str(i)) for i in [1,2,3]])*[1,2*pi,2*pi]
 return dict(nn=nn,lower=lower,dx=(upper-lower)/nn,direction=int(val('direction')),dt=val('dtpar'),native='native_control=T' in text,patch=int(val('reference_patch')))

def snapshots(d,time):
 c=config(d);prefix='transport' if abs(time-.6)<1e-10 else 'sample-'+str(round(time/c['dt']))
 files=sorted(d.glob(prefix+'-rank-*.txt'))
 if not files:return None
 # Fortran ES25.16 omits E for three-digit exponents (e.g. 3.29-100).
 def number(token):
  if isinstance(token,bytes):token=token.decode()
  return float(re.sub(r'(?<=\d)([+-]\d{3})$',r'E\1',token))
 return np.concatenate([np.loadtxt(f,ndmin=2,converters=number) for f in files])

def geometry(data,c):
 patch=data[:,0].astype(int);idx=data[:,1:4].astype(int)
 x=c['lower']+(idx+.5)*c['dx'];r,theta,phi=x.T
 st,ct,sp,cp=np.sin(theta),np.cos(theta),np.sin(phi),np.cos(phi)
 frame=np.empty((len(x),3,3))
 frame[:,:,0]=np.array([st*cp,st*sp,ct]).T
 frame[:,:,1]=np.array([ct*cp,ct*sp,-st]).T
 frame[:,:,2]=np.array([-sp,cp,np.zeros(len(x))]).T
 sel=patch==1
 frame[sel,0,:]*=-1
 temp=frame[sel,1,:].copy();frame[sel,1,:]=frame[sel,2,:];frame[sel,2,:]=temp
 pos=r[:,None]*frame[:,:,0]
 # Analytic FV volume, matching native spherical cell volumes.
 dr,dt,dp=c['dx'];volume=((r+dr/2)**3-(r-dr/2)**3)/3*(np.cos(theta-dt/2)-np.cos(theta+dt/2))*dp
 return idx,pos,frame,volume

def measure(data,c,time,weight=None):
 idx,pos,frame,volume=geometry(data,c)
 if weight is not None:volume*=weight
 ang=math.radians(146.25);normal=np.array([math.cos(ang),math.sin(ang),0]);vel=c['direction']*np.array([-normal[1],normal[0],0])
 q=pos-1.5*normal-(time-.3)*vel;R=.65;a0=.001*R/((8/math.sqrt(7))*(6/7)**3)
 exact=-8*a0/R**2*np.maximum(0,1-np.sum(q*q,axis=1)/R**2)[:,None]**3*np.cross(q,normal)
 cart=np.einsum('nij,nj->ni',frame,data[:,9:12]);dv=np.einsum('nij,nj->ni',frame,data[:,5:8]/data[:,4,None])-vel
 ref=np.sum(exact*exact,axis=1);norm=np.sum(volume*ref)
 error_density=np.sum((cart-exact)**2,axis=1)
 error=np.sum(volume*error_density);gain=np.sum(volume*np.sum(cart*exact,axis=1))/norm
 shape=np.sum(volume*np.sum((cart-gain*exact)**2,axis=1))/norm
 peak=int(np.argmax(np.sum(cart*cart,axis=1)))
 return dict(relative_l2=float(np.sqrt(error/norm)),gain=float(gain),shape_error=float(np.sqrt(shape)),energy_ratio=float(np.sum(volume*np.sum(cart*cart,axis=1))/norm),reference_weighted_velocity_error=float(np.sqrt(np.sum(volume*ref*np.sum(dv*dv,axis=1))/norm)),peak_fraction=float(np.linalg.norm(cart[peak])/.001),peak_position=pos[peak].tolist(),peak_patch=int(data[peak,0]),peak_reference_fraction=float(np.linalg.norm(exact[peak])/.001),peak_vector_error_over_analytic_peak=float(np.sqrt(error_density[peak])/.001)),volume,ref,error_density

def yy_weights(data,c):
 _,pos,_,_=geometry(data,c);other=pos.copy();sel=data[:,0]==0
 other[sel]=np.column_stack([-pos[sel,0],pos[sel,2],pos[sel,1]])
 theta=np.arccos(other[:,2]/np.linalg.norm(other,axis=1));phi=np.arctan2(other[:,1],other[:,0])
 inside=(theta>=.109375*2*pi)&(theta<=.390625*2*pi)&(np.abs(phi)<=.390625*2*pi)
 return np.where(inside,.5,1.)

results={}
for d in sorted(root.glob('*-p*-s*-d*-*')):
 f=d/'transport-history.txt'
 if not f.exists():continue
 h=np.loadtxt(f,ndmin=2)
 if abs(h[-1,1]-.6)>1e-10:continue
 c=config(d);data=snapshots(d,.6)
 v,_,_,_=measure(data,c,.6,None if c['native'] else yy_weights(data,c))
 # Independent norms from dumped states must reproduce the in-solver histories.
 assert abs(v['relative_l2']-math.sqrt(h[-1,15]/h[-1,16]))<1e-10,d.name
 v.update(initial_l2=float(np.sqrt(h[0,15]/h[0,16])),active_divb=float(h[:,2].max()),ghost_divb=float(h[:,4].max()),rho_min=float(h[:,13].min()),p_min=float(h[:,14].min()),velocity_max=float(h[:,7].max()))
 results[d.name]=v
common={}
for direction in [1,-1]:
 yd=root/f'yy-p0-s2-d{direction}-base';nd=[root/f'native-p{i}-s2-d{direction}-base' for i in [0,1]]
 if not all(x.name in results for x in [yd,*nd]):continue
 cy=config(yd)
 for time in [0.,.404,.522,.6]:
  yy=snapshots(yd,time)
  dense=np.full((2,*cy['nn'],8),np.nan)
  for d in nd:
   c=config(d);native=snapshots(d,time);offset=(c['lower']-cy['lower'])/cy['dx']
   assert np.max(np.abs(c['dx']/cy['dx']-1))<1e-13
   assert np.max(np.abs(offset-np.rint(offset)))<1e-12
   idx=native[:,1:4].astype(int)+np.rint(offset).astype(int);valid=np.all((idx>=0)&(idx<cy['nn']),axis=1);ix=idx[valid]
   dense[(native[valid,0].astype(int),ix[:,0],ix[:,1],ix[:,2])]=native[valid,4:12]
  idx=yy[:,:4].astype(int);matched=dense[tuple(idx.T)];mask=np.isfinite(matched[:,0])
  same=yy[mask].copy();uncoupled=same.copy();uncoupled[:,4:12]=matched[mask]
  weight=yy_weights(same,cy);ym,vol,ref,ye=measure(same,cy,time,weight);nm,_,_,ne=measure(uncoupled,cy,time,weight)
  difference=float(np.sqrt(np.sum(vol*np.sum((same[:,9:12]-uncoupled[:,9:12])**2,axis=1))/np.sum(vol*ref)))
  fullvol=geometry(yy,cy)[3]*yy_weights(yy,cy);e=np.sum(yy[:,9:12]**2,axis=1)
  outside=float(np.sum(fullvol[~mask]*e[~mask])/np.sum(fullvol*e))
  ip=int(np.argmax(np.sum(same[:,9:12]**2,axis=1)));cell=same[ip,1:4].astype(int)
  at_same_point=measure(uncoupled[ip:ip+1],cy,time)[0]
  peak_detail=dict(yy_point_vector_error_over_peak=ym['peak_vector_error_over_analytic_peak'],native_same_point_vector_error_over_peak=at_same_point['peak_vector_error_over_analytic_peak'],yy_cell=cell.tolist(),yy_patch=int(same[ip,0]),distance_to_angular_edge_cells=int(min(cell[1],cell[2],cy['nn'][1]-1-cell[1],cy['nn'][2]-1-cell[2])),native_at_same_point_fraction=float(np.linalg.norm(uncoupled[ip,9:12])/.001),delta_local_B_over_peak=((same[ip,9:12]-uncoupled[ip,9:12])/.001).tolist())
  common[f'd{direction}-t{time}']=dict(error_squared_increase=float(np.sum(vol*np.maximum(ye-ne,0))/np.sum(vol*ref)),error_squared_decrease=float(np.sum(vol*np.maximum(ne-ye,0))/np.sum(vol*ref)),peak_detail=peak_detail,yy=ym,uncoupled_native=nm,direct_difference_l2=difference,yy_energy_outside_common_fraction=outside,matched_cells=int(mask.sum()))
summary=dict(cases=results,common_grid=common)
(root/'controls-analysis.json').write_text(json.dumps(summary,indent=2))
for k,v in results.items():print(k,'L2',round(v['relative_l2'],5),'gain',round(v['gain'],5),'shape',round(v['shape_error'],5),'flow',format(v['reference_weighted_velocity_error'],'.3e'))
for k,v in common.items():print('common',k,'YY/native',round(v['yy']['relative_l2'],5),round(v['uncoupled_native']['relative_l2'],5),'difference',round(v['direct_difference_l2'],5),'outside',format(v['yy_energy_outside_common_fraction'],'.2e'))
if a.complete:
 assert len(results)==17
 assert all(v['active_divb']<1e-10 and v['rho_min']>0 and v['p_min']>0 for v in results.values())
 assert all((root/(k+'-memory/command.err')).stat().st_size==0 for k in results)
 assert len(common)==8
 print('Completed controls and independently reconstructed norms verified.')

if a.plot:
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
 variants=['base','mc','hll','hllmc'];labels=['TVDLF + VL','TVDLF + MC','HLL + VL','HLL + MC']
 x=np.arange(4)
 for shift,prefix,label in [(-.18,'native-p1-s2-d1-','Native Yang orientation'),(.18,'yy-p0-s2-d1-','Yin-Yang')]:
  axes[0].bar(x+shift,[100*results[prefix+v]['relative_l2'] for v in variants],width=.36,label=label)
 axes[0].set(xticks=x,xticklabels=labels,ylabel='Final relative L2 error (%)',title='Same grid, different spatial methods')
 axes[0].legend(fontsize=8);axes[0].grid(axis='y',alpha=.2)
 for direction,style in [(1,'-'),(-1,'--')]:
  for prefix,color,label in [('yy-p0','tab:blue','Yin-Yang'),('native-p1','tab:orange','Native Yang orientation')]:
   h=np.loadtxt(root/f'{prefix}-s2-d{direction}-base'/'transport-history.txt')
   axes[1].plot(h[:,1],h[:,11]/.001,style,color=color,label=label+(' forward' if direction==1 else ' reverse'))
 axes[1].set(xlabel='Time',ylabel='Global peak / analytic peak',title='Baseline peak histories')
 axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
 fig.savefig(root/'controls-summary.png',dpi=160)

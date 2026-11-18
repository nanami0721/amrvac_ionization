#!/usr/bin/env python3
"""Equilibrium errors and discrete budgets; process one snapshot at a time."""
from pathlib import Path
import argparse,json,sys,os,subprocess
sys.dont_write_bytecode=True
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('run',type=Path);a=p.parse_args()
if os.environ.get('AMRVAC_MEMORY_GUARD_ACTIVE')!='1':
 guard=Path(__file__).resolve().parents[1]/'interface_transport/convergence/memory_guard.py'
 out=a.run.resolve().with_name(a.run.name+'-analysis-memory')
 raise SystemExit(subprocess.run([sys.executable,str(guard),'--out',str(out),'--',sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]]).returncode)
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));import check_mhd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'interface_transport'));from diagnostics import basis,local_coords
u=check_mhd.reader(a.source)
info=json.loads((a.run/'run.json').read_text());assert info['returncode']==0 and info['stderr']==0
isyy=info['grid']=='yinyang';prefixes=['eq_yin','eq_yang'] if isyy else ['eq'];snapshots=sorted(a.run.glob(prefixes[0]+'[0-9][0-9][0-9][0-9].dat'));rows=[]
axis=np.array([np.sin(np.deg2rad(20)),0,np.cos(np.deg2rad(20))]);regions=['global','common','common_seam','common_off_seam','inner','outer']
for f in snapshots:
 number=f.stem[-4:];acc={key:np.zeros(19) for key in regions};maxima={key:np.zeros(5) for key in regions};rho_min=pmin=beta_min=float('inf');t=None
 for patch,prefix in enumerate(prefixes):
  h,w=check_mhd.load(u,a.run/f'{prefix}{number}.dat');t=h['time'];dx=(h['xmax']-h['xmin'])/h['domain_nx'];axes=[lo+(np.arange(n)+.5)*dd for lo,n,dd in zip(h['xmin'],h['domain_nx'],dx)]
  r,theta,phi=np.meshgrid(*axes,indexing='ij');bb=basis(theta,phi,patch);pos=r[...,None]*bb[..., :,0];normal=pos/r[...,None]
  rhoex=np.ones(r.shape) if info['case'] in [1,3] else np.exp(2*(1/r-1));pex=.5*rhoex
  Bex=np.zeros_like(pos) if info['case']==0 else .5/np.sqrt(info.get('beta',4.))*(3*np.sum(normal*axis,axis=-1)[...,None]*normal-axis)/r[...,None]**3
  if info['case']==3:Bex=np.broadcast_to(axis/np.sqrt(info.get('beta',4.)),pos.shape)
  local=np.einsum('...ji,...j->...i',bb,Bex);rho=w[...,0];vel2=np.sum(w[...,1:4]**2,axis=-1)/rho**2;pr=(h['gamma']-1)*(w[...,4]-.5*rho*vel2-.5*np.sum(w[...,5:8]**2,axis=-1))
  assert rho.min()>0 and pr.min()>0
  rho_min=min(rho_min,float(rho.min()));pmin=min(pmin,float(pr.min()))
  Berr=np.linalg.norm(w[...,5:8]-local,axis=-1);Bn=np.linalg.norm(local,axis=-1);mach=np.sqrt(vel2/(h['gamma']*.5))
  vA2=np.sum(local**2,axis=-1)/rhoex
  alfven_mach2=np.divide(vel2,vA2,out=np.zeros_like(vel2),where=vA2>0)
  Bnum2=np.sum(w[...,5:8]**2,axis=-1)
  if np.any(Bnum2>0):beta_min=min(beta_min,float(np.min(2*pr[Bnum2>0]/Bnum2[Bnum2>0])))
  data=np.full((*h['domain_nx'],9),np.nan);coverage=np.zeros(h['domain_nx'],np.int8)
  files=list(a.run.glob(f'equilibrium_p{patch}_r*_i{h["it"]:07d}.bin'))
  if not files:
   assert info['restart'], f'Missing native diagnostics at time {t}, patch {patch}'
   continue # Only earlier copied restart snapshots may lack diagnostics.
  for dump in files:
   with dump.open('rb') as stream:
    head=np.fromfile(stream,np.int32,5);stamp=np.fromfile(stream,np.float64,1);v=np.fromfile(stream,np.float64).reshape(-1,12)
   assert head[0]==1 and head[1]==patch and head[3]==h['it'] and head[4]==len(v) and abs(stamp[0]-t)<1e-12
   idx=np.rint((v[:,:3]-h['xmin'])/dx-.5).astype(int);assert np.all(idx>=0) and np.all(idx<h['domain_nx']);ix=tuple(idx.T);data[ix]=v[:,3:];np.add.at(coverage,ix,1)
  assert np.all(coverage==1) and np.isfinite(data).all()
  vol=data[...,8];div=abs(data[...,0]);j2=np.sum(data[...,1:4]**2,axis=-1);divdef2=(data[...,0]-data[...,4])**2;jdef2=np.sum((data[...,1:4]-data[...,5:8])**2,axis=-1)
  global_coord=local_coords(pos,0);gt,gp=global_coord[...,1],np.mod(global_coord[...,2],2*np.pi)
  common=(r>1.2)&(r<1.8)&(gt>np.pi/3)&(gt<2*np.pi/3)&(gp>np.deg2rad(100))&(gp<np.deg2rad(190));seam=np.zeros(r.shape,bool)
  for frame in [0,1]:
   q=local_coords(pos,frame);tt,pp=q[...,1],q[...,2];lo=.1875*np.pi;hi=.8125*np.pi;width=np.deg2rad(10)
   seam|=(np.minimum(abs(tt-lo),abs(tt-hi))<width)&(pp>=-hi-width)&(pp<=hi+width)
   seam|=(np.minimum(abs(pp+hi),abs(pp-hi))<width)&(tt>=lo-width)&(tt<=hi+width)
  masks=dict(global_=np.ones(r.shape,bool),common=common,common_seam=common&seam,common_off_seam=common&~seam,inner=r<1.1,outer=r>1.9);masks['global']=masks.pop('global_')
  for key,mask in masks.items():
   if not mask.any():continue
   vv=vol[mask];fields=[np.ones_like(r),mach**2,abs(rho-rhoex),rhoex,abs(pr-pex),pex,Berr,Bn,Berr**2,div**2,j2,rho,w[...,4],divdef2,jdef2,alfven_mach2,vel2,np.divide(div**2,Bn**2,out=np.zeros_like(div),where=Bn>0),np.divide(divdef2,Bn**2,out=np.zeros_like(div),where=Bn>0)]
   acc[key]+=np.array([np.sum(z[mask]*vv) for z in fields]);maxima[key]=np.maximum(maxima[key],[mach[mask].max(),div[mask].max(),np.sqrt(j2[mask]).max(),np.sqrt(alfven_mach2[mask]).max(),np.max(abs(pr[mask]-pex[mask])/pex[mask])])
 for key,z in acc.items():
  if z[0]==0:continue
  rows.append(dict(time=t,crossings=t/info['tau'],region=key,volume=z[0],mach_rms=np.sqrt(z[1]/z[0]),mach_max=maxima[key][0],rho_L1=z[2]/z[3],pressure_L1=z[4]/z[5],B_L1=z[6]/z[7] if z[7]>0 else None,B_error_rms=np.sqrt(z[8]/z[0]),div_rms=np.sqrt(z[9]/z[0]),div_max=maxima[key][1],current_rms=np.sqrt(z[10]/z[0]),current_max=maxima[key][2],mass=z[11],energy=z[12],div_defect_rms=np.sqrt(z[13]/z[0]),current_defect_rms=np.sqrt(z[14]/z[0]),rho_min=rho_min,pressure_min=pmin,beta_min=beta_min if np.isfinite(beta_min) else None,alfven_mach_rms=np.sqrt(z[15]/z[0]) if z[7]>0 else None,alfven_mach_max=maxima[key][3] if z[7]>0 else None,velocity_rms=np.sqrt(z[16]/z[0]),div_relative_rms=np.sqrt(z[17]/z[0]) if z[7]>0 else None,div_defect_relative_rms=np.sqrt(z[18]/z[0]) if z[7]>0 else None,pressure_relative_max=maxima[key][4]))
streams=[np.atleast_2d(np.loadtxt(f)) for f in sorted(a.run.glob('budget_r*.txt'))];assert len(streams)==info['ranks']
assert all(np.array_equal(q[:,0],streams[0][:,0]) for q in streams)
bud=sum(streams);bud[:,0]=streams[0][:,0]
base=bud[0];balance=[]
for q in bud:
 mass_delta=q[1]-base[1];energy_delta=q[2]-base[2]
 balance.append(dict(time=q[0],mass= q[1],energy=q[2],radial_mass_in=q[3]-base[3],radial_energy_in=q[4]-base[4],gravity_work=q[5]-base[5],other_source_energy=(q[6]-base[6])-(q[5]-base[5]),weighted_flux_mass=q[7]-base[7],weighted_flux_energy=q[8]-base[8],mass_closure=(mass_delta-(q[7]-base[7]))/max(abs(base[1]),1e-30),energy_closure=(energy_delta-(q[8]-base[8])-(q[6]-base[6]))/max(abs(base[2]),1e-30),mass_minus_radial=mass_delta-(q[3]-base[3]),energy_minus_radial_gravity=energy_delta-(q[4]-base[4])-(q[5]-base[5])))
last=max(x['time'] for x in rows);assert abs(last-t)<1e-12, 'Final snapshot lacks native diagnostics'
final=next(x for x in rows if x['time']==last and x['region']=='common')
passed=final['mach_rms']<.01 and final['mach_max']<.05 and final['rho_L1']<.02 and final['pressure_L1']<.02 and (final['B_L1']<.02 if final['B_L1'] is not None else final['B_error_rms']<1e-12)
closure=max(max(abs(q['mass_closure']),abs(q['energy_closure'])) for q in balance)
(a.run/'analysis.json').write_text(json.dumps(dict(rows=rows,budget=balance,engineering_pass=bool(passed),budget_closure_max=closure),indent=2)+'\n');print(json.dumps(dict(final=final,engineering_pass=bool(passed),budget_closure_max=closure)),flush=True)
if closure>1e-10:raise SystemExit('discrete budget closure check failed')

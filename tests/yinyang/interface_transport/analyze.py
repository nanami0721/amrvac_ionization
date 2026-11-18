#!/usr/bin/env python3
from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
import numpy as np
from scipy.interpolate import RegularGridInterpolator
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));import check_mhd
from diagnostics import geometry,exact,rotation,divergence
src,root=map(Path,sys.argv[1:3]);u=check_mhd.reader(src);rows=[]
for d in sorted(root.glob('c*-p*-d*-s*')):
 if not (d/'run.json').exists():continue
 info=json.loads((d/'run.json').read_text())
 if info['returncode'] or info['debug']:continue
 case,direction,placement=info['case'],info['direction'],info['placement'];rot=rotation(placement)
 for snap in [0,1,2]:
  states=[];geoms=[];fields=[];exactfields=[]
  for patch,name in enumerate(['yin','yang']):
   h,w=check_mhd.load(u,d/f'adv_{name}{snap:04d}.dat');g=geometry(h,patch);B=np.einsum('...ij,...j->...i',g['basis'],w[...,5:8]);ex=exact(g['pos'],h['time'],case,direction,placement)
   states.append((h,w,ex));geoms.append(g);fields.append(B);exactfields.append(ex[1])
  interps=[RegularGridInterpolator(g['axes'],B,bounds_error=False,fill_value=None) for g,B in zip(geoms,fields)]
  exinterps=[RegularGridInterpolator(g['axes'],B,bounds_error=False,fill_value=None) for g,B in zip(geoms,exactfields)]
  totals={key:0. for key in ['vol','err','den','err2','den2','roierr','roiden','seamerr','seamden','div2','divex2','divdefect2','divvol','seamdiv2','seamvol','Berr2','berrden','v2','rho2','mag2','Bflow2','psi2','p2','sroierr','sroiden']}
  pmin=rhomin=1e99;peak=0.;seammax=0.
  for patch,((h,w,ex),g,B) in enumerate(zip(states,geoms,fields)):
   erho,eB,ev,center=ex;vol=g['vol'];roi=np.linalg.norm(g['pos']-center,axis=-1)<.54;seam=g['seam']
   v=np.einsum('...ij,...j->...i',g['basis'],w[...,1:4])/w[...,0,None];pr=(h['gamma']-1)*(w[...,4]-.5*w[...,0]*np.sum(v*v,-1)-.5*np.sum(B*B,-1))
   if case==0:err=abs(w[...,0]-erho);signal=erho-1;amplitude=w[...,0]-1
   else:err=np.linalg.norm(B-eB,axis=-1);signal=np.linalg.norm(eB-.3*rot[:,1],axis=-1);amplitude=np.linalg.norm(B-.3*rot[:,1],axis=-1)
   for key,value in [('vol',vol.sum()),('err',np.sum(err*vol)),('den',np.sum(signal*vol)),('err2',np.sum(err*err*vol)),('den2',np.sum(signal*signal*vol)),('roierr',np.sum(err[roi]*vol[roi])),('roiden',np.sum(signal[roi]*vol[roi])),('seamerr',np.sum(err[seam]*vol[seam])),('seamden',np.sum(signal[seam]*vol[seam])),('Berr2',np.sum(np.sum((B-eB)**2,-1)*vol)),('berrden',np.sum(np.sum(eB*eB,-1)*vol)),('v2',np.sum(np.sum((v-ev)**2,-1)*vol)),('rho2',np.sum((w[...,0]-erho)**2*vol)),('mag2',np.sum((np.linalg.norm(B,axis=-1)-np.linalg.norm(eB,axis=-1))**2*vol)),('Bflow2',np.sum(np.sum((B-eB)*rot[:,0],axis=-1)**2*vol)),('psi2',np.sum(w[...,8]**2*vol)),('p2',np.sum((pr-.1)**2*vol)),('sroierr',np.sum(err[roi&seam]*vol[roi&seam])),('sroiden',np.sum(signal[roi&seam]*vol[roi&seam]))]:totals[key]+=float(value)
   div=divergence(g,w[...,5:8],interps[1-patch],patch);localexact=np.einsum('...ji,...j->...i',g['basis'],eB);dex=divergence(g,localexact,exinterps[1-patch],patch)
   core=np.zeros_like(roi);core[1:-1]=True;sc=core&seam
   for key,value in [('div2',np.sum(div[core]**2*vol[core])),('divex2',np.sum(dex[core]**2*vol[core])),('divdefect2',np.sum((div[core]-dex[core])**2*vol[core])),('divvol',vol[core].sum()),('seamdiv2',np.sum(div[sc]**2*vol[sc])),('seamvol',vol[sc].sum())]:totals[key]+=float(value)
   seammax=max(seammax,float(np.max(abs(div[sc]))));peak=max(peak,float(amplitude.max()));pmin=min(pmin,float(pr.min()));rhomin=min(rhomin,float(w[...,0].min()))
  t=totals
  result=dict(name=d.name,**info,snapshot=snap,time=h['time'],steps=h['it'],relative_L1=t['err']/t['den'],relative_L2=np.sqrt(t['err2']/t['den2']),roi_L1=t['roierr']/t['roiden'],seam_L1_contribution=t['seamerr']/t['den'],seam_signal_fraction=t['seamden']/t['den'],seam_roi_L1=(t['sroierr']/t['sroiden'] if t['sroiden']>1e-8*t['den'] else None),div_rms=np.sqrt(t['div2']/t['divvol']),div_exact_sample_rms=np.sqrt(t['divex2']/t['divvol']),div_defect_rms=np.sqrt(t['divdefect2']/t['divvol']),seam_div_rms=np.sqrt(t['seamdiv2']/t['seamvol']),seam_div_max=seammax,B_relative_rms=np.sqrt(t['Berr2']/t['berrden']),velocity_rms=np.sqrt(t['v2']/t['vol']),rho_rms=np.sqrt(t['rho2']/t['vol']),B_magnitude_rms=np.sqrt(t['mag2']/t['vol']),B_transport_component_error_rms=np.sqrt(t['Bflow2']/t['vol']),psi_rms=np.sqrt(t['psi2']/t['vol']),pressure_rms_error=np.sqrt(t['p2']/t['vol']),peak=peak,pressure_min=pmin,rho_min=rhomin,volume_bias=t['vol']/(4*np.pi/3*(h['xmax'][0]**3-h['xmin'][0]**3))-1)
  assert pmin>0 and rhomin>0
  rows.append(result)
(root.parent/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
for x in rows:
 if x['snapshot']==2:print(x['name'],f"L1={x['relative_L1']:.5f} ROI={x['roi_L1']:.5f} seam={x['seam_L1_contribution']:.5f} div={x['div_rms']:.3e} defect={x['div_defect_rms']:.3e}")

#!/usr/bin/env python3
"""Compare get_divb in the same physical regions; no full state arrays required."""
from pathlib import Path
import argparse,json,sys
sys.dont_write_bytecode=True
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import check_mhd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from diagnostics import geometry,basis,rotation,local_coords
from analyze_native import exact_fv

def seam_mask(pos):
    # Original Yin-Yang patch edges, also used for the spherical interior control.
    lo=np.array([1.,.1875*np.pi,-.8125*np.pi]);hi=np.array([2.,.8125*np.pi,.8125*np.pi])
    width=np.deg2rad(10.);mask=np.zeros(pos.shape[:-1],bool)
    for patch in [0,1]:
        q=local_coords(pos,patch);t,p=q[...,1],q[...,2]
        mask|=(np.minimum(abs(t-lo[1]),abs(t-hi[1]))<width)&(p>=lo[2]-width)&(p<=hi[2]+width)
        mask|=(np.minimum(abs(p-lo[2]),abs(p-hi[2]))<width)&(t>=lo[1]-width)&(t<=hi[1]+width)
    return mask

def patches(util,directory,grid,time):
    for patch in ([0] if grid=='spherical' else [0,1]):
        prefix='adv' if grid=='spherical' else ['adv_yin','adv_yang'][patch]
        match=[]
        for f in sorted(directory.glob(prefix+'[0-9][0-9][0-9][0-9].dat')):
            with f.open('rb') as stream:h=util.get_header(stream)
            if abs(h['time']-time)<1e-12:match.append(h)
        assert len(match)==1,(directory,time)
        h=match[0]
        if grid=='yinyang':g=geometry(h,patch)
        else:
            dx=(h['xmax']-h['xmin'])/h['domain_nx']
            axes=[a+(np.arange(n)+.5)*dd for a,n,dd in zip(h['xmin'],h['domain_nx'],dx)]
            r,t,p=np.meshgrid(*axes,indexing='ij');bb=basis(t,p,0)
            vol=(r*r+dx[0]**2/12)*dx[0]*2*np.sin(t)*np.sin(dx[1]/2)*dx[2]
            g=dict(axes=axes,dx=dx,r=r,t=t,p=p,basis=bb,pos=r[...,None]*bb[..., :,0],vol=vol)
        g['patch']=patch
        data=np.full((*h['domain_nx'],10),np.nan);coverage=np.zeros(h['domain_nx'],np.int8)
        for f in directory.glob(f'native_p{patch}_r*_i{h["it"]:07d}.bin'):
            with f.open('rb') as stream:
                head=np.fromfile(stream,np.int32,5);stamp=np.fromfile(stream,np.float64,2);raw=np.fromfile(stream,np.float64)
            assert tuple(head[:2])==(1,patch) and head[3]==h['it']
            assert abs(stamp[0]-time)<1e-12 and raw.size==head[4]*10
            v=raw.reshape(-1,10);idx=np.rint((v[:,:3]-h['xmin'])/g['dx']-.5).astype(int)
            assert np.all(idx>=0) and np.all(idx<h['domain_nx'])
            assert np.max(abs((idx+.5)*g['dx']+h['xmin']-v[:,:3]))<1e-12
            ix=tuple(idx.T);data[ix]=v;np.add.at(coverage,ix,1)
        assert np.all(coverage==1) and np.isfinite(data).all()
        ref,vol=exact_fv(g,time,dict(case=1,direction=1,placement=1))
        difference=float(abs(ref-data[...,4]).max());assert difference<2e-12
        assert np.allclose(vol,data[...,5],rtol=2e-13,atol=1e-16)
        yield h,g,data,difference

def main():
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('spherical_root',type=Path);p.add_argument('yinyang_root',type=Path);p.add_argument('--scales',type=int,nargs='+',default=[1,2]);a=p.parse_args()
    util=check_mhd.reader(a.source);rows=[]
    regions=['roi','radial_first_layer','inner_fixed','outer_fixed','roi_minus_first','interior','interior_seam','interior_off_seam']
    for scale in a.scales:
        for grid,root in [('spherical',a.spherical_root),('yinyang',a.yinyang_root)]:
            directory=root/(f's{scale}' if grid=='spherical' else f'c1-p1-d1-s{scale}')
            info=json.loads((directory/'run.json').read_text());assert info['returncode']==0 and info['stderr']==0
            for time in [0.,.3,.6]:
                accum={key:dict(volume=0.,num2=0.,ref2=0.,defect2=0.,psi2=0.,Berr2=0.,max=-1.,reference_max=0.,defect_max=0.,cells=0) for key in regions}
                check=0.
                for h,g,v,roundtrip in patches(util,directory,grid,time):
                    check=max(check,roundtrip);rad=g['r'];pos=g['pos'];center=rotation(1)@np.array([time-.3,0,1.5])
                    roi=np.linalg.norm(pos-center,axis=-1)<.54
                    edge=np.zeros(rad.shape,bool);edge[0]=True;edge[-1]=True
                    bulk=(rad>1.2)&(rad<1.8);seam=seam_mask(pos)
                    masks=dict(roi=roi,radial_first_layer=roi&edge,inner_fixed=roi&(rad<1.1),outer_fixed=roi&(rad>1.9),roi_minus_first=roi&~edge,interior=roi&bulk,interior_seam=roi&bulk&seam,interior_off_seam=roi&bulk&~seam)
                    for key,mask in masks.items():
                        if not mask.any():continue
                        z=accum[key];vol=g['vol'][mask];num=v[...,3][mask];ref=v[...,4][mask];defect=num-ref
                        z['volume']+=float(vol.sum());z['cells']+=int(mask.sum())
                        for name,field in [('num2',num),('ref2',ref),('defect2',defect),('psi2',v[...,6][mask])]:z[name]+=float(np.sum(field*field*vol))
                        z['Berr2']+=float(np.sum(np.sum(v[...,7:10][mask]**2,axis=-1)*vol))
                        z['reference_max']=max(z['reference_max'],float(abs(ref).max()));z['defect_max']=max(z['defect_max'],float(abs(defect).max()))
                        mx=float(abs(num).max())
                        if mx>z['max']:
                            i=np.argmax(np.where(mask,abs(v[...,3]),-1));idx=np.unravel_index(i,rad.shape)
                            z['max']=mx;z['peak']=dict(patch=g['patch'],r=float(rad[idx]),position=pos[idx].tolist())
                for key,z in accum.items():
                    assert z['volume']>0,(grid,scale,time,key)
                    rows.append(dict(grid=grid,scale=scale,time=time,region=key,volume=z['volume'],cells=z['cells'],div_rms=np.sqrt(z['num2']/z['volume']),div_max=z['max'],reference_rms=np.sqrt(z['ref2']/z['volume']),reference_max=z['reference_max'],reference_defect_correlation=(z['num2']-z['ref2']-z['defect2'])/(2*np.sqrt(z['ref2']*z['defect2'])) if z['ref2']*z['defect2']>0 else None,defect_rms=np.sqrt(z['defect2']/z['volume']),defect_max=z['defect_max'],psi_rms=np.sqrt(z['psi2']/z['volume']),B_error_rms=np.sqrt(z['Berr2']/z['volume']),peak=z['peak'],operator_check=check))
    (a.spherical_root/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
    for x in rows:
        if x['time']==.6:print(x['grid'],x['scale'],x['region'],f"RMS={x['div_rms']:.7f} max={x['div_max']:.7f} defect={x['defect_rms']:.7f}")
if __name__=='__main__':main()

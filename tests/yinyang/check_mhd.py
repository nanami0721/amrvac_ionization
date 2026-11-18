#!/usr/bin/env python3
"""Read native patch snapshots; report uniform-flow error and overlap-weighted totals.

Usage: python check_mhd.py SOURCE RUN_DIR [SNAPSHOT_INDEX]
Overlap weights use 16x16 angular subcell quadrature (1 outside, 1/2 within).
Their quadrature is approximate: report volume bias separately. This is not a
conservative overset correction or a discrete solver-divergence diagnostic.
"""
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np


def reader(source):
    spec = importlib.util.spec_from_file_location('datutil', Path(source)/'tools/python/amrvac_pytools/datfiles/reading/datfile_utilities.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load(util, path):
    with open(path, 'rb') as f:
        h = util.get_header(f)
        levels, indices, offsets = util.get_tree_info(f)
        assert np.all(levels == 1)
        w = np.empty((*h['domain_nx'], h['nw']))
        for index, offset in zip(indices, offsets):
            start = (index-1)*h['block_nx']
            sl = tuple(slice(a, a+n) for a,n in zip(start,h['block_nx']))
            w[sl] = util.get_single_block_data(f, offset, (*h['block_nx'],h['nw']))
    assert np.isfinite(w).all(), path
    return h,w


def diagnose(util, directory, snapshot, initial_case=0):
    sums = np.zeros(3)
    errors=[]; divs=[]; interface=[]; minrho=1e99; minp=1e99
    for patch,name in enumerate(('yin','yang')):
        h,w=load(util,Path(directory)/f'yy_{name}{snapshot:04d}.dat')
        dx=(h['xmax']-h['xmin'])/h['domain_nx']
        axes=[a+(np.arange(n)+.5)*d for a,n,d in zip(h['xmin'],h['domain_nx'],dx)]
        r,t,p=np.meshgrid(*axes,indexing='ij'); st,ct,sp,cp=np.sin(t),np.cos(t),np.sin(p),np.cos(p)
        basis=np.stack([np.stack([st*cp,st*sp,ct],-1),np.stack([ct*cp,ct*sp,-st],-1),np.stack([-sp,cp,np.zeros_like(p)],-1)],-1)
        if patch:
            basis=basis[..., [0,2,1],:];basis[...,0,:]*=-1
        rho=w[...,0];v=np.einsum('...ij,...j->...i',basis,w[...,1:4])/rho[...,None]
        B=np.einsum('...ij,...j->...i',basis,w[...,5:8])
        pressure=(h.get('gamma',5/3)-1)*(w[...,4]-.5*rho*np.sum(v*v,-1)-.5*np.sum(B*B,-1))
        minrho=min(minrho,float(rho.min()));minp=min(minp,float(pressure.min()))
        expected_rho=1.
        if initial_case==1:
            position=r[...,None]*basis[..., :,0]
            offset=position-np.array([.1,.2,-.1])*h['time']-np.array([-1.2,.8,0.])
            expected_rho=1.+.1*np.exp(-np.sum(offset*offset,-1)/.04)
        err=np.maximum.reduce([abs(rho-expected_rho),abs(pressure-1),np.max(abs(v-[.1,.2,-.1]),-1),np.max(abs(B-[.2,-.1,.3]),-1)])
        # Integrate the partition of unity over each angular cell. The two
        # rotated patches have weights 1/2 in their common physical region.
        angular=np.zeros((len(axes[1]),len(axes[2])))
        subdivisions=16
        for jt in range(subdivisions):
            tt=axes[1][:,None]+((jt+.5)/subdivisions-.5)*dx[1]
            ds=np.cos(tt-dx[1]/(2*subdivisions))-np.cos(tt+dx[1]/(2*subdivisions))
            for jp in range(subdivisions):
                pp=axes[2][None,:]+((jp+.5)/subdivisions-.5)*dx[2]
                ot=np.arctan2(np.sqrt((np.sin(tt)*np.cos(pp))**2+np.cos(tt)**2),np.sin(tt)*np.sin(pp))
                op=np.arctan2(np.cos(tt),-np.sin(tt)*np.cos(pp))
                overlap=(ot>=h['xmin'][1])&(ot<=h['xmax'][1])&(op>=h['xmin'][2])&(op<=h['xmax'][2])
                angular+=ds*dx[2]/subdivisions*np.where(overlap,.5,1.)
        vol=((r+dx[0]/2)**3-(r-dx[0]/2)**3)/3*angular[None,:,:]
        sums += [np.sum(vol),np.sum(rho*vol),np.sum(w[...,4]*vol)]
        errors.append(float(np.sum(err*vol)))
        band=np.zeros(r.shape,bool);band[:,:2,:]=True;band[:,-2:,:]=True;band[:,:,:2]=True;band[:,:,-2:]=True
        interface.append(float(err[band].max()))
        div=(np.gradient(r*r*w[...,5],dx[0],axis=0)/(r*r)+np.gradient(st*w[...,6],dx[1],axis=1)/(r*st)+np.gradient(w[...,7],dx[2],axis=2)/(r*st))
        divs.append(float(np.max(abs(div[1:-1,1:-1,1:-1]))))
    exact=4*np.pi/3*(h['xmax'][0]**3-h['xmin'][0]**3)
    return dict(time=h['time'],volume=sums[0],volume_relative_bias=sums[0]/exact-1,mass=sums[1],energy=sums[2],uniform_error_l1=sum(errors)/sums[0],interface_error_max=max(interface),divB_fd_max=max(divs),rho_min=minrho,pressure_min=minp)

if __name__=='__main__':
    util=reader(sys.argv[1]);result=diagnose(util,sys.argv[2],int(sys.argv[3]) if len(sys.argv)>3 else 1)
    initial=diagnose(util,sys.argv[2],0)
    result['mass_relative_change']=result['mass']/initial['mass']-1
    result['energy_relative_change']=result['energy']/initial['energy']-1
    print(json.dumps(result,indent=2))
    assert result['rho_min']>0 and result['pressure_min']>0

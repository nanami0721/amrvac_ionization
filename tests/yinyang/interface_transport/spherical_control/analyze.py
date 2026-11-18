#!/usr/bin/env python3
"""Compare identical local magnetic-error norms in spherical and Yin-Yang runs."""
from pathlib import Path
import sys, json, math
sys.dont_write_bytecode = True
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import check_mhd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from diagnostics import geometry, basis, exact, rotation

source, spherical, yinyang = map(Path, sys.argv[1:4])
util = check_mhd.reader(source)
rows = []
for grid, root in [('spherical', spherical), ('yinyang', yinyang)]:
    for scale in [1, 2, 3]:
        directory = root / f'c1-p1-d1-s{scale}'
        info = json.loads((directory / 'run.json').read_text())
        assert info['returncode'] == 0 and info['stderr'] == 0
        for snap in [0, 1, 2]:
            total = np.zeros(6)
            max_error = 0.
            pmin = rhomin = float('inf')
            for patch in ([0] if grid == 'spherical' else [0, 1]):
                prefix = 'adv' if grid == 'spherical' else ['adv_yin', 'adv_yang'][patch]
                h, w = check_mhd.load(util, directory / f'{prefix}{snap:04d}.dat')
                dx = (h['xmax'] - h['xmin']) / h['domain_nx']
                assert np.allclose(dx, [1/(16*scale), np.pi/(32*scale), np.pi/(32*scale)], rtol=0, atol=1e-14)
                if grid == 'yinyang':
                    g = geometry(h, patch)
                    b, pos, vol = g['basis'], g['pos'], g['vol']
                else:
                    axes = [a + (np.arange(n)+.5)*d for a,n,d in zip(h['xmin'], h['domain_nx'], dx)]
                    r,t,p = np.meshgrid(*axes, indexing='ij')
                    b = basis(t,p,0)
                    pos = r[...,None]*b[..., :,0]
                    vol = ((r+dx[0]/2)**3-(r-dx[0]/2)**3)/3 * (np.cos(t-dx[1]/2)-np.cos(t+dx[1]/2))*dx[2]
                erho, eB, ev, center = exact(pos, h['time'], 1, 1, 1)
                B = np.einsum('...ij,...j->...i', b, w[...,5:8])
                error = np.linalg.norm(B-eB, axis=-1)
                signal = np.linalg.norm(eB-.3*rotation(1)[:,1], axis=-1)
                roi = np.linalg.norm(pos-center, axis=-1) < .54
                total += [np.sum(error[roi]*vol[roi]), np.sum(signal[roi]*vol[roi]),
                          np.sum(error[roi]**2*vol[roi]), np.sum(signal[roi]**2*vol[roi]),
                          np.sum(vol[roi]), np.sum(roi)]
                max_error = max(max_error, float(error[roi].max()))
                pressure = (h['gamma']-1)*(w[...,4]-.5*np.sum(w[...,1:4]**2, axis=-1)/w[...,0]-.5*np.sum(w[...,5:8]**2, axis=-1))
                pmin = min(pmin, float(pressure.min()))
                rhomin = min(rhomin, float(w[...,0].min()))
            assert pmin > 0 and rhomin > 0
            rows.append(dict(grid=grid, scale=scale, snapshot=snap, time=h['time'], steps=h['it'],
                             roi_L1=total[0]/total[1], roi_L2=np.sqrt(total[2]/total[3]),
                             roi_B_max=max_error, roi_volume=total[4], roi_signal_integral=total[1],
                             rho_min=rhomin, pressure_min=pmin))
orders = []
for grid in ['spherical','yinyang']:
    for snap in [1,2]:
        q = sorted([x for x in rows if x['grid']==grid and x['snapshot']==snap], key=lambda x:x['scale'])
        for a,b in zip(q,q[1:]):
            orders.append(dict(grid=grid, time=b['time'], scales=[a['scale'],b['scale']],
                **{key:math.log(a[key]/b[key])/math.log(b['scale']/a['scale']) for key in ['roi_L1','roi_L2','roi_B_max']}))
assert max(x['roi_L1'] for x in rows if x['snapshot']==0) < 1e-12
out = spherical.parent / 'comparison.json'
out.write_text(json.dumps(dict(rows=rows, orders=orders), indent=2)+'\n')
for x in rows:
    if x['snapshot']==2: print(x['grid'], x['scale'], 'L1', x['roi_L1'], 'L2', x['roi_L2'], 'steps', x['steps'])
for x in orders:
    if abs(x['time']-.6)<1e-12: print(x)

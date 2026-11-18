#!/usr/bin/env python3
"""Compare two patch reconstructions at identical physical sphere quadrature points.

This measures overlap consistency, NOT a conservative-interface flux. The common
sampling mask excludes cells lacking an all-interior trilinear stencil. Sample
counts and covered areas are reported so changing coverage cannot be hidden.
"""
import numpy as np


def basis(theta, phi, patch):
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    q = np.stack((np.stack((st*cp, st*sp, ct), -1),
                  np.stack((ct*cp, ct*sp, -st), -1),
                  np.stack((-sp, cp, np.zeros_like(phi)), -1)), -1)
    if patch:
        q = q[..., [0, 2, 1], :].copy()
        q[..., 0, :] *= -1
    return q


def reconstruct(header, state, xyz, patch):
    local = xyz.copy()
    if patch:
        local = local[:, [0, 2, 1]]
        local[:, 0] *= -1
    radius = np.linalg.norm(local, axis=1)
    coords = np.column_stack((radius, np.arccos(np.clip(local[:, 2]/radius, -1, 1)),
                              np.arctan2(local[:, 1], local[:, 0])))
    dx = (header['xmax']-header['xmin'])/header['domain_nx']
    q = (coords-header['xmin'])/dx-.5
    lo = np.floor(q).astype(int)
    valid = np.all((lo >= 0) & (lo+1 < header['domain_nx']), axis=1)
    field = np.full(xyz.shape, np.nan)
    # Reconstruct Cartesian B, not spherical components, at each shared point.
    indices = lo[valid]
    fraction = q[valid]-indices
    values = np.zeros((len(indices), 3))
    for a in range(2):
        for b in range(2):
            for c in range(2):
                shift = np.array([a, b, c]); k = indices+shift
                coord = header['xmin']+(k+.5)*dx
                weight = np.prod(np.where(shift, fraction, 1-fraction), axis=1)
                sample = state[k[:, 0], k[:, 1], k[:, 2], 5:8]
                values += weight[:, None]*np.einsum('nij,nj->ni', basis(coord[:, 1], coord[:, 2], patch), sample)
    field[valid] = values
    return field, valid


def compare(util, directory, snapshot, prefix='yy', nmu=64, nphi=128):
    from pathlib import Path
    import check_mhd
    headers, states = [], []
    for patch in ('yin', 'yang'):
        h, w = check_mhd.load(util, Path(directory)/f'{prefix}_{patch}{snapshot:04d}.dat')
        headers.append(h); states.append(w)
    if headers[0]['time'] != headers[1]['time']:
        raise ValueError('patch timestamps disagree')
    if not np.array_equal(headers[0]['domain_nx'], headers[1]['domain_nx']):
        raise ValueError('patch meshes disagree')
    mu = -1+(np.arange(nmu)+.5)*2/nmu
    phi = (np.arange(nphi)+.5)*2*np.pi/nphi
    mu, phi = np.meshgrid(mu, phi, indexing='ij')
    normal = np.column_stack((np.sqrt(1-mu.ravel()**2)*np.cos(phi.ravel()),
                              np.sqrt(1-mu.ravel()**2)*np.sin(phi.ravel()), mu.ravel()))
    h = headers[0]
    radii = h['xmin'][0]+np.array([.25, .5, .75])*(h['xmax'][0]-h['xmin'][0])
    rows = []
    for radius in radii:
        first, ok0 = reconstruct(headers[0], states[0], radius*normal, 0)
        second, ok1 = reconstruct(headers[1], states[1], radius*normal, 1)
        valid = ok0 & ok1
        if not valid.any():raise ValueError('no common all-interior samples')
        a, b, n = first[valid], second[valid], normal[valid]
        difference = a-b
        scale = float(np.mean(.5*(np.linalg.norm(a, axis=1)+np.linalg.norm(b, axis=1))))
        absolute = float(np.mean(np.linalg.norm(difference, axis=1)))
        area = float(4*np.pi*radius**2/(nmu*nphi))
        rows.append(dict(radius=float(radius),samples=int(valid.sum()),area=float(area*valid.sum()),
                         B_difference_l1=absolute,B_difference_relative_l1=absolute/scale if scale else None,
                         signed_flux_difference=float(area*np.sum(difference*n)),
                         absolute_normal_flux_difference=float(area*np.sum(np.abs(np.sum(difference*n, axis=1)))),
                         unsigned_flux_scale=float(area*np.sum(.5*(np.abs(np.sum(a*n, axis=1))+np.abs(np.sum(b*n, axis=1)))))))
    return dict(time=h['time'],nmu=nmu,nphi=nphi,normalization='mean magnitude of both Cartesian reconstructions',
                scope='common interior sphere samples; not a closed-sphere net flux or strict conservation test',spheres=rows)

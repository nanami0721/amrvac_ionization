"""Independent geometry preparation probe; does not execute AMRVAC/VACPP.

Run from an external task directory. Cell-centred angular ghosts only;
radial ghosts, staggered locations and AMR donors require separate checks.
"""
import json
import math


def rotate(v):
    return (-v[0], v[2], v[1])


def xyz(t, p):
    return (math.sin(t)*math.cos(p), math.sin(t)*math.sin(p), math.cos(t))


def angles(v):
    return math.atan2(math.hypot(v[0], v[1]), v[2]), math.atan2(v[1], v[0])


def coverage(n, ng, width):
    # Fixed interior cell counts, symmetric extension of physical patch bounds.
    delta = width * math.pi / (2*n)
    tmin, tmax = math.pi/4-delta, 3*math.pi/4+delta
    pmin, pmax = -3*math.pi/4-delta, 3*math.pi/4+delta
    dt, dp = (tmax-tmin)/n, (pmax-pmin)/(3*n)
    count = missing = 0
    roundtrip = 0.0
    sample = None
    for j in range(-ng, n+ng):
        for k in range(-ng, 3*n+ng):
            if 0 <= j < n and 0 <= k < 3*n:
                continue
            v = xyz(tmin+(j+0.5)*dt, pmin+(k+0.5)*dp)
            t, p = angles(rotate(v))
            recovered = rotate(xyz(t, p))
            roundtrip = max(roundtrip, max(abs(a-b) for a,b in zip(v,recovered)))
            a, b = (t-tmin)/dt-0.5, (p-pmin)/dp-0.5
            # Strict four-point bilinear stencil consisting of donor interior cells.
            good = 0 <= math.floor(a) < n-1 and 0 <= math.floor(b) < 3*n-1
            count += 1
            if not good:
                missing += 1
                if sample is None:
                    sample = [j, k, a, b]
    return dict(ntheta=n, nphi=3*n, ghosts=ng, delta_over_nominal_dtheta=width,
                ghost_points=count, missing_interior_stencils=missing,
                first_missing=sample, coordinate_roundtrip_max=roundtrip)


def main():
    # Orthogonal, proper, involutory rotation: explicit independent arithmetic.
    for v in [(1.,0.,0.), (0.,1.,0.), (0.,0.,1.), (1.3,-2.1,4.7)]:
        assert rotate(rotate(v)) == v
        assert sum(x*x for x in v) == sum(x*x for x in rotate(v))
    cases = [coverage(n,g,w) for n in (32,64) for g in (2,3,4)
             for w in (0.51,1.,2.,3.,4.)]
    assert max(c['coordinate_roundtrip_max'] for c in cases) < 3e-15
    print(json.dumps(dict(scope='geometric centre-ghost survey, not MHD validation',
                          cases=cases), indent=2))


if __name__ == '__main__':
    main()

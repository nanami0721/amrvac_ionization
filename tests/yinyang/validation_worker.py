#!/usr/bin/env python3
"""One bounded-memory analysis job; invoked by run_validation under its watchdog."""
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
if os.environ.get('AMRVAC_MEMORY_GUARD_ACTIVE') != '1':
    raise SystemExit('analysis requires the process-tree memory guard')
import numpy as np
import check_mhd as check

job = json.loads(Path(sys.argv[1]).read_text())
u = check.reader(Path(job['source']))
number = job['snapshot']
result = {}
if job['operation'] == 'compare':
    maximum = 0.0
    for patch in ('yin', 'yang'):
        name = f'yy_{patch}{number:04d}.dat'
        ha, a = check.load(u, Path(job['left']) / name)
        hb, b = check.load(u, Path(job['right']) / name)
        assert ha['time'] == hb['time']
        assert np.isfinite(a).all() and np.isfinite(b).all()
        maximum = max(maximum, float(np.max(np.abs(a-b))))
        assert np.allclose(a, b, atol=1e-13, rtol=0)
        del a, b
    result = dict(max_difference=maximum)
elif job['operation'] == 'overlap':
    import overlap_diagnostics
    result = overlap_diagnostics.compare(u, Path(job['run']), number)
elif job['operation'] == 'diagnose':
    result = check.diagnose(u, Path(job['run']), number)
elif job['operation'] == 'positive':
    for patch in ('yin', 'yang'):
        h, w = check.load(u, Path(job['run']) / f'yy_{patch}{number:04d}.dat')
        assert np.isfinite(w).all() and w[..., 0].min() > 0
        pressure = (h.get('gamma', 5/3)-1)*(w[..., 4] - .5*np.sum(w[..., 1:4]**2, -1)/w[..., 0] - .5*np.sum(w[..., 5:8]**2, -1))
        assert pressure.min() > 0
        result[patch] = dict(min_density=float(w[..., 0].min()), min_pressure=float(pressure.min()))
        del w, pressure
else:
    raise ValueError('unknown analysis operation')
Path(job['output']).write_text(json.dumps(result, indent=2)+'\n')

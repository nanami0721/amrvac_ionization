#!/usr/bin/env python3
"""Compare matched forward-seam runs with different radii on the same shell."""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('narrow',type=Path);p.add_argument('wide',type=Path);a=p.parse_args()
n=json.loads((a.narrow/'analysis.json').read_text());w=json.loads((a.wide/'analysis.json').read_text())
assert all(n['geometry'][k]==w['geometry'][k] for k in ['radial_min','radial_max'])
assert n['geometry']['radius']<w['geometry']['radius']
rows=[]
for scale in [1,2]:
 key=f'debug-d1-p1-s{scale}-np2-a1'
 x=n['cases'][key];y=w['cases'][key]
 # Parameter files must differ only in radius, keeping dt, grid and all methods.
 def without_radius(path):return '\n'.join(line for line in path.read_text().splitlines() if not line.strip().startswith('radius='))
 assert without_radius(a.narrow/key/'amrvac.par')==without_radius(a.wide/key/'amrvac.par')
 ratio=y['final']['relative_l2']/x['final']['relative_l2']
 rows.append(dict(scale=scale,narrow_initial_l2=x['initial']['relative_l2'],wide_initial_l2=y['initial']['relative_l2'],narrow_final_l2=x['final']['relative_l2'],wide_final_l2=y['final']['relative_l2'],wide_over_narrow_final_l2=ratio,narrow_peak=x['final']['peak_fraction'],wide_peak=y['final']['peak_fraction']))
 assert ratio<1
out=dict(narrow_geometry=n['geometry'],wide_geometry=w['geometry'],rows=rows,interpretation='Fixed travel time and physical distance; width changes both spatial sampling and travel distance in units of loop width. This is not a convergence-order measurement.')
(a.wide/'width-comparison.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

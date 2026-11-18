#!/usr/bin/env python3
"""Run sequentially to avoid overlapping timing jobs. Source/builds must be external."""
from pathlib import Path
import argparse,subprocess,shutil,time,json
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--executable',type=Path,required=True);p.add_argument('--fine',action='store_true');p.add_argument('--coarse',action='store_true');p.add_argument('--accuracy',action='store_true');p.add_argument('--seam',action='store_true');a=p.parse_args()
a.root=a.root.resolve()
source_root=Path(__file__).resolve().parents[3]
if a.root==source_root or source_root in a.root.parents:p.error('run directory must be outside source repository')
a.root.mkdir(parents=True,exist_ok=True)
for scale in ([3] if a.accuracy else ([2] if a.fine else [1])):
 for tr in ([0] if a.accuracy else ([2] if a.seam else [0,1])):
  for geom in (['yinyang'] if a.accuracy else ['spherical','yinyang']):
   name=f'{geom}-s{scale}-t{tr}';d=a.root/name;d.mkdir(exist_ok=False);shutil.copy2(a.executable,d/'amrvac')
   nt,np_,lo,hi,pl,ph=(32*scale,64*scale,0,.5,0,1) if geom=='spherical' else (20*scale,52*scale,.09375,.40625,-.40625,.40625)
   nr=16*scale
   if a.accuracy:nr,nt,np_=40,48,124
   theta='pole' if geom=='spherical' else 'cont';phi='periodic' if geom=='spherical' else 'cont'
   par=f'''&filelist
 base_filename='adv'
 autoconvert=.false.
/
&savelist
 itsave(1,1)=0
 itsave(1,2)=0
 ditsave_log=10
 dtsave_dat=0.3d0
/
&stoplist
 it_max=100000
 time_max=0.6d0
/
&methodlist
 time_stepper='threestep'
 flux_scheme=20*'tvdlf'
 limiter=20*'minmod'
/
&boundlist
 typeboundary_min1=8*'special'
 typeboundary_max1=8*'special'
 typeboundary_min2=8*'{theta}'
 typeboundary_max2=8*'{theta}'
 typeboundary_min3=8*'{phi}'
 typeboundary_max3=8*'{phi}'
/
&meshlist
 refine_max_level=1
 block_nx1=8
 block_nx2=4
 block_nx3=4
 domain_nx1={nr}
 domain_nx2={nt}
 domain_nx3={np_}
 xprobmin1=1.d0
 xprobmax1=2.d0
 xprobmin2={lo}
 xprobmax2={hi}
 xprobmin3={pl}
 xprobmax3={ph}
/
&paramlist
 courantpar=0.3d0
/
&mhd_list
 typedivbfix='glm'
 boundary_divbfix=6*.false.
/
&advection_test
 test_geometry='{geom}'
 trajectory={tr}
/
'''
   (d/'amrvac.par').write_text(par)
   t=time.monotonic()
   with (d/'run.log').open('w') as out,(d/'run.err').open('w') as err:
    r=subprocess.run(['/usr/local/openmpi/bin/mpirun','--oversubscribe','-np','4','./amrvac'],cwd=d,stdout=out,stderr=err,timeout=1800)
   result=dict(case=name,returncode=r.returncode,wall=time.monotonic()-t,stderr_bytes=(d/'run.err').stat().st_size)
   (d/'run.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
   if r.returncode:raise SystemExit(r.returncode)

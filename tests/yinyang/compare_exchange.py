#!/usr/bin/env python3
"""Compare every stored cell before/after a communication-only change.

python compare_exchange.py BASELINE_RUN CANDIDATE_RUN --source SOURCE
Both runs must have identical snapshot names, times and physical configurations.
This checks native field arrays exactly (not file headers containing run counters).
"""
import argparse,json,re,sys
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
import check_mhd
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('baseline',type=Path);p.add_argument('candidate',type=Path)
p.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[2]);a=p.parse_args()
u=check_mhd.reader(a.source)
old={f.name:f for f in a.baseline.glob('yy_*.dat')};new={f.name:f for f in a.candidate.glob('yy_*.dat')}
assert old and old.keys()==new.keys(),'snapshot inventory differs'
for name in sorted(old):
 h,w=check_mhd.load(u,old[name]);hh,ww=check_mhd.load(u,new[name])
 assert h['time']==hh['time'] and h['it']==hh['it'],name
 assert w.shape==ww.shape and np.array_equal(w,ww),name
log=(a.candidate/'run.log').read_text()
m=re.search(r'Yin-Yang donor samples raw / unique per exchange:\s*(\d+)\s+(\d+)',log)
assert m,'missing unique-donor diagnostic'
raw,unique=map(int,m.groups());assert 0<=unique<=raw
calls=re.search(r'Yin-Yang exchange calls min / max:\s*(\d+)\s+(\d+)',log)
assert calls and calls[1]==calls[2],'rank exchange counts differ'
assert (a.candidate/'run.err').stat().st_size==0,'runtime stderr is not empty'
print(json.dumps(dict(snapshots=len(old),max_abs_difference=0,raw_samples=raw,unique_samples=unique,removed_fraction=1-unique/raw if raw else 0),indent=2))

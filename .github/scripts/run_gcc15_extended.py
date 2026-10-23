#!/usr/bin/env python3
"""Two additive controls on one runner and common MPI; untouched solver source."""
from pathlib import Path
import os,sys,json,hashlib,shutil,subprocess,re,math
ROOT=Path(os.environ['GITHUB_WORKSPACE']);OUT=ROOT/'results';MPI=Path(os.environ['RUNNER_TEMP'])/'comparison-openmpi'
SHA='91f5091258cd8948f3a189a53c1c19ee15eb1974'
COMP=str(Path(shutil.which('gfortran-15')).resolve())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def manifest(p):return {str(f.relative_to(p)):sha(f) for f in p.rglob('*') if f.is_file()}
def run(cmd,cwd,env,log):
 with log.open('w') as f:return subprocess.run([str(v) for v in cmd],cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
mpi_before=manifest(MPI)
comp_before={p:sha(p) for p in [COMP,subprocess.check_output([COMP,'-print-prog-name=f951'],text=True).strip()]}
write(OUT/'shared-mpi-sha256.json',mpi_before);write(OUT/'compiler-sha256.json',comp_before)
historical=json.loads((ROOT/'.github/scripts/gcc15_historical_hashes.json').read_text())['passing']
summary={};all_audits={};common_libs=None
for arm,flags,full in [('gcc15-march-x86-64',['-march=x86-64'],False),('gcc15-AB-full-suite',['-ffp-contract=off','-fno-tree-loop-vectorize'],True)]:
 out=OUT/arm;out.mkdir();source=ROOT/('source-'+arm)
 subprocess.run(['git','clone','--shared','--no-checkout',str(ROOT/'source-seed'),str(source)],check=True)
 subprocess.run(['git','checkout','--detach',SHA],cwd=source,check=True)
 wrapper=out/'compiler-wrapper.py'
 wrapper.write_text('#!/usr/bin/python3\nimport os,sys,json\n'+f'compiler={COMP!r}\nflags={flags!r}\nlog={str(out/"compiler-invocations.jsonl")!r}\n'+'cmd=[compiler]+sys.argv[1:]+flags\nwith open(log,"a") as f:f.write(json.dumps(cmd)+"\\n")\nos.execv(compiler,cmd)\n');wrapper.chmod(0o755)
 env=dict(os.environ,FC=str(wrapper),OMPI_FC=str(wrapper),EXPECTED_GCC_MAJOR='15',EXPECTED_GCC_VERSION='15.2.0',COMPARISON_ARM=arm,CAPTURE_INSTRUCTION_AUDIT='1',AMRVAC_ALL_CASES='1' if full else '0',AMRVAC_BUILD_DIMS='1d1,1d2,1d3,2d2,2d3,3d3' if full else '1d1,1d3,2d2,3d3',AMRVAC_BUILD_JOBS='4',LD_LIBRARY_PATH=str(MPI)+'/lib:/usr/lib/x86_64-linux-gnu')
 for k in ['LIBRARY_PATH','LD_PRELOAD','AMRVAC_CASE_FILTER','MAKEFLAGS','MFLAGS']:env.pop(k,None)
 write(out/'settings.json',dict(flags=flags,full_suite=full,source=SHA,compiler=COMP,mpi_prefix=str(MPI),environment_controls={k:env.get(k) for k in ['ARCH','NUM_PROCS','OMP_NUM_THREADS','HWLOC_LIBXML','LD_PRELOAD']}))
 probe=out/'fp-probe.f90';probe.write_text('subroutine probe(a,b,c,y)\nreal(8)::a,b,c,y\ny=a*b+c\nend subroutine\n')
 for kind in ['target','optimizers','fortran']:
  assert run(['mpif90','-O2','-Q','--help='+kind,'-c',probe,'-o',out/'probe.o'],out,env,out/('compiler-'+kind+'.txt'))==0
 print('START_ARM',arm,flush=True)
 rc=run([sys.executable,ROOT/'.github/scripts/run_gcc15_extended_cases.py',source,out],ROOT,env,out/'execution.log')
 if not (out/'results.json').exists():raise RuntimeError('No case results for '+arm)
 records=json.loads((out/'results.json').read_text());assert len(records)==(535 if full else 12) and rc in (0,1)
 assert not (out/'source-diff.txt').read_bytes()
 commands=[json.loads(l) for l in (out/'compiler-invocations.jsonl').read_text().splitlines()]
 assert all(c[0]==COMP and c[-len(flags):]==flags for c in commands)
 rows=[];libs={}
 for r in records:
  p=out/'cases'/r['directory']/r['log'][:-4];key=r['directory']+'/'+r['log']
  stats={'byte_identical':False,'shape_equal':False,'finite':False,'within_tolerance':False}
  if (p/'actual.log').exists():
   assert sha(p/'actual.log')==r['actual_sha256'] and sha(p/'reference.log')==r['reference_sha256']
   def table(f):return [[float(x.replace('D','E').replace('d','e')) for x in l.split()] for l in f.read_text().splitlines() if l.strip() and not l.lstrip().startswith(('#','!'))]
   try:
    a=table(p/'actual.log');b=table(p/'reference.log');shape=len(a)==len(b) and all(len(x)==len(y) for x,y in zip(a,b));finite=all(math.isfinite(v) for row in a+b for v in row)
    diffs=[(abs(x-y),1e-8+0.5e-5*(abs(x)+abs(y))) for ar,br in zip(a,b) for x,y in zip(ar,br)] if shape else []
    stats=dict(byte_identical=(p/'actual.log').read_bytes()==(p/'reference.log').read_bytes(),shape_equal=shape,finite=finite,within_tolerance=shape and finite and all(d<=t for d,t in diffs),mismatched_values=sum(d>t for d,t in diffs),max_tolerance_ratio=max((d/t for d,t in diffs),default=None))
   except ValueError as e:stats['parse_error']=str(e)
  valid=r['make_returncode']==0 and r['finished_amrvac'] and r['compare_returncode'] is not None
  passed=valid and r['compare_returncode']==0 and stats['within_tolerance']
  if valid:assert stats['within_tolerance']==(r['compare_returncode']==0), 'Comparator disagreement '+key
  ins={}
  if (p/'instruction-audit.json').exists():
   ins=json.loads((p/'instruction-audit.json').read_text());assert ins['instruction_count']>0 and ins['fma_instruction_count']==0
  if (p/'linked-libraries.txt').exists():
   txt=(p/'linked-libraries.txt').read_text();assert 'not found' not in txt
   paths=re.findall(r'(?:=>\s+|^\s*)(/\S+)',txt,re.M);assert len(paths)>4
   hashes={str(Path(v).resolve()):sha(v) for v in paths}
   if common_libs is None:common_libs=hashes
   assert hashes==common_libs
   libs[key]=hashes
  rows.append(dict(case=key,status=r['status'],passed=passed,reference=stats,historical_passing_byte_identical=(r['actual_sha256']==historical[key]) if key in historical else None,fma_instruction_count=ins.get('fma_instruction_count')))
 write(out/'crosscheck.json',rows);write(out/'runtime-libraries.json',libs)
 summary[arm]=dict(cases=len(rows),passed=sum(r['passed'] for r in rows),failed=[r['case'] for r in rows if not r['passed']],reported12_historical_byte_identical=sum(r['historical_passing_byte_identical'] is True for r in rows),flags=flags,runner_returncode=rc)
 all_audits[arm]=rows;write(OUT/'progress.json',summary)
 print('ARM_COMPLETE',arm,json.dumps(summary[arm]),flush=True)
 assert manifest(MPI)==mpi_before and all(sha(p)==v for p,v in comp_before.items())
result=dict(source=SHA,arms=summary,same_compiler_and_mpi_throughout=True,same_direct_runtime_libraries=common_libs,solver_sources_unchanged=True)
write(OUT/'extended-results.json',result)
with open(os.environ['GITHUB_STEP_SUMMARY'],'a') as f:f.write('\n## GCC 15 additional controls\n\n```json\n'+json.dumps(result,indent=2)+'\n```\n')
print('EXTENDED_COMPLETE',json.dumps(result),flush=True)
sys.exit(0 if all(v['passed']==v['cases'] for v in summary.values()) else 1)

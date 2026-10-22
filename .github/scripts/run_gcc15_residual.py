#!/usr/bin/env python3
"""Controlled GCC 15 residual diagnosis; all solver checkouts remain unchanged."""
from pathlib import Path
import os,sys,json,subprocess,hashlib,shutil,re,gzip,importlib.util,time
ROOT=Path(os.environ['GITHUB_WORKSPACE']);OUT=ROOT/'results';TEMP=Path(os.environ['RUNNER_TEMP']);MPI=TEMP/'comparison-openmpi'
SHA='91f5091258cd8948f3a189a53c1c19ee15eb1974'
CASE='hd/Coleman_1D/coleman4_LTE_HHe.log'
CASEPATH=Path('cases/hd/Coleman_1D/coleman4_LTE_HHe')
OFF=['-ffp-contract=off']
COMP={v:str(Path(shutil.which('gfortran-'+v)).resolve()) for v in ['14','15']}
VERS={'14':'14.3.0','15':'15.2.0'}
HISTORICAL=json.loads((ROOT/'.github/scripts/gcc15_historical_hashes.json').read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def command(args,cwd,env,log,check=True):
 with log.open('w') as f:return subprocess.run([str(x) for x in args],cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT,check=check).returncode
def manifest(root):return {str(p.relative_to(root)):sha(p) for p in root.rglob('*') if p.is_file()}
initial_mpi=manifest(MPI)
initial_comp={}
for c in COMP.values():
 for p in [c,subprocess.check_output([c,'-print-prog-name=f951'],text=True).strip()]:initial_comp[p]=sha(p)
write(OUT/'shared-mpi-sha256.json',initial_mpi);write(OUT/'compilers-sha256.json',initial_comp)
summary={};arms={};common_runtime=None

def run_arm(name,version='15',global_flags=None,scoped_flags=None,modules=None,full=False):
 global common_runtime
 flags=OFF+(global_flags or []);scoped_flags=scoped_flags or [];modules=modules or []
 p=OUT/name;p.mkdir();source=ROOT/('source-'+name)
 subprocess.run(['git','clone','--shared','--no-checkout',str(ROOT/'source-seed'),str(source)],check=True)
 subprocess.run(['git','checkout','--detach',SHA],cwd=source,check=True)
 wrapper=p/'compiler-wrapper.py';(p/'vectorization').mkdir()
 wrapper.write_text('#!/usr/bin/python3\nimport os,sys,json\nfrom pathlib import Path\n'+
  f'compiler={COMP[version]!r}\nflags={flags!r}\nscoped={scoped_flags!r}\nmodules={modules!r}\nout={str(p)!r}\n'+
  'args=sys.argv[1:]\nsources=[Path(x).name.lower() for x in args if x.lower().endswith((".f",".f90"))]\n'+
  'selected=bool(sources and any((m=="eos" and (s.startswith("mod_eos_") or s=="mod_eos.f")) or s==m for s in sources for m in modules))\n'+
  'extra=flags+(scoped if selected else [])\n'+
  'if "-c" in args and sources: extra += ["-fopt-info-vec-optimized="+out+"/vectorization/"+sources[0]+".txt"]\n'+
  'cmd=[compiler]+args+extra\nwith open(out+"/compiler-invocations.jsonl","a") as f: f.write(json.dumps({"command":cmd,"sources":sources,"scoped":selected})+"\\n")\nos.execv(compiler,cmd)\n')
 wrapper.chmod(0o755)
 env=dict(os.environ,FC=str(wrapper),OMPI_FC=str(wrapper),EXPECTED_GCC_MAJOR=version,EXPECTED_GCC_VERSION=VERS[version],CAPTURE_INSTRUCTION_AUDIT='1',COMPARISON_ARM=name,AMRVAC_BUILD_JOBS='4',LD_LIBRARY_PATH=str(MPI)+'/lib:/usr/lib/x86_64-linux-gnu')
 env.pop('LIBRARY_PATH',None);env.pop('LD_PRELOAD',None)
 if not full:env.update(AMRVAC_CASE_FILTER=CASE,AMRVAC_BUILD_DIMS='1d1')
 else:env.pop('AMRVAC_CASE_FILTER',None);env.pop('AMRVAC_BUILD_DIMS',None)
 write(p/'settings.json',dict(compiler=COMP[version],version=VERS[version],flags=flags,scoped_flags=scoped_flags,modules=modules,full=full,source=SHA))
 probe=p/'fp-probe.f90';probe.write_text('subroutine probe(a,b,c,y)\nreal(8)::a,b,c,y\ny=a*b+c\nend subroutine\n')
 for kind in ['target','optimizers','fortran','params']:
  command(['mpif90','-O2','-Q','--help='+kind,'-c',probe,'-o',p/'probe.o'],p,env,p/('compiler-'+kind+'.txt'))
 print('START_ARM',name,flush=True)
 rc=command([sys.executable,ROOT/'.github/scripts/ubuntu_comparison.py',source,p],ROOT,env,p/'execution.log',False)
 if not (p/'results.json').exists():raise RuntimeError('Missing results for '+name)
 results=json.loads((p/'results.json').read_text());assert rc in (0,1) and len(results)==(12 if full else 1)
 assert all(r['make_returncode']==0 and r['finished_amrvac'] and r['status'] in ('PASS','NUMERICAL_MISMATCH') for r in results)
 assert not (p/'source-diff.txt').read_bytes()
 libs={}
 for f in p.glob('cases/**/linked-libraries.txt'):
  text=f.read_text();assert 'not found' not in text
  paths=re.findall(r'(?:=>\s+|^\s*)(/\S+)',text,re.M);assert len(paths)>4
  d={str(Path(x).resolve()):sha(x) for x in paths}
  if common_runtime is None:common_runtime=d
  assert d==common_runtime,'Runtime libraries differ: '+name
  libs[str(f.relative_to(p))]=d
 write(p/'runtime-libraries.json',libs)
 for f in p.glob('cases/**/instruction-audit.json'):
  audit=json.loads(f.read_text());assert audit['instruction_count']>0 and audit['fma_instruction_count']==0
 objects={}
 for f in sorted((source/'lib/1d1_debug').glob('*.o')):
  text=subprocess.check_output(['nm','--undefined-only',str(f)],text=True)
  vectors=[l.strip().split()[-1] for l in text.splitlines() if '_ZGV' in l]
  if vectors:objects[f.name]=vectors
 write(p/'object-vector-symbols.json',objects)
 binary=source/'tests/hd/Coleman_1D/amrvac'
 # Save the Coleman binary even if full validation finishes in another directory.
 with binary.open('rb') as src,gzip.open(p/'coleman-amrvac.gz','wb') as dst:shutil.copyfileobj(src,dst)
 for f in (source/'lib/1d1_debug').glob('mod_eos*.f'):shutil.copy2(f,p/f.name)
 h={r['directory']+'/'+r['log']:r['actual_sha256'] for r in results}
 summary[name]=dict(passed=sum(r['status']=='PASS' for r in results),cases=len(results),matches_historical_passing=sum(h[k]==HISTORICAL['passing'][k] for k in h),flags=flags,scoped_flags=scoped_flags,modules=modules)
 arms[name]=(source,p,env,results)
 write(OUT/'progress.json',summary)
 print('ARM_COMPLETE',name,json.dumps(summary[name]),flush=True)
 return summary[name]

# Same host, same MPI binary and same dynamic-library hashes.
run_arm('gcc14-off','14')
assert summary['gcc14-off']['matches_historical_passing']==1
run_arm('gcc15-off')
assert arms['gcc15-off'][3][0]['actual_sha256']==HISTORICAL['gcc15_off'][CASE],'GCC15 baseline changed; investigate environment before causal interpretation'

# Runtime probes change only the chosen vector math implementation.
source,p,env,_=arms['gcc15-off'];cwd=source/'tests/hd/Coleman_1D';binary=cwd/'amrvac';binary_hash=sha(binary)
probe_results={}
for name,functions in [('new-vector-math',['EXP','LOG','POW','SIN','COS']),('exp',['EXP']),('log',['LOG']),('pow',['POW']),('log10',['LOG10']),('sin-cos',['SIN','COS']),('all-vector-math',['EXP','LOG','POW','SIN','COS','LOG10'])]:
 d=OUT/('scalar-probe-'+name);d.mkdir();lib=d/'scalar-math.so'
 command(['gcc-15','-O2','-shared','-fPIC','-mavx2','-fno-builtin','-fno-tree-vectorize','-ffp-contract=off',*['-DSELECT_'+f for f in functions],ROOT/'.github/scripts/vector_math_scalar_probe.c','-lm','-o',lib],d,env,d/'build.log')
 actual=cwd/'coleman4_LTE_HHe.log';actual.unlink(missing_ok=True)
 # Pass LD_PRELOAD only to AMRVAC ranks, not to mpirun itself.
 rc=command(['mpirun','-np','2','env','LD_PRELOAD='+str(lib),'./amrvac','-i','coleman4_LTE_HHe.par'],cwd,env,d/'run.log',False)
 assert rc==0 and 'Finished AMRVAC' in (d/'run.log').read_text() and sha(binary)==binary_hash
 shutil.copy2(actual,d/'actual.log');shutil.copy2(cwd/'correct_output/coleman4_LTE_HHe.log',d/'reference.log')
 cmp=command([source/'tools/fortran/compare_logs','1e-5','1e-8',actual,cwd/'correct_output/coleman4_LTE_HHe.log'],cwd,env,d/'comparison.log',False)
 counts=re.findall(r'SCALAR_VECTOR_CALL (\w+) (\d+)',(d/'run.log').read_text());totals={}
 for fn,n in counts:totals[fn]=totals.get(fn,0)+int(n)
 probe_results[name]=dict(functions=functions,passed=cmp==0,actual_sha256=sha(actual),matches_historical_passing=sha(actual)==HISTORICAL['passing'][CASE],scalarized_vector_call_counts=totals,binary_sha256=binary_hash,preload_sha256=sha(lib))
 write(OUT/'scalar-probes.json',probe_results);print('SCALAR_PROBE',name,json.dumps(probe_results[name]),flush=True)

run_arm('gcc15-no-loop',global_flags=['-fno-tree-loop-vectorize'])
run_arm('gcc15-no-slp',global_flags=['-fno-tree-slp-vectorize'])
# Prefer a module-scoped change over a global change when it restores the baseline.
candidates=[]
for suffix,flag in [('no-loop','-fno-tree-loop-vectorize'),('no-slp','-fno-tree-slp-vectorize')]:
 if summary['gcc15-'+suffix]['matches_historical_passing']:
  name='gcc15-eos-'+suffix;run_arm(name,scoped_flags=[flag],modules=['eos'])
  if summary[name]['matches_historical_passing']:
   candidates.append((name,[],[flag],['eos']))
   objects=json.loads((OUT/'gcc15-off/object-vector-symbols.json').read_text())
   successful=[v['functions'][0].lower() for k,v in probe_results.items() if len(v['functions'])==1 and v['matches_historical_passing'] and v['scalarized_vector_call_counts']]
   priority=['mod_eos_lte_entropy.o','mod_eos_lte.o','mod_eos_interp.o','mod_eos_shared_functions.o']
   ranked=[o for o in priority if o in objects and (not successful or any(any(sym.endswith('_'+fn) for fn in successful) for sym in objects[o]))]
   for obj in ranked[:2]:
    module=obj[:-2]+'.f';tag='gcc15-'+obj[:-2]+'-'+suffix
    run_arm(tag,scoped_flags=[flag],modules=[module])
    if summary[tag]['matches_historical_passing']:
     candidates.insert(0,(tag,[],[flag],[module]));break
  candidates.append(('gcc15-'+suffix,[flag],[],[]))
  break
if not candidates:
 run_arm('gcc15-no-vector',global_flags=['-fno-tree-loop-vectorize','-fno-tree-slp-vectorize'])
 if summary['gcc15-no-vector']['matches_historical_passing']:candidates.append(('gcc15-no-vector',['-fno-tree-loop-vectorize','-fno-tree-slp-vectorize'],[],[]))
validation=[]
for origin,global_flags,scoped_flags,modules in candidates:
 name='validate-'+origin;run_arm(name,global_flags=global_flags,scoped_flags=scoped_flags,modules=modules,full=True)
 validation.append(name)
 if summary[name]['passed']==12:break
assert manifest(MPI)==initial_mpi and all(sha(p)==h for p,h in initial_comp.items())
result=dict(source=SHA,arms=summary,scalar_probes=probe_results,validation=validation,common_runtime_libraries=common_runtime,compiler_and_mpi_unchanged=True,all_solver_sources_unchanged=True,validated_fix=next((n for n in validation if summary[n]['passed']==12),None))
write(OUT/'residual-diagnosis.json',result)
with open(os.environ['GITHUB_STEP_SUMMARY'],'a') as f:f.write('\n## GCC 15 residual diagnosis\n\n```json\n'+json.dumps(result,indent=2)+'\n```\n')
print('DIAGNOSIS_COMPLETE',json.dumps(result),flush=True)
sys.exit(0 if result['validated_fix'] else 1)

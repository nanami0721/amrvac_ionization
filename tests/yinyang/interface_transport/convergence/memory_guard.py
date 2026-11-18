#!/usr/bin/env python3
"""macOS process-tree RSS/physical-footprint polling guard; not an OS hard limit."""
import argparse,ctypes,json,os,signal,subprocess,time,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--limit-mib',type=float,default=1024);p.add_argument('--out',type=Path,required=True);p.add_argument('--system-start-percent',type=float,default=65.);p.add_argument('--system-stop-percent',type=float,default=68.);p.add_argument('--timeout-seconds',type=float,default=3600.);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args();cmd=a.command[1:] if a.command[:1]==['--'] else a.command
if not cmd:p.error('command required')
if not 0<a.timeout_seconds<float('inf'):p.error('finite positive timeout required')
if not 0<a.system_start_percent<a.system_stop_percent<=68.:p.error('system thresholds must satisfy 0 < start < stop <= 68')
a.out.mkdir(parents=True,exist_ok=False)
class Usage(ctypes.Structure):
 _fields_=[('uuid',ctypes.c_ubyte*16)]+[(k,ctypes.c_uint64) for k in ['user','system','idle_wakeups','interrupt_wakeups','pageins','wired','resident','footprint','start','exit']]
lib=ctypes.CDLL('/usr/lib/libproc.dylib');lib.proc_pid_rusage.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_void_p];lib.proc_pid_rusage.restype=ctypes.c_int

def usage(pid):
 u=Usage()
 if lib.proc_pid_rusage(pid,0,ctypes.byref(u)):return None
 return u

def scan(root,known):
 raw=subprocess.check_output(['/bin/ps','-axo','pid=,ppid=,rss='],text=True,timeout=3)
 table={int(x[0]):(int(x[1]),int(x[2])*1024) for line in raw.splitlines() if len(x:=line.split())==3}
 ids={root}
 # Include tracked surviving descendants even if they have become orphaned.
 for pid,birth in known.items():
  u=usage(pid)
  if u and u.start==birth:ids.add(pid)
 while True:
  more={pid for pid,(parent,_) in table.items() if parent in ids}
  if more<=ids:break
  ids|=more
 rss=foot=0;live={}
 for pid in ids:
  if pid not in table:continue
  u=usage(pid)
  if u is None:
   try:os.kill(pid,0)
   except ProcessLookupError:continue
   raise RuntimeError(f'cannot read physical footprint for PID {pid}')
  live[pid]=u.start;rss+=max(table[pid][1],u.resident);foot+=u.footprint
 return rss,foot,live

def stop(proc,known):
 # Kill only this newly created session and verified descendants.
 try:os.killpg(proc.pid,signal.SIGTERM)
 except ProcessLookupError:pass
 for sig in [signal.SIGTERM,signal.SIGKILL]:
  for pid,birth in list(known.items()):
   u=usage(pid)
   if u and u.start==birth:
    try:os.kill(pid,sig)
    except ProcessLookupError:pass
  if sig==signal.SIGTERM:time.sleep(.3)
 try:proc.wait(timeout=5)
 except subprocess.TimeoutExpired:
  proc.kill();proc.wait()

def system_memory():
 total=int(subprocess.check_output(['/usr/sbin/sysctl','-n','hw.memsize'],text=True,timeout=3))
 raw=subprocess.check_output(['/usr/bin/vm_stat'],text=True,timeout=3)
 page=int(re.search(r'page size of (\d+) bytes',raw).group(1))
 free=int(re.search(r'Pages free:\s+(\d+)',raw).group(1))
 speculative=int(re.search(r'Pages speculative:\s+(\d+)',raw).group(1))
 def pages(name):return int(re.search(re.escape(name)+r':\s+(\d+)',raw).group(1))
 anonymous=pages('Anonymous pages')*page
 wired=pages('Pages wired down')*page
 compressor=pages('Pages occupied by compressor')*page
 cache=pages('File-backed pages')*page
 occupied=anonymous+wired+compressor
 if not 0<=occupied<=total:raise RuntimeError('invalid system memory counters')
 swap=subprocess.check_output(['/usr/sbin/sysctl','-n','vm.swapusage'],text=True,timeout=3).strip()
 return dict(total_bytes=total,occupied_bytes=occupied,occupied_percent=100*occupied/total,anonymous_bytes=anonymous,wired_bytes=wired,compressor_bytes=compressor,file_cache_bytes=cache,cache_inclusive_bytes=total-(free+speculative)*page,accounting='anonymous+wired+physical_compressor',swap=swap)

known={};reason=None;peak_rss=peak_foot=0;started=time.monotonic();samples=0;proc=None;peak_system=0.;initial_system=None
with (a.out/'command.log').open('w') as out,(a.out/'command.err').open('w') as err,(a.out/'samples.jsonl').open('w') as log:
 try:
  initial_system=system_memory();peak_system=initial_system['occupied_percent']
  if peak_system>=a.system_start_percent:
   reason='system_start_limit'
   raise RuntimeError('system memory too high to start task')
  proc=subprocess.Popen(cmd,stdout=out,stderr=err,start_new_session=True,env=dict(os.environ,AMRVAC_MEMORY_GUARD_ACTIVE='1',AMRVAC_MEMORY_LIMIT_MIB=str(a.limit_mib)))
  while True:
   rss,foot,live=scan(proc.pid,known);known.update(live);samples+=1;peak_rss=max(peak_rss,rss);peak_foot=max(peak_foot,foot)
   system=system_memory();peak_system=max(peak_system,system['occupied_percent'])
   log.write(json.dumps(dict(system=system,seconds=time.monotonic()-started,rss_bytes=rss,footprint_bytes=foot,pids=list(live)))+'\n');log.flush()
   if time.monotonic()-started>=a.timeout_seconds:
    reason='time_limit';stop(proc,known);break
   if system['occupied_percent']>=a.system_stop_percent:
    reason='system_memory_limit';stop(proc,known);break
   if max(rss,foot)>a.limit_mib*2**20:
    reason='memory_limit';stop(proc,known);break
   if proc.poll() is not None:
    _,_,remaining=scan(proc.pid,known)
    if any(pid!=proc.pid for pid in remaining):
     reason='surviving_descendants';stop(proc,known)
    break
   time.sleep(.2)
 except BaseException as exc:
  if reason is None:reason=f'monitor_error: {exc}'
  if proc is not None:stop(proc,known)
 finally:
  code=proc.returncode if proc is not None else -1
  report=dict(initial_system=initial_system,peak_system_percent=peak_system,system_start_percent=a.system_start_percent,system_stop_percent=a.system_stop_percent,command=cmd,limit_mib=a.limit_mib,sampling_interval_seconds=.2,timeout_seconds=a.timeout_seconds,peak_rss_mib=peak_rss/2**20,peak_footprint_mib=peak_foot/2**20,wall=time.monotonic()-started,samples=samples,returncode=code,stop_reason=reason)
  (a.out/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
raise SystemExit(125 if reason else code)

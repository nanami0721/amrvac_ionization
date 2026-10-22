#!/usr/bin/env python3
"""Single-host, single-compiler tests of target ISA and FP contraction."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(os.environ['GITHUB_WORKSPACE'])
TEMP = Path(os.environ['RUNNER_TEMP'])
OUTPUT = ROOT / 'results'
MPI = TEMP / 'comparison-openmpi'
ARMS = {'default': [], 'x86-64': ['-march=x86-64'],
        'contract-off': ['-ffp-contract=off']}
COMPILER = str(Path(shutil.which('gfortran-14')).resolve())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mpi_manifest():
    return {str(p.relative_to(MPI)): sha(p) for p in sorted(MPI.rglob('*')) if p.is_file()}


def compiler_manifest():
    frontend = subprocess.check_output([COMPILER, '-print-prog-name=f951'], text=True).strip()
    return {p: sha(p) for p in (COMPILER, frontend)}


def linked_manifest(log):
    text = log.read_text()
    if 'not found' in text:
        raise RuntimeError(f'Unresolved library: {log}')
    paths = re.findall(r'(?:=>\s+|^\s*)(/\S+)', text, re.M)
    if not paths:
        raise RuntimeError(f'No libraries recorded: {log}')
    return {str(Path(p).resolve()): sha(p) for p in paths}


def record_options(destination, env):
    probe = destination / 'fp-probe.f90'
    probe.write_text('''subroutine multiply_add(a,b,c,y)
  implicit none
  real(8), intent(in) :: a,b,c
  real(8), intent(out) :: y
  y=a*b+c
end subroutine multiply_add
''')
    for name, args in [('target', ['-Q', '--help=target']),
                       ('optimizer', ['-Q', '--help=optimizers'])]:
        with (destination / f'compiler-{name}.txt').open('w') as stream:
            subprocess.run(['mpif90', '-O2', '-c', str(probe), '-o', str(destination / 'fp-probe.o'), *args],
                           env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
    subprocess.run(['mpif90', '-O2', '-S', str(probe), '-o', str(destination / 'fp-probe.s')], env=env, check=True)
    return {'target': (destination / 'compiler-target.txt').read_text(),
            'optimizer': (destination / 'compiler-optimizer.txt').read_text()}


initial_mpi, initial_compiler = mpi_manifest(), compiler_manifest()
(OUTPUT / 'shared-mpi-sha256.json').write_text(json.dumps(initial_mpi, indent=2) + '\n')
(OUTPUT / 'shared-compiler-sha256.json').write_text(json.dumps(initial_compiler, indent=2) + '\n')
summary, common_runtime = {}, None
for name, flags in ARMS.items():
    destination = OUTPUT / name
    destination.mkdir(exist_ok=True)
    wrapper = destination / 'fortran-wrapper.py'
    # Each arm executes the exact same packaged compiler. Append only the
    # designated option, without editing any AMRVAC architecture/source file.
    wrapper.write_text('#!/usr/bin/python3\nimport json, os, sys\n'
                       + f'compiler = {COMPILER!r}\nextra = {flags!r}\n'
                       + f'log = {str(destination / "compiler-invocations.jsonl")!r}\n'
                       + 'args = sys.argv[1:] + extra\n'
                       + 'with open(log, "a") as stream:\n'
                       + '    stream.write(json.dumps([compiler] + args) + "\\n")\n'
                       + 'os.execv(compiler, [compiler] + args)\n')
    wrapper.chmod(0o755)
    env = dict(os.environ, FC=str(wrapper), OMPI_FC=str(wrapper),
               EXPECTED_GCC_MAJOR='14', EXPECTED_GCC_VERSION='14.3.0',
               LD_LIBRARY_PATH=f'{MPI}/lib:/usr/lib/x86_64-linux-gnu',
               CAPTURE_INSTRUCTION_AUDIT='1', COMPARISON_ARM=name)
    env.pop('LIBRARY_PATH', None)
    opts = record_options(destination, env)
    target = re.search(r'^\s*-march=\s+(\S+)', opts['target'], re.M)
    contract = re.search(r'^\s*-ffp-contract=\S*\s+(\S+)', opts['optimizer'], re.M)
    if not target or not contract:
        raise SystemExit('Unable to record effective target/contraction options')
    if name == 'x86-64' and target[1] != 'x86-64':
        raise SystemExit('Target override ineffective')
    if name == 'contract-off' and contract[1] != 'off':
        raise SystemExit('Contraction override ineffective')
    rc = subprocess.run([sys.executable, str(ROOT / '.github/scripts/ubuntu_comparison.py'),
                         str(ROOT / f'source-{name}'), str(destination)], env=env).returncode
    result_file = destination / 'results.json'
    records = json.loads(result_file.read_text()) if result_file.exists() else []
    complete = len(records) == 12 and all(r['status'] in ('PASS', 'NUMERICAL_MISMATCH') for r in records)
    if not complete or rc not in (0, 1):
        raise SystemExit(f'Incomplete arm {name}; inspect build/runtime logs')
    if (destination / 'source-diff.txt').read_bytes():
        raise SystemExit('AMRVAC source changed')
    manifests = {}
    for log in sorted(destination.glob('cases/**/linked-libraries.txt')):
        manifest = linked_manifest(log)
        if not any('/libgfortran.so.' in p for p in manifest) or not any(str(MPI) in p for p in manifest):
            raise SystemExit(f'Unexpected compiler/MPI runtime: {log}')
        manifests[str(log.relative_to(destination))] = manifest
        if common_runtime is None:
            common_runtime = manifest
        if manifest != common_runtime:
            raise SystemExit(f'Runtime library paths or hashes changed: {log}')
    if len(manifests) != 12:
        raise SystemExit('Missing linkage evidence')
    (destination / 'runtime-library-sha256.json').write_text(json.dumps(manifests, indent=2) + '\n')
    if mpi_manifest() != initial_mpi or compiler_manifest() != initial_compiler:
        raise SystemExit('MPI or compiler installation changed')
    invocations = [json.loads(l) for l in (destination / 'compiler-invocations.jsonl').read_text().splitlines()]
    if any(a[0] != COMPILER or (flags and a[-len(flags):] != flags) for a in invocations):
        raise SystemExit('Compiler invocation provenance mismatch')
    summary[name] = dict(flags=flags, effective_march=target[1], effective_fp_contract=contract[1],
                         returncode=rc, cases=len(records),
                         passed=sum(r['status'] == 'PASS' for r in records),
                         numerical_mismatch=sum(r['status'] == 'NUMERICAL_MISMATCH' for r in records))
    (OUTPUT / 'arm-progress.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('FLAG_ARM_COMPLETE ' + json.dumps({name: summary[name]}), flush=True)

result = dict(arms=summary, shared_mpi_unchanged=True, shared_compiler_unchanged=True,
              identical_runtime_paths_and_hashes=True, runtime_libraries=common_runtime)
(OUTPUT / 'flag-comparison.json').write_text(json.dumps(result, indent=2) + '\n')
with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as stream:
    stream.write('\n## Controlled compiler-flag comparison\n\n```json\n' + json.dumps(result, indent=2) + '\n```\n')
print(json.dumps(result, indent=2), flush=True)
sys.exit(1 if any(r['numerical_mismatch'] for r in summary.values()) else 0)

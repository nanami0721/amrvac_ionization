#!/usr/bin/env python3
"""Vary only the exact Fortran compiler on one host; audit shared binaries."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

root = Path(os.environ['GITHUB_WORKSPACE'])
temp = Path(os.environ['RUNNER_TEMP'])
output = root / 'results'
mpi = temp / 'comparison-openmpi'
versions = ('14.2.0', '14.3.0')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mpi_manifest():
    return {str(p.relative_to(mpi)): sha(p) for p in sorted(mpi.rglob('*')) if p.is_file()}


def linked_manifest(log):
    text = log.read_text()
    if 'not found' in text:
        raise RuntimeError(f'Unresolved dynamic library: {log}')
    paths = re.findall(r'(?:=>\s+|^\s*)(/\S+)', text, re.M)
    if not paths:
        raise RuntimeError(f'No dynamic libraries captured: {log}')
    return {str(Path(s).resolve()): sha(Path(s)) for s in paths}


baseline_mpi = mpi_manifest()
(output / 'shared-mpi-sha256.json').write_text(json.dumps(baseline_mpi, indent=2) + '\n')
all_results = {}
shared_manifest = None
for version in versions:
    compiler = temp / f'exact-gcc-{version}' / 'bin/gfortran'
    env = dict(os.environ, FC=str(compiler), OMPI_FC=str(compiler),
               EXPECTED_GCC_MAJOR='14', EXPECTED_GCC_VERSION=version,
               LD_LIBRARY_PATH=f'/usr/lib/x86_64-linux-gnu:{mpi}/lib',
               LIBRARY_PATH='/usr/lib/x86_64-linux-gnu')
    destination = output / f'gcc-{version}'
    destination.mkdir(exist_ok=True)
    rc = subprocess.run([sys.executable, str(root / '.github/scripts/ubuntu_comparison.py'),
                         str(root / f'source-{version}'), str(destination)], env=env).returncode
    result_file = destination / 'results.json'
    records = json.loads(result_file.read_text()) if result_file.exists() else []
    # A numerical mismatch returns 1 but must not skip the other compiler.
    complete = len(records) == 12 and all(r['status'] in ('PASS', 'NUMERICAL_MISMATCH') for r in records)
    if not complete or rc not in (0, 1):
        raise SystemExit(f'Incomplete compiler arm {version}; inspect build/runtime logs')
    if (destination / 'source-diff.txt').read_bytes():
        raise SystemExit('Tracked AMRVAC source changed')
    manifests = {}
    for log in sorted(destination.glob('cases/**/linked-libraries.txt')):
        manifest = linked_manifest(log)
        if any('/exact-gcc-' in p for p in manifest):
            raise SystemExit(f'Private compiler runtime loaded: {log}')
        if not any('/libgfortran.so.' in p for p in manifest):
            raise SystemExit(f'Missing dynamic Fortran runtime: {log}')
        if not any(str(mpi) in p for p in manifest):
            raise SystemExit(f'Wrong MPI linkage: {log}')
        manifests[str(log.relative_to(destination))] = manifest
        if shared_manifest is None:
            shared_manifest = manifest
        if manifest != shared_manifest:
            raise SystemExit(f'Dynamic library paths/hashes changed: {log}')
    if len(manifests) != 12:
        raise SystemExit('Missing per-case linkage records')
    (destination / 'runtime-library-sha256.json').write_text(json.dumps(manifests, indent=2) + '\n')
    if mpi_manifest() != baseline_mpi:
        raise SystemExit('Shared MPI installation changed')
    all_results[version] = dict(returncode=rc, cases=len(records),
                               passed=sum(r['status'] == 'PASS' for r in records),
                               numerical_mismatch=sum(r['status'] == 'NUMERICAL_MISMATCH' for r in records))
    print('EXACT_ARM_COMPLETE ' + json.dumps(all_results[version]), flush=True)

summary = dict(arms=all_results, shared_mpi_unchanged=True,
               identical_runtime_paths_and_hashes=True, runtime_libraries=shared_manifest)
(output / 'exact-comparison.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2), flush=True)
with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as stream:
    stream.write('\n## Exact GCC comparison audit\n\n```json\n' + json.dumps(summary, indent=2) + '\n```\n')
sys.exit(1 if any(r['numerical_mismatch'] for r in all_results.values()) else 0)

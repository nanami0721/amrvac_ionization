#!/usr/bin/env python3
"""Run the twelve reported regressions against an unchanged source checkout."""
import hashlib
import json
import os
import re
from collections import Counter
from pathlib import Path
import shutil
import subprocess
import sys
import time

SOURCE_SHA = "91f5091258cd8948f3a189a53c1c19ee15eb1974"
CASES = [
    ("hd/Coleman_1D", "coleman3_FI_HHe.log"),
    ("hd/Coleman_1D", "coleman3_FI_pureH.log"),
    ("hd/Coleman_1D", "coleman3_state_HHe.log"),
    ("hd/Coleman_1D", "coleman4_FI_HHe.log"),
    ("hd/Coleman_1D", "coleman4_FI_pureH.log"),
    ("hd/Coleman_1D", "coleman4_LTE_HHe.log"),
    ("mhd/convection_2D", "cv_2d_2step_tvdmu_al.log"),
    ("mhd/ring_thermal_conduction_3D", "rtc_3d.log"),
    ("mhd/solar_atmosphere_3D", "solar_atm_3D_3step_hll_cada_ct.log"),
    ("mhd/solar_atmosphere_3D", "solar_atm_3D_3step_hlld_cada_ct.log"),
    ("mhd/Ambipolar1DWave", "ambi.log"),
    ("mhd/Ambipolar1DWave-B0split", "ambi.log"),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def classify(make_rc, compare_rc, finished):
    if make_rc != 0 or not finished:
        return "BUILD_OR_RUNTIME_ERROR"
    if compare_rc is None:
        return "MISSING_OUTPUT_OR_COMPARATOR"
    return "PASS" if compare_rc == 0 else "NUMERICAL_MISMATCH"


def instruction_audit(binary, saved):
    """Record emitted instructions, not speculative optimizer settings."""
    counts = Counter()
    samples = []
    function = None
    fma_count = 0
    process = subprocess.Popen(["objdump", "-d", "--no-show-raw-insn", str(binary)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in process.stdout:
        label = re.match(r"^[0-9a-f]+ <(.+)>:", line)
        if label:
            function = label[1]
        instruction = re.match(r"^\s*[0-9a-f]+:\s+(\S+)", line)
        if not instruction:
            continue
        opcode = instruction[1]
        counts[opcode] += 1
        if re.match(r"^v(?:fmadd|fmsub|fnmadd|fnmsub)", opcode):
            fma_count += 1
            if len(samples) < 80:
                samples.append(dict(function=function, instruction=line.strip()))
    error = process.stderr.read()
    if process.wait() != 0:
        raise RuntimeError(f"objdump failed: {error}")
    result = dict(binary_sha256=digest(binary), instruction_count=sum(counts.values()),
                  fma_instruction_count=fma_count, opcode_counts=dict(sorted(counts.items())),
                  fma_samples=samples)
    (saved / "instruction-audit.json").write_text(json.dumps(result, indent=2) + "\n")


def main():
    source, output = [Path(s).resolve() for s in sys.argv[1:3]]
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, AMRVAC_DIR=str(source), ARCH="debug",
               NUM_PROCS="2", OMP_NUM_THREADS="1")
    env.pop("MAKEFLAGS", None)
    env.pop("MFLAGS", None)
    sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    if sha != SOURCE_SHA:
        raise SystemExit(f"Wrong source revision: {sha}")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=source):
        raise SystemExit("Source checkout is not clean")
    for directory, name in CASES:
        if not (source / "tests" / directory / "correct_output" / name).is_file():
            raise SystemExit(f"Missing reference: {directory}/{name}")

    def run(args, cwd, logfile):
        with logfile.open("w") as stream:
            return subprocess.run(args, cwd=cwd, env=env, stdout=stream,
                                  stderr=subprocess.STDOUT).returncode

    commands = [
        ["uname", "-a"], ["cat", "/etc/os-release"], ["lscpu"],
        ["gfortran", "--version"], [env.get("FC", "gfortran"), "--version"],
        [env.get("CC", "gcc"), "--version"], ["mpif90", "--version"],
        ["mpif90", "--showme"], ["mpif90", "--showme:command"],
        ["mpirun", "--version"], ["ompi_info", "--version"],
        ["ldd", "--version"],
        ["dpkg-query", "-W", "gfortran*", "gcc-13*", "gcc-14*", "libgfortran*",
         "openmpi*", "libopenmpi*", "libc6*", "libpmix*", "libprrte*"],
    ]
    with (output / "environment.txt").open("w") as stream:
        stream.write(f"SOURCE_SHA={sha}\n")
        for key in ("ARCH", "NUM_PROCS", "OMP_NUM_THREADS", "FC", "CC", "OMPI_FC",
                    "COMPARISON_ARM", "EXPECTED_GCC_MAJOR", "EXPECTED_GCC_VERSION", "LD_LIBRARY_PATH", "LIBRARY_PATH", "HWLOC_LIBXML", "RUNNER_OS", "RUNNER_ARCH", "ImageOS", "ImageVersion"):
            stream.write(f"{key}={env.get(key, '<unset>')}\n")
        for command in commands:
            stream.write("\n$ " + " ".join(command) + "\n")
            stream.flush()
            try:
                rc = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT).returncode
            except FileNotFoundError:
                rc = 127
            stream.write(f"exit_code={rc}\n")

    # Check the backend actually selected by mpif90, not merely the FC variable.
    backend = subprocess.check_output(["mpif90", "--showme:command"], env=env, text=True).strip()
    version = subprocess.check_output(["mpif90", "-dumpfullversion"], env=env, text=True).strip()
    expected = env.get("EXPECTED_GCC_MAJOR")
    exact = env.get("EXPECTED_GCC_VERSION")
    check = dict(backend=backend, version=version, expected_major=expected, expected_version=exact,
                 compiler_matches=(expected is None or version.split(".")[0] == expected)
                 and (exact is None or version == exact))
    (output / "compiler-check.json").write_text(json.dumps(check, indent=2) + "\n")
    print("ACTUAL_COMPILER " + json.dumps(check), flush=True)
    if not check["compiler_matches"]:
        raise SystemExit("MPI wrapper selected the wrong compiler; no regression cases run")

    # A different compiler may be unable to read Ubuntu's prebuilt mpi.mod.
    # Validate the module and two-rank runtime before interpreting numerical tests.
    probe = output / "mpi-probe.f90"
    probe.write_text("""program amrvac_compat_check
  use mpi
  implicit none
  integer :: ierr, rank, nproc, local_value, total
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  local_value = rank + 1
  call MPI_Allreduce(local_value, total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (nproc /= 2 .or. total /= 3) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  if (rank == 0) print *, 'MPI_PROBE_OK', nproc, total
  call MPI_Finalize(ierr)
end program amrvac_compat_check
""")
    binary = output / "mpi-probe"
    probe_compile = run(["mpif90", "-O2", "-fcheck=all", "-ffpe-trap=invalid,zero", str(probe), "-o", str(binary)], output, output / "mpi-probe-build.log")
    probe_run = None
    if probe_compile == 0:
        probe_run = run(["mpirun", "-np", "2", str(binary)], output, output / "mpi-probe-run.log")
    check.update(mpi_probe_compile=probe_compile, mpi_probe_run=probe_run)
    (output / "compiler-check.json").write_text(json.dumps(check, indent=2) + "\n")
    if probe_compile != 0 or probe_run != 0:
        for log in ("mpi-probe-build.log", "mpi-probe-run.log"):
            if (output / log).exists():
                print((output / log).read_text(), flush=True)
        raise SystemExit("MPI/compiler compatibility probe failed; no regression cases run")
    print("MPI_COMPILER_PROBE_OK", flush=True)

    # Build each required dimensional library separately to avoid header races.
    for dim in ("1d1", "1d3", "2d2", "3d3"):
        rc = run(["make", "-j", "2", dim, "ARCH=debug"], source / "lib", output / f"build-{dim}.log")
        if rc:
            # Preserve the failed parallel attempt; retry serially once.
            rc = run(["make", "-j", "1", dim, "ARCH=debug"], source / "lib", output / f"build-{dim}-retry.log")
        if rc:
            raise SystemExit(f"Library build failed: {dim}")

    comparator = source / "tools/fortran/compare_logs"
    records = []
    for directory, name in CASES:
        cwd = source / "tests" / directory
        saved = output / "cases" / directory / name.removesuffix(".log")
        saved.mkdir(parents=True, exist_ok=True)
        actual, reference = cwd / name, cwd / "correct_output" / name
        # A failed build must never inherit an earlier case's run.log.
        actual.unlink(missing_ok=True)
        (cwd / "run.log").unlink(missing_ok=True)
        print(f"START {directory}/{name}", flush=True)
        started = time.monotonic()
        rc = run(["make", "-j", "1", "-f", "test.make", "hdr"], cwd, saved / "header.log")
        if rc == 0:
            rc = run(["make", "-j", "1", "-f", "test.make", name], cwd, saved / "make.log")
        compare_rc = None
        if actual.is_file() and comparator.is_file():
            compare_rc = run([str(comparator), "1.0e-5", "1.0e-8", str(actual), str(reference)],
                             cwd, saved / "comparison.log")
        finished = (cwd / "run.log").is_file() and "Finished AMRVAC" in (cwd / "run.log").read_text(errors="replace")
        for path, destination in ((actual, "actual.log"), (reference, "reference.log"),
                                  (cwd / "run.log", "run.log"), (cwd / "setup.log", "setup.log")):
            if path.is_file():
                shutil.copy2(path, saved / destination)
        if (cwd / "amrvac").is_file():
            run(["ldd", str(cwd / "amrvac")], cwd, saved / "linked-libraries.txt")
            if env.get("CAPTURE_INSTRUCTION_AUDIT") == "1":
                instruction_audit(cwd / "amrvac", saved)
                run(["nm", "-D", "--undefined-only", str(cwd / "amrvac")], cwd, saved / "external-symbols.txt")
        record = dict(directory=directory, log=name, make_returncode=rc,
                      compare_returncode=compare_rc, finished_amrvac=finished,
                      status=classify(rc, compare_rc, finished),
                      elapsed_seconds=round(time.monotonic() - started, 2),
                      actual_sha256=digest(actual), reference_sha256=digest(reference))
        records.append(record)
        (output / "results.json").write_text(json.dumps(records, indent=2) + "\n")
        print(json.dumps(record), flush=True)
    diff_rc = run(["git", "diff", "--exit-code", SOURCE_SHA], source, output / "source-diff.txt")
    summary = [f"### {env.get('COMPARISON_ARM', 'Regression cases')}", "", f"Source: `{sha}`", "", "| Case | Result |", "|---|---|"]
    summary += [f"| {r['directory']}/{r['log']} | {r['status']} |" for r in records]
    summary.append(f"\nTracked source unchanged: {diff_rc == 0}\n")
    (output / "summary.md").write_text("\n".join(summary))
    if env.get("GITHUB_STEP_SUMMARY"):
        with open(env["GITHUB_STEP_SUMMARY"], "a") as stream:
            stream.write("\n".join(summary))
    raise SystemExit(0 if diff_rc == 0 and all(r["status"] == "PASS" for r in records) else 1)


if __name__ == "__main__":
    main()

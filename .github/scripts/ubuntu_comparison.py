#!/usr/bin/env python3
"""Run the twelve reported regressions against an unchanged source checkout."""
import hashlib
import json
import os
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
        ["gfortran", "--version"], ["gfortran-14", "--version"],
        ["gcc-14", "--version"], ["mpif90", "--version"],
        ["mpif90", "--showme"], ["mpif90", "--showme:command"],
        ["mpirun", "--version"], ["ompi_info", "--version"],
        ["ldd", "--version"],
        ["dpkg-query", "-W", "gfortran*", "gcc-14*", "libgfortran*",
         "openmpi*", "libopenmpi*", "libc6*", "libpmix*", "libprrte*"],
    ]
    with (output / "environment.txt").open("w") as stream:
        stream.write(f"SOURCE_SHA={sha}\n")
        for key in ("ARCH", "NUM_PROCS", "OMP_NUM_THREADS", "FC", "CC", "OMPI_FC",
                    "RUNNER_OS", "RUNNER_ARCH", "ImageOS", "ImageVersion"):
            stream.write(f"{key}={env.get(key, '<unset>')}\n")
        for command in commands:
            stream.write("\n$ " + " ".join(command) + "\n")
            stream.flush()
            try:
                rc = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT).returncode
            except FileNotFoundError:
                rc = 127
            stream.write(f"exit_code={rc}\n")

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
        record = dict(directory=directory, log=name, make_returncode=rc,
                      compare_returncode=compare_rc, finished_amrvac=finished,
                      status=classify(rc, compare_rc, finished),
                      elapsed_seconds=round(time.monotonic() - started, 2),
                      actual_sha256=digest(actual), reference_sha256=digest(reference))
        records.append(record)
        (output / "results.json").write_text(json.dumps(records, indent=2) + "\n")
        print(json.dumps(record), flush=True)
    diff_rc = run(["git", "diff", "--exit-code", SOURCE_SHA], source, output / "source-diff.txt")
    summary = [f"Source: `{sha}`", "", "| Case | Result |", "|---|---|"]
    summary += [f"| {r['directory']}/{r['log']} | {r['status']} |" for r in records]
    summary.append(f"\nTracked source unchanged: {diff_rc == 0}\n")
    (output / "summary.md").write_text("\n".join(summary))
    if env.get("GITHUB_STEP_SUMMARY"):
        with open(env["GITHUB_STEP_SUMMARY"], "a") as stream:
            stream.write("\n".join(summary))
    raise SystemExit(0 if diff_rc == 0 and all(r["status"] == "PASS" for r in records) else 1)


if __name__ == "__main__":
    main()

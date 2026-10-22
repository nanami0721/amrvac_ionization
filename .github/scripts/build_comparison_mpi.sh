#!/usr/bin/env bash
set -euo pipefail
result_dir="$GITHUB_WORKSPACE/results"
mpi_prefix="$RUNNER_TEMP/comparison-openmpi"
mkdir -p "$result_dir" "$RUNNER_TEMP/comparison-mpi-build"
cd "$RUNNER_TEMP/comparison-mpi-build"
curl --fail --location --retry 3 \
  https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.10.tar.bz2 \
  -o openmpi-5.0.10.tar.bz2
echo '0acecc4fc218e5debdbcb8a41d182c6b0f1d29393015ed763b2a91d5d7374cc6  openmpi-5.0.10.tar.bz2' | sha256sum -c - | tee "$result_dir/mpi-source-checksum.log"
tar -xf openmpi-5.0.10.tar.bz2
cd openmpi-5.0.10
trap 'rc=$?; if [ "$rc" -ne 0 ]; then tail -n 70 "$result_dir"/mpi-*.log; fi' EXIT
./configure --prefix="$mpi_prefix" \
  --with-pmix=internal --with-prrte=internal \
  --with-hwloc=internal --with-libevent=internal \
  --without-ucx --without-libfabric --without-verbs --without-cuda \
  > "$result_dir/mpi-configure.log" 2>&1
make -j 4 > "$result_dir/mpi-build.log" 2>&1
make install > "$result_dir/mpi-install.log" 2>&1
echo "$mpi_prefix/bin" >> "$GITHUB_PATH"
echo "LD_LIBRARY_PATH=$mpi_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" >> "$GITHUB_ENV"
"$mpi_prefix/bin/ompi_info" --all > "$result_dir/mpi-build-info.txt"

#!/usr/bin/env bash
set -euo pipefail
result_dir="$GITHUB_WORKSPACE/results/compilers"
mkdir -p "$result_dir"
for version in 14.2.0 14.3.0; do
  case "$version" in
    14.2.0) checksum=932bdef0cda94bacedf452ab17f103c0cb511ff2cec55e9112fc0328cbf1d803b42595728ea7b200e0a057c03e85626f937012e49a7515bc5dd256b2bf4bc396 ;;
    14.3.0) checksum=cb4e3259640721bbd275c723fe4df53d12f9b1673afb3db274c22c6aa457865dccf2d6ea20b4fd4c591f6152e6d4b87516c402015900f06ce9d43af66d3b7a93 ;;
  esac
  build_root="$RUNNER_TEMP/build-exact-gcc-$version"
  prefix="$RUNNER_TEMP/exact-gcc-$version"
  logs="$result_dir/$version"
  mkdir -p "$build_root" "$logs"
  cd "$build_root"
  curl --fail --location --retry 3 "https://gcc.gnu.org/pub/gcc/releases/gcc-$version/gcc-$version.tar.xz" -o gcc.tar.xz
  echo "$checksum  gcc.tar.xz" | sha512sum -c - | tee "$logs/source-checksum.txt"
  printf '%s  %s\n' "$checksum" "gcc-$version.tar.xz" > "$logs/SHA512SUMS"
  tar -xf gcc.tar.xz
  mkdir build
  cd build
  trap 'rc=$?; if [ "$rc" -ne 0 ]; then tail -n 70 "$logs"/*.log; fi' EXIT
  # Both releases use the same host compiler, dependency packages and flags.
  CC=gcc-14 CXX=g++-14 CFLAGS=-O2 CXXFLAGS=-O2 \
    ../gcc-$version/configure --prefix="$prefix" \
    --enable-languages=c,fortran --disable-bootstrap --disable-multilib \
    --disable-libsanitizer --with-system-zlib > "$logs/configure.log" 2>&1
  cp config.log "$logs/config.log"
  make -j 4 > "$logs/build.log" 2>&1
  make install > "$logs/install.log" 2>&1
  "$prefix/bin/gfortran" -v > "$logs/version.txt" 2>&1
  test "$("$prefix/bin/gfortran" -dumpfullversion)" = "$version"
  sha256sum "$prefix/bin/gfortran" "$("$prefix/bin/gfortran" -print-prog-name=f951)" > "$logs/compiler-sha256.txt"
  # Discard only this script's disposable source/build tree to preserve disk.
  cd "$RUNNER_TEMP"
  rm -rf "$build_root"
  df -h > "$logs/disk-after.txt"
  echo "EXACT_GCC_READY $version"
done

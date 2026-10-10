#!/bin/bash
# Build stock, fixed, bug-2-only and bug-1-only EPA-ng v0.3.8 from ONE source tree / ONE cmake
# configuration (identical compiler flags; only the patched translation units are rebuilt).
# Needs flex+bison. Output: /opt/tools/epa/epa-ng-{stock,fix,bug2only,bug1only}
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
S=/opt/src/epa-ng; O=/opt/tools/epa; mkdir -p $O
[ -d $S ] || git clone -q --recursive --branch v0.3.8 https://github.com/pierrebarbera/epa-ng.git $S
cd $S
# split the fix into its two hunks
filterdiff -i '*pll_util.cpp' $HERE/epa-ng-fix.patch > /tmp/bug2.patch 2>/dev/null || \
  awk '/^diff --git/{p=/pll_util/} p' $HERE/epa-ng-fix.patch > /tmp/bug2.patch
awk '/^diff --git/{p=/file_io/} p' $HERE/epa-ng-fix.patch > /tmp/bug1.patch
build() { git checkout -q -- src; [ -n "$1" ] && for p in $1; do git apply $p; done; make -j4 > /tmp/epa_build_$2.log 2>&1; cp bin/epa-ng $O/epa-ng-$2; }
build "" stock
build "/tmp/bug2.patch" bug2only     # scaler-shift fix (accuracy)
build "/tmp/bug1.patch" bug1only     # keep SIMD attribute bits (speed)
build "/tmp/bug2.patch /tmp/bug1.patch" fix
git checkout -q -- src
md5sum $O/epa-ng-*

#!/bin/bash
# Build nuflow (incremental) and generate the closed-loop prediction sample.
#   scripts/run_loop_prediction.sh <output.root> [config]   (default config/minerva_me_fhc_numu_CH_loop.yaml)
# Run from anywhere; uses the repo this script lives in.  Fails on build or test failure.
set -euo pipefail
OUT=$(realpath -m "$1")
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CFG=${2:-config/minerva_me_fhc_numu_CH_loop.yaml}
cd "$DIR"
set +u; source setup.sh >/dev/null; set -u
[ -f build/CMakeCache.txt ] || cmake -S . -B build -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build build -j16 2>&1 | tail -2
./build/nuflow-tests >build/tests.log 2>&1 || { tail -30 build/tests.log; echo "nuflow-tests FAILED"; exit 1; }
tail -1 build/tests.log | grep -q "All tests passed" || true
mkdir -p "$(dirname "$OUT")"
./build/nuflow-gen --config "$CFG" -o "$OUT" 2>&1 | grep -v "events$" | tee "${OUT%.root}.log"

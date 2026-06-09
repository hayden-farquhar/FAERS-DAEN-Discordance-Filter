#!/usr/bin/env bash
# Finish the confirmatory sweep locally: 6 parallel shard-workers (one per
# physical core), single-threaded BLAS for determinism + no oversubscription.
# Each worker skips cells whose shard already exists, so only the 25 missing
# cells get computed. Workers do NOT finalise; combine is a separate step.
set -u
cd "$(dirname "$0")/src" || exit 1
RESULTS="$(dirname "$0")/results"
mkdir -p "$RESULTS"
pids=()
for k in 0 1 2 3 4 5; do
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    python3 -m simulation.sweep --shard-index "$k" --num-shards 6 \
    > "$RESULTS/local_worker_$k.log" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do
  wait "$p" || fail=1
done
echo "ALL 6 WORKERS EXITED (fail=$fail)"
exit "$fail"

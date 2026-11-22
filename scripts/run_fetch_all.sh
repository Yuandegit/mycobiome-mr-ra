#!/bin/bash
# Process each chromosome in a fully isolated Python process so any htslib
# buffer growth is released when the process exits. Retry each chromosome up
# to 4 times; skip chromosomes already fetched.
set -u
ROOT=/home/user/mycobiome_joint_mr
PY=/opt/python3.12/bin/python3
SCRIPT=$ROOT/scripts/04_fetch_eas_genotypes.py
OUTDIR=$ROOT/data/gwas_microbe/eas_chr

for c in $(seq 1 22); do
  out=$OUTDIR/chr$c.pkl
  if [ -f "$out" ]; then
    echo "skip chr$c (exists)"
    continue
  fi
  ok=0
  for attempt in 1 2 3 4; do
    echo "=== chr$c attempt $attempt ==="
    if $PY "$SCRIPT" chr$c; then
      if [ -f "$out" ]; then ok=1; break; fi
    fi
    sleep 3
  done
  if [ "$ok" -ne 1 ]; then echo "FAILED chr$c after retries"; fi
done
echo "ALL DONE"

#!/usr/bin/env bash

set -o errexit
set -o nounset

INDIR=$1
OUTDIR=$2
OUTLOG=$3
PCT_DEPTH=$4
MIN_PCT_ALIGNED=$5

NPROC=12


for path in $INDIR/*; do
  base=`basename $path`;
  filter_align.pl -in $path -out $OUTDIR/$base \
    -pct_depth $PCT_DEPTH -min_pct_aligned $MIN_PCT_ALIGNED -require_inform \
    -log $OUTLOG/$base.d$PCT_DEPTH.a$MIN_PCT_ALIGNED &

  # allow to execute up to $NPROC in parallel
  if [[ $(jobs -r -p | wc -l) -ge $NPROC ]]; then
    wait -n # wait for the next job to terminate
  fi
done
wait


#!/usr/bin/env bash

set -o errexit
set -o nounset

INDIR=$1
OUTDIR=$2
NPROC=$3

for filepath in $INDIR/*; do
  file=`basename $filepath`;
  echo "  Computing alignment, using program famsa, for file $file"
  famsa -t $NPROC $INDIR/$file $OUTDIR/$file
done
wait

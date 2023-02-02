#!/usr/bin/env bash

set -o errexit
set -o nounset

INDIR=$1
OUTDIR=$2

NPROC=$3

for filepath in $INDIR/*; do
  file=`basename $filepath`;
  muscle -in $filepath -out $OUTDIR/$file -quiet &
  #clustalw2 -INFILE=$filepath -ALIGN -OUTPUT=FASTA -OUTFILE=$OUTDIR/$file -QUIET &

  # allow to execute up to $NPROC in parallel
  if [[ $(jobs -r -p | wc -l) -ge $NPROC ]]; then
    wait -n # waits for the next job to terminate
  fi
done
wait



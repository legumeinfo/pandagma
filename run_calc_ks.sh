#!/usr/bin/env bash

# NAME:  run_calc_ks.sh - Call calc_ks_from_dag.pl, looping around provided list of 
#            DAGChainer .aligncoords files
#
# SYNOPSIS:   ./run_calc_ks.sh WORKDIR THREADS

set -o errexit -o nounset 

WORKDIR=$1
THREADS=$2

if [ ! -d $WORKDIR/04_dag ]; then
  echo "Expected DAGChainer dir doesn't exist at $WORKDIR/04_dag/ ; exiting."
  exit
fi

if [ ! -f $WORKDIR/02_all_main_cds.fna ]; then
  echo "Expected fasta files don't exist at $WORKDIR/*_cds.fna ; exiting."
  exit
fi

if [ ! -d $WORKDIR/05_kaksout ]; then
  echo "creating output directory $WORKDIR/05_kaksout"
  mkdir -p $WORKDIR/05_kaksout
fi

for DAGFILE in $WORKDIR/04_dag/*aligncoords; do
  base=`basename $DAGFILE _matches.tsv.aligncoords`
  echo "WORKING ON $base"
  bin/calc_ks_from_dag.pl $WORKDIR/*_cds.fna -dagin $DAGFILE -report_out $WORKDIR/05_kaksout/$base.rptout 1> /dev/null 2> /dev/null & 
  if [[ $(jobs -r -p | wc -l) -ge ${THREADS} ]]; then wait -n; fi
done



#!/usr/bin/env bash

#!/bin/sh
# NAME:   run_some_fam.sh - Call pandagma-fam.sh with some sub-commands.
#
# SYNOPSIS:   ./run_some_fam.sh config/YOUR_CONFIG.conf
#     
# Just a convenience if a few commands are being run (rather than the default "all").

set -o errexit -o nounset 

CONFIG=$1
#./pandagma-fam.sh -c $CONFIG -s run_ingest
#./pandagma-fam.sh -c $CONFIG -s run_mmseqs
#./pandagma-fam.sh -c $CONFIG -s run_filter
#./pandagma-fam.sh -c $CONFIG -s run_dagchainer
./pandagma-fam.sh -c $CONFIG -s run_mcl
./pandagma-fam.sh -c $CONFIG -s run_consense
./pandagma-fam.sh -c $CONFIG -s run_cluster_rest
./pandagma-fam.sh -c $CONFIG -s run_add_extra
./pandagma-fam.sh -c $CONFIG -s run_align
./pandagma-fam.sh -c $CONFIG -s run_model_and_trim
./pandagma-fam.sh -c $CONFIG -s run_calc_trees
./pandagma-fam.sh -c $CONFIG -s run_summarize
#./pandagma-fam.sh -c $CONFIG -s run_clean
#./pandagma-fam.sh -c $CONFIG -s run_ReallyClean


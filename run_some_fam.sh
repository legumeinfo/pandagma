#!/usr/bin/env bash

# NAME:   run_some_fam.sh - Call pandagma-fam.sh with some sub-commands.
#
# SYNOPSIS:   ./run_some_fam.sh config/YOUR_CONFIG.conf
#     
# Just a convenience if a few commands are being run (rather than the default "all").

set -o errexit -o nounset 

CONFIG=$1

# Stage 1: ingest through ks_calc
./pandagma-fam.sh -c $CONFIG -s ingest
./pandagma-fam.sh -c $CONFIG -s mmseqs
./pandagma-fam.sh -c $CONFIG -s filter
./pandagma-fam.sh -c $CONFIG -s dagchainer
./pandagma-fam.sh -c $CONFIG -s ks_calc  # Optional but recommended. Time-consuming.

# Stage 2: ks_filter through summarize
./pandagma-fam.sh -c $CONFIG -s ks_filter  # Optional but recommended. Provide the file ks_cutoffs.tsv .
./pandagma-fam.sh -c $CONFIG -s mcl
./pandagma-fam.sh -c $CONFIG -s consense
./pandagma-fam.sh -c $CONFIG -s cluster_rest
./pandagma-fam.sh -c $CONFIG -s add_extra
./pandagma-fam.sh -c $CONFIG -s align
./pandagma-fam.sh -c $CONFIG -s model_and_trim
./pandagma-fam.sh -c $CONFIG -s calc_trees
./pandagma-fam.sh -c $CONFIG -s summarize

# Optional steps to cleanup the work directory
#./pandagma-fam.sh -c $CONFIG -s clean
#./pandagma-fam.sh -c $CONFIG -s ReallyClean


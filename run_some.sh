#!/usr/bin/env bash

#!/bin/sh
# NAME:   run_some.sh - Call pandagma.sh with some sub-commands.
#
# SYNOPSIS:   ./run_some.sh config/YOUR_CONFIG.conf
#     
# Just a convenience if a few commands are being run (rather than the default "all").
# To add further "extra" sequences to previously-generated clusters, use the following commands:
#   ingest, add_extra, pick_exemplars, filter_to_pctile, order_and_name, calc_chr_pairs, summarize, clean

set -o errexit -o nounset 

CONFIG=$1

./pandagma.sh -c $CONFIG -s ingest
#./pandagma.sh -c $CONFIG -s mmseqs
#./pandagma.sh -c $CONFIG -s filter
#./pandagma.sh -c $CONFIG -s dagchainer
#./pandagma.sh -c $CONFIG -s mcl
#./pandagma.sh -c $CONFIG -s consense
#./pandagma.sh -c $CONFIG -s cluster_rest
./pandagma.sh -c $CONFIG -s add_extra
./pandagma.sh -c $CONFIG -s pick_exemplars
./pandagma.sh -c $CONFIG -s filter_to_pctile
./pandagma.sh -c $CONFIG -s order_and_name
./pandagma.sh -c $CONFIG -s calc_chr_pairs
./pandagma.sh -c $CONFIG -s summarize
./pandagma.sh -c $CONFIG -s clean
#./pandagma.sh -c $CONFIG -s ReallyClean


#!/usr/bin/env bash

#!/bin/sh
# NAME:   run_some_pan.sh - Call pandagma-pan.sh with some sub-commands.
#
# SYNOPSIS:   ./run_some_pan.sh config/YOUR_CONFIG.conf
#     
# Just a convenience if a few commands are being run (rather than the default "all").
# To add further "extra" sequences to previously-generated clusters, use the following commands:
#   ingest, add_extra, pick_exemplars, filter_to_pctile, order_and_name, calc_chr_pairs, summarize, clean

set -o errexit -o nounset 

CONFIG=$1

./pandagma-pan.sh -c $CONFIG -s ingest
./pandagma-pan.sh -c $CONFIG -s mmseqs
./pandagma-pan.sh -c $CONFIG -s filter
./pandagma-pan.sh -c $CONFIG -s dagchainer
./pandagma-pan.sh -c $CONFIG -s mcl
./pandagma-pan.sh -c $CONFIG -s consense
./pandagma-pan.sh -c $CONFIG -s cluster_rest
./pandagma-pan.sh -c $CONFIG -s add_extra
./pandagma-pan.sh -c $CONFIG -s pick_exemplars
./pandagma-pan.sh -c $CONFIG -s filter_to_pctile
./pandagma-pan.sh -c $CONFIG -s order_and_name
./pandagma-pan.sh -c $CONFIG -s calc_chr_pairs
./pandagma-pan.sh -c $CONFIG -s summarize
#./pandagma-pan.sh -c $CONFIG -s clean
#./pandagma-pan.sh -c $CONFIG -s ReallyClean


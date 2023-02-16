#!/usr/bin/env bash

#!/bin/sh
# NAME
#     run_some.sh - Call pandagma.sh with some sub-commands.
#     Just a convenience if a few commands are being run (rather than the default "all").
#
# SYNOPSIS
#    ./run_some.sh config/YOUR_CONFIG.conf

set -o errexit -o nounset 

CONFIG=$1

#./pandagma.sh -c $CONFIG -s ingest
#./pandagma.sh -c $CONFIG -s add_extra
#./pandagma.sh -c $CONFIG -s pick_exemplars
#./pandagma.sh -c $CONFIG -s filter_to_pctile
./pandagma.sh -c $CONFIG -s name_pangenes
./pandagma.sh -c $CONFIG -s calc_chr_pairs
./pandagma.sh -c $CONFIG -s summarize
./pandagma.sh -c $CONFIG -s clean


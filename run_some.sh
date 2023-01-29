#!/usr/bin/env bash

#./pandagma.sh -c config/Vigna_7_5.conf -s ingest
./pandagma.sh -c config/Vigna_7_5.conf -s add_extra
./pandagma.sh -c config/Vigna_7_5.conf -s pick_exemplars
./pandagma.sh -c config/Vigna_7_5.conf -s filter_to_core
./pandagma.sh -c config/Vigna_7_5.conf -s name_pangenes
./pandagma.sh -c config/Vigna_7_5.conf -s calc_chr_pairs
./pandagma.sh -c config/Vigna_7_5.conf -s summarize


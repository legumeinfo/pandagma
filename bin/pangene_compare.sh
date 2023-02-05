#!/usr/bin/env bash
# NAME
#     pangene_compare.sh  -- Compare composition of two pangene sets. 
#        Operates on files with a two-column hash format, with the set ID in the first column and gene IDs in the second.
#
# SYNOPSIS
#    pangene_compare.sh FILE1.hsh.tsv FILE2.hsh.tsv
#
# OUTPUT
#    is to STDOUT, and consists of counts of corresponding panIDs from the two sets.
#
# AUTHOR
#     Steven Cannon <steven.cannon@usda.gov>

set -o errexit
set -o nounset

PG1=$1
PG2=$2

LC_ALL=C join -1 2 -2 2 <(LC_ALL=C sort -k2,2 $PG1) \
                        <(LC_ALL=C sort -k2,2 $PG2) |
           awk '{print $2 "\t" $3}' | sort -k1,1 -k2,2 | uniq -c |
           awk -v OFS="\t" '{print $1, $2, $3}' | sort -k2,2 -k1nr,1nr |
           awk -v ONE=$PG1 -v TWO=$PG2 'BEGIN{print "#count" "\t" ONE "\t" TWO} {print}' 


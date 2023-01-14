#!/usr/bin/awk -f
#
# NAME
#   calc_avg_pangene_order.awk - Given input consisting of two columns with the
#      first column having identifiers that occur one or multiple times, and a
#      second column having numbers (e.g. ordinal positions), calculate
#      the average of the numbers per identifier.
# 
# SYNOPSIS
#   calc_avg_pangene_order.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing data in two columns and sorted by col 1.

BEGIN{OFS="\t"}

NR==1 { ct = 1; prev = $1; sum=$2 }
NR>1 && $1 == prev { ct++; sum+=$2 }
NR>1 && $1 != prev { printf "%s\t%d\t%.2f\n", $1, ct, sum/ct; 
                     ct = 1; sum=$2; prev = $1
                   }

END{ printf "%s\t%d\t%.2f\n", $1, ct, sum/ct }


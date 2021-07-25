#!/usr/bin/awk -f
#
# NAME
#   hash_to_rows_by_1st_col.awk - Given input consisting of two columns with the
#      first column having identifiers that occur one or multiple times, and a
#      second column having unique identifiers, pull the IDs in the second column
#      into a single row anchored by the ID in the first column.
# 
# SYNOPSIS
#   ./hash_to_rows_by_1st_col.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing data in two columns and sorted by col 1.

BEGIN{ORS=""; OFS="\t"}

NR==1 { print $1 "\t" $2; count = 1; prev = $1 }
NR>1 && $1 == prev { print "\t" $2; count++ }
NR>1 && $1 != prev { print "\n" $1 "\t" $2; count = 1; prev = $1 }

END{print "\n"}


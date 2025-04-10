#!/usr/bin/awk -f
#
# NAME
#   calc_avg_bitscore_per_fam_hit.awk - Given input consisting of tabular BLAST-like output, with the
#      first column having gene IDs, the second having gene family IDs, and the twelfth having bit scores,
#      sorted by gene ID and then (numeric reverse) by bit score,
#      calculate the average bit scores per gene ID and gene family, reporting one line per gene ID and gene family
#      ... with the values of the first 11 columns for a row coming from the first query-target line.
# 
# SYNOPSIS
#   calc_avg_bitscore_per_fam_hit.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular, blast-like data
# NOTE  The construct "NF--; line=$0" stores all but the last column (bitscore). The average bitscore is swapped in later.

BEGIN{ FS=OFS="\t" }

NR==1 { ct = 1; prev1 = $1; prev2 = $2; sum=$12; 
        NF--; line=$0 }
NR>1 && $1 == prev1 && $2 == prev2 { ct++; sum+=$12 }
NR>1 && ($1 != prev1 || $2 != prev2) { printf "%s\t%.2f\n", line, sum/ct; 
                     ct = 1; sum=$12; prev1 = $1; prev2 = $2;
                     NF--; line=$0
                   }

END{ printf "%s\t%.2f\n", line, sum/ct }


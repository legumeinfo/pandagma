#!/usr/bin/awk -f
#
# NAME
#   top_hmmscan_line.awk - Filter tabular hmmscan output to top hit per query.
# 
# SYNOPSIS
#   ./top_hmmscan_line.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular hmmscan output (or similar)

BEGIN { MAX = 1 }

$3 == prev && count < MAX { print; count++ }
$3 != prev { print; count = 1; prev = $3 }


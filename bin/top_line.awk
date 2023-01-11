#!/usr/bin/awk -f
#
# NAME
#   top_line.awk - Filter tabular blast output to top hit per query.
# 
# SYNOPSIS
#   ./top_line.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular blast output (or similar); the important
#         field is 0 (query ID in blast output).
#
# EQUIVALENT TO
#   awk '$1==prev && ct<2 {print; ct++} $1!=prev {print; ct=1; prev=$1}' FILE(s)

BEGIN { MAX = 1 } 

$1 == prev && count < MAX { print; count++ }
$1 != prev { print; count = 1; prev = $1 }

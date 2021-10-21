#!/usr/bin/awk -f
#
# NAME
#   pile_against_col1.awk - Process a two-column file, pulling elements in the second column
#                           into one row when they have a common first element. Example:
#        from    ABC  01
#                ABC  02
#                DEF  03
#        to      ABC  01  02
#                DEF  03
# 
# SYNOPSIS
#   cat INPUT | pile_against_col1.awk
#     
# OPERANDS
#     STDIN or input file; input consisting of two columns, sorted on the first column.

NR==1 { prev=$1; cat=$0 } 
NR>1 && $1==prev { prev = $1; cat = cat "\t" $2 } 
NR>1 && $1!=prev { print cat; cat=$0; prev=$1 }
END{ print cat }


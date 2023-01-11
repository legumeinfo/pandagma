#!/usr/bin/awk -f
#
# NAME
#   name_ordered_genes.awk 
#     Process input consisting of four columns:
#         ID chr start end     ... sorted by chr then start (-k2,2 -k3n,3n)
#     Generate new gene IDs with the form chr_ord00_ID
#          
# SYNOPSIS
#   cat INPUT | name_ordered_genes.awk
#     
# OPERANDS
#     STDIN or input file; input consisting of four columns, sorted on col2 then col3-numeric

$1~/^#/ {}
$1!~/^#/ && $2==prevchr { ct++; prevchr=$2; printf("%s\t%s_%04d00\t%d\t%d\n", $1, $2, ct, $3, $4) } 
$1!~/^#/ && $2!=prevchr { ct=1; prevchr=$2; printf("%s\t%s_%04d00\t%d\t%d\n", $1, $2, ct, $3, $4) }


#Glycine.pan3.x20_256700
#Medicago.pan1.x8_000100
#Arachis.pan1.x10_001200
#
#Glycine.pan3.chr20_256700
#Medicago.pan1.chr8_000100
#Arachis.pan1.chr10_001200
#


#!/bin/sh
# NAME
#     n50stats.sh - print some basic stats about sequences in a multifasta file
#
# SYNOPSIS
#    cat INPUT_FILE | n50stats.sh 
#
# AUTHOR
#     Steven Cannon <steven.cannon@ars.usda.gov>

awk 'BEGIN { ORS="" } \
     /^>/ && NR==1 { print $0"\n" } \
     /^>/ && NR!=1 { print "\n"$0"\n" } \
     /^[^>]/ { print } \
     END { print "\n" }' |
awk '/^[^>]/ {print length($1); tot_bp+=length($1) } END { print "bases: " tot_bp "\n" }' \
 | sort -n \
 | awk 'BEGIN { min=99999999999999 }
        /bases:/ { N50_ct = $2/2; bases=$2 } 
        /^[0-9]/ { 
           ct++; sum+=$1; if ( sum >= N50_ct && !printed ) { N50=$1; printed = 1 } 
           min = (min < $1) ? min : $1
        } 
        END { 
          ave=sum/ct;
          print "  seqs  min  max  N50  ave  sum";
          printf("  %d  %d  %d  %d  %.1f  %d ", ct, min, $1, N50, ave, bases);
          print "";
        }'



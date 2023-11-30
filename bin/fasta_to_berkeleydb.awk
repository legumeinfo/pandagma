#!/usr/bin/awk -f
# NAME
#     fasta_to_berkeleydb.awk - remove line returns from fasta files, and report in alternating format
#        to permit creation of a berkeleydb. The output format is:
#          seqID1
#          seq1
#          seqID2
#          seq2
#
# SYNOPSIS
#     fasta_to_berkeleydb.awk INPUT_FILE [INPUT FILE...] | db_load -T -t hash FASTA.db

BEGIN { ORS="" }

/^>/ && NR==1 { print substr($1,2) "\n" }

/^>/ && NR!=1 { print "\n" substr($1,2) "\n" }

/^[^>]/ { print }

END { print "\n" }


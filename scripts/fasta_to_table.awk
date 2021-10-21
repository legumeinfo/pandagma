#!/usr/bin/awk -f
# NAME
#     fasta_to_table.awk - remove line returns from fasta files, and report in table format: ID, seq
#
# SYNOPSIS
#     fasta_to_table.awk INPUT_FILE [INPUT FILE...]

BEGIN { ORS="" }

/^>/ && NR==1 { print substr($0,2) "\t" }

/^>/ && NR!=1 { print "\n" substr($0,2) "\t" }

/^[^>]/ { print }

END { print "\n" }


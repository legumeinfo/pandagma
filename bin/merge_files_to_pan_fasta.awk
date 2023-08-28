#!/usr/bin/awk -f
#
# NAME
#   merge_files_to_pan_fasta.awk - Combine sets of fasta files from one directory into a new
#                                  merged fasta file with the filename prefixed to the sequence ID, e.g.
#                                    >pan00001__gene_ID
# SYNOPSIS
#   merge_files_to_pan_fasta.awk FASTA_DIR/*

function add_file() {
  fasta_file = path[nf]
}

{
  nf = split(FILENAME, path, "/")
  add_file()
}

/^>/ { printf(">%s__%s\n", fasta_file, substr($1,2)) }

/^[^>]/ { print $1 }

END { add_file() } 


# Functionally equivalent to the following, but faster:
#   for path in 07_pan_fasta/*; do
#     pan_file=`basename $path`
#     cat $path | awk -v panID=$pan_file ' $1~/^>/ {print ">" panID "__" substr($0,2) }
#                       $1!~/^>/ {print $1} ' >> 07_pan_fasta_cds.fna
#   done


#!/usr/bin/awk -f
#
# NAME
#   augment_cluster_sets.awk - Combine sets of files from two directories,
#                              EXTRA_DIR and FASTA_DIR, to
#                              produce a new gene-family "cluster" file.
# SYNOPSIS
#   augment_cluster_sets.awk leftovers_dir=EXTRA_DIR FASTA_DIR/*

function add_leftovers() {
  leftovers_file = leftovers_dir "/" path[nf]
  while(getline fasta_line < leftovers_file == 1)
      if (fasta_line ~ /^>/) printf("\t%s", substr(fasta_line,2))
  close(leftovers_file)
  printf("\n")
}
FNR == 1 {
  if (NR != 1) add_leftovers() # previous pan cluster
  nf = split(FILENAME, path, "/")
  printf("%s", path[nf])
}
/^>/ { printf("\t%s", substr($1,2)) }
END { add_leftovers() } 

## Functionally equivalent to the following (but much faster):
#  for path in 07_pan_fasta/*; do
#    file=`basename $path`
#    if [[ -f 11_pan_leftovers/$file ]]; then
#      cat <(awk -v ORS="" '$1~/^>/ {print substr($1,2) "\t"}' $path) \
#          <(awk -v ORS="" '$1~/^>/ {print substr($1,2) "\t"} END{print "\n"}' 11_pan_leftovers/$file) |
#            awk -v CLUST=$file '{print CLUST "\t" $0}'
#    else
#      awk -v ORS="" '$1~/^>/ {print substr($1,2) "\t"} END{print "\n"}' $path |
#        awk -v CLUST=$file '{print CLUST "\t" $0}'
#    fi | sed 's/\t$//'
#  done > 12_syn_pan_aug.clust.tsv


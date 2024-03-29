#!/usr/bin/env bash

set -o errexit
set -o nounset
# NAME
#   gene_trans_db_load.sh -- A little script to create berkeleydb files for use by gene_translate.pl
#
# SYNOPSIS
#   gene_trans_db_load.sh $pan_hash $pan_clust
#
#   This script takes two positional arguments, with paths to two files generated by
#   pandagma-pan.sh: one for a (gzipped) hash file, suffixes hsh.tsv.gz;
#   and one for a (gzipped) "cluster" file, suffixes clust.tsv.gz.
#   The BerkeleyDB utility "db_load" must be in $PATH.
#
# AUTHOR
#   Steven Cannon <steven.cannon@usda.gov>

pan_hash=$1   # Glycine.pan4.RK4P.hsh.tsv.gz
pan_clust=$2  # Glycine.pan4.RK4P.clust.tsv.gz

# For db_load, print the key-value pairs on alternate lines (key odd, value even)
zcat < $pan_hash |
  awk '{print $1 "\t" tolower($2)}' |
  perl -pe 's/(\S+)\t(\S+)\.\D*\d+/$2\n$1/' |
  db_load -T -t hash gn_to_pan.db &

zcat < $pan_clust |
  perl -lane 'print $F[0];
              @genes = ();
              for $gene (@F[1..scalar(@F)-1]){$gene =~ s/(\S+)\.\D*\d+$/$1/; push(@genes, $gene)};
              print join (",", @genes)' |
    db_load -T -t hash pan_to_gn.db &

wait
echo "Done generating BDB files gn_to_pan.db and pan_to_gn.db"


#!/usr/bin/env bash
# NAME
#     compare_pangene_runs.sh  
#
# SYNOPSIS
#    compare_pangene_runs.sh out_dir1 out_dir2 out_dir3
#
# OUTPUT
#    
#
# AUTHOR
#     Steven Cannon <steven.cannon@usda.gov>

set -o errexit
set -o nounset

PG1=$1
PG2=$2
PG3=$3

LC_ALL=C join -1 2 -2 2 <(LC_ALL=C sort -k2,2 $PG1/18_syn_pan_aug_extra.hsh.tsv) \
                        <(LC_ALL=C sort -k2,2 $PG2/18_syn_pan_aug_extra.hsh.tsv) |
         awk '{print $2 "\t" $3}' | sort -k1,1 -k2,2 | uniq -c |
         awk -v OFS="\t" '{print $1, $2, $3, $1}' | sort -k2,2 -k1nr,1nr |
         awk -v -OFS="\t" '{print $2, $3, $1}' | top_line.awk | 
         awk '{print $1 "\t" $2}' > $PG3/map_PG1_PG2.tsv

# Replace the new pan-set IDs with the previous ones in the tsv files
for filepath in $PG1/*.hsh.tsv $PG2/*.clust.tsv $PG1/*counts.tsv; do
  base=`basename $filepath`
  echo $base
  echo "hash_into_table_id.pl $PG3/map_PG1_PG2.tsv -table $filepath > $PG3/$base"
  hash_into_table_id.pl $PG3/map_PG1_PG2.tsv -table $filepath > $PG3/$base
  echo
done

# Replace the new pan-set IDs with the previous ones in the fasta files that use the panID as sequence ID
for filepath in $PG1/21_*.f?a $PG1/22_*.f?a ; do
  base=`basename $filepath`
  echo $base
  echo "hash_into_fasta_id.pl -hash $PG3/map_PG1_PG2.tsv -fasta $filepath > $PG3/$base"
  hash_into_fasta_id.pl -hash $PG3/map_PG1_PG2.tsv -fasta $filepath > $PG3/$base
  echo
done

# Replace the new pan-set IDs with the previous ones in the fasta files with position_panID as sequence ID
#      >Glycine.pan4.chr01_000100_pan47277
# Make a temporary file first that has IDs of the form
#      >pan47277 Glycine.pan4.chr01_000100
for filepath in $PG1/23_*_posn_*.f?a; do
  base=`basename $filepath`
  echo $base
  cat $filepath | perl -pe 's/>(\S+)_(pan\d+) (.+)/>$2 $1 $3/' > $PG3/tmp.$base
  hash_into_fasta_id.pl -hash $PG3/map_PG1_PG2.tsv -fasta $PG3/tmp.$base |
    perl -pe 's/>(\S+) (\S+) (.+)/>$2_$1 $3/' > $PG3/$base
  rm $PG3/tmp.$base
  echo
done




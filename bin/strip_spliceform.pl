#!/usr/bin/env perl

# In multicolumn input in which the third column is a mRNA ID with splice variant,
# return the data with the splice variant string stripped from the third column.
# Intended use: in tabular "gene family assignment" data derived from hmmsearch,
# produce a gene ID in the third column.
# The regex below should handle geneID.1, geneID.m1, geneID.mRNA1, geneID-T1

while (<>){
  my @F = split(/\t/, $_);
  if ( $F[2] =~ /(.+)\.\w\d+$/ ||
       $F[2] =~ /(.+\w\w+)\.\d{1,2}\.\d{1,2}$/ ||
       $F[2] =~ /(.+)\.\d+$/ ||
       $F[2] =~ /(.+)-T\d+$/ ||
       $F[2] =~ /(.+)\.mRNA\d+$/ ||
       $F[2] =~ /^([^.]+\.[^.]+\.[^.]+\.[^.]+\.[^.]+)$/ ){
         $F[2]=$1; print join("\t",@F)
       }
  else { $F[2]="XX"; print join ("\t",@F) } # print a flag; this data will need special handling
}


#!/usr/bin/env perl

use warnings;
use strict;

my $usage = <<EOS;
  Synopsis: cat FILE.clust.tsv | reshape_mcl_to_hash.pl

  Reshape from mcl output format (clustered IDs on one line, with first item being a cluster ID) ...
    ID_1 geneA geneB geneC
    ID_2 geneD
  to a hash format ...
    ID_1 geneA
    ID_1 geneB
    ID_1 geneC
    ID_2 geneD
EOS

while (<>){
  chomp;
  my @items = split(/\t/, $_);
  my $ct = scalar(@items); 
  for my $i ( 1..$ct-1 ){ print "$items[0]\t$items[$i]\n" }
}


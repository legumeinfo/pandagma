#!/usr/bin/env perl
use strict;
use warnings;

my %type_collate = (
  region => 0.00,
  gene => 0.01,
  pseudogene => 0.02,
  transposable_element => 0.021,
  transposable_element_gene => 0.022,
  mRNA => 0.03,
  transcript => 0.04,
  ncRNA => 0.05,
  lnc_RNA => 0.06,
  snoRNA => 0.07,
  snRNA => 0.08,
  rRNA => 0.09,
  tRNA => 0.10,
  exon => 0.11,
  three_prime_UTR => 0.12,
  CDS => 0.13,
  five_prime_UTR => 0.14,
  protein_match => 0.15,
  cDNA_match => 0.16,
  match => 0.17,
  match_part => 0.18,
  protein_match => 0.19,
  expressed_sequence_match => 0.20,
  translated_nucleotide_match => 0.21,
  expressed_sequence_match => 0.22,
  contig => 0.23,
);
my $line;
while ($line = <>) {
   last unless $line =~ /^#/;
   print $line;
}
my @lines = ($line, <>);
@lines = map {my @a = split /\t/; \@a;} @lines;
@lines = sort {
    $a->[0] cmp $b->[0] || $a->[3] <=> $b->[3] || $type_collate{$a->[2]} cmp $type_collate{$b->[2]}
} @lines;
foreach my $line (@lines) {
    print join("\t", @$line);
}

# Changes
#2023-07-08 Add some types for collate sort. With extra items, switch to floats (0.00..0.23) for the collate sort.

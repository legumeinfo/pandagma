#!/usr/bin/env perl
use strict;
use warnings;
use feature "say";
use Getopt::Long;

my ($help, $verbose);
GetOptions ("verbose"  => \$verbose,
            "help"     => \$help
            );

my $usage = <<EOS;
Sort gff3 by feature type and position, such that child features occur below their parents. 
The sorted output should be suitable for indexing using Samtools tabix.

Usage: histplot.pl < file.gff3 

    -verbose   Some debug info, printed as comments to top of file (not for production)
    -help      This message
EOS

die "$usage\n" if $help;

my %type_collate = (
  region => 10,
  scaffold => 20,
  contig => 30,
  gene => 40,
  pseudogene => 50,
  transposable_element_gene => 60,
  mRNA => 70,
  transcript => 80,
  ncRNA => 90,
  lnc_RNA => 100,
  snoRNA => 110,
  snRNA => 120,
  rRNA => 130,
  tRNA => 140,
  exon => 150,
  three_prime_UTR => 160,
  CDS => 170,
  five_prime_UTR => 180,
  protein_match => 190,
  cDNA_match => 200,
  match => 210,
  match_part => 220,
  expressed_sequence_match => 230,
  translated_nucleotide_match => 240,
);

my $last_collate_number = 10*(keys %type_collate);

# Print comment lines. If a feature type in column 3 is not in %type_collate, 
# then add the new item at the end (with a higher value, to be sorted lower).
my $line;
my @lines;
while ($line = <>) {
  chomp $line;
  my @parts = split("\t", $line);
  next if ($line =~ /^\s*$/);
  if ($line =~ /^#/){
    say $line;
    next;
  }
  else {
    push @lines, $line;
    my $type = $parts[2];
    unless ($type_collate{$type}){ 
      $last_collate_number += 10;
      $type_collate{$type} = $last_collate_number;
    }
  }
}

# Debugging, with -verbose : print all types (column 3) and the count of that type
if ($verbose){
  foreach my $key (sort { $type_collate{$a} <=> $type_collate{$b} } keys %type_collate){
    say "# Collate position of type $key: $type_collate{$key}";
  }
}

@lines = map {my @a = split /\t/; \@a;} @lines;
@lines = sort {
  $a->[0] cmp $b->[0] || $a->[3] <=> $b->[3] || $type_collate{$a->[2]} cmp $type_collate{$b->[2]}
} @lines;
foreach my $line (@lines) {
  say join("\t", @$line);
}

__END__
Versions
2023-07-08 Add some types for collate sort. With extra items
2023-12-18 Handle arbitrary types (col 3), adding nonstandard ones to the end of the type_collate hash.


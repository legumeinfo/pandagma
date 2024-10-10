#!/usr/bin/env perl

use strict; 
use warnings;
use feature "say";
use Getopt::Long;

my ($outfile, $verbose, $help);
my $format = "list";

GetOptions ( 
  "format:s"  => \$format,
  "outfile:s" => \$outfile,
  "verbose"   => \$verbose,
  "help"      => \$help,
);

my $usage = <<EOS;
From a gene annotation file in gff format, report the ID for the longest mRNA for each gene, 
and the parent (gene) ID.

Usage: cat FILE.gff | longest_variant_from_gff.pl
  -format   Output format. Options: list (default) or bed. The bed format is an extension of bed6,
             with mRNA length in the fifth-column score, and parent ID in the seventh column:
               seqID, start, end, mRNAID, length, orient, geneID
  -outfile  Filename for output; otherwise output is to STDOUT.
  -verbose  Some intermediate output
  -help     This message
EOS

die "$usage\n" if $help || (-t STDIN);
die "Option -format must be either \"list\" (default) or \"bed\"" if $format !~ /list|bed/;

my $OUT;
if ($outfile){
  open($OUT, ">", $outfile) or die "Can't open out $outfile: $!\n";
}

my @gene_mrna_coords_AoA;
my @F;
while ( my $line = <> ) {
  chomp $line;
  next if ( $line =~ /^#/ );
  my @F = split(/\t+/, $line);
  my ($seqID, $type, $start, $end, $orient, $ninth) = ($F[0], $F[2], $F[3], $F[4], $F[6], $F[8]);
  my @attrs = split(/;/, $ninth);
  if ($type =~ /mRNA/){
    my $mrnaID = $ninth;
    my $parent = $ninth;
    $mrnaID =~ m/ID=([^;]+)/; $mrnaID = $1;
    $parent =~ m/.*arent=([^;]+)/; $parent = $1;
    #push @gene_mrna_coords_AoA, [$seqID, $parent, $mrnaID, $start, $end, $end-$start];
    push @gene_mrna_coords_AoA, [$seqID, $start, $end, $mrnaID, $end-$start, $orient, $parent];
  }
}

# Sort the array of arrays by seqID, then by parent, then by mRNA length (reverse numerically)
my @sorted = sort {
  $a->[0] cmp $b->[0]  # seqID first
           ||
  $a->[6] cmp $b->[6]  # then gene ID
           ||
  $b->[4] <=> $a->[4]  # mRNA length (longest to shortest)
           ||
  $a->[1] <=> $b->[1]  # mRNA start (as tiebreaker if needed)
} @gene_mrna_coords_AoA;

my %seen_gene;
for my $aref ( @sorted ) {
  my ($seqID, $start, $end, $mrnaID, $length, $orient, $parent) = @$aref;
  if ($seen_gene{$parent}){ # We've seen this gene, so it's not first (and longest)
    if ($verbose){
      if ($format =~ /bed/){
        &printstr( join("\t", "#SKIP", $start, $end, $mrnaID, $length, $orient, $parent ) );
      }
      else {
        &printstr( "#SKIP $mrnaID" );
      }
    }
    next;
  }
  else {
    $seen_gene{$parent}++;
    if ($format =~ /bed/){
      &printstr( join("\t", $seqID, $start, $end, $mrnaID, $length, $orient, $parent ) );
    }
    else {
      &printstr( $mrnaID );
    }
  }
}

#####################
sub printstr {
  my $str_to_print = join("", @_);
  if ($outfile) {
    say $OUT $str_to_print;
  }
  else {
    say $str_to_print;
  }
}

__END__
VERSIONS
S. Cannon
2023-03-15 Initial version. Tested on a GenBank RefSeq GFF.
2023-03-06 Add options for printing to list or bed6+ format; and to file or stdout


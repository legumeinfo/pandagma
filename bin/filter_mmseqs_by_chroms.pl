#!/usr/bin/env perl

# PROGRAM: filter_mmseqs_by_chroms.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2021
# see description under Usage

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($chr_pat, $verbose, $help);
my $keepself = "";

GetOptions (
  "chr_pat=s" =>  \$chr_pat,   # required
  "keepself"  =>  \$keepself,
  "verbose"  =>   \$verbose,
  "help" =>       \$help
);

my $scriptname = basename($0);

my $usage = <<EOS;
  Usage: cat HOMOLOGY_FILE[S] | $scriptname -chr_pat FILE [-options] 
  
  Given homology data (on STDIN) with the form:
      chromosome__gene__start__end   chromosome__gene__start__end
    or
      chromosome__gene__start__end__+   chromosome__gene__start__end__-

  ... filter on patterns in the two chromosomes, provided in an input file
  via the flag -chr_pat .
  
  Example1 from chr_pat file ...
    11 11
    13 11
    20 20
  ... to match these gene ID pairs (genes, starts, and ends have placeholders here):
  glyma.Wm82.gnm2.Gm11__gene__start__end  glyma.Wm82.gnm2.Gm11__gene__start__end
  glyma.Lee.gnm1.Gm20__gene__start__end  glyso.PI483463.gnm1.Gs20__gene__start__end
  glyso.W05.gnm1.Chr11__gene__start__end  glyma.Wm82.gnm2.Gm13__gene__start__end
  glyma.Wm82.gnm2.Gm13__gene__start__end  glyso.W05.gnm1.Chr11__gene__start__end

  Example2 from chr_pat file ...
    1 1
    9 10
    10 10
  ... to match these gene ID pairs (genes, starts, and ends have placeholders here):
  Zm-CML103_NAM-1.chr1__gene_start_end  Zm-Oh7B_NAM-1.chr1__gene_start_end
  Zm-CML103_NAM-1.chr10__gene_start_end  Zm-Oh7B_NAM-1.chr10__gene_start_end
  Zm-Oh7B_NAM-1.chr9__gene_start_end  Zm-CML103_NAM-1.chr10__gene_start_end  

  Flags and parameters:
   -chr_pat  -- file with regex for filtering chromosomes in two fields, as indicated above **
   -keepself -- (boolean) Include chromosome self-matches (from the same genome) [false]
   -verbose  -- (boolean) for some intermediate output to STDERR
   -help     -- (boolean) for this info
   
   ** = required
EOS

die "\n$usage\n" 
  if ($help or !defined($chr_pat) );

# read chr patterns in
open( my $PAT_IN, '<', $chr_pat ) or die "can't open chr_pat $chr_pat: $!";
my %rexen;
while (<$PAT_IN>){ # two fields, e.g. "01 01"
  chomp;
  next if (/^#/ || /^$/);
  my $line = $_;
  my ($left, $right) = split(/\s+/);
  $left =~ s/^0+(\d+)/$1/;
  $right =~ s/^0+(\d+)/$1/;
  
  my $chr_rex;
  # Store the pattern for later regex use
  $chr_rex = qr(\S+\D0*$left \S+\D0*$right$);
  if ($verbose) {print "$left $right; $chr_rex\n"}
  $rexen{$chr_rex}++;

  # Store the reverse pattern, if the left and right chroms are not the same
  if ($left ne $right) {
    $chr_rex = qr(\S+\D0*$right \S+\D0*$left$);
    if ($verbose) {print "$right $left; $chr_rex\n"}
    $rexen{$chr_rex}++;
  }
}

# Process homology data, comparing to the stored chr patterns
while (my $line = <>) {
  $line =~ s/>//g; # data shouldn't have ">", but do this to make sure.
  my ($gene1, $gene2) = split(/\t/, $line);

  $gene1 =~ s/__[+-]$//; # strip orientation, since that is not used by subsequent DAGChainer
  $gene2 =~ s/__[+-]$//; # strip orientation, since that is not used by subsequent DAGChainer

  my @parts1 = split(/__/, $gene1);
  my @parts2 = split(/__/, $gene2);
  my $chr1 = $parts1[0];
  my $geneID1 = $parts1[1];
  my $chr2 = $parts2[0];
  my $geneID2 = $parts2[1];
  unless ($keepself){ next if ($chr1 eq $chr2 && $geneID1 eq $geneID2) };
  my $matches = 0;
  foreach my $chr_rex (keys %rexen){
    if ($verbose){print "TEST: $chr1 $chr2 =~ /$chr_rex/\n"}
    if ("$chr1 $chr2" =~ /$chr_rex/){
      print join("\t", @parts1), "\t", join("\t", @parts2);
      $matches++;
    }
  }
  if ($matches == 0) {
    next;
    #print "NO MATCH:"$parts1[1]\t$parts2[1]\n";
  }
  if ($verbose){print "\n"}
}

__END__
VERSIONS
2021-07-11 S. Cannon. 
2021-10-15 strip ">" from data
2021-10-19 For matches between different chromosomes, handle both directions, e.g. 11 13 and 13 11
           and terminate regexes with "$" to ensure match of chromosome number, not preceding text
2021-10-19 Add "keepself" flag
2022-12-31 In regex, ignore leading zeroes in e.g. chr01
2023-03-10 Change -noself to -keepself, and allow orientation code in chromosome__gene fields.

#!/usr/bin/env perl

# PROGRAM: filter_mmseqs_by_chroms.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2021
# see description under Usage

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use feature "say";

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
    optionally with a preceding cluster-id column:
      pan00010  chromosome__gene__start__end   chromosome__gene__start__end

    or an equivalent 12-column -m8 format, with the first two fields as above:
      chromosome__gene__start__end__+   chromosome__gene__start__end__- [blast -m8 format]

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
  if ($verbose) {say "$left $right; $chr_rex"}
  $rexen{$chr_rex}++;

  # Store the reverse pattern, if the left and right chroms are not the same
  if ($left ne $right) {
    $chr_rex = qr(\S+\D0*$right \S+\D0*$left$);
    if ($verbose) {say "$right $left; $chr_rex"}
    $rexen{$chr_rex}++;
  }
}

# Process homology data, comparing to the stored chr patterns
my @rest;
while (my $line = <>) {
  chomp $line;
  $line =~ s/>//g; # data shouldn't have ">", but do this to make sure.
  my @fields = split(/\t/, $line);
  my ($panID, $gene1, $gene2);
  my $format; # to note either two- or three-column form; values 2 or 3

  if (scalar(@fields) == 2 && $fields[0] =~ /__/ && $fields[1] =~ /__/){
    ($gene1, $gene2) = @fields;
    $format = 2;
  }
  elsif (scalar(@fields) == 3 && $fields[1] =~ /__/ && $fields[2] =~ /__/){
    ($panID, $gene1, $gene2) = @fields;
    $format = 3;
  }
  elsif (scalar(@fields) == 12 && $fields[0] =~ /__/ && $fields[1] =~ /__/){
    ($gene1, $gene2, @rest) = @fields;
    $format = 12;
  }
  else {
    warn "WARN| Unexpected format: $line\n";
    next;
  }

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
    if ($verbose){say "TEST: $chr1 $chr2 =~ /$chr_rex/"}
    if ("$chr1 $chr2" =~ /$chr_rex/){
      if ($format == 2){
        say join("\t", @parts1), "\t", join("\t", @parts2);
      }
      elsif ($format == 3) {
        say $panID, "\t", join("\t", @parts1), "\t", join("\t", @parts2);
      }
      elsif ($format == 12) {
        say join("__", @parts1), "\t", join("__", @parts2), "\t", join("\t", @rest)
      }
      $matches++;
    }
  }
  if ($matches == 0) {
    next;
    #say "NO MATCH:"$parts1[1]\t$parts2[1]";
  }
  if ($verbose){ say "" }
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
2023-08-27 Handle a three-column form, to allow match data with leading pangene ID field
2023-12-08 Add option to handle 12-column BLAST -m8 format

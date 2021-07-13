#!/usr/bin/env perl

# PROGRAM: filter_mmseqs_by_chroms.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2021
# see description under Usage

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($chr_pat, $help);

GetOptions (
  "chr_pat=s" =>  \$chr_pat,   # required
  "help" =>       \$help
);

my $scriptname = basename($0);

my $usage = <<EOS;
  Usage: cat HOMOLOGY_FILE[S] | $scriptname -chr_pat FILE [-options] 
  
  Given homology data (on STDIN) with the form:
    chromosome__gene__start__end   chromosome__gene__start__end
  filter on patterns in the two chromosomes, provided in an input file
  via the flag -chr_pat .
  
  Example from chr_pat ...
    11 11
    13 11
    11 13
    20 20
  ... to match these gene ID pairs (genes, starts, and ends have placeholders here):
  glyma.Wm82.gnm2.Gm11__gene__start__end  glyma.Wm82.gnm2.Gm11__gene__start__end
  glyma.Lee.gnm1.Gm20__gene__start__end  glyso.PI483463.gnm1.Gs20__gene__start__end
  glyso.W05.gnm1.Chr11__gene__start__end  glyma.Wm82.gnm2.Gm13__gene__start__end
  glyma.Wm82.gnm2.Gm13__gene__start__end  glyso.W05.gnm1.Chr11__gene__start__end

  Flags and parameters:
   -chr_pat -- regex for filtering chromosomes in two fields, as indicated above **
   -help       for more info
   
   ** = required
EOS

die "\n$usage\n" 
  if ($help or !defined($chr_pat) );

# read chr patterns in
open( my $PAT_IN, '<', $chr_pat ) or die "can't open chr_pat $chr_pat: $!";
my %rexen;
while (<$PAT_IN>){ # two fields, e.g. "01 01"
  chomp;
  next if (/^#/);
  my $line = $_;
  $line =~ s/\s+/ /; # separate the two fields by a single space, if not already
  $line =~ s/(\S+)\s+(\S+)/\\S+$1 \\S+$2/;
  my $chr_rex = qr($line);
  $rexen{$chr_rex}++;
}

# Process homology data, comparing to the stored chr patterns
while (my $line = <>) {
  my ($gene1, $gene2) = split(/\t/, $line);
  my @parts1 = split(/__/, $gene1);
  my @parts2 = split(/__/, $gene2);
  my $matches = 0;
  foreach my $chr_rex (keys %rexen){
    if ("$parts1[0] $parts2[0]" =~ /$chr_rex/){
      #print $line;
      print join("\t", @parts1), "\t", join("\t", @parts2);
      $matches++;
    }
  }
  if ($matches == 0) {
    next;
    #print "NO MATCH:"$parts1[1]\t$parts2[1]\n";
  }
}

__END__
VERSIONS
v0.01 2021-07-11 S. Cannon. 


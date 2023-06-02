#!/usr/bin/env perl

# PROGRAM: filter_mmseqs_by_quotas.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2023
# see description under Usage

use strict;
use warnings;
use feature "say";
use Getopt::Long;
use File::Basename;
use Data::Dumper;

my ($quotas, $verbose, $help);
my $species_regex = '^([^.]+)\..+';

GetOptions (
  "quotas=s" =>        \$quotas,    # required
  "species_regex:s" => \$species_regex,   
  "verbose" =>         \$verbose,
  "help" =>            \$help
);

my $scriptname = basename($0);

my $usage = <<EOS;
  Usage: cat HOMOLOGY_FILE[S] | $scriptname -quotas FILE [-options] 

  Given homology data (on STDIN) with the form:
      qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    where qseqid and sseqid have one of the following two forms:
      chromosome__gene__start__end   chromosome__gene__start__end 
    or
      chromosome__gene__start__end__+   chromosome__gene__start__end__-

  ... filter to the number of sseqid matches per qseqid indicated in the quotas file,
  such that up to the quota number of subject chromosomes are seen per query.
  Output is to stdout.

  Flags and parameters:
    -quotas   -- file with regex for filtering chromosomes in two fields, as indicated above **
    -species_regex  -- quoted regular expression for recognizing and capturing the species ID from 
                       the query and subject fields. Default: '^([^.]+)\..+'
    -verbose  -- (boolean) for some intermediate output to STDERR
    -help     -- (boolean) for this info
  
    ** = required

  The quotas file should have two fields:
    - field 1 matching the first part of the subject seqID (sseqid), sufficient to identify the 
              species or unique taxon in gene families -- for example, "Phaseolus" or "Zea".
    - field 2 being the number of expected paralogs in a gene family of sseqid sequences, 
              relative to the taxonomic depth of interest. For example, at the level of the legume family, 
              two paralogs of Phaseolus are expected relative to the family root, but four of Glycine.

  Example of a quotas file for four species, with self-matches included, corresponding 
  with the evolutionary depth and WGD histories for the legume family (originating ~70 mya):
    Arachis	4
    Cercis	1
    Glycine	4
    Phaseolus	2
EOS

die "\n$usage\n" if ( $help or !defined($quotas) );

my $SP_REX = qr($species_regex);

# read quotas file
open( my $QUOTA_FH, '<', $quotas ) or die "can't open quotas $quotas: $!";
my %quotas;
while (<$QUOTA_FH>){ # two fields, e.g. Glycine 4
  chomp;
  next if (/^#/ || /^$/);
  my $line = $_;
  my ($sbj, $quota) = split(/\s+/);
  #say "AA: ($sbj, $quota)";
  $quotas{$sbj} = $quota;
}

# Process homology data, with default sort by query & E-value
my %ct_Qgn_Schr;
my %seen_Qgn_Schr;
my @homology_by_qry;
my $ct1 = 0;
while (my $line = <>) {
  chomp $line;

  # The query and subject should have 5 sub-fields, separated by "__":
  # species_chr  geneID  start  end  orientation 
  # The remainder, in @rest, should have 10 elements (of 12 fields in the BLAST -m8 format input)
  my ($Q_id, $S_id, @rest) = split(/\t/, $line);
  my @Qparts = split(/__/, $Q_id);
  my @Sparts = split(/__/, $S_id);
  
  # Push all elements back into an array ... swapping query and subject (Q S ==> S Q)
  # and containing the split sub-arrays:
  # 5 in @Sparts, 5 in @Qparts, 10 in @rest
  #    0     1    2       3     4       5     6    7       8     9
  #   S_chr S_id S_start S_end S_or    Q_chr Q_id Q_start Q_end Q_or   
  #    10    11  12   13  14     15   16     17    18     19
  #   piden len mism gap qstart qend sstart ssend evalue bitscore
  push(@homology_by_qry, [@Sparts, @Qparts, @rest]);

  my $Q_sp = $Qparts[0];
  $Q_sp =~ s/$SP_REX/$1/;
  my $S_sp = $Sparts[0];
  $S_sp =~ s/$SP_REX/$1/;

  # say "Q_sp: $Q_sp; S_sp: $S_sp"; 

  my ($Qchr, $Qgn, $Qsta, $Qend, $Qor) = @Qparts;
  my ($Schr, $Sgn, $Ssta, $Send, $Sor) = @Sparts;
  my $Qgn_Schr = "$Qgn.x.$Schr";
  
  #say "Qgn: $Qgn   Schr: $Schr   Qgn_Schr: $Qgn_Schr";
  $ct_Qgn_Schr{$Qgn}++ unless ($seen_Qgn_Schr{$Qgn_Schr});
  $seen_Qgn_Schr{$Qgn_Schr}++;
  if ( $ct_Qgn_Schr{$Qgn} <= $quotas{"$S_sp"} ){ # direction: Q ==> S b/c sorting is by Q+score and we're looking to S
    $ct1++;
    if ($verbose){
      say "BB: $Qchr\t[$Qgn]\t[$Schr]\t$Sgn\t$ct_Qgn_Schr{$Qgn}\t$ct1\t$quotas{$S_sp}";  # for debugging 
    }
    # Print Qchr Qgn Qsta Qend Sor   Schr Sgn Ssta Send Sor
    unless ($verbose){
      say join("\t", @Qparts, @Sparts);
    }
  }
}

#say "\n", Dumper(@homology_by_qry), "\n";

# Sort the homology data by subject_identity (S_id; field 2), then by bitscore
my @homology_by_sbj =
   sort {
     $a->[1] cmp $b->[1] # sort 1st (alpha) by 2nd field (S_id)
            ||
     $b->[19] <=> $a->[19] # reverse sort 2nd by 20th field (bitscore)
   } @homology_by_qry;

#say "\n", Dumper(@homology_by_sbj), "\n";

# Process homology data, sorted by subject
undef %ct_Qgn_Schr;
undef %seen_Qgn_Schr;
my $ct2 = 0;
foreach my $row (@homology_by_sbj) {
  my @parts = @{$row};
  #say "CC: @parts[0..10]";

  # Note: S and Q are swapped relative to input data, since they are swapped in @homology_by_sbj
  my $Q_sp = $parts[0];
  $Q_sp =~ s/$SP_REX/$1/;
  my $S_sp = $parts[5];
  $S_sp =~ s/$SP_REX/$1/;

  # say "Q_sp: $Q_sp; S_sp: $S_sp"; 

  # Now revert to Q for first group of fields and S for second group of fields
  my ($Qchr, $Qgn) = ($parts[0], $parts[1]);
  my ($Schr, $Sgn) = ($parts[5], $parts[6]);
  my $Qgn_Schr = "$Qgn.$Schr";
  
  #say "Qgn: $Qgn   Schr: $Schr   Qgn_Schr: $Qgn_Schr";
  $ct_Qgn_Schr{$Qgn}++ unless ($seen_Qgn_Schr{$Qgn_Schr});
  $seen_Qgn_Schr{$Qgn_Schr}++;
  if ( $ct_Qgn_Schr{$Qgn} <= $quotas{$S_sp}){ 
    $ct2++;
    if ($verbose){
      say "DD: [$Qchr]\t$Qgn\t$Schr\t[$Sgn]\t$ct_Qgn_Schr{$Qgn}\t$ct2\t$quotas{$S_sp}";  # for debugging
    }
    # Print Qchr Qgn Qsta Qor   Qend Schr Sgn Ssta Send Sor 
    unless ($verbose){
      say join("\t", @parts[0..4], @parts[5..9]);
    }
  }
}

__END__
VERSIONS
2023-05-16 S. Cannon. Initial version.
2023-05-19 Add back orientation field to the output, and swap Q and S in second-pass processing
2023-06-01 Change quota input from Q & S fields to just one field - giving the expected duplication count wrt the family root.

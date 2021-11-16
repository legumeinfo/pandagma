#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum min max);
use Bio::SeqIO;

my $usage = <<EOS;
Given a fasta file of sequences in pan-genes sets, identify a representative sequence to 
for each pan-gene, with representative chosen from the median(ish) sequence length.
In cases where the median would be the average of two central values, 
the larger of the two is selected - e.g. of (100,150), 150 is returned. 

The deflines should have the form
  >pan00001__gene_ID
The defline should be structured with as: pan-geneID and gene ID, separated by "__".

Usage: pick_family_rep.pl -infile family_fasta.fa  [options]
   Input has two fields, separated by "__"
     panID__geneID

  OPTIONS:
    -infile   Input fasta file, with pangene IDs + gene IDs
    -outfile  Specify OUT_FH; otherwise, default to STDOUT.
    -verbose  For some debugging info, to STDOUT
    -help     This message. 
EOS

my ($infile, $outfile, $verbose, $help);

GetOptions (
  "infile=s" =>  \$infile,
  "outfile:s" => \$outfile,
  "verbose" =>   \$verbose,
  "help" =>      \$help,
);

die "\n$usage\n" if ( $help or ! $infile );

my $logstr;

my $OUT_FH;
if ($outfile) { open ($OUT_FH, ">", $outfile) or die "\nUnable to open output file for writing: $!\n\n"; }

# Put file into an array, and get count of each scaffold
my (%HoA_PI_lengths, %geneID_seq, %geneID_len);

my $seqin = Bio::SeqIO -> new(-file => "$infile", -format => "Fasta");
while (my $seqobj = $seqin->next_seq() ) {
  my $display_id = $seqobj->display_id();
  if ($display_id !~ /__/){
    die "Sequence ID should have two fields, separated by '__'. Problematic ID: $display_id\n"
  }
  my ($panID, $gene_ID) = split(/__/, $display_id);
  my $len = $seqobj->length;
  my $seq = $seqobj->seq;

  push @{$HoA_PI_lengths{$panID}}, $len;
  $geneID_seq{$display_id} = $seq;
  $geneID_len{$display_id} = $len;
}

my (%count_PI, %min_PI, %median_PI, %max_PI);
for my $panID ( keys %HoA_PI_lengths ) {

  my $count = scalar(@{ $HoA_PI_lengths{$panID} });
  $count_PI{$panID} = $count;

  my $min = min(@{ $HoA_PI_lengths{$panID} });
  $min_PI{$panID} = $min;

  my $median = sprintf("%i", calc_median(\@{ $HoA_PI_lengths{$panID} }));
  $median_PI{$panID} = $median;

  my $max = max(@{ $HoA_PI_lengths{$panID} });
  $max_PI{$panID} = $max; 
}

# Traverse sequences again, reporting the first sequence per panID with median length
if ($verbose){print "panID__geneID\tcount\tmin\tmedian\tmax\n" }
$seqin = Bio::SeqIO -> new(-file => "$infile", -format => "Fasta");
my %seen_PI;
while (my $seqobj = $seqin->next_seq() ) {
  my $display_id = $seqobj->display_id();
  my ($panID, $gene_ID) = split(/__/, $display_id);
  my $len = $seqobj->length;
  my $seq = $seqobj->seq;

  if($len == $median_PI{$panID}){
    if ($seen_PI{$panID}){ next }
    if ($verbose) { 
      print "$display_id\t$count_PI{$panID}\t$min_PI{$panID}\t$median_PI{$panID}\t$max_PI{$panID}\n";
    }
    if ($outfile){
      print $OUT_FH ">$display_id\n$seq\n";
    }
    else {
      print ">$display_id\n$seq\n";
    }
    $seen_PI{$panID}++;
  }
}

##################################################
# SUBRUTINES
#
# Return median value for an array of numbers
# See http://stackoverflow.com/questions/5119034/using-perl-to-find-median-mode-standard-deviation
sub calc_median {
  my $value_ref = shift;
  my @values = @$value_ref;
  my $median;
  my $mid = int ((scalar @values)/2);
  my @sorted_values = sort {$a <=> $b} @values;
  #if (@values % 2) {
  $median = $sorted_values[ $mid ];
  #} else {
  #  $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2; ## Not doing this; want larger of the two
  #} 
  return $median;
}

__END__
2021

v01 11-15 Initial version.


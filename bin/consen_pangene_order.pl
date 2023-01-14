#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use File::Basename;
my $PROGRAM_NAME = basename($0);

my $usage = <<EOS;
Given five-column input with pan-IDs, chromosomes, and order along each chromosome,
and genomic start and end coordinates, report consensus chromosome for each pan-gene, 
and average position for genes with the consensus chromosome for that pan-gene.

Usage: cat PANGENE_TABLE | consen_pangene_order.pl [options]
   Input has five columns:
     panID  chr  order-per-chr, start, end

  OPTIONS:
    -outfile  Specify OUT_FH; otherwise, default to STDOUT.
    -prefix   Characters used to prefix the consensus chromosome names; default "chr"
                Examples: Glycine.pan3.chr, Medicago.pan1.chr
    -verbose  For some debugging info, to STDOUT
    -help     This message. 
EOS

my ($outfile, $verbose, $help);
my $prefix="chr";

GetOptions (
  "outfile:s" => \$outfile,
  "prefix:s" =>  \$prefix,
  "verbose" =>   \$verbose,
  "help" =>      \$help,
);

#die "\n$usage\n" if ( $help or not -t STDIN );
die "\n$usage\n" if ( $help );

my $logstr;

my $OUT_FH;
if ($outfile) { open ($OUT_FH, ">", $outfile) or die "\nUnable to open output file for writing: $!\n\n"; }

# Put file into an array, and get count of each scaffold
my %count_panID_chr;
my %HoH_panID_chr;
my @pangene_table;
while (<>) {
  chomp;
  next unless (/^\S+/);
  my $line = $_;
  my ($panID, $chr, $order, $start, $end) = split(/\t/, $line);
  next if ($chr =~ /chloro|chl|CP|mito|MT|Pt$|ctg|contig|tig|pilon|scaff|sc\d+|un\w+\d+/i);
  $HoH_panID_chr{$panID}{$chr}++;
  push(@pangene_table, $line);
}

my %top_chr;
my %chr_ct_top_chr;
if ($verbose) {print "#pangeneID\tchr:count ...\n"}
foreach my $panID (sort keys %HoH_panID_chr) {
  if ($verbose){ $logstr .= "$panID\t" }
  # Sort chromosomes by count of chromosomes seen for this panID
  foreach my $chr ( sort { $HoH_panID_chr{$panID}{$b} <=> $HoH_panID_chr{$panID}{$a} } 
                            keys %{$HoH_panID_chr{$panID}} ) {
    my $chr_ct = $HoH_panID_chr{$panID}{$chr};
    # Store first chromosome (most frequent), and the number of times it was seen.
    unless ($top_chr{$panID}){ 
      $top_chr{$panID} = $chr; 
      $chr_ct_top_chr{$panID} = $chr_ct;
    }
    if ($verbose){ $logstr .= "$chr:$chr_ct " }
  }
  if ($verbose){ $logstr .= "\n" }
}
if ($verbose){ print "$logstr\n" }

if ($verbose) {print "#pangeneID\ttop_chr:count\n"}
foreach my $panID (sort keys %top_chr) {
  if ($verbose){ printf "%s\t%s:%s\n", $panID, $top_chr{$panID}, $chr_ct_top_chr{$panID} }
}
if ($verbose) {print "\n"}

#if ($OUT_FH){
#  print $OUT_FH "#pangeneID\ttop_chr\tmedian_order\n";
#} else {
#  print "#pangeneID\ttop_chr\tmedian_order\n";
#}

my %orders_per_top_chr_HoA;
my %starts_per_top_chr_HoA;
my %ends_per_top_chr_HoA;
foreach (@pangene_table) {
  my ($panID, $chr, $order, $start, $end) = split(/\t/);
  next if ($chr =~ /chloro|chl|CP|mito|MT|Pt$|ctg|contig|tig|pilon|scaff|sc\d+|un\w+\d+/i);
  $chr =~ s/\S+\.\D+(\d+)/$1/;
  if ($chr eq $top_chr{$panID}){
    push @{ $orders_per_top_chr_HoA{$panID} }, $order;
    push @{ $starts_per_top_chr_HoA{$panID} }, $start;
    push @{ $ends_per_top_chr_HoA{$panID} }, $end;
  }
}

for my $panID ( keys %starts_per_top_chr_HoA ) {
  my $median_order = sprintf("%i", calc_median(\@{ $orders_per_top_chr_HoA{$panID} }));
  my $median_start = sprintf("%i", calc_median(\@{ $starts_per_top_chr_HoA{$panID} }));
  my $median_end = sprintf("%i", calc_median(\@{ $ends_per_top_chr_HoA{$panID} }));

  my $average_order = sprintf("%i", calc_average(\@{ $starts_per_top_chr_HoA{$panID} }));

  if ($OUT_FH){
    #print $OUT_FH "$panID\t$prefix$top_chr{$panID}\t$median_order\t$average_order\n"
    print $OUT_FH "$panID\t$prefix$top_chr{$panID}\t$median_order\t$median_start\t$median_end\n"
  } else {
    #print "$panID\t$prefix$top_chr{$panID}\t$median_order\t$average_order\n";
    print "$panID\t$prefix$top_chr{$panID}\t$median_order\t$median_start\t$median_end\n";
  }
  #my $average_order = sprintf("%i", calc_average(\@{ $starts_per_top_chr_HoA{$panID} }));
  #print "$panID: MED: $median_order\n";
  #print "$panID: AVE: $average_order\n";
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
  if (@values % 2) {
    $median = $sorted_values[ $mid ];
  } else {
    $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
  } 
  return $median;
}

sub calc_average {
  my $value_ref = shift;
  my @values = @$value_ref;
  my $average = sum(@values)/scalar(@values);
  return $average;
}

__END__
2022
10-29 Initial version.
11-04 Read from STDIN rather than file
2023
01-13 Rename from consen_pangene_order.pl to consen_pangene_order.pl. This version operates
       on three-column ordinal data, rather than five-column with genomic position.


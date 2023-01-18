#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use File::Basename;
use Data::Dumper;
use feature "say";
my $PROGRAM_NAME = basename($0);

my $usage = <<EOS;
Given five-column input (pan-IDs, chromosomes, annot.chr#, start, end), 
calculate gene order per annotation, then report
  (1) pan-ID;
  (2) consensus (median) chromosome for each pan-gene; 
  (3) consensus ordering on the chromosome, based on median of orders of the genes in each pan-gene;
  (4) median start position from the pan-genes;
  (5) median end position from the pan-genes.

Usage: cat PANGENE_TABLE | consen_pangene_order.pl [options]
   Input has five columns:
     panID  geneID  chr  annot.chr  start  end

  OPTIONS:
    -outfile  Specify OUT_FH; otherwise, default to STDOUT.
    -prefix   Characters used to prefix the consensus chromosome names; default "chr"
                Examples: Glycine.pan3.chr, Medicago.pan1.chr
    -make_new_id  Make new gene IDs based on chromosome + order, e.g. Prefix01_000100
    -verbose  For some debugging info, to STDOUT
    -help     This message. 
EOS

my ($outfile, $make_new, $verbose, $help);
my $prefix="chr";

GetOptions (
  "outfile:s" => \$outfile,
  "prefix:s" =>  \$prefix,
  "make_new" =>  \$make_new,
  "verbose" =>   \$verbose,
  "help" =>      \$help,
);

#die "\n$usage\n" if ( $help or not -t STDIN );
die "\n$usage\n" if ( $help );

my $logstr;

my $OUT_FH;
if ($outfile) { open ($OUT_FH, ">", $outfile) or die "\nUnable to open output file for writing: $!\n\n"; }

# Put file into an array, and get count of each molecule (chromosome or scaffold.)
my %count_panID_chr;
my %HoH_panID_chr;
my @pangene_table;
while (<>) {
  chomp;
  next unless (/^\S+/);
  my $line = $_;
  my ($panID, $gene, $ann_chr, $start, $end) = split(/\t/, $line);
  # Split third field, a annot.chr string, into annot and chr
  $ann_chr =~ /(\S+)\.(\D+\w+\D+)(\d+)$/;
  my ($ann, $chr_pre, $chr) = ($1, $2, $3);
  $chr_pre =~ s/[_.]$//;
  # Next: skip genes on scaffolds and other non-chromosome molecules
  if ($chr_pre =~ /chloro|chl|CP|mito|MT|ctg|contig|tig|pilon|scaff|sc|super|un\w+\d+/i){
    if ($verbose){ say "Skipping [$chr_pre $chr]\t$gene" }
  }
  else {
    $chr =~ s/^0*([^0]+)/$1/;
    #say "A:\t($panID, $gene, $ann, $chr_pre, $chr)";
    $HoH_panID_chr{$panID}{$chr}++;
    my @six_elts = ( $panID, $gene, $ann, $chr, $start, $end );
    push( @pangene_table, \@six_elts );
  }
}

# Sort @pangene_table to determine order per (1) annot; (2) chromosome; (3) posn
my @sorted_table = sort { 
     $a->[2] cmp $b->[2] || $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] 
   } @pangene_table; 

# Calculate gene order per annot-and-chromosome
my ($prAnn, $prChr) = ("", "");
my $ord=0;
my @elts_with_order;
my @pangene_table_ordered;
my %seen_chr;
foreach my $row ( @sorted_table ) {
  my ($panID, $gene, $ann, $chr, $start, $end) = @$row;
  #say "$panID, $gene, $ann, $chr, $start, $end";
  unless ( $seen_chr{$chr} ){ $seen_chr{$chr}++ }
  if ($ann eq $prAnn && $chr == $prChr){
    $ord++;
    @elts_with_order = ( $panID, $ann, $chr, $ord, $start, $end );
    push ( @pangene_table_ordered, [@elts_with_order] );
    #say join("\t", @elts_with_order);
    ($prAnn, $prChr) = ($ann, $chr);
  }
  elsif ($ann ne $prAnn || $chr != $prChr){
    $ord=1;
    @elts_with_order = ( $panID, $ann, $chr, $ord, $start, $end );
    push ( @pangene_table_ordered, [@elts_with_order] );
    #say join("\t", @elts_with_order);
    ($prAnn, $prChr) = ($ann, $chr);
    $ord++;
  }
}
my $num_chrs = keys %seen_chr;
#say "CHROMOSOMES: $num_chrs";

# Find the most frequent chromosome for each pan-gene set
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

#say Dumper(\@pangene_table_ordered);

# Store order, start, and end for each panID
my %orders_per_top_chr_HoA;
my %starts_per_top_chr_HoA;
my %ends_per_top_chr_HoA;
foreach my $row ( @pangene_table_ordered ) {
  my ( $panID, $ann, $chr, $order, $start, $end ) = @$row;
  #say "$panID, $ann, $chr, $order, $start, $end";
  if ($chr eq $top_chr{$panID}){
    push @{ $orders_per_top_chr_HoA{$panID} }, $order;
    push @{ $starts_per_top_chr_HoA{$panID} }, $start;
    push @{ $ends_per_top_chr_HoA{$panID} }, $end;
  }
}

my %consen_pan_coords;
my $most_freq_chr;
for my $panID ( keys %starts_per_top_chr_HoA ) {
  my $median_order = sprintf("%i", calc_median(\@{ $orders_per_top_chr_HoA{$panID} }));
  my $median_start = sprintf("%i", calc_median(\@{ $starts_per_top_chr_HoA{$panID} }));
  my $median_end = sprintf("%i", calc_median(\@{ $ends_per_top_chr_HoA{$panID} }));
  if ( $num_chrs >= 10 ){ # num_chrs has two digits, so pad the chr numbers
    $most_freq_chr = sprintf("$prefix%02d", $top_chr{$panID});
  }
  else { # num_chrs <=9, so don't pad the chr numbers
    $most_freq_chr = sprintf("$prefix%d", $top_chr{$panID});
  }
  my @pan_coords = ( $most_freq_chr, $median_order, $median_start, $median_end );
  #say "($panID, $most_freq_chr, $median_order, $median_start, $median_end)";
  $consen_pan_coords{$panID} = [ @pan_coords ];
}

# Build array of consensus pangene positions, 
# sorted by (1) chromosome, (2) median order, (3) genomic position.
my @sorted_consen_pan_coords;
for my $key (sort { 
    $consen_pan_coords{$a}[0] cmp $consen_pan_coords{$b}[0] or 
    $consen_pan_coords{$a}[1] <=> $consen_pan_coords{$b}[1] or 
    $consen_pan_coords{$a}[2] <=> $consen_pan_coords{$b}[2] 
  } keys %consen_pan_coords ) {
  #say "$key => ", join(", ", @{$consen_pan_coords{$key}});
  push @sorted_consen_pan_coords, [$key, @{$consen_pan_coords{$key}}];
}

#say Dumper(\@sorted_consen_pan_coords);

# Calculate final pangene order per chromosome.
# Use the consensus ordering based on (1) chromosome, (2) median order, (3) genomic position,
# but generate a new monotonic increasing ordering, with ties broken by genomic position.
# Don't print genomic position ($start, $end), which may conflict with consensus order.
$prChr = "";
my $final_ord=1;
foreach my $row ( @sorted_consen_pan_coords ) {
  my ($key, $chr, $ord, $start, $end) = @$row;
  #say join("\t", ($key, $chr, $ord, $start, $end));
  if ($chr eq $prChr){ $final_ord++ }
  elsif ($chr ne $prChr){ $final_ord=1 }
  if ($make_new){ # Construct a new pseudo-ID from the chrom & position
    &printstr( sprintf("%s\t%s_%06d\n", $key, $chr, $final_ord*100 ) );
  }
  else { # print bare numeric order
    &printstr( sprintf("%s\t%s\t%d\n", $key, $chr, $final_ord) );
  }
  $prChr = $chr;
}

##################################################
# SUBRUTINES

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

# Print to outfile or to stdout
sub printstr {
  my $str_to_print = join("", @_);
  if ($outfile) {
    print $OUT_FH $str_to_print;
  }
  else {
    print $str_to_print;
  }
}


__END__
2022
10-29 Initial version.
11-04 Read from STDIN rather than file
2023
01-13 Rename from consen_pangene_order.pl to consen_pangene_order.pl. This version operates
       on the same five-column input data with genomic position, but calculates modal 
       chromosome and gene orders per annotation set; then calculates modal order for the pangene.


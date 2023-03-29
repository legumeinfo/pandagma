#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum min max);
use Bio::SeqIO;
use feature "say";

my $usage = <<EOS;
Given fasta sequences in pan-genes sets (on STDIN), identify a representative sequence to 
for each pan-gene, with representative chosen from the median(ish) sequence length.
In cases where the median would be the average of two central values, 
the larger of the two is selected - e.g. of (100,150), 150 is returned. 

The deflines should have the form
  >pan00001__gene_ID
The defline should be structured with as: pan-geneID and gene ID, separated by "__".

Usage: cat FAMILY_FILE.fa | pick_family_rep.pl [options]
   Input has two fields, separated by "__"
     panID__geneID

  OPTIONS:
    -prefer   Pattern for matching ID(s) to be selected as representative, if that ID is available
              with near-modal length.
    -nostop   Ignore sequences that contain internal stop codons (assumes protein sequence).
    -outfile  Specify OUT_FH; otherwise, default to STDOUT.
    -verbose  For some debugging info, to STDOUT
    -help     This message. 
EOS

my ($prefer, $outfile, $nostop, $verbose, $help);

GetOptions (
  "prefer:s" =>  \$prefer,
  "nostop"    => \$nostop,
  "outfile:s" => \$outfile,
  "verbose" =>   \$verbose,
  "help" =>      \$help,
);

die "\n$usage\n" if ( $help );

my $OUT_FH;
if ($outfile) { open ($OUT_FH, ">", $outfile) or die "\nUnable to open output file for writing: $!\n\n"; }

# Put file into an array, and get count of each scaffold
my (%HoA_PI_lengths, %geneID_seq, %geneID_len);

my $seqin = Bio::SeqIO -> new(-fh => \*STDIN, -format => "Fasta");
while (my $seqobj = $seqin->next_seq() ) {
  my $display_id = $seqobj->display_id();
  if ($display_id !~ /__/){
    die "Sequence ID should have two fields, separated by '__'. Problematic ID: $display_id\n"
  }
  my ($panID, $gene_ID) = split(/__/, $display_id);
  my $len = $seqobj->length;
  my $seq = $seqobj->seq;
  my $type = $seqobj->alphabet;

  # Remove terminal stop character if present
  $seq =~ s/\*$//; 
  $seq =~ s/\.$//; 

  if ($nostop){
    if ($type eq 'dna'){
      die "Option -nostop was selected, but the sequence of $gene_ID appears to be nucleotide\n";
    }
    else {
      if ($seq =~ /\*.+/ || $seq =~ /\..+/){
        warn "Skipping sequence $gene_ID due to an internal stop codon\n";
      }
      else { # protein sequence looks OK
        push @{$HoA_PI_lengths{$panID}}, $len;
        $geneID_seq{$display_id} = $seq;
        $geneID_len{$display_id} = $len;
      }
    }
  }
  else { # Option -nostop was not selected, so handle the sequence without evaluating stops
    push @{$HoA_PI_lengths{$panID}}, $len;
    $geneID_seq{$display_id} = $seq;
    $geneID_len{$display_id} = $len;
  }
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
  #say "$panID\t$count\t$min\t$median\t$max";
}

my $PREX;
if ($prefer){ $PREX = qr/$prefer/ }

# Traverse sequences again, reporting the first sequence per panID with median length
if ($verbose){say "count\tmin\tmedian\tmax\tlen\tverdict\tpanID__geneID" }
my %seen_PI;
my %verdict;
my %combo_IDs_verdict;
for my $combo_ID (sort keys %geneID_seq){

  my ($panID, $gene_ID) = split(/__/, $combo_ID);
  my $seq = $geneID_seq{$combo_ID};
  my $len = length($seq);

  if ($prefer){
    if ( $gene_ID =~ /$PREX/){ 
      if ( $len == $median_PI{$panID} ){
        $verdict{$combo_ID} = 1; 
      }
      else { # not equal to the median, so assign to 0; won't be used
        $verdict{$combo_ID} = 0;
      }
    }
    else { # $gene_ID not matching $PREX. May be used if a better one isn't seen.
      if ( $len == $median_PI{$panID} ){
        $verdict{$combo_ID} = rand(); # Not top-rated, but acceptable. 
      }
      else { # not equal to the median, so assign to 0; won't be used
        $verdict{$combo_ID} = 0;
      }
    }
  }
  else { # no $prefer string was provided
    if ( $len == $median_PI{$panID} ){
      $verdict{$combo_ID} = rand(); # Not top-rated, but acceptable. 
    }
    else { # not equal to the median, so assign to 0; won't be used
      $verdict{$combo_ID} = 0;
    }
  }
  #say "$panID\t$verdict{$combo_ID}\t$combo_ID";
  $combo_IDs_verdict{$combo_ID} = { verdict => $verdict{$combo_ID}, panID => $panID };
}

my $hash_ref = \%combo_IDs_verdict;

my @sorted = sort {
    $hash_ref->{$a}{panID} cmp $hash_ref->{$b}{panID}
    or
    $hash_ref->{$b}{verdict} <=>  $hash_ref->{$a}{verdict}
  } keys %$hash_ref;

my %seen_panID;
for my $combo_ID ( @sorted ){
  my ( $panID, $geneID ) = split(/__/, $combo_ID);
  if ( $seen_panID{$panID} ){ next }
  #say "$verdict{$combo_ID}\t$panID\t$combo_ID";
  my @seq_chunks = ( $geneID_seq{$combo_ID} =~ m/(.{1,100})/g);
  if ($outfile){
    say $OUT_FH ">$combo_ID\n", join("\n", @seq_chunks);
  }
  else { #print to stdout
    say ">$combo_ID\n", join("\n", @seq_chunks);
  }
  $seen_panID{$panID}++;
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
  $median = $sorted_values[ $mid ];
  return $median;
}

__END__
2021-11-15 Initial version.
2022-12-23 Rewrite: Take fasta input on STDIN, and take in an optional "preferred" ID regex
2022-12-25 Print to either stdout or out-file
2023-01-27 Add option to check for stop codons in protein sequence
2023-03-28 Avoid internal stops represented by '.', and strip terminal stop character (* and .)


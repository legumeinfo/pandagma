#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = <<EOS;
Generates a list of unique strings consisting of amino acid sequences, e.g. IVWF, PDIS, ...
The intended purpose is to generate codes that will be interpreted as amino acid
sequences by alignment programs -- for example, to encode a list of elements as peptide strings.
Glutamine and asparagine (Q and N) are not used; they have been reserved to indicate orientation.

Usage: make_peptide_hash.pl [options]

  OPTIONS:
    -number    (pos integer; default 100000) Number of unique strings to return.
    -width     (pos integer; default 4) Number of random AA characters to generate in each peptide.
    -unstable  (boolean; default unset) Generate peptides that are novel between runs. 
                 By default this is false (unset), and the random numbers are generated using the same 
                 random number seed each time, to give random-but-reproducible output.
    -seed      (integer; default 7) Seed used to generate random-but-stable numbers. Doesn't apply if -unstable.
    -help      This message. 
EOS

my $number = 100000;
my $seed = 7;  # Provided to srand unless -unstable is set
my $width = 4;
my $unstable;
my $help;

GetOptions (
  "number:i" => \$number,
  "seed:i" =>   \$seed,
  "width:i" =>  \$width,
  "unstable" => \$unstable,
  "help" =>     \$help,
);

die "\n$usage\n" if ( $help );

unless ($unstable){
  srand($seed);
}

# Coding groups
my @ZZ = split("", "AVLIFWPMSCEYDHKRGT");
my @FF = split("", "QN");

my %seen_motif;
my $code = "";
my $count = 0;
my $prefix="pan";

foreach my $i (1 .. $number){
  while (1){
    $code = make_motif($width);
    if ($seen_motif{$code}){
      next;
    }
    else {
      $count++;
      my $panstr = sprintf("$prefix%05d", $count);
      print "$panstr\t$code\n";
      $seen_motif{$code}++;
      last;
    }
  }
}

sub make_motif {
  my $width = shift;
  my @AAs = ();
  for my $posn ( 0 .. $width-1 ){
    my $RN = int(rand(18));
    $AAs[$posn] = $ZZ[$RN];
  }
  return join("", @AAs);
}

__END__
2023 S. Cannon
02-09 Initial version.
02-10 Drop coding groups. Use 18 AAs for the motif, reserving Q/N for punctuation & orientation.
02-11 Add srand and width params. Restructure make_motif subroutine.

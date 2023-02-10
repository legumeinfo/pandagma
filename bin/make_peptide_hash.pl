#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = <<EOS;
Generates a list of unique strings consisting of amino acid sequences, e.g. IVWFS, PDISE, ...
The intended purpose is to generate codes that will be interpreted as amino acid
sequences by alignment programs -- for example, to encode a list of elements as peptide strings.
Glutamine and asparagine (Q and N) are not used; they have been reserved to indicate orientation.

Usage: make_peptide_hash.pl [options]

  OPTIONS:
    -number   Number of unique strings to return [100000]
    -help     This message. 
EOS

my $number = 100000;
my $help;

GetOptions (
  "number:i" => \$number,
  "help" =>     \$help,
);

die "\n$usage\n" if ( $help );

# Coding groups
my @ZZ = split("", "AVLIFWPMSCEYDHKRGT");
my @FF = split("", "QN");

my %seen_motif;
my $code = "";
my $count = 0;
my $prefix="pan";

foreach my $i (1 .. $number){
  while (1){
    $code = make_motif();
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
  my $rA = int(rand(11));
  my $rB = int(rand(12));
  my $rC = int(rand(11));
  my $rD = int(rand(12));
  my $rE = int(rand(11));
  #my $rF = int(rand(2));
  my $code = join("", $ZZ[$rA], $ZZ[$rB], $ZZ[$rC], $ZZ[$rD], $ZZ[$rE]);
  return $code;
}

__END__
2023 S. Cannon
02-09 Initial version.
02-10 Drop idea of coding groups for the motif positions. Use 18 AAs for the motif, reserving Q/N for punctuation & orientation.

#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = <<EOS;
Generates a list of unique strings consisting of amino acid sequences, e.g. IVWFS, PDISE, ...
The strings are peptide-like, with alternating nonpolar (mostly) and polar (mostly) residues.
The intended purpose is to generate codes that will be interpreted as amino acid
sequences by alignment programs -- for example, to encode a list of elements as peptide strings.
Glutamine and tyrosine (Q and Y) are not used; they have been reserved to indicate orientation.

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
my @AA = split("", "AVLIFWPMSCE");
my @BB = split("", "SYEDHKRTFLVG");
my @CC = split("", "ICLAWEFMSVP");
my @DD = split("", "EDFLHYRVKTSG");
my @EE = split("", "IFEVCWLMPSA");
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
  my $code = join("", $AA[$rA], $BB[$rB], $CC[$rC], $DD[$rD], $EE[$rE]);
  return $code;
}

__END__
2023 S. Cannon
02-09 Initial version.


	Coding groups
count
1	Alanine	Ala	A	1	nonpolar
2	Valine	Val	V	2	nonpolar
2	Leucine	Leu	L	3	nonpolar
1	Isoleucine	Ile	I	4	nonpolar
2	Phenylalanine	Phe	F	5	nonpolar, aromatic
1	Tryptophan	Trp	W	6	nonpolar, aromatic
1	Proline	Pro	P	7	nonpolar
1	Methionine	Met	M	8	nonpolar
2	Serine	Ser	S	9	polar
2	Cysteine	Cys	C	10	polar
	Glutamic acid	Glu	E	11	acidic

2	Valine	Val	V	1	nonpolar
2	Leucine	Leu	L	2	nonpolar
2	Serine	Ser	S	3	polar
1	Tyrosine	Tyr	Y	4	polar, aromatic
2	Glutamic acid	Glu	E	5	acidic
1	Aspartic acid	Asp	D	6	acidic
1	Histidine	His	H	7	basic
1	Lysine	Lys	K	8	basic
1	Arginine	Arg	R	9	basic
1	Glycine	Gln	Q	10	polar
1	Threonine	Thr	T	11	polar
2	Cysteine	Cys	C	12	polar
2	Phenylalanine	Phe	F	13	nonpolar, aromatic

1	Glutamine	Gln	Q	1	polar
1	Asparagine	Asn	N	2	polar


#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use feature "say";

my ($help, $total);
my $divisor = 1;
my $min_bin_ct = 1;
GetOptions ("divisor:i"    => \$divisor, 
            "min_bin_ct:i" => \$min_bin_ct,
            "total"        => \$total,
            "help"         => \$help,
            );

my $usage = <<EOS;
Converts a two-column histogram (bins, counts) to an ascii plot of count values.

Usage: cat stats* | pan_histplot.pl -divisor 150

    -divisor     Optional divisor to apply to histogram values, to lower the plot amplitude [1]
    -min_bin_ct  Minimum count for printing a bin [1]
    -total       (boolean; false) Flag indicates to display the bin * bin_count -- i.e. the gene count at that bin
    -help        (boolean) This message
EOS

die "$usage\n" if $help; 

my $header = "bin\tabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRST";

print "\nHistograms of cluster sizes, with amplitude divisor = $divisor and minimum bin count = $min_bin_ct\n";

while (my $line = <>) {
  chomp $line;
  if ($line =~ /^Counts/){ print "\n$line\n$header\n"; }
  elsif ($line =~ /^(\d+)\t(\d+)/){
    my ($bin, $value) = ($1, $2);
    my $printbin = sprintf("%i", $bin);
    if ($value >= $min_bin_ct){
      if ($total){
        say "$printbin\t", "." x ($value*$bin/$divisor);
      }
      else {
        say "$printbin\t", "." x ($value/$divisor);
      }
    }
  }
  else {}
}

__END__
VERSIONS
2023-02-03 Initial version for pandagma based on histplot and the one-liner perl -ane 'print "$F[0]\t", "." x $F[1], "\n"'
2023-08-08 Add option to print total per bin -- i.e. the gene count at that bin
2023-09-09 Extend header row to 100 characters


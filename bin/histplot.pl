#!/usr/bin/perl
## #!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $help;
my $divisor = 1;
GetOptions ("divisor:s"  => \$divisor, 
            "help"       => \$help
            );

my $usage = <<EOS;
Converts a two-column histogram (bins, counts) to an ascii plot of count values.

Usage: cat STRING_OF_NUMBERS | histogram -n -s1 | histplot.pl -divisor 10

    -divisor   Optional divisor to apply to histogram values, to lower the plot amplitude
                 Default: 1 (no modification of values)
    -help      This message
EOS

die "$usage\n" if $help; 

BEGIN{print "bin\tabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRST\n"}

while (<>) {
  my ($bin, $value) = split(/\t/, $_);
  my $printbin = sprintf("%.2f", $bin);
  print "$printbin\t", "." x ($value/$divisor), "\n";
}

__END__
VERSIONS
sc = Steven Cannon
v01 2017-08-15 initial, based on the one-liner perl -ane 'print "$F[0]\t", "." x $F[1], "\n"'
v02 2017-01-29 use printf to format the bin output


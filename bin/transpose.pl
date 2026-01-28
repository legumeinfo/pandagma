#!/usr/bin/env perl
use strict; 
use warnings;

# NAME
#   transpose.pl - Given tabular data on STDIN, return transposed table on STDOUT
#
# SYNOPSIS
#   cat TABLE | transpose.pl
#
# OPERANDS
#   Tabular data on STDIN

my(%data);  
my($maxcol, $rownum) = (0, 0);
while (<>) {
  my(@row) = split /\s+/;
  my($colnum) = 0;
  foreach my $val (@row) {
    $data{$rownum}{$colnum++} = $val;
  }
  $rownum++;
  $maxcol = $colnum if $colnum > $maxcol;
}

my $maxrow = $rownum;
for (my $col = 0; $col < $maxcol; $col++) {
  for (my $row = 0; $row < $maxrow; $row++) {
    printf "%s%s", ($row == 0) ? "" : "\t",
        defined $data{$row}{$col} ? $data{$row}{$col} : "";
  }
  print "\n";
}


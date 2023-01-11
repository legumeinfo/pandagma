#!/usr/bin/env perl
# rename.pl -- change filenames
use strict;
use warnings;

my $usage = <<EOS;
Usage: perl $0 pattern file(s)
e.g. 
  perl $0 's/\.orig//' *.orig
  perl $0 'y/A-Z/a-z/ unless /^Make/' *
  perl $0 's/MtChr/Mt/' *
EOS

my $op = shift or die "\n$usage\n";

for (@ARGV) {
  my $was = $_;
  eval $op;
  die if $@;
  rename($was,$_) unless $was eq $_;
}

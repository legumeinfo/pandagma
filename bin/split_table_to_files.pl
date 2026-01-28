#!/usr/bin/env perl 

use strict;
use warnings;
use Getopt::Long;

my ($out_dir, $keep, $help);
my $field_num = 1;

GetOptions (
  "out_dir=s"     => \$out_dir,
  "field_num:i"   => \$field_num,
  "keep"          => \$keep,
  "help"          => \$help
);

my $usage = <<EOS;
Usage: cat TABLE | split_table_to_files.pl -out out_dir  [options]
  Input is tabular data, passed via STDIN, to be split into sub-tables, 
  with the sub-tables named by values in the sorted input table.
  
  Required:
  [data on STDIN]
  -out_dir   Name of directory for output

  Optional
  -field_name   Number of field to use for naming output files (one-indexed) Default: 1
                Values in this field should be suitable for file-naming, 
                  and the file should be sorted by this column.
  -keep         Keep field_name in the output files
  -help         Displays this info
EOS

die "\n$usage\n" unless ($out_dir || $help);

opendir (DIR, "$out_dir") or mkdir($out_dir); close DIR;

my %seen;

while (my $line = <>) {
  chomp $line;  
  next if $line =~ m/^#|^$/;
  
  my @elements = split /\t/, $line;
  my $file_out = $elements[$field_num-1];
  my @rest;
  if ($keep){ # Keep the $field_num column in the output
    @rest = @elements;
  }
  else { # Remove the $field_num column
    @rest = splice(@elements,$field_num);
  }
  
  if ($seen{$file_out}) {
    print OUT join("\t", @rest),"\n";
  }
  else {
    $seen{$file_out}++;
    my $out_file = "$out_dir/$file_out";
    open (OUT, "> $out_file") or die "can't open out $out_file: $!";
    #print OUT "$line\n";
    print OUT join("\t", @rest),"\n";
  }
}

__END__

# Versions
2008-01-08 basic; works OK.
2008-02-25 Change @ARGV input test to allow $field_num = 0
2008-04-04 Change from zero-indexing in the field number to one-indexing
2008-11-06 Add chomp, in case file name is taken from last field
2019-12-14 Add Getopt, and don't print filename in output tables
2023-11-07 Add option -keep to allow keeping the field_name in the output files
2026-01-27 Take tabular data in via STDIN rather than as a named file

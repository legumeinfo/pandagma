#!/usr/bin/env perl

use strict;
use warnings;
use feature "say";
use Getopt::Long;
use JSON;

my $usage = <<EOS;
  Synopsis: pangene_search.pl -pan_file PANGENE_FILE [-format plain]  QUERY_STRING

  Given a search string (either a gene name or pangene ID), return the pangene ID 
  that contains that string, and the constituent gene IDs in that pangene set.
  The output is to stdout, in JSON format, as a hash (with ID as key) of an array (holding the pangene set).
  The file containing the pangene sets should be tab-separated, with the pangene ID 
  in the first column and one or more gene names in second and subsequent fields.

  Examples:
      pangene_search.pl -pan_file glysp.mixed.pan2.TV81.clust.tsv glyma.Lee.gnm1.ann1.GlymaLee.06G127600

    A default filename may be provided, so this should work if the file is available:
      pangene_search.pl glyma.Lee.gnm1.ann1.GlymaLee.06G127600

    A "non-fully-qualified" gene name or pangene ID should work, but may result in multiple matching records.
      pangene_search.pl GlymaLee.06G127600
      
  Required: 
  query_string  Search string, passed in via ARGV[0]
  
  Options:
  -pan_file  Name of pangene file  (a file is required, but a default may be set in the program).
  -format    clust, hash, json (default "clust"). 
               "clust" is a list: panID geneID geneID ...
               "hash" is a two-column format: panID geneID / panID geneID / panID geneID / ...
               "json" is similar to "clust", but in a json object.
  -help      This message.
  
EOS

my ($pan_file, $query, $help);
my $format = "clust";

$pan_file="glysp.mixed.pan2.TV81.clust.tsv"; # Default pangene set

GetOptions (
  "pan_file:s" => \$pan_file,   
  "format:s"   => \$format,
  "help"       => \$help,
);

die "$usage\n" unless ($ARGV[0]);
die "$usage\n" if ($help);

if ($format !~ /clust|hash|json/){
  die "Unrecognized format \"$format\". Allowed formats are: clust, hash, json (default \"clust\").\n";
}

$query=$ARGV[0];

# Trim off splice variant suffix, if present
$query =~ s/\.\d+$//;

open(my $FH, '<', $pan_file) or die "can't open pangene pan_file $pan_file: $!";
my $regex = qr/$query/i;
my $found_it = 0;
while (my $line = <$FH>) {
  if ($line =~ /$regex/) {
    chomp $line;
    my @set=split(/\t/, $line);
    my $ID=$set[0];
    my $size=scalar(@set)-1;
    my @genes=@set[1..$size];
    if ($format =~ /json/){
      my $gene_set_HoA = { $ID => \@genes };
      say encode_json $gene_set_HoA;
    }
    elsif ($format =~ /clust/) { 
      print $ID;
      for my $geneID (@genes){
        print "\t$geneID";
      }
      print "\n";
    }
    elsif ($format =~ /hash/) {
      for my $geneID (@genes){
        say "$ID\t$geneID";
      }
    }
    else {
      die "Unrecognized format: $format\n";
    }
  }
  $found_it = 1;
}
unless ($found_it == 1){ say "WARNING: $query not found" }

__END__

# Steven Cannon 

Versions
2020-02-21 Initial version.
2020-02-29 Report multiple families if found.
2022-02-04 Offer three output formats. Add as a utility to pandagma.

A partial name may find multiple pangene sets, e.g. 06G127600 (from Lee, Zh13, GlysoPI, Wm82.gnm2) 
... generally because unqualified names (i.e. name fragments) may identify different genes 
in different annotations.


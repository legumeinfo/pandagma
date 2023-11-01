#!/usr/bin/env perl

# Program: pangene_tabularize.pl
# Program description: see usage message.
# Steven Cannon 2023

use strict;
use warnings;
use DB_File;
use Getopt::Long;
use feature "say";

my $pan_to_gn = "pan_to_gn.db";
my $hash_method = "bdb";
my $to_annot = "all";
my $json = 0;
my $annot_list;
my $help;

my $usage = <<EOS;
  Usage:  pangene_tabularize.pl -annot_list FILE [-options]
  
  Given a stream of accession identifiers and a pangene file, return a table of genes in the order of the accession identifiers.
  
  The pangene file should be structured thus, as a Berkeley DB file**:
    panID 
    gene1,gene2,gene3,... 
    panID
    gene1,gene2,gene3,... 

  ** The other supported alternative is plain-text, tab-separated:
    panID gene1 gene2 gene3 ... 
    panID gene1 gene2 gene3 ... 

  Required:
    -annot_list  List of annotation prefixes, e.g. 
                   glyma.Wm82.gnm1.ann1
                   glyma.Wm82.gnm2.ann1
                   glyma.Wm82.gnm4.ann1

  Options: 
    -pan_to_gn   [pan_to_gn.db] Path to hash file of pangeneIDs and the list of all geneIDs in 
                   that pangene, either in BerkeleyDB or tsv hash format. 
                   In either case, the gene list should be comma-separated -
                   the panID as the key and the value as the collection of genes.
                   For the tsv option:
                     panID  gene1,gene2,gene3,...
                   For the BerkeleyDB option, the key and value alternate:
                     panID
                     gene1,gene2,gene3,...
    -hash_method [bdb] bdb or tsv
    -help        This message.
EOS

GetOptions(
  'annot_list=s'  => \$annot_list,  # required
  'pan_to_gn:s'   => \$pan_to_gn,   # required - but default is set as pan_to_gn.db
  'hash_method:s' => \$hash_method,
  'h|help'        => \$help,
);

die "$usage\n" if $help;
die "$usage\n" unless $annot_list;
if ( $hash_method =~ /bdb|BerkeleyDB/i && $pan_to_gn !~ /\.db$/ ){
  die "The -hash_method was set as $hash_method but the -pan_to_gn file doesn't have a suffix of .db." .
      "Please check that the BerkeleyDB file exists. If it needs to be created, try gene_trans_db_load.sh\n";
}
elsif ( $hash_method =~ /tsv/ && $pan_to_gn !~ /\.tsv.gz$|\.tsv$/ ){
  die "The -hash_method was set as BerkeleyDB but the -pan_to_gn file, [$pan_to_gn],\n" .
      "doesn't have a suffix of .tsv or .tsv.gz. If it needs to be created, try gene_trans_db_load.sh\n";
}

open my $ANN_FH, "<", $annot_list or die "Can't open in annot_list $annot_list: $!\n";

my %pan_to_gn_hsh;
my $P2G_FH;
if ($hash_method =~ /bdb|BerkeleyDB/i){
  tie %pan_to_gn_hsh, 'DB_File', $pan_to_gn, O_RDWR|O_CREAT, 0666, $DB_HASH
      or die "Cannot open pan_to_gn $pan_to_gn: $!\n";
}
elsif ($hash_method =~ /tsv/){
  if ( $pan_to_gn =~ /gz$/ ){
    open $P2G_FH, "zcat < $pan_to_gn |" or die "Can't open pan_to_gn $pan_to_gn:$!";
  }
  else {
    open $P2G_FH, "<", $pan_to_gn or die "Can't open pan_to_gn $pan_to_gn:$!";
  }
  while (<$P2G_FH>){
    chomp;
    my ($pan, @genes_in_pangene) = split("\t", $_);

    # Next: strip splice variant from mRNA IDs to get gene IDs
    my @genes_in_pangene_stripped = map { local $_ = $_; s/\.\D*\d+$//; $_ } @genes_in_pangene;
    my $pangene_string = join(",", @genes_in_pangene_stripped);
    $pan_to_gn_hsh{$pan} = $pangene_string;
    # say join("\t", $pan, scalar(@genes_in_pangene), $genes_in_pangene[0], $genes_in_pangene_stripped[1]);
  }
  close $P2G_FH;
}
else { die "Flag -hash_method must be either BerkeleyDB or tsv.\n" }

my @annots;
while (<$ANN_FH>){
  chomp;
  my $annot = $_;
  push @annots, $annot;
}
print join("\t", "#pangene", @annots), "\n";

my $ct = 0;
my $test_limit = 9999999 ; # Set to a small positive integer to examine just that number of pangenes
while ( my ($pan, $pangenes) = each( %pan_to_gn_hsh ) ) {
  my @pangene_ary = split(/,/, $pangenes);
  my @target_genes_in_pangene;
  push @target_genes_in_pangene, $pan;
  for my $annot (@annots){
    my @annot_matches = grep(/$annot/ig, @pangene_ary);
    my @matches_stripped;
    if (scalar(@annot_matches) == 0){ $annot_matches[0] = "NONE" }
    for my $gene (@annot_matches){
      $gene =~ s/$annot\.//i;
      push @matches_stripped, $gene;
    }
    push @target_genes_in_pangene, join(",", @matches_stripped);
  }
  print join("\t", @target_genes_in_pangene);
  $ct++;
  last if $ct >= $test_limit;
  say "";
}
say "";

__END__

Versions
2023-11-01 Initial functioning version



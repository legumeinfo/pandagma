#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use feature "say";

my $usage = <<EOS;
  Synopsis:  hash_into_table_id.pl HASH_FILE(s) -table TABLE_FILE [options] 
  
  Read key-value pairs from one or more hash files as first argument(s), and a file with tabular data.
  Swap the IDs with the values in the two indicated columns of the tabular data with values from the hash.
  The hash can be created to match the full ID, including splice variants,
  OR created to match just the gene IDs, excluding splice variants.
  In the latter case, provide -splice_regex to give the form of the splice variant.
  
  Required:
  [ARGV list]  One or more files containing key-value hash, where first column has IDs from tabular file.
  -table  Name of file with tabular data.
  
  Options:
  -idx1     (integer) Column index for first column to swap (zero-indexed). Default 1 (col 2)
  -idx2     (integer) Column index for second column to swap (zero-indexed). Default 2 (col 3)
  -splice_regex   (string) regular expression to use to exclude the 
                  splice variant suffix of a feature name during the match. 
         Example 1: For transcripts like Gene1234.1,      use "\\.\\d+"  
         Example 2: For transcripts like Gene1234-mRNA-1, use "-mRNA-\\d+" 
         Example 3: For proteins like    Gene1234.1.p,    use "\\.\\d+\\.p"
  -help   (boolean) This message.
  
EOS

my ($table, $splice_regex, $help, $SPL_RX);
my ($idx1, $idx2) = (1, 2);

GetOptions (
  "table:s" =>         \$table,
  "idx1:i" =>          \$idx1,
  "idx2:i" =>          \$idx2,
  "splice_regex:s" =>  \$splice_regex,   
  "help" =>            \$help,
);

die "$usage\n" if ($help || !defined($table));
if (scalar(@ARGV)==0){ die "\nPlease provide one or more files with hash of IDs. See -help for guidance.\n\n" }

if ( $splice_regex ){ $SPL_RX=qr/$splice_regex/ }
else { $SPL_RX=qr/$/ }

# read in the hash
my %hash;

foreach my $hash_file ( @ARGV ) {
  open ( my $HSH, '<', $hash_file ) or die "Can't open in input_hash, $hash_file: $!";
  my $ct=0;
  while (<$HSH>) {
    chomp;
    /(\S+)\s+(.+)/;
    my ($id, $hash_val) = ($1,$2);
    $hash{$id} = $hash_val;   # print "$id, $hash{$id}\n";
    $ct++;
  }
  #say $ct;
  close $HSH;
}

# Read in the tabular data
open ( my $TABLE_FH, "<", $table ) or die "Can't open in table, $table: $!";
while ( <$TABLE_FH> ){
  chomp;
  my $line = $_;
  my @fields = split(/\t/, $line);

  my $ID1 = $fields[$idx1];
  my $ID2 = $fields[$idx2];

  my ($first_OK, $second_OK);

  # strip off splice variant for ID1
  $ID1 =~ m/(.+)($SPL_RX)$/;
  my ($base_id1, $suffix1) = ($1, $2);
  if ($splice_regex){
    $suffix1 =~ s/$SPL_RX//; 
    $suffix1 =~ s/-/./; 
  }
  
  if (defined $hash{$base_id1}){
    $fields[$idx1] = $hash{$base_id1};
    $first_OK = 1;
  }
  else {
    $hash{$base_id1} = "$base_id1 HASH UNDEFINED";
    warn "WARNING: HASH UNDEFINED for $base_id1\n";
    $first_OK = 0;
  }

  # strip off splice variant for ID2
  $ID2 =~ m/(.+)($SPL_RX)$/;
  my ($base_id2, $suffix2) = ($1, $2);
  if ($splice_regex){
    $suffix2 =~ s/$SPL_RX//; 
    $suffix2 =~ s/-/./; 
  }
  
  if (defined $hash{$base_id2}){
    $fields[$idx2] = $hash{$base_id2};
    $second_OK = 1;
  }
  else {
    $hash{$base_id2} = "$base_id2 HASH UNDEFINED";
    warn "WARNING: HASH UNDEFINED for $base_id2\n";
    $second_OK = 0;
  }

  if ($first_OK && $second_OK){
    say join("\t", @fields);
  }
  else {
    warn "WARN| New fields not OK. Field 1: $first_OK; Field 2: $second_OK\n"
  }
}

__END__

# Steven Cannon 

Versions
2023-08-27 New script, derived from hash_into_fasta_id.pl


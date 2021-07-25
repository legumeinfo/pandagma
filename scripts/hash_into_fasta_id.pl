#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $usage = <<EOS;
  Synopsis: hash_into_fasta_id.pl [options] -fasta FILE -hash FILE
  
  Read a key-value hash file and a fasta file. 
  Swap the IDs with the values from the hash file.
  
  Required:
  -fasta     (string) Fasta file
  -hash      (string) Key-value hash file, where first column has IDs from fasta file
  
  Options:
  -suff_regex (string) quoted portion of regular expression to use to exclude the last 
               part of a feature name during the match. 
               For example, if a gene "ASDF" has a splice variant "-mRNA-1" :
                 ASDF-mRNA-1    
               then supply "\\-mRNA-\\d+" ; then ASDF will be replaced by the hashed gene name.
  -keep_definition  (boolean) the defline in the fasta file will be retained in the output.
  -out_file  (string) write to this file; otherwise, to stdout.
  
  -help   (boolean) This message.
  
EOS

my ($fasta_file, $hash_file, $suff_regex, $out_file, $keep_definition, $help);

GetOptions (
  "fasta_file=s" =>    \$fasta_file,   # required
  "hash_file=s" =>     \$hash_file,   # required
  "out_file:s" =>      \$out_file,   
  "suff_regex:s" =>    \$suff_regex,   
  "keep_definition" => \$keep_definition,
  "help" =>            \$help,
);

die "$usage\n" unless (defined($fasta_file) and defined($hash_file));
die "$usage\n" if ($help);

my $REX = "$suff_regex";
 #print "[$suff_regex] [$REX]\n";

# read hash in
open(my $HSH, '<', $hash_file) or die "can't open in input_hash, $hash_file: $!";
my %hash;
while (<$HSH>) {
  chomp;
  /(\S+)\s+(.+)/;
  my ($id, $hash_val) = ($1,$2);
  $hash{$id} = $hash_val;   # print "$id, $hash{$id}\n";
}

# Read in the sequence using the Bioperl SeqIO;
my $in;
if ($fasta_file =~ /.gz$/) {
  open my $fh, "gzip -dc $fasta_file |" or die $!;
  $in = Bio::SeqIO->new(-fh => $fh, -format => 'Fasta');
} else {
  $in = Bio::SeqIO->new(-file => $fasta_file, -format => 'Fasta');
}

my $OUT;
if (defined($out_file)){
  open ($OUT, '>', $out_file) or die "can't open out out_file, $hash_file: $!";
}

# Load the sequence into a Bio::Seq object
while ( my $seq = $in->next_seq ) {

  # first part of the fasta def line
  my $display_id = $seq->display_id();

  # strip off splice variant
  my $base_id = $display_id;
  my $suffix = "";
  if (defined($suff_regex)) {
    $display_id =~ m/(.+)($REX)$/;
    $base_id = $1;
    $suffix = $2;
    #print "\n[$base_id] [$suffix]\n"
  }
  else {} # do nothing; $base_id = $display_id
  
  # rest of fasta def line
  my $desc = $seq->desc();
  $desc = "" unless defined ($desc);
  
  # sequence
  my $sequence = $seq->seq();
  
  if ( $keep_definition ) { # keep the preexisting fasta description
      $hash{$base_id} = "$base_id HASH UNDEFINED" unless defined ($hash{$base_id});
      if (defined($out_file)) {
        print $OUT ">$hash{$base_id}$suffix $desc\n$sequence\n";
      }
      else {
        print ">$hash{$base_id}$suffix $desc\n$sequence\n"
      }
  }
  else { # NOT $keep_definition; don't keep the preexisting fasta description
      $hash{$base_id} = "$base_id HASH UNDEFINED" unless defined ($hash{$base_id});
      if (defined($out_file)) {
        print $OUT ">$hash{$base_id}$suffix\n$sequence\n";
      }
      else {
        print ">$hash{$base_id}$suffix\n$sequence\n";
      }
  }
}
print "\n";

__END__

# Steven Cannon 

Versions
v01 2014-05-21 New script, derived from hash_into_fasta_description.pl
v02 2018-02-09 Handle suffixes (e.g. for splice variants)
v03 2019-05-07 Print original ID if no hash is found

#!/usr/bin/env perl

# PROGRAM: get_superfam_from_family_file.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2024; derived from get_fasta_from_family_file.pl
# see description under Usage

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($family_file, $in_dir, $out_dir, $add_IDs, $help);

GetOptions (
  "family_file=s" => \$family_file, # required
  "in_dir=s" =>      \$in_dir,      # required
  "out_dir=s" =>     \$out_dir,     # required
  "add_IDs:s" =>     \$add_IDs,  
  "help" =>          \$help
);

my $scriptname = basename($0);

my $usage = <<EOS;
  Usage: $scriptname -family_file FILENAME -in_dir DIRNAME -out_dir DIRNAME [-options]
  
  Read a file of superfamily membership, consisting of rows each containing a list of gene families;
  and optionally, an initial value containing a name for the superfamily in a family.
  The row can begin with a superfamily ID, or unique IDs can be assigned as an option by this program.
  Print sequences from each family to new fasta files, one per superfamily.
  Format of family file is like this:
    #Fam_ID gene_IDs
    ...
    Superfam_998  Legume.fam3.19708 Legume.fam3.04638 Legume.fam3.17076 Legume.fam3.19540
    Superfam_999  Legume.fam3.19735 Legume.fam3.04436 Legume.fam3.07027 Legume.fam3.12077
    Superfam_1000 Legume.fam3.19745 Legume.fam3.14543 Legume.fam3.14843 Legume.fam3.18466
   
   -in_dir      DIRNAME   ** Directory with family fasta files (to be merged into the superfamily files)
   -out_dir     DIRNAME   ** Directory for output (to contain one file per line of family_file)
   -family_file FILENAME  ** File with family ID and list of sequence IDs, described above
   -add_IDs     (string)  -- Supply name for each family, as stringNUM, where NUM is [1..family_count]
   -help        (boolean) -- Display this help
   
                          ** = required
EOS

die "\n$usage\n" if ($help or !defined($in_dir) or !defined($out_dir) or !defined($family_file) );

# Read in family file
open( my $FAM_FH, '<', $family_file ) or die "can't open family_file $family_file: $!";

# Traverse family file and write one fasta file per family
while (<$FAM_FH>) {
  chomp;
  next if ( $_ =~ /^$/ or $_ =~ /^#/);

  my ($superfam_ID, @seq_IDs);
  my (@items) = split(/\s+/, $_);
  if ($add_IDs){
    my $fam_num = sprintf("%05d", $.);
    $superfam_ID = "$add_IDs$fam_num";
    @seq_IDs = @items;
  }
  else {
    $superfam_ID = shift(@items);
    @seq_IDs = @items;
  }

  # Open file for writing new fasta out
  system("cat /dev/null > $out_dir/$superfam_ID");
  
  print "Writing $superfam_ID\n";
  for my $seq_ID (@seq_IDs){
    system("cat $in_dir/$seq_ID >> $out_dir/$superfam_ID");
  }
}

__END__
VERSIONS

2024-12-28 Initial version, based on get_fasta_from_family_file.pl

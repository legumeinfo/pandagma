#!/usr/bin/env perl

# PROGRAM: get_fasta_from_family_file.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2017
# see description under Usage

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

my ($family_file, $out_dir, $add_IDs, $help);

GetOptions (
  "family_file=s" => \$family_file,   # required
  "out_dir=s" =>     \$out_dir,   # required
  "add_IDs:s" =>     \$add_IDs,  
  "help" =>          \$help
);

my $scriptname = basename($0);

my $usage = <<EOS;
  Usage: $scriptname FASTA_FILE[S] -family_file FILENAME -out_dir DIRNAME [-options]
  
  Read a file of gene family membership, consisting of rows each containing a list of genes in a family.
  The row can begin with a family ID, or unique IDs can be assigned as an option by this program.
  Print sequences from each family to new fasta files, one per family.
  Format of family file is like this:
    #Fam_ID gene_IDs
    L_v7lGl cicar.Ca_26226  cicar.Ca_27089  lotja.Lj2g3v2291660.1 lotja.Lj0g3v0113769.1
    L_V7ll5 glyma.Glyma.01G021500 glyma.Glyma.09G200100 phavu.Phvul.002G145500  vigun.Viungv11005429m
    L_V7mqj glyma.Glyma.01G098900 glyma.Glyma.02G175600 medtr.Medtr8g079750.1 medtr.Medtr4g062590.1
    L_V7P9b medtr.Medtr3g449470.1 medtr.Medtr6g029430.1 medtr.Medtr6g012600.1 medtr.Medtr7g045640.1
   
   -out_dir     DIRNAME   ** Directory for output (to contain one file per line of family_file)
   -family_file FILENAME  ** File with family ID and list of sequence IDs, described above
   -add_IDs     (string)  -- Supply name for each family, as stringNUM, where NUM is [1..family_count]
   -help        (boolean) -- Display this help
   
                          ** = required
EOS

die "\n$usage\n" 
  if ($help or scalar(@ARGV)==0 or !defined($out_dir) or !defined($family_file) );

# Read fasta file(s) into ID-sequence hash
my %seq_hsh;
foreach my $input_fas ( @ARGV ) {
  chomp($input_fas);
  print "    $input_fas\n";
  
  my $seqio_obj;
  if ($input_fas =~ /.gz$/) {
    open my $fh, "gzip -dc $input_fas |" or die $!;
    $seqio_obj = Bio::SeqIO->new(-fh => $fh, -format => 'Fasta');
  } else {
    $seqio_obj = Bio::SeqIO->new(-file => $input_fas, -format => 'Fasta');
  }

  while ( my $seq_obj = $seqio_obj->next_seq ) {
    my $display_id = $seq_obj->display_id();
    $seq_hsh{$display_id} = $seq_obj->seq();
  }
}

# Read in family file
open( my $FAM_FH, '<', $family_file ) or die "can't open family_file $family_file: $!";

# Traverse family file and write one fasta file per family
while (<$FAM_FH>) {
  chomp;
  next if ( $_ =~ /^$/ or $_ =~ /^#/);

  my ($fam_ID, @seq_IDs);
  my (@items) = split(/\s+/, $_);
  if ($add_IDs){
    my $fam_num = sprintf("%05d", $.);
    $fam_ID = "$add_IDs$fam_num";
    @seq_IDs = @items;
  }
  else {
    $fam_ID = shift(@items);
    @seq_IDs = @items;
  }

  # Open file for writing new fasta out
  open( my $OUT_FH, '>', "$out_dir/$fam_ID" ) or die "can't open out $out_dir/$fam_ID: $!"; 
  
  for my $seq_ID (@seq_IDs){
    if (defined($seq_hsh{$seq_ID})) {
      print $OUT_FH ">$seq_ID\n$seq_hsh{$seq_ID}\n";
    }
    else {
      warn "WARNING: can't find sequence for $seq_ID\n";
    }
  }
}

__END__
VERSIONS

v01 2017-07-02 Initial version
v02 2021-07-09 Takes fasta FILES in via stdin. Add option for program to assign fam IDs.


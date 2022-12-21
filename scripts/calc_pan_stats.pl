#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum min max);

my $usage = <<EOS;
Given pan-gene results in "clust.tsv" format (pan-gene ID in first column, followed by
gene IDs, tab-separated), and a regex for capturing assembly/annotation string from
the gene IDs, report modal orthogroup (OG) size and composition per annotation-set.

Usage: calc_pan_stats.pl -pan FILE.clust.tsv [options]

  REQUIRED: 
    -panfile      Pan-gene cluster file

  OPTIONS:
    -annot_regex  Regular expression for capturing annotation name from prefixed gene ID. Default is
                  \"([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+\" for four dot-separated fields (gensp.accn.gnm1.ann1.GENEID)
                  Or, for Zea assembly+annot string (Zm00032ab.GENEID): \"(\\D+\\d+\\D+)\\d+.+\" 
    -out_matrix   Output filename for matrix of counts per pan-ID and annotation
    -help         This message. 
EOS

my ($panfile, $out_matrix, $help);
my $annot_regex = '([^.]+\.[^.]+\.[^.]+\.[^.]+)\.\S+';

GetOptions (
  "panfile=s" => \$panfile,
  "annot_regex:s" => \$annot_regex,
  "out_matrix:s" => \$out_matrix,
  "help" =>      \$help,
);

die "\n$usage\n" if ( $help or ! $panfile );

my $REX = qr/$annot_regex/;

open (my $PAN, "<", $panfile) or die "Can't open in $panfile: $!\n";

my $OUTMTX;
if ($out_matrix){
  open ($OUTMTX, ">", $out_matrix) or die "Can't open out $out_matrix: $!\n";
}

my %annots;
my $ct_annot_global;
while (<$PAN>){
  chomp;
  my ($first, @rest) = split(/\t/, $_);
  my $ct_this_og=0;
  my %seen; 
  foreach my $gene (@rest){
    $gene =~ $REX;
    my $ann=$1;
    unless($seen{$ann}){$seen{$ann}++; $ct_this_og++;}
    unless($annots{$ann}){$annots{$ann}++; $ct_annot_global++;}
  }
}

my %annot_cts;
my @annot_names = ( sort(keys %annots) );

seek $PAN,0,0;
&printstr("#pan_ID\t", join("\t", @annot_names), "\n");
while (my $line = <$PAN>){
  chomp $line;

  my @this_line = split(/\t/, $line);
  &printstr($this_line[0], "\t");
  for my $annot (@annot_names){
    my @found = grep(/$annot/, @this_line);
    &printstr(scalar(@found), "\t");
    $annot_cts{$annot}++ if (scalar(@found)>0);
  }
  &printstr("\n");
}

for my $annot (sort keys %annot_cts){
  print $annot, "\t", $annot_cts{$annot}, "\n";
}

#####################
sub printstr {
  my $str_to_print = join("", @_);
  if ($out_matrix) {
    print $OUTMTX $str_to_print;
  }
  else {
    print $str_to_print;
  }
}


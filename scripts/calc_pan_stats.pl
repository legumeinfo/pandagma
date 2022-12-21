#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

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

# Get list of all annotation names
my %annots; # Hash (list) of all annotation names
while (<$PAN>){
  chomp;
  my ($pan_id, @genes) = split(/\t/, $_);
  my $ct_this_og=0;
  foreach my $gene (@genes){
    $gene =~ /$REX/;
    my $ann=$1;  # annotation name, captured with $annot_regex
    unless($annots{$ann}){$annots{$ann}++}
  }
}

my %annot_cts;
my %gene_cts_per_annot;
my @annot_names = ( sort(keys %annots) );

# Get gene counts per annotation per OG
seek $PAN,0,0;
&printstr("#pan_ID\tanns_in_og\tgenes_in_og\t", join("\t", @annot_names), "\n");
while (my $line = <$PAN>){
  chomp $line;

  my @parts = split(/\t/, $line);
  my $anns_in_og = 0;
  my $genes_in_og = 0;
  my $annot_cts_in_og = "";
  for my $annot (@annot_names){
    my $ct_per_gn_per_og = scalar(grep(/$annot/, @parts));
    if ($ct_per_gn_per_og>0){ $anns_in_og++; $genes_in_og += $ct_per_gn_per_og }
    $annot_cts_in_og = $annot_cts_in_og . $ct_per_gn_per_og . "\t";
    if ( $ct_per_gn_per_og>0 ){
      $annot_cts{$annot}++;  # increments for this annot if there are any genes in this og
      $gene_cts_per_annot{$annot} += $ct_per_gn_per_og; # accumulates total gene count 
    }
  }
  &printstr("$parts[0]\t$anns_in_og\t$genes_in_og\t$annot_cts_in_og\n");
}

# Print OGs per annotation set
print "#Annotation\tannot_cts\ttotal_genes\n";
for my $annot (sort keys %annot_cts){
  print "$annot\t$annot_cts{$annot}\t$gene_cts_per_annot{$annot}\n";
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

__END__
S. Cannon
2022-12-21  Initial version




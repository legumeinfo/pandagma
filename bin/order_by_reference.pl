#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Parallel::ForkManager;
use feature "say";

my $usage = <<EOS;
Given a file with panIDs, geneIDs, annotation names, and positional information,
and the name of a prefered annotation to use as a reference, return the order
of panIDs following the order of genes in the preferred annotation.

Usage: order_by_reference.pl -pan_table PANGENE_TABLE -pref_annot STRING

  PANGENE_TABLE, e.g. e.g. 22*_posn.hsh.tsv, has six columns:
    panID  geneID  annot.chr  start  end  orientation

  REQUIRED:
    -pan_table  PANGENE_TABLE
    -pref_annot  STRING, e.g. 'G19833.gnm2.ann1' (bean) or 'Zm00001eb' (maize)

  OPTIONS:
    -consen_out  File with panID ordering, based on gene order from the consensus annotation.
    -unplaced_out File with panIDs missing from the consensus annotation.
    -annot_regex  Regular expression for capturing annotation name from gene ID, e.g.
                    \"([^.]+\\.[^.]+\\.[^.]+\\.[^.]+)\\..+\"
                      for four dot-separated fields, e.g. vigan.Shumari.gnm1.ann1 (default)
                    or \"(\\D+\\d+\\D+)\\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    -verbose   (boolean) For some debugging info, to STDOUT. Use -v -v to give more output.
    -help      (boolean) This message. 
EOS

my $pan_table;
my $pref_annot;
my $help;
my $consen_out;
my $unplaced_out;
my $verbose=0;
my $logstr="";
my $annot_regex = "([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+"; 

GetOptions (
  "pan_table=s" =>    \$pan_table,
  "pref_annot=s" =>   \$pref_annot,
  "consen_out:s" =>   \$consen_out,
  "unplaced_out:s" => \$unplaced_out,
  "annot_regex:s" =>  \$annot_regex,
  "verbose+" =>       \$verbose,
  "help" =>           \$help,
);

die "\n$usage\n" if ( $help || ! $pan_table || ! $pref_annot );

$pref_annot =~ s/['"]//g;
my $PREF_REX = qr/$pref_annot/;

# Do a first-pass reading of the pan_table to get counts per molecule, to help later
# bypass probable scaffolds.
# Read pan_table into an array, and get count of each molecule (chromosome or scaffold.)
open (my $PAN_FH, "<", $pan_table) or die "Can't open in pan_table: $pan_table\n";
my $chr_gene_count = 0;
my %chr_hsh;
while (<$PAN_FH>) {
  chomp;
  next unless (/^\S+/);
  my $line = $_;
  my ($panID, $gene, $ann_chr, $start, $end, $orient) = split(/\t/, $line);

  # From the third field, a annot.chr string, extract chr
  $ann_chr =~ /\S+\.(\D+\w+\D+)(\d+)$/;
  my ($chr_pre, $chr) = ($1, $2, $3);
  if ( !defined $chr_pre || !defined $chr ){
    next;
  }
  $chr_pre =~ s/[_.]$//;
  # Next: skip genes on scaffolds and other non-chromosome molecules
  unless ( $chr_pre =~ /chloro|chl|CP|mito|MT|ctg|contig|tig|pilon|scaff|sc|super|un\w+\d+/i ){
    $chr_gene_count++;
    $chr =~ s/^0*([^0]+)/$1/;
    $chr_hsh{$chr}++;
  }
}

# Do a second-pass reading of the pan_table and read pan_table into an array.
open ($PAN_FH, "<", $pan_table) or die "Can't open in pan_table: $pan_table\n";
my @pangene_table;
my %HoH_panID_chr;
my %seen_mol;
my %skipped_mols;
while (<$PAN_FH>) {
  chomp;
  next unless (/^\S+/);
  my $line = $_;
  my ($panID, $gene, $ann_chr, $start, $end, $orient) = split(/\t/, $line);

  # From the second field (prefixed genes), extract the annot name
  my $ANN_REX = $annot_regex;
  my $ann = $gene;
  $ann =~ s/$ANN_REX/$1/;

  # From the third field, a annot.chr string, extract chr
  $ann_chr =~ /\S+\.(\D+\w+\D+)(\d+)$/;
  my ($chr_pre, $chr) = ($1, $2, $3);
  if ( !defined $chr_pre || !defined $chr ){
    if ($verbose){ say "For pan-gene consensus, skipping unrecognized annotation-prefix-chr pattern: $ann_chr" }
    next;
  }
  $chr_pre =~ s/[_.]$//;
  # Next: skip genes on scaffolds and other non-chromosome molecules
  if ( $chr_pre =~ /chloro|chl|CP|mito|MT|ctg|contig|tig|pilon|scaff|sc|super|un\w+\d+/i ){
    if ($verbose>1){ say "For pan-gene consensus, skipping non-chromosome gene [$chr_pre $chr]\t$gene" }
  }
  else {
    $chr_gene_count++;
    $chr =~ s/^0*([^0]+)/$1/;
    if ($chr_hsh{$chr} < $chr_gene_count/100){
      if ($seen_mol{"$chr_pre$chr"}){ next }
      else {
        $seen_mol{"$chr_pre$chr"}++;
        $skipped_mols{"$chr_pre$chr"}++; 
      }
    }
    else {
      $chr_hsh{$chr}++;
      #say "A:\t($panID, $gene, $ann, $chr_pre, $chr)";
      $HoH_panID_chr{$panID}{$chr}++;
      my @seven_elts = ( $panID, $gene, $ann, $chr, $start, $end, $orient );
      push( @pangene_table, \@seven_elts );
    }
  }
}
if (%skipped_mols){
  say "Skipped the folowing molecules because of low gene count:";
  for my $short ( keys %skipped_mols ){
    say "  $short: $skipped_mols{$short}";
  }
}

# Sort @pangene_table to determine order per (1) annot; (2) chromosome; (3) posn
my @sorted_table = sort { 
     $a->[2] cmp $b->[2] || $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] 
   } @pangene_table; 

# Calculate gene order per annot-and-chromosome.
my ($prev_ann, $prev_chr) = ("", "");
my $ord=0;
my @elts_with_order;
my %seen_panIDs;
my @pangene_table_ordered;
my %pangene_elts_per_ann;
my $pref_annot_label;
my %annots;
my %seen_chr;
my $prev_panID="";
foreach my $row ( @sorted_table ) {
  my ($panID, $gene, $ann, $chr, $start, $end, $orient) = @$row;
  #say join("\t", "AA: ", $panID, $gene, $ann, $chr, $start, $end, $orient);
  unless ( $seen_chr{$chr} ){ $seen_chr{$chr}++ }
  unless ( $annots{$ann} ){ $annots{$ann}++ }
  unless ( $seen_panIDs{$panID} ){ $seen_panIDs{$panID}++ }
  if ($ann =~ /$PREF_REX/){ $pref_annot_label = $ann }
  if ($panID eq $prev_panID){ # Use just one of a tandem duplicated panID
    next;
  }
  else { # not a tandem duplicate
    if ($ann eq $prev_ann && $chr == $prev_chr){
      $ord++;
      @elts_with_order = ( $panID, $ann, $chr, $ord, $start, $end, $orient );
      push ( @pangene_table_ordered, [@elts_with_order] );
      $pangene_elts_per_ann{$ann}{$panID} = [ $panID, $ann, $chr, $ord, $start, $end, $orient ];
      #say join("\t", @elts_with_order);
      ($prev_ann, $prev_chr) = ($ann, $chr);
    }
    elsif ($ann ne $prev_ann || $chr != $prev_chr){
      $ord=1;
      @elts_with_order = ( $panID, $ann, $chr, $ord, $start, $end, $orient );
      push ( @pangene_table_ordered, [@elts_with_order] );
      $pangene_elts_per_ann{$ann}{$panID} = [ $panID, $ann, $chr, $ord, $start, $end, $orient ];
      #say join("\t", @elts_with_order);
      ($prev_ann, $prev_chr) = ($ann, $chr);
      $ord++;
    }
  }
  $prev_panID = $panID;
}
my $num_chrs = keys %seen_chr;

#say Dumper(\%pangene_elts_per_ann);

# Traverse $pangene_elts_per_ann and report panIDs in order for the preferred annot
#  The data structure is HoHoA:    annot => panID => array
# Sort hash of hash of arrays: first by chromosome, then by position
my $CONS_FH;
if ($consen_out){
  open $CONS_FH, ">", $consen_out or die "Can't open out consen_out file $consen_out:$!\n";
}
my %seen_pref_annot_panID;
foreach my $panID ( sort { ${$pangene_elts_per_ann{$pref_annot_label}}{$a}[2] <=>
                           ${$pangene_elts_per_ann{$pref_annot_label}}{$b}[2] 
                              ||
                           ${$pangene_elts_per_ann{$pref_annot_label}}{$a}[3] <=>
                           ${$pangene_elts_per_ann{$pref_annot_label}}{$b}[3]
                         } keys %{$pangene_elts_per_ann{$pref_annot_label}} ){
  my ($panID, $ann, $chr, $ord, $start, $end, $orient) = @{$pangene_elts_per_ann{$pref_annot_label}{$panID}};
  unless ($seen_pref_annot_panID{$panID}){ $seen_pref_annot_panID{$panID}++ }
  if ($consen_out){ say $CONS_FH join ("\t", $panID, $chr, $ord, $orient); }
  else { say join ("\t", $panID, $chr, $ord, $orient); }
}

# Check and report panIDs missing from the preferred annotation
my $UN_FH;
my %not_seen_panIDs;
if ($unplaced_out){
  open $UN_FH, ">", $unplaced_out or die "Can't open out unplaced_out file $unplaced_out:$!\n";
}
if ($unplaced_out){ # Check and repport only if an output file was indicated.
  foreach my $panID (keys %seen_panIDs){
    unless (exists $seen_pref_annot_panID{$panID}){ say $UN_FH $panID; }
  }
}

__END__


2023
S. Cannon
02-23 Initial version, based on order_gapfill.pl
02-25 Make summary report of %skipped_mols


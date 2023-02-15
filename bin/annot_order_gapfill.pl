#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use feature "say";

my $usage = <<EOS;
Given alignment of gene order with pangene IDs determined by alignment of gene orders 
(by annot_order_encode.pl and annot_order_decode.pl), place leftover pangenes
relative to the pangenes that have established, alignment-based placements.

Usage: annot_order_gapfill.pl -consen_table CONSEN_TABLE  -unplaced UNPLACED_LIST \
                               -pan_table PANGENE_TABLE
                               
  CONSEN_TABLE, e.g. consen_gene_order.tsv, has four columns:
    panID  chr#  order  orientation

  UNPLACED_LIST is a list of panIDs, e.g. consen_pan_unplaced.tsv

  PANGENE_TABLE, e.g. e.g. 22*_posn.hsh.tsv, has six columns:
    panID  geneID  annot.chr  start  end  orientation

Strategy: For each pangene without placement in CONSEN_GENE_ORDER (but with placement 
information in PANGENE_TABLE), score every gene on the chromosome to which it belongs as 
being "before" or "after" the target gene, with respect to each annotation set. This 
will give a hash of genes on a chromosome, with values indicating "beforeness" and 
"afterness". Using that information, place the target gene in the PANGENE_TABLE.

  REQUIRED:
    -consen     CONSEN_GENE_ORDER
    -unplaced   UNPLACED_LIST
    -pan_table  PANGENE_TABLE

  OPTIONS:
    -outfile    File with previously unplaced genes added into the consen_table.
    -annot_regex  Regular expression for capturing annotation name from gene ID, e.g.
                    \"([^.]+\\.[^.]+\\.[^.]+\\.[^.]+)\\..+\"
                      for four dot-separated fields, e.g. vigan.Shumari.gnm1.ann1 (default)
                    or \"(\\D+\\d+\\D+)\\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    -verbose   (boolean) For some debugging info, to STDOUT. Use -v -v to give more output.
    -help      (boolean) This message. 
EOS

my ($consen, $unplaced, $pan_table);
my $help;
my $outfile;
my $verbose=0;
my $logstr="";
my $annot_regex = "([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+"; 

GetOptions (
  "consen=s" => \$consen,
  "unplaced=s" => \$unplaced,
  "pan_table=s" =>  \$pan_table,
  "outfile:s" =>   \$outfile,
  "annot_regex:s" => \$annot_regex,
  "verbose+" =>     \$verbose,
  "help" =>         \$help,
);

die "\n$usage\n" if ( $help || ! $consen || ! $unplaced || ! $pan_table );

# Read list of unplaced panIDs into a hash
open (my $UN_FH, "<", $unplaced) or die "Can't open in unplaced file: $unplaced\n";
my %unplaced;
while (<$UN_FH>){
  chomp;
  my $panID = $_;
  $unplaced{$panID}++;
}

# Read table of panIDs with alignment-based placements
open (my $CONSEN_FH, "<", $consen) or die "Can't open in consensus file: $consen\n";
my %consen_table; # Working version; will be augmented with unplaced panIDs
my %consen_table_entire; # Version with merged "uni-key", to avoid collapse of repeated panIDs
while (<$CONSEN_FH>){
  chomp;
  my ( $panID, $chr, $order, $orient ) = split(/\t/, $_);
  $chr =~ s/chr0*//;
  $consen_table{$chr}{$panID} = [ $panID, $chr, $order, $orient ];
  my $merged_key = join("__", $panID, $chr, $order, $orient );
  $consen_table_entire{$merged_key} = [ $panID, $chr, $order, $orient ];
}


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
      say "Skipping molecule $chr_pre$chr because of low gene count: $chr_hsh{$chr}";
      say " ... suggesting a scaffold or other atypical chromosome.";
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

# Sort @pangene_table to determine order per (1) annot; (2) chromosome; (3) posn
my @sorted_table = sort { 
     $a->[2] cmp $b->[2] || $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] 
   } @pangene_table; 

say "Calculating gene order per annot-and-chromosome";
my ($prAnn, $prChr) = ("", "");
my $ord=0;
my @elts_with_order;
my @pangene_table_ordered;
my %pangene_elts_per_ann;
my %annots;
my %seen_chr;
foreach my $row ( @sorted_table ) {
  my ($panID, $gene, $ann, $chr, $start, $end, $orient) = @$row;
  #say join("\t", "AA: ", $panID, $gene, $ann, $chr, $start, $end, $orient);
  unless ( $seen_chr{$chr} ){ $seen_chr{$chr}++ }
  unless ( $annots{$ann} ){ $annots{$ann}++ }
  if ($ann eq $prAnn && $chr == $prChr){
    $ord++;
    @elts_with_order = ( $panID, $ann, $chr, $ord, $start, $end, $orient );
    push ( @pangene_table_ordered, [@elts_with_order] );
    $pangene_elts_per_ann{$ann}{$panID} = [ $panID, $ann, $chr, $ord, $start, $end, $orient ];
    #say join("\t", @elts_with_order);
    ($prAnn, $prChr) = ($ann, $chr);
  }
  elsif ($ann ne $prAnn || $chr != $prChr){
    $ord=1;
    @elts_with_order = ( $panID, $ann, $chr, $ord, $start, $end, $orient );
    push ( @pangene_table_ordered, [@elts_with_order] );
    $pangene_elts_per_ann{$ann}{$panID} = [ $panID, $ann, $chr, $ord, $start, $end, $orient ];
    #say join("\t", @elts_with_order);
    ($prAnn, $prChr) = ($ann, $chr);
    $ord++;
  }
}
my $num_chrs = keys %seen_chr;

say "Finding the most frequent chromosome for each pan-gene set";
my %top_chr;
my %chr_ct_top_chr;
if ($verbose>1) {print "#pangeneID\tchr:count ...\n"}
foreach my $panID (sort keys %HoH_panID_chr) {
  if ($verbose>1){ $logstr .= "$panID\t" }
  # Sort chromosomes by count of chromosomes seen for this panID
  foreach my $chr ( sort { $HoH_panID_chr{$panID}{$b} <=> $HoH_panID_chr{$panID}{$a} } 
                            keys %{$HoH_panID_chr{$panID}} ) {
    my $chr_ct = $HoH_panID_chr{$panID}{$chr};
    # Store first chromosome (most frequent), and the number of times it was seen.
    unless ($top_chr{$panID}){ 
      $top_chr{$panID} = $chr; 
      $chr_ct_top_chr{$panID} = $chr_ct;
    }
    if ($verbose>1){ $logstr .= "$chr:$chr_ct " }
  }
  if ($verbose>1){ $logstr .= "\n" }
}
if ($verbose>1){ print "$logstr\n" }

if ($verbose>1) {print "#pangeneID\ttop_chr:count\n"}
foreach my $panID (sort keys %top_chr) {
  if ($verbose>1){ printf "%s\t%s:%s\n", $panID, $top_chr{$panID}, $chr_ct_top_chr{$panID} }
}
if ($verbose>1) {print "\n"}

#say Dumper(\@pangene_table_ordered);

# For each unplaced target gene (on a chromosome), score every gene (on that chromosome)
# as being either before the target gene or after it, for each annotation set.
# The verdict for each gene will be stored in $target_gene_scores_HoH{$target_panID}{$panID},
# containing small negative integer values for genes before the target and small postive values after it.
my %target_gene_scores_HoH;
my %orient_panID; # to hold predominant orientation, as + or - integer, keyed on panID
my $count;
say "Scoring each gene relative to unplaced target genes";
foreach my $target_panID (keys %unplaced){
  $count++;
  my $main_chr = $top_chr{$target_panID};
  if ($verbose){
    say "$count\tDominant chr of $target_panID is $main_chr";
  }
  foreach my $ann (keys %annots){
    if ( defined $pangene_elts_per_ann{$ann}{$target_panID}[2] ){
      my $target_panID_chr = $pangene_elts_per_ann{$ann}{$target_panID}[2];
      my $target_panID_order = $pangene_elts_per_ann{$ann}{$target_panID}[3];
      #say "CHECK: $target_panID\t$target_panID_chr\t$target_panID_order\t$ann";
      #say Dumper($pangene_elts_per_ann{$ann}{$target_panID});
      my $compare;
      foreach my $panID ( keys %{$pangene_elts_per_ann{$ann}} ){
        my ($panID, $ann, $chr, $ord, $start, $end, $orient) = @{$pangene_elts_per_ann{$ann}{$panID}};
        if (defined $chr && defined $target_panID_chr && $chr == $target_panID_chr){
          if ($ord < $target_panID_order){ # This panID comes before the target
            $target_gene_scores_HoH{$target_panID}{$panID}--;
          }
          else { # This panID comes at or after the target
            $target_gene_scores_HoH{$target_panID}{$panID}++;
          }
          # Accumulate and store a consensus orientation for this panID# 
          $orient =~ /-/ ?  $orient_panID{$panID}-- : $orient_panID{$panID}++;
          #say "$panID\t$chr\t$ann\t$ord\t$orient\t$target_gene_scores_HoH{$target_panID}{$panID}";
        }
      }
    }
  }
}

#say Dumper($target_gene_scores_HoH{$target_panID});
#say Dumper(%consen_table);
# Now traverse the table of panIDs that have been positioned by alignment, and look up
# whether each ID occurs before or after the target ID, relative to available annotations.
$count = 0;
foreach my $target_panID (keys %unplaced){
  $count++;
  if ($verbose>1){ say "$count\t$target_panID"; }
  my $main_chr = $top_chr{$target_panID};
  my ($before_ID, $after_ID);
  my @fore_aft_score;
  $fore_aft_score[0] = -9;
  my @consen_IDs_ordered;
  my $idx=1;
  #say "MAIN CHR: $main_chr";
  foreach my $panID ( sort { $consen_table{$main_chr}{$a}[2] <=> $consen_table{$main_chr}{$b}[2] }
                      keys %{$consen_table{$main_chr}} ){
    if (defined $target_gene_scores_HoH{$target_panID}{$panID}){
      my $placement_score = $target_gene_scores_HoH{$target_panID}{$panID};
      $fore_aft_score[$idx] = $placement_score;
      $consen_IDs_ordered[$idx] = $panID;
      #say join("\t", @{$consen_table{$main_chr}{$panID}}, $placement_score);

      if ($fore_aft_score[$idx-2] < 0 && $fore_aft_score[$idx-1] < 0 && $fore_aft_score[$idx] >=0){
        my $FORE_ID = $consen_IDs_ordered[$idx-1];
        my $FORE_POS = ${$consen_table{$main_chr}{$FORE_ID}}[2];
        my $HERE = $target_panID;
        my $AFT_ID =  $panID;
        my $AFT_POS = ${$consen_table{$main_chr}{$AFT_ID}}[2];
        my $here_orient;
        $orient_panID{$target_panID} < 0 ? $here_orient = "-" : $here_orient = "+";
        my $new_pos = ($AFT_POS + $FORE_POS)/2;
        
        if ($verbose>1){
          say "FORE: ", join("\t", @{$consen_table{$main_chr}{$FORE_ID}}, $fore_aft_score[$idx-1]);
          say "HERE: ", join("\t", $target_panID, $main_chr, ($AFT_POS + $FORE_POS)/2, $here_orient, 0);
          say "AFT:  ", join("\t", @{$consen_table{$main_chr}{$AFT_ID}}, $fore_aft_score[$idx]);
        }
        $consen_table{$main_chr}{$target_panID} = [ $target_panID, $main_chr, $new_pos, $here_orient ];
        my $merged_key = join("__", $target_panID, $main_chr, $new_pos, $here_orient );
        $consen_table_entire{$merged_key}++;
        last;
      }
      $idx++;
    }
  }
  #say Dumper(%consen_table);
}

# Final traversal of %consen_table to print the results - now with added elements from %unplaced.
# First put it into an array of arrays, before sorting.
# Renumber the panIDs, since the initial numbers may not have sufficient space to permit
# addition of unplaced genes within the available ordered integers, given the placement scheme.
my @consen_table_AoA;
my $chr_count = keys %chr_hsh;
my $OUT_FH;
if ($outfile){
  open ($OUT_FH, ">", $outfile) or die "Can't open outfile: $outfile\n";
}
for my $record (keys %consen_table_entire){
  my @parts = split(/__/, $record);
  my ($panID, $chr, $order, $orient) = @parts;
  my $printchr;
  if ($chr_count <= 9){ # don't zero-pad the chromosome number
    $printchr = sprintf("chr%d", $chr);
  }
  else { # pad chromosome number for two digits
    $printchr = sprintf("chr%02d", $chr);
  }
  my @new_record = ($panID, $printchr, $order, $orient);
  push @consen_table_AoA, \@new_record;
}
my @sorted_consen_table_AoA = sort {
  $a->[1] cmp $b->[1] ||   # chromosome string
  $a->[2] <=> $b->[2]      # position
} @consen_table_AoA;

my $new_ordinals = 0;
%seen_chr = ();
for my $row (@sorted_consen_table_AoA){
  my ($panID, $chr, $orig_order, $orient) = @$row;
  if ($seen_chr{$chr}){
    $new_ordinals += 100;
  }
  else { # not seen_chr
    $seen_chr{$chr}++;
    $new_ordinals = 100;
  }
  
  if ($outfile){
    say $OUT_FH join ("\t", $panID, $chr, $new_ordinals, $orient, $orig_order);
  }
  else {
    say join ("\t", $panID, $chr, $new_ordinals, $orient, $orig_order);
  }
}

__END__
2023
S. Cannon
02-14 Initial version, based on consen_pangene_order.pl


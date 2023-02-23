#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Graph;
use Graph::Directed;
use feature "say";

my $usage = <<EOS;
Given six-column input (pan-IDs, gene_IDs, annot.chr#, start, end, orientation), 
calculate consensus gene order. (That's the objective, anyway).

NOTE: This script is experimental (i.e. it doesn't work with real genomic data!).
The problem is that the pruning redundant panIDs leaves gaps; and witout pruning, there are cycles.


Usage: pan_graph_order.pl -pan_table PANGENE_TABLE.tsv [options]
   PANGENE_TABLE.tsv has six columns:
     panID  geneID  chr  annot.chr  start  end  orientation
     pan00001	phalu.G27455.gnm1.ann1.Pl11G0000359800.1	phalu.G27455.gnm1.Pl11	44272147	44274071	-
     pan00001	phalu.G27455.gnm1.ann1.Pl11G0000365800.1	phalu.G27455.gnm1.Pl11	44964635	44966334	+
     pan00001	phavu.5-593.gnm1.ann1.Pv5-593.11G183000.1	phavu.5-593.gnm1.Chr11	53201894	53202521	+
     pan00001	phavu.5-593.gnm1.ann1.Pv5-593.11G183100.1	phavu.5-593.gnm1.Chr11	53221588	53222946	+

  REQUIRED:
    -pan_table  PANGENE_TABLE

  OPTIONS:
    -outfile  Specify OUT_FH; otherwise, default to STDOUT.
    -prefix   Characters used to prefix the consensus chromosome names; default "chr"
                Examples: Glycine.pan3.chr, Medicago.pan1.chr
    -annot_regex  Regular expression for capturing annotation name from gene ID, e.g.
               \"([^.]+\\.[^.]+\\.[^.]+\\.[^.]+)\\..+\"
                 for four dot-separated fields, e.g. vigan.Shumari.gnm1.ann1 (default)
               or \"(\\D+\\d+\\D+)\\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    -verbose  For some debugging info, to STDOUT. For even more output, call again: -v -v
    -help     This message. 
EOS

my ($outfile, $pan_table, $help);
my $prefix="chr";
my $verbose=0;
my $annot_regex = "([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+";

GetOptions (
  "pan_table=s" => \$pan_table,
  "outfile:s" => \$outfile,
  "prefix:s" =>  \$prefix,
  "verbose+" =>  \$verbose,
  "help" =>      \$help,
);

my $dontuse = <<NOPE;
NOTE: This script is experimental (i.e. it doesn't work with real genomic data!). The problem
is that the pruning redundant panIDs leaves gaps; and witout pruning, there are cycles.
The script is left here in case parts are useful for some other case where a graph is appropriate.
Or maybe as a skeleton, left in warning.
NOPE

die "$dontuse\n";

#die "\n$usage\n" if ( $help or not -t STDIN );
die "\n$usage\n" if ( $help || !$pan_table );

my $logstr;

my $OUT_FH;
if ($outfile) { open ($OUT_FH, ">", $outfile) or die "\nUnable to open output file for writing: $!\n\n"; }

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
        say "Skipping molecule $chr_pre$chr because of low gene count: $chr_hsh{$chr}";
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

# Sort @pangene_table to determine order per (1) annot; (2) chromosome; (3) posn
my @sorted_table = sort {
     $a->[2] cmp $b->[2] || $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]
   } @pangene_table;

say "# Calculate gene order per annot-and-chromosome";
say "# Also filter to the first panID occurrance per chr & annot.\n";
my ($prAnn, $prChr) = ("", "");
my $ord=0;
my @elts_with_order;
my %seen_annot;
my @pangene_table_ordered;
my %pangene_elts_per_ann;
my %seen_panID_per_annot_and_chr;
foreach my $row ( @sorted_table ) {
  my ($panID, $gene, $ann, $chrnum, $start, $end, $orient) = @$row;
  #say join("\t", "AA: ", $panID, $gene, $ann, $chrnum, $start, $end, $orient);
  unless ($seen_annot{$ann}){ $seen_annot{$ann}++ }
  if ($seen_panID_per_annot_and_chr{$ann}{$chrnum}{$panID}){
    #say "Skip $panID from $chrnum, $ann";
    next;
  }
  else {
    $seen_panID_per_annot_and_chr{$ann}{$chrnum}{$panID}++;
    if ($ann eq $prAnn && $chrnum == $prChr){
      $ord++;
      @elts_with_order = ( $panID, $ann, $chrnum, $ord, $start, $end, $orient );
      push ( @pangene_table_ordered, [@elts_with_order] );
      $pangene_elts_per_ann{$ann}{$panID} = [ $panID, $ann, $chrnum, $ord, $start, $end, $orient ];
      #say join("\t", @elts_with_order);
      ($prAnn, $prChr) = ($ann, $chrnum);
    }
    elsif ($ann ne $prAnn || $chrnum != $prChr){
      $ord=1;
      @elts_with_order = ( $panID, $ann, $chrnum, $ord, $start, $end, $orient );
      push ( @pangene_table_ordered, [@elts_with_order] );
      $pangene_elts_per_ann{$ann}{$panID} = [ $panID, $ann, $chrnum, $ord, $start, $end, $orient ];
      #say join("\t", @elts_with_order);
      ($prAnn, $prChr) = ($ann, $chrnum);
      $ord++;
    }
  }
}

my $num_annots = keys %seen_annot;
say "Number of annotations: $num_annots";

say "# Find the most frequent chromosome for each pan-gene set\n";
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

say "# Store order, start, and end for each panID.";
my %orders_per_top_chr_HoA;
my %starts_per_top_chr_HoA;
my %ends_per_top_chr_HoA;
my %orients_per_top_chr_HoA;
my $orient_num;
my %chr_panID_ct; # keys will be the chrIDs
foreach my $row ( @pangene_table_ordered ) {
  my ( $panID, $ann, $chr, $order, $start, $end, $orient ) = @$row;
  #say "BB: $panID, $ann, $chr, $order, $start, $end $orient";
  $chr_panID_ct{$chr}++;
  if ($orient eq "+"){ $orient_num = 1 }
  elsif ($orient = "-"){ $orient_num = -1 }
  else { $orient_num = 0 }
  if ($chr eq $top_chr{$panID}){
    push @{ $orders_per_top_chr_HoA{$panID} }, $order;
    push @{ $starts_per_top_chr_HoA{$panID} }, $start;
    push @{ $ends_per_top_chr_HoA{$panID} }, $end;
    push @{ $orients_per_top_chr_HoA{$panID} }, $orient_num;
  }
}
my $num_chrs = keys %chr_panID_ct;

#say Dumper(\%chr_panID_ct);
#say Dumper(\%orients_per_top_chr_HoA);

my %consen_pan_coords;
my $most_freq_chr;
for my $panID ( keys %starts_per_top_chr_HoA ) {
  my $median_order = sprintf("%i", calc_median(\@{ $orders_per_top_chr_HoA{$panID} }));
  my $median_start = sprintf("%i", calc_median(\@{ $starts_per_top_chr_HoA{$panID} }));
  my $median_end = sprintf("%i", calc_median(\@{ $ends_per_top_chr_HoA{$panID} }));
  my $median_orient = sprintf("%i", calc_median(\@{ $orients_per_top_chr_HoA{$panID} }));
  $most_freq_chr = $top_chr{$panID};
  my @pan_coords = ( $most_freq_chr, $median_order, $median_start, $median_end, $median_orient );
  #say "($panID, $most_freq_chr, $median_order, $median_start, $median_end, $median_orient)";
  $consen_pan_coords{$panID} = [ @pan_coords ];
}

say "# Build array of consensus pangene positions, sorted by";
say "#   (1) chromosome, (2) median order, (3) genomic position.";
my @sorted_consen_pan_coords;
for my $panID (sort { 
    $consen_pan_coords{$a}[0] cmp $consen_pan_coords{$b}[0] or 
    $consen_pan_coords{$a}[1] <=> $consen_pan_coords{$b}[1] or 
    $consen_pan_coords{$a}[2] <=> $consen_pan_coords{$b}[2] 
  } keys %consen_pan_coords ) {
  my ($chrnum, $posn, $start, $end) = @{$consen_pan_coords{$panID}};
  #say "($panID => $chrnum, $posn, $start, $end)";
  push @sorted_consen_pan_coords, [$panID, @{$consen_pan_coords{$panID}}];
}

#say Dumper(\@sorted_consen_pan_coords);

say "# Calculate initial pangene order per chromosome, based on ";
say "#  (1) chromosome, (2) median order, (3) genomic position.";
say "# The main purpose is to find a first and a last gene on each chromosome, for use";
say "# in deriving a graph-based path between the two.\n";
my %seen_chr;
my %first_panIDs_per_chr;
my %last_panIDs_per_chr;
my $ordering=1;
my @prev_row;
my $idx = 0;
foreach my $row ( @sorted_consen_pan_coords ) {
  my ($panID, $chrnum, $ord, $start, $end, $orient_num) = @$row;
  #say "($panID, $chrnum, $ord, $start, $end, $orient_num)";
  if ( $seen_chr{$chrnum} ){ 
    $seen_chr{$chrnum}++; 
    @prev_row = @$row;
    $idx++;
    next;
  }
  else { # haven't seen this chrnum yet
    $seen_chr{$chrnum}++;
    if ($idx > 0){ # If $idx == 0, there is no previous record
      my ($panID, $chrnum, $ord, $start, $end, $orient_num) = @prev_row;
      if ($verbose){say join("\t", $idx-1, @prev_row), "\n";}
      $last_panIDs_per_chr{$chrnum} = $panID;
    }
    if ($verbose){say join("\t", $idx, @$row);}
    $first_panIDs_per_chr{$chrnum} = $panID;
    $idx++;
  }
}
{# Handle the last element of @sorted_consen_pan_coords
  my ($panID, $chrnum, $ord, $start, $end, $orient_num) = @{$sorted_consen_pan_coords[$idx-1]};
  $last_panIDs_per_chr{$chrnum} = $panID;
  if ($verbose){say join ("\t", $idx-1, @{$sorted_consen_pan_coords[$idx-1]} );}
}

#say Dumper(\%first_panIDs_per_chr);
#say Dumper(\%last_panIDs_per_chr);


say "# Traverse \@pangene_table_ordered again and store provisional panID pairs";
say "# and weights,  based on panID adjacencies in each annotation set.\n";
# Below, in preparing edges (pairs) for the graph, a high weight is "costly" 
# in the graph and low is "inexpensive" or "low-friction".
# The weight ($wt = $seen_panID_pair{"$prev_panID,$panID"}) is started high, at $num_annots,
# and lowers (toward 0) as the panID pair is seen in other annotations.weighted 
%seen_panID_per_annot_and_chr = ();
my %seen_panID_pair;
my %panID_pairs_per_chr;
my ($prev_panID, $panID) = ("", "");
my ( $ann, $chrnum, $order, $start, $end, $orient );
my $prev_chr;
foreach my $row ( @pangene_table_ordered ) {
  ( $panID, $ann, $chrnum, $order, $start, $end, $orient ) = @$row;
  my $wt;

  if ($seen_panID_per_annot_and_chr{$chrnum}{$ann}{$panID}){
    $seen_panID_per_annot_and_chr{$chrnum}{$ann}{$panID}++;
    next; # Store only one instance of a panID in a chromosome, to avoid short-circuits in the graphs
  }
  else {
    $seen_panID_per_annot_and_chr{$chrnum}{$ann}{$panID}++;
    # Also Remember that we (plan to) have seen $last_panIDs_per_chr{$prev_chr}
    if ($chrnum > 1){ # Set weight to 1 (low/good), because we want this path traversed.
      $seen_panID_per_annot_and_chr{$prev_chr}{$ann}{$last_panIDs_per_chr{$prev_chr}} = 1; 
    }
    #say join("\t", $chrnum, $panID, $ann, $seen_panID_per_annot_and_chr{$chrnum}{$ann}{$panID});
  
    #say "CC: $panID, $ann, $chrnum, $order, $start, $end $orient";
    if ($seen_chr{$chrnum}){
      if ($seen_panID_pair{"$prev_panID,$panID"}){ 
        $seen_panID_pair{"$prev_panID,$panID"}--; 
      } else {
        $seen_panID_pair{"$prev_panID,$panID"} = $num_annots; # Set high initially
      }
      $wt = $seen_panID_pair{"$prev_panID,$panID"};
      #say join("\t", ("Y:", $chrnum, $prev_panID, $panID, $wt));
      #$g1->add_weighted_edge($prev_panID,$panID, $wt);
      $panID_pairs_per_chr{$chrnum}{"$prev_panID,$panID"} = $wt;
      # Remember that we have seen $first_panIDs_per_chr{$chrnum}
    }
    else { # not $seen_chr{$chrnum}; 
      $wt = 1; # A low ("good") value - as opposed to the high of $num_annots
      if ($chrnum > 1){ # Handle the last panID of the previous chromosome
        # Set last panID to the one determined to be the last globally.
        $panID = $last_panIDs_per_chr{$prev_chr};
        #say join("\t", ("Z:", $chrnum-1, $prev_panID, $panID, $wt)), "\n";
        $panID_pairs_per_chr{$chrnum-1}{"$prev_panID,$panID"} = $wt;
      }
      # Set the first panID to the one determined to be first globally.
      $prev_panID = $first_panIDs_per_chr{$chrnum};
      # Remember that we have seen $first_panIDs_per_chr{$chrnum}
      $seen_panID_per_annot_and_chr{$chrnum}{$ann}{$first_panIDs_per_chr{$chrnum}}++;
  
      #say join("\t", ("X:", $chrnum, $prev_panID, $panID, $wt));
      $panID_pairs_per_chr{$chrnum}{"$prev_panID,$panID"} = $wt;
      $seen_chr{$chrnum}++;
    }
  }
  $prev_panID = $panID;
  $prev_chr = $chrnum;
}
{ # Handle the last panID of the previous chromosome
  # Set last panID to the one determined to be the last globally.
  my $panID = $last_panIDs_per_chr{$prev_chr};
  my $wt = 1; # A low value - as opposed to the high of $num_annots
  #say join("\t", ("Z:", $chrnum, $prev_panID, $panID, $wt)), "\n";
  $panID_pairs_per_chr{$chrnum}{"$prev_panID,$panID"} = $wt;
}

# Report panID pairs, filtering out panID pairs in which the left (FROM) or right (TO)
# of the pair has been seen on this chr, in a lower-weighted path.
my %seen_FROM_panID_per_chr;
my %seen_TO_panID_per_chr;
my %filtered_panID_pairs_per_chr;
foreach my $chrnum ( sort { $a <=> $b } keys %panID_pairs_per_chr ){
  foreach my $pair ( sort { $panID_pairs_per_chr{$chrnum}{$a} 
                        <=> $panID_pairs_per_chr{$chrnum}{$b} } 
                     keys %{$panID_pairs_per_chr{$chrnum}}){
    my $wt = $panID_pairs_per_chr{$chrnum}{$pair};
    my ($first, $second) = split(/,/, $pair);
    #say "($first, $second), $wt";

    #if ($seen_FROM_panID_per_chr{$chrnum}{$first} ){ 
    #  #say "Skip first for $chrnum: ($first, $second), $wt";
    #  next; 
    #}
    #else { $seen_FROM_panID_per_chr{$chrnum}{$first}++; }

    #if ($seen_TO_panID_per_chr{$chrnum}{$second} ){ 
    #  #say "Skip second for $chrnum: ($first, $second), $wt";
    #  next; 
    #}
    #else { $seen_TO_panID_per_chr{$chrnum}{$second}++; }
    
    #say "Keep for $chrnum: ($first, $second), $wt";
    $filtered_panID_pairs_per_chr{$chrnum}{"$first,$second"} = $wt;
    #say join("\t", $chrnum, $first, $second, $wt);
  }
}

# Populate a graph for each chromosome, based on lowest-weight panID pairs
# from the annotation sets. Due to the pruning to eliminate duplicates, the
# chromosome graphs may be fragmented.
my %graph_per_chr;
foreach my $chrnum ( sort { $a <=> $b } keys %filtered_panID_pairs_per_chr ){
  $graph_per_chr{$chrnum} = Graph::Directed->new;
  foreach my $pair ( sort { $filtered_panID_pairs_per_chr{$chrnum}{$a} 
                        <=> $filtered_panID_pairs_per_chr{$chrnum}{$b} } 
                     keys %{$filtered_panID_pairs_per_chr{$chrnum}}){
    my $wt = $filtered_panID_pairs_per_chr{$chrnum}{$pair};
    my ($first, $second) = split(/,/, $pair);

    # Add edge to the graph for this chr
    $graph_per_chr{$chrnum}->add_weighted_edge($first,$second, $wt);
  }
}

# Traverse the chromosome graphs; and also report the count of unique panIDs per chr
say "# chr\tpan_count\tstart\tend";
foreach $chrnum (sort {$a <=> $b} keys %filtered_panID_pairs_per_chr){
  my $count_panIDs_per_chr = keys %{$filtered_panID_pairs_per_chr{$chrnum}};
  my $start = $first_panIDs_per_chr{$chrnum};
  my $end = $last_panIDs_per_chr{$chrnum};
  say join("\t", ("$chrnum:", $count_panIDs_per_chr, $start, $end));
  my @successors = $graph_per_chr{$chrnum}->all_successors($start);
  say "  Number of successors for $start: ", scalar(@successors);
  #say join(", ", @successors), "\n";
  #my @path = $graph_per_chr{$chrnum}->SP_Dijkstra($start, $end);
  #say "Shortest path has ", scalar(@path), " elements, vs. ", $count_panIDs_per_chr, " expected.";
  #say "SP_Dijkstra($start, $end): \t@path\n";
}

## Print all pairs in each chromosome graph
#foreach $chrnum (sort {$a <=> $b} keys %filtered_panID_pairs_per_chr){
#  my $start = $first_panIDs_per_chr{$chrnum};
#  my $end = $last_panIDs_per_chr{$chrnum};
#
#  #say "The graph for chr $chrnum is: ";
#  #say $graph_per_chr{$chrnum}, "\n";
#  
#  #my $sptg = $graph_per_chr{$chrnum}->SPT_Dijkstra($start);
#  #say "The SSSP for chr $chrnum is: ";
#  #say $sptg;
#}


##################################################
# SUBRUTINES

# Return median value for an array of numbers
# See http://stackoverflow.com/questions/5119034/using-perl-to-find-median-mode-standard-deviation
sub calc_median {
  my $value_ref = shift;
  my @values = @$value_ref;
  my $median;
  my $annot = int ((scalar @values)/2);
  my @sorted_values = sort {$a <=> $b} @values;
  if (@values % 2) {
    $median = $sorted_values[ $annot ];
  } else {
    $median = ($sorted_values[$annot-1] + $sorted_values[$annot])/2;
  } 
  return $median;
}

# Print to outfile or to stdout
sub printstr {
  my $str_to_print = join("", @_);
  if ($outfile) {
    print $OUT_FH $str_to_print;
  }
  else {
    print $str_to_print;
  }
}

__END__
2023
02-20 New script "pan_graph_order.pl", derived from order_by_consensus.pl


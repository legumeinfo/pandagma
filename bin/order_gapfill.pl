#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Parallel::ForkManager;
use Symbol;
use IO::Handle;
use Scalar::Util qw< openhandle >;
use List::Util qw(sum);
use feature "say";

$! = 1;  # disable buffering, to get finer-grained real-time process & debugging output

my $usage = <<EOS;
Given alignment of gene order with pangene IDs determined previously (by alignment, with
order_encode.pl and order_decode.pl or by reference, with order_by_reference.pl), place the
leftover pangenes relative to the pangenes that have established, alignment-based placements.

Usage: order_gapfill.pl -consen_table CONSEN_TABLE  -unplaced UNPLACED_LIST \
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
    -nproc      Maximum number of processes to be created. 0 means no forking. (default 0).
    -retain     Retain tempfiles from search of unplaced genes against all others (default false)
    -annot_regex  Regular expression for capturing annotation name from gene ID, e.g.
                    \"([^.]+\\.[^.]+\\.[^.]+\\.[^.]+)\\..+\"
                      for four dot-separated fields, e.g. vigan.Shumari.gnm1.ann1 (default)
                    or \"(\\D+\\d+\\D+)\\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    -verbose   (boolean) For some debugging info, to STDOUT. Use -v -v to give more output.
    -help      (boolean) This message. 
EOS

my ($consen, $unplaced, $pan_table);
my $help;
my $nproc=0;
my $outfile;
my $verbose=0;
my $retain;
my $logstr="";
my $annot_regex = "([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+"; 

GetOptions (
  "consen=s" =>      \$consen,
  "unplaced=s" =>    \$unplaced,
  "pan_table=s" =>   \$pan_table,
  "outfile:s" =>     \$outfile,
  "annot_regex:s" => \$annot_regex,
  "nproc:i" =>       \$nproc,
  "retain" =>        \$retain,
  "verbose+" =>      \$verbose,
  "help" =>          \$help,
);

die "\n$usage\n" if ( $help || ! $consen || ! $unplaced || ! $pan_table );

# Read list of unplaced panIDs into a hash
open (my $UN_FH, "<", $unplaced) or die "Can't open in unplaced file: $unplaced\n";
my %unplaced;
my $ct_unplaced;
while (<$UN_FH>){
  chomp;
  my $panID = $_;
  $unplaced{$panID}++;
  $ct_unplaced++;
}

# Read table of panIDs with alignment- or reference-based placements
open (my $CONSEN_FH, "<", $consen) or die "Can't open in consensus file: $consen\n";
my %consen_table; # Working version; will be augmented with unplaced panIDs
my %consen_table_entire; # Version with merged "uni-key", to avoid collapse of repeated panIDs
my %orient_known_panIDs_by_chr;
my %posn_known_panIDs_by_chr;
while (<$CONSEN_FH>){
  chomp;
  my ( $panID, $chr, $order, $orient ) = split(/\t/, $_);
  $order = 1*$order;
  $chr =~ s/chr0*//;
  $consen_table{$chr}{$panID} = [ $panID, $chr, $order, $orient ];
  my $merged_key = join("__", $panID, $chr, $order, $orient );
  $consen_table_entire{$chr}{$merged_key} = [ $panID, $chr, $order, $orient ];
  $orient_known_panIDs_by_chr{$chr}{$panID} = $orient; # Assumes one panID per chr, (may not be true)
  $posn_known_panIDs_by_chr{$chr}{$panID} = $order; # Assumes one panID per chr, (may not be true)
}

# Do a first-pass reading of the pan_table to get counts per molecule, to help later
# bypass probable scaffolds.
# Read pan_table into an array, and get count of each molecule (chromosome or scaffold.)
open (my $PAN_FH, "<", $pan_table) or die "Can't open in pan_table: $pan_table\n";
my $chr_gene_count=0;
my %panIDs;
my %chr_hsh;
while (<$PAN_FH>) {
  chomp;
  next unless (/^\S+/);
  my $line = $_;
  my ($panID, $gene, $ann_chr, $start, $end, $orient) = split(/\t/, $line);

  unless ($panIDs{$panID}){$panIDs{$panID}++};

  # From the third field, a annot.chr string, extract chr
  $ann_chr =~ /\S+\.(\D+\w+\D+)(\d+)$/;
  my ($chr_pre, $chr) = ($1, $2, $3);
  if ( !defined $chr_pre || !defined $chr ){
    next;
  }
  $chr_pre =~ s/[_.]$//;
  # Next: skip genes on scaffolds and other non-chromosome molecules
  unless ( $chr_pre =~ /chloro|chl|^CP|mito|ctg|contig|tig|pilon|^scaff|^sc|^super|^un\w+\d+/i ){
    $chr_gene_count++;
    $chr =~ s/^0*([^0]+)/$1/;
    $chr_hsh{$chr}++;
  }
}

$annot_regex =~ s/['"]//g;
my $ANN_REX = qr/$annot_regex/;
say "# Annotation regex is: $ANN_REX";

my $ct_panIDs = keys %panIDs;
if ($verbose){ say "unplaced: $ct_unplaced; panIDs: $ct_panIDs; ratio: ", $ct_unplaced/$ct_panIDs; }
if ($ct_unplaced/$ct_panIDs > 0.5){ 
  warn "\nWARNING: The number of unplaced genes is high, at $ct_unplaced vs. $ct_panIDs total panIDs.";
  warn "This suggests a problem in the previous step(s). Check that -annot_regex is set appropriately,";
  warn "so that the second column of the PANGENE_TABLE can be parsed into annotation names.";
  warn "  annot_regex =  $ANN_REX";
  warn "Other debugging steps: turn on printing at commented \"say\" statements containing AA: BB: CC: etc.\n";
  die;
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
  #say $line;

  # From the second field (prefixed genes), extract the annot name
  my $ann = $gene;
  $ann =~ s/$ANN_REX/$1/;

  # From the third field, a annot.chr string, extract chr
  $ann_chr =~ /^\S+\.(\D+\w+\D+)(\d+)$/;
  my ($chr_pre, $chr) = ($1, $2, $3);
  if ( !defined $chr_pre || !defined $chr ){
    if ($verbose>2){ say "For pan-gene consensus, skipping unrecognized annotation-prefix-chr pattern: $ann_chr" }
    next;
  }
  $chr_pre =~ s/[_.]$//;
  # Next: skip genes on scaffolds and other non-chromosome molecules
  if ( $chr_pre =~ /chloro|chl|^CP|mito|ctg|contig|tig|pilon|^scaff|^sc|^super|^un\w+\d+/i ){
    if ($verbose>2){ say "For pan-gene consensus, skipping non-chromosome gene [$chr_pre $chr]\t$gene" }
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
      #say "AA:\t($panID, $gene, $ann, $chr_pre, $chr)";
      #say "HoH_panID_chr: $panID, $gene, $ann, $chr, $start, $end, $orient";
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

# Sort @pangene_table to determine order per annot [2]; chromosome [3];  posn [4]
my @sorted_table = sort { 
     $a->[2] cmp $b->[2] || $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] 
   } @pangene_table; 

say "# Calculating gene order per annot-and-chromosome";
my ($prAnn, $prChr) = ("", 0);
my $ord=0;
my @elts_with_order;
my @pangene_table_ordered;
my %pangene_elts_per_ann;
my %annots;
my %seen_chr;
foreach my $row ( @sorted_table ) {
  my ($panID, $gene, $ann, $chr, $start, $end, $orient) = @$row;
  #say join("\t", "BB: ", $panID, $gene, $ann, $chr, $start, $end, $orient);
  unless ( $seen_chr{$chr} ){ $seen_chr{$chr}++ }
  unless ( $annots{$ann} ){ $annots{$ann}++ }
  if ($ann eq $prAnn && $chr == $prChr){
    $ord++;
    @elts_with_order = ( $panID, $ann, $chr, $ord, $start, $end, $orient );
    push ( @pangene_table_ordered, [@elts_with_order] );
    $pangene_elts_per_ann{$ann}{$panID} = [ $panID, $ann, $chr, $ord, $start, $end, $orient ];
    ##say join("\t", @elts_with_order);
    ($prAnn, $prChr) = ($ann, $chr);
  }
  elsif ($ann ne $prAnn || $chr != $prChr){
    $ord=1;
    @elts_with_order = ( $panID, $ann, $chr, $ord, $start, $end, $orient );
    push ( @pangene_table_ordered, [@elts_with_order] );
    $pangene_elts_per_ann{$ann}{$panID} = [ $panID, $ann, $chr, $ord, $start, $end, $orient ];
    ##say join("\t", @elts_with_order);
    ($prAnn, $prChr) = ($ann, $chr);
    $ord++;
  }
}
my $num_chrs = keys %seen_chr;
my $ct_annots = keys %annots;
say "#  $num_chrs chrs and $ct_annots annotations";

if ($ct_annots > 100){ 
  warn "\nWARNING: The number of annotations is high, at $ct_annots. Check that -annot_regex is set,";
  warn "appropriately, so that the second column of the PANGENE_TABLE can be parsed into annotation names.";
  warn "  annot_regex =  $ANN_REX";
  warn "Other debugging steps: turn on printing at commented \"say\" statements containing AA: BB: CC: etc.\n";
  die;
}

say "# Finding the most frequent chromosome for each panID";
my %top_chr;
my %chr_ct_top_chr;
if ($verbose>2) {print "#pangeneID\tchr:count ...\n"}
foreach my $panID (sort keys %HoH_panID_chr) {
  if ($verbose>2){ $logstr .= "$panID\t" }
  # Sort chromosomes by count of chromosomes seen for this panID
  foreach my $chr ( sort { $HoH_panID_chr{$panID}{$b} <=> $HoH_panID_chr{$panID}{$a} } 
                            keys %{$HoH_panID_chr{$panID}} ) {
    my $chr_ct = $HoH_panID_chr{$panID}{$chr};
    # Store first chromosome (most frequent), and the number of times it was seen.
    unless ($top_chr{$panID}){ 
      $top_chr{$panID} = $chr; 
      $chr_ct_top_chr{$panID} = $chr_ct;
    }
    if ($verbose>2){ $logstr .= "$chr:$chr_ct " }
  }
  if ($verbose>2){ $logstr .= "\n" }
}
if ($verbose>2){ print "$logstr\n" }

if ($verbose>2) {print "#pangeneID\ttop_chr:count\n"}
foreach my $panID (sort keys %top_chr) {
  if ($verbose>2){ printf "%s\t%s:%s\n", $panID, $top_chr{$panID}, $chr_ct_top_chr{$panID} }
}
if ($verbose>2) {print "\n"}


##########
# For each unplaced target gene (on a chromosome), score every gene (on that chromosome)
# as being either before the target gene or after it, for each annotation set.
# The verdict for each gene will be stored in $target_gene_scores_by_annot{$ann}{$chr}{$target_panID}{$panID},
# containing small negative integer values for genes before the target and small postive values after it.

my $tmpdir = "consen_tmp";
unless ( -d $tmpdir ){ mkdir $tmpdir or die "Can't mkdir $tmpdir $!\n" };
say "# Scoring each gene relative to $ct_unplaced unplaced target panIDs.";
say "# This step is done in parallel, writing to temp files at $tmpdir/";

my $pm1 = Parallel::ForkManager->new($nproc);
DATA_LOOP:
foreach my $ann (keys %annots){
  my $pid = $pm1->start and next DATA_LOOP;
  say "# Starting $ann";
  my $ann_file = "$tmpdir/target_gene_scores_$ann";
  my $SYM = gensym; # method from Symbol;  creates an anonymous glob and returns a reference to it.
  open $SYM, ">", $ann_file or die "Can't open out $ann_file: $!\n";
  foreach my $target_panID (keys %unplaced){
    #print ".";
    select()->flush();
    my $main_chr = $top_chr{$target_panID};
    #say "$count\tDominant chr of $target_panID is $main_chr\n  ";
    #say "ann: $ann; target_panID: $target_panID";
    if ( defined $pangene_elts_per_ann{$ann}{$target_panID}[2] ){
      my $target_panID_chr = $pangene_elts_per_ann{$ann}{$target_panID}[2];
      my $target_panID_order = $pangene_elts_per_ann{$ann}{$target_panID}[3];
      #say "CHECK: $target_panID\t$target_panID_chr\t$target_panID_order\t$ann";
      foreach my $panID ( keys %{$pangene_elts_per_ann{$ann}} ){
        my ($panID, $ann, $chr, $ord, $start, $end, $orient) = @{$pangene_elts_per_ann{$ann}{$panID}};
        #if (defined ($ann) ){ say "CC: ($panID, $ann, $chr, $ord, $start, $end, $orient)" }
        if (defined $chr && defined $target_panID_chr && $chr == $target_panID_chr){
          if ($ord < $target_panID_order){ # This panID comes before the target
            say $SYM join ("\t", $chr, $target_panID, $panID, -1);
          }
          else { # This panID comes at or after the target
            say $SYM join ("\t", $chr, $target_panID, $panID, 1);
          }
        }
      }
    }
  }
  say "\n#  Finished $ann\n";
  $pm1->finish; # No structures to return; stored data written to files.
  if (openhandle($SYM)){close $SYM};
}
$pm1->wait_all_children;
say "\n# Finished scoring genes relative to unplaced target genes.";


say "# Combine the scores in the target_gene_scores_by_annot hashes, combining across annotations,";
say "# and record which panIDs were used in merged_target_gene_scores.";
my $chr_count = keys %chr_hsh;
my %merged_target_gene_scores;
my %used_target_panIDs;
foreach my $ann ( keys %annots ){
  if ($verbose){ say "  Processing tempfile for annotation $ann" }
  my $ann_file = "$tmpdir/target_gene_scores_$ann";
  my $SYM = gensym; # method from Symbol;  creates an anonymous glob and returns a reference to it.
  open $SYM, "<", $ann_file or die "Can't open in $ann_file: $!\n";
  while (<$SYM>){
    chomp;
    my ($chr, $target_panID, $panID, $orient) = split("\t", $_);
    $merged_target_gene_scores{$chr}{$target_panID}{$panID} += $orient;
    $used_target_panIDs{$chr}{$panID} = $target_panID unless $used_target_panIDs{$chr}{$panID};
  }
  close $SYM;
  # Clean up large temp files
  unless ($retain){
    if ($verbose){ say "  Removing tempfile |$SYM| $ann_file" }
    unlink $ann_file or die "Can't unlink file $ann_file: $!\n";
  }
}

#say Dumper(\%merged_target_gene_scores);

# Calculate consensus panID orientations (in %orient_target_panID), 
my %orient_target_panID;
foreach my $ann (keys %annots){
  #say "ANNOT: $ann";
  foreach my $target_panID (keys %unplaced){
    my $main_chr = $top_chr{$target_panID};
    #print "Dominant chr of $target_panID is $main_chr\n  ";
    if ( defined $pangene_elts_per_ann{$ann}{$target_panID}[2] ){
      my $target_panID_chr = $pangene_elts_per_ann{$ann}{$target_panID}[2];
      my ($panID, $ann, $chr, $ord, $start, $end, $orient) = @{$pangene_elts_per_ann{$ann}{$target_panID}};
      if (defined $chr && defined $target_panID_chr && $chr == $target_panID_chr){
        # Accumulate and store a consensus orientation for this panID
        if ( $orient =~ /-/ ){ $orient_target_panID{$panID} = "-" }
        else { $orient_target_panID{$panID} = "+" }
        #say "DD: ($panID, $ann, $chr, $ord, $start, $end, $orient), $target_panID_chr, $orient_target_panID{$panID}";
      }
    }
  }
}

# Make revised hashes to hold panID info for panIDs included ("used") in %merged_target_gene_scores
my %consen_table_entire_used; 
my %orient_known_panIDs_by_chr_used;
my %posn_known_panIDs_by_chr_used;
foreach my $chr ( keys %consen_table_entire ){
  foreach my $merged_key ( keys %{$consen_table_entire{$chr}} ){
    my ( $panID, $chr, $order, $orient ) = @{$consen_table_entire{$chr}{$merged_key}};
    foreach my $target_panID ( $used_target_panIDs{$chr}{$panID} ){
      if (defined $target_panID && defined $merged_target_gene_scores{$chr}{$target_panID}{$panID} ){
        #say "target_panID: $target_panID";
        #say "EE: $chr, $target_panID, $panID, $merged_target_gene_scores{$chr}{$target_panID}{$panID}";
        my $merged_key = join("__", $panID, $chr, $order, $orient );
        $consen_table_entire_used{$chr}{$merged_key} = [ $panID, $chr, $order, $orient ];
        $orient_known_panIDs_by_chr_used{$chr}{$panID} = $orient;
        $posn_known_panIDs_by_chr_used{$chr}{$panID} = $order;
      }
    }
  }
}

# Now traverse the table of panIDs that have been positioned by alignment, and look up
# whether each ID occurs before or after the target ID, relative to available annotations.
my $count = 0;
my %consen_IDs_ordered_per_chr;
my %signs_by_chr;
my %missed_panIDs_by_chr;
foreach my $target_panID (keys %unplaced){
  $count++;
  if ($verbose >2){ say "$count\t$target_panID"; }
  my $main_chr = $top_chr{$target_panID};
  my @fore_aft_score;
  $fore_aft_score[0] = -9;
  my $idx=0;
  if ($verbose >1){ say "MAIN CHR: $main_chr; panID: $target_panID" }
  
  # Sort consen_table by chromosome, then position.
  my $ct_missed;
  my $prev_orient = -1;
  foreach my $merged_key ( sort { $consen_table_entire_used{$main_chr}{$a}[2] <=> 
                                  $consen_table_entire_used{$main_chr}{$b}[2] }
                      keys %{$consen_table_entire_used{$main_chr}} ){
    my ($panID, $main_chr, $order, $orient) = split(/__/, $merged_key);
    #say "Merged key has elements ($panID, $main_chr, $order, $orient)";
    # Compare the target panID against all others in this chr, to get a composite before/after placement score.
    if (defined $merged_target_gene_scores{$main_chr}{$target_panID}{$panID}){
      # fore_aft_score holds the count of target_panIDs that occur before (-) or after (+) the present panID
      $fore_aft_score[$idx] = $merged_target_gene_scores{$main_chr}{$target_panID}{$panID};
      unless (exists $consen_IDs_ordered_per_chr{$main_chr}[$idx]){
        $consen_IDs_ordered_per_chr{$main_chr}[$idx] = $panID;
      }

      if ( $fore_aft_score[$idx] < 0 ){ $signs_by_chr{$main_chr}{$target_panID}{$idx} = -1 } 
      elsif ( $fore_aft_score[$idx] == 0 ){ $signs_by_chr{$main_chr}{$target_panID}{$idx} = 0 }
      else { $signs_by_chr{$main_chr}{$target_panID}{$idx} = 1 }
      $prev_orient = $signs_by_chr{$main_chr}{$target_panID}{$idx};
      $idx++;
    }
    else {
      $ct_missed++;
      if ($verbose>2){ say "$ct_missed -- No target gene score for $target_panID vs. $panID"; }
      $signs_by_chr{$main_chr}{$target_panID}{$idx} = $prev_orient;
      $idx++;  # Increment the index, to keep in sync with positions from sorted merged_target_gene_scores
      $missed_panIDs_by_chr{$main_chr}{$merged_key} = [ $panID, $main_chr, $order, $orient ];
    }
    #say join( "\t", "FF:", @{$consen_table_entire_used{$main_chr}{$merged_key}} );
  }
  
  # Determine consensus placement of this target panID
  if ($verbose){say "Do find_placement of target_panID $target_panID on chr $main_chr";}
  my $transition_idx = find_placement( $main_chr, $target_panID, \%signs_by_chr );

  if (defined $transition_idx){
    my $AFT_ID =  $consen_IDs_ordered_per_chr{$main_chr}[$transition_idx];
    my $AFT_ORIENT = $orient_known_panIDs_by_chr_used{$main_chr}{$AFT_ID}; 
    my $AFT_POS = $posn_known_panIDs_by_chr_used{$main_chr}{$AFT_ID};
    my $MERGED_AFT_ID =  join ("__", $AFT_ID, $main_chr, $AFT_POS, $AFT_ORIENT);
  
    my ($FORE_ID, $FORE_ORIENT, $FORE_POS, $MERGED_FORE_ID) = ("", "", "", "");
    if (exists $consen_IDs_ordered_per_chr{$main_chr}[$transition_idx-1]){
      $FORE_ID = $consen_IDs_ordered_per_chr{$main_chr}[$transition_idx-1];
      $FORE_ORIENT = $orient_known_panIDs_by_chr_used{$main_chr}{$FORE_ID};
      $FORE_POS = $posn_known_panIDs_by_chr_used{$main_chr}{$FORE_ID};
      $MERGED_FORE_ID = join ("__", $FORE_ID, $main_chr, $FORE_POS, $FORE_ORIENT);
  
      my $NEW_POS = ($AFT_POS + $FORE_POS)/2;
    
      my $HERE_ID = $target_panID;
      my $HERE_ORIENT = $orient_target_panID{$HERE_ID};
      my $MERGED_HERE_ID = join ("__", $target_panID, $main_chr, $NEW_POS, $HERE_ORIENT); 
    
      # Update the data structures that are indexed by position, to reflect the added panID
      splice @{$consen_IDs_ordered_per_chr{$main_chr}}, $transition_idx, 0, $target_panID;
      $signs_by_chr{$main_chr}{$target_panID}{$NEW_POS} = 0;
      $orient_known_panIDs_by_chr_used{$main_chr}{$target_panID} = $HERE_ORIENT;
      $posn_known_panIDs_by_chr_used{$main_chr}{$target_panID} = $NEW_POS;
    
      my $merged_key = join("__", $target_panID, $main_chr, $NEW_POS, $HERE_ORIENT );
      $consen_table_entire_used{$main_chr}{$merged_key} = [ $target_panID, $main_chr, $NEW_POS, $HERE_ORIENT ];
    
      if ($verbose>1){
        say join("\t", "FORE:", $FORE_ID, $main_chr, $FORE_POS, $FORE_ORIENT, $FORE_POS );
        say join("\t", "HERE:", $target_panID, $main_chr, $NEW_POS, $HERE_ORIENT, $NEW_POS );
        say join("\t", "AFT:", $AFT_ID, $main_chr, $AFT_POS, $AFT_ORIENT, $AFT_POS );
        say "Size of consen_IDs_ordered on chr $main_chr: ", scalar(@{$consen_IDs_ordered_per_chr{$main_chr}});
        say "";
      }
    }
    else {
      say "Transition index failure for ", $transition_idx-1;
    }
  }
  else {
    say "Transition index (panID placement) undefined for $target_panID";
  }

  if ($verbose>1){ say "=================\n"; }
}

# Merge %consen_table_entire_used and %missed_panIDs_by_chr
#my %merged_consen_table_entire = (%missed_panIDs_by_chr, %consen_table_entire_used);
my %merged_consen_table_entire;
my ($ct_cte, $ct_cteu, $ct_miss);
foreach my $chr  (1 .. $num_chrs ){
  while ( my ($k,$v) = each( %{$consen_table_entire{$chr}} ) ) {
    $merged_consen_table_entire{$chr}{$k} = $v;
    $ct_cte++;
  }
  while ( my ($k,$v) = each( %{$consen_table_entire_used{$chr}} ) ) {
    $merged_consen_table_entire{$chr}{$k} = $v;
    $ct_cteu++;
  }
  while ( my ($k,$v) = each( %{$missed_panIDs_by_chr{$chr}} ) ) {
    $merged_consen_table_entire{$chr}{$k} = $v;
    $ct_miss++;
  }
}
my $ct_mcte;
foreach my $chr (1 .. $num_chrs ){
  while ( my ($k,$v) = each( %{$merged_consen_table_entire{$chr}} ) ) {
    $ct_mcte++;
  }
}
#say "cte: $ct_cte\ncte: $ct_cteu\nmiss: $ct_miss\nmerged: $ct_mcte";

# Final traversal of %consen_table to print the results - now with added elements from %unplaced.
# First put it into an array of arrays, before sorting.
# Renumber the panIDs, since the initial numbers may not have sufficient space to permit
# addition of unplaced genes within the available ordered integers, given the placement scheme.
my @consen_table_AoA;
my $OUT_FH;
if ($outfile){
  open ($OUT_FH, ">", $outfile) or die "Can't open outfile: $outfile\n";
}
for my $chr (keys %merged_consen_table_entire){
  for my $merged_key (keys %{$merged_consen_table_entire{$chr}} ){
    my @parts = split(/__/, $merged_key);
    my ($panID, $chr, $order, $orient) = @parts;
    my $printchr;
    if ($chr_count <= 9){ # don't zero-pad the chromosome number
      $printchr = sprintf("chr%d", $chr);
    }
    else { # pad chromosome number for two digits
      $printchr = sprintf("chr%02d", $chr);
    }
    my @new_record = ($panID, $printchr, $order, $orient);
    #say "($panID, $printchr, $order, $orient)";
    push @consen_table_AoA, \@new_record;
  }
}
my @sorted_consen_table_AoA = sort {
  $a->[1] cmp $b->[1] ||   # chromosome string
  $a->[2] <=> $b->[2]      # position
} @consen_table_AoA;

my $new_ordinals = 0;
%seen_chr = ();
my $prev_panID;
$count = 0;
for my $row (@sorted_consen_table_AoA ){
  my ($panID, $chr, $orig_order, $orient ) = @$row;
  if ($count == 0 ){
    if ($seen_chr{$chr}){
      $new_ordinals += 100;
    }
    else { # not seen_chr
      $seen_chr{$chr}++;
      $new_ordinals = 100;
    }
  }
  elsif ($count > 0 && $panID ne $prev_panID ){
    if ($seen_chr{$chr}){
      $new_ordinals += 100;
    }
    else { # not seen_chr
      $seen_chr{$chr}++;
      $new_ordinals = 100;
    }
  }
  else { next }

  if ($outfile){
    #say $OUT_FH join ("\t", $panID, $chr, $new_ordinals, $orient, $orig_order);
    say $OUT_FH join ("\t", $panID, $chr, $new_ordinals, $orient);
  }
  else {
    #say join ("\t", $panID, $chr, $new_ordinals, $orient, $orig_order);
    say join ("\t", $panID, $chr, $new_ordinals, $orient);
  }
  $prev_panID = $panID;
  $count++;
}

#################### Subroutines

# Find consensus placement of unplaced panID
sub find_placement {
  my ($chr, $panID, $signs_by_chr_ref) = @_;
  my %signs_by_chr = %$signs_by_chr_ref;
  my $transition_idx;
  my $width = 5; # width of slices before and after to evaluate for signs, around position index $i
  my $seen_transition = 0;
  my @elements; 
  # Put hash contents (signs) into array by index(position)
  foreach my $idx ( sort { $a <=> $b } keys %{$signs_by_chr{$chr}{$panID}} ){
    my $sign = ${$signs_by_chr{$chr}{$panID}}{$idx};
    $elements[$idx] = $sign;
  }
  if ($verbose>2){ say join ("", @elements); }
  foreach my $i ( $width-1 .. @elements-($width+1) ){
    my @fore_elts = @elements[ $i-($width-1) .. $i ];
    my @aft_elts = @elements[ $i+1 .. $i+$width ];
    if ($verbose>2){
      say sum(@fore_elts), " ", sum(@aft_elts), "\t@fore_elts   @aft_elts\t$i"; 
    }
    if ( sum(@fore_elts) < 0 && sum(@aft_elts) == $width &&
         abs(sum(@aft_elts)) > abs(sum(@fore_elts)) && 
         $seen_transition == 0 ){
      if ($verbose>2){ say join ("\t", "  target:", $panID, $chr,  "transition_idx:", $i); }
      $transition_idx = $i;
      $seen_transition = 1;
      last;
    }
  }
  return $transition_idx;
}

__END__
2023
S. Cannon
02-14 Initial version, based on order_by_consensus.pl
02-22 Retrieve data structures from Parallel::ForkManager with run_on_finish callback
02-23 Yank Parallel::ForkManager because of inconsistency in retrieval of data structure.
02-25 More testing. Merge original consen_gene_order table with missed and formerly unplaced panIDs.
02-26 Add back ForkManager after restructuring loop and writing to tmp files.
02-27 Fix REGEX for chromosome prefix, removing patterns that can match Mtrun
02-28 Fix bug in determining the consensus orientation.
03-05 Improve runtime feedback and debugging info.
03-11 More checking for uninitialized transition index / panID
12-11 Stop reporting per-gene progress 

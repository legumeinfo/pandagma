#!/usr/bin/env perl

use strict;
use warnings;
use Bio::AlignIO;
use Getopt::Long;
use List::Util qw(max);
use List::MoreUtils qw(uniq);

my ($in_align, $out_align, $log, $require_inform, $help, $verbose);
my $depth = 1; # minimum number of characters for reporting column
my $pct_depth = 0; # minimum percentage of sequence count to use for keeping a column
my $min_pct_aligned = 10; # minimum proportion of aligned characters after trimming to $depth
GetOptions ("in_align=s"        => \$in_align,  # required
            "out_align=s"       => \$out_align,  # required
            "log:s"             => \$log,  
            "depth:i"           => \$depth,
            "pct_depth:i"       => \$pct_depth,
            "min_pct_aligned:i" => \$min_pct_aligned,
            "require_inform"    => \$require_inform,
            "verbose"           => \$verbose,
            "help"              => \$help
            );

my $usage = <<EOS;
Usage: filter_align.pl -in ALIGN.fa -out ALIGN_OUT.fa

  Eliminate columns in a fasta-format alignment that just contain gaps.

  Required: 
    -in_align         filename
    -out_align        filename
  Optional:
    -log              filename
    -depth            minimum number of characters for reporting column; default 1
    -pct_depth        minimum percentage of seq count for reporting column; default 0
                        If both depth and pct_depth, the higher will be used - e.g. for 
                        -depth 4 -pct 20  and a 50-seq alignment, 10 will be in effect
    -min_pct_aligned  minimum percent of aligned characters after trimming to 
                        depth for retaining sequence; default 10
    -require_inform   Require that all columns have difference information, 
                        i.e. have more than one residue or base
    -verbose          print some debugging info
    -help             this message
EOS

die "$usage\n" if ($help or !defined $in_align or !defined $out_align);

my $in = new Bio::AlignIO ( -file => $in_align, -format => 'fasta' );

open my $OUT_FH, ">", $out_align or die "cant open out $out_align: $!\n";

my $LOG_FH;
if ($log) {
  open $LOG_FH, ">", $log or die "cant open out $log: $!\n";
}

my $aln = $in->next_aln();

# Create an array containing all columns as strings
my @aln_cols;
my @aln_col_scores;
my $seq_count;
my %seen_col; # to record if this $aln_col_scores[] has been initialized

eval {
  foreach my $seqobj ( $aln->each_seq() ) {
    my $col_idx = 0;
    $seq_count++;
    foreach my $char ( split("", $seqobj->seq()) ) {
      if ($seen_col{$col_idx}){} # otherwise initialize this $aln_col_scores[$col_idx]
      else {$seen_col{$col_idx}++; $aln_col_scores[$col_idx] = 0} 
      $char =~ s/X/-/; # Replace X character with gap character (dash)
      $aln_cols[$col_idx] .= $char;
      $aln_col_scores[$col_idx]++ unless ($char =~ /-|N/i);
      $col_idx++;
    }
  } 
  if ($@) { # Report trapped error
    print "Trapped error in filter_align.pl, input file $in_align\n";
    die "$@";
  }
};

# Skip columns that have too many gap characters
my $gapchar = $aln->gap_char();
my @no_gap_cols = ();
my $col_idx = 0;
foreach my $col ( @aln_cols ) {
  my $columnlength = length($col);

  my @col_chars = ();
  my $ct_uniq_chars = "";
  my $skip_low_info_col = 0;
  if ($require_inform) {
    @col_chars = split("", $col);
    $ct_uniq_chars = uniq @col_chars;
    if ( $ct_uniq_chars == 1 ){ $skip_low_info_col = 1 }
    else { $skip_low_info_col = 0 }
  }

  # test if chars in this colmn are below $depth or below $pct_depth
  my $min_effective_depth = max( $depth, $seq_count*($pct_depth/100) );
  #print "min_effective_depth: max( $depth, $seq_count*($pct_depth/100) ) = $min_effective_depth\n";
  if ($verbose) { print "$col\t$aln_col_scores[$col_idx]\t$col_idx\t$min_effective_depth\n" }
  if ( $skip_low_info_col == 1 ){
    #print "uniq=$ct_uniq_chars; ??=$skip_low_info_col; depth=$aln_col_scores[$col_idx]; SkipAA $col_idx\n";
    $col_idx++;
    next
  }
  elsif ( $aln_col_scores[$col_idx]<$min_effective_depth ) {
    #print "uniq=$ct_uniq_chars; ??=$skip_low_info_col; depth=$aln_col_scores[$col_idx]; SkipBB $col_idx\n";
    $col_idx++;
    next
  }
  else {
    #print "uniq=$ct_uniq_chars; ??=$skip_low_info_col; depth=$aln_col_scores[$col_idx]; AddZZZ $col_idx\n";
    push @no_gap_cols, $col;
    $col_idx++;
  }
}

# Reconstruct the sequence strings from columnar data
my @seq_strs = ();
foreach my $col ( @no_gap_cols ) {
  my $col_idx = 0;
  foreach my $char ( split("", $col) ) {
    $seq_strs[$col_idx] .= $char;
    $col_idx++; 
  }
}

# replace the old sequences strings with the new ones
my $i=0;
eval {
  foreach my $seqobj ( $aln->each_seq ) {
    $seqobj->seq($seq_strs[$i++]);
    my $display_id = $seqobj->display_id();
    my $sequence = $seqobj->seq();
    if (!defined($sequence)){
      warn "WARNING -- In $in_align, No sequence found for $display_id; skipping\n";
      next;
    }
  
    my $desc = $seqobj->desc();
    my $seq_len = length($sequence);
    if (!defined $seq_len || $seq_len == 0){
      warn "Got a sequence with no length for alignment $in_align, ID $display_id\n";
      next;
    }
    else {
      my $ct_dashes += () = $sequence =~ /-/g;
      my $ct_residues = $seq_len - $ct_dashes;
      my $pct_aligned_residues = int(10000*$ct_residues/$seq_len)/100;
      if ($verbose) { print "$min_pct_aligned\t$ct_residues\t$ct_dashes\t$pct_aligned_residues\n"; }
      if ($pct_aligned_residues < $min_pct_aligned) {
        if (defined $log) {
          print $LOG_FH "For alignment $in_align, cutting $display_id: only " .
            "$ct_residues/$seq_len = $pct_aligned_residues% of $min_pct_aligned% " .
            "required after trimming to align depth $depth and percent depth $pct_depth\n";
        }
        else {
          print "For alignment $in_align, cutting $display_id: only " .
            "$ct_residues/$seq_len = $pct_aligned_residues% of $min_pct_aligned% " .  
            "required after trimming to align depth $depth and percent depth $pct_depth\n";
        }
      }
      else {
        if (defined $desc) {
          print $OUT_FH ">$display_id $desc\n$sequence\n";
        }
        else {
          print $OUT_FH ">$display_id\n$sequence\n";
        }
      }
    }
  }
  if ($@) { # Report trapped error
    print "Trapped error in filter_align.pl, input file $in_align\n";
    die "$@";
  }
};

__END__
VERSIONS
sc = Steven Cannon
picked up this script on-line somewhere, ca. 2010
17-07-17 sc Report usage message, and print to stdout
17-08-01 sc Rewrite. Add getopts. Add filtering by alignment depth.
17-08-20 sc Add warning about gap-only sequences after trimming to given align depth
            Add min_pct_aligned to cut short sequences
            Add option to remove short sequences. Write out to specified outfile.
            Add log file.
18-02-23 sc Add option to filter by depth by percent of sequences
19-09-30 sc Also filter against Ns (in addition to dashes)
23-01-22 sc Remove columns that have no information (are all the same residue or base)
23-10-09 sc Test if seq_len is defined and >0
23-11-15 sc Handle case of missing sequence
24-09-19 sc Trap and report error when traversing $seqobj->seq()

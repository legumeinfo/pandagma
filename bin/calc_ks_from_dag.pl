#!/usr/bin/env perl

# Program: calc_ks_from_dag.pl
# Program description: see usage message.
# Steven Cannon 2023

use strict;
use warnings;
use vars qw($CODONSIZE $KS_CUTOFF $KA_CUTOFF);
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::MAFFT;
use Bio::Tools::CodonTable;
use Bio::Align::DNAStatistics;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Phylo::PAML::Yn00;
use Getopt::Long;
use File::Basename;
use feature "say";
#use Data::Dumper;
#use Carp;

my $scriptname = basename($0);

my $usage = <<EOS;
  Usage: $scriptname FASTA_FILE[S] -dagin FILE.aligncoords -report_out FILE [options]
  
  Calculates Ks for genes in DAGchainer synteny gene pairs.
  
  Required:
    [on ARGV]   -- One or more fasta files with CDS sequence containing genes from FILE.aligncoords
    -dagin      -- An .aligncoords file, from DAGchainer
                    IDchrA idA start stop IDchrB idB start stop Eval
    -report_out -- filename for output report
    
  Options:
    -aln_out              name for optional CDS alignment file, e.g. chrA.chrB.outcds
    -h/--help             See this information
    -v/--verbose          Run in verbose mode
EOS

my ($CODONSIZE, $KA_CUTOFF, $KS_CUTOFF, $frame, $codontable, $verbose) = (3, 4, 4, 0, 1, 0);
my ($dagin, $chrA_fas, $chrB_fas, $report_out, $aln_out, $help);
my $align_method = "clustalw";

GetOptions( 
  'dagin=s'           => \$dagin,      # required
  'report_out=s'      => \$report_out, # required
  'aln_out:s'         => \$aln_out,
  'h|help'            => \$help,
  'am|align_method:s' => \$align_method,
  'v|verbose+'        => \$verbose,
);

if ($help) { die "\n$usage" }
if (not defined($dagin))      { die "\n  -dagin is required. See -help for guidance.\n\n" }
if (not defined($report_out)) { die "\n  -report_out is required. See -help for guidance.\n\n" }
if (scalar(@ARGV)==0){ die "\nPlease provide FASTA_FILE(S) as the first argument(s). See -help for guidance.\n\n" }

##################################################
# Read fasta file(s) into ID-sequence hash
my %seq_hsh;
foreach my $input_fas ( @ARGV ) {
  chomp($input_fas);
  say "    Processing FASTA: $input_fas";

  my $seqio_obj;
  if ($input_fas =~ /.gz$/) {
    open my $fh, "gzip -dc $input_fas |" or die $!;
    $seqio_obj = Bio::SeqIO->new(-fh => $fh, -format => 'Fasta');
  } else {
    $seqio_obj = Bio::SeqIO->new(-file => $input_fas, -format => 'Fasta');
  }

  while ( my $seq_obj = $seqio_obj->next_seq ) {
    my $display_id = $seq_obj->display_id();
    #$seq_hsh{$display_id} = $seq_obj->seq();
    $seq_hsh{$display_id} = $seq_obj;
  }
}

##################################################
# Process the DAGchainer output (FILE.aligncoords) 
open (my $DAG, "< $dagin") or die "can't open $dagin: $!";
open (my $RPTOUT, "> $report_out") or die "cannot open out $report_out:$!";

my $align_obj;
if ($aln_out) { $align_obj = new Bio::AlignIO('-format' => "fasta", '-file' => ">$aln_out") }

my $table = new Bio::Tools::CodonTable();
my $alignengine;
if ($align_method =~ /clustalw/i){
  $alignengine = Bio::Tools::Run::Alignment::Clustalw->new('ktuple' => 2, 'matrix' => 'BLOSUM');
}
elsif ($align_method =~ /mafft/i){
  $alignengine = Bio::Tools::Run::Alignment::MAFFT->new();
}
else {
  die "Alignment method is not recognized. Available methods: clustalw and mafft";
}

my $count = 0;
my $dag_header = "nullstring";
my @dag_header_bits;
my (@out_KaKs_ary_all, $out_KaKs_aryref_all);
my (@out_Ka_ary, @out_Ks_ary, @out_KaKs_ary, $out_Ka_aryref, $out_Ks_aryref, $out_KaKs_aryref);
my (%out_KaKsSum_hsh, $out_KaKsSum_hshref);
my ($dag_startA, $dag_stopA, $dag_startB, $dag_stopB) = (999999999999999999,0,999999999999999999,0);

while (<$DAG>) {
  chomp;
  my $line = $_;
  
  if ($line =~ /^#/) {
    # Process and print summary header lines from previous DAGchainer block
    my ($ave_Ka, $ave_Ks, $ave_KaKs, $median_Ka, $median_Ks, $median_KaKs);
    if ($count > 0) { # means $ave_Ka, $ave_Ks were calculable, so $count incremented; prepare for printing.
      $ave_Ka =   sprintf ("%.4f", ($out_KaKsSum_hshref->{"sum_Ka"})/$count );
      $ave_Ks =   sprintf ("%.4f", ($out_KaKsSum_hshref->{"sum_Ks"})/$count );
      $ave_KaKs = sprintf ("%.4f", ($out_KaKsSum_hshref->{"sum_KaKs"})/$count );
      
      $median_Ka =   sprintf ("%.4f", calc_median(\@$out_Ka_aryref));
      $median_Ks =   sprintf ("%.4f", calc_median(\@$out_Ks_aryref));
      $median_KaKs =   sprintf ("%.4f", calc_median(\@$out_KaKs_aryref));
    }
    else { # means $ave_Ka, $ave_Ks weren't calculable, so $count didn't increment. Flag with 999s.  
      ($ave_Ka, $ave_Ks, $ave_KaKs, $median_Ka, $median_Ks, $median_KaKs) = (999, 999, 999, 999, 999, 999);
    }
    
    my ($chrA, $chrB, $diag_orient, $diag_num, $diag_score, $diag_ct) = @dag_header_bits;
    # GFF:  seqid  source  type  start  end  score  strand  phase  attributes
    unless ($dag_header eq "nullstring") {
      say $RPTOUT  "##GFF_A\t$chrA\tDAGchainer\tsynteny\t$dag_startA\t$dag_stopA\t$diag_score\t$diag_orient\t." . 
                       "\tID=$chrA.$chrB.$diag_num.$diag_orient;median_Ks=$median_Ks;ave_Ks=$ave_Ks;matches=$chrB:$dag_startB..$dag_stopB";
      say $RPTOUT  "##GFF_B\t$chrB\tDAGchainer\tsynteny\t$dag_startB\t$dag_stopB\t$diag_score\t$diag_orient\t." .
                       "\tID=$chrB.$chrA.$diag_num.$diag_orient;median_Ks=$median_Ks;ave_Ks=$ave_Ks;matches=$chrA:$dag_startA..$dag_stopA";
      say $RPTOUT  "#%SUM\tcount\tdag_startA\tdag_stopA\tdag_startB\tdag_stopB\t" . 
                       "ave_Ka\tave_Ks\tave_KaKs\tmedian_Ka\tmedian_Ks\tmedian_KaKs";
      say $RPTOUT  "##SUM\t$count\t$dag_startA\t$dag_stopA\t$dag_startB\t$dag_stopB\t" . 
                       "$ave_Ka\t$ave_Ks\t$ave_KaKs\t$median_Ka\t$median_Ks\t$median_KaKs";
      say $RPTOUT  "#%DATA\tA_id\tB_id\taln_len\tKa\tKs\tKaKs\tblock_Ks";
    }

    # Print a line for each gene-pair-data-line from the aligncoords file. 
    # At the end of each line, add median_Ks for the block, to allow filtering in post-processing.
    foreach my $data_line (@$out_KaKs_aryref_all){ 
      say $RPTOUT  "$data_line\t$median_Ks";
    }
    
    $out_KaKsSum_hshref->{"sum_Ka"} = 0;
    $out_KaKsSum_hshref->{"sum_Ks"} = 0;
    $out_KaKsSum_hshref->{"sum_KaKs"} = 0;
    
    @$out_KaKs_aryref = ();
    @$out_Ka_aryref = ();
    @$out_Ks_aryref = ();
    @$out_KaKs_aryref_all = ();
    
    $count = 0;
    ($dag_startA, $dag_stopA, $dag_startB, $dag_stopB) = (999999999999999999,0,999999999999999999,0);
    
    # DAG headers look like one of the two following:
    ## alignment Gm10 vs. Gm20 (reverse) Alignment #3  score = 348.9 (num aligned pairs: 10):
    ## alignment Gm10 vs. Gm20 Alignment #3  score = 1458.8 (num aligned pairs: 32):
    
    $dag_header = $line;
    if ($dag_header =~ /reverse/) { $diag_orient = "-" } else { $diag_orient = "+" }
    $dag_header =~ /alignment ([^ ]+) vs. ([^ ]+) .*Alignment #(\d+) +score = (\d+.\d+) .num aligned pairs: (\d+).:/;
    ($chrA, $chrB, $diag_num, $diag_score, $diag_ct) = ($1, $2, "$3", $4, $5);
    @dag_header_bits = ($chrA, $chrB, $diag_orient, $diag_num, $diag_score, $diag_ct);
    
    next
  }
  else { # Process data line and print Ka and Ks results.
    #say "AA: $line";
    
    my ($IDchrA, $idA, $startA, $stopA, $IDchrB, $idB, $startB, $stopB, $Eval) = split /\s+/, $line;
  
    $dag_startA = $startA unless ($dag_startA < $startA);
    $dag_stopA  = $stopA unless  ($dag_stopA >  $stopA);    
    $dag_startB = $startB unless ($dag_startB < $startB);
    $dag_stopB  = $stopB unless  ($dag_stopB >  $stopB);    
    
    my ($nuc_objA, $nuc_objB) = ($seq_hsh{$idA}, $seq_hsh{$idB});
    
    if ($verbose) {
      say "\ncount:\t$count";
      say "\tIDchrA\t[idA]\tstartA\t[chrAhsh]";
      say "\t$IDchrA, [$idA], $startA, [$nuc_objA]";
      say "\tIDchrB\t[idB]\tstartB\t[chrBhsh]";
      say "\t$IDchrB, [$idB], $startB, [$nuc_objB]";
    }
    
    my $prot_objA = $nuc_objA->translate(-frame => 0, -codontable_id => $codontable);
    my $prot_objB = $nuc_objB->translate(-frame => 0, -codontable_id => $codontable);
   
    ## Align sequences
    my $dna_aln;
    eval {
      if ($verbose) { say "ALIGNING $nuc_objA, $nuc_objB, $prot_objA, $prot_objB"; }
      $dna_aln = align_pair($nuc_objA, $nuc_objB, $prot_objA, $prot_objB);
    }; 
    warn $@ if $@;
    
    ## Calculate and report Ka and Ks
    eval {
      $count++;
      if ($verbose) { say "CALCULATING KA & KS for $nuc_objA, $nuc_objB, $prot_objA, $prot_objB"; }
      ($out_KaKs_aryref_all, $out_KaKsSum_hshref, $out_Ka_aryref, $out_Ks_aryref, $out_KaKs_aryref, $count) = 
        KaKs_report($dna_aln, $nuc_objA, $nuc_objB, $count);
    }; 
    warn $@ if $@;
  }
}

# flush last results
say $RPTOUT "$dag_header"; 

my ($ave_Ka, $ave_Ks, $ave_KaKs, $median_Ka, $median_Ks, $median_KaKs);
if ($count > 0) { # means $ave_Ka, $ave_Ks were calculable, so $count incremented; prepare for printing.
  $ave_Ka =   sprintf ("%.4f", ($out_KaKsSum_hshref->{"sum_Ka"})/$count );
  $ave_Ks =   sprintf ("%.4f", ($out_KaKsSum_hshref->{"sum_Ks"})/$count );
  $ave_KaKs = sprintf ("%.4f", ($out_KaKsSum_hshref->{"sum_KaKs"})/$count );
  
  $median_Ka =   sprintf ("%.4f", calc_median(\@$out_Ka_aryref));
  $median_Ks =   sprintf ("%.4f", calc_median(\@$out_Ks_aryref));
  $median_KaKs = sprintf ("%.4f", calc_median(\@$out_KaKs_aryref));
}
else { # means $ave_Ka, $ave_Ks weren't calculable, so $count didn't increment. Flag with 999s.  
  ($ave_Ka, $ave_Ks, $ave_KaKs, $median_Ka, $median_Ks, $median_KaKs) = (999, 999, 999, 999, 999, 999);
}

my ($chrA, $chrB, $diag_orient, $diag_num, $diag_score, $diag_ct) = @dag_header_bits;
# GFF:  seqid  source  type  start  end  score  strand  phase  attributes
say $RPTOUT  "##GFF_A\t$chrA\tDAGchainer\tsynteny\t$dag_startA\t$dag_stopA\t$diag_score\t$diag_orient\t." . 
                 "\tID=$chrA.$chrB.$diag_num.$diag_orient;median_Ks=$median_Ks;ave_Ks=$ave_Ks;matches=$chrB:$dag_startB..$dag_stopB";
say $RPTOUT  "##GFF_B\t$chrB\tDAGchainer\tsynteny\t$dag_startB\t$dag_stopB\t$diag_score\t$diag_orient\t." .
                 "\tID=$chrB.$chrA.$diag_num.$diag_orient;median_Ks=$median_Ks;ave_Ks=$ave_Ks;matches=$chrA:$dag_startA..$dag_stopA";
say $RPTOUT  "#%SUM\tcount\tdag_startA\tdag_stopA\tdag_startB\tdag_stopB\t" . 
                 "ave_Ka\tave_Ks\tave_KaKs\tmedian_Ka\tmedian_Ks\tmedian_KaKs";
say $RPTOUT  "##SUM\t$count\t$dag_startA\t$dag_stopA\t$dag_startB\t$dag_stopB\t" . 
                 "$ave_Ka\t$ave_Ks\t$ave_KaKs\t$median_Ka\t$median_Ks\t$median_KaKs";
say $RPTOUT  "#%DATA\tA_id\tB_id\taln_len\tKa\tKs\tKaKs\tblock_Ks";

# Print a line for each gene-pair-data-line from the aligncoords file. 
# At the end of each line, add median_Ks for the block, to allow filtering in post-processing.
foreach my $data_line (@$out_KaKs_aryref_all){ 
  say $RPTOUT  "$data_line\t$median_Ks"; 
}

####################################################################################################
######### subroutines ##############################################################################

sub align_pair {
  my ($nuc_objA, $nuc_objB, $prot_objA, $prot_objB) = @_;
  
  my @nucpair = ($nuc_objA, $nuc_objB);
  my @protpair = ($prot_objA, $prot_objB);
      
  my $aa_aln = $alignengine->align(\@protpair);
  
  my $dna_aln = new Bio::SimpleAlign;
  my $seqorder = 0;
  my $aa_aln_len = $aa_aln->length;
  foreach my $seq ( $aa_aln->each_seq ) {    
    my $newseq;
    
    foreach my $pos ( 1..$aa_aln_len ) { 
      my $loc = $seq->location_from_column($pos);
      my $dna = ''; 
      if( !defined $loc || $loc->location_type ne 'EXACT' ) {
        $dna = '---';
      } 
      else {
        my $aa_start = $loc->start;
        my $aa_end = $loc->end;
        
        # to readjust to codon boundaries, end needs to be +1 so we can just multiply by CODONSIZE to get this:
        my $nt_start = (($aa_start - 1)*$CODONSIZE) +1;
        my $nt_end = ($aa_end)*$CODONSIZE;
        
        if ( $nt_start <=0 || $nt_end > $nuc_objB->length() ) {
          if ( $verbose>1 ) {
              say "\tcodons don't match for ($aa_start,$aa_end) ($nt_start,$nt_end)";
          }
          $dna = '---';
        } 
        else {
            $dna = $nucpair[$seqorder]->subseq($nt_start,$nt_end);
        }
      }
      $newseq .= $dna;
    }
    $seqorder++;
    
    # Readjust to codon boundaries (note sequences start with 1)
    my $new_dna = new Bio::LocatableSeq(
                 -display_id  => $seq->id(),
                 -start => (($seq->start - 1) * $CODONSIZE) + 1, 
                 -end   => ($seq->end * $CODONSIZE),
                 -strand => $seq->strand,
                 -seq   => $newseq
         );
                     
    $dna_aln->add_seq($new_dna);
  }
  
  if ($verbose) { say "WRITING $dna_aln" }
  if ($aln_out) { $align_obj->write_aln($dna_aln) }
  
  return $dna_aln;
}

sub KaKs_report {
  my ($dna_aln, $nuc_objA, $nuc_objB, $count) = @_;
  #say "BB: $dna_aln, $nuc_objA, $nuc_objB, $count";
  
  # accumulate string for the report
  my $A_id = $nuc_objA->display_id;
  my $B_id = $nuc_objB->display_id;
  my $aln_len = $dna_aln->length;

  my ($kaks_factory, $Ka, $Ks, $KaKs);
  
  # Use codeml method (yn00 is't handled in this script yet; probably no reason to.)
  $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => -2, 'seqtype' => 1 } );
  $kaks_factory->alignment($dna_aln);
  my ($rc, $parser) = $kaks_factory->run();
  my $result = $parser->next_result;
  my $MLmatrix = $result->get_MLmatrix();
  
  $Ka = $MLmatrix->[0]->[1]->{'dN'};
  $Ks = $MLmatrix->[0]->[1]->{'dS'};
  $KaKs = $MLmatrix->[0]->[1]->{'omega'};
  #say "DD: $A_id\t$B_id\t$aln_len\t$Ka\t$Ks\t$KaKs";
  
  if ($Ka < $KA_CUTOFF) {
    $out_KaKsSum_hsh{"sum_Ka"} += $Ka;
    push @out_Ka_ary, $Ka;
  }
  if ($Ks < $KS_CUTOFF) {
    $out_KaKsSum_hsh{"sum_Ks"} += $Ks;
    push @out_Ks_ary, $Ks;
  }
  if ($Ka < $KA_CUTOFF and $Ks < $KS_CUTOFF) {
    $out_KaKsSum_hsh{"sum_KaKs"} += $KaKs;
    push @out_KaKs_ary, $KaKs;
  }
  if ($Ka >= $KA_CUTOFF or $Ks >= $KS_CUTOFF) {
    $count--; # high Ks isn't included in out_KaKsSum_hsh, so de-increment
  }
  
  push @out_KaKs_ary_all, "\t$A_id\t$B_id\t$aln_len\t$Ka\t$Ks\t$KaKs";
  #say "EE: \t$A_id\t$B_id\t$aln_len\t$Ka\t$Ks\t$KaKs";
  
  return (\@out_KaKs_ary_all, \%out_KaKsSum_hsh, \@out_Ka_ary, \@out_Ks_ary, \@out_KaKs_ary, $count);
}

# Return median value for an array of numbers
# See http://stackoverflow.com/questions/5119034/using-perl-to-find-median-mode-standard-deviation
sub calc_median {
  my ($array_ref) = @_;
  my @values = @$array_ref;
  my $median;
  my $mid = int((scalar @values)/2);
  my $len = scalar(@values);
  my @sorted_values = sort {$a <=> $b} @values;
  #say Dumper(@sorted_values);
  if (@values % 2) {
    $median = $sorted_values[ $mid ];
  } else {
    $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
  }
  #say "FF: Median is $median";
  return $median;
}

__END__

Versions
2023-08-07 Initial, based on older dag_ks.pl script.
2023-08-26 Tweak help text: STDIN --> ARGV

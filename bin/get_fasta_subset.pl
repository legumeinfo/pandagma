#!/usr/bin/env perl

# PROGRAM: get_fasta_subset.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2005
# see description under Usage

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use feature "say";

####### File IO ########

my ($input_fas, $list_IDs, $help, $verbose, $output_fas);
my $xclude = 0;
my $fasta_order = 0;
my $clobber = 0;

GetOptions (
  "input_fas=s" =>  \$input_fas,   # required
  "output_fas=s" => \$output_fas,   # required
  "list_IDs=s" =>   \$list_IDs,   # required
  "xclude" =>       \$xclude, 
  "fasta_order" =>  \$fasta_order,
  "clobber" =>      \$clobber,
  "verbose+" =>     \$verbose,
  "help" =>         \$help
);

my $scriptname = basename($0);

my $usage = <<EOS;
  Usage: $scriptname [-options]
  
  Read a list of fasta IDs and a fasta file and either exclude or include those
  sequences into a new fasta file.
   
   -input_fas: input fasta file **
   -output_fas: output fasta file **
   -list_IDs: list of IDs to include or exclude **
   -xclude: boolean 
      -x (True) means exclude seqs in list 
         (False; default) means include sequences in the list 
   -fasta_order: boolean
      -f (True) means retain the order in the input fasta file
         (False; default) means report sequences in list order.
           Has no effect if -xclude is true.
   -clobber    (boolean) clobber (overwrite) existing output file
   -verbose    (boolean) for more stuff to terminal; may call multiple times
   -help       for more info
   
   ** = required
EOS

die "\n$usage\n" 
  if ($help or !defined($input_fas) or !defined($output_fas) or !defined($list_IDs) );

####### The program ########
# read hash in
open( my $LIST_FH, '<', $list_IDs ) or die "can't open list_IDs $list_IDs: $!";

# Open file for writing new fasta out
if (-e $output_fas && not($clobber) ) {die "Output file exists.\n"} # don't clobber existing file
open( my $OUT_FH, '>', $output_fas ) or die "can't open output_fas $output_fas: $!"; 

# Put elements of LIST into array in case we need to retain list order
my @list_ary;
# Put elements of LIST into hash; they may be ordered by array later.
my %hash;

my $ct=0;
while (<$LIST_FH>) {
  chomp;
  next if ( $_ =~ /^\s*$/ );
  my $id = $_;
  $id =~ s/\s+$//; # strip trailing whitespace
  #say "II id: [$id]";
  $hash{$id} = $id; 
  $list_ary[$ct] = $id;
  $ct++;
}

# Read in the sequence using the Bioperl SeqIO;
my $in  = Bio::SeqIO->new(-file => $input_fas , '-format' => 'Fasta');

my (%seq_hsh, %desc_hsh);

# Load the sequence into a Bio::Seq object
while ( my $seq = $in->next_seq ) {

  # get parts of the fasta seq
  my $display_id = $seq->display_id();
  my $desc = $seq->desc();
  my $sequence = $seq->seq();
  
  if ($xclude) {
    unless (defined($hash{$display_id})) {
      #say "XX display_id: {$display_id}";
      say $OUT_FH ">$display_id $desc\n$sequence";
    }
  }
  else { # $xclude is false; therefore, include seqs in the list
    if (defined($hash{$display_id})) {
      if (defined($desc)) { # there IS a $desc
        if ($fasta_order) {
          #say "AA display_id: {$display_id}";
          say $OUT_FH ">$display_id $desc\n$sequence";
        }
        else { # not in $fasta_order but in list order
          #say "BB display_id: {$display_id}";
          $seq_hsh{$display_id} = $sequence;
          $desc_hsh{$display_id} = $desc;
        }
      }
      else { # there is NOT a $desc
        if ($fasta_order) {
          #say "CC display_id: {$display_id}";
          say $OUT_FH ">$display_id\n$sequence";
        }
        else { # not in $fasta_order but in list order
          #say "DD display_id: {$display_id}";
          $seq_hsh{$display_id} = $sequence;
          $desc_hsh{$display_id} = " ";
        }
      }
    }
  }
}

if ((not $fasta_order) and (not $xclude)) { # Report sequences in list order
  my $ct=0;
  foreach my $id (@list_ary) {
    #say "YY id $ct: [$id]";
    if (exists $desc_hsh{$id}) {
      my $defline = ">$id $desc_hsh{$id}";
      $defline =~ s/\s+$//; # strip whitespace
      say $OUT_FH "$defline\n$seq_hsh{$id}";
    }
    else { # no desc
      say $OUT_FH ">$id\n$seq_hsh{$id}";
    }
    $ct++;
  }
}

__END__
VERSIONS

2005 S. Cannon. 
2007-05-05 Add getopts; clean up a bit
2015-02-18 Add option to report sequences either in fasta or list order
2015-07-13 don't clobber existing fasta file
2017-03-28 When reading list into file, skip blank lines. Allow clobbering, with -c
2023-03-12 Strip trailing whitespace in list and in fasta output. Add some intermediate/debugging statements.


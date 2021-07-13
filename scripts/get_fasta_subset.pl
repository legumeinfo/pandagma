#!/usr/bin/env perl
# #!/usr/bin/env perl

# PROGRAM: get_fasta_subset.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2005
# see description under Usage

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

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

while (<$LIST_FH>) {
  chomp;
  next if ( $_ =~ /^$/ );
  if (defined $_) {
    $hash{$_} = $_; 
    $list_ary[$.] = $_;
  }
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
    # print "{$display_id} ";
    unless (defined($hash{$display_id})) {
      print $OUT_FH ">$display_id $desc\n$sequence\n";
    }
  }
  else { # $xclude is false; therefore, include seqs in the list
    if (defined($hash{$display_id})) {
      if (defined($desc)) { # there IS a $desc
        if ($fasta_order) {
          print $OUT_FH ">$display_id $desc\n$sequence\n";
        }
        else { # not in $fasta_order but in list order
          $seq_hsh{$display_id} = $sequence;
          $desc_hsh{$display_id} = $desc;
        }
      }
      else { # there is NOT a $desc
        if ($fasta_order) {
          print $OUT_FH ">$display_id\n$sequence\n";
        }
        else { # not in $fasta_order but in list order
          $seq_hsh{$display_id} = $sequence;
          $desc_hsh{$display_id} = " ";
        }
      }
    }
  }
}

if ((not $fasta_order) and (not $xclude)) { # Report sequences in list order
  foreach my $id (@list_ary) {
    if (defined $id) {
      print $OUT_FH ">$id $desc_hsh{$id}\n$seq_hsh{$id}\n";
    }
  }
}

__END__
VERSIONS

v0.01 2005 S. Cannon. 
v0.02 May05'07 SC Add getopts; clean up a bit
v0.03 2015-02-18 SC Add option to report sequences either in fasta or list order
v0.04 2015-07-13 SC don't clobber existing fasta file
v0.05 2017-03-28 SC When reading list into file, skip blank lines. Allow clobbering, with -c

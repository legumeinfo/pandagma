#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = <<EOS;
 Synopsis: hash_into_fasta_id.pl [options] -hash FILE -fasta FILE
  
  Read a key-value hash file, and a fasta file (may be compressed or not).
  Swap the IDs with the values from the hash file.
  The hash file can be created to match the full ID, including splice variants,
  OR created to match just the gene IDs, excluding splice variants.
  In the latter case, provide -splice_regex to give the form of the splice variant.
  
  Required:
  -hash      Key-value hash file, where first column has IDs from fasta file.
  -fasta     Fasta file; may be compressed or not.
  
  Options:
  -outfile   Output filename. If omitted, output is to STDOUT.
  -nodef          (boolean) Don't print the def-line description (print only the ID).
  -splice_regex   (string) regular expression to use to exclude the 
                  splice variant suffix of a feature name during the match. 
         Example 1: For transcripts like Gene1234.1,      use "\\.\\d+"  
         Example 2: For transcripts like Gene1234-mRNA-1, use "-mRNA-\\d+" 
         Example 3: For proteins like    Gene1234.1.p,    use "\\.\\d+\\.p"
  -strip_regex    (string) regular expression to use for stripping unwanted characters
                  from splice variant.
         Example 1: For transcripts like Gene1234.1,      Omit. Nothing to strip.
         Example 2: For transcripts like Gene1234-mRNA-1, use "-mRNA" (note "-" left as delimiter)
         Example 3: For proteins like    Gene1234.1.p,    use "\\.p"
  
  -help   (boolean) This message.
  
EOS

my ($hash_file, $fasta_file, $out_file, $splice_regex, $strip_regex, $nodef, $help, $SPL_RX, $STR_RX);

GetOptions (
  "hash_file=s" =>     \$hash_file,   # required
  "fasta_file=s" =>    \$fasta_file,  # required
  "splice_regex:s" =>  \$splice_regex,   
  "strip_regex:s" =>   \$strip_regex,   
  "out_file:s" =>      \$out_file,   
  "nodef" =>           \$nodef,   
  "help" =>            \$help,
);

die "$usage\n" unless (defined($hash_file) && defined($fasta_file));
die "$usage\n" if ($help);

if ( $splice_regex ){ $SPL_RX=qr/$splice_regex/ }
else { $SPL_RX=qr/$/ }

if ( $strip_regex ){ $STR_RX=qr/$strip_regex/ }
else { $STR_RX=qr/$/ }

# read hash in
open ( my $HSH, '<', $hash_file ) or die "can't open in input_hash, $hash_file: $!";
my %hash;
while (<$HSH>) {
  chomp;
  /(\S+)\s+(.+)/;
  my ($id, $hash_val) = ($1,$2);
  $hash{$id} = $hash_val;   # print "$id, $hash{$id}\n";
}

# Read in the sequence 
my ($FASTA_FH, $OUT_FH);
if ( $fasta_file =~ /gz$/ ){
  open( $FASTA_FH, "zcat $fasta_file|" ) or die "Can't do zcat $fasta_file| : $!";
}
else {
  open ( $FASTA_FH, "<", $fasta_file ) or die "Can't open in $fasta_file: $!\n";
}
if ( $out_file ){
  open ( $OUT_FH, ">", $out_file ) or die "Can't open out $out_file: $!\n";
}
else {
  open ( $OUT_FH, ">&", \*STDOUT) or die;
}
while ( <$FASTA_FH> ){
  chomp;
  my $line = $_;
  my ($display_id, $desc, $seq, $base_id, $suffix);
  if ($line =~ /^>(\S+) +(\S.+)/){
    $display_id = $1;
    $desc = $2;
    $seq = "";
    #print "1:[$display_id] {$desc}\n";

    # strip off splice variant
    $display_id =~ m/(.+)($SPL_RX)$/;
    ($base_id, $suffix) = ($1, $2);
    if ($strip_regex){
      $suffix =~ s/$STR_RX//; 
      $suffix =~ s/-/./; 
    }
    #print "[$base_id] [$suffix]\n";
    
    $hash{$base_id} = "$base_id HASH UNDEFINED" unless defined ($hash{$base_id});
    if ($nodef){ # DON'T print the defline description
      print $OUT_FH ">$hash{$base_id}$suffix\n";
    }
    else { # DO print the defline description
      print $OUT_FH ">$hash{$base_id}$suffix $desc\n";
    }
  }
  elsif ($line =~ /^>(\S+) *$/){
    $display_id = $1;
    $seq = "";
    #print "2:[$display_id]\n";

    # strip off splice variant
    $display_id =~ m/(.+)($SPL_RX)$/;
    ($base_id, $suffix) = ($1, $2);
    #print "[$base_id] [$suffix]\n";
    
    $hash{$base_id} = "$base_id HASH UNDEFINED" unless defined ($hash{$base_id});
    print $OUT_FH ">$hash{$base_id}$suffix\n";
  }
  elsif ($line !~ /^>/){
    $seq .= $line;
    print $OUT_FH "$seq\n";
  }
}

__END__

# Steven Cannon 

Versions
v01 2014-05-21 New script, derived from hash_into_fasta_description.pl
v02 2018-02-09 Handle suffixes (e.g. for splice variants)
v03 2019-05-07 Print original ID if no hash is found
v04 2021-11-01 Don't print final ">" without ID or sequence!
v05 2021-11-04 Add warning for undefined hash
v06 2022-10-04 Remove BioPerl dependency, and take fasta in via STDIN.
                Change handling of the splice variant matching.
v07 2022-10-11 Take sequence file as parameter, handling compressed and uncompressed files.
                Print to named file or to STDOUT.
                Add flag "-strip_regex" to remove e.g. ".p" from protein gene IDs: Gene123.1.p --> Gene123.1

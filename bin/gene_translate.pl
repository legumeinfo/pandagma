#!/usr/bin/env perl

# Program: gene_translate.pl
# Program description: see usage message.
# Steven Cannon 2023

use strict;
use warnings;
use DB_File;
use Getopt::Long;
use feature "say";
# Load module JSON below if -json flag is set

my $gn_to_pan = "gn_to_pan.db";
my $pan_to_gn = "pan_to_gn.db";
my $hash_method = "BerkeleyDB";
my $to_annot = "all";
my $json = 0;
my ($gensp, $strain, $gnm_ver, $ann_ver);
my $help;

my $usage = <<EOS;
  Usage: cat gene-list | gene_translate.pl [-options]
  
  Given a stream of "FROM" genes, return a table of "FROM - TO" correspondences and pangene IDs.
  
  In the default mode, the FROM list is prefixedof the form gensp.strain.gnm#.ann# :
    glyma.Zh13.gnm2.ann1.SoyZH13_14G044301
    glyma.Zh13.gnm2.ann1.SoyZH13_09G142600

  Alternatively, provide a list of the unprefixed base identifiers:
    SoyZH13_14G044301
    SoyZH13_09G142600
  and, with flags, provide all other components that specify the annotation version:
    -gensp glyma     (first three leters of the genus + first two letters of the species)
    -strain Zh13
    -gnm_ver 2 
    -ann_ver 1

  Also required: paths to two hash files, preferably as Berkeley DB files**:
    one, a hash of   gene  panID    (default gn_to_pan.db)
    another, of      panID gene1,gene2,gene3,...  (default pan_to_gn.db)

  ** (The other supported alternative is plain-text hash files of the same format.
      The plain-text option is about 100 times slower than the BerkelyDB option.)

  Options: 
    -gn_to_pan   [gn_to_pan.db] Path to hash file of geneIDs and pangeneIDs, 
                   either in BerkeleyDB or tsv hash format. 
    -pan_to_gn   [pan_to_gn.db] Path to hash file of pangeneIDs and the list of all geneIDs in 
                   that pangene, either in BerkeleyDB or tsv hash format. 
                   In either case, the gene list should be comma-separated -
                   the panID as the key and the value as the collection of genes.
                   For the tsv option:
                     panID  gene1,gene2,gene3,...
                   For the BerkeleyDB option, the key and value alternate:
                     panID
                     gene1,gene2,gene3,...
    -hash_method [BerkeleyDB] BerkeleyDB or tsv
    -to_annot    [all] Annotation string for genes to report - for example, "glyma.Wm82.gnm4.ann1"
                   If not provided or if specified as "all", 
                   all genes in the matching pangenes will be reported.
    -gensp       If input isn't prefixed, supply first three leters of the genus + first two of species
    -strain      If input isn't prefixed, supply strain, e.g. Wm82
    -gnm_ver     If input isn't prefixed, supply genome version, e.g. 2
    -ann_ver     If input isn't prefixed, supply annotation version, e.g. 1
    -json        boolean; default false. If set, encode output in json.
    -help        This message.
EOS

GetOptions(
  'gn_to_pan:s'   => \$gn_to_pan,  
  'pan_to_gn:s'   => \$pan_to_gn,  
  'hash_method:s' => \$hash_method,
  'to_annot:s'    => \$to_annot,
  'gensp:s'       => \$gensp,
  'strain:s'      => \$strain,
  'gnm_ver:s'     => \$gnm_ver,
  'ann_ver:s'     => \$ann_ver,
  'json'          => \$json,
  'h|help'        => \$help,
);

die "$usage\n" if ($help);

if ($json){
  use JSON; 
}

my %gn_to_pan_hsh;
my %pan_to_gn_hsh;

if ($hash_method =~ /BerkeleyDB/){
  tie %gn_to_pan_hsh, 'DB_File', $gn_to_pan, O_RDWR|O_CREAT, 0666, $DB_HASH
      or die "Cannot open gn_to_pan $gn_to_pan: $!\n";
  
  tie %pan_to_gn_hsh, 'DB_File', $pan_to_gn, O_RDWR|O_CREAT, 0666, $DB_HASH
      or die "Cannot open pan_to_gn $pan_to_gn: $!\n";
}
elsif ($hash_method =~ /tsv/){
  open my $G2P_FH, "<", $gn_to_pan or die "Can't open gn_to_pan $gn_to_pan:$!";
  while (<$G2P_FH>){
    chomp;
    my ($gn, $pan) = split("\t");
    $gn_to_pan_hsh{$gn} = $pan;
    $pan_to_gn_hsh{$pan} = $gn;
  }
}
else { die "Flag -hash_method must be either BerkeleyDB or tsv.\n" }

# Read gene list from stdin and look corresponding genes
my %gene_set_hsh;
while (<>){
  chomp;
  my $from_gn;
  if ($_ =~ /^[^.]+\.[^.]+\.[^.]+\.[^.]+/){
    $from_gn = $_;
  }
  else { # Infer that the gene isn't prefixed
    if (length($gensp) && length($strain) && length($gnm_ver) && length($ann_ver)){
      $from_gn = "$gensp.$strain.gnm$gnm_ver.ann$ann_ver.$_";
    }
    else {
      say "\nInput gene $_ lacks the four dot-spearated fields that specify assembly version";
      say "and at least one of the flags -gensp, -strain, -gnm_ver, -ann_ver was not set.\n";
      die "Please either provide a list of prefixed genes OR bare IDs with the required flags.\n\n";
    }
  }

  # Trim splice variant suffix if present
  $from_gn =~ s/\.\D*\d+$//;

  # For case-insensitive matching, lowercase IDs in the query list
  my $from_gn_lc = lc($from_gn);

  my $pan = $gn_to_pan_hsh{$from_gn_lc};

  # Put the pangene vector into an array
  my @to_gn_ary = split(/,/, $pan_to_gn_hsh{$pan});

  if ($to_annot =~ /all/i){
    say "$from_gn\t$pan\t", join ("\t", @to_gn_ary);
  }
  elsif ($to_annot =~ /^[^.]+\.[^.]+\.[^.]+\.[^.]+/){ # TO DO: Genralize. Assumes LIS DataStore prefix form.
    my @to_matches = grep(/$to_annot/ig, @to_gn_ary);
    if (scalar(@to_matches) == 0){ @to_matches = "NONE" }
    $gene_set_hsh{$from_gn} = join("\t", $pan, @to_matches);
  }
  else {
    die "Option -to_annot doesn't have the expected four dot-separated fields: [$to_annot]\n";
  }
}

# Report results
if ($json){
  say encode_json \%gene_set_hsh;
}
else {
  while ( my($from_gn, $pan_and_to_gns) = each(%gene_set_hsh)) {
    say "$from_gn\t$pan_and_to_gns";
  }
}

__END__

Versions
2023-10-29 Initial functioning version
2023-10-30 Add -json option
2023-11-01 Simplify handling of case for matched genes

#!/usr/bin/env perl

# PROGRAM: root_tree_by_species.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2017
# see description under Usage

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

####### File IO ########

#my $raxml_executable = "/home/scannon/.conda/envs/hybseq/bin/raxml-ng"; # c2s2 via conda
my $raxml_executable = "raxml-ng"; 
my ($align_in, $tree_in, $roots, $use_first_root, $full, $help, $verbose, $execute);
my $model="LG+G8+F";

GetOptions (
  "align_in=s" =>     \$align_in,   # required
  "tree_in:s" =>      \$tree_in,  
  "roots:s" =>        \$roots,   
  "use_first_root" => \$use_first_root,
  "full" =>           \$full,
  "verbose+" =>       \$verbose,
  "execute" =>        \$execute,
  "help" =>           \$help
);

my $usage = <<EOS;
  Usage: perl $0 [-options]
  
  Calculate a RAxML tree and root it by specified sequence patterns, if those sequences are present.
  Rooting is done by a sequence name or substring in the name, if the name is provided and present; 
  otherwise, the rooting uses the tree-height minimizing model of RAxML.

  For example, if there are outgroup species
    arath.AT1G77750.1
    arath.AT5G14320.1
    orysa.LOC_Os03g49710.1
  and the outgroup strings are supplied as follows:
    -roots "arath,orysa"
  then rooting will be done on arath and orysa ... unless -use_first_root is specified,
  in which case the first found of the list (e.g. arath) will be used as the single sequence for rooting.  
  If none is found, then the tree will be rooted using the RAxML model "-f I":
    "It roots the tree by rooting it at the branch that best balances the subtree lengths
     (sum over branches in the subtrees) of the left and right subtree."

  Flags and parameters:
   -align_in:   Input alignment (only used for identifying sequence IDs) **
   -tree_in:    User-provided tree (optional)
   -roots:      List of strings found in sequence names to be used as an outgroup,
                   comma-separated, within double quotes, e.g. -roots="orysa,arath"
                   If missing, then the tree will be rooted using the RAxML model "-f I"
   -use_first_root:  Use the first outgroup sequence encountered from the list to do the rooting
                     ... rather than all indicated outgroups together.
   -full:       Do raxml-ng "--all" option, i.e. all-in-one (ML search + bootstrapping); use 100 bootstraps
   -execute:    Execute the command
   -verbose:    For more stuff to terminal
   -help:       For more info
   
   ** = required
EOS

die "\n$usage\n" if ($help or !defined($align_in));

warn "Note: -execute wasn't set, so the command will be printed but not executed.\n" if (!defined($execute));

# read alignment in
open( my $ALIGN_IN, '<', $align_in ) or die "can't open align_in $align_in: $!";

# read tree if provided
my $TREE_IN;
if ( $tree_in ){
  open ( my $TREE_IN, '<', $tree_in ) or die "can't open tree_in $tree_in: $!";
}

my %seq_IDs;
while (my $line = <$ALIGN_IN>) {
  next unless ($line =~ /^>(\S+) *.*/);
  my $ID=$1;
  #print "$ID\n";
  $seq_IDs{$ID} = $ID;
}

my $raxml_command;
my $seed=int(100*rand());

my $search_mode;
if ($full){ $search_mode = "--all --tree pars{10} --bs-trees 100" } 
else { $search_mode = "--search1" }

if (defined($roots)) { # List of roots is provided
  # Put elements of LIST into array in case we need to retain list order
  my @outgrp_strs = split(/,/, $roots);
  
  # Given the ougroup root-strings (e.g. Arath, Oryza), find the corresponding sequence names 
  my %seen_seq_ID;
  my $first_outgrp="";
  my $set_of_out_IDs="";
  for my $root_str (@outgrp_strs) {
    #print "===== ", $root_str, " =====\n";
    for my $seq_ID (%seq_IDs) {
      #print "  $root_str $seq_ID\n";
      if ($seq_ID =~ /$root_str/) {
        if (length($first_outgrp)==0){ $first_outgrp = $seq_ID }
        if ($seen_seq_ID{$seq_ID}) { next }
        else {
          #print "$seq_ID\n";
          $set_of_out_IDs .= "$seq_ID,";
          $seen_seq_ID{$seq_ID}++;
        }
      }
    }
  }

  if (length($first_outgrp) > 0) { # Do outgroup rooting
    my $outgrp_seq_IDs;
    if ($use_first_root){
      $outgrp_seq_IDs = $first_outgrp;
    }
    else {
      $outgrp_seq_IDs = $set_of_out_IDs;
    }
    print "OUT: $first_outgrp\n";
    if ($verbose) { print "1: tree calc with outgroup rooting\n" }
    if ($tree_in){
      my $exec_str = \
        "$raxml_executable $search_mode --model $model --outgroup $outgrp_seq_IDs --msa $align_in --tree $tree_in --redo";
      if ( $execute ){ system($exec_str) }
      else { print "$exec_str\n" } 
    }
    else { # no starting tree was provided
      my $exec_str = "$raxml_executable $search_mode --model $model --outgroup $outgrp_seq_IDs --msa $align_in --redo";
      if ( $execute ){ system($exec_str) }
      else { print "$exec_str\n" } 
    }
  }
  else {
    if ($verbose) { print "2: tree calc; no outgroup rooting \n" }
    if ($tree_in){
      my $exec_str = "$raxml_executable $search_mode --model $model --tree $tree_in --msa $align_in --redo";
      if ( $execute ){ system($exec_str) }
      else { print "$exec_str\n" } 
    }
    else { # no starting tree was provided
      my $exec_str = "$raxml_executable $search_mode --model $model --msa $align_in --redo";
      if ( $execute ){ system($exec_str) }
      else { print "$exec_str\n" } 
    }
  }
} 
else { # No list of roots provided, so do tree calc without rooting
  if ($verbose) { print "3: tree calc without rooting\n" }
  if ($tree_in){
    my $exec_str = "$raxml_executable $search_mode --model $model --msa $align_in --tree $tree_in --redo";
    if ( $execute ){ system($exec_str) }
    else { print "$exec_str\n" } 
  }
  else { # no starting tree was provided
    my $exec_str = "$raxml_executable $search_mode --model $model --msa $align_in --redo";
    if ( $execute ){ system($exec_str) }
    else { print "$exec_str\n" } 
  }
}


__END__
VERSIONS
SC = Steven Cannon
v0.01 17-03-28 SC Initial version
v0.02 17-07-20 SC Add warning if -execute isn't set. 
                    Don't print outgroup $seq_ID again if already seen. 
                    Add -threads
                    Implement midpoint rooting, if outgroup list isn't provided or no outgroup is found in tree
v0.03 17-07-20 SC Remove option for taking in a user tree. 
v0.04 18-07-31 SC Fix midpoint rooting by providing a tree to root
v0.05 21-07-06 SC Update for raxml-ng
v0.06 21-09-19 SC Make input tree optional
v0.07 22-08-28 SC Add options to root on just one outgroup sequence, and to do raxml-ng --all

root_tree_by_species.pl -tree 09_trees/59026833_a -align 08_hmmalign_trimmed/59026833_a -out $PWD/09_trees_rt_by_outgrp -roots "ambtr,zeama,orysa,vitvi,solly,arath,prupe" -v -e


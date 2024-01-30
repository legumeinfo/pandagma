#!/bin/bash
#SBATCH   # Your ...
#SBATCH   # sbatch ...
#SBATCH   # commands ...
#SBATCH   # here ...
#SBATCH   # ...

set -o errexit
set -o nounset
set -o xtrace

date   # print timestamp

# If using conda environment for dependencies:
conda activate pandagma

PDGPATH=/path/to/your/pandagma
CONFIG=$PDGPATH/config/family3_sup_ur.conf

echo "Config: $CONFIG"

export PATH=$PATH:$PDGPATH:$PDGPATH/bin
echo "PATH: $PATH"
printf "Testing path: hash_into_fasta_id.pl is "; 
  which hash_into_fasta_id.pl
echo

##########
# Run all main steps, assuming ./out_pandagma exists from previous run of pandagma fam
#pandagma fsup -c $CONFIG

# Below: run indicated steps.

##########
# Run individual steps; equivalent to pandagma fsup -c $CONFIG
#pandagma fsup -c $CONFIG -s ingest
#pandagma fsup -c $CONFIG -s fam_consen
#pandagma fsup -c $CONFIG -s search_families
#pandagma fsup -c $CONFIG -s tabularize
#pandagma fsup -c $CONFIG -s realign_and_trim
#pandagma fsup -c $CONFIG -s calc_trees
pandagma fsup -c $CONFIG -s summarize

date   # print timestamp



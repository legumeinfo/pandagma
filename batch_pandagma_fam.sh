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
CONFIG=$PDGPATH/config/family3_18_3.conf 

echo "Config: $CONFIG"

export PATH=$PATH:$PDGPATH:$PDGPATH/bin
echo "PATH: $PATH"
printf "Testing path: hash_into_fasta_id.pl is "; 
  which hash_into_fasta_id.pl
echo

##########
# Run all main steps
#pandagma fam -c $CONFIG

# Below: run indicated steps.

##########
# Run individual steps; equivalent to pandagma fam -c $CONFIG
pandagma fam -c $CONFIG -s ingest
pandagma fam -c $CONFIG -s mmseqs
pandagma fam -c $CONFIG -s filter
pandagma fam -c $CONFIG -s dagchainer
pandagma fam -c $CONFIG -s ks_calc
pandagma fam -c $CONFIG -s ks_filter
pandagma fam -c $CONFIG -s mcl
pandagma fam -c $CONFIG -s consense
pandagma fam -c $CONFIG -s cluster_rest
pandagma fam -c $CONFIG -s add_extra
pandagma fam -c $CONFIG -s tabularize
pandagma fam -c $CONFIG -s summarize

##########
# Optional alignment, modeling, and tree-calculation steps
#pandagma fam -c $CONFIG -s align
#pandagma fam -c $CONFIG -s model_and_trim
#pandagma fam -c $CONFIG -s calc_trees

##########
# Optional work-directory cleanup steps
#pandagma fam -c $CONFIG -s clean
#rm -rf ./pandagma_work

date   # print timestamp



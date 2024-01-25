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

PDGPATH=$PWD
CONFIG=$PWD/config/XXX.conf

echo "Config: $CONFIG"

export PATH=$PATH:$PDGPATH:$PDGPATH/bin
echo "PATH: $PATH"

##########
# Test PATH
which pandagma
which calc_ks_from_dag.pl

##########
## Run all main steps
#pandagma pan -c $CONFIG

##########
## Run individual steps
#pandagma pan -c $CONFIG -s ingest
#pandagma pan -c $CONFIG -s mmseqs
#pandagma pan -c $CONFIG -s filter
#pandagma pan -c $CONFIG -s dagchainer
#pandagma pan -c $CONFIG -s mcl
#pandagma pan -c $CONFIG -s consense
#pandagma pan -c $CONFIG -s cluster_rest
#pandagma pan -c $CONFIG -s add_extra
#pandagma pan -c $CONFIG -s pick_exemplars
#pandagma pan -c $CONFIG -s filter_to_pctile
#pandagma pan -c $CONFIG -s order_and_name
#pandagma pan -c $CONFIG -s calc_chr_pairs
#pandagma pan -c $CONFIG -s summarize

##########
# Optional alignment and tree-construction steps
#pandagma pan -c $CONFIG -s align
#pandagma pan -c $CONFIG -s model_and_trim
#pandagma pan -c $CONFIG -s calc_trees
pandagma pan -c $CONFIG -s xfr_aligns_trees

##########
## Optional work-directory cleanup steps
#pandagma pan -c $CONFIG -s clean
#rm -rf pandagma_work

date   # print timestamp


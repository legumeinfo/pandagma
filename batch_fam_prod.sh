#!/bin/bash
#SBATCH -A m4440
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -n 30    #  number of cores/tasks in this job
#SBATCH -t 23:00:00
#SBATCH -C cpu
#SBATCH -J pand-fam2
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -o errexit
set -o nounset
set -o xtrace

date   # print timestamp

# If using conda environment for dependencies:
module load conda
conda activate pandagma

PDGPATH=$PWD
CONFIG=$PWD/config/family3_22_3.conf

echo "Config: $CONFIG"

export PATH=$PATH:$PDGPATH/bin
echo "PATH: $PATH"

##########
# Test PATH
which pandagma
which calc_ks_from_dag.pl

##########
## Fetch relevant data files; e.g.
#mkdir -p data
#make -C data -f $PWD/get_data/family3_22_3.mk

##########
## Filter transposable elements
#pandagma TEfilter -c $CONFIG

##########
## Run all main steps, assuming input data files exist in ./data
## Work directory will be ./work_pandagma
## Output will go to ./out_pandagma
#pandagma fam -c $CONFIG -d data_TEfilter

##########
## Run individual steps
#pandagma fam -c $CONFIG -s ingest -d data_TEfilter
#pandagma fam -c $CONFIG -s mmseqs
#pandagma fam -c $CONFIG -s filter
#pandagma fam -c $CONFIG -s dagchainer
#pandagma fam -c $CONFIG -s ks_calc
#pandagma fam -c $CONFIG -s ks_filter
#pandagma fam -c $CONFIG -s mcl
#pandagma fam -c $CONFIG -s consense
#pandagma fam -c $CONFIG -s cluster_rest
#pandagma fam -c $CONFIG -s add_extra
#pandagma fam -c $CONFIG -s tabularize
#pandagma fam -c $CONFIG -s align_protein
#pandagma fam -c $CONFIG -s model_and_trim
#pandagma fam -c $CONFIG -s calc_trees
pandagma fam -c $CONFIG -s xfr_aligns_trees
pandagma fam -c $CONFIG -s summarize

##########
## Optional work-directory cleanup steps
#pandagma fam -c $CONFIG -s clean
#rm -rf ./work_pandagma

date   # print timestamp


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

CONFIG=$PWD/config/XXX.conf
IMAGE=XXX_PATH_TO_YOUR_SINGULARITY_IMAGE.sif

export SINGULARITY_CLEANENV=TRUE

##########
## Fetch relevant data files; e.g.
#mkdir -p data
#singularity exec $IMAGE make -C data -f /usr/local/pandagma/get_data/family3_4_3.mk

##########
## Run all main steps, assuming input data files exist in ./data
## Work directory will be ./work_pandagma
## Output will go to ./out_pandagma
singularity exec $IMAGE pandagma fam -c $CONFIG

##########
## Run individual steps
#singularity exec $IMAGE pandagma fam -c $CONFIG -s ingest
#singularity exec $IMAGE pandagma fam -c $CONFIG -s mmseqs
#singularity exec $IMAGE pandagma fam -c $CONFIG -s filter
#singularity exec $IMAGE pandagma fam -c $CONFIG -s dagchainer
#singularity exec $IMAGE pandagma fam -c $CONFIG -s ks_calc
#singularity exec $IMAGE pandagma fam -c $CONFIG -s ks_filter
#singularity exec $IMAGE pandagma fam -c $CONFIG -s mcl
#singularity exec $IMAGE pandagma fam -c $CONFIG -s consense
#singularity exec $IMAGE pandagma fam -c $CONFIG -s cluster_rest
#singularity exec $IMAGE pandagma fam -c $CONFIG -s add_extra
#singularity exec $IMAGE pandagma fam -c $CONFIG -s tabularize
#singularity exec $IMAGE pandagma fam -c $CONFIG -s align_protein
#singularity exec $IMAGE pandagma fam -c $CONFIG -s model_and_trim
#singularity exec $IMAGE pandagma fam -c $CONFIG -s calc_trees
#singularity exec $IMAGE pandagma fam -c $CONFIG -s summarize

##########
## Optional work-directory cleanup steps
#singularity exec $IMAGE pandagma fam -c $CONFIG -s clean
#rm -rf ./work_pandagma

date   # print timestamp


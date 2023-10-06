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

ml singularityCE

##########
# Test PATH
SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE which pandagma-pan.sh
SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE which calc_ks_from_dag.pl

##########
## Run all main steps
SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG

##########
## Run individual steps
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s ingest
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s mmseqs
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s filter
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s dagchainer
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s ks_calc
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s ks_filter
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s mcl
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s consense
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s cluster_rest
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-pan.sh -c $CONFIG -s add_extra
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-fam.sh -c $CONFIG -s align
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-fam.sh -c $CONFIG -s model_and_trim
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-fam.sh -c $CONFIG -s calc_trees
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-fam.sh -c $CONFIG -s summarize

##########
## Optional work-directory cleanup steps
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-fam.sh -c $CONFIG -s clean
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma-fam.sh -c $CONFIG -s ReallyClean

date   # print timestamp


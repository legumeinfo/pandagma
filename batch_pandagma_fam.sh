#!/bin/bash
#SBATCH --job-name="pandagma_fam"     #  name of the job submitted
#SBATCH -p short                  #  name of the queue you are submitting job to
#SBATCH -N 1     #  number of nodes in this job
#SBATCH -n 38    #  number of cores/tasks in this job; 38 is about half a node
#SBATCH -t 20:00:00     #  time allocated for this job hours:mins:seconds
#SBATCH --mail-user=steven.cannon@usda.gov     # enter your email address to receive emails
#SBATCH --mail-type=END,FAIL     # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"     #  standard out %j adds job number to outputfile
## #SBATCH -q cicgru-mem      #  name of the queue you are submitting job to  QOS: cicgru-mem
## #SBATCH --mem=350G
## #SBATCH -e "stderr.%j.%N"     # optional but it prints out standard error

set -o errexit
set -o nounset
set -o xtrace

PDGPATH=$1  # /project/legume_project/steven.cannon/pandagma/fam
IMAGE=$2    # /project/legume_project/singularity_images/pandagma-2023-05-15.sif
CONFIG=$3   # config/family_7_0.conf

date   # print timestamp

ml singularity

# If using conda environment of this name:
##source activate pandagma

##singularity run $IMAGE -c $CONFIG
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE which pandagma_fam.sh
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec --cleanenv $IMAGE pandagma_fam.sh -c $CONFIG

# Run indicated steps
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec $IMAGE pandagma_fam.sh -c $CONFIG -s mmseqs
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec $IMAGE pandagma_fam.sh -c $CONFIG -s filter
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec $IMAGE pandagma_fam.sh -c $CONFIG -s dagchainer
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec $IMAGE pandagma_fam.sh -c $CONFIG -s mcl
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec $IMAGE pandagma_fam.sh -c $CONFIG -s consense
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec $IMAGE pandagma_fam.sh -c $CONFIG -s cluster_rest
#SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec $IMAGE pandagma_fam.sh -c $CONFIG -s add_extra
SINGULARITYENV_PREPEND_PATH=$PDGPATH:$PDGPATH/bin singularity exec $IMAGE pandagma_fam.sh -c $CONFIG -s align_and_trim

date   # print timestamp

# On Atlas:
##SBATCH -A legume_project
##SBATCH -p atlas
##SBATCH --mem=0  
##SBATCH --exclusive


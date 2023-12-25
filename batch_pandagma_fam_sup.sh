#!/bin/bash
#SBATCH --time=23:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=30   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="panfsup"
#SBATCH --mail-user=steven.cannon@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## #SBATCH --qos=legume-mem
## #SBATCH -q cicgru-mem      #  name of the queue you are submitting job to  QOS: cicgru-mem
## #SBATCH --mem=350G
## #SBATCH -e "stderr.%j.%N"     # optional but it prints out standard error

set -o errexit
set -o nounset
set -o xtrace

date   # print timestamp

# If using conda environment for dependencies:
ml miniconda
source activate pandagma

PDGPATH=/project/legume_project/steven.cannon/pandagma/fam
CONFIG=/project/legume_project/steven.cannon/pandagma/fam/config/family3_sup_ur.conf

echo "Config: $CONFIG"

export PATH=$PATH:$PDGPATH:$PDGPATH/bin
echo "PATH: $PATH"
printf "Testing path: hash_into_fasta_id.pl is "; 
  which hash_into_fasta_id.pl
echo

##########
# Run all main steps
#pandagma-fsup.sh -c $CONFIG

# Below: run indicated steps.

##########
# Run individual steps; equivalent to pandagma-fsup.sh -c $CONFIG
#pandagma-fsup.sh -c $CONFIG -s ingest
#pandagma-fsup.sh -c $CONFIG -s fam_consen
#pandagma-fsup.sh -c $CONFIG -s search_families
#pandagma-fsup.sh -c $CONFIG -s tabularize
#pandagma-fsup.sh -c $CONFIG -s realign_and_trim
#pandagma-fsup.sh -c $CONFIG -s calc_trees
pandagma-fsup.sh -c $CONFIG -s summarize

date   # print timestamp



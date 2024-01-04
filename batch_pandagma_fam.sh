#!/bin/bash
#SBATCH --time=23:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=30   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="panfam"
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

PDGPATH=/project/legume_project/steven.cannon/pandagma/pandagma
CONFIG=/project/legume_project/steven.cannon/pandagma/pandagma/config/family3_18_3.conf 

echo "Config: $CONFIG"

export PATH=$PATH:$PDGPATH:$PDGPATH/bin
echo "PATH: $PATH"
printf "Testing path: hash_into_fasta_id.pl is "; 
  which hash_into_fasta_id.pl
echo

##########
# Run all main steps
#pandagma-fam.sh -c $CONFIG

# Below: run indicated steps.

##########
# Run individual steps; equivalent to pandagma-fam.sh -c $CONFIG
pandagma-fam.sh -c $CONFIG -s ingest
pandagma-fam.sh -c $CONFIG -s mmseqs
pandagma-fam.sh -c $CONFIG -s filter
pandagma-fam.sh -c $CONFIG -s dagchainer
pandagma-fam.sh -c $CONFIG -s ks_calc
pandagma-fam.sh -c $CONFIG -s ks_filter
pandagma-fam.sh -c $CONFIG -s mcl
pandagma-fam.sh -c $CONFIG -s consense
pandagma-fam.sh -c $CONFIG -s cluster_rest
pandagma-fam.sh -c $CONFIG -s add_extra
pandagma-fam.sh -c $CONFIG -s tabularize
pandagma-fam.sh -c $CONFIG -s summarize

##########
# Optional alignment, modeling, and tree-calculation steps
#pandagma-fam.sh -c $CONFIG -s align
#pandagma-fam.sh -c $CONFIG -s model_and_trim
#pandagma-fam.sh -c $CONFIG -s calc_trees

##########
# Optional work-directory cleanup steps
#pandagma-fam.sh -c $CONFIG -s clean
#pandagma-fam.sh -c $CONFIG -s ReallyClean

date   # print timestamp



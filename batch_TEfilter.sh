#!/bin/bash
#SBATCH --time=23:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=30   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="TEfilter"
#SBATCH --mail-user=steven.cannon@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -o errexit
set -o nounset
set -o xtrace

date   # print timestamp

# If using conda environment for dependencies:
#conda activate pandagma
ml miniconda
source activate pandagma

PDGPATH=$PWD
CONFIG=$PWD/config/Glycine_7_3_2.conf
#CONFIG=$PWD/config/Glycine5_21_29_6.conf
#CONFIG=$PWD/config/family_ur.conf

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
pandagma TEfilter -c $CONFIG -s TEfilter

date   # print timestamp


#!/bin/bash
#SBATCH --time=23:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=30   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="pandagma-pan"
#SBATCH --mail-user=XXX  # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -o errexit
set -o nounset
#set -o xtrace

date   # print timestamp

module load miniconda
source activate pandagma

PDGPATH=$PWD
CONFIG=$PWD/config/Zea_7_2_0.conf

#########
# Fetch relevant data files; e.g.
mkdir -p data
get_data/get_Zea.sh

echo "Config: $CONFIG"

export PATH=$PATH:$PDGPATH:$PDGPATH/bin
echo "PATH: $PATH"

##########
# Test PATH
which pandagma
which calc_ks_from_dag.pl


#########
# Run all main steps
pandagma pan -c $CONFIG

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
#pandagma pan -c $CONFIG -s align_cds
#pandagma pan -c $CONFIG -s align_protein
#pandagma pan -c $CONFIG -s model_and_trim
#pandagma pan -c $CONFIG -s calc_trees
#pandagma pan -c $CONFIG -s xfr_aligns_trees

##########
## Optional work-directory cleanup steps
#pandagma pan -c $CONFIG -s clean
#rm -rf work_pandagma

date   # print timestamp


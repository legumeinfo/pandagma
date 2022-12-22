
pandagma.sh version

# NOTE 1: All of the following commands depend on parameters being set in
# pandagma.conf (or equivalent file, pointed to as below), before continuing.
# Example:  export PANDAGMA_CONF=conf/Zea.conf

# Note 2: All of the commands below will be run in order if the program is called
# without specifying any subcommand, e.g.
#   nohup pandagma.sh run &

# get info from matching GFF/BED and FNA files
pandagma.sh run ingest

# do mmseqs on all pairings of annotation sets
pandagma.sh run mmseqs

# filter based on expected chromosome pairings
pandagma.sh run filter

# Identify syntenic blocks
pandagma.sh run dagchainer

# Calculate clusters using Markov clustering
pandagma.sh run mcl

# Calculate a consensus sequence for each pan-gene set, using vsearch.
pandagma.sh run consense

# Add extra annotation sets, by homology to prior clusters
pandagma.sh run add_extra

# Calculate orthogroup composition and filter fasta files to core orthogroups.
pandagma.sh run filter_to_core

# Assign pan-gene names with consensus chromosome and position information
pandagma.sh run name_pangenes

# Calculate observed chromosome pairs
pandagma.sh run calc_chr_pairs

# Move results into output directory, and report some summary statistics
pandagma.sh run summarize



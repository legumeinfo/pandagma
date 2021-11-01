
pandagma.sh version

# Set initial default values
pandagma.sh init

# NOTE: modify pandagma.conf before continuing

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

# Assign pan-gene names with consensus chromosome and position information
pandagma.sh run name_pangenes

# Calculate observed chromosome pairs
pandagma.sh run calc_chr_pairs

# Move results into output directory, and report some summary statistics
pandagma.sh run summarize



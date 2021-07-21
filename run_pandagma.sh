
pandagma.sh version

# Set initial default values
pandagma.sh init

# Override default values for some variables
#pandagma.sh config clust_iden 0.98
#pandagma.sh config mcl_inflation 2

# get info from matching GFF and FNA files
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

# Move results into output directory, and report some summary statistics
pandagma.sh run summarize



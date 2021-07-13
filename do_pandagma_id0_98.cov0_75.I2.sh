
time pandagma.sh version

# Set initial default values
time pandagma.sh init

# Override default values for some variables
#time pandagma.sh config clust_iden 0.98
#time pandagma.sh config mcl_inflation 2

# get info from matching GFF and FNA files
time pandagma.sh run ingest

# do mmseqs on all pairings of annotation sets
time pandagma.sh run mmseqs

# filter based on expected chromosome pairings
time pandagma.sh run filter

# Identify syntenic blocks
time pandagma.sh run dagchainer

# Calculate clusters using Markov clustering
time pandagma.sh run mcl

# Calculate a consensus sequence for each pan-gene set, using vsearch.
time pandagma.sh run consense

# Move results into output directory, and report some summary statistics
time pandagma.sh run summarize



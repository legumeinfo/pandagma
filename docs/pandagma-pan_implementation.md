## Implementation outline for the pangene workflow of Pandagma

The steps below indicate the processes carried out for the `pandagma pan` workflow, 
called for example as `pandagma pan -c CONFIG_FILE.conf`

See the README.md file for further details and examples, and `pandagma pan -h` for help, or
`pandagma pan -h` for more information - for example, about setting parameters in the config file.

### Step: ingest

In this step, the annotation files are prepared for analysis. The main task here is
to add positional information from GFF3 or 6-column BED to FASTA IDs.
A hash file is made that records the existing gene IDs and a new set of IDs that 
include positional information. The format of the IDs with positional information is:
```
  chromosome__geneID__start__end__orient
  vigan.Shumari.gnm1.chr01__vigan.Shumari.gnm1.ann1.Vigan.01G000100.01__44431__47510__+
```

The chromosome and gene prefixing information is used by convention within this analysis, 
but is not necessary as long as a regular expression is provided that enables identification
of the annotation name and the chromosome number.

### Step: mmseqs

In this step, the CDS (nucleotide) sequences each annotation set are compared against
each other annotation set. This effectively results in an all-by-all comparison of all genes
from all annotation sets (though for the upper triangle of comparisons -- i.e., 
avoiding calculating B x A if A x B has been calculated already).

This homology search step is handled by mmseqs, subprogram easy-cluster, on pairs of annotation sets.
The output format is tabular, similar to the m8 output from blastn.

The parameters for this run are:
```
    --min-seq-id "$clust_iden"   # set in config file; typically 0.95 for pangene CDS search
    -c "$clust_cov"       # set in the config file; typically 0.5 for pangene CDS search
    --cov-mode 0          # coverage value applies to both the query and target sequences.
    --cluster-reassign
```
The parameters "clust_iden" and "clust_cov" are set in the config file. 
Reasonable values for these parameters are `clust_iden=0.95` and `$clust_cov=0.5`. 
The argument `--cov-mode 0` means that the indicated coverage.
value applies to both the query and target sequences. If clust_cov=0.5, 
then at least 50% of the query sequence and 50% of the target sequence must align to be reported.
         
### Step: filter

In this step, the homology results from mmseqs are evaluated relative to the `expected_chr_matches`,
if that constraint data has been provided in the config file. This data filters the results, 
for example, to retain matches between Chr01 and Chr01 but not Chr01 and Chr02 (if the latter is not
expected in this data). If no `expected_chr_matches` vector is provided, then this filtering is not applied.

### Step: dagchainer

In this step, the filtered matches are further evaluated for synteny between the annotation sets.
This evaluation is handled with DAGchainer, using the following parameters:
```
    -M 50    # Maximum match score 
    -E 1e-5  # Maximum E-value
    -A 6     # Minimum number of aligned pairs in a synteny block
    -s       # Include tandem alignments
    -g       # Length of a gap in bp; calculated for each annotation set, as the average distance 
               between genes from both annotation sets; stored in the variable $ave_gene_gap
    -D       # maximum distance allowed between two matches in basepairs; calculated as 
               20*$ave_gene_gap
```

### Step: mcl

In this step, genes are clustered using Markov clustering (mcl) based on the gene pairs from the previous step.
The key parameter in this step is the inflation parameter, which is provided in the
configuration file. Based on evaluation of results for a range of inflation parameters,
we recommend `mcl_inflation="1.6"`.

### Step: consense

Step "consense" adds previously unclustered sequences into an "augmented" pan-gene set, by homology.
The implementation is as follows:

For each pan-gene set, retrieve sequences into a multifasta file, also prefixing the gene IDs with panID

Pick a representative seq. for each orthogroup - as a sequence with the median length for that orthogroup.

Identify which genes are absent from the clusters generated in the mcl step. These
may have avoided being clustered due to being outside of a synteny block.

Retrieve the non-clustered genes, and search them against genes already clustered. This search
is conducted using `mmseqs easy-search`, with parameters
```
    --search-type 3    # nucleotide
    --cov-mode 5       # short sequence needs to be at least x% of the other seq.
    -c "${clust_cov}"  # set in the config file; 0.5 is a reasonable value.
```

Then ...
  - Filter the mmseqs easy-search result by `$clust_iden` and report panID, qry_gene, sbj_gene.
  - Filter based on list of expected chromosome pairings if provided.
  - Make augmented cluster sets. 
  - Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)

### Step: cluster_rest

In this step, the remaining unclustered sequences are clustered, filtered on chromosome 
patterns defined in expected_chr_matches, then combined with clusters derived 
from synteny or restricted homology (12_syn_pan_aug.clust.tsv). Unclustered singletons
are skipped. 

This search is conducted using `mmseqs easy-cluster`, with parameters
```
    --min-seq-id "$clust_iden"   # set in config file; typically 0.95
    -c "$clust_cov"       # set in the config file; typically 0.5 
    --cov-mode 0          # coverage value applies to both the query and target sequences.
    --cluster-reassign
```

### Step: add_extra

Add extra annotation sets (if provided) to the augmented clusters, by homology.
For each pan-gene set, retrieve sequences into a multifasta file.
Search non-clustered genes against pan-gene consensus sequences

This search
  is conducted using `mmseqs easy-search`, with parameters
```
    --search-type 3    # nucleotide
    --cov-mode 5       # short sequence needs to be at least x% of the other seq.
    -c "${clust_cov}"  # set in the config file; typically 0.5
```

Reasonable values for these parameters are `clust_iden=0.95` and `$clust_cov=0.5`. 

Next ...
  - Filter on chromosome patterns defined in expected_chr_matches and by identity and top match.
  - Join panIDs to unclustered genes. This restrictive filtering requires that genes have a chromosome match.
  - Derive a cluster-format file from the hash of panIDs and extra genes.
  - Make augmented cluster sets. 
  - Reshape from hash to mcl output format (clustered IDs on one line)
  - For each pan-gene set, retrieve sequences into a multifasta file.
  - The directory and multifasta file are called 19_pan_aug_leftover_merged_cds
  - Merge files in 19_pan_aug_leftover_merged_cds, prefixing IDs with panID__

### Step: pick_exemplars

  - Pick representative (exemplar) sequence for each pan-gene set (protein and CDS)
  - Get all protein sequences corresponding with 18_syn_pan_aug_extra.clust.tsv
  - Get protein sequences into pan-gene sets, corresponding with 18_syn_pan_aug_extra.clust.tsv
  - The directory and multifasta file are called 19_pan_aug_leftover_merged_prot 
  - Get all protein sequences from files in 19_pan_aug_leftover_merged_prot
  - Pick a representative sequence for each pangene set - as a sequence with the median length for that set.
  - Do this for both the CDS and protein sets.

### Step: tabularize

  - Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv
  - Field for sorting 18_syn_pan_aug_extra.table.tsv is sort_col: 8
  - Sort, putting header row at top.

### Step: filter_to_pctile

  - Calculate matrix of gene counts per orthogroup and annotation set
  - Select orthogroups with genes from selected percentiles of max_annot_ct annotation sets

### Step: order_and_name

Two options are provided for ordering genes, specified in the configuration file with
order_method="reference" or order_method="alignment". Within pandagma-pan, the ordering is
handled (whether by reference or alignment method) with the order_and_name subcommand. 

For the “reference” method, the ordering of the pangenes is determined relative to the
specified "preferred annotation". A pangene that contains a gene from the preferred
annotation takes its chromosome ordinal position from a sorting of the genes in the
preferred annotation (or from the first gene in the case of local paralogs within a
pangene). A subsequent gap-filling step (described below) places the remaining pangenes,
when possible, relative to the ordered reference-based pangenes.

For the “alignment” method, the consensus (most frequent) chromosome is first identified
for the genes in a pangene set (relying on the provided annotation regular expression to
parse and identify chromosome numbers from annotation data; genes on scaffolds are set
aside for the purpose of placement determination). Then, the gene ordering from each
annotation set is encoded into a pseudo-protein “pseudo-chromosome” sequence, which is
aligned with the corresponding pseudo-chromosome from all other encoded annotations, using
standard multiple sequence alignment methods (we are using ‘famsa’ for this alignment).
Then the gene placements are then decoded and used to calculate a consensus order for all
genes and their pangenes.

In greater detail: The genes from each annotation (for all pangenes) have an order within
that annotation set. That ordering information, along with the panID, is hashed to a
unique synthetic five-residue peptide string. This peptide is used to encode and later
decode the following information: gene-id, annotation-name, chromosome,
order-in-annotation, start, end, orientation. Then, for each annotation, the peptide
strings are concatenated into a peptide-sequence “pseudo-chromosome.” A typical such
encoded “chromosomal” sequence tends to be about 10k residues. Then, each such
chromosome-representing sequence from all annotation sets are aligned in protein space, to
generate a protein multiple sequence alignment. Finally, each sequence is from these
alignments is decoded to recover the gene-id, annotation-name, chromosome,
order-in-annotation, start, end, orientation. The result is a table that provides a
provisional ordering of most of the pangenes.

Finally, for both the reference and alignment ordering methods, some pangenes may not have
a placement following the steps above – either because a gene from the “preferred
annotation” is missing from a given pangene, or because no consensus ordinal position
could be determined based on pseudo-chromosomal multiple sequence alignment. To place the
leftover genes (when possible), a script “order_gapfill.pl” is used. This script compares
the positions of genes in an unplaced pangene with the positions of the corresponding
annotations from genes that have a consensus ordinal position. More specifically: for each
pangene without placement, score every gene on the chromosome to which it belongs as being
"before" or "after" the target gene, with respect to each annotation set. This gives a
hash of genes on a chromosome, with values indicating "beforeness" and "afterness". The
consensus placement puts the pangene at an inflection point where genes scored
predominantly as “before” transition to predominantly “after”. 

This implementation is not a formal HMM, but is similar in concept, as the algorithm
evaluates sequential data (a vector of signed integers that represents the “before” states
and “after” states for all genes in a chromosome relative to an unplaced gene) and
identifies a state-transition point in a that vector. The gene with unknown position is
placed at that transition point. The gap-filling step is among the most time-consuming
part of the pangene workflow, so has been parallelized using Parallel::ForkManager.

The gapfilling process may still leave some genes unplaced. Those genes are reported as
such, without ordinal position.

### Step: calc_chr_pairs

Generate a report of observed chromosome pairs. This is reported in the file observed_chr_pairs.tsv.
This information can be used to identify gene associations between other chromosomes 
than those expected. This information can be used to identify if there are translocations in some 
of the annotation sets, or if there are nontrivial chromosome correspondences between species. 
This is most informative if expected_chr_matches is not set in the config.

### Step: align_cds

This is an optional step (not run as part of "all", but runnable as e.g. 
`pandagma fam -c CONFIG.conf -s align`).

  - Retrieve sequences for each family, preparatory to aligning them.
  - For each pan-gene set, retrieve sequences into a multifasta file.
  - Align the gene families, using `famsa` with default parameters (two threads,
    but running (NPROC/2) jobs at a time).

### Step: model_and_trim

This is an optional step (not run as part of "all", but runnable as e.g. 
`pandagma fam -c CONFIG.conf -s model_and_trim`).

  - For each alignment calculated in the alignment step, build an HMM using hmmbuild.
  - Realign the fasta sequences to the corresponding HMM.
  - Trim HMM alignments to match-states.
  - Filter alignments prior to tree calculation, using parameters
```
    filter_align.pl -in "$filepath" -out 23_hmmalign_trim2/"$file" \
                    -depth 3 # minimum number of characters for reporting a column.
                    -pct_depth 20 # At least 20% of annotations need to have a residue at an
                        alignment position in order to retain that position in the alignment.
                    -min_pct_aligned 20  # minimum percent of aligned characters after trimming to 
                        depth for retaining sequence.
```

### Step: calc_trees

This is an optional step (not run as part of "all", but runnable as e.g. 
`pandagma fam -c CONFIG.conf -s calc_trees`).

For each cleaned alignment calculated in the previous step, calculate a phylogenetic tree.
```
    fasttree -quiet "$filepath"
```

### Step: xfr_aligns_trees

This is an optional step (not run as part of "all", but runnable as e.g. 
`pandagma fam -c CONFIG.conf -s xfr_aligns_trees`).

Transfer the relevant directories of alignments, HMMs, cleaned alignments, and trees,
to the output directory.

### Step: summarize

Summarize: Copy results into output directory, and report some summary statistics.
This step creates an output directory (default `./out_pandagma`).

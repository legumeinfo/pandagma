## Implementation outline for the gene family workflow of Pandagma

The steps below indicate the processes carried out for the `pandagma fam` workflow, 
called for example as `pandagma fam -c CONFIG_FILE.conf`

See the README.md file for further details and examples, and `pandagma fam -h` for help, or
`pandagma fam -h` for more information - for example, about setting parameters in the config file.

### Step: ingest

In this step, the annotation files are prepared for analysis. The main task here is
to add positional information from GFF3 or 6-column BED to FASTA IDs.
A hash file is made that records the existing gene IDs and a new set of IDs that 
include positional information. The format of the IDs with positional information is:
```
  chromosome__geneID__start__end__orient
  
  e.g.
  Phaseolus.pan3.chr01__Phaseolus.pan3.chr01_000100_pan21140__10001__11000__+
  or 
  vicfa.Hedin2.gnm1.chr1S__vicfa.Hedin2.gnm1.ann1.1g000080.1__165357__167950__+
```

The chromosome and gene prefixing information is used by convention within this analysis, 
but is not necessary as long as a regular expression is provided that enables identification
of the annotation name and the chromosome number.

### Step: mmseqs

In this step, the protein sequences each annotation set are compared against
each other annotation set. This effectively results in an all-by-all comparison of all genes
from all annotation sets (though for the upper triangle of comparisons -- i.e., 
avoiding calculating B x A if A x B has been calculated already).

This homology search step is handled by mmseqs, subprogram easy-cluster, on pairs of annotation sets.
The output format is tabular, similar to the m8 output from blastn.

The parameters for this run are:
```
    --min-seq-id "$clust_iden"   # set in config file; typically 0.40 for protein search
    -c "$clust_cov"       # set in the config file; typically 0.40 for protein search 
    --cov-mode 0          # coverage value applies to both the query and target sequences.
    --cluster-reassign
    --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
```

Note the output format, which uses fields similar to the BLAST .m8 output, plus aligned sequences
that are used later for Ks calculation, if ks_calc is run.
         
### Step: filter

In this step, the homology results from mmseqs are evaluated relative to the
`expected_quotas`, if that constraint data has been provided in the config file. This file
sets the number of gene matches retained from a specified species. For example, `arath 2`
indicates that up to two _Arabidopsis thaliana_ matches should be retained relative to
queries from any gene, while `glyma 4` indicates that up to four _Glycine max_ matches
should be retained relative to queries from any gene (4 rather than 2 because the
papilionoid subfamily to which _Glycine_ belongs experienced a whole-genome duplication,
and there was a second whole-genome duplication within the _Glycine_ genus itself). Note
that these "expected quotas" do not necessarily constrain a gene family to be limited to
this number of genes from any given annotation, because the family (aka cluster or
component) can grow due to other matches, transitively; but quotas are useful to limit the
file sizes of the (potentially very large) homology results, and quotas can also help
limit the growth of (potentially very large) clusters during the initial clustering step.

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

### Step: ks_calc

In this step, Ks values are calculated based on gene paris from processed DAGChainer output, 
generating gene-pair and block-median Ks values for subsequent filtering. 
For this step, alignments from mmseqs are converted to berkeleydb key-value hashes, 
where the key is composed of `query subject` `subject query`, in lexical order -- 
since the output of DAGChainer has query and subject in lexical order.

This step is managed by the script calc_ks_from_dag.pl. Because of the large number of 
Ks calculations that need to be made, this step is slow -- though it has been substantially
speeded up by using the alignments produced by mmseqs in the preceding mmseqs step.

### Step: ks_filter

In this step, 
From optional provided ${WORK_DIR}/stats/ks_peaks.tsv file and from processed DAGChainer output,i 
(with gene-pair and block-median Ks values added by calc_ks_from_dag.pl), filter on 
block-median Ks values, and combine the homology data into a file with gene pairs to be clustered.

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
    --search-type 1    # protein
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
    --search-type 1    # protein
    --cov-mode 5       # short sequence needs to be at least x% of the other seq.
    -c "${clust_cov}"  # set in the config file; typically 0.5
```

Next ...
  - Filter on chromosome patterns defined in expected_chr_matches and by identity and top match.
  - Join panIDs to unclustered genes. This restrictive filtering requires that genes have a chromosome match.
  - Derive a cluster-format file from the hash of panIDs and extra genes.
  - Make augmented cluster sets. 
  - Reshape from hash to mcl output format (clustered IDs on one line)
  - For each pan-gene set, retrieve sequences into a multifasta file.
  - The directory and multifasta file are called 19_pan_aug_leftover_merged_cds
  - Merge files in 19_pan_aug_leftover_merged_prot, prefixing IDs with panID__

### Step: tabularize

  - Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv
  - Field for sorting 18_syn_pan_aug_extra.table.tsv is sort_col: 8
  - Sort, putting header row at top.

### Step: align_protein

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

# pandagma
Pandagma is a collection of tools for calculating pangene sets and gene families. 

Authors: Steven Cannon, Hyunoh Lee, Nathan Weeks, Joel Berendzen, 2020-2025.

### Overview

The Pandagma software package provides methods for efficiently and sensitively identifying
pangene and gene family sets for annotation sets from eukaryotic genomes, with methods for
handling polyploidy and for targeting family construction at specified taxonomic depths.
The pangene and gene family workflows are called from the main driver program,
respectively, with pandagma pan and pandagma fam. The pangene workflow is suited for
identifying sets of genes that correspond by synteny (sequence similarity and gene order)
across genomes (in particular, sequences and positions) at the level of a species or
across species within a genus. A filter for expected chromosome matches can be provided in
configuration to limit the merger of gene pairs; this is particularly important for
handling polyploid genomes (to limit the collapse of paralogs between homoeologous
chromosomes). The gene family workflow is suited for identifying gene families among more
distant taxonomic groupings, from taxonomic families to kingdoms. Expected duplication
quotas and synonymous-site changes (Ks) can be set via configuration to help restrict
mergers beyond a desired taxonomic depth. Additional workflows allow the pre-process
removal of sequences with homology to transposable elements, and the addition of
annotation sets to pangene sets or gene families previously calculated, permitting
maintenance and continuity of analyses over time.

There are two main workflows: `pandagma pan` for pangene sets, and `pandagma fam` for gene
families. The pangene workflow is designed to operate primarily on CDS sequences, at the
level of a species or genus. The family workflow is designed to operate on protein
sequences, at the level of an organismal family or order.

Tandem arrays and paralogs are both generally incorporated into clusters in both the -pan
and -fam workflows, but components that span unexpected chromosome pairs are disallowed if
an expected_chr_matches array is provided. For example, in Glycine, the chromosome
pairings are 01-01, 02-02 … 20-20 could be specified, and a known translocation (11-13)
could also be specified and allowed. The aspects of the process that affect the clustering
behavior are (1) the initial homology search of all against all genes, which should find
available paralogs (including tandem duplicates); (2) the filtering for synteny, mediated
by DAGChainer; (3) restriction to allowed chromosome pairs, if configured, for
pandagma-pan; (4) clustering (handled by MCL); (5) addition of leftover genes to initial
clusters based on homology – either with or without chromosome-pair restrictions, as
configured.

### Output

The pangene workflow produces the following files (described in the `MANIFEST_output_pan.yml` file):
  - Cluster files in two formats
  - CDS and protein fasta files, each containing a representative sequence from each cluster
  - A counts table, giving counts of genes from each annotation-set in each pangene
  - A table of genes in each family or pangene, in columns by accession
  - CDS and protein files for clusters with larger than a specified number of annotation 
    sets represented (default >=25th percentile), with identifiers that contain the cluster's 
    modal chromosome and an ordinal position, of the form `Glycine.pan5.chr01_000100_pan46447`
  - A file of summary statistics, reporting characteristics of each input file and of the pangene clusters.

The gene family workflow produces the following files (described in the `MANIFEST_output_fam.yml` file):
  - Cluster files in two formats, and from three steps in the clustering process (corresponding with
    the initial clusters incorporating designated main annotations and strict synteny constraints, 
    then clusters of the main annotations with leftover genes added by homology, 
    then clusters incorporating extra annotation sets if provided).
  - CDS and protein fasta files, each containing a representative sequence from each cluster
  - A counts table, giving counts of genes from each annotation-set in each pangene
  - A table of genes in each family or pangene, in columns by accession
  - CDS and protein files for clusters with larger than a specified number of annotation 
    sets represented (default >=25th percentile), with identifiers that contain the cluster's 
    modal chromosome and an ordinal position, of the form `Glycine.pan5.chr01_000100_pan46447`
  - A file of summary statistics, reporting characteristics of each input file and of the pangene clusters.
  - Directory of fasta files for all families corresponding with the main cluster result

Also, for both the pangene and family workflows, the following set of directories is produced
if the option alignment, modeling, and tree-calculation steps are run:
  - Directory of alignments (if align is run)
  - Directory of Hidden Markov Models of all alignments (if model_and_trim is run)
  - Directory of alignments of family files back to the HMMs (if model_and_trim is run)
  - Directory of alignments, trimmed of non-match states (if model_and_trim is run)
  - Directory of phylogenetic trees (if calc_trees is run)

### Outline of the pangene workflow

The pangene workflow (`pandagma pan`) is essentially as follows:
* Add positional information to the gene IDs
* Find gene homologies among pairs of annotation sets, using mmseqs2
* Filter by synteny (DAGChainer) -- and, optionally, by a list of allowed chromosome pairings.
* Cluster (mcl)
* Add back the genes that were excluded to this point 
* Add "extra" annotation sets by homology 
* Generate a tabular form of the pangenes
* Identify a representative sequence for each pangene
* Calculate consensus order of the pangenes, based on whole-chromosome multiple alignments of ordered panIDs
* Calculate and report statistics

Optionally, alignments and trees can also be calculated.

See additional implementation details at docs/docs/pandagma-pan_implementation.md 

### Outline of the gene family workflow

The gene family workflow (`pandagma fam`) is similar to the pangene workflow, but
operates on proteins rather than CDS, and uses an optional "quotas" file to determine the
initial expected gene duplication relationships between species, considering known or
suspected whole-genome duplication histories at the evolutionary depth of interest. The
workflow can also filter using synonymous-distance (Ks) parameters calculated from the
data and thresholds provided by the user.
* Add positional information to the gene IDs
* Find gene homologies among pairs of annotation sets, using mmseqs2
* Filter by synteny (DAGChainer)
* Filter by additional criteria (either or both can be used):
   - by a list of provided "quotas", with whole-genome duplication counts per species.
   - by Ks peak values, calculated by the ks_calc step
* Cluster (mcl)
* Add back the genes that were excluded to this point 
* Add "extra" annotation sets by homology 
* Generate a tabular form of the pangenes
* Align sequences in each family
* Build HMMs for each alignment
* Realign to the HMMs and trim to HMM match-states
* Calculate gene trees
* Calculate and report statistics

See additional implementation details at docs/docs/pandagma-fam_implementation.md 

### Additional, optional workflows and utilities

Several additional optional worflows are available:
* `pandagma TEfilter` - to compare CDS and protein data sets against a user-provided set of transposable elements 
    or other sequences to be filtered against. The intended use is to remove such sequences prior to calculating
    pangenes and gene families.
* `pandagma pupdate` - to compare and generate a mapping between two pangene sets.
* `pandagma fsup` - to "supplement" or add species/annotations into gene families calculated previously. 
    The selected new annotation sets are compared against results of a prior full run of pandagma-fam.sh, using one of
    two methods: align_sample to compare new proteomes against a sample of protein sequences from multifasta files, OR
    hmmsearch, to compare against a hmmpress-ed hmm database of compressed-and-indexed HMMs.
    Of the two methods, hmmsearch is generally more accurate than hmm_sample (aligns more sequences), but takes longer.

Additionally, a large number of utilities are available in the bin directory.
In general, each of these has full usage information, accessible with `program_name -h`

## Installation methods <a name="installation"></a>

First, clone this repository and change your working directory:

    git clone https://github.com/legumeinfo/pandagma.git
    cd pandagma

Next, proceed to follow one of the two supported installation paths (conda or Apptainer/SingularityCE).

### Installation method 1: installation of scripts and dependencies with conda

Create a conda environment called `pandagma` from the environment.yml in this repository: 

    conda env create

Set a shell variable (e.g., PANDAGMA_ROOT) that contains the path to the pandagma git repository, add the pandagma/bin directory to your PATH, and activate the conda environment:
```
    export PANDAGMA_ROOT=/path/to/repo/pandagma
    export PATH=$PANDAGMA_ROOT/bin:$PATH
    conda activate pandagma
      #or
    source activate pandagma
```


### Installation method 2: create a SingularityCE/Apptainer container image

To build a [SingularityCE](https://sylabs.io/singularity/) / [Apptainer](https://apptainer.org/) container image from the provided 
[definition file](https://apptainer.org/docs/user/latest/definition_files.html) (`pandagma.def`):
after cloning this repository:
```
  singularity build pandagma.sif singularity.def
```

To run pandamga using singularity, use `singularity exec --cleanenv pandagma.sif <command> [options]`, e.g.:
```
  singularity exec --cleanenv pandagma.sif pandagma <subcommand> -c CONFIG_FILE <options>
```

## Usage for the main `pandagma pan` <a name="usage_pan"></a>

```
Usage: 
       pandagma pan -c CONFIG_FILE [options]
   or
       pandagma pan -c CONFIG_FILE -s SUBCOMMAND [options]

  Required:
           -c (path to the config file)

  Options: -s (subcommand to run. If "all" or omitted, all steps will be run; otherwise, run specified step)
           -w (working directory, for temporary and intermediate files [default: './work_pandagma'].)
           -d (data directory, for annotation files [default: './data'; or set to data_TEfilter following 'pandagma TEfilter']
           -o OUTPUT_DIR (name for the output directory [default: './out_pandagma'].
                Applicable only to "all" and "summarize" steps.)
           -O (ordering method, for placing pan-genes. Options: 
                "reference" (default; uses preferred_annot to order, then gap-filling for missing panIDs.)
                "alignment" (uses whole-chromosome alignment of ordered panIDs from all annotations)
           -n (number of processors to use. Defaults to number of processors available to the parent process)
           -r (retain. Don't do subcommand "clean" after running "all".)
           -v (version)
           -h (help)
           -m (more information)

The bin/ directory should be added to your PATH, to make the scripts there accessible.

Primary coding and protein sequences (both fasta) and annotation (GFF3 or BED) files must be listed in the
config file, in the arrays cds_files, annotation_files, and protein_files. See example files.

Note that the annotation and CDS files need to be listed in CORRESPONDING ORDER in the config.

FASTA deflines are assumed to be formatted with the ID separated from any additional fields by a space:

    >id1 ...optional fields...

BED or GFF annotation files must have corresponding IDs specified in mRNA feature ID attributes.
CDS coordinates are derived from CDS features:

    Chr1	.	mRNA	1111	6666	.	.	.	ID=id1
    Chr1	.	CDS	2222	3333	.	.	.	Parent=id1
    Chr1	.	CDS	4444	5555	.	.	.	Parent=id1

BED files are assumed to contain a single feature for each primary coding sequence specified, 
where the coordinates represent the minimum start (0-based) and maximum end positions of the 
primary coding sequence in the reference:

  molecule, feature-start, feature-end, mRNA-ID, score(0), strand, gene-ID

Optionally, an expected_chr_matches array variable can be specified in GENUS.conf,
which provides anticipated chromosome pairings, e.g.

expected_chr_matches=(
  01 01
  02 02
  ...
# allows for translocation between 11 and 13
  11 13
# allows for translocation between 13 and 11
  13 11
)

These pairings are used in a regular expression to identify terminal portions of molecule IDs, e.g.
  glyma.Wm82.gnm2.Gm01  glyso.PI483463.gnm1.Gs01
  glyma.Wm82.gnm2.Gm13  glyso.W05.gnm1.Chr11
If an expected_chr_matches array variable is not defined, then no such filtering will be done.

At the end of the process, remaining genes will be added to initial clusters, based on homology.
Remaining genes may be those falling on unanchored scaffolds, or on chromosomes by not part of
synteny blocks and so not making it into the synteny-based clusters.
```
Subcommands for the **pangene** workflow, `pandagma pan`, in order they are usually run:
```
                all - All of the steps below, except for clean
                        (Or equivalently: omit the -s flag; "all" is default)
             ingest - Prepare the assembly and annotation files for analysis
             mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies
             filter - Filter the synteny results for chromosome pairings, returning gene pairs.
         dagchainer - Run DAGchainer to filter for syntenic blocks
                mcl - Derive clusters, with Markov clustering
           consense - Calculate a consensus sequences from each pangene set,
                      adding sequences missed in the first clustering round.
       cluster_rest - Retrieve unclustered sequences and cluster those that can be.
          add_extra - Add other gene model sets to the primary clusters. Useful for adding
                      annotation sets that may be of lower or uncertain quality.
         tabularize - Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv
     pick_exemplars - Pick representative sequence for each pangene
     filter_to_core - Calculate orthogroup composition and filter fasta files to core orthogroups.
     order_and_name - Assign pangene names with consensus chromosomes and ordinal positions.
     calc_chr_pairs - Report observed chromosome pairs; useful for preparing expected_chr_matches
          summarize - Move results into output directory, and report summary statistics.
```

Subcommands for the **gene family** workflow, `pandagma fam`, in order they are usually run:

```
  Run these first (if using the ks_peaks.tsv file; otherwise, run all main steps and 
  ks filtering will be done using parameters ks_block_wgd_cutoff and max_pair_ks)
                all - All of the steps below, except for ks_filter and clean
                        (Or equivalently: omit the -s flag; \"all\" is default).
             ingest - Prepare the assembly and annotation files for analysis.
             mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies.
             filter - Filter the synteny results for chromosome pairings, returning gene pairs.
         dagchainer - Run DAGchainer to filter for syntenic blocks.
            ks_calc - Calculation of Ks values on gene pairs from DAGchainer output.

  Evaluate the stats/ks_histplots.tsv and stats/ks_peaks_auto.tsv files and
  put ks_peaks.tsv into the \${WORK_DIR}/stats directory, then run the following commands:
          ks_filter - Filtering based on provided ks_peaks.tsv file (assumes prior ks_calc step)
                mcl - Derive clusters, with Markov clustering.
           consense - Calculate a consensus sequences from each pan-gene set,
                      adding sequences missed in the first clustering round.
       cluster_rest - Retrieve unclustered sequences and cluster those that can be.
          add_extra - Add other gene model sets to the primary clusters. Useful for adding
                      annotation sets that may be of lower or uncertain quality.
         tabularize - Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv
          summarize - Move results into output directory, and report summary statistics.

  If generating alignments, models, and trees, run the following steps:
              align - Align families.
     model_and_trim - Build HMMs and trim the alignments, preparatory to calculating trees.
         calc_trees - Calculate gene trees.

  For both pandagma-pan and pandagma-fam, run either of the following subcommands separately if you wish:
              clean - Clean (delete) files in the working directory that are not needed
                        for later addition of data using add_extra and subsequent run commands.
                        By default, "clean" is run as part of "all" unless the -r flag is set.
```

Variables in the config file for the **pangene workflow**, `pandagma pan`:

```
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.50]
        consen_iden - Minimum identity threshold for consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of "extra" annotations [80]
      mcl_inflation - Inflation parameter, for Markov clustering [default: 2]
        strict_synt - For clustering of the "main" annotations, use only syntenic pairs (1)
                        The alternative (0) is to use all homologous pairs that satisfy expected_chr_matches
      consen_prefix - Prefix to use in names for genomic ordered consensus IDs [Genus.pan1]
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g. 
                        "([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+" 
                          for four dot-separated fields, e.g. vigan.Shumari.gnm1.ann1
                        or "(\D+\d+\D+)\d+.+" for Zea assembly+annot string, e.g. Zm00032ab
    preferred_annot - String to match and select an annotation set, from a gene ID.
                        This is used for picking representative IDs+sequence from an orthogroup, when
                        this annotation is among those with the median length for the orthogroup.
                        Otherwise, one is selected at random from those with median length.
    min_align_count - Minimum number of sequences in a family to trigger alignments, modeling, and trees [4]
min_annots_in_align - Minimum number of distinct annotation groups in an alignment to retain it [2]
```

Variables in the config file for the **family workflow**, `pandagma fam`:
```
         clust_iden - Minimum identity threshold for mmseqs clustering [0.40]
          clust_cov - Minimum coverage for mmseqs clustering [0.40]
        consen_iden - Minimum identity threshold for consensus generation [0.30]
         extra_iden - Minimum identity threshold for mmseqs addition of "extra" annotations [0.30]
      mcl_inflation - Inflation parameter, for Markov clustering [1.6]
        strict_synt - For clustering of the "main" annotations, use only syntenic pairs [1]
                        The alternative (0) is to use all homologous pairs that satisfy expected_quotas
      ks_low_cutoff - For inferring Ks peak per species pair. Don't consider Ks block-median values less than this. [0.5]
       ks_hi_cutoff - For inferring Ks peak per species pair. Don't consider Ks block-median values greater than this. [2.0]
         ks_binsize - For calculating and displaying histograms. [0.05]
ks_block_wgd_cutoff - Fallback, if a ks_peaks.tsv file is not provided. [1.75]
        max_pair_ks - Fallback value for excluding gene pairs, if a ks_peaks.tsv file is not provided. [4.0]

      consen_prefix - Prefix to use in orthogroup names
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g.
                        \"([^.]+\.[^.]+)\..+\"
                          for two dot-separated fields, e.g. vigan.Shumari
                        or \"(\D+\d+\D+)\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    preferred_annot - String to match and select an annotation set, from a gene ID.
                        This is used for picking representative IDs+sequence from an orthogroup, when
                        this annotation is among those with the median length for the orthogroup.
                        Otherwise, one is selected at random from those with median length.
    min_align_count - Minimum number of sequences in a family to trigger alignments, modeling, and trees
min_annots_in_align - Minimum number of distinct annotation groups in an alignment to retain it
```

## Example run for the pangene workflow <a name="pangene_example"></a>

1. Download the annotation data from a remote source (CDS, protein, and GFF or BED), and transform if 
    needed. You can do this with a simple shell script that executes curl commands and then 
    applies some transformations. See the files in get_data/ for examples. There are \"get_data\" 
    shell scripts and Makefiles.


```
       mkdir data
       # using a conda environment
       make -C data -j 2 -f $PANDAGMA_ROOT/get_data/Glycine_7_3_2.mk
       # using SingularityCE/Apptainer
       singularity exec pandagma.sif make -C data -j 2 -f /usr/local/pandagma/get_data/Glycine_7_3_2.mk
```
2. Create a config file to provide program parameters and indicate sequence and coordinate files to be analyzed.
    The config file sets program parameters and constraints, and then lists annotation files and fasta files.

    The best way to start a config file for a new run is to copy and modify an existing config file in the repository.
    As of mid-2024, sample config files are provided for eight pangene runs (seven for legume genera and
    one for maize (Zea)). Also see parameter descriptions above for both workflows.

    In addition to the main parameters a vector of chromosome-match constraints
    can be provided in the `expected_chr_matches` vector. This indicates anticipated chromosome pairings
    (see description above, and examples in the pangene config files). Here are values for this parameter,
    for the example pangene configuration file Glycine_7_3_2.conf (the optional file naming convention
    here indicating 7 `main` annotations, 3 `extra_constr` annotations, and 2 `extra_free` annotations):

```
      expected_chr_matches=(
        01 01
        02 02
        03 03
        ...
        20 20
```

3. Prepare the annotation and fasta files, and check correspondence with listings in config file.
    The annotation and fasta files need to be listed in corresponding order.
    The nth listed GFF file corresponds to the nth listed FASTA file. Also note that 
    there are blocks for annotation_files and fasta_files, and for "extra" files of two types. 
    The "main" set should list annotation sets that you consider trustworthy and "central" in some 
    respect. The latter (extra) may be problematic in some way - for example, with questionable assemblies 
    or annotations, or perhaps annotations for other species in the genus. For the main set, both 
    homology and synteny information is used in the clustering. Annotations in the extra set are placed
    into clusters from the main set, using only homology information - and either with or without
    constraints on chromosome matches to the pangene clusters 
    (e.g., chr1 genes only going to chr1 pangene clusters). 
```
        annotation_files  -- "main" set, to establish clusters based on homology + synteny
        cds_files  -- "main" set, to establish clusters based on homology + synteny
        
        annotation_files_extra_constr -- optional extra annotations, constrained by chromosome match
        cds_files_extra_constr -- optional extra annotations, constrained by chromosome match
        
        annotation_files_extra_free -- optional extra annotations, not constrained by chromosome match
        cds_files_extra_free -- optional extra annotations, not constrained by chromosome match
```

4. Start the run. The examples below assume a run using `$PANDAGMA_ROOT/config/Glycine_7_3_2.conf`
       
   This workflow is best run in an HPC environment. If your environment uses job scheduling such as slurm, 
   then you will modify a batch submission script to submit and control the job. Examples are provided for
   calling the pandagma workflows using singularity and conda.

5. Examine the output, and adjust parameters and possibly the initial chromosome correspondences.
    Output will go into a directory specified by the `-o OUT_DIR` option (default "./out_pandagma").
    
    After a successful run, the main output files will be copied from the work directory 
    (default `./work_pandagma`) to the output directory (default: `./out`), along with
    a manifest file with the following descriptions:

```
      cat MANIFEST_output_pan.yml
      ---
      18_syn_pan_aug_extra.clust.tsv: Pan-gene sets, in cluster format - i.e., ID in
        first column, followed by tab-separated gene list.
      18_syn_pan_aug_extra.counts.tsv: Matrix of counts of genes per annotation set for each pan-gene set.
      18_syn_pan_aug_extra.hsh.tsv: Pan-gene sets, in a two-column hash format, with
        the set ID in the first column and genes in the second.
      18_syn_pan_aug_extra_complement.fna: Complement of genes in this pan-gene set;
        i.e. not clustered, presumed to be singletons.
      18_syn_pan_aug_extra_table.tsv.gz: Table of genes in each family or pangene, in columns by accession
      21_pan_fasta_clust_rep_cds.fna: CDS pan-gene sequence, inclusive (not filtered
        by minimum cluster size or annotation-set representation).
      21_pan_fasta_clust_rep_prot.faa: Protein pan-gene sequence, inclusive (not
        filtered by minimum cluster size or annotation-set representation).
      23_syn_pan_pctl25_posn_cds.fna: CDS pan-gene sequence, omitting pan-genes
        smaller than 25% of the mode, with derived pan-gene IDs corresponding with
        consensus chromosome and ordinal position.
      23_syn_pan_pctl25_posn_prot.faa: Protein pan-gene sequence, omitting pan-genes
        smaller than 25% of the mode, with derived pan-gene IDs corresponding with
        consensus chromosome and ordinal position.
      stats.txt: Descriptive statistics about program parameters, input
        sequences, and pan-gene products.
```

5. (continued; examining the output ...)
   The summary of the run is given in the file stats.txt. Look at the modal values
   in the histograms, the report of proportion of each assembly with matches, etc.
   One of the output files that may be of use in a subsequent run is observed_chr_pairs.tsv .
   This can indicate possible translocations among genomes in the input data. These values can be 
   used to constrain matches in a subsequent run, via an array variable expected_chr_matches
   defined in the pandagma config file.
   For example, in Zea, the values
   for one run are:
```
           1  1  39612   <== Main chromosome correspondences
           5  5  30156
           2  2  30074
           3  3  27673
           4  4  25903
           8  8  23763
           6  6  20986
           7  7  19743
           9  9  19090
          10 10  16841
           9 10    748  <== corresponding with a known translocation in accession Oh7B
           1  3    527
           1  2    515
```


## Example run for the family workflow <a name="family_example"></a>

1. Download the annotation data from a remote source (CDS, protein, and GFF or BED), and transform if 
   needed. You can do this with a simple shell script that executes curl commands and then applies some 
   transformations. See the files in get_data/ for examples. There are \"get_data\" scripts for Glycine, 
   Medicago, Phaseolus, Vigna, and Zea.
    
       mkdir data
       # using a conda environment
       make -C data -j 2 -f $PANDAGMA_ROOT/get_data/family_7_3.mk
       # using SingularityCE/Apptainer
       singularity exec pandagma.sif make -C data -j 2 -f /usr/local/pandagma/get_data/family_7_3.mk
    
2. Create a config file to provide program parameters and indicate sequence and coordinate files to be analyzed.
    As with the pangene workflow (above), the config file sets program parameters and constraints, 
    and then lists annotation files and fasta files.

    The best way to start a config file for a new run is to copy and modify an existing config file in the repository.
    As of mid-2025, sample config files are provided for six gene family runs, and for four family "supplement" (`pandagma fsup`) runs.
    These config files give typical parameters, and contain comments where needed about parameter choices.
    Also see parameter descriptions above for both workflows.

    In addition to the main parameters, an optional set of homology constraints can be set in
    a vector of `expected_quotas`, which sets the number of gene matches from a specified species 
    that should be considered from a query gene, during the initial homology search step. 
    Here are the values used in the example config file, `family_7_3.conf`, for seven legume species 
    and three outgroup species. Here, `arath 2` indicates that up to two _Arabidopsis thaliana_ 
    matches should be retained relative to queries from any gene; `glyma 4` indicates that up to four
    _Glycine max_ matches should be retained relative to queries from any gene (4 rather than 2 because 
    the papilionoid subfamily to which _Glycine_ belongs experienced a whole-genome duplication, and
    there was a second whole-genome duplication within the _Glycine_ genus itself). Note that 
    these "expected quotas" do not necessarily constrain a gene family to be limited to this 
    number of genes from any given annotation, because the family (aka cluster or component) can
    grow due to other matches, transitively; but quotas are useful to limit the file sizes of the
    (potentially very large) homology results, and quotas can also help limit the growth of 
    (potentially very large) clusters during the initial clustering step.

```
      expected_quotas=(
        arath 2
        glyma 4
        medtr 2
        phavu 2
        pissa 2
        prupe 1
        sento 2
        singl 2
        vigun 2
        vitvi 1
      )
```

3. Prepare the annotation and fasta files, and check correspondence with listings in config file.
    The annotation and fasta files need to be listed in corresponding order.
    The nth listed GFF file corresponds to the nth listed FASTA file.
    Also note that there are blocks for annotation_files and fasta_files, 
    and for annotation_files_extra and fasta_files_extra.
    The former, main set should list annotation sets that you consider trustworthy and "central" in 
    some respect. The latter (data_extra) may be problematic in some way, or may be from outgroup species.
    For the main set, both homology and synteny information are used in the clustering. 
    Annotations in the extra set are placed into clusters from the main set, using 
    only homology information. 

4. Start the run. The examples below assume a run using `config/family_7_3.conf`
       
   This workflow is best run in an HPC environment. If your environment uses job scheduling such as slurm, 
   then you will modify a batch submission script to submit and control the job. Examples are provided for 
   calling the pandagma workflows using singularity and conda.

5. The typical usage for the family workflow is to run all steps from `ingest` through `ks_calc`; 
   then evalute the Ks peaks and add a `ks_peaks.tsv` file to work directory; and then run the  
   remaining steps (see **8** below).

    An intermediate output file, `stats/ks_peaks_auto.tsv`, is written to the work directory
    This should be examined for biological plausibility, along with the other 
    Ks results (histograms) in the work_pandagma/stats subdirectory.
    The `ks_peaks_auto.tsv` file can be examined and used to create a file named `ks_peaks.tsv`
    with changes relative to `ks_peaks_auto.tsv` if necessary to reflect known or suspected WGD histories. 
    If stats/ks_peaks.tsv is not provided, then Ks filtering will be done using values provided
    in the config file for ks_block_wgd_cutoff and max_pair_ks.

6. Run steps `ks_filter` through `summarize`.
    The family workflow can be run straight through, without providing a `ks_peaks.tsv` file; 
    but the file does permit fine-grained control of Ks thresholds, appropriate for each species pair. 
    The typical usage for the family workflow is to run all steps from `ingest` through `ks_calc`; 
    then evalute the Ks peaks and add a `ks_peaks.tsv` file to the data_fam directory; 
    then run the steps `ks_filter` through `summarize`. In step 8 here, we run that last set of steps.

7. Examine the output, and adjust parameters and possibly the initial chromosome correspondences.
    After a successful run, the main output files will be copied from the work directory 
    (default: `./work_pandagma`) to the output directory (default: `./out_pandagma`).
    Output will go into a directory specified by the `-o OUT_DIR` option (default "./out_pandagma").

    The summary of the run is given in the file stats.txt .
    Look at the modal values in the histograms, the report of proportion of each assembly 
    with matches, etc.
    
    The output directory will include a manifest file, `MANIFEST_output_fam.yml`, with the following descriptions:
```
      cat MANIFEST_output_fam.yml
      ---
      06_syn_pan.clust.tsv: Gene family sets of "main" annotations only, clustered by homology and synteny, 
        in cluster format - i.e., ID in first column, followed by tab-separated gene list.
      06_syn_pan_ge3.hsh.tsv: Gene family sets of "main" annotations only, clustered by homology and synteny,
        in a two-column hash format, with the set ID in the first column and genes in the second. 
        Filtered to include three or more members.
      12_syn_pan_aug.clust.tsv: Gene family sets of "main" annotations only, clustered by homology and synteny
        and augmented by leftover genes into the initial clusters plus clusters of remaining genes, 
        in cluster format - i.e., ID in first column, followed by tab-separated gene list.
      12_syn_pan_aug.hsh.tsv: Gene family sets of "main" annotations only, clustered by homology and synteny
        and augmented by leftover genes into the initial clusters plus clusters of remaining genes, 
        in a two-column hash format, with the set ID in the first column and genes in the second.
      18_syn_pan_aug_extra.clust.tsv: Gene family sets, in cluster format - i.e., ID in
        first column, followed by tab-separated gene list.
      18_syn_pan_aug_extra.counts.tsv: Matrix of counts of genes per annotation set for each gene family.
      18_syn_pan_aug_extra.hsh.tsv: Gene family file, in a two-column hash format, with
        the set ID in the first column and genes in the second.
      18_syn_pan_aug_extra.table.tsv: Tabular format of all gene families, in columns by species. Multiple paralogs
        from a species in a gene family are comma-separated within the column of the species.
      19_pan_aug_leftover_merged: Directory of fasta files for all families 
        corresponding with 18_syn_pan_aug_extra.clust.tsv
      stats.txt: Descriptive statistics about program parameters, input
        sequences, and Gene family products.
```

8. Optionally (for both the pangene and family workflows), run steps `align`, `model_and_trim`, 
    and `calc_trees`.
    These steps align all families (or pangenes), then calculate HMMs of each alignment; 
    then realign to the HMMs and trim out the non-match-state characters, and then calculate 
    gene trees for each trimmed alignment. 
    The output will include corresponding directories, which will be copied to the 
    indicated output directory, 
    with the other results.


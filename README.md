# pandagma
Pandagma is a collection of tools for calculating pangene sets and gene families. 
There are two main workflows: pandagma-pan.sh for pangene sets, and pandagma-fam.sh for gene families.
The pangene workflow is designed to operate primarily on CDS sequences, at the level of a species or genus.
The family workflow is designed to operate on protein sequences, at the level of an organismal family or order.

Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2023.

The pangene workflow (pandagma-pan.sh) is essentially as follows:
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

Alignments and trees can also be calculated with these optional steps:
* align
* model_and_trim
* calc_trees
* xfr_aligns_trees


The gene family workflow (pandagma-fam.sh) is similar to the pangene workflow, but operates on proteins rather than CDS,
and uses an optional "quotas" file to determine the initial expected gene duplication relationships between species,
considering known or suspected whole-genome duplication histories at the evolutionary depth of interest.
* Add positional information to the gene IDs
* Find gene homologies among pairs of annotation sets, using mmseqs2
* Filter by synteny (DAGChainer) 
*   Additional filter options -- either or both can be used:
*     - by a list of provided "quotas", with whole-genome duplication counts per species.
*     - by Ks peak values, calculated by the ks_calc step
* Cluster (mcl)
* Add back the genes that were excluded to this point 
* Add "extra" annotation sets by homology 
* Generate a tabular form of the pangenes
* Align sequences in each family
* Build HMMs for each alignment
* Realign to the HMMs and trim to HMM match-states
* Calculate gene trees
* Calculate and report statistics

## Installation methods <a name="installation"></a>

### Installation method 1: by creating a Singularity container image

To build a [Singularity](https://singularity.hpcng.org/) container image from the provided [Singularity definition file](https://singularity.hpcng.org/user-docs/master/definition_files.html) (`singularity.def`), after cloning this repository:

    singularity build [--remote] pandagma.sif singularity.def

Where the `--remote` option is used if building the image as an unprivileged user using (by default) the [Sylabs Cloud Remote Builder](https://cloud.sylabs.io/builder).

To run pandamga using singularity, use `singularity run --cleanenv pandagma.sif [options]`, e.g.:

    singularity run --cleanenv pandagma.sif -c CONFIG_FILE

### Installation method 2: manual installation of scripts and dependencies

These dependencies are required: 
```  
  bioperl, bioperl-run, perl-parallel-forkmanager, perl-list-moreutils 
  mmseqs, dagchainer, mcl, EMBOSS, famsa, fasttree, hmmer
```
These need to be installed and available in the script's environment.

Once those are available on your PATH, the program can be called directly with its options.
(see Usage below).

Installing the dependencies is up to you. They can be installed via a suitable package manager. 
For example, using conda: 
~~~
  conda create -n pandagma
  conda install -n pandagma -c conda-forge -c bioconda perl-bioperl-core perl-parallel-forkmanager \
    perl-bioperl-run perl-list-moreutils dagchainer mcl mmseqs2 emboss famsa fasttree hmmer
  conda activate pandagama
~~~
Or on some systems, activate conda thus:
~~~
  source activate pandagma
~~~

## Usage for the main pandagma-pan.sh <a name="usage"></a>

~~~
Usage: 
       nohup ./pandagma-pan.sh -c CONFIG_FILE [options] &
   or
       nohup ./pandagma-pan.sh -c CONFIG_FILE -s SUBCOMMAND [options] &

   ... or using a suitable job submission script if using a workload manager such as slurm.

  Required:
           -c (path to the config file)

  Options: -s (subcommand to run. If \"all\" or omitted, all steps will be run; otherwise, run specified step)
           -w (working directory, for temporary and intermediate files.
                Must be specified in config file if not specified here.)
           -n (number of processors to use. Defaults to number of processors available to the parent process)
           -r (retain. Don't do subcommand \"clean\" after running \"all\".)
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

BED files are assumed to contain a single feature for each primary coding sequence specified, where the coordinates represent the minimum start (0-based) and maximum end positions of the the primary coding sequence in the reference:

  molecule, feature-start, feature-end, mRNA-ID, score(0), strand, gene-ID

Optionally, a file specified in the expected_chr_matches variable can be specified in GENUS.conf,
which provides anticipated chromosome pairings, e.g.
  01 01
  02 02
  ...
  11 13  # allows for translocation between 11 and 13
  13 11  # allows for translocation between 13 and 11
These pairings are used in a regular expression to identify terminal portions of molecule IDs, e.g.
  glyma.Wm82.gnm2.Gm01  glyso.PI483463.gnm1.Gs01
  glyma.Wm82.gnm2.Gm13  glyso.W05.gnm1.Chr11
If an expected_chr_matches file is not provided, then no such filtering will be done.

At the end of the process, remaining genes will be added to initial clusters, based on homology.
Remaining genes may be those falling on unanchored scaffolds, or on chromosomes by not part of
synteny blocks and so not making it into the synteny-based clusters.

Subcommands (in order they are usually run):
                all - All of the steps below, except for clean and ReallyClean
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
     calc_chr_pairs - Report observed chromosome pairs; useful for preparing expected_chr_matches.tsv
          summarize - Move results into output directory, and report summary statistics.

  Run either of the following subcommands separately if you wish:
              clean - Clean (delete) files in the working directory that are not needed
                        for later addition of data using add_extra and subsequent run commands.
                        By default, \"clean\" is run as part of \"all\" unless the -r flag is set.
        ReallyClean - Do complete clean-up of files in the working directory.
                        Use this if you want to start over, OR if you are satisified with the results and
                        don't anticipate adding other annotation sets to this pangene set.

Variables in pandagma config file:
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.50]
        consen_iden - Minimum identity threshold for consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of "extra" annotations [80]
      mcl_inflation - Inflation parameter, for Markov clustering [default: 2]
        strict_synt - For clustering of the \"main\" annotations, use only syntenic pairs (1)
                        The alternative (0) is to use all homologous pairs that satisfy expected_chr_matches.tsv
      consen_prefix - Prefix to use in names for genomic ordered consensus IDs [Genus.pan1]
       out_dir_base - Base name for the output directory [default: './out']
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g. 
                        "([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+" 
                          for four dot-separated fields, e.g. vigan.Shumari.gnm1.ann1
                        or "(\D+\d+\D+)\d+.+" for Zea assembly+annot string, e.g. Zm00032ab
    preferred_annot - String to match and select an annotation set, from a gene ID.
                        This is used for picking representative IDs+sequence from an orthogroup, when
                        this annotation is among those with the median length for the orthogroup.
                        Otherwise, one is selected at random from those with median length.
           work_dir - Working directory, for temporary and intermediate files. 

~~~

## Example run <a name="example"></a>

1. Clone the program and associated files from github:
      clone https://github.com/legumeinfo/pandagma
      I typically rename the downloaded repository to the genus name that I am working on, e.g.
    
            git clone https://github.com/legumeinfo/pandagma.git
            mv pandagma Zea
            cd Zea

2. Make a work directory. I typically make this in the same area as the pandagma directory, e.g.

            mkdir ../work_Zea
            ls ..
              work_Zea  Zea

3. Get into a suitable work environment (computation node), and load dependencies.
    The script has these third-party dependencies: **bioperl, mcl, mmseqs2, DAGchainer, famsa, EMBOSS**

    These can be loaded using a module-loading system, or with a package manager such as conda, or
    via a Singularity image. 
      
          pandagma_sing_img=$YOURPATH/pandagma.sif
    
    OR:
    
          conda activate pandagma    
            # pandagma is the name of the conda environment where I've installed the dependencies above.
    
4. Download the annotation data from a remote source (CDS, protein, and GFF or BED), and transform if needed.
    I do this with a simple shell script that executes curl commands and then applies some transformations.
    See the files in get_data/ for examples. There are \"get_data\" scripts for Glycine, Medicago, Phaseolus, Vigna, and Zea.
    
          ./get_data/get_Zea.sh
            # This puts the data first into data_orig/, 
            # and puts transformed data into data/
    
5. Create a config file to provide program parameters and indicate sequence and coordinate files to be analyzed.
    The config file sets nine program parameters, and then lists annotation files and fasta files.
    The annotation and fasta files need to be listed in corresponding order.
    The nth listed GFF file corresponds to the nth listed FASTA file.
    Also note that there are blocks for annotation_files and fasta_files, 
    and for annotation_files_extra and fasta_files_extra.
    The former, main set should list annotation sets that you consider trustworthy and "central" in some respect.
    The latter (data_extra) may be problematic in some way - for example, with questionable assemblies or annotations, 
    or perhaps annotations for other species in the genus. For the main set, both homology and synteny information
    is used in the clustering. Annotations in the extra set are placed into clusters from the main set, using 
    only homology information. Currently, the program is constructed to use at least one "extra" annotation.
    Future versions of the program may remove this requirement.

6. Start the run. The examples below assume a run using pandagma_Zea_7_2_nuc.conf.
       
    Calling the program directly:

      conda activate pandagma

      nohup ./pandagma-pan.sh -c config/Zea_7_2.conf &


    Using a Singularity image:

         pandagma_sing_img=$YOURPATH/pandagma-v_YOUR_VERSION
         nohup singularity exec --cleanenv $pandagma_sing_img ./pandagma-pan.sh -c config/Zea_7_2.conf &

7. Examine the output, and adjust parameters and possibly the initial chromosome correspondences.
    Output will go into a directory composed from a provided prefix name (default "out") and
    information about key parameter values, e.g.

          out_Zea_7_2/

    The summary of the run is given in the file stats.[parameters].txt .
    Look at the modal values in the histograms, the report of proportion of each assembly with matches, etc.
    One of the output files that may be of use in a subsequent run is observed_chr_pairs.tsv .
    This can indicate possible translocations among genomes in the input data. For example, in Zea, the values
    for one run are:

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
           9 10    748  <== corresponding with a known translocation in Oh7B
           1  3    527
           1  2    515

    These values can be used to constrain matches in a subsequent run, via the file data/expected_chr_matches.tsv
    (For the Zea run, these values were already provided, generated by get_data/get_Zea.sh 



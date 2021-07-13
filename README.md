# pandagma
Generate pan-gene sets, given a collection of genome assemblies and corresponding gene models.

The pandagma.sh script generates pan-gene clusters, given input consisting of genome assemblies and
gene models (GFFs) for a collection of accessions and/or closely related species.
The script drives several programs that do the primary work of clustering, identifying synteny, and 
additional searches and refinements. The other programs are: mmseqs, dagchainer, mcl, and vsearch. 
These need to be installed and available in the script's environment.
Authors: Steven Cannon, Joel Berendzen, 2020-2021

Create conda environment with required dependencies
~~~
  conda create -n pandagma
  conda install -n pandagma -c conda-forge -c bioconda perl-bioperl-core
  conda install -n pandagma -c conda-forge -c bioconda dagchainer
  conda install -n pandagma -c conda-forge -c bioconda mcl
  conda install -n pandagma -c conda-forge -c bioconda vsearch
  conda install -n pandagma -c conda-forge -c bioconda mmseqs2
  conda activate pandagama
~~~

The pandagma.sh has several "run commands," which are run in a particular order to do the clustering
and refinement. The several commands are best called via an additional calling script. 
The following commands should produce a reasonable collection of pan-genes, given suitably
named assembly and gene-model files:

~~~
pandagma.sh version

# Set initial default values
pandagma.sh init

# get info from matching GFF and FNA files
pandagma.sh run ingest

# do mmseqs on all genome pairings
pandagma.sh run mmseqs

# filter based on expected chromosome pairings
pandagma.sh run filter

# Identify syntenic blocks
pandagma.sh run dagchainer

# Calculate clusters using Markov clustering
pandagma.sh run mcl

# Calculate a consensus sequence for each pan-gene set, using vsearch.
pandagma.sh run consense

# Move results into output directory, and report some summary statistics
pandagma.sh run summarize
~~~

~~~
Usage: ${pkg} SUBCOMMAND [SUBCOMMAND_OPTIONS]

By default, name-matched assembly (fasta) and annotation (GFF) files are expected
in the data/ directory, within the working directory from where this script is called.
(To set a different data directory name, see discussion of environment variables below).
Example of name-matched files within the data/ directory:
  genome1.fna genome1.gff3
  genome2.fna genome2.gff3
  genomeXYZ.fna genomeXYZ.gff3

Optionally, a file \"expected_chr_matches.tsv\" can be provided (also in the data/ directory),
which provides anticipated chromosome pairings, e.g.
  01 01
  02 02
  ...
  11 13  # allows for translocation between 11 and 13
  13 11  # allows for translocation between 13 and 11
These pairings are used in a regular expression to identify terminal portions of molecule IDs, e.g.
  glyma.Wm82.gnm2.Gm01  glyso.PI483463.gnm1.Gs01
  glyma.Wm82.gnm2.Gm13  glyso.W05.gnm1.Chr11
If \"expected_chr_matches.tsv\" is not provided, then no such filtering will be done.

At the end of the process, remaining genes will be added to initial clusters, based on homology.
Remaining genes may be those falling on unanchored scaffolds, or on chromosomes by not part of
synteny blocks and so not making it into the synteny-based clusters.

Subommands (in order they are usually run):
            version - Get installed package version
               init - Initialize parameters required for run
             config - View/set run parameters, '-h' for help
         run ingest - Prepare the assembly and annotation files for analysis
         run mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies
         run filter - Filter the synteny results for chromosome pairings, returning gene pairs.
     run dagchainer - Run DAGchainer to filter for syntenic blocks
            run mcl - Derive clusters, with Markov clustering
       run consense - calculate a consensus sequences from each pan-gene set,
                       If possible add sequences missed in the first clustering round.
      run summarize - Move results into output directory, and report summary statistics.
       clear_config - Clear all config variables
              clean - Delete work directory

Variables (accessed by \"config\" subcommand):
      max_main_jobs - max number of significant jobs to run concurrently [default: processors/5]
      max_lite_jobs - max number of lightweight jobs to run concurrently [default: processors/2]
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.98]
          clust_cov - Minimum coverage for mmseqs clustering [0.75]
        consen_iden - Minimum identity threshold for vsearch consensus generation [0.80]
          fasta_ext - Extension of FASTA files
            gff_ext - Extension of GFF files
         pan_prefix - Prefix to use as a prefix for pangene clusters [default: pan]
       out_dir_base - base name for the output directory [default: './${pkg}_out']
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]
        mcl_threads - Threads to use in Markov clustering [default: processors/5]
            version - version of this script at config time

Environmental variables (may be set externally):
         ${PKG}_DIR - Location of the config directory, currently
                       \"${root_dir}\"
    ${PKG}_WORK_DIR - Location of working files, currently
                       \"${work_dir}\"
    ${PKG}_DATA_DIR - Location of intput data files (assemblies and GFs), currently
                       \"$PWD/data\"
              NPROC - Number of processors to use to set the max_main_jobs and max_lite_jobs
                       configuration variable. If not present, max_main_jobs is set
                       to a fifth of the number of processors on the system
                       and max_lite_jobs is set to half the number of processors.
~~~

 



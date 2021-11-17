# pandagma
Generate pan-gene sets, given coding sequences and corresponding gene models, from a set of related accessions or species.

Pandagma is a workflow that derives pan-gene sets, given input of sequences (CDS or peptide) and
annotations (BED or GFF) from several annotations. It is designed to be suitable for constructing
pan-genes for annotations across several different species within a genus, or annotations within
one species.

Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2021.

The workflow is essentially as follows:
* Add positional information to the gene IDs
* Find gene homologies among pairs of annotation sets, using mmseqs2
* Filter by synteny (DAGChainer) -- and, optionally, by a list of allowed chromosome pairings.
* Cluster (mcl)
* Add back the genes that were excluded to this point 
* Add "extra" annotation sets by homology 
* Calculate and report statistics


## Installation methods <a name="installation"></a>

### Installation method 1 (recommended): by creating a Singularity container image


To build a [Singularity](https://singularity.hpcng.org/) container image from the provided [Singularity definition file](https://singularity.hpcng.org/user-docs/master/definition_files.html) (`singularity.def`):

    singularity build [--remote] pandagma.sif singularity.def

Where the `--remote` option is used if building the image as an unprivileged user using (by default) the [Sylabs Cloud Remote Builder](https://cloud.sylabs.io/builder).

To run pandamga using singularity, use `singularity run pandagma.sif [run command]`, optionally setting environment variables, e.g.:

    singularity run pandagma.sif init
    ... modify pandagma.conf ...
    mkdir /path/to/pandagma/work_dir # create work dir on fast, node-local storage
    singularity run --env NPROC=<number of cpus> --env PANDAGMA_WORK_DIR=/path/to/work/dir pandagma.sif run

### Installation method 2: manual installation of scripts and dependencies

The scripts should be added to your path - either by copying them to a directory that is in your path 
(e.g. `cp scripts/* ~/bin/`), or by adding the script directory to your path 
(e.g. `export PATH=$PATH:~/src/pandagma/scripts` ... if that is where you have them).

Also, these dependencies are required: 
  bioperl, mmseqs, dagchainer, and mcl.
These need to be installed and available in the script's environment.

Installing the dependencies is up to you. They can be installed via a suitable package manager. 
For example, using conda: 
~~~
  conda create -n pandagma
  conda install -n pandagma -c conda-forge -c bioconda perl-bioperl-core
  conda install -n pandagma -c conda-forge -c bioconda dagchainer
  conda install -n pandagma -c conda-forge -c bioconda mcl
  conda install -n pandagma -c conda-forge -c bioconda mmseqs2
  conda activate pandagama
~~~

## Usage for the main pandagma.sh <a name="usage"></a>

~~~
Usage: 
       pandagma.sh run
   or
       pandagma.sh run SUBCOMMAND

Add the scripts directory to your PATH
Export these environment variables:
  PANDAGMA_WORK_DIR
  PANDAGMA_CONF
  NPROC

Primary coding sequence (fasta) and annotation (GFF3 or BED) files must be listed in the
fasta_files and annotation_files variables defined in pandagma.conf, which by default must exist
within the working directory from where this script is called.

FASTA deflines are assumed to be formatted with the ID separated from any additional fields by a space:

    >id1 ...optional fields...

GFF annotation files must have corresponding IDs specified in mRNA feature ID attributes.
CDS coordinates are derived from CDS features:

    Chr1	.	mRNA	1111	6666	.	.	.	ID=id1
    Chr1	.	CDS	2222	3333	.	.	.	Parent=id1
    Chr1	.	CDS	4444	5555	.	.	.	Parent=id1

GFF files are not required to be sorted.

BED files are assumed to contain a single feature for each primary coding sequence specified, where the coordinates represent the minimum start (0-based) and maximum end positions of the the primary coding sequence in the reference:

	Chr1	2221	5555	id1

Optionally, a file specified in the expected_chr_matches variable can be specified in pandagma.conf,
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

Subommands (in order they are usually run):
            version - Get installed package version
         run ingest - Prepare the assembly and annotation files for analysis
         run mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies
         run filter - Filter the synteny results for chromosome pairings, returning gene pairs.
     run dagchainer - Run DAGchainer to filter for syntenic blocks
            run mcl - Derive clusters, with Markov clustering
       run consense - calculate a consensus sequences from each pan-gene set,
                       If possible add sequences missed in the first clustering round.
      run add_extra - Add other gene model sets to the primary clusters. Useful for adding
                       annotation sets that may be of lower or uncertain quality.
      run summarize - Move results into output directory, and report summary statistics.

Variables in pandagma config file:
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.98]
          clust_cov - Minimum coverage for mmseqs clustering [0.75]
        consen_iden - Minimum identity threshold for vsearch consensus generation [0.80]
         pan_prefix - Prefix to use as a prefix for pangene clusters [default: pan]
       out_dir_base - base name for the output directory [default: './${pkg}_out']
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]
        mcl_threads - Threads to use in Markov clustering [default: processors/5]
            version - version of this script at config time

Environment variables:
         PANDAGMA_CONF - Path of the pandagma config file (default ${PWD})
     PANDAGMA_WORK_DIR - Location of working files (default ${PWD}/work)
                 NPROC - Number of processors to use (default 1)
~~~

## Detailed instructions <a name="details"></a>

1. Clone the program and associated files from github:
      clone https://github.com/legumeinfo/pandagma
      I typically rename the downloaded repository to the genus name that I am working on, e.g.
    
            clone https://github.com/legumeinfo/pandagma.git
            mv pandagma zea
            cd zea

2. Make a work directory. I typically make this in the same area as the pandagma directory, e.g.

            mkdir ../work_zea
            ls ..
              work_zea  zea

3. Get into a suitable work environment (computation node), and load dependencies.
    The script has these third-party dependencies: **bioperl, mcl, mmseqs2**

    These can be loaded using a module-loading system, or with a package manager such as conda, or
    via a Singularity image. 
      
          pandagma_sing_img=$YOURPATH/pandagma-v2021-11-16.sif
    
    OR:
    
          conda activate pandagma    
             # pandagma is the conda environment where I've installed bioperl, mcl, mmseqs2
    
4. Download the annotation data from a remote source (CDS or peptide, and GFF or BED), and transform if needed.
    I do this with a simple shell script that executes curl commands and then applies some transformations.
    See the files in get_samples/ for examples. There are scripts for Glycine, Medicago, Vigna, and Zea.
    
          ./get_samples/get_Zea.sh
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

6. Set environment variables.
    The script learns about the work directory, the number of processors, the config file, and the helper scripts,
    via environment variables.

    If dependencies are manually installed and the script is called directly:

          export NPROC=10
          export PANDAGMA_WORK_DIR=$PWD/../work_zea
          export PANDAGMA_CONF=config/pandagma_Zea_7_2_nuc.conf
          export PATH=$PWD/scripts:$PATH
        
    If the dependencies and script are called via Singularity image: set the environment variables as 
    part of the singularity invocation (see next step)
    
7. Start the run. The examples below assume a run using pandagma_Zea_7_2_nuc.conf.
       
    Calling the program directly:

         nohup pandagma.sh > nohup_7_2.out &

    Using a Singularity image:

         pandagma_sing_img=$YOURPATH/pandagma-v2021-11-16.sif
         nohup singularity exec --env NPROC=$(( ${SLURM_JOB_CPUS_PER_NODE}/5 ))  \
                      --env PANDAGMA_WORK_DIR=$PWD/../work_zea \
                      --env PANDAGMA_CONF=config/pandagma_Zea_7_2_nuc.conf \
                      --cleanenv $pandagma_sing_img pandagma.sh run > nohup_7_2.out &

8. Examine the output, and adjust parameters and possibly the initial chromosome correspondences.
    Output will go into a directory composed from a provided prefix name (default "out") and
    information about key parameter values, e.g.

          out_7_2.NUC.id95.cov60.cns80.ext80.I2/

    The summary of the run is given in the file stats.[parameters].txt .
    Look at the modal values in the histograms, the report of proportion of each assembly with matches, etc.
    One of the output files that may be of use in a subsequent run is observed_chr_pairs.tsv .
    This can indicate possible translocations among genomes in the input data. For example, in Zea, the values
    for one run are:

           1   1 39596   <== Main chromosome correspondences
           5   5 30158
           2   2 30038
           3   3 27658
           4   4 25902
           8   8 23759
           6   6 20969
           7   7 19769
           9   9 19112
          10  10 16833
           9  10   740  <== corresponding with a known translocation in Oh7B
           1   5   544
           1   3   522

    These values can be used to constrain matches in a subsequent run, via the file data/expected_chr_matches.tsv
    (For the Zea run, these values were already provided, generated by get_samples/get_Zea.sh 



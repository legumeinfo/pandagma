# pandagma
Pandagma is a workflow that derives pan-gene sets, given input of gene sequences (CDS and peptide) and
annotations (BED or GFF) from several annotations. It is designed to be suitable for constructing
pan-genes for annotations across several different species within a genus, or annotations within
one species.

Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2023.

The workflow is essentially as follows:
* Add positional information to the gene IDs
* Find gene homologies among pairs of annotation sets, using mmseqs2
* Filter by synteny (DAGChainer) -- and, optionally, by a list of allowed chromosome pairings.
* Cluster (mcl)
* Add back the genes that were excluded to this point 
* Add "extra" annotation sets by homology 
* Identify a representative sequence for each pan-gene
* Calculate consensus order of the pan-genes
* Calculate and report statistics

## Installation methods <a name="installation"></a>

### Installation method 1 (recommended): by creating a Singularity container image

To build a [Singularity](https://singularity.hpcng.org/) container image from the provided [Singularity definition file](https://singularity.hpcng.org/user-docs/master/definition_files.html) (`singularity.def`):

    singularity build [--remote] pandagma.sif singularity.def

Where the `--remote` option is used if building the image as an unprivileged user using (by default) the [Sylabs Cloud Remote Builder](https://cloud.sylabs.io/builder).

To run pandamga using singularity, use `singularity run pandagma.sif [run command]`, e.g.:

  NOTE: THE FOLLOWING STEPS NEED TO BE TESTED AND MODIFIED IN 2023
    singularity pandagma.sif init
    ... modify config/Example.conf ...
    mkdir /path/to/pandagma/work_dir # create work dir on fast, node-local storage
    singularity pandagma.sif -c CONFIG_FILE &

### Installation method 2: manual installation of bin and dependencies

The bin should be added to your path - either by copying them to a directory that is in your path 
(e.g. `cp bin/* ~/bin/`), or by adding the script directory to your path 
(e.g. `export PATH=$PATH:~/src/pandagma/bin` ... if that is where you have them).

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
       nohup ./pandagma.sh -c CONFIG_FILE [options] &
   or
       nohup ./pandagma.sh -c CONFIG_FILE -s SUBCOMMAND [options] &

Primary coding and protein sequences (both fasta) and annotation (GFF3 or BED) files must be listed in the
config file, in the arrays fasta_files, annotation_files, and protein_files. See example files.

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

Subommands (in order they are usually run):
            version - Get installed package version
             ingest - Prepare the assembly and annotation files for analysis
             mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies
             filter - Filter the synteny results for chromosome pairings, returning gene pairs.
         dagchainer - Run DAGchainer to filter for syntenic blocks
                mcl - Derive clusters, with Markov clustering
           consense - Calculate a consensus sequences from each pan-gene set,
                       If possible add sequences missed in the first clustering round.
          add_extra - Add other gene model sets to the primary clusters. Useful for adding
                       annotation sets that may be of lower or uncertain quality.
     filter_to_core - Calculate orthogroup composition and filter fasta files to core orthogroups.
      name_pangenes - Assign pan-gene names with consensus chromosomes and ordinal positions.
     calc_chr_pairs - Report observed chromosome pairs; useful for preparing expected_chr_matches.tsv
          summarize - Move results into output directory, and report summary statistics.

Variables in pandagma config file:
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.60]
        consen_iden - Minimum identity threshold for consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of "extra" annotations [90]
      min_core_prop - Minimum fraction of annotation sets for an orthogroup to be "core" [0.333333333333333]
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]
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

## Detailed instructions <a name="details"></a>

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
    The script has these third-party dependencies: **bioperl, mcl, mmseqs2, DAGchainer**

    These can be loaded using a module-loading system, or with a package manager such as conda, or
    via a Singularity image. 
      
          pandagma_sing_img=$YOURPATH/pandagma-v2021-11-16.sif
    
    OR:
    
          conda activate pandagma    
            # pandagma is the conda environment where I've installed the dependencies above.
    
4. Download the annotation data from a remote source (CDS or peptide, and GFF or BED), and transform if needed.
    I do this with a simple shell script that executes curl commands and then applies some transformations.
    See the files in get_data/ for examples. There are bin for Glycine, Medicago, Vigna, and Zea.
    
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

      nohup ./pandagma.sh -c config/Zea_7_2.conf &


    Using a Singularity image:

  NOTE: THE FOLLOWING STEPS NEED TO BE TESTED AND MODIFIED IN 2023
         pandagma_sing_img=$YOURPATH/pandagma-v_YOUR_VERSION
         nohup singularity exec --cleanenv $pandagma_sing_img ./pandagma.sh -c config/Zea_7_2.conf &

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



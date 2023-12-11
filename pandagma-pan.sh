#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates pan-gene 
# clusters using the programs mmseqs, dagchainer, and mcl. 
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2023
#
scriptname=`basename "$0"`
version="2023-12-10"
set -o errexit -o errtrace -o nounset -o pipefail

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

HELP_DOC="""Compute pan-gene clusters using a combination of synteny and homology,
using the programs mmseqs, dagchainer, and mcl, and additional pre- and post-refinement steps.

Usage:
  ./$scriptname -c CONFIG_FILE [options]

  Required:
           -c (path to the config file)

  Options: -s (subcommand to run. If \"all\" or omitted, all steps will be run; otherwise, run specified step)
           -w (working directory, for temporary and intermediate files. 
                Must be specified in config file if not specified here.)
           -n (number of processors to use. Defaults to number of processors available to the parent process)
           -o (ordering method, for placing pan-genes. Options: 
                \"reference\" (default; uses preferred_annot to order, then gap-filling for missing panIDs.)
                \"alignment\" (uses whole-chromosome alignment of ordered panIDs from all annotations)
           -r (retain. Don't do subcommand \"clean\" after running \"all\".)
           -v (version)
           -h (help)
           -m (more information)

Environment requirements: The following packages need to be available in your PATH:
    mmseqs dagchainer mcl cons famsa hmmalign hmmbuild fasttree

Also, please add the pandagma utility programs in the bin directory adjacent to pandagma-pan.sh, e.g.
    PATH=$PWD/bin:\$PATH

Primary coding and protein sequences (both fasta) and annotation (GFF3 or BED) files must be listed in the
config file, in the arrays cds_files, annotation_files, and protein_files. See example files.

Note that the annotation and CDS files need to be listed in CORRESPONDING ORDER in the config.

Subcommands (in order they are usually run):
                all - All of the steps below, except for clean and ReallyClean
                        (Or equivalently: omit the -s flag; \"all\" is default)
             ingest - Prepare the assembly and annotation files for analysis
             mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies
             filter - Filter the synteny results for chromosome pairings, returning gene pairs.
         dagchainer - Run DAGchainer to filter for syntenic blocks
                mcl - Derive clusters, with Markov clustering
           consense - Calculate a consensus sequences from each pan-gene set, 
                      adding sequences missed in the first clustering round.
       cluster_rest - Retrieve unclustered sequences and cluster those that can be.
          add_extra - Add other gene model sets to the primary clusters. Useful for adding
                      annotation sets that may be of lower or uncertain quality.
     pick_exemplars - Pick representative sequence for each pan-gene
   filter_to_pctile - Calculate orthogroup composition and filter fasta files by selected percentiles.
         tabularize - Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv
     order_and_name - Assign pan-gene names with consensus chromosomes and ordinal positions.
     calc_chr_pairs - Report observed chromosome pairs; useful for preparing expected_chr_matches.tsv
          summarize - Copy results into output directory, and report summary statistics.

  Run the following subcommands separately if you wish:
              align - Align families.
     model_and_trim - Build HMMs and trim the alignments, preparatory to calculating trees.
         calc_trees - Calculate gene trees.
   xfr_aligns_trees - Transfer alignments, HMMs, and trees to output directory

              clean - Clean (delete) files in the working directory that are not needed 
                        for later addition of data using add_extra and subsequent run commands.
                        By default, \"clean\" is run as part of \"all\" unless the -r flag is set.
        ReallyClean - Do complete clean-up of files in the working directory.
                        Use this if you want to start over, OR if you are satisified with the results and
                        don't anticipate adding other annotation sets to this pan-gene set.
"""

MORE_INFO="""
Optionally, a file specified in the expected_chr_matches variable can be specified in pandagma.conf,
which provides anticipated chromosome pairings, e.g.
  01 01
  02 02
  ...
  11 13  # allows for translocation between 11 and 13
  12 12
These pairings are used in a regular expression to identify terminal portions of molecule IDs, e.g.
  glyma.Wm82.gnm2.Gm01  glyso.PI483463.gnm1.Gs01
  glyma.Wm82.gnm2.Gm13  glyso.W05.gnm1.Chr11
If an expected_chr_matches file is not provided, then no such filtering will be done.

Variables in pandagma config file (Set the config with the CONF environment variable)
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.60]
        consen_iden - Minimum identity threshold for consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of \"extra\" annotations [90]
      mcl_inflation - Inflation parameter, for Markov clustering [1.6]
        strict_synt - For clustering of the \"main\" annotations, use only syntenic pairs [1]
                        The alternative (0) is to use all homologous pairs that satisfy expected_chr_matches.tsv
      consen_prefix - Prefix to use in names for genomic ordered consensus IDs [Genus.pan1]
       out_dir_base - Base name for the output directory [default: './out']
           pctl_low - Lower percentile bin size cutoff relative to modal count (usu. the count of annots) [25]
           pctl_med - Medium percentile bin size cutoff relative to modal count (usu. the count of annots) [50]
            pctl_hi - High percentile bin size cutoff relative to modal count (usu. the count of annots) [75]
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g. 
                        \"([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+\" 
                          for four dot-separated fields, e.g. vigan.Shumari.gnm1.ann1
                        or \"(\D+\d+\D+)\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    preferred_annot - String to match and select an annotation set, from a gene ID.
                        This is used for picking representative IDs+sequence from an orthogroup, when
                        this annotation is among those with the median length for the orthogroup.
                        Otherwise, one is selected at random from those with median length.
       order_method - Method to determine consensus panID order. reference or alignment [default reference]
           work_dir - Working directory, for temporary and intermediate files. 
"""

########################################
# Helper functions begin here

version() {
  echo $scriptname $version
}

##########
canonicalize_paths() {
  echo "Entering canonicalize_paths. Fasta files: ${cds_files[@]}"

  cds_files=($(realpath --canonicalize-existing "${cds_files[@]}"))
  annotation_files=($(realpath --canonicalize-existing "${annotation_files[@]}"))
  protein_files=($(realpath --canonicalize-existing "${protein_files[@]}"))
  if (( ${#cds_files_extra_constr[@]} > 0 ))
  then
    cds_files_extra_constr=($(realpath --canonicalize-existing "${cds_files_extra_constr[@]}"))
    annotation_files_extra_constr=($(realpath --canonicalize-existing "${annotation_files_extra_constr[@]}"))
  fi
  if (( ${#cds_files_extra_free[@]} > 0 ))
  then
    cds_files_extra_free=($(realpath --canonicalize-existing "${cds_files_extra_free[@]}"))
    annotation_files_extra_free=($(realpath --canonicalize-existing "${annotation_files_extra_free[@]}"))
  fi
  readonly chr_match_list=${expected_chr_matches:+$(realpath "${expected_chr_matches}")}
  readonly submit_dir=${PWD}

  fasta_file=$(basename "${cds_files[0]}" .gz)
  fna="${fasta_file##*.}"

  export ANN_REX=${annot_str_regex}
}

##########
cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}

##########
check_seq_type () {
  someseq=${1}
  export proportion_nuc=$(echo $someseq | fold -w1 | 
                          awk '$1~/[ATCGN]/ {nuc++} $1!~/[ATCGN]/ {not++} END{print nuc/(nuc+not)}')
  perl -le '$PN=$ENV{"proportion_nuc"}; if ($PN>0.9){print 3} else {print 1}'
}

##########
calc_seq_stats () {
  # Given fasta file on STDIN and an environment variable "annot_name", report:
  # seqs  min  max  N50  ave  annotation_set
  awk 'BEGIN { ORS="" } 
       /^>/ && NR==1 { print $0"\n" } 
       /^>/ && NR!=1 { print "\n"$0"\n" } 
       /^[^>]/ { print } 
       END { print "\n" }' |
  awk '/^[^>]/ {print length($1); tot_bp+=length($1) } END { print "bases: " tot_bp "\n" }' \
   | sort -n \
   | awk -v ANN=$annot_name 'BEGIN { min=99999999999999 }
          /bases:/ { N50_ct = $2/2; bases=$2 } 
          /^[0-9]/ { 
             ct++; sum+=$1; if ( sum >= N50_ct && !printed ) { N50=$1; printed = 1 } 
             min = (min < $1) ? min : $1
          } 
          END { 
            ave=sum/ct;
            printf(" %4d  %4d  %4d  %4d  %7.1f  %4s\n", ct, min, $1, N50, ave, ANN);
          }'
}

########################################
# run functions 

##########
run_ingest() {
# Add positional information from GFF3 or 4- or 6-column BED to FASTA IDs
# BED start coordinate converted to 1-based
  cd "${WORK_DIR}"
  echo; echo "Run ingest: from fasta and gff or bed data, create fasta with IDs containing positional info."
  
  mkdir -p 02_fasta_nuc 02_fasta_prot 01_posn_hsh stats

    # Prepare the tmp.gene_count_start to be joined, in run_summarize, with tmp.gene_count_end_pctl??_end.
    # This is captured from the gene IDs using the annot_str_regex set in the config file.
    cat /dev/null > stats/tmp.gene_count_start
    cat /dev/null > stats/tmp.fasta_seqstats
    start_time=`date`
    printf "Run started at: $start_time\n" > stats/tmp.timing

  export ANN_REX=${annot_str_regex}

  echo "  Get position information from the main annotation sets."
  cat /dev/null > 02_all_main_cds.fna # Collect all starting sequences, for later comparisons
  for (( file_num = 0; file_num < ${#cds_files[@]} ; file_num++ )); do
    file_base=$(basename ${cds_files[file_num]%.*})
    cat_or_zcat "${cds_files[file_num]}" >> 02_all_main_cds.fna # Collect original seqs for later comparisons
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | 
      gff_or_bed_to_hash5.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files[file_num]}" \
                          -hash 01_posn_hsh/$file_base.hsh \
                          -out 02_fasta_nuc/$file_base
    # calc basic sequence stats
    annot_name=$(basename 02_fasta_nuc/$file_base | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
    printf "  Main:  " >> stats/tmp.fasta_seqstats
    cat_or_zcat "${cds_files[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
  done

  echo "  Get position information from the extra annotation sets, if any,"
  echo "  for the annotations that will be added chromosome constraints (_constr)"
  if (( ${#cds_files_extra_constr[@]} > 0 ))
  then
    cat /dev/null > 02_all_extra_cds.fna # Collect all starting sequences, for later comparisons
    for (( file_num = 0; file_num < ${#cds_files_extra_constr[@]} ; file_num++ )); do
      file_base=$(basename ${cds_files_extra_constr[file_num]%.*})
      cat_or_zcat "${cds_files_extra_constr[file_num]}" >> 02_all_extra_cds.fna  # Collect original seqs for later comparisons
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra_constr[file_num]}" | 
        gff_or_bed_to_hash5.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files_extra_constr[file_num]}" \
                            -hash 01_posn_hsh/$file_base.hsh \
                            -out 02_fasta_nuc/$file_base
      # calc basic sequence stats
      annot_name=$(basename 02_fasta_nuc/$file_base | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )

      printf "  Extra: " >> stats/tmp.fasta_seqstats
      cat_or_zcat "${cds_files_extra_constr[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
    done
  fi

  echo "  Get position information from the extra annotation sets, if any,"
  echo "  for the annotations that will be added without chromosome constraints (_free)"
  if (( ${#cds_files_extra_free[@]} > 0 ))
  then
    cat /dev/null > 02_all_extra_cds.fna # Collect all starting sequences, for later comparisons
    for (( file_num = 0; file_num < ${#cds_files_extra_free[@]} ; file_num++ )); do
      file_base=$(basename ${cds_files_extra_free[file_num]%.*})
      cat_or_zcat "${cds_files_extra_free[file_num]}" >> 02_all_extra_cds.fna  # Collect original seqs for later comparisons
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra_free[file_num]}" | 
        gff_or_bed_to_hash5.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files_extra_free[file_num]}" \
                            -hash 01_posn_hsh/$file_base.hsh \
                            -out 02_fasta_nuc/$file_base
      # calc basic sequence stats
      annot_name=$(basename 02_fasta_nuc/$file_base | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )

      printf "  Extra: " >> stats/tmp.fasta_seqstats
      cat_or_zcat "${cds_files_extra_free[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
    done
  fi

  cat /dev/null > 01_combined_posn.hsh # Collect combined file of position hashes
    echo "  Make a combined file of position hashes for later use in step add_extra"
      cat 01_posn_hsh/*.hsh >> 01_combined_posn.hsh

  echo "  Count starting sequences, for later comparisons"
  for file in 02_fasta_nuc/*.$fna; do
    awk '$0~/UNDEFINED/ {ct++} 
      END{if (ct>0){print "Warning: " FILENAME " has " ct " genes without position (HASH UNDEFINED)" } }' $file
    cat $file | grep '>' | perl -pe 's/__/\t/g' | cut -f2 | # extracts the GeneName from combined genspchr__GeneName__start__end__orient
      perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex.+/$1/' |
      grep -v UNDEFINED | sort | uniq -c | awk '{print $2 "\t" $1}' >> stats/tmp.gene_count_start
  done

  echo "  Also get protein files"
  if (( ${#protein_files[@]} > 0 ))
  then
    for (( file_num = 0; file_num < ${#protein_files[@]} ; file_num++ )); do
      echo "  Copying protein file ${protein_files[file_num]}"
      cp "${protein_files[file_num]}" 02_fasta_prot/
    done
  fi

  sort -o stats/tmp.gene_count_start stats/tmp.gene_count_start
}

##########
run_mmseqs() {
  # Do mmseqs clustering on all pairings of the main annotation sets (not the extra ones though)
  cd "${WORK_DIR}"
  echo; echo "Run mmseqs -- at ${clust_iden} percent identity and minimum of ${clust_cov}% coverage."
  #

  if [ -d 03_mmseqs ]; then rm -rf 03_mmseqs ; fi
  mkdir -p 03_mmseqs 03_mmseqs_tmp
  for (( file1_num = 0; file1_num < ${#cds_files[@]} ; file1_num++ )); do
    qry_base=$(basename ${cds_files[file1_num]%.*} .$fna)
    for (( file2_num = $( expr $file1_num + 1 ); file2_num < ${#cds_files[@]} ; file2_num++ )); do
      sbj_base=$(basename ${cds_files[file2_num]%.*} .$fna)
      echo "  Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      { cat 02_fasta_nuc/$qry_base.$fna 02_fasta_nuc/$sbj_base.$fna ; } |
        mmseqs easy-cluster stdin 03_mmseqs/${qry_base}.x.${sbj_base} $MMTEMP \
         --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null & # background
        # allow to execute up to $MMSEQSTHREADS in parallel
        if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    echo
  done
  wait # wait for last jobs to finish
}

##########
run_filter() {
  echo; echo "From mmseqs cluster output, split out the following fields: molecule, gene, start, stop."
  cd "${WORK_DIR}"
  mkdir -p 04_dag
  if [[ -f ${chr_match_list} ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns from file ${chr_match_list}"
    for mmseqs_path in 03_mmseqs/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      echo "  $outfilebase"
      cat ${mmseqs_path} | filter_mmseqs_by_chroms.pl -chr_pat ${chr_match_list} |
        awk 'NF==8' |  # matches for genes with coordinates. The case of <8 can happen for seqs with UNDEF position.
        cat > 04_dag/${outfilebase}_matches.tsv &

      # allow to execute up to $NPROC in parallel
      if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
    done
    wait # wait for last jobs to finish
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches.tsv file was provided, so proceeding without chromosome-pair filtering."
    for mmseqs_path in 03_mmseqs/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      cat ${mmseqs_path} | perl -pe 's/__/\t/g; s/\t[\+-]//g' |
        awk 'NF==8' |  # matches for genes with coordinates. The case of <8 can happen for seqs with UNDEF position.
        cat > 04_dag/${outfilebase}_matches.tsv 
    done
  fi
}

##########
run_dagchainer() {
  # Identify syntenic blocks, using DAGchainer
  cd "${WORK_DIR}"
  dagchainer_args='-M 50 -E 1e-5 -A 6 -s'  # -g and -D are calculated from the data
  echo; echo "Run DAGchainer using arguments \"${dagchainer_args}\" (-g and -D are calculated from the data)"
  # Check and preemptively remove malformed \*_matches.file, which can result from an aborted run
  if [ -f 04_dag/\*_matches.tsv ]; then rm 04_dag/\*_matches.tsv; fi
  for match_path in 04_dag/*_matches.tsv; do
    align_file=`basename $match_path _matches.tsv`
    qryfile=$(echo $align_file | perl -pe 's/(\S+)\.x\..+/$1/')
    sbjfile=$(echo $align_file | perl -pe 's/\S+\.x\.(\S+)/$1/')

    echo "Find average distance between genes for the query and subject files: "
    echo "  $qryfile and $sbjfile"
    ave_gene_gap=$(cat 02_fasta_nuc/$qryfile.$fna 02_fasta_nuc/$sbjfile.$fna | 
                     awk '$1~/^>/ {print substr($1,2)}' | perl -pe 's/__/\t/g' | sort -k1,1 -k3n,3n |
                     awk '$1 == prev1 && $3 > prev4 {sum+=$3-prev4; ct++; prev1=$1; prev3=$3; prev4=$4};
                          $1 != prev1 || $3 <= prev4 {prev1=$1; prev3=$3; prev4=$4}; 
                          END{print 100*int(sum/ct/100)}')
    let "max_gene_gap = ave_gene_gap * 20"

    echo "Running DAGchainer on comparison: $align_file"
    echo "  Calculated DAGchainer parameters: -g (ave_gene_gap): $ave_gene_gap -D (max_gene_gap): $max_gene_gap"; echo
    echo " run_DAG_chainer.pl $dagchainer_args  -g $ave_gene_gap -D $max_gene_gap -i \"${OLDPWD}/${match_path}\""

    # run_DAG_chainer.pl writes temp files to cwd;
    # use per-process temp directory to avoid any data race
    (
      tmpdir=$(mktemp -d)
      cd "${tmpdir}"
      run_DAG_chainer.pl $dagchainer_args  -g $ave_gene_gap -D $max_gene_gap -i "${OLDPWD}/${match_path}" 1>/dev/null
      rmdir ${tmpdir}
    ) &
    # allow to execute up to $NPROC in parallel
    [ $(jobs -r -p | wc -l) -ge ${NPROC} ] && wait -n
  done
  wait # wait for last jobs to finish

  if [ "$strict_synt" -eq 1 ]; then
    # Combine the synteny pairs
    cat 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' | 
      awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
  else 
    # Combine the homology pairs (filtered by chromosome matches if provided) and the synteny pairs
    cat 04_dag/*_matches.tsv 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' | 
      awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
  fi
}

##########
run_mcl() {
  # Calculate clusters using Markov clustering
  cd "${WORK_DIR}"
  printf "\nDo Markov clustering with inflation parameter $mcl_inflation and ${NPROC} threads\n"
  echo "MCL COMMAND: mcl 05_filtered_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv"
  mcl 05_filtered_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv \
    1>/dev/null
 
  echo "  Add cluster IDs"
  awk '{padnum=sprintf("%05d", NR); print "pan" padnum "\t" $0}' tmp.syn_pan.clust.tsv > 06_syn_pan.clust.tsv
  rm tmp.syn_pan.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 06_syn_pan.clust.tsv > 06_syn_pan.hsh.tsv
}

##########
run_consense() {
  echo; 
  printf "\nStep \"consense\" will add previously unclustered"
  printf "\nsequences into an \"augmented\" pan-gene set, by homology.\n"
  cd "${WORK_DIR}"
  if [ -d 07_pan_fasta ]; then rm -rf 07_pan_fasta; fi
  mkdir -p 07_pan_fasta lists

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  echo "    Fasta file:" "${cds_files[@]}"
  get_fasta_from_family_file.pl "${cds_files[@]}" -fam 06_syn_pan.clust.tsv -out 07_pan_fasta

  echo "  Merge fasta files in 07_pan_fasta, prefixing IDs with panID"
  merge_files_to_pan_fasta.awk 07_pan_fasta/* > 07_pan_fasta_cds.fna

  echo "  Pick a representative seq. for each orthogroup - as a sequence with the median length for that OG."
  echo "  cat 07_pan_fasta_cds.fna | 
    pick_family_rep.pl -prefer $preferred_annot -out 08_pan_fasta_clust_rep_cds.fna"
  cat 07_pan_fasta_cds.fna | pick_family_rep.pl -prefer $preferred_annot -out 08_pan_fasta_clust_rep_cds.fna

  echo "  Get sorted list of all genes, from the original fasta files"
  cat_or_zcat "${cds_files[@]}" | awk '/^>/ {print substr($1,2)}' | sort > lists/09_all_genes

  echo "  Get sorted list of all clustered genes"
  awk '$1~/^>/ {print $1}' 07_pan_fasta/* | sed 's/>//' | sort > lists/09_all_clustered_genes

  echo "  Get list of genes not in clusters"
  comm -13 lists/09_all_clustered_genes lists/09_all_genes > lists/09_genes_not_in_clusters

  echo "  Retrieve the non-clustered genes"
  cat_or_zcat "${cds_files[@]}" |
    get_fasta_subset.pl -in /dev/stdin -clobber -lis lists/09_genes_not_in_clusters \
                        -out 09_genes_not_in_clusters.fna 

  echo "  Search non-clustered genes against genes already clustered."

  # Check sequence type (in case this run function is called separately from the usually-prior ones)
  someseq=$(head 07_pan_fasta_cds.fna | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
  SEQTYPE=$(check_seq_type "${someseq}") # 3=nuc; 1=pep
  echo "SEQTYPE is: $SEQTYPE"

  if [ -d 10_place_leftovers ]; then rm -rf 10_place_leftovers; fi
  mkdir -p 10_place_leftovers
  mmseqs easy-search 09_genes_not_in_clusters.fna \
                     07_pan_fasta_cds.fna \
                     10_place_leftovers/unclust.x.07_pan_fasta.m8 \
                     03_mmseqs_tmp \
                     --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null 

  echo "  Filter unclust.x.07_pan_fasta.m8 by clust_iden and report: panID, qry_gene, sbj_gene"
  cat 10_place_leftovers/unclust.x.07_pan_fasta.m8 | top_line.awk |
    awk -v IDEN=${clust_iden} '$3>=IDEN {print $2 "\t" $1}' | 
    perl -pe 's/^(pan\d+)__(\S+)\t(\S)/$1\t$2\t$3/' > 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj.tsv

  echo "  Add gene-pair information"
  hash_into_table_2cols.pl 01_posn_hsh/*hsh -table 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj.tsv |
    cat > 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_posn.tsv

  echo; echo "  Filter based on list of expected chromosome pairings if provided"
  if [[ -f ${chr_match_list} ]]; then  
    echo "Filtering on chromosome patterns from file ${chr_match_list}"
    cat 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_posn.tsv | filter_mmseqs_by_chroms.pl -chr_pat ${chr_match_list} |
      awk 'NF==9' |  # matches for genes with coordinates. The case of <9 can happen for seqs with UNDEF position.
      cat > 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_matches.tsv 
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches.tsv file was provided, so proceeding without chromosome-pair filtering."
    cat 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_posn.tsv | perl -pe 's/__/\t/g; s/\t[\+-]//g' |
      awk 'NF==9' |  # matches for genes with coordinates. The case of <9 can happen for seqs with UNDEF position.
      cat > 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_matches.tsv
  fi

  echo "  Place unclustered genes into their respective pan-gene sets"
  cat 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_matches.tsv | cut -f1,3 | 
  sort -u -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > 11_syn_pan_leftovers.clust.tsv

  echo "  Retrieve sequences for the leftover genes"
  mkdir -p 11_pan_leftovers
  get_fasta_from_family_file.pl "${cds_files[@]}" \
    -fam 11_syn_pan_leftovers.clust.tsv -out 11_pan_leftovers/

  echo "  Make augmented cluster sets. Uniqify on the back end, converting from mcl clust format to hsh (clust_ID gene)."
  cat /dev/null > 12_syn_pan_aug_pre.hsh.tsv
  augment_cluster_sets.awk leftovers_dir=11_pan_leftovers 07_pan_fasta/* |
    perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' |
    sort -k1,1 -k2,2 | uniq > 12_syn_pan_aug_pre.hsh.tsv 

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  cat 12_syn_pan_aug_pre.hsh.tsv | hash_to_rows_by_1st_col.awk > 12_syn_pan_aug_pre.clust.tsv
}

##########
run_cluster_rest() {
  echo
  echo; echo "== Retrieve unclustered sequences and cluster those that can be =="
  cd "${WORK_DIR}"

  echo "  Retrieve genes present in the original CDS files but absent from 12_syn_pan_aug.hsh"
  cut -f2 12_syn_pan_aug_pre.hsh.tsv | LC_ALL=C sort > lists/lis.12_syn_pan_aug_complement
  get_fasta_subset.pl -in 02_all_main_cds.fna -out 12_syn_pan_aug_complement.fna \
    -lis lists/lis.12_syn_pan_aug_complement -xclude -clobber

  MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
  complmt_self_compare="12_syn_pan_aug_complement.x.12_syn_pan_aug_complement"
  cat 12_syn_pan_aug_complement.fna |
    mmseqs easy-cluster stdin 03_mmseqs/$complmt_self_compare $MMTEMP \
    --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null

  echo "  Add gene-pair information"
  hash_into_table_2cols.pl 01_posn_hsh/*hsh -idx1 0 -idx2 1 -table 03_mmseqs/${complmt_self_compare}_cluster.tsv |
    cat > 03_mmseqs/${complmt_self_compare}_cluster_posn.tsv

  if [[ -f ${chr_match_list} ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns from file ${chr_match_list}"
    cat 03_mmseqs/${complmt_self_compare}_cluster_posn.tsv |
      filter_mmseqs_by_chroms.pl -chr_pat ${chr_match_list} > 04_dag/${complmt_self_compare}_matches.tsv
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches.tsv file was provided, so proceeding without chromosome-pair filtering."
    cat 03_mmseqs/${complmt_self_compare}_cluster.tsv |
      perl -pe 's/__/\t/g; s/\t[\+-]$//' > 04_dag/${complmt_self_compare}_matches.tsv 
  fi

  echo "  Extract query and subject matches to be clustered"
  cat 04_dag/${complmt_self_compare}_matches.tsv | cut -f2,6 > 04_dag/${complmt_self_compare}_pairs.tsv

  echo "  Cluster the remaining sequences that have matches"
  mcl  04_dag/${complmt_self_compare}_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan_aug_complement.clust.tsv

  echo "  Find number of clusters in initial (06) results"
  clust_count_06=`wc -l 06_syn_pan.clust.tsv | awk '{print $1}'`
 
  echo "  Add cluster IDs"
  cat tmp.syn_pan_aug_complement.clust.tsv | 
    awk -v START=$clust_count_06 'NF>1 {padnum=sprintf("%05d", NR+START); print "pan" padnum "\t" $0}' |
    cat > 12_syn_pan_aug_complement.clust.tsv
  rm tmp.syn_pan_aug_complement.clust.tsv

  echo "  Combine clusters derived from synteny or restricted homology (12_syn_pan_aug.clust.tsv)"
  echo "  and clusters from the complement of that set. Skip singletons (\"clusters\" of size 1)."
  cat 12_syn_pan_aug_pre.clust.tsv 12_syn_pan_aug_complement.clust.tsv | awk 'NF>1' > 12_syn_pan_aug.clust.tsv

  echo "  Reshape from hash to mcl output format (clustered IDs on one line)"
  cat 12_syn_pan_aug.clust.tsv | hash_to_rows_by_1st_col.awk > 12_syn_pan_aug.hsh.tsv
}

##########
run_add_extra() {
  echo; echo "== Add extra annotation sets (if provided) to the augmented clusters, by homology =="
  cd "${WORK_DIR}"

  if [ -d 13_pan_aug_fasta ]; then rm -rf 13_pan_aug_fasta; fi
  if [ -d 13_extra_out_constr_dir ]; then rm -rf 13_extra_out_constr_dir; fi
  if [ -d 13_extra_out_free_dir ]; then rm -rf 13_extra_out_free_dir; fi
  mkdir -p 13_extra_out_constr_dir 13_extra_out_free_dir 13_pan_aug_fasta

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  get_fasta_from_family_file.pl "${cds_files[@]}" -fam 12_syn_pan_aug.clust.tsv -out 13_pan_aug_fasta
  
  echo "  Merge fasta files in 13_pan_aug_fasta, prefixing IDs with panID__"
  merge_files_to_pan_fasta.awk 13_pan_aug_fasta/* > 13_pan_aug_fasta.fna

  echo "  Store the panID - geneID in a hash to be retrieved later, after filtering by chromosome"
  awk '$1~/^>/ {print substr($1,2)}' 13_pan_aug_fasta.fna | 
    perl -pe 's/(pan\d+)__(.+)/$2\t$1/' | sort -k1,1 -k2,2 > 13_pan_gene.hsh
  perl -pi -e 's/>pan\d+__/>/' 13_pan_aug_fasta.fna

  echo "  Add position information to 13_pan_aug_fasta.fna"
  hash_into_fasta_id.pl -nodef -fasta 13_pan_aug_fasta.fna \
                        -hash 01_combined_posn.hsh \
                        -out 13_pan_aug_fasta_posn.fna

  if (( ${#cds_files_extra_constr[@]} > 0 )) || (( ${#cds_files_extra_free[@]} > 0 ))
  then # handle the "extra" annotation files
    echo "  Search non-clustered genes against pan-gene consensus sequences"
    # Check sequence type (in case this run function is called separately from the usually-prior ones)
    someseq=$(head 07_pan_fasta_cds.fna | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
    SEQTYPE=$(check_seq_type "${someseq}") # 3=nuc; 1=pep
    echo "SEQTYPE is: $SEQTYPE"

    for file in "${cds_files_extra_constr[@]}"; do
      file_base=$(basename $file .gz)
      echo "Extra: 02_fasta_nuc/$file_base"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      mmseqs easy-search "02_fasta_nuc/$file_base" 13_pan_aug_fasta_posn.fna 13_extra_out_constr_dir/${file_base}.x.all_cons.m8 \
                   $MMTEMP --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null & # background

      if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    wait # wait for jobs to finish

    for file in "${cds_files_extra_free[@]}"; do
      file_base=$(basename $file .gz)
      echo "Extra: 02_fasta_nuc/$file_base"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      mmseqs easy-search "02_fasta_nuc/$file_base" 13_pan_aug_fasta_posn.fna 13_extra_out_free_dir/${file_base}.x.all_cons.m8 \
                   $MMTEMP --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null & # background

      if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    wait # wait for jobs to finish

    if [[ -f ${chr_match_list} ]]; then  # filter based on list of expected chromosome pairings if provided
      echo "Filtering on chromosome patterns from file ${chr_match_list} and by identity and top match."
      echo "First find top match without respect to chromosome match, then with chromosome matches."
      echo "A top hit with chromosome match trumps a top hit on the wrong chromosome."
      for m8file in 13_extra_out_constr_dir/*.m8; do
        base=`basename $m8file .m8`
        echo "  Processing file $m8file"
        cat $m8file | 
            top_line.awk | awk -v IDEN=${extra_iden} '$3>=IDEN {print $1 "\t" $2}' |
            perl -pe 's/^\S+__(\S+)__\d+__\d+\t\S+__(\S+)__\d+__\d+$/$1\t$2/' > 13_extra_out_constr_dir/$base.top
        cat $m8file | filter_mmseqs_by_chroms.pl -chr_pat ${chr_match_list} |
            top_line.awk | awk -v IDEN=${extra_iden} '$3>=IDEN {print $1 "\t" $2}' | 
            perl -pe 's/^\S+__(\S+)__\d+__\d+\t\S+__(\S+)__\d+__\d+$/$1\t$2/' > 13_extra_out_constr_dir/$base.top
      done
    else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
      echo "WARNING: No expected_chr_matches.tsv file was provided, but annotations were indicated in the "
      echo "file sets annotation_files_extra_constr and cds_files_extra_constr. Please check that "
      echo "expected_chr_matches.tsv is in the data directory OR that extra annotations are in the "
      echo "annotation_files_extra_free and cds_files_extra_free file collections."
      exit 1;
    fi

    for m8file in 13_extra_out_free_dir/*.m8; do
      base=`basename $m8file .m8`
      echo "  Processing file $m8file"
      cat $m8file | 
          top_line.awk | awk -v IDEN=${extra_iden} '$3>=IDEN {print $1 "\t" $2}' |
          perl -pe 's/^\S+__(\S+)__\d+__\d+__[+-]\t\S+__(\S+)__\d+__\d+__[+-]$/$1\t$2/' > 13_extra_out_free_dir/$base.top
    done

    echo "Join panIDs to unclustered genes"
    cat 13_extra_out_*_dir/*top | sort -k2,2 -k1,1 | join -1 2 -2 1 - 13_pan_gene.hsh | 
      awk 'NF==3 {print $3 "\t" $2}' | sort -k1,1 -k2,2 > 14_syn_pan_extra.hsh.tsv

    echo "Derive a cluster-format file from the hash of panIDs and extra genes"
    cat 14_syn_pan_extra.hsh.tsv | hash_to_rows_by_1st_col.awk > 14_syn_pan_extra.clust.tsv

    echo "  Retrieve sequences for the extra genes"
    if [ -d 16_pan_leftovers_extra ]; then rm -rf 16_pan_leftovers_extra; fi
    mkdir -p 16_pan_leftovers_extra
    get_fasta_from_family_file.pl "${cds_files_extra_constr[@]}" "${cds_files_extra_free[@]}" \
       -fam 14_syn_pan_extra.clust.tsv -out 16_pan_leftovers_extra/
  
    echo "  Make augmented cluster sets. Uniqify on the back end, converting from mcl clust format to hsh (clust_ID gene)."
    augment_cluster_sets.awk leftovers_dir=16_pan_leftovers_extra 13_pan_aug_fasta/* |
      perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' |
      sort -k1,1 -k2,2 | uniq > 18_syn_pan_aug_extra.hsh.tsv 

    echo "  Reshape from hash to mcl output format (clustered IDs on one line)"
    cat 18_syn_pan_aug_extra.hsh.tsv | hash_to_rows_by_1st_col.awk > 18_syn_pan_aug_extra.clust.tsv

    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    echo "    Fasta file:" "${protein_files[@]}"
    if [ -d 19_pan_aug_leftover_merged_cds ]; then rm -rf 19_pan_aug_leftover_merged_cds; fi
    mkdir -p 19_pan_aug_leftover_merged_cds
    get_fasta_from_family_file.pl "${cds_files[@]}" "${cds_files_extra_constr[@]}" "${cds_files_extra_free[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_cds
  
    echo "  Merge files in 19_pan_aug_leftover_merged_cds, prefixing IDs with panID__"
    merge_files_to_pan_fasta.awk 19_pan_aug_leftover_merged_cds/* > 19_pan_aug_leftover_merged_cds.fna

  else  
    echo "== No annotations were designated as \"extra\", so just promote the syn_pan_aug files as syn_pan_aug_extra. ==" 
    cp 07_pan_fasta_cds.fna 19_pan_aug_leftover_merged_cds.fna
    cp 12_syn_pan_aug.clust.tsv 18_syn_pan_aug_extra.clust.tsv
    cp 12_syn_pan_aug.hsh.tsv 18_syn_pan_aug_extra.hsh.tsv
  fi
}

##########
run_pick_exemplars() {
  echo; echo "== Pick representative (exemplar) sequence for each pan-gene set (protein and CDS) =="
  cd "${WORK_DIR}"

  echo "  Get all protein sequences corresponding with 18_syn_pan_aug_extra.clust.tsv"
  cat /dev/null > 20_pan_fasta_prot.faa
  for filepath in 02_fasta_prot/*.gz; do 
    zcat $filepath >> 20_pan_fasta_prot.faa
  done

  echo "  Get protein sequences into pan-gene sets, corresponding with 18_syn_pan_aug_extra.clust.tsv"
  if [ -d 19_pan_aug_leftover_merged_prot ]; then rm -rf 19_pan_aug_leftover_merged_prot; fi
  mkdir -p 19_pan_aug_leftover_merged_prot
  get_fasta_from_family_file.pl 20_pan_fasta_prot.faa \
              -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_prot

  echo "  Get all protein sequences from files in 19_pan_aug_leftover_merged_prot"
  merge_files_to_pan_fasta.awk 19_pan_aug_leftover_merged_prot/* > 19_pan_aug_leftover_merged_prot.faa

  echo "  Pick a representative sequence for each pangene set - as a sequence with the median length for that set."
  echo "    == first proteins:"
  cat 19_pan_aug_leftover_merged_prot.faa | pick_family_rep.pl \
    -nostop -prefer $preferred_annot -out 21_pan_fasta_clust_rep_prot.faa

  echo "    == then CDS sequences, corresponding with 21_pan_fasta_clust_rep_prot.faa"
  cat 21_pan_fasta_clust_rep_prot.faa | awk '$1~/^>/ {print substr($1,2)}' > lists/lis.21_pan_fasta_clust_rep
  get_fasta_subset.pl -in 19_pan_aug_leftover_merged_cds.fna -list lists/lis.21_pan_fasta_clust_rep \
                    -clobber -out 21_pan_fasta_clust_rep_cds.fna

  perl -pi -e 's/__/  /' 21_pan_fasta_clust_rep_cds.fna
  perl -pi -e 's/__/  /' 21_pan_fasta_clust_rep_prot.faa

  echo "  Retrieve genes present in the original CDS files but absent from 18_syn_pan_aug_extra"
  cut -f2 18_syn_pan_aug_extra.hsh.tsv | LC_ALL=C sort > lists/lis.18_syn_pan_aug_extra
  cat 02_all_main_cds.fna 02_all_extra_cds.fna > 02_all_cds.fna
  get_fasta_subset.pl -in 02_all_cds.fna -out 18_syn_pan_aug_extra_complement.fna \
    -lis lists/lis.18_syn_pan_aug_extra -xclude -clobber
}

##########
run_tabularize() {
  echo; echo "== Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv"
  cd "${WORK_DIR}"

  # Get table header 
  pangene_tabularize.pl -pan 18_syn_pan_aug_extra.clust.tsv -annot_str_regex $ANN_REX > tmp.18_syn_pan_aug_extra.clust.tsv
  head -1 tmp.18_syn_pan_aug_extra.clust.tsv > tmp.table_header

  # Find column on which to sort first
  export PR_ANN=$preferred_annot
  sort_col=$( cat tmp.table_header | 
     perl -ane 'BEGIN{use List::Util qw(first); $PA=$ENV{"PR_ANN"}}; $idx = first { $F[$_] =~ /$PA/ } 0..$#F; print ($idx+1) ' )

  echo "  Field for sorting 18_syn_pan_aug_extra.table.tsv is sort_col: $sort_col"

  echo "  Sort, putting header row at top, and don't print pangenes that are all NONE"
    cat tmp.18_syn_pan_aug_extra.clust.tsv |
    sort -k$sort_col,$sort_col -k2,2 | 
    sed '/^$/d; /^#pangene/d' |
    perl -lane '$ct=0; for $gn (@F){if ($gn=~/NONE/){$ct++}}; if ($ct<(scalar(@F)-1)){print $_}' |
    cat tmp.table_header - > 18_syn_pan_aug_extra.table.tsv

    rm tmp.18_syn_pan_aug_extra.clust.tsv
    rm tmp.table_header
}

##########
run_filter_to_pctile() {
  cd "${WORK_DIR}"

  echo "  Calculate matrix of gene counts per orthogroup and annotation set"
  calc_pan_stats.pl -annot_regex $ANN_REX -pan 18_syn_pan_aug_extra.clust.tsv -out 18_syn_pan_aug_extra.counts.tsv
  max_annot_ct=$(cat 18_syn_pan_aug_extra.counts.tsv | 
                       awk '$1!~/^#/ {print $2}' | sort -n | uniq | tail -1)

  echo "  Select orthogroups with genes from selected percentiles of max_annot_ct annotation sets"
  for percentile in $pctl_low $pctl_med $pctl_hi; do
    cat 18_syn_pan_aug_extra.counts.tsv | 
      awk -v PCTL=$percentile -v ANNCT=$max_annot_ct '$2>=ANNCT*(PCTL/100)&& $1!~/^#/ {print $1}' |
        cat > lists/lis.18_syn_pan_aug_extra.pctl${percentile}

    echo "  Get a fasta subset with only genes from at least $(($max_annot_ct*$pctl_low/100)) annotation sets"
    get_fasta_subset.pl -in 21_pan_fasta_clust_rep_cds.fna \
                        -list lists/lis.18_syn_pan_aug_extra.pctl${percentile} \
                        -clobber -out 22_pan_fasta_rep_pctl${percentile}_cds.fna

    echo "  Get a clust.tsv file with orthogroups with at least min_core_prop*max_annot_ct annotation sets"
    join <(LC_ALL=C sort -k1,1 lists/lis.18_syn_pan_aug_extra.pctl${percentile}) \
         <(LC_ALL=C sort -k1,1 18_syn_pan_aug_extra.clust.tsv) |
            cat > 22_syn_pan_aug_extra_pctl${percentile}.clust.tsv
  done
}

##########
run_order_and_name() {
  cd "${WORK_DIR}"

  echo "  Reshape from mcl output format, clustered IDs on one line, to a hash format"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 22_syn_pan_aug_extra_pctl${pctl_low}.clust.tsv |
    cat > 22_syn_pan_aug_extra_pctl${pctl_low}.hsh.tsv

  echo "  Add positional information to the hash output."
  join -a 1 -1 2 -2 1 <(sort -k2,2 22_syn_pan_aug_extra_pctl${pctl_low}.hsh.tsv) <(cat 01_posn_hsh/*hsh | sort -k1,1) | 
    perl -pe 's/__/\t/g; s/ /\t/g' | awk -v OFS="\t" '{print $2, $1, $3, $5, $6, $7}' |
    sort -k1,1 -k2,2 > 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv

  echo "Order method: $order_method"
  if [[ $order_method =~ "align" ]]; then
    echo "==  The next several steps will determine panID ordering using the \"alignment\" ordering option,"
    echo "    based on alignment of panID orders in each annotation, for each chromosome."
    echo "  If code_table/pan_to_peptide.tsv doesn't exist, generate it."

    mkdir -p code_table
    if [ ! -f code_table/pan_to_peptide.tsv ]; then
      echo "   Generating hash file code_table/pan_to_peptide.tsv"
      make_peptide_hash.pl > code_table/pan_to_peptide.tsv
    fi
  
    echo "  Encode pan-genes as unique peptide strings, and extract annotation sets with encoded"
    echo "  IDs ordered along chromosomes, permitting whole-chromosome alignments of gene order."
    echo "  Calling script order_encode.pl on file 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv"
    mkdir -p 23_encoded_chroms
    order_encode.pl -pan_table 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv \
                          -code_table code_table/pan_to_peptide.tsv -utilized code_table/utilized.tsv \
                          -annot_regex "${ANN_REX}" -outdir 23_encoded_chroms 
  
    echo "  Align chromosome sequences with peptide-encoded-gene-orders."
    mkdir -p 23_encoded_chroms_aligned
    do_align.sh 23_encoded_chroms 23_encoded_chroms_aligned $(( $NPROC/2 ))
  
    echo "  Filter chromosome alignments to consensus motifs"
    mkdir -p 23_encoded_chroms_filt1
    for filepath in 23_encoded_chroms_aligned/*; do 
      base=`basename $filepath`
      cons -sequence $filepath -outseq 23_encoded_chroms_filt1/$base -name $base &
      if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    wait
  
    echo "  Decode motifs to recover pangene IDs in order along each chromosome."
    cat /dev/null > consen_gene_order.tsv
    for filepath in 23_encoded_chroms_filt1/*; do
      chr=`basename $filepath`
      order_decode.pl -align $filepath -code code_table/utilized.tsv -verbose |
        top_line.awk >> consen_gene_order.tsv 
    done
  
    echo "  Find panIDs that weren't placed"
    comm -23 <( cut -f1 code_table/utilized.tsv ) <( cut -f1 consen_gene_order.tsv | sort -u ) |
      awk '{print $1}' > consen_pan_unplaced.txt

  elif [[ $order_method =~ "reference" ]]; then
    echo "==  The next several steps will determine panID ordering using the \"reference\" ordering option,"
    echo "    based on the panID order from the supplied \"preferred annotation\", $preferred_annot"
    echo " order_by_reference.pl -pan_table 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv \\ "
    echo "   -annot_regex \"${ANN_REX}\" -pref_annot $preferred_annot \\ "
    echo "   -consen_out consen_gene_order.tsv -unplaced_out consen_pan_unplaced.txt -verbose "
    order_by_reference.pl -pan_table 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv \
                          -annot_regex "${ANN_REX}" -pref_annot $preferred_annot \
                          -consen_out consen_gene_order.tsv -unplaced_out consen_pan_unplaced.txt -verbose 
  else 
    echo "Ordering method (-o $order_method ) is not one of \"alignment\" or \"reference\". Aborting."
    exit 1
  fi # End of $order_method =~ /align/

  gapfill_threads=$(($NPROC/2))
  echo "  Fill gaps in the alignment-based pangene ordering. This step is time-consuming."
  echo "  so is run in parallel, using $gapfill_threads threads."
  echo "   order_gapfill.pl -verbose -annot_regex \"${ANN_REX}\" -consen consen_gene_order.tsv \\";
  echo "     -pan_table 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv \\";
  echo "     -unplaced consen_pan_unplaced.txt -out consen_gene_order_gapfilled.tsv -nproc $gapfill_threads";
  order_gapfill.pl -verbose -annot_regex "${ANN_REX}" -consen consen_gene_order.tsv \
    -pan_table 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv \
    -unplaced consen_pan_unplaced.txt -out consen_gene_order_gapfilled.tsv -nproc $gapfill_threads

  echo "  Reshape defline into a hash, e.g. pan20175 Phaseolus.pan2.chr11_222300_pan20175 222300 223300 -"
  echo "  Note: these \"positions\" and sizes are artificial, representing only inferred relative positions."
  cat consen_gene_order_gapfilled.tsv | awk -v PRE=$consen_prefix '
    {printf("%s\t%s.%s_%06d_%s\t%d\t%d\t%s\n", $1, PRE, $2, $3, $1, $3*100, $3*100+1000, $4)}' > consen_posn.hsh

  echo "  Extract a version of the consen_gene_order file with immediate tandem duplicates collapsed"
  cat consen_posn.hsh | top_line.awk > consen_posn_trim.hsh

  echo "  Hash position information into fasta file (CDS)"
  hash_into_fasta_id.pl -fasta 22_pan_fasta_rep_pctl${pctl_low}_cds.fna -hash consen_posn_trim.hsh |
    grep -v "HASH UNDEFINED" > 22_pan_fasta_rep_pctl${pctl_low}_posn_cdsTMP.fna

  echo "  Reshape defline, and sort by position (CDS)"
  fasta_to_table.awk 22_pan_fasta_rep_pctl${pctl_low}_posn_cdsTMP.fna | sed 's/__/\t/g; s/ /\t/' | 
    sort | awk '{print ">" $1, $2, $3, $4, $5; print $6}' > 23_syn_pan_pctl${pctl_low}_posn_cds.fna

  echo "  Also get corresponding protein sequences for genes in lists/lis.18_syn_pan_aug_extra.pctl${pctl_low}"
  if [ ! -f 20_pan_fasta_prot.faa ]; then
    for filepath in 02_fasta_prot/*.gz; do 
      zcat $filepath >> 20_pan_fasta_prot.faa
    done
  fi

  cat 23_syn_pan_pctl${pctl_low}_posn_cds.fna | awk '$1~/^>/ {print $5}' > lists/lis.23_syn_pan_pctl${pctl_low}_posn
  get_fasta_subset.pl -in 20_pan_fasta_prot.faa -list lists/lis.23_syn_pan_pctl${pctl_low}_posn \
                      -clobber -out 23_syn_pan_pctl${pctl_low}_posn_proteinTMP.faa

  echo "  Get directory of protein multifasta sequences for each pangene, for (separate) protein alignments"

  if [ -d 22_syn_pan_aug_extra_pctl${pctl_low} ]; then rm -rf 22_syn_pan_aug_extra_pctl${pctl_low} ; fi
  mkdir -p 22_syn_pan_aug_extra_pctl${pctl_low}
  get_fasta_from_family_file.pl 20_pan_fasta_prot.faa \
    -family_file 22_syn_pan_aug_extra_pctl${pctl_low}.clust.tsv -out_dir 22_syn_pan_aug_extra_pctl${pctl_low}

  echo "  Hash pan-ID into protein fasta file"
  cat 23_syn_pan_pctl${pctl_low}_posn_cds.fna |
    awk '$1~/^>/ { print $5 "\t" substr($1,2) "__" $2 "__" $3 "__" $4 }' > lists/23_syn_pan_pctl${pctl_low}_posn.hsh
  hash_into_fasta_id.pl -hash lists/23_syn_pan_pctl${pctl_low}_posn.hsh -swap_IDs -nodef \
    -fasta 23_syn_pan_pctl${pctl_low}_posn_proteinTMP.faa |
    perl -pe 's/__/ /g' > 23_syn_pan_pctl${pctl_low}_posn_prot.faa
}

##########
run_align_cds() {
  cd "${WORK_DIR}"

  echo "== Align CDS sequences =="

  # If not already present, retrieve sequences for each family, preparatory to aligning them
  if [[ -d 19_pan_aug_leftover_merged_cds ]]; then
    : # do nothing; the directory and file(s) exist
  else
    mkdir -p 19_pan_aug_leftover_merged_cds
    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    get_fasta_from_family_file.pl "${cds_files[@]}" "${cds_files_extra_constr[@]}" "${cds_files_extra_free[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_cds
  fi

  echo; echo "== Align nucleotide sequence for each gene family =="
  mkdir -p 20_aligns_cds
  for filepath in 19_pan_aug_leftover_merged_cds/*; do
    file=`basename $filepath`;
    echo "  Computing alignment, using program famsa, for file $file"
    famsa -t 2 19_pan_aug_leftover_merged_cds/$file 20_aligns_cds/$file 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_align_protein() {
  cd "${WORK_DIR}"

  echo "== Align protein sequences =="

  # If not already present, retrieve sequences for each family, preparatory to aligning them
  if [[ -d 19_pan_aug_leftover_merged_prot ]]; then
    : # do nothing; the directory and file(s) exist
  else
    mkdir -p 19_pan_aug_leftover_merged_prot
    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    get_fasta_from_family_file.pl "${protein_files[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_prot
  fi

  echo; echo "== Align protein sequence for the each gene family =="
  mkdir -p 20_aligns_prot
  for filepath in 19_pan_aug_leftover_merged_prot/*; do
    file=`basename $filepath`;
    echo "  Computing alignment, using program famsa, for file $file"
    famsa -t 2 19_pan_aug_leftover_merged_prot/$file 20_aligns_prot/$file 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_model_and_trim() {
  echo; echo "== Build HMMs =="
  cd "${WORK_DIR}"
  mkdir -p 21_hmm
  for filepath in 20_aligns_cds/*; do
    file=`basename $filepath`;
    hmmbuild -n $file 21_hmm/$file $filepath 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Realign to HMMs =="
  mkdir -p 22_hmmalign
  for filepath in 21_hmm/*; do
    file=`basename $filepath`;
    printf "$file "
    hmmalign --trim --outformat A2M -o 22_hmmalign/$file 21_hmm/$file 19_pan_aug_leftover_merged_cds/$file &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo

  echo; echo "== Trim HMM alignments to match-states =="
  mkdir -p 23_hmmalign_trim1
  for filepath in 22_hmmalign/*; do
    file=`basename $filepath`;
    printf "$file "
    cat $filepath |
      perl -ne 'if ($_ =~ />/) {print $_} else {$line = $_; $line =~ s/[a-z]//g; print $line}' |
      sed '/^$/d' > 23_hmmalign_trim1/$file &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo

  echo; echo "== Filter alignments prior to tree calculation =="
  mkdir -p 23_hmmalign_trim2 23_hmmalign_trim2_log
  min_depth=3
  min_pct_depth=20
  min_pct_aligned=20
  for filepath in 23_hmmalign_trim1/*; do
    file=`basename $filepath`
    printf "$file "
    filter_align.pl -in $filepath -out 23_hmmalign_trim2/$file -log 23_hmmalign_trim2_log/$file \
                    -depth $min_depth -pct_depth $min_pct_depth -min_pct_aligned $min_pct_aligned &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo
}

##########
run_calc_trees() {
  echo; echo "== Calculate trees =="
  cd "${WORK_DIR}"

  mkdir -p 24_trees

  echo; echo "== Move small (<4) and very low-entropy families (sequences are all identical) to the side =="
  mkdir -p 23_pan_aug_small_or_identical
  min_seq_count=4
  # Below, "count" is the number of unique sequences in the alignment.
  for filepath in 23_hmmalign_trim2/*; do
    file=`basename $filepath`
    count=$(awk '$1!~/>/ {print FILENAME "\t" $1}' $filepath | sort -u | wc -l);
    if [[ $count -lt $min_seq_count ]]; then
      echo "Set aside small or low-entropy family $file";
      mv $filepath 23_pan_aug_small_or_identical/
    fi;
  done

  # By default, FastTreeMP uses all available threads.
  # It is more efficient to run more jobs on one core each by setting an environment variable.
  OMP_NUM_THREADS=1
  export OMP_NUM_THREADS
  for filepath in 23_hmmalign_trim2/*; do
    file=`basename $filepath`
    echo "  Calculating tree for $file"
    fasttree -nt -quiet $filepath > 24_trees/$file &
  done
  wait
}

##########
run_xfr_aligns_trees() {
  echo; echo "Copy alignment and tree results into output directory"

  conf_base=`basename $CONF .conf`
  full_out_dir="${out_dir_base}_$conf_base"

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  for dir in 20_aligns_cds 20_aligns_prot 21_hmm 22_hmmalign 23_hmmalign_trim2 24_trees; do
    if [ -d ${WORK_DIR}/$dir ]; then
      echo "Copying directory $dir to output directory"
      cp -r ${WORK_DIR}/$dir ${full_out_dir}/
    else
      echo "Warning: couldn't find dir ${WORK_DIR}/$dir; skipping"
    fi
  done
}

##########
run_calc_chr_pairs() {
  cd "${WORK_DIR}"
  echo "Generate a report of observed chromosome pairs"

  echo "  Identify gene pairs, using mmseqs --easy_cluster"
    mmseqs easy-cluster 19_pan_aug_leftover_merged_cds.fna 24_pan_fasta_clust 03_mmseqs_tmp \
    --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null

  echo "   Extract chromosome-chromosome correspondences"
  cut -f2,3 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv | sort -k1,1 > 22_syn_pan_aug_extra_pctl${pctl_low}_posn.gene_chr.hsh

  echo "   From mmseqs cluster table, prune gene pairs, keeping only those in the same pangene cluster"
  cat 24_pan_fasta_clust_cluster.tsv | perl -pe 's/__/\t/g' | awk '$1==$3' |
    perl -pe 's/(pan\d+)\t/$1__/g' > 24_pan_fasta_cluster_pruned.tsv

  echo "   Join chromosome numbers to gene pairs and count chromosome correspondences among the pairs."
  join <(perl -pe 's/pan\d+__//g' 24_pan_fasta_cluster_pruned.tsv | sort -k1,1) \
       22_syn_pan_aug_extra_pctl${pctl_low}_posn.gene_chr.hsh | perl -pe 's/ /\t/g' | sort -k2,2 | 
     join -1 2 -2 1 - 22_syn_pan_aug_extra_pctl${pctl_low}_posn.gene_chr.hsh | 
     awk 'BEGIN{IGNORECASE=1; OFS="\t"} 
          $3!~/cont|scaff|sc|pilon|ctg|contig|tig|mito|mt$|cp$|pt$|chl|unanchor|unkn/ && \
          $4!~/cont|scaff|sc|pilon|ctg|contig|tig|mito|mt$|cp$|pt$|chl|unanchor|unkn/ \
          {print $3, $4}' | perl -pe 's/^\S+\.\D+(\d+)\s+\S+\.\D+(\d+)/$1\t$2/' |
          perl -pe 's/^0(\d+)/$1/; s/\t0(\d+)/\t$1/' |
     awk -v OFS="\t" '$1<=$2 {print $1, $2} $1>$2 {print $2, $1}' |
     sort | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' | sort -k3nr,3nr > observed_chr_pairs.tsv
}

##########
run_summarize() {
  echo; echo "Summarize: Copy results into output directory, and report some summary statistics"
 
  # Determine if the sequence looks like nucleotide or like protein.
  someseq=$(head ${WORK_DIR}/07_pan_fasta_cds.fna | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
  SEQTYPE=$(check_seq_type "${someseq}")
  case "$SEQTYPE" in 
    3) ST="NUC" ;;
    1) ST="PEP" ;;
  esac

  conf_base=`basename $CONF .conf`
  full_out_dir="${out_dir_base}_$conf_base"
  stats_file=${full_out_dir}/stats.$conf_base.txt
  export ANN_REX=${annot_str_regex}

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  for file in 06_syn_pan.clust.tsv 06_syn_pan.hsh.tsv \
              12_syn_pan_aug.clust.tsv 12_syn_pan_aug.hsh.tsv \
              18_syn_pan_aug_extra.clust.tsv  18_syn_pan_aug_extra.hsh.tsv 18_syn_pan_aug_extra.counts.tsv \
              18_syn_pan_aug_extra_complement.fna 18_syn_pan_aug_extra.table.tsv \
              21_pan_fasta_clust_rep_cds.fna 21_pan_fasta_clust_rep_prot.faa \
              22_pan_fasta_rep_pctl${pctl_low}_cds.fna \
              22_pan_fasta_rep_pctl${pctl_med}_cds.fna \
              22_pan_fasta_rep_pctl${pctl_hi}_cds.fna \
              22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv \
              23_syn_pan_pctl${pctl_low}_posn_cds.fna 23_syn_pan_pctl${pctl_low}_posn_prot.faa \
              observed_chr_pairs.tsv ; do
    if [ -f ${WORK_DIR}/$file ]; then
      cp ${WORK_DIR}/$file ${full_out_dir}/
    else 
      echo "Warning: couldn't find file ${WORK_DIR}/$file; skipping"
    fi
  done

  echo "Copy manifest file into the output directory"
  if [ -f "${submit_dir}/manifests/MANIFEST_output_pan.yml" ]; then
    cp "${submit_dir}/manifests/MANIFEST_output_pan.yml" $full_out_dir/
  else 
    echo "Couldn't find file manifests/MANIFEST_output_pan.yml"
  fi

  printf "Run of program $scriptname, version $version\n" > ${stats_file}

  end_time=`date`
  cat ${WORK_DIR}/stats/tmp.timing >> ${stats_file}
  printf "Run ended at:   $end_time\n\n" >> ${stats_file}

  # Report sequence type 
  if [[ "$SEQTYPE" == 3 ]]; then
    echo "Sequence type: nucleotide" >> ${stats_file}
  else
    echo "Sequence type: protein" >> ${stats_file}
  fi

  echo "  Report parameters from config file"
  printf "Parameter  \tvalue\n" >> ${stats_file}
  for key in ${pandagma_conf_params}; do
    printf '%-15s\t%s\n' ${key} "${!key}" >> ${stats_file}
  done

  printf "\nOutput directory for this run:\t${full_out_dir}\n" >> ${stats_file}

  echo "  Report threshold for inclusion in \"pctl${pctl_low}\""
  max_annot_ct=$(cat ${full_out_dir}/18_syn_pan_aug_extra.counts.tsv |
                       awk '$1!~/^#/ {print $2}' | sort -n | uniq | tail -1)
  pctl_low_threshold=$(awk -v MCP=$pctl_low -v MAC=$max_annot_ct 'BEGIN{print MCP*MAC/100}')

  echo "  Report orthogroup composition statistics for the three main cluster-calculation steps"

  printf '\n  %-20s\t%s\n' "Statistic" "value" >> ${stats_file}
  printf "The global mode may be for a smaller OG size. Modes below are greater than the specified core threshold.\n" \

  clustcount=1
  # When the number of main annotation sets is small (<4), the core-threshold-ceiling may be lower than 
  # the largest number of clusters. In that case, set $CTceil=2 to the smallest cluster size ($CTceil=2)
  for clustering in 06_syn_pan 12_syn_pan_aug 18_syn_pan_aug_extra; do
    clustfile=${full_out_dir}/$clustering.clust.tsv
    if [[ -f $clustfile ]]; then
      if [[ $clustcount == 1 ]]; then
        printf "\n== Initial clusters (containing only genes within synteny blocks)\n" >> ${stats_file}
        (( largest=$(awk '{print NF-1}' $clustfile | sort -n | tail -1) ))
        (( mode=$(awk "{print NF-1}" $clustfile |
          sort -n | uniq -c | awk '{print $1 "\t" $2}' | sort -n | tail -1 | awk '{print $2}') ))
      elif [[ $clustcount == 2 ]]; then
        printf "\n== Augmented clusters (unanchored sequences added to the initial clusters)\n" >> ${stats_file}
        (( largest=$(awk '{print NF-1}' $clustfile | sort -n | tail -1) ))
        (( mode=$(awk "{print NF-1}" $clustfile |
          sort -n | uniq -c | awk '{print $1 "\t" $2}' | sort -n | tail -1 | awk '{print $2}') ))
      elif [[ $clustcount == 3 ]]; then
        printf "\n== Augmented-extra clusters (with sequences from extra annotation sets)\n" >> ${stats_file}
        (( largest=$(awk '{print NF-1}' $clustfile | sort -n | tail -1) ))
        CTceil=$(echo $pctl_low_threshold | awk '{print int($1+1)}')
        if (( $CTceil>$largest )); then let "CTceil=2"; fi; export CTceil
        (( mode=$(awk -v CT=$CTceil "(NF-1)>=CT {print NF-1}" $clustfile |
          sort -n | uniq -c | awk '{print $1 "\t" $2}' | sort -n | tail -1 | awk '{print $2}') ))
        printf "    The pctl${pctl_low} set consists of orthogroups with at least %.0f genes per OG (>= %d * %d/100 sets).\n" \
           $CTceil $max_annot_ct $pctl_low >> ${stats_file} >> ${stats_file}
      fi

      export mode
      (( clusters=$(wc -l < $clustfile) ))
      (( num_at_mode=$(awk -v MODE=$mode '(NF-1)==MODE {ct++} END{print ct}' $clustfile) ))
      (( seqs_clustered=$(awk '{sum+=NF-1} END{print sum}' $clustfile) ))

      printf "  %-20s\t%s\n" "Cluster file" "$clustering.clust.tsv" >> ${stats_file}
      printf "  %-20s\t%s\n" "num_of_clusters" $clusters >> ${stats_file}
      printf "  %-20s\t%s\n" "largest_cluster" $largest >> ${stats_file}
      if (( $clustcount<3 )); then
        printf "  %-20s\t%d\n" "modal_clst_size" $mode >> ${stats_file}
        printf "  %-20s\t%d\n" "num_at_mode" $num_at_mode >> ${stats_file}
      else
        printf "  %-20s\t%d\n" "modal_clst_size>=$CTceil" $mode >> ${stats_file}
        printf "  %-20s\t%d\n" "num_at_mode>=$CTceil" $num_at_mode >> ${stats_file}
      fi
      printf "  %-20s\t%d\n" "seqs_clustered" $seqs_clustered >> ${stats_file}
  
      (( clustcount=$((clustcount+1)) ))
    else
      printf "File $clustfile is not available; skipping\n"
      printf "File $clustfile is not available; skipping\n" >> ${stats_file}
    fi
  done

  echo "  Print sequence composition statistics for each annotation set"
  printf "\n== Sequence stats for CDS files\n" >> ${stats_file}
  printf "  Class:  seqs     min max    N50    ave     annotation_name\n" >> ${stats_file} 
  if [ -f ${WORK_DIR}/stats/tmp.fasta_seqstats ]; then
    cat ${WORK_DIR}/stats/tmp.fasta_seqstats >> ${stats_file}

  printf "\n  Avg:   " >> ${stats_file} 
    cat ${WORK_DIR}/stats/tmp.fasta_seqstats | transpose.pl |
      perl -ane 'BEGIN{use List::Util qw(sum)}; 
                 if ($F[0]=~/^\d+/){
                   $sum=sum @F; $ct=scalar(@F); $avg=$sum/$ct;
                   printf " %4d ", $avg;
                 END{print "   all_annot_sets\n"}
                 }' >> ${stats_file}
  fi

  printf "\n== Sequence stats for final pangene CDS files -- pctl${pctl_low} and trimmed\n" >> ${stats_file}
  printf "  Class:   seqs     min max    N50    ave     annotation_name\n" >> ${stats_file} 
  annot_name=23_syn_pan_pctl${pctl_low}_posn_cds.fna
    printf "  pctl${pctl_low}: " >> ${stats_file}
    cat_or_zcat "${WORK_DIR}/23_syn_pan_pctl${pctl_low}_posn_cds.fna" | calc_seq_stats >> ${stats_file}

  echo "  Print per-annotation-set coverage stats (sequence counts, sequences retained)"
  #   tmp.gene_count_start was generated during run_ingest
  printf "\n== Proportion of initial genes retained in the \"aug_extra\" and \"pctl${pctl_low}\" sets:\n" \
    >> ${stats_file}

  cut -f2 ${WORK_DIR}/18_syn_pan_aug_extra.hsh.tsv | 
    perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' |
    sort | uniq -c | awk '{print $2 "\t" $1}' > ${WORK_DIR}/stats/tmp.gene_count_all_end

  cut -f2 ${WORK_DIR}/22_syn_pan_aug_extra_pctl${pctl_low}.hsh.tsv | 
    perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' |
    sort | uniq -c | awk '{print $2 "\t" $1}' > ${WORK_DIR}/stats/tmp.gene_count_pctl${pctl_low}_end

  paste ${WORK_DIR}/stats/tmp.gene_count_start \
        ${WORK_DIR}/stats/tmp.gene_count_all_end \
        ${WORK_DIR}/stats/tmp.gene_count_pctl${pctl_low}_end | 
    awk 'BEGIN{print "  Start\tEnd_all\tEnd_core\tPct_kept_all\tPct_kept_core\tAnnotation_name"} 
        { printf "  %i\t%i\t%i\t%2.1f\t%2.1f\t%s\n", $2, $4, $6, 100*($4/$2), 100*($6/$2), $1 }'  >> ${stats_file}

  echo "  Print counts per accession"
  if [ -f ${full_out_dir}/18_syn_pan_aug_extra.counts.tsv ]; then
    printf "\n== For all annotation sets, counts of genes-in-orthogroups and counts of orthogroups-with-genes:\n" \
      >> ${stats_file}
    printf "  gns-in-OGs  OGs-w-gns  OGs-w-gns/gns  pct-non-null-OGs  pct-null-OGs  annot-set\n" \
      >> ${stats_file}
    cat ${full_out_dir}/18_syn_pan_aug_extra.counts.tsv | transpose.pl | 
      perl -lane 'next if ($.<=3); 
        $ct=0; $sum=0; $nulls=0; $OGs=0;
        for $i (@F[1..(@F-1)]){
          $OGs++;
          if ($i>0){$ct++; $sum+=$i}
          if ($i==0){$nulls++}
        }; 
        printf("  %d\t%d\t%.2f\t%.2f\t%.2f\t%s\n", $sum, $ct, 100*$ct/$sum, 100*($OGs-$nulls)/$OGs, 100*$nulls/$OGs, $F[0])' \
        >> ${stats_file}
  fi

  echo "  Print histograms"
  if [ -f ${full_out_dir}/06_syn_pan.clust.tsv ]; then
    printf "\nCounts of initial clusters by cluster size, file 06_syn_pan.clust.tsv:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/06_syn_pan.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi

  if [ -f ${full_out_dir}/12_syn_pan_aug.clust.tsv ]; then
    printf "\nCounts of augmented clusters by cluster size, file 12_syn_pan_aug.clust.tsv:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/12_syn_pan_aug.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi

  if [ -f ${full_out_dir}/18_syn_pan_aug_extra.clust.tsv ]; then
    printf "\nCounts of augmented-extra clusters by cluster size, file 18_syn_pan_aug_extra.clust.tsv:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/18_syn_pan_aug_extra.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi
}

##########
run_clean() {
  echo "Clean (delete) files in the working directory that are not needed with subsequent add_extra"
  cd "${WORK_DIR}"
  echo "  work_dir: $PWD"
  if [ -d MMTEMP ]; then rm -rf MMTEMP/*; 
  fi
  for dir in 11_pan_leftovers 13_extra_out_dir 16_pan_leftovers_extra 19_pan_aug_leftover_merged_prot \
    22_syn_pan_aug_extra_pctl${pctl_low}; do
    if [ -d $dir ]; then echo "  Removing directory $dir"; rm -rf $dir &
    fi
  done
  #for file in 10* 11* 14* 20* 21* 23* 24* consen*; do
  for file in 10* 11* 14* 20*; do
    if [ -f $file ]; then echo "  Removing file $file"; rm $file; 
    fi
  done
  wait
  cd $OLDPWD
}

##########
run_ReallyClean() {
  echo "Doing complete clean-up of files in the working directory (remove all files and directories)"
  if [ -d ${WORK_DIR} ]; then
    cd "${WORK_DIR}"
    echo "  work_dir: $PWD"
    if [ -d 01_posn_hsh ]; then
      echo "Expected directory 01_posn_hsh exists in the work_dir (${WORK_DIR}),"
      echo "so proceeding with complete clean-up from that location."
      for dir in *; do
        rm -rf $dir &
      done
    else 
      echo "Expected directory 01_posn_hsh is not present in the work_dir (${WORK_DIR}),"
      echo "so aborting the -K ReallyClean step. Please do this manually if you wish."
      cd $OLDPWD
      exit 1;
    fi
  fi
  wait
  cd $OLDPWD
}

########################################
# Main program

NPROC=$(command -v nproc > /dev/null && nproc || getconf _NPROCESSORS_ONLN)
CONFIG="null"
optarg_work_dir="null"
optarg_order_method="null"
step="all"
retain="no"

export NPROC=${NPROC:-1}
export MMSEQS_NUM_THREADS=${NPROC} # mmseqs otherwise uses all cores by default

# mmseqs uses significant number of threads on its own. Set a maximum, which may be below NPROC.
MMSEQSTHREADS=$(( 4 < ${NPROC} ? 4 : ${NPROC} ))

pandagma_conf_params='clust_iden clust_cov consen_iden extra_iden mcl_inflation 
  strict_synt out_dir_base pctl_low pctl_med pctl_hi 
  consen_prefix annot_str_regex order_method preferred_annot work_dir'

##########
# Command-line interpreter

while getopts "c:w:s:n:o:rvhm" opt
do
  case $opt in
    c) CONFIG=$OPTARG; echo "Config: $CONFIG" ;;
    w) optarg_work_dir=$OPTARG; echo "Work dir: $optarg_work_dir" ;;
    s) step=$OPTARG; echo "step(s): $step" ;;
    n) NPROC=$OPTARG; echo "processors: $NPROC" ;;
    o) optarg_order_method=$OPTARG; echo "order method: $optarg_order_method" ;;
    r) retain="yes" ;;
    v) version ;;
    h) printf >&2 "$HELP_DOC\n" && exit 0 ;;
    m) printf >&2 "$HELP_DOC\n$MORE_INFO\n" && exit 0 ;;
    *) printf >&2 "$HELP_DOC\n" && exit 1 ;;
  esac
done

if [ $CONFIG == "null" ]; then
  printf "\nPlease provide the path to a config file: -c CONFIG\n" >&2
  printf "\nRun \"$scriptname -h\" for help.\n\n" >&2
  exit 1;
else
  export CONF=${CONFIG}
fi

# Add shell variables from config file
. "${CONF}"

if [ "$#" -eq 0 ]; then
  printf >&2 "$HELP_DOC\n" && exit 0;
fi

shift $(expr $OPTIND - 1)

if [ "$order_method" != "$optarg_order_method" ] && [ "$optarg_order_method" != "null" ]; then
  echo "Command-line option for order_method was \"$optarg_order_method\", overriding"
  echo "the setting of \"$order_method\" from the config."
  order_method=$optarg_order_method
fi

##########
# Check for existence of work directory.
# work_dir provided via cli argument takes precedence over the variable specified in the conf file
if [ $optarg_work_dir == "null" ] && [ -z ${work_dir+x} ]; then
  echo "\nPlease provide the path to a work directory, for intermediate files," >&2
  echo "either in the config file (work_dir=) or with option -w" >&2
  printf "\nRun \"$scriptname -h\" for help.\n\n" >&2
  exit 1;
elif [ $optarg_work_dir != "null" ] && [ -d $optarg_work_dir ]; then
  work_dir=$optarg_work_dir
  export WORK_DIR=${work_dir}
elif [ $optarg_work_dir == "null" ] && [ -d $work_dir ]; then
  export WORK_DIR=${work_dir}
elif [ ! -d $work_dir ] && [ ! -d $optarg_work_dir ]; then
  echo "Neither work directory $work_dir nor $optarg_work_dir was not found. Please specify path in the config file or with option -w." >&2
  exit 1;
elif [ ! -d $work_dir ]; then
  echo "Work directory $work_dir was not found. Please specify path in the config file or with option -w." >&2
  exit 1;
else 
  echo "Check odd condition: optarg_work_dir: $optarg_work_dir; work_dir: $work_dir"
  exit 1;
fi

if [[ $step == "clean" ]] || [[ $step == "ReallyClean" ]] ; then
  run_$step
  echo "Command \"$step\" was run for cleanup in the work directory: $WORK_DIR"
  exit 0;
fi

# Get paths to the fasta and annotation files
canonicalize_paths

# Check for existence of third-party executables
missing_req=0
for program in mmseqs dagchainer mcl cons famsa hmmalign hmmbuild run_DAG_chainer.pl; do
  if ! type $program &> /dev/null; then
    echo "Warning: executable $program is not on your PATH."
    missing_req=$((missing_req+1))
  fi
done
if [ "$missing_req" -gt 0 ]; then 
  printf "\nPlease add the programs above to your environment and try again.\n\n"
  exit 1; 
fi

# Check that the bin directory is in the PATH
if ! type hash_into_fasta_id.pl &> /dev/null; then
  printf "\nPlease add the pandagma bin directory to your PATH and try again. Try the following:\n"
  printf "\n  PATH=$PWD/bin:\$PATH\n\n"
  exit 1; 
fi

# Run all specified steps (except clean -- see below; and  ReallyClean, which can be run separately).
# Also, the steps align_cds, align_protein, model_and_trim, calc_trees, and xfr_aligns_trees may be run separately.
commandlist="ingest mmseqs filter dagchainer mcl consense cluster_rest add_extra pick_exemplars \
             filter_to_pctile tabularize order_and_name  calc_chr_pairs summarize"

if [[ $step =~ "all" ]]; then
  for command in $commandlist; do
    echo "RUNNING STEP run_$command"
    run_$command
    echo
  done
else
  run_${step};
fi

if [[ $step =~ "all" ]] && [[ $retain == "yes" ]]; then
  echo "Flag -r (retain) was set, so skipping clean-up of the work directory."
  echo "If you wish to do a cleanupt separately, you can call "
  echo "  .pandagma-pan.sh -c $CONF -s clean";
elif [[ $step =~ "all" ]] && [[ $retain == "no" ]]; then
  echo "Calling the medium cleanup function \"run_clean\" ..."
  run_clean  
fi


#!/bin/bash
#
# Configuration and run script for "pandagma", which generates pan-gene clusters using the programs 
# mmseqs, dagchainer, mcl, and vsearch. These are used to do the initial clustering, 
# synteny-finding, and re-clustering.
# Authors: Steven Cannon, Joel Berendzen, 2020-2021
#
scriptname=`basename "$0"`
version="0.9.0"
set -e # stop on errors
scriptstart=$(date +%s)
pkg="pandagma"
PKG="$(echo ${pkg} | tr /a-z/ /A-Z/)"
PKG_DIR="${PKG}_DIR"
PKG_WORK_DIR="${PKG}_WORK_DIR"
if [ -z "${!PKG_DIR}" ]; then
  root_dir=~/.local/share/${pkg}
else
  root_dir="${!PKG_DIR}"
fi
if [ -z "${!PKG_WORK_DIR}" ]; then
  work_dir="/tmp/${pkg}-work"
else
  work_dir="${!PKG_WORK_DIR}"
fi
if [ -z "${!PKG_DATA_DIR}" ]; then
  data_dir="$PWD/data"
else
  data_dir="${!PKG_DATA_DIR}"
fi
etc_dir="${root_dir}/etc"
dag_dir="${work_dir}/dag"
mmseqs_dir="${work_dir}/mmseqs"
pan_fasta_dir="${work_dir}/pan_fasta"
pan_consen_dir="${work_dir}/pan_consen"
leftovers_dir="${work_dir}/pan_leftovers"
error_exit() {
  echo >&2 "ERROR -- unexpected exit from ${BASH_SOURCE} script at line:"
  echo >&2 "   $BASH_COMMAND"
}
trap error_exit EXIT

TOP_DOC="""Compute pan-gene clusters using the programs mmseqs, dagchainer, and mcl, and
additional pre- and post-refinement steps.

Usage:
        ${pkg} SUBCOMMAND [SUBCOMMAND_OPTIONS]

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
"""
#
# Helper functions begin here
#
set_value() {
  if [ "$2" == "-d" ]; then
    rm -f "${etc_dir}/${1}"
  else
    echo "$2" >"${etc_dir}/${1}"
  fi
}
get_value() {
  if [ -e ${etc_dir}/${1} ]; then
    cat ${etc_dir}/${1}
  else
    trap - EXIT
    echo >&2 "ERROR -- value for $1 variable not found."
    exit 1
  fi
}
#
perl_defs() {
  #
  # Perl local::lib settings, where Bioperl::SeqIO is installed
  #
  perlbase=${SYN_BIN}/perl5
  PATH="${perlbase}/bin${PATH:+:${PATH}}"; export PATH
  PERL5LIB="${perlbase}/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB
  PERL_LOCAL_LIB_ROOT="${perlbase}${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"
  export PERL_LOCAL_LIB_ROOT
  PERL_MB_OPT="--install_base \"${perlbase}\""; export PERL_MB_OPT
  PERL_MM_OPT="INSTALL_BASE=${perlbase}"; export PERL_MM_OPT
}
#
# run functions
#
run_ingest() {
  # Prepare the assembly and annotation files for analysis. Add positional info to gene IDs.
  gff_ext=$(get_value gff_ext)
  fasta_ext=$(get_value fasta_ext)
  gff_files="$(ls ${data_dir}/*.${gff_ext})"
  fasta_files="$(ls ${data_dir}/*.${fasta_ext})"
  n_fasta=`ls ${data_dir}/*.${fasta_ext} | wc -l`
  echo "ingest -- combine data from ${n_fasta} ${gff_ext} and ${fasta_ext} files"
  for path in $gff_files; do
    base=$(basename $path .${gff_ext})
    base_no_ann=$(echo $base | perl -pe 's/\.ann\d+\.\w+//')
    cat $path | awk -v OFS="\t" '$3=="mRNA" {print $1, $4, $5, $9}' |
      perl -pe 's/ID=([^;]+);\S+/$1/' >${work_dir}/${base_no_ann}.bed
  done
  # add positional information to FASTA ids
  for path in ${work_dir}/*.bed; do
    base=$(basename $path .bed)
    cat $path | awk '{print $4 "\t" $1 "__" $4 "__" $2 "__" $3}' \
      >${work_dir}/${base}.hsh
    hash_into_fasta_id.pl\
      -fasta ${data_dir}/${base}.${fasta_ext} \
      -hash ${work_dir}/${base}.hsh \
      -suff_regex \
      >${work_dir}/${base}.${fasta_ext}
  done
  echo
}
#
run_mmseqs() {
  # Do mmseqs clustering on all genome pairings
  mm_clust_iden=$(get_value clust_iden)
  mm_clust_cov=$(get_value clust_cov)
  n_jobs=$(get_value max_main_jobs)
  fasta_ext=$(get_value fasta_ext)
  echo "run mmseqs -- at ${mm_clust_iden} percent identity and minimum of ${mm_clust_cov}% coverage."
  start_time=$(date +%s)
  #
  for qry_path in ${work_dir}/*.${fasta_ext}; do
    qry_base=$(basename $qry_path .${fasta_ext})
    for sbj_path in ${work_dir}/*.${fasta_ext}; do
      sbj_base=$(basename $sbj_path .${fasta_ext})
      if [[ "$qry_base" > "$sbj_base" ]]; then
         echo "Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
         mmseqs easy-cluster $qry_path $sbj_path ${mmseqs_dir}/${qry_base}.x.${sbj_base} tmp \
           --min-seq-id $mm_clust_iden \
           -c $mm_clust_cov \
           --cov-mode 0 \
           --cluster-reassign 2>/dev/null 1>/dev/null &  ## in background, for parallelization  

        # allow to execute up to $n_jobs in parallel
        if [[ $(jobs -r -p | wc -l) -ge $n_jobs ]]; then wait -n; fi
      fi
    done
  done
  wait # wait for last jobs to finish
  end_time=$(date +%s)
  set_value mmseqs_time_s $((end_time-start_time))
}
#
run_filter() {
  echo "From mmseqs cluster output, split out the following fields: molecule, gene, start, stop."
  n_jobs=$(get_value max_main_jobs)
  chr_match_list=${data_dir}/expected_chr_matches.tsv
  if [[ -f ${chr_match_list} ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns from file ${chr_match_list}"
    for mmseqs_path in ${mmseqs_dir}/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      cat $mmseqs_path |
        filter_mmseqs_by_chroms.pl -chr ${chr_match_list} > ${dag_dir}/${outfilebase}_matches.tsv &

      # allow to execute up to $n_jobs in parallel
      if [[ $(jobs -r -p | wc -l) -ge $n_jobs ]]; then wait -n; fi
    done
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches.tsv file was provided, so proceeding without chromosome-pair filtering."
    for mmseqs_path in ${mmseqs_dir}/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      cat $mmseqs_path |
        perl -pe 's/__/\t/g' > ${dag_dir}/${outfilebase}_matches.tsv
    done
  fi
}
#
run_dagchainer() {
  # Identify syntenic blocks, using DAGchainer
  dag_args=$(get_value dagchainer_args)
  echo; echo "Run DAGchainer, using args \"${dag_args}\""
  start_time=$(date +%s)
  for match_path in ${dag_dir}/*_matches.tsv; do
    align_file=`basename $match_path _matches.tsv`
    echo "Running DAGchainer on comparison: $align_file"
    run_DAG_chainer.pl $dag_args \
      -i $match_path 2>/dev/null 1>/dev/null 
  done
  wait
  end_time=$(date +%s)
  set_value dag_time_s $((end_time-start_time))
  #Extract single-linkage synteny anchors
  cat /dev/null > ${work_dir}/synteny_blocks.tsv
  printf "matches\tscore\trev\tid1\tid2\n" >${work_dir}/synteny_blocks.tsv

  cat /dev/null > ${work_dir}/homology_pairs.tsv
  for path in ${dag_dir}/*_matches.tsv; do
    cat $path | awk '$1!~/^#/ {print $2 "\t" $6}' \
      >>${work_dir}/homology_pairs.tsv
  done

  cat /dev/null > ${work_dir}/synteny_pairs.tsv
  for path in ${dag_dir}/*.aligncoords; do
    cat $path | awk '$1!~/^#/ {print $2 "\t" $6}' \
      >>${work_dir}/synteny_pairs.tsv
    cat $path | grep \#\# | grep -v reverse |
      awk '{print substr($14,0,length($14)-2) "\t" $10 "\t" 1 "\t" $3 "\t" $5}' \
        >>${work_dir}/synteny_blocks.tsv
    cat $path | grep \#\# | grep reverse |
      awk '{print substr($15,0,length($15)-2) "\t" $11 "\t" 1 "\t" $3 "\t" $5}' \
        >>${work_dir}/synteny_blocks.tsv
  done
}
#
run_mcl() {
  # Calculate clusters using Markov clustering
  mcl_I=$(get_value mcl_inflation)
  mcl_te=$(get_value mcl_threads)
  prefix=$(get_value pan_prefix)
  printf "\nCalculate clusters. use Markov clustering with inflation parameter $mcl_I and $mcl_te threads\n"
  echo "MCL COMMAND: mcl ${work_dir}/synteny_pairs.tsv -I $mcl_I -te $mcl_te --abc -o ${work_dir}/tmp.syn_pan.clust.tsv"
  mcl ${work_dir}/synteny_pairs.tsv -I $mcl_I -te $mcl_te --abc -o ${work_dir}/tmp.syn_pan.clust.tsv \
    2>/dev/null 1>/dev/null
 
  # Add cluster IDs
  cat ${work_dir}/tmp.syn_pan.clust.tsv |
    awk -v PRE=$prefix '{padnum=sprintf("%05d", NR); print PRE padnum "\t" $0}' > ${work_dir}/syn_pan.clust.tsv

  # Reshape from mcl output form (clustered IDs on one line) to a hash format (clust_ID gene)
  cat ${work_dir}/syn_pan.clust.tsv | 
    perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' > ${work_dir}/syn_pan.hsh.tsv
}
#
run_consense() {
  echo "Calculate a consensus sequence for each pan-gene set, using vsearch."
  echo "Then add previously unclustered sequences into an \"augmented\" pan-gene set, by homology."
  fasta_ext=$(get_value fasta_ext)
  fasta_files="$(ls ${data_dir}/*.${fasta_ext})"
  vs_consen_iden=$(get_value consen_iden)
  mm_clust_iden=$(get_value clust_iden)
  n_jobs=$(get_value max_lite_jobs)

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  get_fasta_from_family_file.pl ${fasta_files} -fam ${work_dir}/syn_pan.clust.tsv -out ${pan_fasta_dir} 

  echo "  Calculate consensus sequences for each pan-gene set."
  ls ${pan_fasta_dir} | xargs -I{} -n 1 -P $n_jobs \
    vsearch --cluster_fast ${pan_fasta_dir}/{} --id ${vs_consen_iden} --fasta_width 0 \
            --consout ${pan_consen_dir}/{} 2>/dev/null 1>/dev/null 

  echo "  Combine consensus sequnces into one multifasta file"
  cat /dev/null > ${work_dir}/syn_pan_consen.fna
  for path in ${pan_consen_dir}/*; do
    file=`basename $path`;
    awk -v FN=$file '$1~/^>/ {print ">" FN " " substr($1,2)} 
                     $1!~/^>/ {print $1}' $path >> ${work_dir}/syn_pan_consen.fna
  done 
  rm ${pan_consen_dir}/*

  echo "  Get sorted list of all genes, from the original fasta files"
  cat ${fasta_files} | awk '$1~/^>/ {print $1}' | sed 's/>//' | sort > ${work_dir}/lis.all_genes

  echo "  Get sorted list of all clustered genes"
  cat ${pan_fasta_dir}/* | awk '$1~/^>/ {print $1}' | sed 's/>//' | sort > ${work_dir}/lis.all_clustered_genes

  echo "  Get list of genes not in clusters"
  comm -13 ${work_dir}/lis.all_clustered_genes ${work_dir}/lis.all_genes > ${work_dir}/lis.genes_not_in_clusters

  echo "  Retrieve the non-clustered genes"
  cat ${fasta_files} > ${work_dir}/all_genes.fna
  get_fasta_subset.pl -in ${work_dir}/all_genes.fna \
                      -out ${work_dir}/genes_not_in_clusters.fna \
                      -lis ${work_dir}/lis.genes_not_in_clusters -clobber

  echo "  Search non-clustered genes against pan-gene consensus sequences"
  mmseqs easy-search ${work_dir}/genes_not_in_clusters.fna \
                     ${work_dir}/syn_pan_consen.fna \
                     ${work_dir}/unclust.x.all_cons.m8 tmp --search-type 3 --cov-mode 5 -c 0.5

  echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
  top_line.awk ${work_dir}/unclust.x.all_cons.m8 | 
    awk -v IDEN=${mm_clust_iden} '$3>=IDEN {print $2 "\t" $1}' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk >  ${work_dir}/syn_pan_leftovers.clust.tsv

  echo "  Retrieve sequences for the "leftover" genes"
  get_fasta_from_family_file.pl ${fasta_files} \
    -fam ${work_dir}/syn_pan_leftovers.clust.tsv -out ${leftovers_dir}

  echo "  Make augmented cluster sets"
  cat /dev/null > ${work_dir}/syn_pan_augmented.clust.tsv
  for path in ${work_dir}/pan_fasta/*; do
    file=`basename $path` 
    if [[ -f ${leftovers_dir}/$file ]]; then
      cat <(awk -v ORS="" -v ID=$file 'BEGIN{print ID "\t"} $1~/^>/ {print substr($1,2) "\t"}' $path) \
          <(awk -v ORS="" '$1~/^>/ {print substr($1,2) "\t"} END{print "\n"}' ${leftovers_dir}/$file)
    else
      awk -v ORS=""  -v ID=$file 'BEGIN{print ID "\t"} $1~/^>/ {print substr($1,2) "\t"} END{print "\n"}' $path
    fi | sed 's/\t$//'
  done > ${work_dir}/syn_pan_augmented.clust.tsv

  echo "  Reshape from mcl output form (clustered IDs on one line) to a hash format (clust_ID gene)"
  cat ${work_dir}/syn_pan_augmented.clust.tsv | 
    perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' > ${work_dir}/syn_pan_augmented.hsh.tsv
}
#
run_summarize() {
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"
  mcl_I=$(get_value mcl_inflation)
  mm_clust_iden=$(get_value clust_iden)
  mm_clust_cov=$(get_value clust_cov)
  mmseqs_time=$(get_value mmseqs_time_s)
  dag_time=$(get_value dag_time_s)
  out_dir="$(get_value out_dir_base)"
  full_out_dir=`echo "$out_dir.id${mm_clust_iden}.cov${mm_clust_cov}.I${mcl_I}" | perl -pe 's/(\d)\.(\d+)/$1_$2/g'`
  stats_file=${full_out_dir}/stats.tsv

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  cp ${work_dir}/synteny_blocks.tsv ${full_out_dir}/synteny_blocks.tsv
  cp ${work_dir}/syn_pan.clust.tsv ${full_out_dir}/syn_pan.clust.tsv
  cp ${work_dir}/syn_pan.hsh.tsv ${full_out_dir}/syn_pan.hsh.tsv
  cp ${work_dir}/syn_pan_augmented.clust.tsv ${full_out_dir}/syn_pan_augmented.clust.tsv
  cp ${work_dir}/syn_pan_augmented.hsh.tsv ${full_out_dir}/syn_pan_augmented.hsh.tsv

  printf "Run of program $scriptname, version $version\n\n" > ${stats_file}

  printf "Parameter  \tvalue\n" >> ${stats_file}
  for key in $(ls ${etc_dir}); do
    value="$(get_value ${key})"
    if [[ ${key} != +(dag|mmseqs)_time_s ]]; then
      printf '%-15s\t%s\n' ${key} ${value} >> ${stats_file}
    fi
  done

  printf "\nOutput directory for this run:\t${full_out_dir}\n" >> ${stats_file}

  printf '%-20s\t%s\n' "Statistic" "value" >> ${stats_file}

  let "n_blocks=$(cat ${full_out_dir}/synteny_blocks.tsv | wc -l)-1"
  printf '%-20s\t%s\n' $synteny_blocks $n_blocks >> ${stats_file}

#
  printf "\n== Initial clusters (containing only genes within synteny blocks)\n" >> ${stats_file}
  let "clusters=$(cat ${full_out_dir}/syn_pan.clust.tsv | wc -l)"
  printf '%-20s\t%s\n' "num_of_clusters" $clusters >> ${stats_file}

  let "largest=$(cat ${full_out_dir}/syn_pan.clust.tsv | awk "{print NF-1}" | head -1)"
  printf '%-20s\t%s\n' "largest_cluster" $largest >> ${stats_file}

  let "mode=$(cat ${full_out_dir}/syn_pan.clust.tsv | awk "{print NF-1}" | \
    uniq -c | sort -n | tail -1 | awk '{print $2}')"
  printf '%-20s\t%s\n' "modal_clst_size" $mode >> ${stats_file}

  let "num_at_mode=$(cat ${full_out_dir}/syn_pan.clust.tsv | awk "{print NF-1}" | \
    uniq -c | sort -n | tail -1 | awk '{print $1}')"
  printf '%-20s\t%s\n' "num_at_mode" $num_at_mode >> ${stats_file}
  
#
  printf "\n== Augmented clusters (unanchored sequences added to the initial clusters)\n" >> ${stats_file}
  let "clustersA=$(cat ${full_out_dir}/syn_pan_augmented.clust.tsv | wc -l)"
  printf '%-20s\t%s\n' "num_of_clusters" $clustersA >> ${stats_file}

  let "largestA=$(cat ${full_out_dir}/syn_pan_augmented.clust.tsv | awk "{print NF-1}" | sort -n | tail -1)"
  printf '%-20s\t%s\n' "largest_cluster" $largestA >> ${stats_file}

  let "modeA=$(cat ${full_out_dir}/syn_pan_augmented.clust.tsv | awk "{print NF-1}" | \
    uniq -c | sort -n | tail -1 | awk '{print $2}')"
  printf '%-20s\t%s\n' "modal_clst_size" $modeA >> ${stats_file}

  let "numA_at_mode=$(cat ${full_out_dir}/syn_pan_augmented.clust.tsv | awk "{print NF-1}" | \
    sort -n | uniq -c | sort -n | tail -1 | awk '{print $1}')"
  printf '%-20s\t%s\n' "num_at_mode" $numA_at_mode >> ${stats_file}

  # # run times
  # printf "\nRun times (seconds):" >> ${stats_file}
  # printf "\nmmseqs\t${mmseqs_time}" >> ${stats_file}
  # printf "\nDAGchainer\t${dag_time}\n" >> ${stats_file}

  # histograms
  printf "\nCounts of initial clusters by cluster size:\n" >> ${stats_file}
  cut -f1 ${full_out_dir}/syn_pan.hsh.tsv | sort | uniq -c |
    awk '{print $1}' | uniq -c | awk '{print $2 "\t" $1}' | sort -k1n,1n >> ${stats_file}

  printf "\nCounts of augmented clusters by cluster size:\n" >> ${stats_file}
  cut -f1 ${full_out_dir}/syn_pan_augmented.hsh.tsv | sort | uniq -c |
    awk '{print $1}' | sort -n | uniq -c | awk '{print $2 "\t" $1}' | sort -k1n,1n >> ${stats_file}

  echo
  cat ${stats_file}
}
#
# top-level command functions
#
config() {
  CONFIG_DOC="""Sets/displays key/value pairs for the $pkg build system.

Usage:
   $scriptname config [-h] [KEY] [VALUE | -d]

Options:
   -h   Prints this help message and exits
   -d   Deletes the setting of KEY

Arguments:
   If KEY is absent, all values will be displayed
   If KEY is present but VALUE is absent, the value will be displayed
   If KEY and VALUE are present, the value will be set
"""
  if [ "$#" -eq 0 ]; then
    echo >&2 "$CONFIG_DOC"
    param="all"
  elif [ "$1" == "-h" ]; then
    echo >&2 "$CONFIG_DOC"
    param=all
  else
    param="$1"
  fi
  if [ "$param" == "all" ]; then
      trap - EXIT
      for key in $(ls ${etc_dir}); do
        value="$(get_value ${key})"
        if [[ ${key} != +(dag|mmseqs)_time_s ]]; then
          printf '%-20s\t%s\n' ${key} ${value} >&1
        fi
      done
      exit 0
  fi
  if [ "$#" -eq 1 ]; then
    if [ -e ${etc_dir}/${param} ]; then
      echo "$(get_value $param)"
    else
      trap - EXIT
      echo >&2 "ERROR -- \"${1}\" has not been set"
      exit 1
    fi
  elif [ "$#" -eq 2 ]; then # set
    set_value $param $2
  else
    trap - EXIT
    echo >&2 "$CONFIG_DOC"
    echo >&2 "ERROR -- too many arguments (${#})."
    exit 1
  fi
}
#
clear_config() {
  echo "clearing configuration directory"
  rm -f ${etc_dir}/*
}
#
init() {
  # Initialize parameters required for run. Write these to files, for persistent access through the program.
  echo "setting run configuration parameters"
  [ -z "$NPROC" ] && NPROC=$(nproc)
  set_value clust_iden "0.98"
  set_value clust_cov "0.75"
  set_value consen_iden "0.80"
  set_value max_main_jobs $(($NPROC/5))
  set_value max_litejobs $(($NPROC/2))
  set_value mcl_inflation 2
  set_value mcl_threads $(($NPROC/5))
  set_value dagchainer_args ""
  set_value fasta_ext "fna"
  set_value gff_ext "gff3"
  set_value pan_prefix "pan"
  set_value out_dir_base "${pkg}_out"
  config all

  echo "Data directory (for assembly and GFF files): $data_dir"
  echo "Work directory (for temporary files): $work_dir"; echo
}
#
run() {
  RUN_DOC="""Run an analysis step

Usage:
   $scriptname run [STEP]

Steps:
   If STEP is not set, the following steps will be run in order,
   otherwise the step is run by itself:
                init - initialize parameters required for run
              ingest - get info from matching GFF and FNA files
              mmseqs - run mmseqs for all gene sets
              filter - select gene matches from indicated chromosome pairings
          dagchainer - compute Directed Acyclic Graphs
                 mcl - calculate Markov clusters
            consense - calculate a consensus sequences from each pan-gene set, 
                         If possible add sequences missed in the first clustering round.
           summarize - compute synteny stats
"""
  commandlist="init ingest mmseqs filter dagchainer mcl consense summarize"
  if [ "$#" -eq 0 ]; then # run the whole list
    for package in $commandlist; do
      echo "RUNNING PACKAGE run_$package"
      run_$package
      echo
    done
  else
    command="$1"
    shift 1
    case $commandlist in
    *"$command"*)
      run_$command $@
      ;;
    $commandlist)
      trap - EXIT
      echo >&2 "$RUN_DOC"
      if [ "$command" == "-h" ]; then
        exit 0
      fi
      echo >&2 "ERROR -- unrecognized run step \"$1\""
      exit 1
      ;;
    esac
  fi
}
#
version() {
  echo $scriptname $version
}
#
clean() {
  echo "cleaning work directory and tmpOut files"
  rm -rf $work_dir 
  rm -f .*.tmpOut
}
#
# Command-line interpreter.
#
if [ "$#" -eq 0 ]; then
  trap - EXIT
  echo >&2 "$TOP_DOC"
  exit 0
fi
#
#
# Create directories if needed
dirlist="root_dir work_dir etc_dir mmseqs_dir dag_dir pan_fasta_dir pan_consen_dir leftovers_dir"
for dirvar in $dirlist; do
    dirname="${!dirvar}"
    if [ ! -d "$dirname" ]; then
      echo "creating directory \"${dirname}\" as $dirvar"
      mkdir -p $dirname
    fi
done
#
command="$1"
shift 1
case $command in
"config")
  config $@
  ;;
"clean")
  clean $@
  ;;
"clear_config")
  clear_config $@
  ;;
"init")
  init $@
  ;;
"run")
  run $@
  ;;
"version")
  version $@
  ;;
*)
  trap - EXIT
  echo >&2 "ERROR -- command \"$command\" not recognized."
  exit 1
  ;;
esac
trap - EXIT
exit 0

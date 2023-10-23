#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates gene-family 
# orthogroups using the programs mmseqs, dagchainer, and mcl. 
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2023
#
scriptname=`basename "$0"`
version="2023-10-18"
set -o errexit -o errtrace -o nounset -o pipefail

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

HELP_DOC="""Place annotation sets (CDS and protein) into pangene or gene family sets, using hmmsearch 
to compare the indicated annotations against HMMs calculated in a prior run of pandagma-pan or pandagma-fam.

Usage:
  ./$scriptname -c CONFIG_FILE [options]

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

Environment requirements: The following packages need to be available in your PATH:
    hmmsearch hmmalign fasttree

Also, please add the pandagma utility programs in the bin directory adjacent to pandagma-fam.sh, e.g.
    PATH=$PWD/bin:\$PATH

Subcommands (in order they are usually run):
                all - All of the steps below, except for ks_filter, clean and ReallyClean
                        (Or equivalently: omit the -s flag; \"all\" is default).
             ingest - Prepare the assembly and annotation files for analysis.
             mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies.
             filter - Filter the synteny results for chromosome pairings, returning gene pairs.
         dagchainer - Run DAGchainer to filter for syntenic blocks.
            ks_calc - Calculation of Ks values on gene pairs from DAGchainer output.
          ks_filter - Filtering based on provided ks_peaks.tsv file (assumes prior ks_calc step)
                mcl - Derive clusters, with Markov clustering.
           consense - Calculate a consensus sequences from each pan-gene set, 
                      adding sequences missed in the first clustering round.
       cluster_rest - Retrieve unclustered sequences and cluster those that can be.
          add_extra - Add other gene model sets to the primary clusters. Useful for adding
                      annotation sets that may be of lower or uncertain quality.
              align - Align families.
     model_and_trim - Build HMMs and trim the alignments, preparatory to calculating trees.
         calc_trees - Calculate gene trees.
          summarize - Move results into output directory, and report summary statistics.

  Run the following subcommands separately if you wish:
              clean - Clean (delete) files in the working directory that are not needed 
                        for later addition of data using add_extra and subsequent run commands.
                        By default, \"clean\" is run as part of \"all\" unless the -r flag is set.
        ReallyClean - Do complete clean-up of files in the working directory.
                        Use this if you want to start over, OR if you are satisified with the results and
                        don't anticipate adding other annotation sets to this pan-gene set.
"""

MORE_INFO="""
Variables in pandagma config file (Set the config with the CONF environment variable)
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.60]
        consen_iden - Minimum identity threshold for consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of \"extra\" annotations [90]
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]

      consen_prefix - Prefix to use in orthogroup names
       out_dir_base - Base name for the output directory [default: './out']
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g. 
                        \"([^.]+\.[^.]+)\..+\"
                          for two dot-separated fields, e.g. vigan.Shumari
                        or \"(\D+\d+\D+)\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    preferred_annot - String to match and select an annotation set, from a gene ID.
                        This is used for picking representative IDs+sequence from an orthogroup, when
                        this annotation is among those with the median length for the orthogroup.
                        Otherwise, one is selected at random from those with median length.
           work_dir - Working directory, for temporary and intermediate files. 

File sets (arrays):
      protein_files
          cds_files
 """

########################################
# Helper functions begin here

version() {
  echo $scriptname $version
}

##########
canonicalize_paths() {
  echo "Entering canonicalize_paths."

  protein_files=($(realpath --canonicalize-existing "${protein_files[@]}"))
  cds_files=($(realpath --canonicalize-existing "${cds_files[@]}"))

  readonly submit_dir=${PWD}

  prot_file=$(basename "${protein_files[0]}" .gz)
  faa="${prot_file##*.}"

  cds_file=$(basename "${cds_files[0]}" .gz)
  fna="${cds_file##*.}"

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
# Retrieve nucleotide and protein data sets
  cd "${WORK_DIR}"
  echo; echo "Run ingest: pull nucleotide and protein data sets into the work directory."
  echo "Note that this is a simpler ingest than for the primary gene family or pangene construction,"
  echo "as the identifiers aren't modified with positional information."
  
  mkdir -p 02_fasta_nuc_add 02_fasta_prot_add 

    # Prepare the tmp.gene_count_start to be joined, in run_summarize, with tmp.gene_count_end_pctl??_end.
    # This is captured from the gene IDs using the annot_str_regex set in the config file.
    cat /dev/null > stats/tmp.gene_count_start_add
    cat /dev/null > stats/tmp.fasta_seqstats_add
    start_time=`date`
    printf "Run started at: $start_time\n" > stats/tmp.timing

  export ANN_REX=${annot_str_regex}

  echo "  Pull the protein files locally"
  cat /dev/null > 02_all_main_prot_add.faa # Collect all starting protein sequences, for later comparisons
  for (( file_num = 0; file_num < ${#protein_files[@]} ; file_num++ )); do
    file_base=$(basename ${protein_files[file_num]%.*})
    zcat "${protein_files[file_num]}" >> 02_all_main_prot_add.faa # Collect original seqs for later comparisons
    zcat "${protein_files[file_num]}" > 02_fasta_prot_add/$file_base
    # calc basic sequence stats
    annot_name=$(basename 02_fasta_prot/$file_base | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
    printf "  Added with hmmsearch:  " >> stats/tmp.fasta_seqstats_add
    cat_or_zcat "${protein_files[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats_add
  done

  echo "  Pull the CDS files locally"
  cat /dev/null > 02_all_main_cds_add.fna # Collect all starting cds sequences, for later comparisons
  for (( file_num = 0; file_num < ${#cds_files[@]} ; file_num++ )); do
    file_base=$(basename ${cds_files[file_num]%.*})
    zcat "${cds_files[file_num]}" >> 02_all_main_cds_add.fna # Collect original seqs for later comparisons
    zcat "${cds_files[file_num]}" > 02_fasta_nuc_add/$file_base
  done

  echo "  Count starting sequences, for later comparisons"
  for file in 02_fasta_prot_add/*.$faa; do
    awk '$0~/UNDEFINED/ {ct++} 
      END{if (ct>0){print "Warning: " FILENAME " has " ct " genes without position (HASH UNDEFINED)" } }' $file
    cat $file | grep '>' | perl -pe 's/__/\t/g' | cut -f2 | # extracts the GeneName from combined genspchr__GeneName__start__end__orient
      perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex.+/$1/' |
      grep -v UNDEFINED | sort | uniq -c | awk '{print $2 "\t" $1}' >> stats/tmp.gene_count_start_add
  done

  sort -o stats/tmp.gene_count_start_add stats/tmp.gene_count_start_add
}

run_fam_consen() {
  echo
  echo; echo "== Generate a consensus sequence for each family =="
  cd "${WORK_DIR}"

  if [ $consen_method == "hmmemit" ]; then
    if [ -d 21_hmmemit ]; then rm -rf 21_hmmemit; fi
    mkdir -p 21_hmmemit
  
    MINL=0.25 # show consensus as 'any' (X/N) unless >= this fraction
    MINU=0.5  # show consensus as upper case if >= this fraction
    for filepath in 21_hmm/* ; do
      base=`basename $filepath`
      hmmemit -C --minl $MINL --minu $MINU -o 21_hmmemit/$base $filepath &
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
    done
    wait
  
    cat 21_hmmemit/* | perl -pe 's/-consensus//' > 30_consen.faa

  elif [ $consen_method == "cons" ]; then
    if [ -d 30_hmmalign_cons ]; then rm -rf 30_hmmalign_cons; fi
    mkdir -p 30_hmmalign_cons
  
    for filepath in 23_hmmalign_trim2/* ; do
      base=`basename $filepath`
      cons -sequence $filepath -outseq 30_hmmalign_cons/$base -name $base &
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
    done
    wait
  
    cat 30_hmmalign_cons/* > 30_consen.faa

  else 
    echo "Unrecognized consen_method: $consen_method."
    echo "Expected values are \"hmmemit\" or \"cons\"."
    exit
  fi
}

##########
run_search_consen() {
  cd "${WORK_DIR}"
  echo "Search provided annotation sets (protein) against family consensus sequences (from hmmemit)"

  if [ -d 33_mmseqs_fam_match ]; then rm -rf 33_mmseqs_fam_match; fi
  mkdir -p 33_mmseqs_fam_match 33_mmseqs_tmp

  SEQTYPE=1 # 3=nuc; 1=pep

  for filepath in 02_fasta_prot_add/*; do 
    base=`basename $filepath .$faa`
    echo "  Search $base.$faa against family hmmemit consensus sequences"
    mmseqs easy-search $filepath \
                       30_consen.faa \
                       33_mmseqs_fam_match/$base.x.consen.m8 \
                       33_mmseqs_tmp \
                       --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null
  done
}

##########
run_place_matches() {
  echo; echo "== Place sequences into families based on top mmsearch hits."
  echo "  Use the consen_iden threshold of $consen_iden."

  cd "${WORK_DIR}"
  export FAM_PRE=${consen_prefix}
  cat 33_mmseqs_fam_match/*.m8 | top_line.awk |
    awk -v IDEN=${consen_iden} '$3>=IDEN {print $2 "\t" $1}' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > 34_add_vs_fam_consen.clust.tsv
}

##########
run_realign_and_trim() {
  echo; echo "== Retrieve sequences for each family, preparatory to aligning them =="
  cd "${WORK_DIR}"
  mkdir -p 35_add_in_fams_prot
  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  get_fasta_from_family_file.pl "${protein_files[@]}" \
    -fam 34_add_vs_fam_consen.clust.tsv -out 35_add_in_fams_prot

  echo; echo "== Realign to HMMs =="
  mkdir -p 42_hmmalign
  for filepath in 35_add_in_fams_prot/*; do 
    file=`basename $filepath`;
    hmmalign --trim --outformat A2M --amino -o 42_hmmalign/$file 21_hmm/$file 35_add_in_fams_prot/$file &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Trim HMM alignments to match-states =="
  mkdir -p 43_hmmalign_trim1
  for filepath in 42_hmmalign/*; do 
    file=`basename $filepath`;
    cat $filepath | 
      perl -ne 'if ($_ =~ />/) {print $_} else {$line = $_; $line =~ s/[a-z]//g; print $line}' |
      sed '/^$/d' > 43_hmmalign_trim1/$file &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Filter alignments prior to tree calculation =="
  mkdir -p 43_hmmalign_trim2 43_hmmalign_trim2_log
  min_depth=3
  min_pct_depth=20
  min_pct_aligned=20
  for filepath in 43_hmmalign_trim1/*; do 
    file=`basename $filepath`
    filter_align.pl -in $filepath -out 43_hmmalign_trim2/$file -log 43_hmmalign_trim2_log/$file \
                    -depth $min_depth -pct_depth $min_pct_depth -min_pct_aligned $min_pct_aligned &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_calc_trees() {
  echo; echo "== Calculate trees =="
  cd "${WORK_DIR}"
  
  mkdir -p 44_trees

  echo; echo "== Move small (<4) and very low-entropy families (sequences are all identical) to the side =="
  mkdir -p 43_pan_aug_small_or_identical
  min_seq_count=4
  # Below, "count" is the number of unique sequences in the alignment.
  for filepath in 43_hmmalign_trim2/*; do
    file=`basename $filepath`
    count=$(awk '$1!~/>/ {print FILENAME "\t" $1}' $filepath | sort -u | wc -l);
    if [[ $count -lt $min_seq_count ]]; then
      echo "Set aside small or low-entropy family $file";
      mv $filepath 43_pan_aug_small_or_identical/
    fi;
  done

  # By default, FastTreeMP uses all available threads.
  # It is more efficient to run more jobs on one core each by setting an environment variable.
  OMP_NUM_THREADS=1
  export OMP_NUM_THREADS
  for filepath in 43_hmmalign_trim2/*; do
    file=`basename $filepath`
    echo "  Calculating tree for $file"
    fasttree -quiet $filepath > 44_trees/$file &
  done
  wait
}

##########
run_summarize() {
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"

  WORK_DIR=$work_dir
  cd "${WORK_DIR}"
  echo "  work_dir: $PWD"

  echo "  Calculate matrix of gene counts per orthogroup and annotation set"
  calc_pan_stats.pl -annot_regex $ANN_REX -pan 18_syn_pan_aug_extra.clust.tsv -out 18_syn_pan_aug_extra.counts.tsv
  max_annot_ct=$(cat 18_syn_pan_aug_extra.counts.tsv | 
                       awk '$1!~/^#/ {print $2}' | sort -n | uniq | tail -1)
 
  conf_base=`basename $CONF .conf`
  full_out_dir="${out_dir_base}_$conf_base"
  stats_file=${full_out_dir}/stats.$conf_base.txt
  export ANN_REX=${annot_str_regex}

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  for file in 06_syn_pan.clust.tsv 06_syn_pan_ge3.hsh.tsv \
              12_syn_pan_aug.clust.tsv 12_syn_pan_aug.hsh.tsv \
              18_syn_pan_aug_extra.clust.tsv  18_syn_pan_aug_extra.hsh.tsv 18_syn_pan_aug_extra.counts.tsv; do
    if [ -f ${WORK_DIR}/$file ]; then
      cp ${WORK_DIR}/$file ${full_out_dir}/
    else 
      echo "Warning: couldn't find file ${WORK_DIR}/$file; skipping"
    fi
  done

  for dir in 19_pan_aug_leftover_merged_prot 21_hmm 42_hmmalign 23_hmmalign_trim2 24_trees; do
    if [ -d ${WORK_DIR}/$dir ]; then
      echo "Copying directory $dir to output directory"
      cp -r ${WORK_DIR}/$dir ${full_out_dir}/
    else 
      echo "Warning: couldn't find dir ${WORK_DIR}/$dir; skipping"
    fi
  done

  printf "Run of program $scriptname, version $version\n" > ${stats_file}

  end_time=`date`
  cat ${WORK_DIR}/stats/tmp.timing >> ${stats_file}
  printf "Run ended at:   $end_time\n\n" >> ${stats_file}

  echo "  Report parameters from config file"
  printf "Parameter  \tvalue\n" >> ${stats_file}
  for key in ${pandagma_conf_params}; do
    printf '%-15s\t%s\n' ${key} "${!key}" >> ${stats_file}
  done

  printf "\nOutput directory for this run:\t${full_out_dir}\n" >> ${stats_file}

  echo "  Report orthogroup composition statistics for the three main cluster-calculation steps"

  echo "  Print sequence composition statistics for each annotation set"
  printf "\n== Sequence stats for protein files\n" >> ${stats_file}
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

  echo "  Print per-annotation-set coverage stats (sequence counts, sequences retained)"
  #   tmp.gene_count_start was generated during run_ingest
  printf "\n== Proportion of initial genes retained in the \"aug_extra\" set:\n" \
    >> ${stats_file}

  cut -f2 ${WORK_DIR}/18_syn_pan_aug_extra.hsh.tsv | 
    perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' |
    sort | uniq -c | awk '{print $2 "\t" $1}' > ${WORK_DIR}/stats/tmp.gene_count_all_end

  paste ${WORK_DIR}/stats/tmp.gene_count_start \
        ${WORK_DIR}/stats/tmp.gene_count_all_end |
    awk 'BEGIN{print "  Start\tEnd\tPct_kept\tAnnotation_name"} 
        { printf "  %i\t%i\t%2.1f\t%s\n", $2, $4, 100*($4/$2), $1 }'  >> ${stats_file}

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
  for dir in 11_pan_leftovers 13_extra_out_dir 16_pan_leftovers_extra 19_pan_aug_leftover_merged_prot; do
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
step="all"
retain="no"

export NPROC=${NPROC:-1}

pandagma_conf_params='consen_iden out_dir_base pctl_low pctl_med pctl_hi consen_prefix annot_str_regex work_dir'

##########
# Command-line interpreter

while getopts "c:w:s:n:o:rvhm" opt
do
  case $opt in
    c) CONFIG=$OPTARG; echo "Config: $CONFIG" ;;
    w) optarg_work_dir=$OPTARG; echo "Work dir: $optarg_work_dir" ;;
    s) step=$OPTARG; echo "step(s): $step" ;;
    n) NPROC=$OPTARG; echo "processors: $NPROC" ;;
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
for program in hmmscan famsa; do
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
commandlist="ingest fam_consen search_consen place_matches realign_and_trim calc_trees summarize"

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
  echo "  .pandagma-fam.sh -c $CONF -s clean";
elif [[ $step =~ "all" ]] && [[ $retain == "no" ]]; then
  echo "Calling the medium cleanup function \"run_clean\" ..."
  run_clean  
fi


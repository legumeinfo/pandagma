#!/usr/bin/env bash

version="2023-02-23"
set -o errexit -o errtrace -o nounset -o pipefail -o posix

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

# to help assign a heredoc value to a variable. The internal line return is intentional.
define(){ o=; while IFS=$'\n' read -r a; do o="$o$a"'
'; done; eval "$1=\$o"; }

declare scriptname annot_name fam_dir annot_str_regex dependencies commandlist

version() {
  echo "$scriptname" $version
}

cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}

check_seq_type () {
  someseq=${1}
  proportion_nuc=$(echo "$someseq" | fold -w1 | 
    awk '$1~/[ATCGN]/ {nuc++} $1!~/[ATCGN]/ {not++} END{print nuc/(nuc+not)}')
  export proportion_nuc
  perl -le '$PN=$ENV{"proportion_nuc"}; if ($PN>0.9){print 3} else {print 1}'
}

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
   | awk -v ANN="$annot_name" 'BEGIN { min=99999999999999 }
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


##########
run_align_cds() {
  cd "${WORK_DIR}" || exit

  echo "== Align CDS sequences =="

  # If not already present, retrieve sequences for each family, preparatory to aligning them
  # 19_pan_aug_leftover_merged_cds
  if [[ -d 19_palmc ]]; then
    : # do nothing; the directory and file(s) exist
  else
    mkdir -p 19_palmc
    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    get_fasta_from_family_file.pl "${cds_files[@]}" "${cds_files_extra_constr[@]}" "${cds_files_extra_free[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_palmc
  fi

  echo; echo "== Move small families to the side =="
  mkdir -p 19_palmc_small
  # Below, count_seqs = number of sequences in the alignment; count_annots = number of unique annotation groups.
  min_annots_in_align=2 # require at least this many distinct annotation groups in an alignment to retain it.
  for filepath in 19_palmc/*; do
    file=$(basename "$filepath")
    count_seqs=$(awk '$1!~/>/ {print FILENAME "\t" $1}' "$filepath" | wc -l);
    count_annots=$( perl -lane 'BEGIN{$REX=qr($ENV{"ANN_REX"})}; if($F[0]=~/$REX/){ print $1 }' $filepath | sort -u | wc -l)
    if [[ $count_seqs -lt $min_align_count ]] || [[ $count_annots -lt $min_annots_in_align ]]; then
      echo "Set aside small family $file; $count_seqs sequences, $count_annots distinct annotations";
      mv "$filepath" 19_palmc_small/
    fi;
  done

  echo; echo "== Align nucleotide sequence for each gene family =="
  mkdir -p 20_aligns_cds
  for filepath in 19_palmc/*; do
    file=$(basename "$filepath");
    echo "  Computing alignment, using program famsa, for file $file"
    famsa -t 2 19_palmc/"$file" 20_aligns_cds/"$file" 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_align_protein() {
  cd "${WORK_DIR}" || exit

  echo "== Align protein sequences =="

  # If not already present, retrieve sequences for each family, preparatory to aligning them
  if [[ -d 19_palmp ]]; then
    : # do nothing; the directory and file(s) exist
  else
    mkdir -p 19_palmp
    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    get_fasta_from_family_file.pl "${protein_files[@]}" "${protein_files_extra[@]}" \
      "${protein_files_extra_constr[@]}" "${protein_files_extra_free[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_palmp
  fi

  echo; echo "== Move small families to the side =="
  mkdir -p 19_palmp_small
  # Below, count_seqs = number of unique sequences in the alignment; count_annots = # of unique annotation groups.
  min_annots_in_align=2 # require at least this many distinct annotation groups in an alignment to retain it.
  for filepath in 19_palmp/*; do
    file=$(basename "$filepath")
    count_seqs=$(awk '$1!~/>/ {print FILENAME "\t" $1}' "$filepath" | sort -u | wc -l);
    count_annots=$( perl -lane 'BEGIN{$REX=qr($ENV{"ANN_REX"})}; if($F[0]=~/$REX/){ print $1 }' $filepath | sort -u | wc -l)
    if [[ $count_seqs -lt $min_align_count ]] || [[ $count_annots -lt $min_annots_in_align ]]; then
      echo "Set aside small family $file; $count_seqs sequences, $count_annots annotations";
      mv "$filepath" 19_palmp_small/
    fi;
  done

  echo; echo "== Align protein sequence for the each gene family =="
  mkdir -p 20_aligns_prot
  for filepath in 19_palmp/*; do
    file=$(basename "$filepath");
    # echo "  Computing alignment, using program famsa, for file $file"
    famsa -t 2 19_palmp/"$file" 20_aligns_prot/"$file" 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_model_and_trim() {
  echo; echo "== Build HMMs =="
  cd "${WORK_DIR}" || exit
  mkdir -p 21_hmm
  for filepath in 20_aligns_prot/*; do
    file=$(basename "$filepath");
    hmmbuild -n "$file" 21_hmm/"$file" "$filepath" 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Realign to HMMs =="
  mkdir -p 22_hmmalign
  for filepath in 21_hmm/*; do
    file=$(basename "$filepath");
    # printf "%s " "$file"
    hmmalign --trim --outformat A2M -o 22_hmmalign/"$file" 21_hmm/"$file" 19_palmp/"$file" &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo

  echo; echo "== Trim HMM alignments to match-states =="
  mkdir -p 23_hmmalign_trim1
  for filepath in 22_hmmalign/*; do
    file=$(basename "$filepath");
    # printf "%s " "$file"
    perl -ne 'if ($_ =~ />/) {print $_} else {$line = $_; $line =~ s/[a-z]//g; print $line}' "$filepath" |
      sed '/^$/d' > 23_hmmalign_trim1/"$file" &
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
    file=$(basename "$filepath")
    # printf "%s " "$file"
    filter_align.pl -in "$filepath" -out 23_hmmalign_trim2/"$file" -log 23_hmmalign_trim2_log/"$file" \
                    -depth $min_depth -pct_depth $min_pct_depth -min_pct_aligned $min_pct_aligned &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo
}

run_calc_trees() {
  case $scriptname in
    'pandagma pan')  dir_prefix=2 ;;
    'pandagma fam')  dir_prefix=2 ;;
    'pandagma fsup') dir_prefix=4
  esac
  echo; echo "== Calculate trees =="
  cd "${WORK_DIR}"

  mkdir -p "${dir_prefix}4_trees"

  # By default, FastTreeMP uses all available threads.
  # It is more efficient to run more jobs on one core each by setting an environment variable.
  OMP_NUM_THREADS=1
  export OMP_NUM_THREADS
  for filepath in "${dir_prefix}3_hmmalign_trim2"/*; do
    file=$(basename "$filepath")
    echo "  Calculating tree for $file"
    if [[ "$scriptname" =~ "pandagma pan" ]]; then  # pan; calculate from nucleotide alignments
      echo "fasttree -nt -quiet $filepath > ${dir_prefix}4_trees/$file"
      fasttree -quiet -nt "$filepath" > "${dir_prefix}4_trees/$file"
    else # fam or fsup; calculate from protein alignments
      echo "fasttree -quiet $filepath > ${dir_prefix}4_trees/$file"
      fasttree -quiet "$filepath" > "${dir_prefix}4_trees/$file" &
    fi
    # allow to execute up to $NPROC concurrent asynchronous processes
    if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
  done
  wait
}

run_xfr_aligns_trees() {
  echo; echo "Copy alignment and tree results into output directory"

  full_out_dir="${out_dir}"

  cd "${submit_dir}" || exit

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p "$full_out_dir"
  fi

  for dir in 20_align* 21_hmm 22_hmmalign 23_hmmalign_trim2 24_trees; do
    if [ -d "${WORK_DIR}"/$dir ]; then
      echo "Copying directory $dir to output directory"
      cp -r "${WORK_DIR}"/$dir "${full_out_dir}"/
    else
      echo "Warning: couldn't find dir ${WORK_DIR}/$dir; skipping"
    fi
  done
}

main_pan_fam() {
  if [ "$#" -eq 0 ]; then
    echo >&2 "$HELP_DOC" && exit 0;
  fi
  
  optarg_order_method="null"
  
  NPROC=$( ( command -v nproc > /dev/null && nproc ) || getconf _NPROCESSORS_ONLN)
  TIMEFORMAT=$'\nreal\t%3lR\nuser\t%3lU\nsys\t%3lS\ncpu\t%P%%' # add % cpu
  CONFIG="null"
  step="all"
  retain="no"
  
  export NPROC=${NPROC:-1}
  export MMSEQS_NUM_THREADS=${NPROC} # mmseqs otherwise uses all cores by default
  
  # mmseqs uses significant number of threads on its own. Set a maximum, which may be below NPROC.
  : ${MMSEQSTHREADS:=$(( 4 < NPROC ? 4 : NPROC ))}
  
  ##########
  # Command-line interpreter
  
  while getopts "c:d:w:s:n:O:o:f:rvhm" opt
  do
    case $opt in
      c) CONFIG=$OPTARG; echo "Config: $CONFIG" ;;
      d) optarg_data_dir=$(realpath -e "$OPTARG"); echo "Data dir: $optarg_data_dir" ;;
      w) optarg_work_dir=$OPTARG; echo "Work dir: $optarg_work_dir" ;;
      s) step=$OPTARG; echo "step(s): $step" ;;
      n) NPROC=$OPTARG; echo "processors: $NPROC" ;;
      O) optarg_order_method=$OPTARG; echo "order method: $optarg_order_method" ;;
      o) out_dir=$OPTARG; echo "out_dir: $out_dir" ;;
      f) fam_dir=$(realpath --canonicalize-existing "$OPTARG") ;;
      r) retain="yes" ;;
      v) version ;;
      h) echo >&2 "$HELP_DOC" && exit 0 ;;
      m) printf >&2  "%s\n%s\n" "$HELP_DOC" "$MORE_INFO" && exit 0 ;;
      *) echo >&2 echo "$HELP_DOC" && exit 1 ;;
    esac
  done
  
  export fam_dir

  case ${step} in 
    all|summarize|xfr_aligns_trees) : "${out_dir:=out_pandagma}" ;; # default to ./out_pandagma 
    *) if [ "${out_dir:-}" ]; then
         printf '\noption -o OUT_DIR not applicable to step %s\n' "$step" >&2
         exit 1
      fi
  esac
  
  if [ "$CONFIG" == "null" ]; then
    printf "\nPlease provide the path to a config file: -c CONFIG\n" >&2
    printf "\nRun \"%s -h\" for help.\n\n" "$scriptname" >&2
    exit 1;
  else
    export CONF=${CONFIG}
  fi
  
  # Add shell variables from config file
  # shellcheck source=/dev/null
  . "${CONF}"
  
  export ANN_REX=${annot_str_regex}

  shift $(( OPTIND - 1 ))
  
  if [ "${scriptname}" != 'pandagma pan' ] && \
     [ "${scriptname}" != 'pandagma TEfilter' ] && \
     [ "${order_method:-}${optarg_order_method:-}" != "null" ]; then
    printf '\norder_method is applicable only to pandagma pan\n' >&2
    exit 1 
  elif [ "${order_method:-}" != "$optarg_order_method" ] && [ "$optarg_order_method" != "null" ]; then
    echo "Command-line option for order_method was \"$optarg_order_method\", overriding"
    echo "the setting of \"$order_method\" from the config."
    order_method=$optarg_order_method
  fi
  
  ##########
  # Set the DATA_DIR, defaulting to ./data if -d DATA_DIR option not specified
  export DATA_DIR=${optarg_data_dir:-${PWD}/data}
  echo "From common: $DATA_DIR"
  # Set the work_dir, defaulting to ./work_pandagma if -w WORK_DIR option not specified
  export WORK_DIR=${optarg_work_dir:-${PWD}/work_pandagma}
  
  if [[ $step == "clean" ]] ; then
    run_"$step"
    echo "Command \"$step\" was run for cleanup in the work directory: $WORK_DIR"
    exit 0;
  fi
  
  # Get paths to the fasta and annotation files
  canonicalize_paths
  
  # Check for existence of third-party executables
  missing_req=0
  for program in $dependencies; do
    if ! type "$program" &> /dev/null; then
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
    printf "\n  PATH=%s/bin:\%s\n\n" "$PWD" "$PATH"
    exit 1; 
  fi
  
  # create work directory if it doesn't already exist
  mkdir -p "${WORK_DIR}"

  # Run all specified steps (except clean -- see below).
  if [[ $step =~ "all" ]]; then
    for command in $commandlist; do
      echo "RUNNING STEP run_$command"
      time run_"$command"
      echo
    done
  else
    run_"${step}";
  fi
  
  if ! command -v run_clean; then
    exit # run_clean not defined for this workflow
  elif [[ $step =~ "all" ]] && [[ $retain == "yes" ]]; then
    echo "Flag -r (retain) was set, so skipping clean-up of the work directory."
    echo "If you wish to do a cleanup separately, you can call "
    echo "  ${scriptname} -c $CONF -s clean";
  elif [[ $step =~ "all" ]] && [[ $retain == "no" ]]; then
    echo "Calling the medium cleanup function \"run_clean\" ..."
    run_clean  
  fi
}

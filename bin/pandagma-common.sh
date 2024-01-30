#!/usr/bin/env bash

version="2023-01-28"
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

run_calc_trees() {
  case $scriptname in
    'pandagma pan') fasttree_opts='-nt'; dir_prefix=2 ;;
    'pandagma fam') fasttree_opts=''; dir_prefix=2 ;;
    'pandagma fsup') fasttree_opts=''; dir_prefix=4
  esac
  echo; echo "== Calculate trees =="
  cd "${WORK_DIR}"

  mkdir -p "${dir_prefix}4_trees"

  echo; echo "== Move small (<4) and very low-entropy families (sequences are all identical) to the side =="
  mkdir -p "${dir_prefix}3_pan_aug_small_or_identical"
  min_seq_count=4
  # Below, "count" is the number of unique sequences in the alignment.
  for filepath in "${dir_prefix}3_hmmalign_trim2"/*; do
    file=$(basename "$filepath")
    count=$(awk '$1!~/>/ {print FILENAME "\t" $1}' "$filepath" | sort -u | wc -l);
    if [[ $count -lt $min_seq_count ]]; then
      echo "Set aside small or low-entropy family $file";
      mv "$filepath" "${dir_prefix}3_pan_aug_small_or_identical"/
    fi;
  done

  # By default, FastTreeMP uses all available threads.
  # It is more efficient to run more jobs on one core each by setting an environment variable.
  OMP_NUM_THREADS=1
  export OMP_NUM_THREADS
  for filepath in "${dir_prefix}3_hmmalign_trim2"/*; do
    file=$(basename "$filepath")
    echo "  Calculating tree for $file"
    fasttree "${fasttree_opts}" -quiet "$filepath" > "${dir_prefix}4_trees/$file" &
    # allow to execute up to $NPROC concurrent asynchronous processes
    if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
  done
  wait
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
      d) optarg_data_dir=$OPTARG; echo "Data dir: $optarg_data_dir" ;;
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
    all|summarize|xfr_aligns_trees) : "${out_dir:=out_pandagma}" ;; # default to ./pandagma_out 
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

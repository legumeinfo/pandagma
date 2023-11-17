#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates gene-family 
# orthogroups using the program mmseqs2.
# This workflow, pandagma-psup.sh, is to be used to add selected annotation sets to
# a previously calculated pangene set
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2023

scriptname=`basename "$0"`
version="2023-11-16"
set -o errexit -o errtrace -o nounset -o pipefail

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

HELP_DOC="""Place annotation sets (CDS and protein) into pangene or gene family sets, using mmseqs2
to compare the indicated supplemental annotations against the output from a previous pangene run.

The following files need to be available in the work directory, from a previous pangene calculation:
  18_syn_pan_aug_extra.clust.tsv and  21_pan_fasta_clust_rep_cds.fna   

Usage:
  ./$scriptname -c CONFIG_FILE [options]

  Required:
           -c (path to the config file)

  Options: -s (subcommand to run. If \"all\" or omitted, all steps will be run; otherwise, run specified step)
           -w (working directory, for temporary and intermediate files. 
                Must be specified in config file if not specified here.)
           -n (number of processors to use. Defaults to number of processors available to the parent process)
           -v (version)
           -h (help)
           -m (more information)

Environment requirements: The following packages need to be available in your PATH:
    mmseqs2 

Also, please add the pandagma utility programs in the bin directory adjacent to the pandagma- scripts, e.g.
    PATH=$PWD/bin:\$PATH

Subcommands (in order they are usually run):
                all - All of the steps below (Or equivalently: omit the -s flag; \"all\" is default).
             ingest - Prepare the assembly and annotation files for analysis.
    search_pangenes - Search supplementary sequence sets against representative pangene sequences.
          add_extra - Place sequences into pangenes based on top mmseqs hits.
   filter_to_pctile - 
          summarize - Move results into output directory, and report summary statistics.
"""

MORE_INFO="""
Variables in pandagma config file (Set the config with the CONF environment variable)
        consen_iden - Minimum identity threshold for consensus generation [0.30]
          clust_cov - Minimum alignment coverage [0.40]
      consen_prefix - Prefix to use in orthogroup names
       out_dir_base - Base name for the output directory [default: './out']
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g. 
                        \"([^.]+\.[^.]+)\..+\"
                          for two dot-separated fields, e.g. vigan.Shumari
                        or \"(\D+\d+\D+)\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
           work_dir - Working directory, for temporary and intermediate files. 

File sets (arrays):
          cds_files
   annotation_files
          cds_files_extra  (optional)
   annotation_files_extra  (optional)
      protein_files
   annotation_files_sup
      protein_files_sup
          cds_files_sup
 """

########################################
# Helper functions begin here

version() {
  echo $scriptname $version
}

##########
canonicalize_paths() {
  echo "Entering canonicalize_paths."

  annotation_files=($(realpath --canonicalize-existing "${annotation_files[@]}"))
  cds_files=($(realpath --canonicalize-existing "${cds_files[@]}"))
  protein_files=($(realpath --canonicalize-existing "${protein_files[@]}"))

  annotation_files_sup=($(realpath --canonicalize-existing "${annotation_files_sup[@]}"))
  cds_files_sup=($(realpath --canonicalize-existing "${cds_files_sup[@]}"))
  protein_files_sup=($(realpath --canonicalize-existing "${protein_files_sup[@]}"))

  if (( ${#cds_files_extra[@]} > 0 ))
  then
    cds_files_extra=($(realpath --canonicalize-existing "${cds_files_extra[@]}"))
    annotation_files_extra=($(realpath --canonicalize-existing "${annotation_files_extra[@]}"))
  fi

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
# Add positional information from GFF3 or 4- or 6-column BED to FASTA IDs
# BED start coordinate converted to 1-based
  cd "${WORK_DIR}"
  echo; echo "Run ingest: from fasta and gff or bed data, create fasta with IDs containing positional info."

  export ANN_REX=${annot_str_regex}
  mkdir -p 02_fasta_nuc_sup 02_fasta_prot_sup 01_posn_hsh_sup stats_sup
  mkdir -p 02_fasta_nuc 02_fasta_prot 01_posn_hsh stats

  # Prepare the tmp.gene_count_start to be joined, in run_summarize, with tmp.gene_count_end_pctl??_end.
  # This is captured from the gene IDs using the annot_str_regex set in the config file.
  cat /dev/null > stats_sup/tmp.gene_count_start_sup
  cat /dev/null > stats_sup/tmp.fasta_seqstats
  start_time=`date`
  printf "Run started at: $start_time\n" > stats_sup/tmp.timing

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
    printf "  Main:  " >> stats_sup/tmp.fasta_seqstats
    cat_or_zcat "${cds_files[file_num]}" | calc_seq_stats >> stats_sup/tmp.fasta_seqstats
  done

  echo "  Get position information from the extra annotation sets, if any."
  if (( ${#cds_files_extra[@]} > 0 ))
  then
    cat /dev/null > 02_all_extra_cds.fna # Collect all starting sequences, for later comparisons
    for (( file_num = 0; file_num < ${#cds_files_extra[@]} ; file_num++ )); do
      file_base=$(basename ${cds_files_extra[file_num]%.*})
      cat_or_zcat "${cds_files_extra[file_num]}" >> 02_all_extra_cds.fna  # Collect original seqs for later comparisons
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra[file_num]}" |
        gff_or_bed_to_hash5.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files_extra[file_num]}" \
                            -hash 01_posn_hsh/$file_base.hsh \
                            -out 02_fasta_nuc/$file_base
      # calc basic sequence stats
      annot_name=$(basename 02_fasta_nuc/$file_base | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
      printf "  Extra: " >> stats_sup/tmp.fasta_seqstats
      cat_or_zcat "${cds_files_extra[file_num]}" | calc_seq_stats >> stats_sup/tmp.fasta_seqstats
    done
  fi

  echo "  Get position information from the supplementary annotation sets, if any."
  if (( ${#cds_files_sup[@]} > 0 ))
  then
    cat /dev/null > 02_all_sup_cds.fna # Collect all starting sequences, for later comparisons
    for (( file_num = 0; file_num < ${#cds_files_sup[@]} ; file_num++ )); do
      file_base=$(basename ${cds_files_sup[file_num]%.*})
      cat_or_zcat "${cds_files_sup[file_num]}" >> 02_all_sup_cds.fna  # Collect original seqs for later comparisons
      echo "  Adding positional information to sup fasta file $file_base"
      cat_or_zcat "${annotation_files_sup[file_num]}" |
        gff_or_bed_to_hash5.awk > 01_posn_hsh_sup/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files_sup[file_num]}" \
                            -hash 01_posn_hsh_sup/$file_base.hsh \
                            -out 02_fasta_nuc_sup/$file_base
      # calc basic sequence stats
      annot_name=$(basename 02_fasta_nuc_sup/$file_base | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
      printf "  sup: " >> stats_sup/tmp.fasta_seqstats
      cat_or_zcat "${cds_files_sup[file_num]}" | calc_seq_stats >> stats_sup/tmp.fasta_seqstats
    done
  fi

  echo "  Count starting sequences, for later comparisons"
  for file in 02_fasta_nuc/*.$fna 02_fasta_nuc_sup/*.$fna; do
    awk '$0~/UNDEFINED/ {ct++} 
      END{if (ct>0){print "Warning: " FILENAME " has " ct " genes without position (HASH UNDEFINED)" } }' $file
    cat $file | grep '>' | perl -pe 's/__/\t/g' | cut -f2 | # extracts the GeneName from id-string-with-coords
      perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex.+/$1/' |
      grep -v UNDEFINED | sort | uniq -c | awk '{print $2 "\t" $1}' >> stats_sup/tmp.gene_count_start
  done

  echo "  Also get protein files"
  if (( ${#protein_files[@]} > 0 ))
  then
    for (( file_num = 0; file_num < ${#protein_files[@]} ; file_num++ )); do
      echo "  Copying protein file ${protein_files[file_num]}"
      cp "${protein_files[file_num]}" 02_fasta_prot/
    done
  fi
  if (( ${#protein_files_sup[@]} > 0 ))
  then
    for (( file_num = 0; file_num < ${#protein_files_sup[@]} ; file_num++ )); do
      echo "  Copying protein file ${protein_files_sup[file_num]}"
      cp "${protein_files_sup[file_num]}" 02_fasta_prot_sup/
    done
  fi

  sort -o stats_sup/tmp.gene_count_start stats_sup/tmp.gene_count_start
}

##########
run_search_pangenes() {
  cd "${WORK_DIR}"
  echo "Search provided supplemental annotation sets against a file of pangene representative sequences"

  if [ -d 33_mmseqs_pan_match ]; then rm -rf 33_mmseqs_pan_match; fi
  mkdir -p 33_mmseqs_pan_match 33_mmseqs_tmp

  SEQTYPE=3 # 3=nuc; 1=pep

  for filepath in 02_fasta_nuc_sup/*; do 
    base=`basename $filepath .$fna`
    echo "  Search $base.$fna against representative sequences to determine best pangene match"
    mmseqs easy-search $filepath \
                       21_pan_fasta_clust_rep_cds.fna \
                       33_mmseqs_pan_match/$base.x.21_pan_fasta_clust_rep_cds.m8 \
                       33_mmseqs_tmp \
                       --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null
  done

  for filepath in 33_mmseqs_pan_match/*; do
    base=`basename $filepath`
    cat $filepath | perl -pe 's/^(\S+)\t([^_]+)__\S+/$1\t$2/' | 
      sort -k1,1 -k12nr,12nr | top_line.awk > 33_mmseqs_pan_match/$base.top &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_add_extra() {
  cd "${WORK_DIR}"
  echo "== Place sequences into pangenes based on top mmseqs hits, using the consen_iden threshold of $consen_iden."

  cat 33_mmseqs_pan_match/*.m8.top | top_line.awk |
    perl -pe 's/^\S+__(\S+)__\d+__\d+__[+-]\t/$1\t/' |
    awk -v IDEN=${consen_iden} '$3>=IDEN {print $2 "\t" $1}' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > 41_sup_in_21_rep_cds.clust.tsv
  
  echo "  Retrieve sequences for the extra genes"
    if [ -d 42_sup_in_21_rep ]; then rm -rf 42_sup_in_21_rep; fi
    mkdir -p 42_sup_in_21_rep
    get_fasta_from_family_file.pl "${cds_files_sup[@]}" \
       -fam 41_sup_in_21_rep_cds.clust.tsv -out 42_sup_in_21_rep/

    echo "  Make augmented cluster sets - starting from 19_pan_aug_leftover_merged_cds/ if available;"
    echo "  else from 18_syn_pan_aug_extra.clust.tsv"
    if [[ -d 19_pan_aug_leftover_merged_cds ]]; then
      : # do nothing; the directory and file(s) exist
    else 
      if [ -f 18_syn_pan_aug_extra.clust.tsv ]; then
        mkdir -p 19_pan_aug_leftover_merged_cds
        get_fasta_from_family_file.pl "${cds_files[@]}" "${cds_files_extra[@]}" \
          -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_cds
      else
        echo "Neither 19_pan_aug_leftover_merged_cds/ nor 18_syn_pan_aug_extra.clust.tsv is available"
        echo "in the work directory. Please provide at least 18_syn_pan_aug_extra.clust.tsv in ${WORK_DIR}."
        exit 1
      fi
    fi
    augment_cluster_sets.awk leftovers_dir=42_sup_in_21_rep 19_pan_aug_leftover_merged_cds/* |
      cat > 48_syn_pan_sup.clust.tsv

    echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
    perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 48_syn_pan_sup.clust.tsv \
      > 48_syn_pan_sup.hsh.tsv

    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    echo "    Fasta file:" "${cds_files[@]}"
    if [ -d 49_pan_sup_merged_cds ]; then rm -rf 49_pan_sup_merged_cds; fi
    mkdir -p 49_pan_sup_merged_cds
    get_fasta_from_family_file.pl "${cds_files[@]}"  "${cds_files_extra[@]}"  "${cds_files_sup[@]}" \
      -fam 48_syn_pan_sup.clust.tsv -out 49_pan_sup_merged_cds

    echo "  Merge files in 49_pan_sup_merged_cds, prefixing IDs with panID__"
    merge_files_to_pan_fasta.awk 49_pan_sup_merged_cds/* > 49_pan_sup_merged_cds.fna
}

##########
run_pick_exemplars() {
  echo; echo "== Pick representative (exemplar) sequence for each pan-gene set (protein and CDS) =="
  cd "${WORK_DIR}"

  echo "  Get all protein sequences corresponding with 48_syn_pan_extra.counts.tsv"
  cat /dev/null > 50_pan_fasta_prot.faa
  for filepath in 02_fasta_prot/*.gz 02_fasta_prot_sup/*.gz; do
    zcat $filepath >> 50_pan_fasta_prot.faa
  done

  echo "  Get protein sequences into pan-gene sets, corresponding with 19_pan_aug_leftover_merged_cds.fna"
  if [ -d 49_pan_sup_merged_cds ]; then rm -rf 49_pan_sup_merged_cds; fi
  mkdir -p 49_pan_sup_merged_cds
  get_fasta_from_family_file.pl 50_pan_fasta_prot.faa \
              -fam 48_syn_pan_sup.clust.tsv -out 49_pan_sup_merged_cds

  echo "  Get all protein sequences from files in 49_pan_sup_merged_cds"
  merge_files_to_pan_fasta.awk 49_pan_sup_merged_cds/* > 49_pan_sup_merged_cds.faa

  echo "  Pick a representative sequence for each pangene set - as a sequence with the median length for that set."
  echo "    == first proteins:"
  cat 49_pan_sup_merged_cds.faa | pick_family_rep.pl \
    -nostop -prefer $preferred_annot -out 51_pan_fasta_clust_rep_prot.faa

  echo "    == then CDS sequences, corresponding with 51_pan_fasta_clust_rep_prot.faa"
  cat 51_pan_fasta_clust_rep_prot.faa | awk '$1~/^>/ {print substr($1,2) "__" $2}' > lists/lis.51_pan_fasta_clust_rep
  get_fasta_subset.pl -in 49_pan_sup_merged_cds.fna -list lists/lis.51_pan_fasta_clust_rep \
                    -clobber -out 51_pan_fasta_clust_rep_cds.fna

  perl -pi -e 's/__/  /' 51_pan_fasta_clust_rep_cds.fna
  perl -pi -e 's/__/  /' 51_pan_fasta_clust_rep_prot.faa

  echo "  Retrieve genes present in the original CDS files but absent from 48_syn_pan_sup.hsh.tsv"
  cut -f2 48_syn_pan_sup.hsh.tsv | LC_ALL=C sort > lists/lis.48_syn_pan_sup
  cat 02_all_main_cds.fna 02_all_extra_cds.fna 02_all_main_cds_sup.fna > 02_all_cds.fna
  get_fasta_subset.pl -in 02_all_cds.fna -out 48_syn_pan_sup_complement.fna \
    -lis lists/lis.48_syn_pan_sup -xclude -clobber
}

##########
run_filter_to_pctile() {
  cd "${WORK_DIR}"
  echo "  Calculate matrix of gene counts per orthogroup and annotation set"
  calc_pan_stats.pl -annot_regex $ANN_REX -pan 48_syn_pan_sup.clust.tsv -out 48_syn_pan_extra.counts.tsv
  max_annot_ct=$(cat 48_syn_pan_extra.counts.tsv |
                       awk '$1!~/^#/ {print $2}' | sort -n | uniq | tail -1)

  echo "  Select orthogroups with genes from selected percentiles of max_annot_ct annotation sets"
  for percentile in $pctl_low $pctl_med $pctl_hi; do
    cat 48_syn_pan_extra.counts.tsv |
      awk -v PCTL=$percentile -v ANNCT=$max_annot_ct '$2>=ANNCT*(PCTL/100)&& $1!~/^#/ {print $1}' |
        cat > lists/lis.48_syn_pan_extra.pctl${percentile}

    echo "  Get a fasta subset with only genes from at least $(($max_annot_ct*$pctl_low/100)) annotation sets"
    get_fasta_subset.pl -in 51_pan_fasta_clust_rep_cds.$fna \
                        -list lists/lis.48_syn_pan_extra.pctl${percentile} \
                        -clobber -out 52_pan_fasta_rep_pctl${percentile}_cds.$fna

    echo "  Get a clust.tsv file with orthogroups with at least min_core_prop*max_annot_ct annotation sets"
    join <(LC_ALL=C sort -k1,1 lists/lis.48_syn_pan_extra.pctl${percentile}) \
         <(LC_ALL=C sort -k1,1 48_syn_pan_sup.clust.tsv) |
            cat > 52_pan_sup_pctl${percentile}.clust.tsv
  done
}


##########
run_summarize() {
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"

  WORK_DIR=$work_dir
  cd "${WORK_DIR}"
  echo "  work_dir: $PWD"

  echo "  Calculate matrix of gene counts per orthogroup and annotation set"
  calc_pan_stats.pl -annot_regex $ANN_REX 
     -pan 21_pan_fasta_clust_rep_cds.clust.tsv -out 21_pan_fasta_clust_rep_cds.counts.tsv
  max_annot_ct=$(cat 21_pan_fasta_clust_rep_cds.counts.tsv | 
                       awk '$1!~/^#/ {print $2}' | sort -n | uniq | tail -1)
 
  conf_base=`basename $CONF .conf`
  full_out_dir="${out_dir_base}_$conf_base"
  stats_file=${full_out_dir}/stats_sup.$conf_base.txt
  export ANN_REX=${annot_str_regex}

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  for file in 48_syn_pan_sup.clust.tsv 48_syn_pan_sup.hsh.tsv \
      52_pan_fasta_rep_pctl${pctl_low}_cds.fna \
      52_pan_fasta_rep_pctl${pctl_med}_cds.fna \
      52_pan_fasta_rep_pctl${pctl_hi}_cds.fna ; do
    if [ -f ${WORK_DIR}/$file ]; then
      cp ${WORK_DIR}/$file ${full_out_dir}/
    else 
      echo "Warning: couldn't find file ${WORK_DIR}/$file; skipping"
    fi
  done

#I AM HERE
  printf "Run of program $scriptname, version $version\n" > ${stats_file}

  end_time=`date`
  cat ${WORK_DIR}/stats_sup/tmp.timing >> ${stats_file}
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
  if [ -f ${WORK_DIR}/stats_sup/tmp.fasta_seqstats ]; then
    cat ${WORK_DIR}/stats_sup/tmp.fasta_seqstats >> ${stats_file}

  printf "\n  Avg:   " >> ${stats_file} 
    cat ${WORK_DIR}/stats_sup/tmp.fasta_seqstats | transpose.pl |
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
    sort | uniq -c | awk '{print $2 "\t" $1}' > ${WORK_DIR}/stats_sup/tmp.gene_count_all_end

  paste ${WORK_DIR}/stats_sup/tmp.gene_count_start \
        ${WORK_DIR}/stats_sup/tmp.gene_count_all_end |
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

########################################
# Main program

NPROC=$(command -v nproc > /dev/null && nproc || getconf _NPROCESSORS_ONLN)
CONFIG="null"
optarg_work_dir="null"
step="all"

export NPROC=${NPROC:-1}

pandagma_conf_params='consen_iden out_dir_base pctl_low pctl_med pctl_hi consen_prefix annot_str_regex work_dir'

##########
# Command-line interpreter

while getopts "c:w:s:n:o:vhm" opt
do
  case $opt in
    c) CONFIG=$OPTARG; echo "Config: $CONFIG" ;;
    w) optarg_work_dir=$OPTARG; echo "Work dir: $optarg_work_dir" ;;
    s) step=$OPTARG; echo "step(s): $step" ;;
    n) NPROC=$OPTARG; echo "processors: $NPROC" ;;
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

# Get paths to the fasta and annotation files
canonicalize_paths

# Check for existence of third-party executables
missing_req=0
for program in mmseqs hmmscan famsa; do
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

# Run all specified steps 
commandlist="ingest search_pangenes add_extra filter_to_pctile summarize"

if [[ $step =~ "all" ]]; then
  for command in $commandlist; do
    echo "RUNNING STEP run_$command"
    run_$command
    echo
  done
else
  run_${step};
fi


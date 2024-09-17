#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates gene-family 
# orthogroups using the programs mmseqs, the hmmer package and others.
# This workflow, pandagma-fsup.sh, is to be used to add selected annotation sets to
# a collection of HMMs calculated as part of a prior full run of pandagma-fam.sh.
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2023
#
scriptname='pandagma fsup'

# shellcheck source=/dev/null
. pandagma-common.sh

define HELP_DOC <<'EOS'
Place annotation sets (CDS and protein) into pangene or gene family sets, using hmmsearch 
to compare the indicated annotations against HMMs calculated in a prior run of pandagma-pan or pandagma-fam.

Usage:
  $scriptname -c CONFIG_FILE [options]

  Required:
           -c (path to the config file)
           -f FAM_DIR (path to directory from previous \"pandagma fam\ -o FAM_DIR\" run)

  Options: -s (subcommand to run. If \"all\" or omitted, all steps will be run; otherwise, run specified step)
           -w (working directory, for temporary and intermediate files [default: './work_pandagma'].)
           -o OUTPUT_DIR (name for the output directory [default: './out_pandagma'].
                Applicable only to "all" and "summarize" steps.)
           -n (number of processors to use. Defaults to number of processors available to the parent process)
           -v (version)
           -h (help)
           -m (more information)

Environment requirements: The following packages need to be available in your PATH:
    hmmsearch hmmalign fasttree

Also, please add the pandagma utility programs in the bin directory adjacent to pandagma-fam.sh, e.g.
    PATH=$PWD/bin:\$PATH

Subcommands (in order they are usually run):
                all - All of the steps below
                        (Or equivalently: omit the -s flag; \"all\" is default).
             ingest - Prepare the assembly and annotation files for analysis.
         fam_consen - Generate a consensus sequence for each family.
    search_families - Search provided annotation sets (protein) against family consensus sequences.
   realign_and_trim - Build HMMs and trim the alignments, preparatory to calculating trees.
         calc_trees - Calculate gene trees.
          summarize - Move results into output directory, and report summary statistics.
EOS

define MORE_INFO <<'EOS'
cat <<'EOS'
Variables in pandagma config file (Set the config with the CONF environment variable)
        consen_iden - Minimum identity threshold for consensus generation [0.30]
          clust_cov - Minimum alignment coverage [0.40]
      consen_method - Method for placing sequences into pangenes. One of:
                        hmmemit cons align_sample [align_sample]
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g. 
                        \"([^.]+\.[^.]+)\..+\"
                          for two dot-separated fields, e.g. vigan.Shumari
                        or \"(\D+\d+\D+)\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    min_align_count - Minimum number of sequences in a family to trigger alignments, modeling, and trees [4]
min_annots_in_align - Minimum number of distinct annotation groups in an alignment to retain it [2]

File sets (arrays):
   annotation_files
          cds_files
      protein_files
EOS

########################################
# Helper functions begin here

canonicalize_paths() {
  echo "Entering canonicalize_paths."

  if ! [[ -v fam_dir ]]; then
    echo "pandagma fsup -f FAM_DIR not specified" >&2
    exit 1
  fi

  cd "${DATA_DIR}" || exit

  export ANN_REX=${annot_str_regex}

  mapfile -t cds_files < <(realpath --canonicalize-existing "${cds_files[@]}")
  mapfile -t protein_files < <(realpath --canonicalize-existing "${protein_files[@]}")

  cd "${OLDPWD}" || exit
  readonly submit_dir=${PWD}

  prot_file=$(basename "${protein_files[0]}" .gz)
  faa="${prot_file##*.}"
}

########################################
# run functions 

##########
run_ingest() {
# Retrieve nucleotide and protein data sets
  cd "${WORK_DIR}" || exit
  echo; echo "Run ingest: pull nucleotide and protein data sets into the work directory."
  echo "Note that this is a simpler ingest than for the primary gene family or pangene construction,"
  echo "as the identifiers aren't modified with positional information."
  
  if [ -d 02_fasta_nuc_sup ]; then rm -rf 02_fasta_nuc_sup; fi
  if [ -d 02_fasta_prot_sup ]; then rm -rf 02_fasta_prot_sup; fi
  mkdir -p 02_fasta_nuc_sup 02_fasta_prot_sup stats

  # Prepare the tmp.gene_count_start to be joined, in run_summarize, with tmp.gene_count_end_sup
  # This is captured from the gene IDs using the annot_str_regex set in the config file.
  cat /dev/null > stats/tmp.gene_count_start_sup
  cat /dev/null > stats/tmp.fasta_seqstats_sup
  start_time=$(date)
  echo "Run started at: $start_time" > stats/tmp.timing

  echo "  Pull the protein files locally"
  cat /dev/null > 02_all_main_prot_sup.faa # Collect all starting protein sequences, for later comparisons
  for (( file_num = 0; file_num < ${#protein_files[@]} ; file_num++ )); do
    file_base=$(basename "${protein_files[file_num]%.*}")
    zcat "${protein_files[file_num]}" >> 02_all_main_prot_sup.faa # Collect original seqs for later comparisons
    zcat "${protein_files[file_num]}" > 02_fasta_prot_sup/"$file_base"
    # calc basic sequence stats
    annot_name=$(basename 02_fasta_prot/"$file_base" | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
    echo "  CHECK: report annot_name to stats file? [$annot_name]"
    printf "  Added with $consen_method:  " >> stats/tmp.fasta_seqstats_sup
    cat_or_zcat "${protein_files[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats_sup
  done

  echo "  Pull the CDS files locally"
  cat /dev/null > 02_all_main_cds_sup.fna # Collect all starting cds sequences, for later comparisons
  for (( file_num = 0; file_num < ${#cds_files[@]} ; file_num++ )); do
    file_base=$(basename "${cds_files[file_num]%.*}")
    zcat "${cds_files[file_num]}" >> 02_all_main_cds_sup.fna # Collect original seqs for later comparisons
    zcat "${cds_files[file_num]}" > 02_fasta_nuc_sup/"$file_base"
  done

  echo "  Count starting sequences, for later comparisons"
  for file in 02_fasta_prot_sup/*."$faa"; do
    grep '>' "$file" | perl -pe 's/>(\S+)/$1/; s/^(\S+)\s+.+/$1/' | 
      perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex.+/$1/' |
      sort | uniq -c | awk '{print $2 "\t" $1}' >> stats/tmp.gene_count_start_sup
  done

  sort -o stats/tmp.gene_count_start_sup stats/tmp.gene_count_start_sup
}

run_fam_consen() {
  echo
  echo; echo "== Generate a consensus sequence for each family =="
  cd "${WORK_DIR}" || exit

  if [ "$consen_method" == "hmmemit" ]; then
    if [ -d 21_hmmemit ]; then rm -rf 21_hmmemit; fi
    mkdir -p 21_hmmemit
  
    MINL=0.25 # show consensus as 'any' (X/N) unless >= this fraction
    MINU=0.5  # show consensus as upper case if >= this fraction
    for filepath in "${fam_dir}"/21_hmm/* ; do
      base=$(basename "$filepath")
      hmmemit -C --minl $MINL --minu $MINU -o 21_hmmemit/"$base" "$filepath" &
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
    done
    wait
  
    cat 21_hmmemit/* | perl -pe 's/-consensus//' > 30_consen.faa

  elif [ "$consen_method" == "cons" ]; then
    if [ -d 30_hmmalign_cons ]; then rm -rf 30_hmmalign_cons; fi
    mkdir -p 30_hmmalign_cons
  
    for filepath in "${fam_dir}"/23_hmmalign_trim2/* ; do
      base=$(basename "$filepath")
      cons -sequence "$filepath" -outseq 30_hmmalign_cons/"$base" -name "$base" &
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
    done
    wait
  
    cat 30_hmmalign_cons/* > 30_consen.faa

  elif [ "$consen_method" == "align_sample" ]; then

    enough_seqs=10  # Don't search against more than this number of sequences in a family
    echo "  Prepare to search against a random sampling of $enough_seqs sequences from each family"
    echo "  (or as many as available if there are fewer than $enough_seqs in the family)."

    cat /dev/null > 23_hmmalign_trim2.fam_ID.faa

    for filepath in "${fam_dir}"/23_hmmalign_trim2/*; do
      base=$(basename "$filepath")
      fasta_to_table.awk "$filepath" | shuf |
        awk -v FAM="$base" -v MAX=$enough_seqs 'ct<MAX {print ">" FAM "__" $1; print $2; ct++}' >> 23_hmmalign_trim2.fam_ID.faa
    done

  else 
    echo "Unrecognized consen_method: $consen_method."
    echo "Expected values are \"hmmemit\", \"cons\", or \"align_sample\"."
    exit 1
  fi
}

##########
run_search_families() {
  cd "${WORK_DIR}" || exit
  echo "Search provided annotation sets (protein) against family consensus sequences, using method $consen_method"

  if [ -d 33_mmseqs_fam_match ]; then rm -rf 33_mmseqs_fam_match; fi
  mkdir -p 33_mmseqs_fam_match 33_mmseqs_tmp

  SEQTYPE=1 # 3=nuc; 1=pep

  if [ "$consen_method" == "hmmemit" ] || 
     [ "$consen_method" == "cons" ]; then
    for filepath in 02_fasta_prot_sup/*; do 
      base=$(basename "$filepath" ."$faa")
      echo "  Search $base.$faa against family hmmemit consensus sequences"
      mmseqs easy-search "$filepath" \
                         30_consen.faa \
                         33_mmseqs_fam_match/"$base".x.consen.m8 \
                         33_mmseqs_tmp \
                         --search-type ${SEQTYPE} --cov-mode 5 -c "${clust_cov}" 1>/dev/null
    done
    wait

    for filepath in 33_mmseqs_fam_match/*; do
      base=$(basename "$filepath")
      top_line.awk "$filepath" > 33_mmseqs_fam_match/"$base".top &
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
    done
    wait

  elif [ "$consen_method" == "align_sample" ]; then
    for filepath in 02_fasta_prot_sup/*; do 
      base=$(basename "$filepath" ."$faa")
      echo "  Search $base.$faa against a sampling of sequences in family alignments in 23_hmmalign_trim2 to determine best family match"
      mmseqs easy-search "$filepath" \
                         23_hmmalign_trim2.fam_ID.faa \
                         33_mmseqs_fam_match/"$base".x.23_hmmalign_trim2.fam_ID.m8 \
                         33_mmseqs_tmp \
                         --search-type ${SEQTYPE} --cov-mode 5 -c "${clust_cov}" 1>/dev/null
    done

    for filepath in 33_mmseqs_fam_match/*; do
      base=$(basename "$filepath")
      < "$filepath" perl -pe 's/^(\S+)\t([^_]+)__\S+/$1\t$2/' | 
        sort -k1,1 -k12nr,12nr | top_line.awk > 33_mmseqs_fam_match/"$base".top &
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
    done
    wait

  else
    echo "Unrecognized consen_method: $consen_method."
    echo "Expected values are \"hmmemit\", \"cons\", or \"align_sample\"."
    exit 1
  fi

  echo "== Place sequences into families based on top mmsearch hits, using the consen_iden threshold of $consen_iden."

  cat 33_mmseqs_fam_match/*.m8.top | top_line.awk |
    awk -v IDEN="${consen_iden}" '$3>=IDEN {print $2 "\t" $1}' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > 34_sup_vs_fam_consen.clust.tsv
}

##########
run_tabularize() {
  echo; echo "== Derive a table-format version of 34_sup_vs_fam_consen.clust.tsv"
  cd "${WORK_DIR}" || exit

  # Get table header
  pangene_tabularize.pl -pan 34_sup_vs_fam_consen.clust.tsv -annot_str_regex "$ANN_REX" > tmp.34_sup_vs_fam_consen.clust.tsv
  head -1 tmp.34_sup_vs_fam_consen.clust.tsv > tmp.table_header

  echo "  Sort, putting header row at top, and don't print pangenes that are all NONE"
    sort -k1,1 tmp.34_sup_vs_fam_consen.clust.tsv |
    sed '/^$/d; /^#pangene/d' |
    perl -lane '$ct=0; for $gn (@F){if ($gn=~/NONE/){$ct++}}; if ($ct<(scalar(@F)-1)){print $_}' |
    cat tmp.table_header - > 34_sup_vs_fam_consen.table.tsv

    rm tmp.34_sup_vs_fam_consen.clust.tsv
    rm tmp.table_header
}

##########
run_realign_and_trim() {
  echo; echo "== Retrieve sequences for each family, preparatory to aligning them =="
  cd "${WORK_DIR}" || exit
  mkdir -p 35_sup_in_fams_prot
  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  get_fasta_from_family_file.pl "${protein_files[@]}" \
    -fam 34_sup_vs_fam_consen.clust.tsv -out 35_sup_in_fams_prot

  echo; echo "== Realign to HMMs =="
  mkdir -p 42_hmmalign
  for filepath in 35_sup_in_fams_prot/*; do 
    file=$(basename "$filepath");
    hmmalign --trim --outformat A2M --amino -o 42_hmmalign/"$file" "${fam_dir}/21_hmm/$file" 35_sup_in_fams_prot/"$file" &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Trim HMM alignments to match-states =="
  mkdir -p 43_hmmalign_trim1
  for filepath in 42_hmmalign/*; do 
    file=$(basename "$filepath");
    < "$filepath" perl -ne 'if ($_ =~ />/) {print $_} else {$line = $_; $line =~ s/[a-z]//g; print $line}' |
      sed '/^$/d' > 43_hmmalign_trim1/"$file" &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Filter alignments prior to tree calculation =="
  mkdir -p 43_hmmalign_trim2 43_hmmalign_trim2_log
  min_depth=3
  min_pct_depth=20
  min_pct_aligned=20
  for filepath in 43_hmmalign_trim1/*; do 
    file=$(basename "$filepath")
    # printf "%s " "$file"
    filter_align.pl -in "$filepath" -out 43_hmmalign_trim2/"$file" -log 43_hmmalign_trim2_log/"$file" \
                    -depth $min_depth -pct_depth $min_pct_depth -min_pct_aligned $min_pct_aligned &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_summarize() {
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"

  cd "${WORK_DIR}" || exit
  echo "  work_dir: $PWD"

  echo "  Calculate matrix of gene counts per orthogroup and annotation set"
  calc_pan_stats.pl -annot_regex "$ANN_REX" -pan 34_sup_vs_fam_consen.clust.tsv -out 34_sup_vs_fam_consen.counts.tsv

  full_out_dir="${out_dir}"
  stats_file=${full_out_dir}/stats.txt

  cd "${submit_dir}" || exit

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p "$full_out_dir"
  fi

  for file in 34_sup_vs_fam_consen.clust.tsv 34_sup_vs_fam_consen.counts.tsv 34_sup_vs_fam_consen.table.tsv; do
    if [ -f "${WORK_DIR}"/$file ]; then
      cp "${WORK_DIR}"/$file "${full_out_dir}"/
    else 
      echo "Warning: couldn't find file ${WORK_DIR}/$file; skipping"
    fi
  done

  # Directories 35_sup_in_fams_prot 43_hmmalign_trim2 44_trees are transferred in xfr_aligns_trees (in pandagma-common.sh)

  end_time=$(date)

  {
    echo "Run of program $scriptname, version $version" > "${stats_file}"
    cat "${WORK_DIR}"/stats/tmp.timing 
    echo "Run ended at:   $end_time"
    echo 
  } >> "${stats_file}"

  echo "  Report parameters from config file"
  printf "Parameter  \tvalue\n" >> "${stats_file}"
  for key in ${pandagma_conf_params}; do
    printf '%-15s\t%s\n' "${key}" "${!key}" >> "${stats_file}"
  done

  {
    printf "Output directory for this run:\t%s" "$full_out_dir"

    echo "  Report orthogroup composition statistics for the three main cluster-calculation steps"

    echo "  Print sequence composition statistics for each annotation set"
    printf "\n== Sequence stats for protein files\n" 
    printf "  Class:  seqs     min max    N50    ave     annotation_name\n" 
    if [ -f "${WORK_DIR}"/stats/tmp.fasta_seqstats_sup ]; then
      cat "${WORK_DIR}"/stats/tmp.fasta_seqstats_sup 

    printf "\n  Avg:   " >> "${stats_file}" 
      < "${WORK_DIR}"/stats/tmp.fasta_seqstats_sup transpose.pl |
        perl -ane 'BEGIN{use List::Util qw(sum)}; 
                 if ($F[0]=~/^\d+/){
                   $sum=sum @F; $ct=scalar(@F); $avg=$sum/$ct;
                   printf " %4d ", $avg;
                 END{print "   all_annot_sets\n"}
                 }' 
  fi
  } >> "${stats_file}"

  echo "  Print per-annotation-set coverage stats (sequence counts, sequences retained)"
  #   tmp.gene_count_start was generated during run_ingest
  printf "\n== Proportion of initial genes retained in the \"aug_extra\" set:\n" \
    >> "${stats_file}"

  transpose.pl "${WORK_DIR}"/34_sup_vs_fam_consen.counts.tsv |
    perl -MList::Util=sum -lane 'if ($.>3){print $F[0], "\t", sum(@F[1..(scalar(@F)-1)]) }' |
      cat > "${WORK_DIR}"/stats/tmp.gene_count_end_sup

  paste "${WORK_DIR}"/stats/tmp.gene_count_start_sup \
        "${WORK_DIR}"/stats/tmp.gene_count_end_sup |
    awk 'BEGIN{print "  Start\tEnd\tPct_kept\tAnnotation_name"} 
        { printf "  %i\t%i\t%2.1f\t%s\n", $2, $4, 100*($4/$2), $1 }'  >> "${stats_file}"

  echo "  Print counts per accession"
  if [ -f "${WORK_DIR}"/_sup_vs_fam_consen.counts.tsv ]; then
    printf "\n== For all annotation sets, counts of genes-in-orthogroups and counts of orthogroups-with-genes:\n"
    printf "  gns-in-OGs  OGs-w-gns  OGs-w-gns/gns  pct-non-null-OGs  pct-null-OGs  annot-set\n" 
    transpose.pl "${WORK_DIR}"/34_sup_vs_fam_consen.counts.tsv | 
      perl -lane 'next if ($.<=3); 
        $ct=0; $sum=0; $nulls=0; $OGs=0;
        for $i (@F[1..(@F-1)]){ $OGs++; if ($i>0){$ct++; $sum+=$i} if ($i==0){$nulls++} }; 
        printf("  %d\t%d\t%.2f\t%.2f\t%.2f\t%s\n", $sum, $ct, 100*$ct/$sum, 100*($OGs-$nulls)/$OGs, 100*$nulls/$OGs, $F[0])' \
      >> "${stats_file}"
  fi

  if [ -f "${full_out_dir}"/34_sup_vs_fam_consen.clust.tsv ]; then
    printf "\nCounts of clusters by cluster size, file 34_sup_vs_fam_consen.clust.tsv:\n" >> "${stats_file}"
    awk '{print NF-1}' "${WORK_DIR}"/34_sup_vs_fam_consen.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> "${stats_file}"
  fi
}

########################################
# Main program

pandagma_conf_params='consen_iden clust_cov consen_method annot_str_regex'

# Run all specified steps.
# Those steps. Steps realign_and_trim calc_trees xfr_aligns_trees are in pandagma-common.sh
export commandlist="ingest fam_consen search_families realign_and_trim calc_trees xfr_aligns_trees summarize"

export dependencies='hmmscan famsa'

declare out_dir version consen_iden clust_cov consen_method annot_str_regex

main_pan_fam "$@"

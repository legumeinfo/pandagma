#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates 
# pan-gene clusters using the programs mmseqs, dagchainer, and mcl. 
# Authors: Steven Cannon, Hyunoh Lee, Joel Berendzen, Nathan Weeks, 2020-2023
#
scriptname='pandagma pan'

. pandagma-common.sh

define HELP_DOC <<'EOS'
Compute pan-gene clusters using a combination of synteny and homology,
using the programs mmseqs, dagchainer, and mcl, and additional pre- and post-refinement steps.

Usage:
  $scriptname -c CONFIG_FILE [options]

  Required:
           -c (path to the config file)

  Options: -s (subcommand to run. If \"all\" or omitted, all steps will be run; otherwise, run specified step)
           -w (working directory, for temporary and intermediate files [default: './pandagma_work'].)
           -o OUTPUT_DIR (name for the output directory [default: './pandagma_out'].
                Applicable only to "all" and "summarize" steps.)
           -n (number of processors to use. Defaults to number of processors available to the parent process)
           -O (ordering method, for placing pan-genes. Options: 
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
                all - All of the steps below, except for clean 
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
     calc_chr_pairs - Report observed chromosome pairs; useful for preparing expected_chr_matches
          summarize - Copy results into output directory, and report summary statistics.

  Run the following subcommands separately if you wish:
              align - Align families.
     model_and_trim - Build HMMs and trim the alignments, preparatory to calculating trees.
         calc_trees - Calculate gene trees.
   xfr_aligns_trees - Transfer alignments, HMMs, and trees to output directory

              clean - Clean (delete) files in the working directory that are not needed 
                        for later addition of data using add_extra and subsequent run commands.
                        By default, \"clean\" is run as part of \"all\" unless the -r flag is set.
''
EOS

define MORE_INFO <<'EOS'
Optionally, an expected_chr_matches array variable can be specified in the config file,
which provides anticipated chromosome pairings, e.g.

expected_chr_matches=(
  01 01
  02 02
  ...
# allows for translocation between 11 and 13
  11 13 
  12 12
)

These pairings are used in a regular expression to identify terminal portions of molecule IDs, e.g.
  glyma.Wm82.gnm2.Gm01  glyso.PI483463.gnm1.Gs01
  glyma.Wm82.gnm2.Gm13  glyso.W05.gnm1.Chr11
If an expected_chr_matches variable is not defined in the pandagma config file, then no such filtering will be done.

Variables in pandagma config file (Set the config with the CONF environment variable)
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.60]
         extra_iden - Minimum identity threshold for mmseqs addition of \"extra\" annotations [90]
      mcl_inflation - Inflation parameter, for Markov clustering [1.6]
        strict_synt - For clustering of the \"main\" annotations, use only syntenic pairs [1]
                        The alternative (0) is to use all homologous pairs that satisfy expected_chr_matches
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
''
EOS

########################################
# Helper functions begin here
canonicalize_paths() {
  echo "Entering canonicalize_paths. Fasta files: "
  echo "${cds_files[@]}"

  cd "${DATA_DIR}" || exit

  mapfile -t cds_files < <(realpath --canonicalize-existing "${cds_files[@]}")
  mapfile -t annotation_files < <(realpath --canonicalize-existing "${annotation_files[@]}")
  mapfile -t protein_files < <(realpath --canonicalize-existing "${protein_files[@]}")

  if [[ -v cds_files_extra_constr ]]
  then
    mapfile -t cds_files_extra_constr < <(realpath --canonicalize-existing "${cds_files_extra_constr[@]}")
    mapfile -t annotation_files_extra_constr < <(realpath --canonicalize-existing "${annotation_files_extra_constr[@]}")
    mapfile -t protein_files_extra_constr < <(realpath --canonicalize-existing "${protein_files_extra_constr[@]}")
  fi  

  if [[ -v cds_files_extra_free ]]
  then
    mapfile -t cds_files_extra_free < <(realpath --canonicalize-existing "${cds_files_extra_free[@]}")
    mapfile -t annotation_files_extra_free < <(realpath --canonicalize-existing "${annotation_files_extra_free[@]}")
    mapfile -t protein_files_extra_free < <(realpath --canonicalize-existing "${protein_files_extra_free[@]}")
  fi  

  cd "${OLDPWD}" || exit
  readonly submit_dir=${PWD}

  fasta_file=$(basename "${cds_files[0]}" .gz)
  fa="${fasta_file##*.}"
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
    start_time=$(date)
    printf "Run started at: %s\n" "$start_time" > stats/tmp.timing

  echo "  Get position information from the main annotation sets."
  cat /dev/null > 02_all_main_cds.fna # Collect all starting sequences, for later comparisons
  for (( file_num = 0; file_num < ${#cds_files[@]} ; file_num++ )); do
    file_base=$(basename "${cds_files[file_num]%.*}")
    cat_or_zcat "${cds_files[file_num]}" >> 02_all_main_cds.fna # Collect original seqs for later comparisons
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | 
      gff_or_bed_to_hash5.awk > 01_posn_hsh/"$file_base".hsh
    hash_into_fasta_id.pl -nodef -fasta "${cds_files[file_num]}" \
                          -hash 01_posn_hsh/"$file_base".hsh \
                          -out 02_fasta_nuc/"$file_base"
    # calc basic sequence stats
    annot_name=$(basename 02_fasta_nuc/"$file_base" | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
    printf "  Main:  " >> stats/tmp.fasta_seqstats
    cat_or_zcat "${cds_files[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
  done

  echo "  Get position information from the extra annotation sets, if any,"
  echo "  for the annotations that will be added chromosome constraints (_constr)"
  cat /dev/null > 02_all_extra_cds.fna # Collect all starting sequences, for later comparisons
  if [[ -v cds_files_extra_constr ]]
  then
    for (( file_num = 0; file_num < ${#cds_files_extra_constr[@]} ; file_num++ )); do
      file_base=$(basename "${cds_files_extra_constr[file_num]%.*}")
      cat_or_zcat "${cds_files_extra_constr[file_num]}" >> 02_all_extra_cds.fna  
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra_constr[file_num]}" | 
        gff_or_bed_to_hash5.awk > 01_posn_hsh/"$file_base".hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files_extra_constr[file_num]}" \
                            -hash 01_posn_hsh/"$file_base".hsh \
                            -out 02_fasta_nuc/"$file_base"
      # calc basic sequence stats
      annot_name=$(basename 02_fasta_nuc/"$file_base" | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )

      printf "  Extra: " >> stats/tmp.fasta_seqstats
      cat_or_zcat "${cds_files_extra_constr[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
    done
  fi

  echo "  Get position information from the extra annotation sets, if any,"
  echo "  for the annotations that will be added without chromosome constraints (_free)"
  if [[ -v cds_files_extra_free ]]
  then
    for (( file_num = 0; file_num < ${#cds_files_extra_free[@]} ; file_num++ )); do
      file_base=$(basename "${cds_files_extra_free[file_num]%.*}")
      cat_or_zcat "${cds_files_extra_free[file_num]}" >> 02_all_extra_cds.fna  
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra_free[file_num]}" | 
        gff_or_bed_to_hash5.awk > 01_posn_hsh/"$file_base".hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files_extra_free[file_num]}" \
                            -hash 01_posn_hsh/"$file_base".hsh \
                            -out 02_fasta_nuc/"$file_base"
      # calc basic sequence stats
      annot_name=$(basename 02_fasta_nuc/"$file_base" | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )

      printf "  Extra: " >> stats/tmp.fasta_seqstats
      cat_or_zcat "${cds_files_extra_free[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
    done
  fi

  cat /dev/null > 01_combined_posn.hsh # Collect combined file of position hashes
    echo "  Make a combined file of position hashes for later use in step add_extra"
      cat 01_posn_hsh/*.hsh >> 01_combined_posn.hsh

  echo "  Count starting sequences, for later comparisons"
  for file in 02_fasta_nuc/*."$fa"; do
    awk '$0~/UNDEFINED/ {ct++} 
      END{if (ct>0){print "Warning: " FILENAME " has " ct " genes without position (HASH UNDEFINED)" } }' "$file"
    grep '>' "$file" | perl -pe 's/__/\t/g' | cut -f2 | # extracts the GeneName from combined genspchr__GeneName__start__end__orient
      perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex.+/$1/' |
      grep -v UNDEFINED | sort | uniq -c | awk '{print $2 "\t" $1}' >> stats/tmp.gene_count_start
  done

  echo "  Also get protein files"
  if [[ -v protein_files ]]
  then
    for (( file_num = 0; file_num < ${#protein_files[@]} ; file_num++ )); do
      echo "  Copying protein file ${protein_files[file_num]}"
      cp "${protein_files[file_num]}" 02_fasta_prot/
    done
  fi
  if [[ -v protein_files_extra_constr ]]
  then
    for (( file_num = 0; file_num < ${#protein_files_extra_constr[@]} ; file_num++ )); do
      echo "  Copying protein file ${protein_files_extra_constr[file_num]}"
      cp "${protein_files_extra_constr[file_num]}" 02_fasta_prot/
    done
  fi
  if [[ -v protein_files_extra_free ]]
  then
    for (( file_num = 0; file_num < ${#protein_files_extra_free[@]} ; file_num++ )); do
      echo "  Copying protein file ${protein_files_extra_free[file_num]}"
      cp "${protein_files_extra_free[file_num]}" 02_fasta_prot/
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
    qry_base=$(basename "${cds_files[file1_num]%.*}" ."$fa")
    for (( file2_num = file1_num + 1; file2_num < ${#cds_files[@]} ; file2_num++ )); do
      sbj_base=$(basename "${cds_files[file2_num]%.*}" ."$fa")
      echo "  Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      { cat 02_fasta_nuc/"$qry_base"."$fa" 02_fasta_nuc/"$sbj_base"."$fa" ; } |
        mmseqs easy-cluster stdin 03_mmseqs/"${qry_base}".x."${sbj_base}" "$MMTEMP" \
         --min-seq-id "$clust_iden" -c "$clust_cov" --cov-mode 0 --cluster-reassign 1>/dev/null
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
  if [[ -v expected_chr_matches ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns defined in expected_chr_matches"
    for mmseqs_path in 03_mmseqs/*_cluster.tsv; do
      outfilebase=$(basename "$mmseqs_path" _cluster.tsv)
      echo "  $outfilebase"
      filter_mmseqs_by_chroms.pl -chr_pat <(printf '%s %s\n' "${expected_chr_matches[@]}") < "${mmseqs_path}" |
        awk 'NF==8' |  # matches for genes with coordinates. The case of <8 can happen for seqs with UNDEF position.
        cat > 04_dag/"${outfilebase}"_matches.tsv &

      # allow to execute up to $NPROC in parallel
      if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
    done
    wait # wait for last jobs to finish
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches variable was defined in the config file, so proceeding without chromosome-pair filtering."
    for mmseqs_path in 03_mmseqs/*_cluster.tsv; do
      outfilebase=$(basename "$mmseqs_path" _cluster.tsv)
      perl -pe 's/__/\t/g; s/\t[\+-]//g' "${mmseqs_path}" |
        awk 'NF==8' |  # matches for genes with coordinates. The case of <8 can happen for seqs with UNDEF position.
        cat > 04_dag/"${outfilebase}"_matches.tsv 
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
    align_file=$(basename "$match_path" _matches.tsv)
    qryfile=$(echo "$align_file" | perl -pe 's/(\S+)\.x\..+/$1/')
    sbjfile=$(echo "$align_file" | perl -pe 's/\S+\.x\.(\S+)/$1/')

    echo "Find average distance between genes for the query and subject files: "
    echo "  $qryfile and $sbjfile"
    ave_gene_gap=$(cat 02_fasta_nuc/"$qryfile"."$fa" 02_fasta_nuc/"$sbjfile"."$fa" | 
                     awk '$1~/^>/ {print substr($1,2)}' | perl -pe 's/__/\t/g' | sort -k1,1 -k3n,3n |
                     awk '$1 == prev1 && $3 > prev4 {sum+=$3-prev4; ct++; prev1=$1; prev3=$3; prev4=$4};
                          $1 != prev1 || $3 <= prev4 {prev1=$1; prev3=$3; prev4=$4}; 
                          END{print 100*int(sum/ct/100)}')
    max_gene_gap=$(( ave_gene_gap * 20 ))

    echo "Running DAGchainer on comparison: $align_file"
    echo "  Calculated DAGchainer parameters: -g (ave_gene_gap): $ave_gene_gap -D (max_gene_gap): $max_gene_gap"; echo
    echo " run_DAG_chainer.pl $dagchainer_args  -g $ave_gene_gap -D $max_gene_gap -i \"${OLDPWD}/${match_path}\""

    # run_DAG_chainer.pl writes temp files to cwd;
    # use per-process temp directory to avoid any data race
    (
      tmpdir=$(mktemp -d)
      cd "${tmpdir}"
      run_DAG_chainer.pl "$dagchainer_args"  -g "$ave_gene_gap" -D "$max_gene_gap" -i "${OLDPWD}/${match_path}" 1>/dev/null
      rmdir "${tmpdir}"
    ) &
    # allow to execute up to $NPROC in parallel
    [ "$(jobs -r -p | wc -l)" -ge "${NPROC}" ] && wait -n
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
  printf "\nDo Markov clustering with inflation parameter %s and %d threads\n" "$mcl_inflation" "${NPROC}"
  echo "MCL COMMAND: mcl 05_filtered_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv"
  mcl 05_filtered_pairs.tsv -I "$mcl_inflation" -te "${NPROC}" --abc -o tmp.syn_pan.clust.tsv \
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
  echo "  < 07_pan_fasta_cds.fna pick_family_rep.pl -prefer $preferred_annot -out 08_pan_fasta_clust_rep_cds.fna"
  < 07_pan_fasta_cds.fna pick_family_rep.pl -prefer "$preferred_annot" -out 08_pan_fasta_clust_rep_cds.fna 

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
                     --search-type "${SEQTYPE}" --cov-mode 5 -c "${clust_cov}" 1>/dev/null 

  echo "  Filter unclust.x.07_pan_fasta.m8 by clust_iden and report: panID, qry_gene, sbj_gene"
  top_line.awk 10_place_leftovers/unclust.x.07_pan_fasta.m8 |
    awk -v IDEN="${clust_iden}" '$3>=IDEN {print $2 "\t" $1}' | 
    perl -pe 's/^(pan\d+)__(\S+)\t(\S)/$1\t$2\t$3/' > 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj.tsv

  echo "  Add gene-pair information"
  hash_into_table_2cols.pl 01_posn_hsh/*hsh -table 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj.tsv |
    cat > 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_posn.tsv

  echo; echo "  Filter based on list of expected chromosome pairings if provided"
  if [[ -v expected_chr_matches ]]; then  
    echo "Filtering on chromosome patterns defined in expected_chr_matches"
    < 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_posn.tsv \
      filter_mmseqs_by_chroms.pl -chr_pat <(printf '%s %s\n' "${expected_chr_matches[@]}") |
      awk 'NF==9' |  # matches for genes with coordinates. The case of <9 can happen for seqs with UNDEF position.
      cat > 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_matches.tsv 
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches variable was defined in the config file, so proceeding without chromosome-pair filtering."
    perl -pe 's/__/\t/g; s/\t[\+-]//g' < 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_posn.tsv |
      awk 'NF==9' |  # matches for genes with coordinates. The case of <9 can happen for seqs with UNDEF position.
      cat > 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_matches.tsv
  fi

  echo "  Place unclustered genes into their respective pan-gene sets"
  cut -f1,3 10_place_leftovers/unclust.x.07_pan_fasta.pan_qry_sbj_matches.tsv | 
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
  hash_to_rows_by_1st_col.awk < 12_syn_pan_aug_pre.hsh.tsv > 12_syn_pan_aug_pre.clust.tsv
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
   < 12_syn_pan_aug_complement.fna mmseqs easy-cluster stdin 03_mmseqs/$complmt_self_compare "$MMTEMP" \
    --min-seq-id "$clust_iden" \
    -c "$clust_cov" \
    --cov-mode 0 \
    --cluster-reassign 1>/dev/null

  echo "  Add gene-pair information"
  hash_into_table_2cols.pl 01_posn_hsh/*hsh -idx1 0 -idx2 1 -table 03_mmseqs/${complmt_self_compare}_cluster.tsv |
    cat > 03_mmseqs/${complmt_self_compare}_cluster_posn.tsv

  if [[ -v expected_chr_matches ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns defined in expected_chr_matches"
    < 03_mmseqs/${complmt_self_compare}_cluster_posn.tsv filter_mmseqs_by_chroms.pl \
      -chr_pat <(printf '%s %s\n' "${expected_chr_matches[@]}") > 04_dag/${complmt_self_compare}_matches.tsv
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches was defined in the config file, so proceeding without chromosome-pair filtering."
    cat 03_mmseqs/${complmt_self_compare}_cluster.tsv |
      perl -pe 's/__/\t/g; s/\t[\+-]$//' > 04_dag/${complmt_self_compare}_matches.tsv 
  fi

  echo "  Extract query and subject matches to be clustered"
  cat 04_dag/${complmt_self_compare}_matches.tsv | cut -f2,6 > 04_dag/${complmt_self_compare}_pairs.tsv

  echo "  Cluster the remaining sequences that have matches"
  mcl  04_dag/${complmt_self_compare}_pairs.tsv -I "$mcl_inflation" -te "${NPROC}" --abc -o tmp.syn_pan_aug_complement.clust.tsv

  echo "  Find number of clusters in initial (06) results"
  clust_count_06=$(wc -l 06_syn_pan.clust.tsv | awk '{print $1}')
 
  echo "  Add cluster IDs"
  awk -v START="$clust_count_06" 'NF>1 {padnum=sprintf("%05d", NR+START); print "pan" padnum "\t" $0}' tmp.syn_pan_aug_complement.clust.tsv |
    cat > 12_syn_pan_aug_complement.clust.tsv
  rm tmp.syn_pan_aug_complement.clust.tsv

  echo "  Combine clusters derived from synteny or restricted homology (12_syn_pan_aug.clust.tsv)"
  echo "  and clusters from the complement of that set. Skip singletons (\"clusters\" of size 1)."
  cat 12_syn_pan_aug_pre.clust.tsv 12_syn_pan_aug_complement.clust.tsv | awk 'NF>1' > 12_syn_pan_aug.clust.tsv

  echo "  Reshape from hash to mcl output format (clustered IDs on one line)"
  hash_to_rows_by_1st_col.awk 12_syn_pan_aug.clust.tsv > 12_syn_pan_aug.hsh.tsv
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

  if [[ -v cds_files_extra_constr ]] || [[ -v cds_files_extra_free ]]
  then # handle the "extra" annotation files
    echo "  Search non-clustered genes against pan-gene consensus sequences"
    # Check sequence type (in case this run function is called separately from the usually-prior ones)
    someseq=$(head 07_pan_fasta_cds.fna | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
    SEQTYPE=$(check_seq_type "${someseq}") # 3=nuc; 1=pep
    echo "SEQTYPE is: $SEQTYPE"

    for file in "${cds_files_extra_constr[@]}"; do
      file_base=$(basename "$file" .gz)
      echo "Extra: 02_fasta_nuc/$file_base"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      mmseqs easy-search "02_fasta_nuc/$file_base" 13_pan_aug_fasta_posn.fna 13_extra_out_constr_dir/"${file_base}".x.all_cons.m8 \
                   "$MMTEMP" --search-type "${SEQTYPE}" --cov-mode 5 -c "${clust_cov}" 1>/dev/null
    done
    wait # wait for jobs to finish

    for file in "${cds_files_extra_free[@]}"; do
      file_base=$(basename "$file" .gz)
      echo "Extra: 02_fasta_nuc/$file_base"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      mmseqs easy-search "02_fasta_nuc/$file_base" 13_pan_aug_fasta_posn.fna 13_extra_out_free_dir/"${file_base}".x.all_cons.m8 \
                   "$MMTEMP" --search-type "${SEQTYPE}" --cov-mode 5 -c "${clust_cov}" 1>/dev/null
    done
    wait # wait for jobs to finish

    if [[ -v expected_chr_matches ]]; then  # filter based on list of expected chromosome pairings if provided
      echo "Filtering on chromosome patterns defined in expected_chr_matches and by identity and top match."
      echo "First find top match without respect to chromosome match, then with chromosome matches."
      echo "A top hit with chromosome match trumps a top hit on the wrong chromosome."
      echo "The placement with chromosome match is prefixed with 1 and that without is prefixed with 2;"
      echo "then the first (top) gene is selected - i.e., 1 geneA before 2 geneA."
      if [[ -v cds_files_extra_constr ]]; then
        for m8file in 13_extra_out_constr_dir/*.m8; do
          base=$(basename "$m8file" .m8)
          echo "  Processing file $m8file"
          top_line.awk "$m8file" | awk -v IDEN="${extra_iden}" '$3>=IDEN {print $1 "\t" $2}' |
              perl -pe 's/^\S+__(\S+)__\d+__\d\S*\t\S+__(\S+)__\d+__\d\S+$/2\t$1\t$2/' > 13_extra_out_constr_dir/"$base".top_2nd
          filter_mmseqs_by_chroms.pl -chr_pat <(printf '%s %s\n' "${expected_chr_matches[@]}") < "$m8file" |
              top_line.awk | awk -v IDEN="${extra_iden}" '$3>=IDEN {print $1 "\t" $2}' | 
              perl -pe 's/^\S+__(\S+)__\d+__\d+\t\S+__(\S+)__\d+__\d+$/1\t$1\t$2/' > 13_extra_out_constr_dir/"$base".top_1st
        done
      fi
    else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
      echo "WARNING: No expected_chr_matches.tsv file was provided, but annotations were indicated in the "
      echo "file sets annotation_files_extra_constr and cds_files_extra_constr. Please check that "
      echo "expected_chr_matches.tsv is in the data directory OR that extra annotations are in the "
      echo "annotation_files_extra_free and cds_files_extra_free file collections."
      exit 1;
    fi

    if [[ -v cds_files_extra_free ]]; then
      for m8file in 13_extra_out_free_dir/*.m8; do
        base=$(basename "$m8file" .m8)
        echo "  Processing file $m8file"
        top_line.awk "$m8file" | awk -v IDEN="${extra_iden}" '$3>=IDEN {print $1 "\t" $2}' |
            perl -pe 's/^\S+__(\S+)__\d+__\d\S*\t\S+__(\S+)__\d+__\d\S+$/1\t$1\t$2/' > 13_extra_out_free_dir/"$base".top_2nd
      done
    fi

    if [ "$strict_synt" -eq 1 ]; then
      echo "Join panIDs to unclustered genes. This RESTRICTIVE filtering requires that genes have a chromosome match."
      cat 13_extra_out_*_dir/*top_* | awk '$1==1' | # restrict to those with chromosome matches
        sort -k3,3 -k2,2 -k1,1 | cut -f2,3 | join -1 2 -2 1 - 13_pan_gene.hsh | 
        awk 'NF==3 {print $3 "\t" $2}' | sort -k1,1 -k2,2 | uniq > 14_syn_pan_extra.hsh.tsv
    else
      echo "Join panIDs to unclustered genes. This PERMISSIVE filtering allows gene placement for genes that lack a chromosome match."
      cat 13_extra_out_*_dir/*top_* | sort -k3,3 -k2,2 -k1,1 | cut -f2,3 | join -1 2 -2 1 - 13_pan_gene.hsh | 
        awk 'NF==3 {print $3 "\t" $2}' | sort -k1,1 -k2,2 | uniq > 14_syn_pan_extra.hsh.tsv
    fi

    echo "Derive a cluster-format file from the hash of panIDs and extra genes"
    hash_to_rows_by_1st_col.awk 14_syn_pan_extra.hsh.tsv  > 14_syn_pan_extra.clust.tsv

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
    hash_to_rows_by_1st_col.awk 18_syn_pan_aug_extra.hsh.tsv > 18_syn_pan_aug_extra.clust.tsv

    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    echo "    Fasta file:" "${protein_files[@]} ${protein_files_extra_constr[@]} ${protein_files_extra_free[@]}"
    if [ -d 19_pan_aug_leftover_merged_cds ]; then rm -rf 19_pan_aug_leftover_merged_cds; fi
    mkdir -p 19_pan_aug_leftover_merged_cds
    get_fasta_from_family_file.pl "${cds_files[@]}" "${cds_files_extra_constr[@]}" "${cds_files_extra_free[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_cds
  
    echo "  Merge files in 19_pan_aug_leftover_merged_cds, prefixing IDs with panID__"
    find 19_pan_aug_leftover_merged_cds -maxdepth 1 -mindepth 1 |
      sort |
        xargs awk '/^>/ {printf(">%s__%s\n", substr(FILENAME,index(FILENAME,"/")+1), substr($0,2)); next} {print}' > 19_pan_aug_leftover_merged_cds.fna

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
    zcat "$filepath" >> 20_pan_fasta_prot.faa
  done

  echo "  Get protein sequences into pan-gene sets, corresponding with 18_syn_pan_aug_extra.clust.tsv"
  if [ -d 19_pan_aug_leftover_merged_prot ]; then rm -rf 19_pan_aug_leftover_merged_prot; fi
  mkdir -p 19_pan_aug_leftover_merged_prot
  get_fasta_from_family_file.pl 20_pan_fasta_prot.faa \
              -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_prot

  echo "  Get all protein sequences from files in 19_pan_aug_leftover_merged_prot"
  find 19_pan_aug_leftover_merged_prot -maxdepth 1 -mindepth 1 |
    sort |
      xargs awk '/^>/ {printf(">%s__%s\n", substr(FILENAME,index(FILENAME,"/")+1), substr($0,2)); next} {print}' > 19_pan_aug_leftover_merged_prot.faa

  echo "  Pick a representative sequence for each pangene set - as a sequence with the median length for that set."
  echo "    == first proteins:"
  pick_family_rep.pl -nostop -prefer "$preferred_annot" \
    -out 21_pan_fasta_clust_rep_prot.faa < 19_pan_aug_leftover_merged_prot.faa

  echo "    == then CDS sequences, corresponding with 21_pan_fasta_clust_rep_prot.faa"
  awk '$1~/^>/ {print substr($1,2)}' 21_pan_fasta_clust_rep_prot.faa > lists/lis.21_pan_fasta_clust_rep
  get_fasta_subset.pl -in 19_pan_aug_leftover_merged_cds.fna -list lists/lis.21_pan_fasta_clust_rep \
                    -clobber -out 21_pan_fasta_clust_rep_cds.fna

  perl -pi -e 's/__/  /' 21_pan_fasta_clust_rep_cds.fna
  perl -pi -e 's/__/  /' 21_pan_fasta_clust_rep_prot.faa

  echo "  Retrieve genes present in the original CDS files but absent from 18_syn_pan_aug_extra"
  cut -f2 18_syn_pan_aug_extra.hsh.tsv | LC_ALL=C sort > lists/lis.18_syn_pan_aug_extra
  cat 02_all_*_cds.fna > 02_all_cds.fna
  get_fasta_subset.pl -in 02_all_cds.fna -out 18_syn_pan_aug_extra_complement.fna \
    -lis lists/lis.18_syn_pan_aug_extra -xclude -clobber
}

##########
run_tabularize() {
  echo; echo "== Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv"
  cd "${WORK_DIR}"

  # Get table header 
  pangene_tabularize.pl -pan 18_syn_pan_aug_extra.clust.tsv -annot_str_regex "$ANN_REX" > tmp.18_syn_pan_aug_extra.clust.tsv
  head -1 tmp.18_syn_pan_aug_extra.clust.tsv > tmp.table_header

  # Find column on which to sort first
  export PR_ANN=$preferred_annot
  sort_col=$( perl -ane 'BEGIN{use List::Util qw(first); 
                $PA=$ENV{"PR_ANN"}}; 
                $idx = first { $F[$_] =~ /$PA/ } 0..$#F; print ($idx+1) ' < tmp.table_header )

  echo "  Field for sorting 18_syn_pan_aug_extra.table.tsv is sort_col: $sort_col"

  echo "  Sort, putting header row at top, and don't print pangenes that are all NONE"
    sort -k"$sort_col","$sort_col" -k2,2 tmp.18_syn_pan_aug_extra.clust.tsv | 
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
  calc_pan_stats.pl -annot_regex "$ANN_REX" -pan 18_syn_pan_aug_extra.clust.tsv -out 18_syn_pan_aug_extra.counts.tsv
  max_annot_ct=$( awk '$1!~/^#/ {print $2}' 18_syn_pan_aug_extra.counts.tsv | sort -n | uniq | tail -1)

  echo "  Select orthogroups with genes from selected percentiles of max_annot_ct annotation sets"
  for percentile in $pctl_low $pctl_med $pctl_hi; do
    awk -v PCTL="$percentile" -v ANNCT="$max_annot_ct" '$2>=ANNCT*(PCTL/100)&& $1!~/^#/ {print $1}' 18_syn_pan_aug_extra.counts.tsv |
        cat > lists/lis.18_syn_pan_aug_extra.pctl"${percentile}"

    echo "  Get a fasta subset with only genes from at least $((max_annot_ct*pctl_low/100)) annotation sets"
    get_fasta_subset.pl -in 21_pan_fasta_clust_rep_cds.fna \
                        -list lists/lis.18_syn_pan_aug_extra.pctl"${percentile}" \
                        -clobber -out 22_pan_fasta_rep_pctl"${percentile}"_cds.fna

    echo "  Get a clust.tsv file with orthogroups with at least min_core_prop*max_annot_ct annotation sets"
    join <(LC_ALL=C sort -k1,1 lists/lis.18_syn_pan_aug_extra.pctl"${percentile}") \
         <(LC_ALL=C sort -k1,1 18_syn_pan_aug_extra.clust.tsv) |
            cat > 22_syn_pan_aug_extra_pctl"${percentile}".clust.tsv
  done
}

##########
run_order_and_name() {
  cd "${WORK_DIR}"

  echo "  Reshape from mcl output format, clustered IDs on one line, to a hash format"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 22_syn_pan_aug_extra_pctl"${pctl_low}".clust.tsv |
    cat > 22_syn_pan_aug_extra_pctl"${pctl_low}".hsh.tsv

  echo "  Add positional information to the hash output."
  join -a 1 -1 2 -2 1 <(sort -k2,2 22_syn_pan_aug_extra_pctl"${pctl_low}".hsh.tsv) <(cat 01_posn_hsh/*hsh | sort -k1,1) | 
    perl -pe 's/__/\t/g; s/ /\t/g' | awk -v OFS="\t" '{print $2, $1, $3, $5, $6, $7}' |
    sort -k1,1 -k2,2 > 22_syn_pan_aug_extra_pctl"${pctl_low}"_posn.hsh.tsv

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
    order_encode.pl -pan_table 22_syn_pan_aug_extra_pctl"${pctl_low}"_posn.hsh.tsv \
                          -code_table code_table/pan_to_peptide.tsv -utilized code_table/utilized.tsv \
                          -annot_regex "${ANN_REX}" -outdir 23_encoded_chroms 
  
    echo "  Align chromosome sequences with peptide-encoded-gene-orders."
    mkdir -p 23_encoded_chroms_aligned
    do_align.sh 23_encoded_chroms 23_encoded_chroms_aligned $(( NPROC/2 ))
  
    echo "  Filter chromosome alignments to consensus motifs"
    mkdir -p 23_encoded_chroms_filt1
    for filepath in 23_encoded_chroms_aligned/*; do 
      base=$(basename "$filepath")
      cons -sequence "$filepath" -outseq 23_encoded_chroms_filt1/"$base" -name "$base" &
      if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
    done
    wait
  
    echo "  Decode motifs to recover pangene IDs in order along each chromosome."
    cat /dev/null > consen_gene_order.tsv
    for filepath in 23_encoded_chroms_filt1/*; do
      order_decode.pl -align "$filepath" -code code_table/utilized.tsv -verbose |
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
    order_by_reference.pl -pan_table 22_syn_pan_aug_extra_pctl"${pctl_low}"_posn.hsh.tsv \
                          -annot_regex "${ANN_REX}" -pref_annot "$preferred_annot" \
                          -consen_out consen_gene_order.tsv -unplaced_out consen_pan_unplaced.txt -verbose 
  else 
    echo "Ordering method (-o $order_method ) is not one of \"alignment\" or \"reference\". Aborting."
    exit 1
  fi # End of $order_method =~ /align/

  gapfill_threads=$((NPROC/2))
  echo "  Fill gaps in the alignment-based pangene ordering. This step is time-consuming."
  echo "  so is run in parallel, using $gapfill_threads threads."
  echo "   order_gapfill.pl -verbose -annot_regex \"${ANN_REX}\" -consen consen_gene_order.tsv \\";
  echo "     -pan_table 22_syn_pan_aug_extra_pctl${pctl_low}_posn.hsh.tsv \\";
  echo "     -unplaced consen_pan_unplaced.txt -out consen_gene_order_gapfilled.tsv -nproc $gapfill_threads";
  order_gapfill.pl -verbose -annot_regex "${ANN_REX}" -consen consen_gene_order.tsv \
    -pan_table 22_syn_pan_aug_extra_pctl"${pctl_low}"_posn.hsh.tsv \
    -unplaced consen_pan_unplaced.txt -out consen_gene_order_gapfilled.tsv -nproc $gapfill_threads

  echo "  Reshape defline into a hash, e.g. pan20175 Phaseolus.pan2.chr11_222300_pan20175 222300 223300 -"
  echo "  Note: these \"positions\" and sizes are artificial, representing only inferred relative positions."
  awk -v PRE="$consen_prefix" '
    {printf("%s\t%s.%s_%06d_%s\t%d\t%d\t%s\n", $1, PRE, $2, $3, $1, $3*100, $3*100+1000, $4)}' consen_gene_order_gapfilled.tsv > consen_posn.hsh

  echo "  Extract a version of the consen_gene_order file with immediate tandem duplicates collapsed"
  top_line.awk consen_posn.hsh > consen_posn_trim.hsh

  echo "  Hash position information into fasta file (CDS)"
  hash_into_fasta_id.pl -fasta 22_pan_fasta_rep_pctl"${pctl_low}"_cds.fna -hash consen_posn_trim.hsh |
    grep -v "HASH UNDEFINED" > 22_pan_fasta_rep_pctl"${pctl_low}"_posn_cdsTMP.fna

  echo "  Reshape defline, and sort by position (CDS)"
  fasta_to_table.awk 22_pan_fasta_rep_pctl"${pctl_low}"_posn_cdsTMP.fna | sed 's/__/\t/g; s/ /\t/' | 
    sort | awk '{print ">" $1, $2, $3, $4, $5; print $6}' > 23_syn_pan_pctl"${pctl_low}"_posn_cds.fna

  echo "  Also get corresponding protein sequences for genes in lists/lis.18_syn_pan_aug_extra.pctl${pctl_low}"
  if [ ! -f 20_pan_fasta_prot.faa ]; then
    for filepath in 02_fasta_prot/*.gz; do 
      zcat "$filepath" >> 20_pan_fasta_prot.faa
    done
  fi

  awk '$1~/^>/ {print $5}' 23_syn_pan_pctl"${pctl_low}"_posn_cds.fna > lists/lis.23_syn_pan_pctl"${pctl_low}"_posn
  get_fasta_subset.pl -in 20_pan_fasta_prot.faa -list lists/lis.23_syn_pan_pctl"${pctl_low}"_posn \
                      -clobber -out 23_syn_pan_pctl"${pctl_low}"_posn_proteinTMP.faa

  echo "  Get directory of protein multifasta sequences for each pangene, for (separate) protein alignments"

  if [ -d 22_syn_pan_aug_extra_pctl"${pctl_low}" ]; then rm -rf 22_syn_pan_aug_extra_pctl"${pctl_low}" ; fi
  mkdir -p 22_syn_pan_aug_extra_pctl"${pctl_low}"
  get_fasta_from_family_file.pl 20_pan_fasta_prot.faa \
    -family_file 22_syn_pan_aug_extra_pctl"${pctl_low}".clust.tsv -out_dir 22_syn_pan_aug_extra_pctl"${pctl_low}"

  echo "  Hash pan-ID into protein fasta file"
  awk '$1~/^>/ { print $5 "\t" substr($1,2) "__" $2 "__" $3 "__" $4 }' 23_syn_pan_pctl"${pctl_low}"_posn_cds.fna |
    cat > lists/23_syn_pan_pctl"${pctl_low}"_posn.hsh
  hash_into_fasta_id.pl -hash lists/23_syn_pan_pctl"${pctl_low}"_posn.hsh -swap_IDs -nodef \
    -fasta 23_syn_pan_pctl"${pctl_low}"_posn_proteinTMP.faa |
    perl -pe 's/__/ /g' > 23_syn_pan_pctl"${pctl_low}"_posn_prot.faa
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
    file=$(basename "$filepath");
    echo "  Computing alignment, using program famsa, for file $file"
    famsa -t 2 19_pan_aug_leftover_merged_cds/"$file" 20_aligns_cds/"$file" 1>/dev/null &
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
    get_fasta_from_family_file.pl "${protein_files[@]}" "${protein_files_extra_constr[@]}" "${protein_files_extra_free[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_prot
  fi

  echo; echo "== Align protein sequence for the each gene family =="
  mkdir -p 20_aligns_prot
  for filepath in 19_pan_aug_leftover_merged_prot/*; do
    file=$(basename "$filepath");
    echo "  Computing alignment, using program famsa, for file $file"
    famsa -t 2 19_pan_aug_leftover_merged_prot/"$file" 20_aligns_prot/"$file" 1>/dev/null &
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
    file=$(basename "$filepath");
    hmmbuild -n "$file" 21_hmm/"$file" "$filepath" 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Realign to HMMs =="
  mkdir -p 22_hmmalign
  for filepath in 21_hmm/*; do
    file=$(basename "$filepath");
    printf "%f " "$file"
    hmmalign --trim --outformat A2M -o 22_hmmalign/"$file" 21_hmm/"$file" 19_pan_aug_leftover_merged_cds/"$file" &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo

  echo; echo "== Trim HMM alignments to match-states =="
  mkdir -p 23_hmmalign_trim1
  for filepath in 22_hmmalign/*; do
    file=$(basename "$filepath");
    printf "%s " "$file"
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
    printf "%f " "$file"
    filter_align.pl -in "$filepath" -out 23_hmmalign_trim2/"$file" -log 23_hmmalign_trim2_log/"$file" \
                    -depth $min_depth -pct_depth $min_pct_depth -min_pct_aligned $min_pct_aligned &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo
}

##########
run_xfr_aligns_trees() {
  echo; echo "Copy alignment and tree results into output directory"

  full_out_dir="${out_dir}"

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p "$full_out_dir"
  fi

  for dir in 20_aligns_cds 20_aligns_prot 21_hmm 22_hmmalign 23_hmmalign_trim2 24_trees; do
    if [ -d "${WORK_DIR}"/$dir ]; then
      echo "Copying directory $dir to output directory"
      cp -r "${WORK_DIR}"/$dir "${full_out_dir}"/
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
    --min-seq-id "$clust_iden" -c "$clust_cov" --cov-mode 0 --cluster-reassign 1>/dev/null

  echo "   Extract chromosome-chromosome correspondences"
  cut -f2,3 22_syn_pan_aug_extra_pctl"${pctl_low}"_posn.hsh.tsv | sort -k1,1 > 22_syn_pan_aug_extra_pctl"${pctl_low}"_posn.gene_chr.hsh

  echo "   From mmseqs cluster table, prune gene pairs, keeping only those in the same pangene cluster"
  perl -pe 's/__/\t/g' 24_pan_fasta_clust_cluster.tsv | awk '$1==$3' |
    perl -pe 's/(pan\d+)\t/$1__/g' > 24_pan_fasta_cluster_pruned.tsv

  echo "   Join chromosome numbers to gene pairs and count chromosome correspondences among the pairs."
  join <(perl -pe 's/pan\d+__//g' 24_pan_fasta_cluster_pruned.tsv | sort -k1,1) \
       22_syn_pan_aug_extra_pctl"${pctl_low}"_posn.gene_chr.hsh | perl -pe 's/ /\t/g' | sort -k2,2 | 
     join -1 2 -2 1 - 22_syn_pan_aug_extra_pctl"${pctl_low}"_posn.gene_chr.hsh | 
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
  someseq=$(head "${WORK_DIR}"/07_pan_fasta_cds.fna | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
  SEQTYPE=$(check_seq_type "${someseq}")

  full_out_dir="${out_dir}"
  stats_file=${full_out_dir}/stats.txt
  export ANN_REX=${annot_str_regex}

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p "$full_out_dir"
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
    if [ -f "${WORK_DIR}"/"$file" ]; then
      cp "${WORK_DIR}"/"$file" "${full_out_dir}"/
    else 
      echo "Warning: couldn't find file ${WORK_DIR}/$file; skipping"
    fi
  done

  echo "Copy manifest file into the output directory"
  if [ -f "${submit_dir}/manifests/MANIFEST_output_pan.yml" ]; then
    cp "${submit_dir}/manifests/MANIFEST_output_pan.yml" "$full_out_dir"/
  else 
    echo "Couldn't find file manifests/MANIFEST_output_pan.yml"
  fi

  echo "Run of program $scriptname, version $version" > "${stats_file}"

  end_time=$(date)
  cat "${WORK_DIR}"/stats/tmp.timing >> "${stats_file}"
  printf "Run ended at:   %s\n\n" "$end_time" >> "${stats_file}"

  # Report sequence type 
  if [[ "$SEQTYPE" == 3 ]]; then
    echo "Sequence type: nucleotide" >> "${stats_file}"
  else
    echo "Sequence type: protein" >> "${stats_file}"
  fi

  echo "  Report parameters from config file"
  printf "Parameter  \tvalue\n" >> "${stats_file}"
  for key in ${pandagma_conf_params}; do
    printf '%-15s\t%s\n' "${key}" "${!key}" >> "${stats_file}"
  done

  printf "\nOutput directory for this run:\t%s\n" "${full_out_dir}" >> "${stats_file}"

  echo "  Report threshold for inclusion in \"pctl${pctl_low}\""
  max_annot_ct=$(awk '$1!~/^#/ {print $2}' "${full_out_dir}"/18_syn_pan_aug_extra.counts.tsv | sort -n | uniq | tail -1)
  pctl_low_threshold=$(awk -v MCP="$pctl_low" -v MAC="$max_annot_ct" 'BEGIN{print MCP*MAC/100}')

  echo "  Report orthogroup composition statistics for the three main cluster-calculation steps"

  printf '\n  %-20s\t%s\n' "Statistic" "value" >> "${stats_file}"
  printf "The global mode may be for a smaller OG size. Modes below are greater than the specified core threshold.\n" \

  clustcount=1
  # When the number of main annotation sets is small (<4), the core-threshold-ceiling may be lower than 
  # the largest number of clusters. In that case, set $CTceil=2 to the smallest cluster size ($CTceil=2)
  for clustering in 06_syn_pan 12_syn_pan_aug 18_syn_pan_aug_extra; do
    clustfile=${full_out_dir}/$clustering.clust.tsv
    if [[ -f $clustfile ]]; then
      if [[ $clustcount == 1 ]]; then
        printf "\n== Initial clusters (containing only genes within synteny blocks)\n" >> "${stats_file}"
        (( largest=$(awk '{print NF-1}' "$clustfile" | sort -n | tail -1) ))
        (( mode=$(awk "{print NF-1}" "$clustfile" |
          sort -n | uniq -c | awk '{print $1 "\t" $2}' | sort -n | tail -1 | awk '{print $2}') ))
      elif [[ $clustcount == 2 ]]; then
        printf "\n== Augmented clusters (unanchored sequences added to the initial clusters)\n" >> "${stats_file}"
        (( largest=$(awk '{print NF-1}' "$clustfile" | sort -n | tail -1) ))
        (( mode=$(awk "{print NF-1}" "$clustfile" |
          sort -n | uniq -c | awk '{print $1 "\t" $2}' | sort -n | tail -1 | awk '{print $2}') ))
      elif [[ $clustcount == 3 ]]; then
        printf "\n== Augmented-extra clusters (with sequences from extra annotation sets)\n" >> "${stats_file}"
        (( largest=$(awk '{print NF-1}' "$clustfile" | sort -n | tail -1) ))
        CTceil=$(echo "$pctl_low_threshold" | awk '{print int($1+1)}')
        if (( CTceil>largest )); then (( CTceil=2 )); fi; export CTceil
        (( mode=$(awk -v CT="$CTceil" "(NF-1)>=CT {print NF-1}" "$clustfile" |
          sort -n | uniq -c | awk '{print $1 "\t" $2}' | sort -n | tail -1 | awk '{print $2}') ))
        printf "    The pctl${pctl_low} set consists of orthogroups with at least %.0f genes per OG (>= %d * %d/100 sets).\n" \
           "$CTceil" "$max_annot_ct" "$pctl_low" >> "${stats_file}"
      fi

      export mode
      (( clusters=$(wc -l < "$clustfile") ))
      (( num_at_mode=$(awk -v MODE="$mode" '(NF-1)==MODE {ct++} END{print ct}' "$clustfile") ))
      (( seqs_clustered=$(awk '{sum+=NF-1} END{print sum}' "$clustfile") ))

      {
        printf "  %-20s\t%s\n" "Cluster file" "$clustering.clust.tsv"
        printf "  %-20s\t%s\n" "num_of_clusters" $clusters
        printf "  %-20s\t%s\n" "largest_cluster" "$largest" 
      } >> "${stats_file}"
      if (( clustcount<3 )); then
        printf "  %-20s\t%d\n" "modal_clst_size" "$mode" >> "${stats_file}"
        printf "  %-20s\t%d\n" "num_at_mode" $num_at_mode >> "${stats_file}"
      else
        printf "  %-20s\t%d\n" "modal_clst_size>=$CTceil" "$mode" >> "${stats_file}"
        printf "  %-20s\t%d\n" "num_at_mode>=$CTceil" $num_at_mode >> "${stats_file}"
      fi
      printf "  %-20s\t%d\n" "seqs_clustered" $seqs_clustered >> "${stats_file}"
  
      (( clustcount=$((clustcount+1)) ))
    else
      printf "File %f is not available; skipping\n" "$clustfile"
      printf "File %f is not available; skipping\n" "$clustfile" >> "${stats_file}"
    fi
  done

  echo "  Print sequence composition statistics for each annotation set"
  printf "\n== Sequence stats for CDS files\n" >> "${stats_file}"
  printf "  Class:  seqs     min max    N50    ave     annotation_name\n" >> "${stats_file}" 
  if [ -f "${WORK_DIR}"/stats/tmp.fasta_seqstats ]; then
    {
      cat "${WORK_DIR}"/stats/tmp.fasta_seqstats
      printf "\n  Avg:   "
    } >> "${stats_file}" 
    transpose.pl "${WORK_DIR}"/stats/tmp.fasta_seqstats |
      perl -ane 'BEGIN{use List::Util qw(sum)}; 
                 if ($F[0]=~/^\d+/){
                   $sum=sum @F; $ct=scalar(@F); $avg=$sum/$ct;
                   printf " %4d ", $avg;
                 END{print "   all_annot_sets\n"}
                 }' >> "${stats_file}"
  fi

  printf "\n== Sequence stats for final pangene CDS files -- pctl%f and trimmed\n" "$pctl_low" >> "${stats_file}"
  printf "  Class:   seqs     min max    N50    ave     annotation_name\n" >> "${stats_file}" 
  annot_name=23_syn_pan_pctl${pctl_low}_posn_cds.fna
    printf "  pctl%f: " "${pctl_low}" >> "${stats_file}"
    cat_or_zcat "${WORK_DIR}/23_syn_pan_pctl${pctl_low}_posn_cds.fna" | calc_seq_stats >> "${stats_file}"

  echo "  Print per-annotation-set coverage stats (sequence counts, sequences retained)"
  #   tmp.gene_count_start was generated during run_ingest
  printf "\n== Proportion of initial genes retained in the \"aug_extra\" and \"pctl%f\" sets:\n" "${pctl_low}" \
    >> "${stats_file}"

  cut -f2 "${WORK_DIR}"/18_syn_pan_aug_extra.hsh.tsv | 
    perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' |
    sort | uniq -c | awk '{print $2 "\t" $1}' > "${WORK_DIR}"/stats/tmp.gene_count_all_end

  cut -f2 "${WORK_DIR}"/22_syn_pan_aug_extra_pctl"${pctl_low}".hsh.tsv | 
    perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' |
    sort | uniq -c | awk '{print $2 "\t" $1}' > "${WORK_DIR}"/stats/tmp.gene_count_pctl"${pctl_low}"_end

  paste "${WORK_DIR}"/stats/tmp.gene_count_start \
        "${WORK_DIR}"/stats/tmp.gene_count_all_end \
        "${WORK_DIR}"/stats/tmp.gene_count_pctl"${pctl_low}"_end | 
    awk 'BEGIN{print "  Start\tEnd_all\tEnd_core\tPct_kept_all\tPct_kept_core\tAnnotation_name"} 
        { printf "  %i\t%i\t%i\t%2.1f\t%2.1f\t%s\n", $2, $4, $6, 100*($4/$2), 100*($6/$2), $1 }'  >> "${stats_file}"

  echo "  Print counts per accession"
  if [ -f "${full_out_dir}"/18_syn_pan_aug_extra.counts.tsv ]; then
    {
    printf "\n== For all annotation sets, counts of genes-in-orthogroups and counts of orthogroups-with-genes:\n" 
    printf "  gns-in-OGs  OGs-w-gns  OGs-w-gns/gns  pct-non-null-OGs  pct-null-OGs  annot-set\n" 
    transpose.pl "${full_out_dir}"/18_syn_pan_aug_extra.counts.tsv | 
      perl -lane 'next if ($.<=3); 
        $ct=0; $sum=0; $nulls=0; $OGs=0;
        for $i (@F[1..(@F-1)]){
          $OGs++;
          if ($i>0){$ct++; $sum+=$i}
          if ($i==0){$nulls++}
        }; 
        printf("  %d\t%d\t%.2f\t%.2f\t%.2f\t%s\n", $sum, $ct, 100*$ct/$sum, 100*($OGs-$nulls)/$OGs, 100*$nulls/$OGs, $F[0])' \
        
    } >> "${stats_file}"
  fi

  echo "  Print histograms"
  if [ -f "${full_out_dir}"/06_syn_pan.clust.tsv ]; then
    printf "\nCounts of initial clusters by cluster size, file 06_syn_pan.clust.tsv:\n" >> "${stats_file}"
    awk '{print NF-1}' "${full_out_dir}"/06_syn_pan.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> "${stats_file}"
  fi

  if [ -f "${full_out_dir}"/12_syn_pan_aug.clust.tsv ]; then
    printf "\nCounts of augmented clusters by cluster size, file 12_syn_pan_aug.clust.tsv:\n" >> "${stats_file}"
    awk '{print NF-1}' "${full_out_dir}"/12_syn_pan_aug.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> "${stats_file}"
  fi

  if [ -f "${full_out_dir}"/18_syn_pan_aug_extra.clust.tsv ]; then
    printf "\nCounts of augmented-extra clusters by cluster size, file 18_syn_pan_aug_extra.clust.tsv:\n" >> "${stats_file}"
    awk '{print NF-1}' "${full_out_dir}"/18_syn_pan_aug_extra.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> "${stats_file}"
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
    if [ -d "$dir" ]; then echo "  Removing directory $dir"; rm -rf "$dir" &
    fi
  done
  #for file in 10* 11* 14* 20* 21* 23* 24* consen*; do
  for file in 10* 11* 14* 20*; do
    if [ -f "$file" ]; then echo "  Removing file $file"; rm "$file"; 
    fi
  done
  wait
  cd "$OLDPWD"
}

########################################
# Main program

pandagma_conf_params='clust_iden clust_cov extra_iden mcl_inflation 
  strict_synt pctl_low pctl_med pctl_hi 
  consen_prefix annot_str_regex order_method preferred_annot'

# The steps align_cds, align_protein, model_and_trim, calc_trees, and xfr_aligns_trees may be run separately.
commandlist="ingest mmseqs filter dagchainer mcl consense cluster_rest add_extra pick_exemplars \
             filter_to_pctile tabularize order_and_name  calc_chr_pairs summarize"

dependencies='mmseqs dagchainer mcl cons famsa hmmalign hmmbuild run_DAG_chainer.pl'

declare clust_iden clust_cov extra_iden mcl_inflation strict_synt \
        pctl_low pctl_med pctl_hi consen_prefix annot_str_regex preferred_annot order_method

main_pan_fam "$@"

#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates pan-gene 
# clusters using the programs mmseqs, dagchainer, and mcl. 
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2022
#
scriptname=`basename "$0"`
version="2022-12-25"
set -o errexit -o errtrace -o nounset -o pipefail

export NPROC=${NPROC:-1}
export MMSEQS_NUM_THREADS=${NPROC} # mmseqs otherwise uses all cores by default
export CONF=${CONF:-${PWD}/pandagma.conf}
export WORK_DIR=${WORK_DIR:-${PWD}/work}

# mmseqs uses a significant number of threads on its own. Set a maximum, which may be below NPROC.
MMSEQSTHREADS=$(( 10 < ${NPROC} ? 10 : ${NPROC} ))

pandagma_conf_params='clust_iden clust_cov consen_iden extra_iden mcl_inflation min_core_prop
                      dagchainer_args out_dir_base consen_prefix annot_str_regex preferred_annot'

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

TOP_DOC="""Compute pan-gene clusters using the programs mmseqs, dagchainer, and mcl, and
additional pre- and post-refinement steps.

Usage:
        pandagma.sh SUBCOMMAND [SUBCOMMAND_OPTIONS]

Primary coding sequence (fasta) and annotation (GFF3 or BED) files must be listed in the
fasta_files and annotation_files variables defined in pandagma.conf, which by default must exist
within the working directory from where this script is called.

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

At the end of the process, remaining genes will be added to initial clusters, based on homology.
Remaining genes may be those falling on unanchored scaffolds, or on chromosomes by not part of
synteny blocks and so not making it into the synteny-based clusters.

Subommands (in order they are usually run):
            version - Report installed package version
         run ingest - Prepare the assembly and annotation files for analysis
         run mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies
         run filter - Filter the synteny results for chromosome pairings, returning gene pairs.
     run dagchainer - Run DAGchainer to filter for syntenic blocks
            run mcl - Derive clusters, with Markov clustering
       run consense - Calculate a consensus sequences from each pan-gene set, 
                      adding sequences missed in the first clustering round.
      run add_extra - Add other gene model sets to the primary clusters. Useful for adding
                      annotation sets that may be of lower or uncertain quality.
  run name_pangenes - Assign pan-gene names with consensus chromosomes and ordinal positions.
 run filter_to_core - Calculate orthogroup composition and filter fasta files to core orthogroups.
 run calc_chr_pairs - Report observed chromosome pairs; useful for preparing expected_chr_matches.tsv
      run summarize - Move results into output directory, and report summary statistics.

Variables in pandagma config file (Set the config with the CONF environment variable)
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.60]
        consen_iden - Minimum identity threshold for consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of \"extra\" annotations [90]
      min_core_prop - Minimum fraction of annotation sets for an orthogroup to be \"core\" [0.333333333333333]
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]
      consen_prefix - Prefix to use in names for genomic ordered consensus IDs [Genus.pan1]
       out_dir_base - Base name for the output directory [default: './out']
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g. 
                        \"([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+\" 
                          for four dot-separated fields, e.g. vigan.Shumari.gnm1.ann1
                        or \"(\D+\d+\D+)\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    preferred_annot - String to match and select an annotation set, from a gene ID.
                        This is used for picking representative IDs+sequence from an orthogroup, when
                        this annotation is among those with the median length for the orthogroup.
                        Otherwise, one is selected at random from those with median length.

Environment variables:
               CONF - Path of the pandagma config file, default: \"${CONF}\"
           WORK_DIR - Location of working files, default: \"${WORK_DIR}\"
              NPROC - Number of processors to use (default 1)
"""
#
# Helper functions begin here
#
canonicalize_paths() {
  fasta_files=($(realpath --canonicalize-existing "${fasta_files[@]}"))
  annotation_files=($(realpath --canonicalize-existing "${annotation_files[@]}"))
  protein_files=($(realpath --canonicalize-existing "${protein_files[@]}"))
  if (( ${#fasta_files_extra[@]} > 0 ))
  then
    fasta_files_extra=($(realpath --canonicalize-existing "${fasta_files_extra[@]}"))
    annotation_files_extra=($(realpath --canonicalize-existing "${annotation_files_extra[@]}"))
  fi
  readonly chr_match_list=${expected_chr_matches:+$(realpath "${expected_chr_matches}")}
  readonly submit_dir=${PWD}

  fasta_file=$(basename "${fasta_files[0]}" .gz)
  fna="${fasta_file##*.}"
}
cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}
check_seq_type () {
  someseq=${1}
  export proportion_nuc=$(echo $someseq | fold -w1 | 
                          awk '$1~/[ATCGN]/ {nuc++} $1!~/[ATCGN]/ {not++} END{print nuc/(nuc+not)}')
  perl -le '$PN=$ENV{"proportion_nuc"}; if ($PN>0.9){print 3} else {print 1}'
}

### run functions ###

run_ingest() {
# Add positional information from GFF3 or 4- or 6-column BED to FASTA IDs
# BED start coordinate converted to 1-based
  cd "${WORK_DIR}"
  echo; echo "Run ingest: from fasta and gff or bed data, create fasta with IDs containing positional info."
  
  mkdir -p 02_fasta_nuc 02_fasta_prot 01_posn_hsh stats

    # Prepare the tmp.gene_count_start to be joined, in run_summarize, with tmp.gene_count_end.
    # This is captured from the gene IDs using the annot_str_regex set in the config file.
    cat /dev/null > stats/tmp.gene_count_start
    cat /dev/null > stats/tmp.fasta_list
    start_time=`date`
    printf "Run started at: $start_time\n" > stats/tmp.timing

  export ANN_REX=${annot_str_regex}

  for (( file_num = 0; file_num < ${#fasta_files[@]} ; file_num++ )); do
    file_base=$(basename ${fasta_files[file_num]%.*})
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | 
      gff_or_bed_to_hash.awk > 01_posn_hsh/$file_base.hsh
    hash_into_fasta_id.pl -nodef -fasta "${fasta_files[file_num]}" \
                          -hash 01_posn_hsh/$file_base.hsh \
                          -out 02_fasta_nuc/$file_base
    annot_name=$(head -1 02_fasta_nuc/$file_base | perl -pe 's/>.+__(.+)__\d+__\d+$/$1/' | 
      perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
    cat_or_zcat "${fasta_files[file_num]}" | 
      awk -v ANNOT=$annot_name '$1~/^>/ {ct++} END{print ANNOT "\t" ct}' >> stats/tmp.gene_count_start
    echo "Main:   $file_base" >> stats/tmp.fasta_list
  done

  # Also get position information from the "extra" annotation sets, if any.
  if (( ${#fasta_files_extra[@]} > 0 ))
  then
    for (( file_num = 0; file_num < ${#fasta_files_extra[@]} ; file_num++ )); do
      file_base=$(basename ${fasta_files_extra[file_num]%.*})
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra[file_num]}" | 
        gff_or_bed_to_hash.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${fasta_files_extra[file_num]}" \
                            -hash 01_posn_hsh/$file_base.hsh \
                            -out 02_fasta_nuc/$file_base
      annot_name=$(head -1 02_fasta_nuc/$file_base | perl -pe 's/>.+__(.+)__\d+__\d+$/$1/' | 
         perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
      cat_or_zcat "${fasta_files_extra[file_num]}" | 
        awk -v ANNOT=$annot_name '$1~/^>/ {ct++} END{print ANNOT "\t" ct}' >> stats/tmp.gene_count_start
      echo "Extra:  $file_base" >> stats/tmp.fasta_list
    done
  fi

  # Also protein files, if available
  if (( ${#protein_files[@]} > 0 ))
  then
    for (( file_num = 0; file_num < ${#protein_files[@]} ; file_num++ )); do
      echo "  Copying protein file ${protein_files[file_num]}"
      cp "${protein_files[file_num]}" 02_fasta_prot/
    done
  fi

  sort -o stats/tmp.gene_count_start stats/tmp.gene_count_start
}

run_mmseqs() {
  # Do mmseqs clustering on all pairings of the main annotation sets (not the extra ones though)
  cd "${WORK_DIR}"
  echo; echo "Run mmseqs -- at ${clust_iden} percent identity and minimum of ${clust_cov}% coverage."
  #

  mkdir -p 03_mmseqs 03_mmseqs_tmp
  for (( file1_num = 0; file1_num < ${#fasta_files[@]} ; file1_num++ )); do
    qry_base=$(basename ${fasta_files[file1_num]%.*} .$fna)
    for (( file2_num = $( expr $file1_num + 1 ); file2_num < ${#fasta_files[@]} ; file2_num++ )); do
      sbj_base=$(basename ${fasta_files[file2_num]%.*} .$fna)
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

run_filter() {
  echo; echo "From mmseqs cluster output, split out the following fields: molecule, gene, start, stop."
  cd "${WORK_DIR}"
  mkdir -p 04_dag
  if [[ -f ${chr_match_list} ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns from file ${chr_match_list}"
    for mmseqs_path in 03_mmseqs/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      echo "  $outfilebase"
      cat ${mmseqs_path} | 
        filter_mmseqs_by_chroms.pl -chr_pat ${chr_match_list} > 04_dag/${outfilebase}_matches.tsv &

      # allow to execute up to $NPROC in parallel
      if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
    done
    wait # wait for last jobs to finish
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches.tsv file was provided, so proceeding without chromosome-pair filtering."
    for mmseqs_path in 03_mmseqs/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      perl -pe 's/__/\t/g' > 04_dag/${outfilebase}_matches.tsv < "${mmseqs_path}"
    done
  fi
}

run_dagchainer() {
  # Identify syntenic blocks, using DAGchainer
  cd "${WORK_DIR}"
  echo; echo "Run DAGchainer, using args \"${dagchainer_args}\""
  for match_path in 04_dag/*_matches.tsv; do
    align_file=`basename $match_path _matches.tsv`
    echo "Running DAGchainer on comparison: $align_file"
    echo "  run_DAG_chainer.pl $dagchainer_args -i $match_path"; echo
    # run_DAG_chainer.pl writes temp files to cwd;
    # use per-process temp directory to avoid any data race
    (
      tmpdir=$(mktemp -d)
      cd "${tmpdir}"
      run_DAG_chainer.pl $dagchainer_args -i "${OLDPWD}/${match_path}" 1>/dev/null
      rmdir ${tmpdir}
    ) &
    # allow to execute up to $NPROC in parallel
    [ $(jobs -r -p | wc -l) -ge ${NPROC} ] && wait -n
  done
  wait # wait for last jobs to finish

  awk '$1!~/^#/ {print $2 "\t" $6}' 04_dag/*_matches.tsv > 05_homology_pairs.tsv
  awk '$1!~/^#/ {print $2 "\t" $6}' 04_dag/*.aligncoords > 05_synteny_pairs.tsv
}

run_mcl() {
  # Calculate clusters using Markov clustering
  cd "${WORK_DIR}"
  printf "\nDo Markov clustering with inflation parameter $mcl_inflation and ${NPROC} threads\n"
  echo "MCL COMMAND: mcl 05_synteny_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv"
  mcl 05_synteny_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv \
    1>/dev/null
 
  # Add cluster IDs
  awk '{padnum=sprintf("%05d", NR); print "pan" padnum "\t" $0}' tmp.syn_pan.clust.tsv > 06_syn_pan.clust.tsv
  rm tmp.syn_pan.clust.tsv

  # Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 06_syn_pan.clust.tsv > 06_syn_pan.hsh.tsv
}

run_consense() {
  echo; 
  echo "Add previously unclustered sequences into an \"augmented\" pan-gene set, by homology."
  cd "${WORK_DIR}"
  mkdir -p 07_pan_fasta lists

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  echo "    Fasta file:" "${fasta_files[@]}"
  get_fasta_from_family_file.pl "${fasta_files[@]}" -fam 06_syn_pan.clust.tsv -out 07_pan_fasta

  cat /dev/null > 07_pan_fasta.$fna
  for path in 07_pan_fasta/*; do
    pan_file=`basename $path`
    cat $path | awk -v panID=$pan_file '
                      $1~/^>/ {print ">" panID "__" substr($0,2) }
                      $1!~/^>/ {print $1}
                    ' >> 07_pan_fasta.$fna
  done

  echo "  Pick a representative seq. for each orthogroup - as a sequence with the median length for that OG."
  echo "  cat 07_pan_fasta.$fna | 
    pick_family_rep.pl -prefer $preferred_annot -out 08_pan_fasta_clust_rep_seq.$fna"
  cat 07_pan_fasta.$fna | pick_family_rep.pl -prefer $preferred_annot -out 08_pan_fasta_clust_rep_seq.$fna

  echo "  Get sorted list of all genes, from the original fasta files"
  cat_or_zcat "${fasta_files[@]}" | awk '/^>/ {print substr($1,2)}' | sort > lists/09_all_genes

  echo "  Get sorted list of all clustered genes"
  awk '$1~/^>/ {print $1}' 07_pan_fasta/* | sed 's/>//' | sort > lists/09_all_clustered_genes

  echo "  Get list of genes not in clusters"
  comm -13 lists/09_all_clustered_genes lists/09_all_genes > lists/09_genes_not_in_clusters

  echo "  Retrieve the non-clustered genes"
  cat_or_zcat "${fasta_files[@]}" |
    get_fasta_subset.pl -in /dev/stdin -clobber -lis lists/09_genes_not_in_clusters \
                        -out 09_genes_not_in_clusters.$fna 

  echo "  Search non-clustered genes against genes already clustered."

  # Check sequence type (in case this run function is called separately from the usually-prior ones)
  someseq=$(head 07_pan_fasta.$fna | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
  SEQTYPE=$(check_seq_type "${someseq}") # 3=nuc; 1=pep
  echo "SEQTYPE is: $SEQTYPE"

  mmseqs easy-search 09_genes_not_in_clusters.$fna \
                     07_pan_fasta.$fna \
                     10_unclust.x.07_pan_fasta.m8 \
                     03_mmseqs_tmp \
                     --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null 

  echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
  echo "  Use the \"main set\" $clust_iden threshold."
  top_line.awk 10_unclust.x.07_pan_fasta.m8 | 
    awk -v IDEN=${clust_iden} '$3>=IDEN {print $2 "\t" $1}' | perl -pe 's/^(pan\d+)__\S+/$1/' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk >  11_syn_pan_leftovers.clust.tsv

  echo "  Retrieve sequences for the leftover genes"
  mkdir -p 11_pan_leftovers
  get_fasta_from_family_file.pl "${fasta_files[@]}" \
    -fam 11_syn_pan_leftovers.clust.tsv -out 11_pan_leftovers/

  echo "  Make augmented cluster sets"
  cat /dev/null > 12_syn_pan_aug.clust.tsv
  augment_cluster_sets.awk leftovers_dir=11_pan_leftovers 07_pan_fasta/* > 12_syn_pan_aug.clust.tsv
}

run_add_extra() {
  if (( ${#fasta_files_extra[@]} > 0 ))
  then # handle the "extra" annotation files
    echo; echo "== Add extra annotation sets to the augmented clusters, by homology =="
    echo "  Search non-clustered genes against pan-gene consensus sequences"
    cd "${WORK_DIR}"
    mkdir -p 13_extra_out_dir 13_pan_aug_fasta

    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    get_fasta_from_family_file.pl "${fasta_files[@]}" -fam 12_syn_pan_aug.clust.tsv -out 13_pan_aug_fasta
    
    cat /dev/null > 13_pan_aug_fasta.$fna
    for path in 13_pan_aug_fasta/*; do
      pan_file=`basename $path`
      cat $path | awk -v panID=$pan_file '
                        $1~/^>/ {print ">" panID "__" substr($0,2) }
                        $1!~/^>/ {print $1}
                      ' >> 13_pan_aug_fasta.$fna
    done

    # Check sequence type (in case this run function is called separately from the usually-prior ones)
    someseq=$(head 07_pan_fasta.$fna | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
    SEQTYPE=$(check_seq_type "${someseq}") # 3=nuc; 1=pep
    echo "SEQTYPE is: $SEQTYPE"

    for path in "${fasta_files_extra[@]}"; do
      fasta_file=`basename ${path%.*}`
      echo "Extra: $fasta_file"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      mmseqs easy-search "${path}" \
                         13_pan_aug_fasta.$fna \
                         13_extra_out_dir/${fasta_file}.x.all_cons.m8 \
                         $MMTEMP \
                         --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null & # background

      if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    wait # wait for jobs to finish
  
    echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
    echo "  Use identity threshold extr_iden: $extra_iden."
    top_line.awk 13_extra_out_dir/*.x.all_cons.m8 |
      awk -v IDEN=${extra_iden} '$3>=IDEN {print $2 "\t" $1}' | perl -pe 's/^(pan\d+)__\S+/$1/' |
      perl -pe 's/^(\w+)\.\d+/$1/' |
      sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > 14_syn_pan_extra.clust.tsv
  
    echo "  Retrieve sequences for the extra genes"
    mkdir -p 16_pan_leftovers_extra
    get_fasta_from_family_file.pl "${fasta_files_extra[@]}" \
       -fam 14_syn_pan_extra.clust.tsv -out 16_pan_leftovers_extra/
  
    echo "  Make augmented cluster sets"
    augment_cluster_sets.awk leftovers_dir=16_pan_leftovers_extra 13_pan_aug_fasta/* |
      cat > 17_syn_pan_aug_extra.clust.tsv

    echo "  Merge fasta sets"
    mkdir -p 19_pan_aug_leftover_merged
    for path in 13_pan_aug_fasta/*; do
      file=`basename $path`
      if [[ -f "16_pan_leftovers_extra/$file" ]]; then
        cat $path 16_pan_leftovers_extra/$file > 19_pan_aug_leftover_merged/$file
      else
        cp $path 19_pan_aug_leftover_merged/
      fi
    done
  
    cat /dev/null > 20_pan_fasta.$fna
    for path in 19_pan_aug_leftover_merged/*; do
      pan_file=`basename $path`
      cat $path | awk -v panID=$pan_file '
                        $1~/^>/ {print ">" panID "__" substr($0,2) }
                        $1!~/^>/ {print $1}
                      ' >> 20_pan_fasta.$fna
    done

    echo "  Pick a representative sequence for each pangene set - as a sequence with the median length for that set."
    echo "  cat 20_pan_fasta.$fna | pick_family_rep.pl -prefer $preferred_annot -out 21_pan_fasta_clust_rep_seq.$fna"
    cat 20_pan_fasta.$fna | pick_family_rep.pl -prefer $preferred_annot -out 21_pan_fasta_clust_rep_seq.$fna
    
    perl -pi -e 's/__/  /' 21_pan_fasta_clust_rep_seq.$fna
  
  else  
    echo "== No annotations were designated as \"extra\", so just promote the syn_pan_aug files as syn_pan_aug_extra. ==" 
    cp 07_pan_fasta.$fna 20_pan_fasta.$fna
    cp 12_syn_pan_aug.clust.tsv 17_syn_pan_aug_extra.clust.tsv
    cp 08_pan_fasta_clust_rep_seq.$fna 21_pan_fasta_clust_rep_seq.$fna
    perl -pi -e 's/__/  /' 21_pan_fasta_clust_rep_seq.$fna
  fi
}

run_filter_to_core() {
  cd "${WORK_DIR}"

  echo "  Calculate matrix of gene counts per orthogroup and annotation set"
  calc_pan_stats.pl -pan 17_syn_pan_aug_extra.clust.tsv -out 17_syn_pan_aug_extra.counts.tsv
  max_annot_ct=$(cat 17_syn_pan_aug_extra.counts.tsv | 
                       awk '$1!~/^#/ {print $2}' | sort -n | uniq | tail -1)

  echo "  Select orthogroups with genes from at least min_core_prop*max_annot_ct annotation sets"
  cat 17_syn_pan_aug_extra.counts.tsv | 
    awk -v MINCORE=$min_core_prop -v ANNCT=$max_annot_ct '$2>=ANNCT*MINCORE && $1!~/^#/ {print $1}' |
      cat > lists/lis.17_syn_pan_aug_extra.core

  echo "  Get a fasta subset with only genes from at least min_core_prop*max_annot_ct annotation sets"
  get_fasta_subset.pl -in 21_pan_fasta_clust_rep_seq.$fna -list lists/lis.17_syn_pan_aug_extra.core \
                      -clobber -out 22_pan_fasta_rep_seq_core.$fna


  echo "  Get a clust.tsv file with orthogroups with at least min_core_prop*max_annot_ct annotation sets"
  join <(LC_ALL=C sort -k1,1 lists/lis.17_syn_pan_aug_extra.core) \
       <(LC_ALL=C sort -k1,1 17_syn_pan_aug_extra.clust.tsv) |
          cat > 22_syn_pan_aug_extra_core.clust.tsv
}

run_name_pangenes() {
  cd "${WORK_DIR}"

  echo "  Reshape from mcl output format, clustered IDs on one line, to a hash format"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 22_syn_pan_aug_extra_core.clust.tsv |
    cat > 22_syn_pan_aug_extra_core.hsh.tsv

  echo "  Add positional information to the hash output."
  join -a 1 -1 2 -2 1 <(sort -k2,2 22_syn_pan_aug_extra_core.hsh.tsv) <(cat 01_posn_hsh/*hsh | sort -k1,1) | 
    perl -pe 's/__/\t/g; s/ /\t/g' | awk -v OFS="\t" '{print $2, $1, $3, $5, $6}' |
    sort -k1,1 -k2,2 > 22_syn_pan_aug_extra_core_posn.hsh.tsv

  echo "  Calculate consensus pan-gene positions"
  cat 22_syn_pan_aug_extra_core_posn.hsh.tsv | consen_pangene_posn.pl -pre ${consen_prefix}.chr |
    sort -k2,2 -k3n,3n | name_ordered_genes.awk > consen_${consen_prefix}.tsv

  echo "  Reshape defline into a hash, e.g. pan47789	Glycine.pan3.chr01__Glycine.pan3.chr01_000100__45224__45786"
  cat consen_${consen_prefix}.tsv | 
    perl -pe 's/^(\S+)\t([^.]+\.[^.]+)\.(chr\d+)_(\d+)\t(\d+)\t(\d+)/$1\t$1__$2.$3_$4__$5__$6/' > consen_posn.hsh

  echo "  Hash position information into fasta file"
  hash_into_fasta_id.pl -fasta 22_pan_fasta_rep_seq_core.$fna -hash consen_posn.hsh |
    grep -v "HASH UNDEFINED" > 22_pan_fasta_rep_seq_core_posnTMP.$fna

  echo "  Reshape defline, and sort by position"
  fasta_to_table.awk 22_pan_fasta_rep_seq_core_posnTMP.$fna | sed 's/__/\t/g; s/ /\t/' | 
    perl -pe 's/^(\w+)\s+(.+)\.(chr\d+)_(\d+)\s+/$1\t$2\t$3\t$4\t/' | sed 's/chr//' | sort -k3n -k4n |
    awk '{print ">" $2 ".chr" $3 "_" $4 " " $1 " " $5 " " $6 " " $7; print $8}' |
      cat > 23_syn_pan_cent_merged_core_posn.$fna
 
  echo "  Re-cluster, to identify neighboring genes that are highly similar"
    mmseqs easy-cluster 23_syn_pan_cent_merged_core_posn.$fna 24_pan_fasta 03_mmseqs_tmp \
    --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null

  echo "  Parse mmseqs clusters into pairs of genes that are similar and ordinally close"
  close=10 # genes within a neighborhood of +-close genes among the consensus-orderd genes
  cat 24_pan_fasta_cluster.tsv | perl -pe 's/\S+\.chr(\d+)_(\d+)/$1\t$2/g' | 
    awk -v CLOSE=$close -v PRE=$consen_prefix '$1==$3 && sqrt(($2/100-$4/100)^2)<=CLOSE \
             {print PRE ".chr" $1 "_" $2 "\t" PRE ".chr" $3 "_" $4}' > 24_pan_fasta_cluster.closepairs

  echo "  Cluster the potential near-duplicate neighboring paralogs"
  mcl  24_pan_fasta_cluster.closepairs --abc -o 24_pan_fasta_cluster.close.clst
 
  echo "  Keep the first gene from the cluster and discard the rest."
  cat 24_pan_fasta_cluster.close.clst | awk '{$1=""}1' | awk '{$1=$1}1' | tr ' ' '\n' | 
    sort -u > lists/lis.clust_genes_remove
  get_fasta_subset.pl -in 23_syn_pan_cent_merged_core_posn.$fna -xclude -clobber \
    -lis lists/lis.clust_genes_remove -out 24_syn_pan_cent_merged_core_posn_trimd.$fna 

  echo "  Also get all corresponding protein sequences for genes in lists/lis.17_syn_pan_aug_extra.core"
  faa="faa"
  cat /dev/null > 20_pan_fasta_prot.$faa
  for filepath in 02_fasta_prot/*.gz; do 
    zcat $filepath >> 20_pan_fasta_prot.$faa
  done
  cat 23_syn_pan_cent_merged_core_posn.$fna | awk '$1~/^>/ {print $5}' > lists/lis.23_syn_pan_cent_merged_core_posn
  get_fasta_subset.pl -in 20_pan_fasta_prot.$faa -list lists/lis.23_syn_pan_cent_merged_core_posn \
                      -clobber -out 23_syn_pan_cent_merged_core_posn_protTMP.$faa

  echo "  Hash pan-ID into protein fasta file"
  cat 23_syn_pan_cent_merged_core_posn.$fna | 
    awk '$1~/^>/ {print $5 "\t" substr($1, 2) "__" $2 "__" $3 "__" $4}' > lists/23_syn_pan_cent_merged_core_posn.hsh
  hash_into_fasta_id.pl -hash lists/23_syn_pan_cent_merged_core_posn.hsh -swap_IDs -nodef \
    -fasta 23_syn_pan_cent_merged_core_posn_protTMP.$faa |
    perl -pe 's/__/ /g' > 23_syn_pan_cent_merged_core_posn_prot.$faa
}

run_calc_chr_pairs() {
  cd "${WORK_DIR}"
  echo "Generate a report of observed chromosome pairs"

  echo "  Identify gene pairs, ussing mmseqs --easy_cluster"
    mmseqs easy-cluster 20_pan_fasta.$fna 24_pan_fasta_clust 03_mmseqs_tmp \
    --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null

  echo "   Extract chromosome-chromosome correspondences"
  cut -f2,3 22_syn_pan_aug_extra_core_posn.hsh.tsv | sort -k1,1 > 22_syn_pan_aug_extra_core_posn.gene_chr.hsh

  echo "   From mmseqs cluster table, prune gene pairs, keeping only those in the same pangene cluster"
  cat 24_pan_fasta_clust_cluster.tsv | perl -pe 's/__/\t/g' | awk '$1==$3' |
    perl -pe 's/(pan\d+)\t/$1__/g' > 24_pan_fasta_cluster_pruned.tsv

  echo "   Join chromosome numbers to gene pairs and count chromosome correspondences among the pairs."
  join <(perl -pe 's/pan\d+__//g' 24_pan_fasta_cluster_pruned.tsv | sort -k1,1) \
       22_syn_pan_aug_extra_core_posn.gene_chr.hsh | perl -pe 's/ /\t/g' | sort -k2,2 | 
     join -1 2 -2 1 - 22_syn_pan_aug_extra_core_posn.gene_chr.hsh | 
     awk 'BEGIN{IGNORECASE=1; OFS="\t"} 
          $3!~/cont|scaff|sc|pilon|ctg|contig|tig|mito|mt$|cp$|pt$|chl|unanchor|unkn/ && \
          $4!~/cont|scaff|sc|pilon|ctg|contig|tig|mito|mt$|cp$|pt$|chl|unanchor|unkn/ \
          {print $3, $4}' | perl -pe 's/^\S+\.\D+(\d+)\s+\S+\.\D+(\d+)/$1\t$2/' |
     awk -v OFS="\t" '$1<=$2 {print $1, $2} $1>$2 {print $2, $1}' |
     sort | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' | sort -k3nr,3nr > observed_chr_pairs.tsv
}

run_summarize() {
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"
 
  # Determine if the sequence looks like nucleotide or like protein.
  someseq=$(head ${WORK_DIR}/07_pan_fasta.$fna | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
  SEQTYPE=$(check_seq_type "${someseq}")
  case "$SEQTYPE" in 
    3) ST="NUC" ;;
    1) ST="PEP" ;;
  esac

  param_string="${ST}.id${clust_iden}.cov${clust_cov}.cns${consen_iden}.ext${extra_iden}.I${mcl_inflation}"
  full_out_dir=$(echo "$out_dir_base.$param_string" | perl -pe 's/0\.(\d+)/$1/g')
  stats_file=${full_out_dir}/stats.$param_string.txt

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  for file in 06_syn_pan.clust.tsv, 06_syn_pan.hsh.tsv, 22_syn_pan_aug_extra_core_posn.hsh.tsv, \
              observed_chr_pairs.tsv, consen_${consen_prefix}.tsv, 21_pan_fasta_clust_rep_seq.fna, \
              23_syn_pan_cent_merged_core_posn.fna, 24_syn_pan_cent_merged_core_posn_trimd.fna, \
              23_syn_pan_cent_merged_core_posn_prot.faa; do
    if [ -f ${WORK_DIR}/$file ]; then
      cp ${WORK_DIR}/$file ${full_out_dir}/
    else 
      echo "Warning: couldn't find file ${WORK_DIR}/$file; skipping"
    fi
  done

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

  printf "Parameter  \tvalue\n" >> ${stats_file}
  for key in ${pandagma_conf_params}; do
    printf '%-15s\t%s\n' ${key} "${!key}" >> ${stats_file}
  done

  printf "\nOutput directory for this run:\t${full_out_dir}\n" >> ${stats_file}

  printf '%-20s\t%s\n' "Statistic" "value" >> ${stats_file}

  printf "\n== Initial clusters (containing only genes within synteny blocks)\n" >> ${stats_file}
  let "clusters=$(wc -l < ${full_out_dir}/06_syn_pan.clust.tsv)"
  printf '%-20s\t%s\n' "num_of_clusters" $clusters >> ${stats_file}

  let "largest=$(awk 'NR == 1 {print NF-1; exit}' ${full_out_dir}/06_syn_pan.clust.tsv)"
  printf '%-20s\t%s\n' "largest_cluster" $largest >> ${stats_file}

  let "mode=$(awk "{print NF-1}" ${full_out_dir}/06_syn_pan.clust.tsv | \
    uniq -c | sort -n | tail -1 | awk '{print $2}')"
  printf '%-20s\t%s\n' "modal_clst_size" $mode >> ${stats_file}

  let "num_at_mode=$(awk "{print NF-1}" ${full_out_dir}/06_syn_pan.clust.tsv | \
    uniq -c | sort -n | tail -1 | awk '{print $1}')"
  printf '%-20s\t%s\n' "num_at_mode" $num_at_mode >> ${stats_file}
  
  let "seqs_clustered=$(wc -l ${full_out_dir}/06_syn_pan.hsh.tsv | awk '{print $1}')"
  printf '%-20s\t%s\n' "seqs_clustered" $seqs_clustered >> ${stats_file}
  
#
  if [ -f ${full_out_dir}/12_syn_pan_aug.clust.tsv ]; then
    printf "\n== Augmented clusters (unanchored sequences added to the initial clusters)\n" >> ${stats_file}
    let "clustersA=$(wc -l < ${full_out_dir}/12_syn_pan_aug.clust.tsv)"
    printf '%-20s\t%s\n' "num_of_clusters" $clustersA >> ${stats_file}

    let "largestA=$(awk "{print NF-1}" ${full_out_dir}/12_syn_pan_aug.clust.tsv | sort -n | tail -1)"
    printf '%-20s\t%s\n' "largest_cluster" $largestA >> ${stats_file}

    let "modeA=$(awk "{print NF-1}" ${full_out_dir}/12_syn_pan_aug.clust.tsv | \
      uniq -c | sort -n | tail -1 | awk '{print $2}')"
    printf '%-20s\t%s\n' "modal_clst_size" $modeA >> ${stats_file}

    let "numA_at_mode=$(awk "{print NF-1}" ${full_out_dir}/12_syn_pan_aug.clust.tsv | \
      sort -n | uniq -c | sort -n | tail -1 | awk '{print $1}')"
    printf '%-20s\t%s\n' "num_at_mode" $numA_at_mode >> ${stats_file}
    
    let "seqs_clustered=$(cat ${full_out_dir}/12_syn_pan_aug.clust.tsv | awk '{ct+=(NF-1)} END{print ct}')"
    printf '%-20s\t%s\n' "seqs_clustered" $seqs_clustered >> ${stats_file}
  fi

#
  if [ -f ${full_out_dir}/17_syn_pan_aug_extra.clust.tsv ]; then
    printf "\n== Augmented-extra clusters (sequences from extra annotation sets have been added)\n" >> ${stats_file}
    let "clustersB=$(wc -l < ${full_out_dir}/17_syn_pan_aug_extra.clust.tsv)"
    printf '%-20s\t%s\n' "num_of_clusters" $clustersB >> ${stats_file}

    let "largestB=$(awk "{print NF-1}" ${full_out_dir}/17_syn_pan_aug_extra.clust.tsv | sort -n | tail -1)"
    printf '%-20s\t%s\n' "largest_cluster" $largestB >> ${stats_file}

    let "modeB=$(awk "{print NF-1}" ${full_out_dir}/17_syn_pan_aug_extra.clust.tsv | \
      sort -n | uniq -c | sort -n | tail -1 | awk '{print $2}')"
    printf '%-20s\t%s\n' "modal_clst_size" $modeB >> ${stats_file}

    let "numB_at_mode=$(awk "{print NF-1}" ${full_out_dir}/17_syn_pan_aug_extra.clust.tsv | \
      sort -n | uniq -c | sort -n | tail -1 | awk '{print $1}')"
    printf '%-20s\t%s\n' "num_at_mode" $numB_at_mode >> ${stats_file}
    
    let "seqs_clustered=$(wc -l ${full_out_dir}/22_syn_pan_aug_extra_core_posn.hsh.tsv | awk '{print $1}')"
    printf '%-20s\t%s\n' "seqs_clustered" $seqs_clustered >> ${stats_file}
  fi

  # Report the starting files:
    printf "\n== Starting fasta files\n" >> ${stats_file}
    cat ${WORK_DIR}/stats/tmp.fasta_list >> ${stats_file}

  # Print per-annotation-set coverage stats (sequence counts, sequences retained)
  cut -f2 ${WORK_DIR}/22_syn_pan_aug_extra_core_posn.hsh.tsv | 
    perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' |
    sort | uniq -c | awk '{print $2 "\t" $1}' > ${WORK_DIR}/stats/tmp.gene_count_end

  # tmp.gene_count_start was generated during run_ingest

  printf "\n== Start and end gene counts per annotation set\n" >> ${stats_file}

  join ${WORK_DIR}/stats/tmp.gene_count_start ${WORK_DIR}/stats/tmp.gene_count_end | 
    awk 'BEGIN{print "\tStart\tEnd\tKept\tAnnotation_name"} 
        { printf "\t%i\t%i\t%2.1f\t%s\n", $2, $3, 100*($3/$2), $1 }'  >> ${stats_file}

  # histograms
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

  if [ -f ${full_out_dir}/17_syn_pan_aug_extra.clust.tsv ]; then
    printf "\nCounts of augmented-extra clusters by cluster size, file 17_syn_pan_aug_extra.clust.tsv:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/17_syn_pan_aug_extra.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi

  # Counts per accession
  if [ -f ${full_out_dir}/17_syn_pan_aug_extra.counts.tsv ]; then
    printf "\nFor all annotation sets, counts of genes-in-orthogroups and counts of orthogroups-with-genes:\n" >> ${stats_file}
    printf "genes-in-OGs\tOGs-w-genes\tOGs-w-genes/genes\tpct-non-null-OGs\tpct-null-OGs\tannotation-set\n" >> ${stats_file}
    cat ${full_out_dir}/17_syn_pan_aug_extra.counts.tsv | transpose.pl | 
      perl -lane 'next if ($.<=3); 
        $ct=0; $sum=0; $nulls=0; $OGs=0;
        for $i (@F[1..(@F-1)]){
          $OGs++;
          if ($i>0){$ct++; $sum+=$i}
          if ($i==0){$nulls++}
        }; 
        printf("%d\t%d\t%.2f\t%.2f\t%.2f\t%s\n", $sum, $ct, 100*$ct/$sum, 100*($OGs-$nulls)/$OGs, 100*$nulls/$OGs, $F[0])' \
        >> ${stats_file}
  fi

  echo
  cat ${stats_file}
}

#
# top-level command functions
#

run() {
  RUN_DOC="""Run an analysis step

Usage:
   $scriptname run [STEP]

Steps:
   If STEP is not set, the following steps will be run in order,
   otherwise the step is run by itself:
              ingest - Get info from matching GFF and FNA files
              mmseqs - Run mmseqs for all gene sets
              filter - Select gene matches from indicated chromosome pairings
          dagchainer - Compute Directed Acyclic Graphs
                 mcl - Calculate Markov clusters
            consense - Calculate a consensus sequences from each pan-gene set, 
                       adding sequences missed in the first clustering round.
           add_extra - Add other gene model sets to the primary clusters. Useful for adding
                       annotation sets that may be of lower or uncertain quality.
       name_pangenes - Assign pan-gene names with consensus chromosome and position information.
      filter_to_core - Calculate orthogroup composition and filter fasta files to core orthogroups.
      calc_chr_pairs - Calculate observed chromosome pairs.
           summarize - Compute synteny stats.
"""
  commandlist="ingest mmseqs filter dagchainer mcl consense add_extra \
               filter_to_core name_pangenes calc_chr_pairs summarize"
  . "${CONF}"
  canonicalize_paths
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
# Command-line interpreter.
#
if [ "$#" -eq 0 ]; then
  echo >&2 "$TOP_DOC"
  exit 0
fi
#
command="$1"
shift 1
case $command in
"run")
  run $@
  ;;
"version")
  version $@
  ;;
*)
  echo >&2 "ERROR -- command \"$command\" not recognized."
  exit 1
  ;;
esac

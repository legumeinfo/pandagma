#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates pan-gene 
# clusters using the programs mmseqs, dagchainer, and mcl. 
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2021
#
scriptname=`basename "$0"`
version="2021-11-16"
set -o errexit -o errtrace -o nounset -o pipefail

export NPROC=${NPROC:-1}
export MMSEQS_NUM_THREADS=${NPROC} # mmseqs otherwise uses all cores by default
export PANDAGMA_CONF=${PANDAGMA_CONF:-${PWD}/pandagma.conf}
export PANDAGMA_WORK_DIR=${PANDAGMA_WORK_DIR:-${PWD}/work}

# mmseqs uses a significant number of threads on its own. Set a maximum, which may be below NPROC.
MMSEQSTHREADS=$(( 10 < ${NPROC} ? 10 : ${NPROC} ))

pandagma_conf_params='clust_iden clust_cov consen_iden extra_iden mcl_inflation dagchainer_args out_dir_base consen_prefix extra_stats'

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
 run calc_chr_pairs - Report observed chromosome pairs; useful for preparing expected_chr_matches.tsv
      run summarize - Move results into output directory, and report summary statistics.

Variables in pandagma config file (Set the config with the PANDAGMA_CONF environment variable)
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.60]
        consen_iden - Minimum identity threshold for consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of \"extra\" annotations [90]
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]
      consen_prefix - Prefix to use in names for genomic ordered consensus IDs [Genus.pan1]
       out_dir_base - Base name for the output directory [default: './out']
        extra_stats - Flag (yes/no) to calculate coverage stats; use only for genes with IDs that have a four-part
                       prefix like glyma.Wm82.gnm1.ann1 [no]


Environment variables:
         PANDAGMA_CONF - Path of the pandagma config file, default:
                       \"${PANDAGMA_CONF}\"
    PANDAGMA_WORK_DIR - Location of working files, default:
                       \"${PANDAGMA_WORK_DIR}\"
              NPROC - Number of processors to use (default 1)
"""
#
# Helper functions begin here
#
canonicalize_paths() {
  fasta_files=($(realpath --canonicalize-existing "${fasta_files[@]}"))
  annotation_files=($(realpath --canonicalize-existing "${annotation_files[@]}"))
  if (( ${#fasta_files_extra[@]} > 0 ))
  then
    fasta_files_extra=($(realpath --canonicalize-existing "${fasta_files_extra[@]}"))
    annotation_files_extra=($(realpath --canonicalize-existing "${annotation_files_extra[@]}"))
  fi
  readonly chr_match_list=${expected_chr_matches:+$(realpath "${expected_chr_matches}")}
  readonly submit_dir=${PWD}

  fasta_file=$(basename "${fasta_files[0]}" .gz)
  faext="${fasta_file##*.}"
}
cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}
check_seq_type() {
  someseq=${1}
  proportion_nuc=$(echo $someseq | fold -w1 | awk '$1~/[ATCGN]/ {nuc++} $1!~/[ATCGN]/ {not++} END{print nuc/(nuc+not)}')
  if (( $(bc <<< "$proportion_nuc > 0.9") )); then
    echo "NUC"
  else
    echo "PEP"
  fi
}

### run functions ###

run_ingest() {
# Add positional information from GFF3 or 4- or 6-column BED to FASTA IDs
# BED start coordinate converted to 1-based
  cd "${PANDAGMA_WORK_DIR}"
  echo; echo "Run ingest: from fasta and gff or bed data, create fasta with IDs containing positional info."
  
  mkdir -p 02_fasta 01_posn_hsh stats

    # Prepare the tmp.gene_count_start to be joined, in run_summarize, with tmp.gene_count_end.
    # Requires chr names to be prefixed with the annotation name, separated by dot from chr/scaff number
    # e.g. glyma.FiskebyIII.gnm1.Gm11  or  Zm-B73_GRAMENE-4.chr1
    cat /dev/null > stats/tmp.gene_count_start
    cat /dev/null > stats/tmp.fasta_list
    start_time=`date`
    printf "Run started at: $start_time\n" > stats/tmp.timing

  for (( file_num = 0; file_num < ${#fasta_files[@]} ; file_num++ )); do
    file_base=$(basename ${fasta_files[file_num]%.*})
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | gff_or_bed_to_hash.awk > 01_posn_hsh/$file_base.hsh
    hash_into_fasta_id.pl -fasta "${fasta_files[file_num]}" \
                          -hash 01_posn_hsh/$file_base.hsh \
                          -out 02_fasta/$file_base
    annot_name=$(head -1 02_fasta/$file_base | 
                 perl -pe 's/>(.+)__.+__\d+__\d+$/$1/' | perl -pe 's/(.+)\.[^.]+$/$1/')
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
      cat_or_zcat "${annotation_files_extra[file_num]}" | gff_or_bed_to_hash.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -fasta "${fasta_files_extra[file_num]}" \
                            -hash 01_posn_hsh/$file_base.hsh \
                            -out 02_fasta/$file_base
      annot_name=$(head -1 02_fasta/$file_base | perl -pe 's/>(.+)__.+__\d+__\d+$/$1/' | perl -pe 's/(.+)\.[^.]+$/$1/')
      cat_or_zcat "${fasta_files_extra[file_num]}" | 
        awk -v ANNOT=$annot_name '$1~/^>/ {ct++} END{print ANNOT "\t" ct}' >> stats/tmp.gene_count_start
      echo "Extra:  $file_base" >> stats/tmp.fasta_list
    done
  fi

  sort -o stats/tmp.gene_count_start stats/tmp.gene_count_start
}

run_mmseqs() {
  # Do mmseqs clustering on all pairings of the main annotation sets (not the extra ones though)
  cd "${PANDAGMA_WORK_DIR}"
  echo; echo "Run mmseqs -- at ${clust_iden} percent identity and minimum of ${clust_cov}% coverage."
  #

  mkdir -p 03_mmseqs 03_mmseqs_tmp
  for (( file1_num = 0; file1_num < ${#fasta_files[@]} ; file1_num++ )); do
    qry_base=$(basename ${fasta_files[file1_num]%.*} .$faext)
    for (( file2_num = $( expr $file1_num + 1 ); file2_num < ${#fasta_files[@]} ; file2_num++ )); do
      sbj_base=$(basename ${fasta_files[file2_num]%.*} .$faext)
      echo "  Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
      { cat 02_fasta/$qry_base.$faext 02_fasta/$sbj_base.$faext ; } |
        mmseqs easy-cluster stdin 03_mmseqs/${qry_base}.x.${sbj_base} 03_mmseqs_tmp \
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
  cd "${PANDAGMA_WORK_DIR}"
  mkdir -p 04_dag
  if [[ -f ${chr_match_list} ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns from file ${chr_match_list}"
    for mmseqs_path in 03_mmseqs/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      echo "  $outfilebase"
      cat ${mmseqs_path} | filter_mmseqs_by_chroms.pl -chr_pat ${chr_match_list} > 04_dag/${outfilebase}_matches.tsv &

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
  cd "${PANDAGMA_WORK_DIR}"
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
  cd "${PANDAGMA_WORK_DIR}"
  printf "\nCalculate clusters. use Markov clustering with inflation parameter $mcl_inflation and ${NPROC} threads\n"
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
  cd "${PANDAGMA_WORK_DIR}"
  mkdir -p 07_pan_fasta lists

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  echo "    Fasta file:" "${fasta_files[@]}"
  get_fasta_from_family_file.pl "${fasta_files[@]}" -fam 06_syn_pan.clust.tsv -out 07_pan_fasta

  cat /dev/null > 07_pan_fasta.$faext
  for path in 07_pan_fasta/*; do
    pan_file=`basename $path`
    cat $path | awk -v panID=$pan_file '
                      $1~/^>/ {print ">" panID "__" substr($0,2) }
                      $1!~/^>/ {print $1}
                    ' >> 07_pan_fasta.$faext
  done

  echo "  Pick a representative sequence for each pangene set - as a sequence with the median length for that set."
  pick_family_rep.pl -in 07_pan_fasta.$faext -out 08_pan_fasta_clust_rep_seq.$faext

  echo "  Get sorted list of all genes, from the original fasta files"
  cat_or_zcat "${fasta_files[@]}" | awk '/^>/ {print substr($1,2)}' | sort > lists/09_all_genes

  echo "  Get sorted list of all clustered genes"
  awk '$1~/^>/ {print $1}' 07_pan_fasta/* | sed 's/>//' | sort > lists/09_all_clustered_genes

  echo "  Get list of genes not in clusters"
  comm -13 lists/09_all_clustered_genes lists/09_all_genes > lists/09_genes_not_in_clusters

  echo "  Retrieve the non-clustered genes"
  cat_or_zcat "${fasta_files[@]}" |
    get_fasta_subset.pl -in /dev/stdin \
                        -out 09_genes_not_in_clusters.$faext \
                        -lis lists/09_genes_not_in_clusters -clobber

  echo "  Search non-clustered genes against genes already clustered."

  # Check sequence type (in case this run function is called separately from the usually-prior ones)
  someseq=$(head 07_pan_fasta.$faext | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
  SEQTYPE=$(check_seq_type "${someseq}")
  if [[ "$SEQTYPE" == "NUC" ]]; then MMST=3; else MMST=1; fi

  mmseqs easy-search 09_genes_not_in_clusters.$faext \
                     07_pan_fasta.$faext \
                     10_unclust.x.07_pan_fasta.m8 \
                     03_mmseqs_tmp \
                     --search-type ${MMST} --cov-mode 5 -c ${clust_cov} 1>/dev/null 

  echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
  echo "  Use the \"main set\" $clust_iden threshold. Re-merge sub-clusters."
  top_line.awk 10_unclust.x.07_pan_fasta.m8 | 
    awk -v IDEN=${clust_iden} '$3>=IDEN {print $2 "\t" $1}' | perl -pe 's/^(pan\d+)__\S+/$1/' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk >  11_syn_pan_leftovers.clust.tsv

  echo "  Retrieve sequences for the leftover genes"
  mkdir -p 11_pan_leftovers
  get_fasta_from_family_file.pl "${fasta_files[@]}" \
    -fam 11_syn_pan_leftovers.clust.tsv -out 11_pan_leftovers/

  echo "  Make augmented cluster sets"
  cat /dev/null > 12_syn_pan_aug.clust.tsv

  echo "  Make augmented cluster sets"
  augment_cluster_sets.awk leftovers_dir=11_pan_leftovers 07_pan_fasta/* > 12_syn_pan_aug.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 12_syn_pan_aug.clust.tsv |
    cat > 12_syn_pan_aug.hsh.tsv
}

run_add_extra() {
  if (( ${#fasta_files_extra[@]} > 0 ))
  then # handle the "extra" annotation files
    echo; echo "Add extra annotation sets to the augmented clusters, by homology"
    echo "  Search non-clustered genes against pan-gene consensus sequences"
    cd "${PANDAGMA_WORK_DIR}"
    mkdir -p 13_extra_out_dir 13_pan_aug_fasta

    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    get_fasta_from_family_file.pl "${fasta_files[@]}" -fam 12_syn_pan_aug.clust.tsv -out 13_pan_aug_fasta
    
    cat /dev/null > 13_pan_aug_fasta.$faext
    for path in 13_pan_aug_fasta/*; do
      pan_file=`basename $path`
      cat $path | awk -v panID=$pan_file '
                        $1~/^>/ {print ">" panID "__" substr($0,2) }
                        $1!~/^>/ {print $1}
                      ' >> 13_pan_aug_fasta.$faext
    done

    # Check sequence type (in case this run function is called separately from the usually-prior ones)
    someseq=$(head 07_pan_fasta.$faext | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
    SEQTYPE=$(check_seq_type "${someseq}")
    if [[ "$SEQTYPE" == "NUC" ]]; then MMST=3; else MMST=1; fi

    for path in "${fasta_files_extra[@]}"; do
      fasta_file=`basename ${path%.*}`
      echo "Extra: $fasta_file"
      mmseqs easy-search "${path}" \
                         13_pan_aug_fasta.$faext \
                         13_extra_out_dir/${fasta_file}.x.all_cons.m8 \
                         03_mmseqs_tmp/ \
                         --search-type ${MMST} --cov-mode 5 -c ${clust_cov} 1>/dev/null & # background

      if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    wait # wait for jobs to finish
  
    echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
    echo "  Use identity threshold extr_iden: $extra_iden. Re-merge sub-clusters."
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

    echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
    perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 17_syn_pan_aug_extra.clust.tsv |
      cat > 17_syn_pan_aug_extra.hsh.tsv

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
  
    cat /dev/null > 20_pan_fasta.$faext
    for path in 19_pan_aug_leftover_merged/*; do
      pan_file=`basename $path`
      cat $path | awk -v panID=$pan_file '
                        $1~/^>/ {print ">" panID "__" substr($0,2) }
                        $1!~/^>/ {print $1}
                      ' >> 20_pan_fasta.$faext
    done

    echo "  Pick a representative sequence for each pangene set - as a sequence with the median length for that set."
    pick_family_rep.pl -in 20_pan_fasta.$faext -out 21_pan_fasta_clust_rep_seq.$faext
    perl -pi -e 's/__/\t/' 21_pan_fasta_clust_rep_seq.$faext
  
  else  # no "extra" fasta files, so just promote the syn_pan_aug files as syn_pan_aug_extra
        # TO DO: handle this better, i.e. report that no "extra" files were provided and skip these "aug" files.
    cp 07_pan_fasta.$faext 20_pan_fasta.$faext
    cp 12_syn_pan_aug.clust.tsv 17_syn_pan_aug_extra.clust.tsv
    cp 12_syn_pan_aug.hsh.tsv 17_syn_pan_aug_extra.hsh.tsv
    cp 08_pan_fasta_clust_rep_seq.$faext 21_pan_fasta_clust_rep_seq.$faext
    #cp syn_pan_cent.$faext syn_pan_cent_merged.$faext
fi
}

run_name_pangenes() {
  cd "${PANDAGMA_WORK_DIR}"

  echo "  Add positional information to the hash output."
  join -a 1 -1 2 -2 1 <(sort -k2,2 17_syn_pan_aug_extra.hsh.tsv) <(cat 01_posn_hsh/*hsh | sort -k1,1) | 
    perl -pe 's/__/\t/g; s/ /\t/g' | awk -v OFS="\t" '{print $2, $1, $3, $5, $6}' |
    sort -k1,1 -k2,2 > 18_syn_pan_aug_extra_posn.hsh.tsv

  echo "  Calculate consensus pan-gene positions"
  cat 18_syn_pan_aug_extra_posn.hsh.tsv | consen_pangene_posn.pl -pre ${consen_prefix}.chr |
    sort -k2,2 -k3n,3n | name_ordered_genes.awk > consen_${consen_prefix}.tsv

  echo "  Reshape defline into a hash, e.g. pan47789	Glycine.pan3.chr01__Glycine.pan3.chr01_000100__45224__45786"
  cat consen_${consen_prefix}.tsv | 
    perl -pe 's/^(\S+)\t([^.]+\.[^.]+)\.(chr\d+)_(\d+)\t(\d+)\t(\d+)/$1\t$1__$2.$3_$4__$5__$6/' > consen_posn.hsh

  echo "  Hash position information into fasta file"
  hash_into_fasta_id.pl -fasta 21_pan_fasta_clust_rep_seq.$faext -hash consen_posn.hsh -keep_definition \
    -out_file 21_pan_fasta_clust_rep_seq_posnTMP.$faext

  # Reshape defline, and sort by position
  fasta_to_table.awk 21_pan_fasta_clust_rep_seq_posnTMP.$faext | sed 's/__/\t/g; s/ /\t/' | 
    perl -pe 's/^(\w+)\s+(.+)\.(chr\d+)_(\d+)\s+/$1\t$2\t$3\t$4\t/' | sed 's/chr//' | sort -k3n -k4n |
    awk '{print ">" $2 ".chr" $3 "_" $4 " " $1 " " $5 " " $6 " " $7; print $8}' > 22_syn_pan_cent_merged_posn.$faext
  rm 21_pan_fasta_clust_rep_seq_posnTMP.fa
 
  echo "  Re-cluster, to identify neighboring genes that are highly similar"
    mmseqs easy-cluster 22_syn_pan_cent_merged_posn.$faext 23_pan_fasta 03_mmseqs_tmp \
    --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null

  # Parse mmseqs clusters into pairs of genes that are similar and ordinally close
  close=10 # genes within a neighborhood of +-close genes among the consensus-orderd genes
  cat 23_pan_fasta_cluster.tsv | perl -pe 's/\S+\.chr(\d+)_(\d+)/$1\t$2/g' | 
    awk -v CLOSE=$close -v PRE=$consen_prefix '$1==$3 && sqrt(($2/100-$4/100)^2)<=CLOSE \
             {print PRE ".chr" $1 "_" $2 "\t" PRE ".chr" $3 "_" $4}' > 23_pan_fasta_cluster.closepairs

  # Cluster the potential near-duplicate neighboring paralogs
  mcl  23_pan_fasta_cluster.closepairs --abc -o 23_pan_fasta_cluster.close.clst
 
  # Keep the first gene from the cluster and discard the rest. Do this by removing the others.
  cat 23_pan_fasta_cluster.close.clst | awk '{$1=""}1' | 
    awk '{$1=$1}1' | tr ' ' '\n' | sort -u > lists/lis.clust_genes_remove

  # Remove neighboring close paralogs, leaving one
  get_fasta_subset.pl -in 22_syn_pan_cent_merged_posn.$faext -out 22_syn_pan_cent_merged_posn_trimd.$faext \
    -lis lists/lis.clust_genes_remove -xclude -clobber
}

run_calc_chr_pairs() {
  cd "${PANDAGMA_WORK_DIR}"
  echo "Generate a report of observed chromosome pairs"

  echo "  Identify gene pairs, ussing mmseqs --easy_cluster"
    mmseqs easy-cluster 20_pan_fasta.$faext 20_pan_fasta_clust 03_mmseqs_tmp \
    --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null

  # Extract chromosome-chromosome correspondences
  cut -f2,3 18_syn_pan_aug_extra_posn.hsh.tsv | sort -k1,1 > 18_syn_pan_aug_extra_posn.gene_chr.hsh

  # From mmseqs cluster table, prune gene pairs, keeping only those in the same pangene cluster
  cat 20_pan_fasta_clust_cluster.tsv | perl -pe 's/__/\t/g' | awk '$1==$3' |
    perl -pe 's/(pan\d+)\t/$1__/g' > 20_pan_fasta_clust_cluster_pruned.tsv

  # Join chromosome numbers to gene pairs and count chromosome correspondences among the pairs.
  join <(perl -pe 's/pan\d+__//g' 20_pan_fasta_clust_cluster_pruned.tsv | sort -k1,1) \
       18_syn_pan_aug_extra_posn.gene_chr.hsh | perl -pe 's/ /\t/g' | sort -k2,2 | 
     join -1 2 -2 1 - 18_syn_pan_aug_extra_posn.gene_chr.hsh | 
     awk 'BEGIN{IGNORECASE=1; OFS="\t"} 
          $3!~/cont|scaff|sc|pilon|ctg|contig|mito|mt$|cp$|pt$|chl|unanchor|unkn/ && \
          $4!~/cont|scaff|sc|pilon|ctg|contig|mito|mt$|cp$|pt$|chl|unanchor|unkn/ \
          {print $3, $4}' | perl -pe 's/^\S+\.\D+(\d+)\s+\S+\.\D+(\d+)/$1\t$2/' |
     awk -v OFS="\t" '$1<=$2 {print $1, $2} $1>$2 {print $2, $1}' |
     sort | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' | sort -k3nr,3nr > observed_chr_pairs.tsv
}

run_summarize() {
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"
 
  # Determine if the sequence looks like nucleotide or like protein. Set default type as NUC
  someseq=$(head ${PANDAGMA_WORK_DIR}/07_pan_fasta.$faext | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
  SEQTYPE=$(check_seq_type "${someseq}")

  param_string="${SEQTYPE}.id${clust_iden}.cov${clust_cov}.cns${consen_iden}.ext${extra_iden}.I${mcl_inflation}"
  full_out_dir=$(echo "$out_dir_base.$param_string" | perl -pe 's/0\.(\d+)/$1/g')
  stats_file=${full_out_dir}/stats.$param_string.txt

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  cp ${PANDAGMA_WORK_DIR}/06_syn_pan.clust.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/06_syn_pan.hsh.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/*_syn_pan_aug*.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/observed_chr_pairs.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/consen_${consen_prefix}.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/22_syn_pan_cent_merged_posn.$faext ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/22_syn_pan_cent_merged_posn_trimd.$faext ${full_out_dir}/

  printf "Run of program $scriptname, version $version\n" > ${stats_file}

  end_time=`date`
  cat ${PANDAGMA_WORK_DIR}/stats/tmp.timing >> ${stats_file}
  printf "Run ended at:   $end_time\n\n" >> ${stats_file}

  # Report sequence type 
  if [[ "$SEQTYPE" == "NUC" ]]; then
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
    
    let "seqs_clustered=$(wc -l ${full_out_dir}/12_syn_pan_aug.hsh.tsv | awk '{print $1}')"
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
    
    let "seqs_clustered=$(wc -l ${full_out_dir}/18_syn_pan_aug_extra_posn.hsh.tsv | awk '{print $1}')"
    printf '%-20s\t%s\n' "seqs_clustered" $seqs_clustered >> ${stats_file}
  fi

  # Report the starting files:
    printf "\n== Starting fasta files\n" >> ${stats_file}
    cat ${PANDAGMA_WORK_DIR}/stats/tmp.fasta_list >> ${stats_file}

  # Print per-annotation-set coverage stats (sequence counts, sequences retained), if stats-extra flag is set
  if [ ${extra_stats,,} = "yes" ]; then
    cut -f3 ${PANDAGMA_WORK_DIR}/18_syn_pan_aug_extra_posn.hsh.tsv | 
      perl -pe 's/(.+)\.([^.]+)$/$1\n/' |
      sort | uniq -c | awk '{print $2 "\t" $1}' > ${PANDAGMA_WORK_DIR}/stats/tmp.gene_count_end

    # tmp.gene_count_start was generated during run_ingest

    printf "\n== Start and end gene counts per annotation set\n" >> ${stats_file}

    join ${PANDAGMA_WORK_DIR}/stats/tmp.gene_count_start ${PANDAGMA_WORK_DIR}/stats/tmp.gene_count_end | 
      awk 'BEGIN{print "\tStart\tEnd\tKept\tAnnotation_name"} 
          { printf "\t%i\t%i\t%2.1f\t%s\n", $2, $3, 100*($3/$2), $1 }'  >> ${stats_file}
  fi

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
       name_pangenes - Assign pan-gene names with consensus chromosome and position information
      calc_chr_pairs - Calculate observed chromosome pairs
           summarize - Compute synteny stats
"""
  commandlist="ingest mmseqs filter dagchainer mcl consense add_extra name_pangenes calc_chr_pairs summarize"
  . "${PANDAGMA_CONF}"
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

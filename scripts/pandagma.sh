#!/bin/bash
#
# Configuration and run script for "pandagma", which generates pan-gene clusters using the programs 
# mmseqs, dagchainer, mcl, and vsearch. These are used to do the initial clustering, 
# synteny-finding, and re-clustering.
# Authors: Steven Cannon, Joel Berendzen, 2020-2021
#
scriptname=`basename "$0"`
version="0.9.1"
set -o errexit -o nounset -o pipefail

export NPROC=${NPROC:-1}
export MMSEQS_NUM_THREADS=${NPROC} # mmseqs otherwise uses all cores by default
export PANDAGMA_CONF=${PANDAGMA_CONF:-${PWD}/pandagma.conf}

readonly work_dir=${PANDAGMA_WORK_DIR:-${PWD}/work}
readonly data_dir=${PANDAGMA_DATA_DIR:-${PWD}/data}
readonly data_extra_dir=${PANDAGMA_DATA_EXTRA_DIR:-${PWD}/data_extra}

dag_dir="${work_dir}/dag"
mmseqs_dir="${work_dir}/mmseqs"
mmseqs_tmp_dir="${work_dir}/mmseqs_tmp"
work_data_dir="${work_dir}/data"
work_data_extra_dir="${work_dir}/data_extra"
work_extra_out_dir="${work_dir}/extra_out_dir"
pan_fasta_dir="${work_dir}/pan_fasta"
pan_consen_dir="${work_dir}/pan_consen"
leftovers_dir="${work_dir}/pan_leftovers"
leftovers_extra_dir="${work_dir}/pan_leftovers_extra"

pandagma_conf_params='clust_iden clust_cov consen_iden mcl_inflation dagchainer_args fasta_ext gff_ext pan_prefix out_dir_base'

error_exit() {
  echo >&2 "ERROR -- unexpected exit from ${BASH_SOURCE} script at line:"
  echo >&2 "   $BASH_COMMAND"
}
trap error_exit EXIT

TOP_DOC="""Compute pan-gene clusters using the programs mmseqs, dagchainer, and mcl, and
additional pre- and post-refinement steps.

Usage:
        pandagma.sh SUBCOMMAND [SUBCOMMAND_OPTIONS]

By default, name-matched primary coding sequence (fasta) and annotation (GFF) files are expected 
in the data/ directory, within the working directory from where this script is called.
(To set a different data directory name, see discussion of environment variables below).
Example of name-matched files within the data/ directory:
  accession1.fna accession1.gff3   
  accession2.fna accession2.gff3  
  accessionXYZ.fna accessionXYZ.gff3

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
         run ingest - Prepare the assembly and annotation files for analysis
         run mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies
         run filter - Filter the synteny results for chromosome pairings, returning gene pairs.
     run dagchainer - Run DAGchainer to filter for syntenic blocks
            run mcl - Derive clusters, with Markov clustering
       run consense - Calculate a consensus sequences from each pan-gene set, 
                      adding sequences missed in the first clustering round.
      run add_extra - Add other gene model sets to the primary clusters. Useful for adding
                      annotation sets that may be of lower or uncertain quality.
      run summarize - Move results into output directory, and report summary statistics.
              clean - Delete work directory

Variables in pandagma config file:
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.98]
          clust_cov - Minimum coverage for mmseqs clustering [0.75]
        consen_iden - Minimum identity threshold for vsearch consensus generation [0.80]
          fasta_ext - Extension of FASTA files
            gff_ext - Extension of GFF files
         pan_prefix - Prefix to use as a prefix for pangene clusters [default: pan]
       out_dir_base - base name for the output directory [default: './out']
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]

Environment variables:
         PANDAGMA_CONF - Path of the pandagma config file, default:
                       \"${PANDAGMA_CONF}\"
    PANDAGMA_WORK_DIR - Location of working files, default:
                       \"${work_dir}\"
    PANDAGMA_DATA_DIR - Location of intput data files (CDS fasta and GFFs), default:
                       \"${data_dir}\"
PANDAGMA_DATA_EXTRA_DIR - Location of data files to be added after primary cluster construction, default:
                       \"${data_extra_dir}\"
              NPROC - Number of processors to use (default 1)
"""
#
# Helper functions begin here
#
cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}
# usage: make_augmented_cluster_sets leftovers_dir=${leftovers_dir} ${pan_fasta_dir}/* > syn_pan_augmented.clust.tsv
make_augmented_cluster_sets() {
  awk 'function add_leftovers() {
        leftovers_file = leftovers_dir "/" path[nf]
        while(getline fasta_line < leftovers_file == 1)
            if (fasta_line ~ /^>/) printf("\t%s", substr(fasta_line,2))
        close(leftovers_file)
        printf("\n")
    }
    FNR == 1 {
        if (NR != 1) add_leftovers() # previous pan cluster
        nf = split(FILENAME, path, "/")
        printf("%s", path[nf])
    }
    /^>/ { printf("\t%s", substr($1,2)) }
    END { add_leftovers() }' "$@"
}
#
# run functions
#
run_ingest() {
  # Prepare the assembly and annotation files for analysis. Add positional info to gene IDs.
  . "${PANDAGMA_CONF}"
  gff_files="$(ls ${data_dir}/*.${gff_ext})"
  n_fasta=`ls ${data_dir}/*.${fasta_ext} | wc -l`
  printf "\nIngest -- combine data from ${n_fasta} ${gff_ext} and ${fasta_ext} files\n"

  if ! command -v hash_into_fasta_id.pl &> /dev/null ; then
    echo >&2 "ERROR: hash_into_fasta_id.pl could not be found. "
    echo >&2 "       Are the necessary computational tools available? The required tools are:"
    echo >&2 "           perl-bioperl-core, mmseqs2, dagchainer, mcl, vsearch"
    exit
  fi

  for path in $gff_files; do
    base=$(basename $path .${gff_ext})
    base_no_ann=$(echo $base | perl -pe 's/\.ann\d+\.\w+//')
    cat_or_zcat "${path}" |
      awk -v OFS="\t" '$3=="mRNA" {print $1, $4, $5, $9}' |
        perl -pe 's/ID=([^;]+);\S+/$1/' >${work_data_dir}/${base_no_ann}.bed
  done
  # add positional information to FASTA ids
  for path in ${work_data_dir}/*.bed; do
    base=$(basename $path .bed)
    awk '{print $4 "\t" $1 "__" $4 "__" $2 "__" $3}' "${path}" \
      >${work_data_dir}/${base}.hsh
    hash_into_fasta_id.pl\
      -fasta ${data_dir}/${base}.${fasta_ext} \
      -hash ${work_data_dir}/${base}.hsh \
      -suff_regex \
      >${work_data_dir}/${base}.fna
  done

  extra_gff_files="$(ls ${data_extra_dir}/*.${gff_ext})"
  printf "\nIf there are extra annotation sets (in data_extra), prepare those for analysis.\n"
  for path in $extra_gff_files; do
    base=$(basename $path .${gff_ext})
    base_no_ann=$(echo $base | perl -pe 's/\.ann\d+\.\w+//')
    cat_or_zcat "${path}" |
      awk -v OFS="\t" '$3=="mRNA" {print $1, $4, $5, $9}' |
        perl -pe 's/ID=([^;]+);\S+/$1/' >${work_data_extra_dir}/${base_no_ann}.bed
  done
  # add positional information to FASTA ids
  for path in ${work_data_extra_dir}/*.bed; do
    base=$(basename $path .bed)
    awk '{print $4 "\t" $1 "__" $4 "__" $2 "__" $3}' "${path}" \
      >${work_data_extra_dir}/${base}.hsh
    hash_into_fasta_id.pl\
      -fasta ${data_extra_dir}/${base}.${fasta_ext} \
      -hash ${work_data_extra_dir}/${base}.hsh \
      -suff_regex \
      >${work_data_extra_dir}/${base}.fna
  done
}
#
run_mmseqs() {
  # Do mmseqs clustering on all pairings of annotation sets.
  . "${PANDAGMA_CONF}"
  echo; echo "Run mmseqs -- at ${clust_iden} percent identity and minimum of ${clust_cov}% coverage."
  #
  for qry_path in ${work_data_dir}/*.fna; do
    qry_base=$(basename $qry_path .fna})
    for sbj_path in ${work_data_dir}/*.fna; do
      sbj_base=$(basename $sbj_path .fna)
      if [[ "$qry_base" > "$sbj_base" ]]; then
         echo "Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
         mmseqs easy-cluster $qry_path $sbj_path ${mmseqs_dir}/${qry_base}.x.${sbj_base} ${mmseqs_tmp_dir} \
           --min-seq-id $clust_iden \
           -c $clust_cov \
           --cov-mode 0 \
           --cluster-reassign 1>/dev/null
      fi
    done
  done
}
#
run_filter() {
  echo; echo "From mmseqs cluster output, split out the following fields: molecule, gene, start, stop."
  chr_match_list=${data_dir}/expected_chr_matches.tsv
  if [[ -f ${chr_match_list} ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns from file ${chr_match_list}"
    for mmseqs_path in ${mmseqs_dir}/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      filter_mmseqs_by_chroms.pl -chr ${chr_match_list} > ${dag_dir}/${outfilebase}_matches.tsv < "${mmseqs_path}" &

      # allow to execute up to $NPROC in parallel
      if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
    done
    wait # wait for last jobs to finish
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches.tsv file was provided, so proceeding without chromosome-pair filtering."
    for mmseqs_path in ${mmseqs_dir}/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      perl -pe 's/__/\t/g' > ${dag_dir}/${outfilebase}_matches.tsv < "${mmseqs_path}"
    done
  fi
}
#
run_dagchainer() {
  # Identify syntenic blocks, using DAGchainer
  . "${PANDAGMA_CONF}"
  echo; echo "Run DAGchainer, using args \"${dag_args}\""
  for match_path in ${dag_dir}/*_matches.tsv; do
    align_file=`basename $match_path _matches.tsv`
    echo "Running DAGchainer on comparison: $align_file"
    echo "  run_DAG_chainer.pl $dag_args -i $match_path"; echo
    # run_DAG_chainer.pl writes temp files to cwd;
    # use per-process temp directory to avoid any data race
    (
      tmpdir=$(mktemp -d)
      cd "${tmpdir}"
      run_DAG_chainer.pl $dag_args -i "${match_path}" 1>/dev/null
      rmdir ${tmpdir}
    ) &
    # allow to execute up to $NPROC in parallel
    [ $(jobs -r -p | wc -l) -ge ${NPROC} ] && wait -n
  done
  wait # wait for last jobs to finish
  #Extract single-linkage synteny anchors
  printf "matches\tscore\trev\tid1\tid2\n" >${work_dir}/synteny_blocks.tsv

  awk '$1!~/^#/ {print $2 "\t" $6}' ${dag_dir}/*_matches.tsv > ${work_dir}/homology_pairs.tsv
  awk '$1!~/^#/ {print $2 "\t" $6}' ${dag_dir}/*.aligncoords > ${work_dir}/synteny_pairs.tsv

  for path in ${dag_dir}/*.aligncoords; do
    awk '/##/ && !/reverse/ {print substr($14,0,length($14)-2) "\t" $10 "\t" 1 "\t" $3 "\t" $5}' "${path}" >> ${work_dir}/synteny_blocks.tsv
    awk '/##/ &&  /reverse/ {print substr($15,0,length($15)-2) "\t" $11 "\t" 1 "\t" $3 "\t" $5}' "${path}" >> ${work_dir}/synteny_blocks.tsv
  done
}
#
run_mcl() {
  # Calculate clusters using Markov clustering
  . "${PANDAGMA_CONF}"
  printf "\nCalculate clusters. use Markov clustering with inflation parameter $mcl_inflation and ${NPROC} threads\n"
  echo "MCL COMMAND: mcl ${work_dir}/synteny_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o ${work_dir}/tmp.syn_pan.clust.tsv"
  mcl ${work_dir}/synteny_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o ${work_dir}/tmp.syn_pan.clust.tsv \
    1>/dev/null
 
  # Add cluster IDs
  awk -v PRE=${pan_prefix} '{padnum=sprintf("%05d", NR); print PRE padnum "\t" $0}' "${work_dir}"/tmp.syn_pan.clust.tsv > ${work_dir}/syn_pan.clust.tsv

  # Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' "${work_dir}"/syn_pan.clust.tsv > ${work_dir}/syn_pan.hsh.tsv
}
#
run_consense() {
  echo; echo "Calculate a consensus sequence for each pan-gene set, using vsearch."
  echo "Then add previously unclustered sequences into an \"augmented\" pan-gene set, by homology."
  . "${PANDAGMA_CONF}"
  fasta_files="$(ls ${data_dir}/*.${fasta_ext})"

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  get_fasta_from_family_file.pl ${fasta_files} -fam ${work_dir}/syn_pan.clust.tsv -out ${pan_fasta_dir} 

  echo "  Calculate consensus sequences for each pan-gene set."
  ls ${pan_fasta_dir} | xargs -I{} -n 1 -P ${NPROC} \
    vsearch --cluster_fast ${pan_fasta_dir}/{} --id ${consen_iden} --fasta_width 0 \
            --consout ${pan_consen_dir}/{} \
            --quiet \
            --threads 1

  echo "  Combine consensus sequences into one multifasta file"

  awk 'FNR==1 { nf=split(FILENAME, FN, "/") }
         /^>/ { print ">" FN[nf] " " substr($1,2); next }
              { print }' ${pan_consen_dir}/* > ${work_dir}/syn_pan_consen.fna
      
  rm ${pan_consen_dir}/*

  echo "  Get sorted list of all genes, from the original fasta files"
  cat_or_zcat ${fasta_files} | awk '/^>/ {print substr($1,2)}' | sort > ${work_dir}/lis.all_genes

  echo "  Get sorted list of all clustered genes"
  awk '$1~/^>/ {print $1}' "${pan_fasta_dir}"/* | sed 's/>//' | sort > ${work_dir}/lis.all_clustered_genes

  echo "  Get list of genes not in clusters"
  comm -13 ${work_dir}/lis.all_clustered_genes ${work_dir}/lis.all_genes > ${work_dir}/lis.genes_not_in_clusters

  echo "  Retrieve the non-clustered genes"
  cat_or_zcat ${fasta_files} |
    get_fasta_subset.pl -in /dev/stdin \
                        -out ${work_dir}/genes_not_in_clusters.fna \
                        -lis ${work_dir}/lis.genes_not_in_clusters -clobber

  echo "  Search non-clustered genes against pan-gene consensus sequences"
  mmseqs easy-search ${work_dir}/genes_not_in_clusters.fna \
                     ${work_dir}/syn_pan_consen.fna \
                     ${work_dir}/unclust.x.all_cons.m8 ${mmseqs_tmp_dir} \
                     --search-type 3 --cov-mode 5 -c 0.5 1>/dev/null

  echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
  top_line.awk ${work_dir}/unclust.x.all_cons.m8 | 
    awk -v IDEN=${clust_iden} '$3>=IDEN {print $2 "\t" $1}' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk >  ${work_dir}/syn_pan_leftovers.clust.tsv

  echo "  Retrieve sequences for the leftover genes"
  get_fasta_from_family_file.pl ${fasta_files} \
    -fam ${work_dir}/syn_pan_leftovers.clust.tsv -out ${leftovers_dir}

  echo "  Make augmented cluster sets"
  make_augmented_cluster_sets leftovers_dir=${leftovers_dir} ${pan_fasta_dir}/* > ${work_dir}/syn_pan_augmented.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' ${work_dir}/syn_pan_augmented.clust.tsv > ${work_dir}/syn_pan_augmented.hsh.tsv
}
#
run_add_extra() {
  echo; echo "Add extra annotation sets to the augmented clusters, by homology"
  . "${PANDAGMA_CONF}"
  fasta_files="$(ls ${data_extra_dir}/*.${fasta_ext})"

  echo "  Search non-clustered genes against pan-gene consensus sequences"
  for path in ${data_extra_dir}/*.${fasta_ext}; do
    fasta_file=`basename $path` 
    echo "Extra: $fasta_file"

    mmseqs easy-search ${path} \
                       ${work_dir}/syn_pan_consen.fna \
                       ${work_extra_out_dir}/${fasta_file}.x.all_cons.m8 ${mmseqs_tmp_dir} \
                       --search-type 3 --cov-mode 5 -c 0.5 1>/dev/null
  done

  echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
  > ${work_dir}/syn_pan_extra.clust.tsv
  top_line.awk "${work_extra_out_dir}"/*.x.all_cons.m8 |
    awk -v IDEN=${clust_iden} '$3>=IDEN {print $2 "\t" $1}' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk >  ${work_dir}/syn_pan_extra.clust.tsv

  echo "  Retrieve sequences for the extra genes"
  get_fasta_from_family_file.pl ${fasta_files} \
    -fam ${work_dir}/syn_pan_extra.clust.tsv -out ${leftovers_extra_dir}

  echo "  Make augmented cluster sets"
  make_augmented_cluster_sets leftovers_dir=${leftovers_dir} ${pan_fasta_dir}/* > ${work_dir}/syn_pan_aug_extra.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' ${work_dir}/syn_pan_aug_extra.clust.tsv > ${work_dir}/syn_pan_aug_extra.hsh.tsv
}
#
run_summarize() {
  . "${PANDAGMA_CONF}"
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"
  full_out_dir=`echo "$out_dir_base.id${clust_iden}.cov${clust_cov}.I${mcl_inflation}" | perl -pe 's/(\d)\.(\d+)/$1_$2/g'`
  stats_file=${full_out_dir}/stats.txt

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  cp ${work_dir}/synteny_blocks.tsv ${full_out_dir}/
  cp ${work_dir}/syn_pan.clust.tsv ${full_out_dir}/
  cp ${work_dir}/syn_pan.hsh.tsv ${full_out_dir}/
  cp ${work_dir}/syn_pan_consen.fna ${full_out_dir}/
  cp ${work_dir}/syn_pan_aug*.tsv ${full_out_dir}/

  printf "Run of program $scriptname, version $version\n\n" > ${stats_file}

  printf "Parameter  \tvalue\n" >> ${stats_file}
  for key in ${pandagma_conf_params}; do
    if [[ ${key} != +(dag|mmseqs)_time_s ]]; then
      printf '%-15s\t%s\n' ${key} "${!key}" >> ${stats_file}
    fi
  done

  printf "\nOutput directory for this run:\t${full_out_dir}\n" >> ${stats_file}

  printf '%-20s\t%s\n' "Statistic" "value" >> ${stats_file}

  let "n_blocks=$(wc -l < ${full_out_dir}/synteny_blocks.tsv)-1"
  printf '%-20s\t%s\n' synteny_blocks $n_blocks >> ${stats_file}

#
  printf "\n== Initial clusters (containing only genes within synteny blocks)\n" >> ${stats_file}
  let "clusters=$(wc -l < ${full_out_dir}/syn_pan.clust.tsv)"
  printf '%-20s\t%s\n' "num_of_clusters" $clusters >> ${stats_file}

  let "largest=$(awk "{print NF-1}" ${full_out_dir}/syn_pan.clust.tsv | head -1)"
  printf '%-20s\t%s\n' "largest_cluster" $largest >> ${stats_file}

  let "mode=$(awk "{print NF-1}" ${full_out_dir}/syn_pan.clust.tsv | \
    uniq -c | sort -n | tail -1 | awk '{print $2}')"
  printf '%-20s\t%s\n' "modal_clst_size" $mode >> ${stats_file}

  let "num_at_mode=$(awk "{print NF-1}" ${full_out_dir}/syn_pan.clust.tsv | \
    uniq -c | sort -n | tail -1 | awk '{print $1}')"
  printf '%-20s\t%s\n' "num_at_mode" $num_at_mode >> ${stats_file}
  
  let "seqs_clustered=$(wc -l ${full_out_dir}/syn_pan.hsh.tsv | awk '{print $1}')"
  printf '%-20s\t%s\n' "seqs_clustered" $seqs_clustered >> ${stats_file}
  
#
  if [ -f ${full_out_dir}/syn_pan_augmented.clust.tsv ]; then
    printf "\n== Augmented clusters (unanchored sequences added to the initial clusters)\n" >> ${stats_file}
    let "clustersA=$(wc -l < ${full_out_dir}/syn_pan_augmented.clust.tsv)"
    printf '%-20s\t%s\n' "num_of_clusters" $clustersA >> ${stats_file}

    let "largestA=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_augmented.clust.tsv | sort -n | tail -1)"
    printf '%-20s\t%s\n' "largest_cluster" $largestA >> ${stats_file}

    let "modeA=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_augmented.clust.tsv | \
      uniq -c | sort -n | tail -1 | awk '{print $2}')"
    printf '%-20s\t%s\n' "modal_clst_size" $modeA >> ${stats_file}

    let "numA_at_mode=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_augmented.clust.tsv | \
      sort -n | uniq -c | sort -n | tail -1 | awk '{print $1}')"
    printf '%-20s\t%s\n' "num_at_mode" $numA_at_mode >> ${stats_file}
    
    let "seqs_clustered=$(wc -l ${full_out_dir}/syn_pan_augmented.hsh.tsv | awk '{print $1}')"
    printf '%-20s\t%s\n' "seqs_clustered" $seqs_clustered >> ${stats_file}
  fi

#
  if [ -f ${full_out_dir}/syn_pan_aug_extra.clust.tsv ]; then
    printf "\n== Augmented-extra clusters (sequences from extra annotation sets have been added)\n" >> ${stats_file}
    let "clustersB=$(wc -l < ${full_out_dir}/syn_pan_aug_extra.clust.tsv)"
    printf '%-20s\t%s\n' "num_of_clusters" $clustersB >> ${stats_file}

    let "largestB=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_aug_extra.clust.tsv | sort -n | tail -1)"
    printf '%-20s\t%s\n' "largest_cluster" $largestB >> ${stats_file}

    let "modeB=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_aug_extra.clust.tsv | \
      uniq -c | sort -n | tail -1 | awk '{print $2}')"
    printf '%-20s\t%s\n' "modal_clst_size" $modeB >> ${stats_file}

    let "numB_at_mode=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_aug_extra.clust.tsv | \
      sort -n | uniq -c | sort -n | tail -1 | awk '{print $1}')"
    printf '%-20s\t%s\n' "num_at_mode" $numB_at_mode >> ${stats_file}
    
    let "seqs_clustered=$(wc -l ${full_out_dir}/syn_pan_aug_extra.hsh.tsv | awk '{print $1}')"
    printf '%-20s\t%s\n' "seqs_clustered" $seqs_clustered >> ${stats_file}
  fi

  # histograms
  if [ -f ${full_out_dir}/syn_pan.clust.tsv ]; then
    printf "\nCounts of initial clusters by cluster size:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/syn_pan.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi

  if [ -f ${full_out_dir}/syn_pan_augmented.clust.tsv ]; then
    printf "\nCounts of augmented clusters by cluster size:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/syn_pan_augmented.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi

  if [ -f ${full_out_dir}/syn_pan_aug_extra.clust.tsv ]; then
    printf "\nCounts of augmented-extra clusters by cluster size:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/syn_pan_aug_extra.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi

  echo
  cat ${stats_file}
}
#
# top-level command functions
#
init() {
  # Initialize parameters required for run. Write these to files, for persistent access through the program.
  echo; echo "Setting run configuration parameters in ${PANDAGMA_CONF}"
  cat <<END > "${PANDAGMA_CONF}"
clust_iden='0.98'
clust_cov='0.75'
consen_iden='0.80'
mcl_inflation='2'
dagchainer_args='-g 10000 -M 50 -D 200000 -E 1e-5 -A 6 -s'
fasta_ext='fna'
gff_ext='gff3'
pan_prefix='pan'
out_dir_base='out'
END

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
                       adding sequences missed in the first clustering round.
           add_extra - Add other gene model sets to the primary clusters. Useful for adding
                       annotation sets that may be of lower or uncertain quality.
           summarize - compute synteny stats
"""
  commandlist="init ingest mmseqs filter dagchainer mcl consense add_extra summarize"
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
dirlist="work_dir mmseqs_dir mmseqs_tmp_dir dag_dir pan_fasta_dir \
         pan_consen_dir leftovers_dir work_data_dir \
         work_data_extra_dir leftovers_extra_dir work_extra_out_dir"
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
"clean")
  clean $@
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

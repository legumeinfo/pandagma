#!/usr/bin/env bash
#
# Configuration and run script for "pandagma", which generates pan-gene clusters using the programs 
# mmseqs, dagchainer, mcl, and vsearch. These are used to do the initial clustering, 
# synteny-finding, and re-clustering.
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2021
#
scriptname=`basename "$0"`
version="2021-11-02"
set -o errexit -o errtrace -o nounset -o pipefail

export NPROC=${NPROC:-1}
export MMSEQS_NUM_THREADS=${NPROC} # mmseqs otherwise uses all cores by default
export PANDAGMA_CONF=${PANDAGMA_CONF:-${PWD}/pandagma.conf}
export PANDAGMA_WORK_DIR=${PANDAGMA_WORK_DIR:-${PWD}/work}

# mmseqs uses a significant number of threads on its own. Set a maximum, which may be below NPROC.
MMSEQSTHREADS=$(( 10 < ${NPROC} ? 10 : ${NPROC} ))

pandagma_conf_params='clust_iden clust_cov consen_iden extra_iden mcl_inflation dagchainer_args pan_prefix out_dir_base consen_prefix extra_stats'

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
  run name_pangenes - Assign pan-gene names with consensus chromosomes and ordinal positions.
 run calc_chr_pairs - Report observed chromosome pairs; useful for preparing expected_chr_matches.tsv
      run summarize - Move results into output directory, and report summary statistics.

Variables in pandagma config file:
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.60]
        consen_iden - Minimum identity threshold for vsearch consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of \"extra\" annotations [90]
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]
         pan_prefix - Prefix to use for pangene clusters [pan]
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
}
cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}
# usage: make_augmented_cluster_sets leftovers_dir=${leftovers_dir} ${pan_fasta_dir}/* > syn_pan_aug.clust.tsv
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

### run functions ###

run_ingest() {
# Add positional information from GFF3 or 4- or 6-column BED to FASTA IDs
# BED start coordinate converted to 1-based
  cd "${PANDAGMA_WORK_DIR}"
  echo; echo "Run ingest: from fasta and gff or bed data, create fasta with IDs containing positional info."
  #
  # For use when calculating stats at the end, count genes in CDS files
    cat /dev/null > tmp.gene_count_start_0

  mkdir -p fasta posn_hsh
  for (( file_num = 0; file_num < ${#fasta_files[@]} ; file_num++ )); do
    file_base=$(basename ${fasta_files[file_num]%.*})
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | gff_or_bed_to_hash.awk > posn_hsh/$file_base.hsh
    hash_into_fasta_id.pl -fasta "${fasta_files[file_num]}" \
                          -hash posn_hsh/$file_base.hsh \
                          -out fasta/$file_base
    cat_or_zcat "${fasta_files[file_num]}" | 
      awk -v FILE=$file_base '$1~/^>/ {ct++} END{print FILE "\t" ct}' >> tmp.gene_count_start_0
  done

  # Also get position information from the "extra" annotation sets, if any.
  if (( ${#fasta_files_extra[@]} > 0 ))
  then
    for (( file_num = 0; file_num < ${#fasta_files_extra[@]} ; file_num++ )); do
      file_base=$(basename ${fasta_files_extra[file_num]%.*})
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra[file_num]}" | gff_or_bed_to_hash.awk > posn_hsh/$file_base.hsh
      cat_or_zcat "${fasta_files_extra[file_num]}" | 
        awk -v FILE=$file_base '$1~/^>/ {ct++} END{print FILE "\t" ct}' >> tmp.gene_count_start_0
    done
  fi

  # Prepare the tmp.gene_count_start to be joined, in run_summarize, with tmp.gene_count_end.
  # Depends on four-part prefix form of the names, e.g. glyma.Wm82.gnm2.ann2.BG1Q.cds_primary.fna
  cat tmp.gene_count_start_0 | perl -pe 's/^(\w+\.\w+\.\w+\.\w+)\.\S+\t(\d+)/$1\t$2/' | sort > tmp.gene_count_start
  rm tmp.gene_count_start_0
}

run_mmseqs() {
  # Do mmseqs clustering on all pairings of annotation sets.
  cd "${PANDAGMA_WORK_DIR}"
  echo; echo "Run mmseqs -- at ${clust_iden} percent identity and minimum of ${clust_cov}% coverage."
  #

  mkdir -p mmseqs mmseqs_tmp
  for qry_path in ${PANDAGMA_WORK_DIR}/fasta/*.fna; do
    qry_base=$(basename $qry_path .fna)
    for sbj_path in ${PANDAGMA_WORK_DIR}/fasta/*.fna; do
      sbj_base=$(basename $sbj_path .fna)
      if [[ "$qry_base" > "$sbj_base" ]]; then
        echo "Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
        { cat $qry_path $sbj_path ; } |
          mmseqs easy-cluster stdin mmseqs/${qry_base}.x.${sbj_base} mmseqs_tmp \
           --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null & # background
          # allow to execute up to $MMSEQSTHREADS in parallel
          if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
      fi
    done
  done
  wait # wait for last jobs to finish
}

run_filter() {
  echo; echo "From mmseqs cluster output, split out the following fields: molecule, gene, start, stop."
  cd "${PANDAGMA_WORK_DIR}"
  mkdir -p dag
  if [[ -f ${chr_match_list} ]]; then  # filter based on list of expected chromosome pairings if provided
    echo "Filtering on chromosome patterns from file ${chr_match_list}"
    for mmseqs_path in mmseqs/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      filter_mmseqs_by_chroms.pl -chr_pat ${chr_match_list} > dag/${outfilebase}_matches.tsv < ${mmseqs_path} &

      # allow to execute up to $NPROC in parallel
      if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
    done
    wait # wait for last jobs to finish
  else   # don't filter, since chromosome pairings aren't provided; just split lines on "__"
    echo "No expected_chr_matches.tsv file was provided, so proceeding without chromosome-pair filtering."
    for mmseqs_path in mmseqs/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path _cluster.tsv`
      perl -pe 's/__/\t/g' > dag/${outfilebase}_matches.tsv < "${mmseqs_path}"
    done
  fi
}

run_dagchainer() {
  # Identify syntenic blocks, using DAGchainer
  cd "${PANDAGMA_WORK_DIR}"
  echo; echo "Run DAGchainer, using args \"${dagchainer_args}\""
  for match_path in dag/*_matches.tsv; do
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

  awk '$1!~/^#/ {print $2 "\t" $6}' dag/*_matches.tsv > homology_pairs.tsv
  awk '$1!~/^#/ {print $2 "\t" $6}' dag/*.aligncoords > synteny_pairs.tsv

  #Extract single-linkage synteny anchors
  printf "matches\tscore\trev\tid1\tid2\n" > synteny_blocks.tsv
  for path in dag/*.aligncoords; do
    awk '/##/ && !/reverse/ {print substr($14,0,length($14)-2) "\t" $10 "\t" 1 "\t" $3 "\t" $5}' "${path}"
    awk '/##/ &&  /reverse/ {print substr($15,0,length($15)-2) "\t" $11 "\t" 1 "\t" $3 "\t" $5}' "${path}"
  done >> synteny_blocks.tsv
}

run_mcl() {
  # Calculate clusters using Markov clustering
  cd "${PANDAGMA_WORK_DIR}"
  printf "\nCalculate clusters. use Markov clustering with inflation parameter $mcl_inflation and ${NPROC} threads\n"
  echo "MCL COMMAND: mcl synteny_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv"
  mcl synteny_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv \
    1>/dev/null
 
  # Add cluster IDs
  awk -v PRE=${pan_prefix} '{padnum=sprintf("%05d", NR); print PRE padnum "\t" $0}' tmp.syn_pan.clust.tsv > syn_pan.clust.tsv

  # Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' syn_pan.clust.tsv > syn_pan.hsh.tsv
}

run_consense() {
  echo; echo "Identify a consensus sequence for each pan-gene set, using vsearch."
  echo "Then add previously unclustered sequences into an \"augmented\" pan-gene set, by homology."
  cd "${PANDAGMA_WORK_DIR}"
  mkdir -p pan_consen pan_cent pan_fasta

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  echo "    Fasta file:" "${fasta_files[@]}"
  get_fasta_from_family_file.pl "${fasta_files[@]}" -fam syn_pan.clust.tsv -out pan_fasta

  echo "  Calculate consensus and centroid sequences for each pan-gene set."
  ls pan_fasta | xargs -I{} -n 1 -P ${NPROC} \
    vsearch --cluster_fast pan_fasta/{} --id ${consen_iden} --fasta_width 0 \
            --consout pan_consen/{} \
            --centroid pan_cent/{} \
            --quiet \
            --threads 1

  echo "  Combine centroid sequences into one multifasta file."
  echo "  If there are multiple such sequences in a cluster, pick the centroid"
  echo "  sequence from the most frequent subcluster for the pangene."
 
  awk 'FNR==1 { nf=split(FILENAME, FN, "/") }
         /^>/ { print ">" FN[nf] " " substr($1,2); next }
              { print }' pan_consen/* | fasta_to_table.awk |
    perl -pe 's/^(\S+) centroid=([^;]+);seqs=(\d+)\t.+/$1\t$2\t$3/' | sort -k1,1 -k3nr,3nr |
    top_line.awk | perl -pe 's/^(\S+)\t(\S+)\t\d+/$1_$2/' > lis.syn_pan_consen_repr
 
  awk 'FNR==1 { nf=split(FILENAME, FN, "/") }
         /^>/ { print ">" FN[nf] "_" substr($1,2); next }
              { print }' pan_cent/* > syn_pan_cent_allTMP.fna
      
    get_fasta_subset.pl -in syn_pan_cent_allTMP.fna \
                        -out syn_pan_cent.fna \
                        -lis lis.syn_pan_consen_repr -clobber

    perl -pi -e 's/^(>\w+)_/$1 /' syn_pan_cent.fna

  # rm pan_consen/* &
  # rm pan_cent/* &

  echo "  Get sorted list of all genes, from the original fasta files"
  cat_or_zcat "${fasta_files[@]}" | awk '/^>/ {print substr($1,2)}' | sort > lis.all_genes

  echo "  Get sorted list of all clustered genes"
  awk '$1~/^>/ {print $1}' pan_fasta/* | sed 's/>//' | sort > lis.all_clustered_genes

  echo "  Get list of genes not in clusters"
  comm -13 lis.all_clustered_genes lis.all_genes > lis.genes_not_in_clusters

  echo "  Retrieve the non-clustered genes"
  cat_or_zcat "${fasta_files[@]}" |
    get_fasta_subset.pl -in /dev/stdin \
                        -out genes_not_in_clusters.fna \
                        -lis lis.genes_not_in_clusters -clobber

  echo "  Search non-clustered genes against pan-gene consensus sequences"
  mmseqs easy-search genes_not_in_clusters.fna \
                     syn_pan_cent.fna \
                     unclust.x.all_cons.m8 mmseqs_tmp \
                     --search-type 3 --cov-mode 5 -c 0.5 1>/dev/null 

  echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
  echo "  Use $extra_iden, which may be set more permissively than clust_iden, used for the initial clusters."
  top_line.awk unclust.x.all_cons.m8 | 
    awk -v IDEN=${extra_iden} '$3>=IDEN {print $2 "\t" $1}' |
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk >  syn_pan_leftovers.clust.tsv

  echo "  Retrieve sequences for the leftover genes"
  mkdir -p pan_leftovers
  get_fasta_from_family_file.pl "${fasta_files[@]}" \
    -fam syn_pan_leftovers.clust.tsv -out pan_leftovers/

  cd "${PANDAGMA_WORK_DIR}"
  echo "  Make augmented cluster sets"
  make_augmented_cluster_sets leftovers_dir=pan_leftovers pan_fasta/* > syn_pan_aug.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' syn_pan_aug.clust.tsv > syn_pan_augmented.hsh.tsv
}

run_add_extra() {
  if (( ${#fasta_files_extra[@]} > 0 ))
  then # handle the "extra" annotation files
    echo; echo "Add extra annotation sets to the augmented clusters, by homology"
    echo "  Search non-clustered genes against pan-gene consensus sequences"
    cd "${PANDAGMA_WORK_DIR}"
    mkdir -p extra_out_dir
    for path in "${fasta_files_extra[@]}"
    do
      fasta_file=`basename ${path%.*}` 
      echo "Extra: $fasta_file"
  
      mmseqs easy-search "${path}" \
                         syn_pan_cent.fna \
                         extra_out_dir/${fasta_file}.x.all_cons.m8 mmseqs_tmp/ \
                         --search-type 3 --cov-mode 5 -c 0.5 1>/dev/null & # background

      if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    wait # wait for jobs to finish
  
    echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
    top_line.awk extra_out_dir/*.x.all_cons.m8 |
      awk -v IDEN=${clust_iden} '$3>=IDEN {print $2 "\t" $1}' |
      sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > syn_pan_extra.clust.tsv
  
    echo "  Retrieve sequences for the syn_pan_augmented genes (from \"run consense\")"
    mkdir -p pan_augmented
    get_fasta_from_family_file.pl "${fasta_files[@]}" -fam syn_pan_aug.clust.tsv -out pan_augmented/
  
    echo "  Retrieve sequences for the extra genes"
    mkdir -p pan_leftovers_extra
    get_fasta_from_family_file.pl "${fasta_files_extra[@]}" -fam syn_pan_extra.clust.tsv -out pan_leftovers_extra/
  
    echo "  Make augmented cluster sets"
    make_augmented_cluster_sets leftovers_dir=pan_leftovers_extra pan_augmented/* > syn_pan_aug_extra.clust.tsv
  
    echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
    perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' syn_pan_aug_extra.clust.tsv > syn_pan_aug_extra.tmp
      # syn_pan_aug_extra.tmp is the main pangene hash file; but it is named "tmp" because 
      # it will be used later, in run_name_pangenes

    ########## 
    echo "  Merge fasta sets"
    mkdir -p pan_aug_leftover_merged
    for path in pan_augmented/*; do
      file=`basename $path`
      if [[ -f "pan_leftovers_extra/$file" ]]; then
        cat $path pan_leftovers_extra/$file > pan_aug_leftover_merged/$file
      else
        cp $path pan_aug_leftover_merged/
      fi
    done

    echo "  Calculate consensus and centroid sequences for each MERGED (augmented/extra) pan-gene set."
    mkdir pan_consen_merged pan_cent_merged
    ls pan_aug_leftover_merged | xargs -I{} -n 1 -P ${NPROC} \
      vsearch --cluster_fast pan_aug_leftover_merged/{} --id ${consen_iden} --fasta_width 0 \
              --consout pan_consen_merged/{} \
              --centroids pan_cent_merged/{} \
              --quiet \
              --threads 1
  
    echo "  Combine MERGED (augmented/extra) centroid sequences into one multifasta file."
    echo "  If there are multiple such sequences in a cluster, pick the centroid"
    echo "  sequence from the most frequent subcluster for the pangene."
  
    awk 'FNR==1 { nf=split(FILENAME, FN, "/") }
           /^>/ { print ">" FN[nf] " " substr($1,2); next }
                { print }' pan_consen_merged/* | fasta_to_table.awk |
      perl -pe 's/^(\S+) centroid=([^;]+);seqs=(\d+)\t.+/$1\t$2\t$3/' | sort -k1,1 -k3nr,3nr |
      top_line.awk | perl -pe 's/^(\S+)\t(\S+)\t\d+/$1_$2/' > lis.syn_pan_consen_merged_repr
  
    awk 'FNR==1 { nf=split(FILENAME, FN, "/") }
           /^>/ { print ">" FN[nf] "_" substr($1,2); next }
                { print }' pan_cent_merged/* > syn_pan_cent_merged_allTMP.fna
        
      get_fasta_subset.pl -in syn_pan_cent_merged_allTMP.fna \
                          -out syn_pan_cent_merged.fna \
                          -lis lis.syn_pan_consen_merged_repr -clobber

    perl -pi -e 's/^(>\w+)_/$1 /' syn_pan_cent_merged.fna

    # rm pan_consen_merged/* &
    # rm pan_cent_merged/* &

  else  # no "extra" fasta files, so just promote the syn_pan_aug files as syn_pan_aug_extra
        # TO DO: handle this better, i.e. report that no "extra" files were provided and skip these "aug" files.
    cp syn_pan_aug.clust.tsv syn_pan_aug_extra.clust.tsv
    cp syn_pan_augmented.hsh.tsv syn_pan_aug_extra.tmp
    cp syn_pan_cent.fna syn_pan_cent_merged.fna
fi
}

run_name_pangenes() {
  cd "${PANDAGMA_WORK_DIR}"

  echo "  Make a version of the hash output that contains positional information"
  join -a 1 -1 2 -2 1 <(sort -k2,2 syn_pan_aug_extra.tmp) <(cat posn_hsh/*hsh | sort -k1,1) | 
    perl -pe 's/__/\t/g; s/ /\t/g' | awk -v OFS="\t" '{print $2, $1, $3, $5, $6}' |
    sort -k1,1 -k2,2 > syn_pan_aug_extra.hsh.tsv

  echo "  Calculate consensus pan-gene positions"
  consen_pangene_posn.pl -pre ${consen_prefix}.chr syn_pan_aug_extra.hsh.tsv | 
    sort -k2,2 -k3n,3n | name_ordered_genes.awk > consen_${consen_prefix}.tsv

  echo "  Reshape defline into a hash, e.g. pan47789	Glycine.pan3.chr01__Glycine.pan3.chr01_000100__45224__45786"
  cat consen_${consen_prefix}.tsv | 
    perl -pe 's/^(\S+)\t([^.]+\.[^.]+)\.(chr\d+)_(\d+)\t(\d+)\t(\d+)/$1\t$2.$3__$2.$3_$4__$5__$6/' > consen_posn.hsh

  echo "  Hash position information into fasta file"
  hash_into_fasta_id.pl -fasta syn_pan_cent_merged.fna -hash consen_posn.hsh -out_file syn_pan_cent_merged_posnTMP.fna

  echo "  Use mmseqs to cluster fasta file, to identify near-duplicate neighboring paralogs"
  # Note: this uses $consen_iden rather than $clust_iden, with $consen_iden potentially more lenient.
  mmseqs easy-cluster syn_pan_cent_merged_posnTMP.fna  syn_pan_cent_merged_posn_mmseqs mmseqs_tmp \
    --min-seq-id $consen_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null 

  # Parse mmseqs clusters into pairs of genes that are similar and ordinally close
  close=10 # genes within a neighborhood of +-close genes among the consensus-orderd genes
  cat syn_pan_cent_merged_posn_mmseqs_cluster.tsv | perl -pe 's/__/\t/g; s/(chr\d+)_(\d+)/$1_$2\t$2/g' |
    awk -v CLOSE=$close '$1==$6 && $2!=$7 && sqrt(($3/100-$8/100)^2)<=CLOSE {print $2 "\t" $7}' |
    cat > syn_pan_cent_merged_posn_mmseqs.closepairs

  # Cluster the potential near-duplicate neighboring paralogs
  mcl  syn_pan_cent_merged_posn_mmseqs.closepairs --abc -o syn_pan_cent_merged_posn_mmseqs.clst
 
  # Keep the first gene from the cluster and discard the rest. Do this by removing the others.
  cat syn_pan_cent_merged_posn_mmseqs.clst | awk '{$1=""}1' | awk '{$1=$1}1' | tr ' ' '\n' > lis.clust_genes_remove

  # Reshape defline
  cat syn_pan_cent_merged_posnTMP.fna | sed 's/__/\t/g' |
    awk '$1~/^>/ {print ">" $2 " " $3 " " $4} $1!~/^>/ {print}' > syn_pan_cent_merged_posn.fna

  # Remove neighboring close paralogs, leaving one
  get_fasta_subset.pl -in syn_pan_cent_merged_posn.fna -out syn_pan_cent_merged_posn_trim.fna \
    -lis lis.clust_genes_remove -xclude -clobber

}

run_calc_chr_pairs() {
  cd "${PANDAGMA_WORK_DIR}"
  echo "  Generate a report of observed chromosome pairs"

  cat syn_pan_aug_extra.hsh.tsv | 
    awk 'BEGIN{IGNORECASE=1} $3!~/cont|scaff|sc|pilon|mito|chl|unanchor/ {print $1 "\t" $3}' |
    perl -pe 's/(^\S+)\t.+\.\D+(\d+)/$1\t$2/' | sort | uniq | pile_against_col1.awk |
    perl -pe 's/^\S+\t//' | sort | uniq -c | perl -pe 's/^ +(\d+)\s/$1\t/' | sort -k1n,1n |
    awk -v OFS="\t" 'NF==4 {print $1, $2, $3 "\n" $1, $2, $4 "\n" $1, $3, $4}
                   NF==5 {print $1, $2, $3 "\n" $1, $2, $4 "\n" $1, $2, $5;
                          print $1, $3, $4 "\n" $1, $3, $5; print $1, $4, $5}
                   NF<4 || NF>5 {print}' |
    awk -v OFS="\t" 'NF==2 {print $1, $2, $2} NF>2 {print}' | sort -k2n,2n -k3n,3n |
    awk -v OFS="\t" 'NR==1 {sum=$1; prev2=$2; prev3=$3}
                     NR>1 && prev2==$2 && prev3==$3 {sum+=$1; prev2=$2; prev3=$3}
                     NR>1 && prev2!=$2 || prev3!=$3 {print sum, prev2, prev3; prev2=$2; prev3=$3; sum=$1}
                     END{print sum, prev2, prev3}' |
    sort -k1nr,1nr -k2n,2n -k3n,3n | sed '/^[[:space:]]*$/d' > observed_chr_pairs.tsv
}

run_summarize() {
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"
  full_out_dir=`echo "$out_dir_base.id${clust_iden}.cov${clust_cov}.cns${consen_iden}.ext${extra_iden}.I${mcl_inflation}" |
                       perl -pe 's/0\.(\d+)/$1/g'`
  stats_file=${full_out_dir}/stats.txt

  cd "${submit_dir}"

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p $full_out_dir
  fi

  cp ${PANDAGMA_WORK_DIR}/synteny_blocks.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/syn_pan.clust.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/syn_pan.hsh.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/syn_pan_aug*.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/observed_chr_pairs.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/consen_${consen_prefix}.tsv ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/syn_pan_cent_merged_posn.fna ${full_out_dir}/
  cp ${PANDAGMA_WORK_DIR}/syn_pan_cent_merged_posn_trim.fna ${full_out_dir}/

  printf "Run of program $scriptname, version $version\n\n" > ${stats_file}

  printf "Parameter  \tvalue\n" >> ${stats_file}
  for key in ${pandagma_conf_params}; do
    printf '%-15s\t%s\n' ${key} "${!key}" >> ${stats_file}
  done

  printf "\nOutput directory for this run:\t${full_out_dir}\n" >> ${stats_file}

  printf '%-20s\t%s\n' "Statistic" "value" >> ${stats_file}

  let "n_blocks=$(wc -l < ${full_out_dir}/synteny_blocks.tsv)-1"
  printf '%-20s\t%s\n' synteny_blocks $n_blocks >> ${stats_file}

#
  printf "\n== Initial clusters (containing only genes within synteny blocks)\n" >> ${stats_file}
  let "clusters=$(wc -l < ${full_out_dir}/syn_pan.clust.tsv)"
  printf '%-20s\t%s\n' "num_of_clusters" $clusters >> ${stats_file}

  let "largest=$(awk 'NR == 1 {print NF-1; exit}' ${full_out_dir}/syn_pan.clust.tsv)"
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
  if [ -f ${full_out_dir}/syn_pan_aug.clust.tsv ]; then
    printf "\n== Augmented clusters (unanchored sequences added to the initial clusters)\n" >> ${stats_file}
    let "clustersA=$(wc -l < ${full_out_dir}/syn_pan_aug.clust.tsv)"
    printf '%-20s\t%s\n' "num_of_clusters" $clustersA >> ${stats_file}

    let "largestA=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_aug.clust.tsv | sort -n | tail -1)"
    printf '%-20s\t%s\n' "largest_cluster" $largestA >> ${stats_file}

    let "modeA=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_aug.clust.tsv | \
      uniq -c | sort -n | tail -1 | awk '{print $2}')"
    printf '%-20s\t%s\n' "modal_clst_size" $modeA >> ${stats_file}

    let "numA_at_mode=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_aug.clust.tsv | \
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
      sort -n | uniq -c | sort -n | tail -1 | awk '{print $2}')"
    printf '%-20s\t%s\n' "modal_clst_size" $modeB >> ${stats_file}

    let "numB_at_mode=$(awk "{print NF-1}" ${full_out_dir}/syn_pan_aug_extra.clust.tsv | \
      sort -n | uniq -c | sort -n | tail -1 | awk '{print $1}')"
    printf '%-20s\t%s\n' "num_at_mode" $numB_at_mode >> ${stats_file}
    
    let "seqs_clustered=$(wc -l ${full_out_dir}/syn_pan_aug_extra.hsh.tsv | awk '{print $1}')"
    printf '%-20s\t%s\n' "seqs_clustered" $seqs_clustered >> ${stats_file}
  fi

  # Print per-annotation-set coverage stats (sequence counts, sequences retained), if stats-extra flag is set
  if [ ${extra_stats,,} = "yes" ]; then
    cut -f2 ${PANDAGMA_WORK_DIR}/syn_pan_aug_extra.hsh.tsv | perl -pe 's/(\w+\.\w+\.\w+\.\w+)\..+/$1/' | 
      sort | uniq -c | awk '{print $2 "\t" $1}' > ${PANDAGMA_WORK_DIR}/tmp.gene_count_end

    # tmp.gene_count_start was generated during run_ingest

    printf "\n== Start and end gene counts per annotation set\n" >> ${stats_file}

    join ${PANDAGMA_WORK_DIR}/tmp.gene_count_start ${PANDAGMA_WORK_DIR}/tmp.gene_count_end | 
      awk 'BEGIN{print "\tStart\tEnd\tKept\tAnnotation_name"} 
          { printf "\t%i\t%i\t%2.1f\t%s\n", $2, $3, 100*($3/$2), $1 }'  >> ${stats_file}
  fi

  # histograms
  if [ -f ${full_out_dir}/syn_pan.clust.tsv ]; then
    printf "\nCounts of initial clusters by cluster size, file syn_pan.clust.tsv:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/syn_pan.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi

  if [ -f ${full_out_dir}/syn_pan_aug.clust.tsv ]; then
    printf "\nCounts of augmented clusters by cluster size, file syn_pan_aug.clust.tsv:\n" >> ${stats_file}
    awk '{print NF-1}' ${full_out_dir}/syn_pan_aug.clust.tsv |
      sort | uniq -c | awk '{print $2 "\t" $1}' | sort -n >> ${stats_file}
  fi

  if [ -f ${full_out_dir}/syn_pan_aug_extra.clust.tsv ]; then
    printf "\nCounts of augmented-extra clusters by cluster size, file syn_pan_aug_extra.clust.tsv:\n" >> ${stats_file}
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
clust_iden='0.95'
clust_cov='0.60'
consen_iden='0.80'
extra_iden='0.90'
mcl_inflation='2'
dagchainer_args='-g 10000 -M 50 -D 200000 -E 1e-5 -A 6 -s'
pan_prefix='pan'
out_dir_base='out'
consen_prefix='Genus.pan1'
extra_stats='no'

##### (required) list of GFF3/4-column BED & FASTA file paths
# Uncomment add file paths to the the annotation_files and fasta_files arrays.
# The nth listed annotation file corresponds to the nth listed FASTA file.
# Files with a ".gz" suffix are uncompressed on-the-fly; otherwise, file suffix
# is ignored.

#annotation_files=(
# file1.gff.gz
# file2.bed.gz
#)

#fasta_files=(
# file1.fna.gz
# file2.fa.gz
#)

#### (optional) Extra GFF3/BED & FASTA files
#annotation_files_extra=(
# file1-extra.gff3.gz
#)

#fasta_files_extra=(
# file1-extra.fa.gz
#)

##### (optional) expected_chr_matches file path
expected_chr_matches=''
END

  echo "Work directory (for temporary files): ${PANDAGMA_WORK_DIR}"; echo
}
#
run() {
  RUN_DOC="""Run an analysis step

Usage:
   $scriptname run [STEP]

Steps:
   If STEP is not set, the following steps will be run in order,
   otherwise the step is run by itself:
                init - Initialize parameters required for run
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
  echo >&2 "ERROR -- command \"$command\" not recognized."
  exit 1
  ;;
esac

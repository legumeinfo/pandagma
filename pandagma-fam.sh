#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates gene-family 
# orthogroups using the programs mmseqs, dagchainer, and mcl. 
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2023
#
scriptname=`basename "$0"`
version="2023-08-08"
set -o errexit -o errtrace -o nounset -o pipefail

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

HELP_DOC="""Compute orthogroups using a combination of synteny and homology,
using the programs mmseqs, dagchainer, and mcl, and additional pre- and post-refinement steps.

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
    mmseqs dagchainer mcl

Also, please add the pandagma utility programs in the bin directory adjacent to pandagma-fam.sh, e.g.
    PATH=$PWD/bin:\$PATH

Primary protein sequences and annotation (GFF3 or BED) files must be listed in the
config file, in the arrays annotation_files and protein_files. See example files.

Subcommands (in order they are usually run):
                all - All of the steps below, except for clean and ReallyClean
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
              align - Align families.
     model_and_trim - Build HMMs and trim the alignments, preparatory to calculating trees.
         calc_trees - Calculate gene trees.
          summarize - Move results into output directory, and report summary statistics.

  Run either of the following subcommands separately if you wish:
              clean - Clean (delete) files in the working directory that are not needed 
                        for later addition of data using add_extra and subsequent run commands.
                        By default, \"clean\" is run as part of \"all\" unless the -r flag is set.
        ReallyClean - Do complete clean-up of files in the working directory.
                        Use this if you want to start over, OR if you are satisified with the results and
                        don't anticipate adding other annotation sets to this pan-gene set.
"""

MORE_INFO="""
Optionally, a file specified in the expected_quotas variable can be specified in pandagma.conf,
which provides the number of expected paralogs for a species in a gene family at the desired evolutionary depth,
considering a known history of whole-genome duplications that affect the included species. For example,
for the legume family, originating ~70 mya:
  Arachis     4
  Cercis      1
  Glycine     4
  Phaseolus   2

These quotas are used in a regular expression to identify the initial portion of the prefixed seqIDs,
for example, \"Cicer\" and \"cerca\" matching the genus and \"gensp\" matching prefixes for these seqIDs:
  Cicer.pan1.chr06
  cerca.ISC453364.gnm3.Chr03

If an expected_quotas.tsv file is not provided, then no such filtering will be done.

Variables in pandagma config file (Set the config with the CONF environment variable)
    dagchainer_args - Argument for DAGchainer command
         clust_iden - Minimum identity threshold for mmseqs clustering [0.95]
          clust_cov - Minimum coverage for mmseqs clustering [0.60]
        consen_iden - Minimum identity threshold for consensus generation [0.80]
         extra_iden - Minimum identity threshold for mmseqs addition of \"extra\" annotations [90]
      mcl_inflation - Inflation parameter, for Markov clustering [default: 1.2]
        strict_synt - For clustering of the \"main\" annotations, use only syntenic pairs (1)
                        The alternative (0) is to use all homologous pairs that satisfy expected_quotas.tsv
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
"""

########################################
# Helper functions begin here

version() {
  echo $scriptname $version
}

##########
canonicalize_paths() {
  echo "Entering canonicalize_paths. Annotation files: ${annotation_files[@]}"

  annotation_files=($(realpath --canonicalize-existing "${annotation_files[@]}"))
  protein_files=($(realpath --canonicalize-existing "${protein_files[@]}"))
  cds_files=($(realpath --canonicalize-existing "${cds_files[@]}"))

  if (( ${#protein_files_extra[@]} > 0 ))
  then
    protein_files_extra=($(realpath --canonicalize-existing "${protein_files_extra[@]}"))
    annotation_files_extra=($(realpath --canonicalize-existing "${annotation_files_extra[@]}"))
  fi

  if (( ${#cds_files_extra[@]} > 0 ))
  then
    cds_files_extra=($(realpath --canonicalize-existing "${cds_files_extra[@]}"))
    #annotation_files_extra=($(realpath --canonicalize-existing "${annotation_files_extra[@]}"))
  fi

  readonly expected_quotas=${expected_quotas:+$(realpath "${expected_quotas}")}
  readonly submit_dir=${PWD}

  fasta_file=$(basename "${protein_files[0]}" .gz)
  faa="${fasta_file##*.}"

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
  
  mkdir -p 02_fasta_nuc 02_fasta_prot 01_posn_hsh stats

    # Prepare the tmp.gene_count_start to be joined, in run_summarize, with tmp.gene_count_end_pctl??_end.
    # This is captured from the gene IDs using the annot_str_regex set in the config file.
    cat /dev/null > stats/tmp.gene_count_start
    cat /dev/null > stats/tmp.fasta_seqstats
    start_time=`date`
    printf "Run started at: $start_time\n" > stats/tmp.timing

  export ANN_REX=${annot_str_regex}

  echo "  Get position information from the main annotation sets (protein)."
  cat /dev/null > 02_all_main_prot.faa # Collect all starting protein sequences, for later comparisons
  for (( file_num = 0; file_num < ${#protein_files[@]} ; file_num++ )); do
    file_base=$(basename ${protein_files[file_num]%.*})
    cat_or_zcat "${protein_files[file_num]}" >> 02_all_main_prot.faa # Collect original seqs for later comparisons
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | 
      gff_or_bed_to_hash5.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${protein_files[file_num]}" \
                          -hash 01_posn_hsh/$file_base.hsh \
                          -out 02_fasta_prot/$file_base
    # calc basic sequence stats
    annot_name=$(basename 02_fasta_prot/$file_base | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
    printf "  Main:  " >> stats/tmp.fasta_seqstats
    cat_or_zcat "${protein_files[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
  done

  echo "  Get position information from the main annotation sets (cds)."
  cat /dev/null > 02_all_main_cds.fna # Collect all starting cds sequences, for later comparisons
  for (( file_num = 0; file_num < ${#cds_files[@]} ; file_num++ )); do
    file_base=$(basename ${cds_files[file_num]%.*})
    cat_or_zcat "${cds_files[file_num]}" >> 02_all_main_cds.fna # Collect original seqs for later comparisons
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | 
      gff_or_bed_to_hash5.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files[file_num]}" \
                          -hash 01_posn_hsh/$file_base.hsh \
                          -out 02_fasta_nuc/$file_base
  done

  echo "  Get position information from the extra annotation sets (protein), if any."
  if (( ${#protein_files_extra[@]} > 0 ))
  then
    cat /dev/null > 02_all_extra_protein.faa # Collect all starting sequences, for later comparisons
    for (( file_num = 0; file_num < ${#protein_files_extra[@]} ; file_num++ )); do
      file_base=$(basename ${protein_files_extra[file_num]%.*})
      cat_or_zcat "${protein_files_extra[file_num]}" >> 02_all_extra_protein.faa  # Collect original seqs for later comparisons
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra[file_num]}" | 
        gff_or_bed_to_hash5.awk > 01_posn_hsh/$file_base.hsh
      hash_into_fasta_id.pl -nodef -fasta "${protein_files_extra[file_num]}" \
                            -hash 01_posn_hsh/$file_base.hsh \
                            -out 02_fasta_prot/$file_base
      # calc basic sequence stats
      annot_name=$(basename 02_fasta_prot/$file_base | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
      printf "  Extra: " >> stats/tmp.fasta_seqstats
      cat_or_zcat "${protein_files_extra[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
    done
  fi

  echo "  Get position information from the extra annotation sets (cds), if any."
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
    done
  fi

  echo "  Count starting sequences, for later comparisons"
  for file in 02_fasta_prot/*.$faa; do
    awk '$0~/UNDEFINED/ {ct++} 
      END{if (ct>0){print "Warning: " FILENAME " has " ct " genes without position (HASH UNDEFINED)" } }' $file
    cat $file | grep '>' | perl -pe 's/__/\t/g' | cut -f2 | # extracts the GeneName from combined genspchr__GeneName__start__end__orient
      perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex.+/$1/' |
      grep -v UNDEFINED | sort | uniq -c | awk '{print $2 "\t" $1}' >> stats/tmp.gene_count_start
  done

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
  SEQTYPE=1;  # 1=pep; 3=nuc
  for (( file1_num = 0; file1_num < ${#protein_files[@]} ; file1_num++ )); do
    qry_base=$(basename ${protein_files[file1_num]%.*} .$faa)
    for (( file2_num = $file1_num; file2_num < ${#protein_files[@]} ; file2_num++ )); do  # file2_num = $file1_num includes self-comparisons
      sbj_base=$(basename ${protein_files[file2_num]%.*} .$faa)
      echo "  Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)

      if [[ ${qry_base} == ${sbj_base} ]]; then # self-comparison, so use flag --add-self-matches
        mmseqs easy-search \
          02_fasta_prot/$qry_base.$faa 02_fasta_prot/$sbj_base.$faa 03_mmseqs/${qry_base}.x.${sbj_base}.m8 $MMTEMP \
          --add-self-matches --search-type ${SEQTYPE} --cov-mode 0 -c ${clust_cov} --min-seq-id ${clust_iden} 1>/dev/null 
      else # not a self-comparison, do omit flag --add-self-matches
        mmseqs easy-search \
          02_fasta_prot/$qry_base.$faa 02_fasta_prot/$sbj_base.$faa 03_mmseqs/${qry_base}.x.${sbj_base}.m8 $MMTEMP \
          --search-type ${SEQTYPE} --cov-mode 0 -c ${clust_cov} --min-seq-id ${clust_iden} 1>/dev/null 
      fi
    done
    echo
  done
  wait # wait for last jobs to finish
}

##########
run_filter() {
  echo; echo "From homology (mmseqs -m8) output, split out the following fields: molecule, gene, start, stop."
  cd "${WORK_DIR}"
  if [ -d 04_dag ]; then rm -rf 04_dag ; fi
  mkdir -p 04_dag
  if [[ -f ${expected_quotas} ]]; then  # filter based on list of match quotas if provided
    echo "Filtering on quotas from file ${expected_quotas}"
    for mmseqs_path in 03_mmseqs/*.m8; do
      outfilebase=`basename $mmseqs_path .m8`
      echo "  Filtering $outfilebase based on expected quotas"
      cat ${mmseqs_path} | filter_mmseqs_by_quotas.pl -quotas ${expected_quotas} |
        perl -pe 's/\t[\+-]//g' |  # strip orientation, which isn't used by DAGChainer 
        cat | # Next: for self-comparisons, suppress same chromosome && different gene ID (local dups)
        perl -lane 'print $_ unless($F[0] eq $F[4] && $F[1] ne $F[5])' |
        awk 'NF==8' |  # matches for genes with coordinates. The case of <8 can happen for seqs with UNDEF position.
        cat > 04_dag/${outfilebase}_matches.tsv &

      # allow to execute up to $NPROC in parallel
      if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
    done
    wait # wait for last jobs to finish
  else   # don't filter, since quotas file isn't provided; just remove orientation, which isn't used by DAGChainer
    echo "No expected_quotas.tsv file was provided, so proceeding without quota (expected gene-count) filtering."
    for mmseqs_path in 03_mmseqs/*_cluster.tsv; do
      outfilebase=`basename $mmseqs_path .m8`
      cat ${mmseqs_path} | cut -f1,2 | 
        perl -pe 's/__/\t/g; s/\t[\+-]//g' | 
        cat | # Next: or self-comparisons, suppress same chromosome && different gene ID (local dups)
        perl -lane 'print $_ unless($F[0] eq $F[4] && $F[1] ne $F[5])' |
        awk 'NF==8' |  # matches for genes with coordinates. The case of <8 can happen for seqs with UNDEF position.
        cat > 04_dag/${outfilebase}_matches.tsv 
    done
  fi
}

##########
run_dagchainer() {
  # Identify syntenic blocks, using DAGchainer
  cd "${WORK_DIR}"
  echo; echo "Run DAGchainer, using args \"${dagchainer_args}\""
  # Check and preemptively remove malformed \*_matches.file, which can result from an aborted run
  if [ -f 04_dag/\*_matches.tsv ]; then rm 04_dag/\*_matches.tsv; fi
  for match_path in 04_dag/*_matches.tsv; do
    #echo "basename $match_path _matches.tsv"
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
    if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
  done
  wait # wait for last jobs to finish
}

##########
run_mcl() {
  # Calculate clusters using Markov clustering
  cd "${WORK_DIR}"
  mkdir -p lists

  printf "\nPreparatory to clustering, combine the homology data into a file with gene pairs to be clustered.\n"

  if [ "$strict_synt" -eq 1 ]; then
    # https://stackoverflow.com/questions/3601515/how-to-check-if-a-variable-is-set-in-bash
    # https://pubs.opengroup.org/onlinepubs/9699919799/utilities/V3_chap02.html#tag_18_06_02
    if [ -z ${ks_block_wgd_cutoff+x} ] || [ -z ${ks_pair_cutoff+x} ] || 
        [ ! -d "${WORK_DIR}/05_kaksout" ] && [ ! "$(ls -A 05_kaksout/*.rptout)" ]; then
      echo "## One or both of ks_block_wgd_cutoff ks_pair_cutoff are unset or 05_kaksout doesn't exist; no Ks filtering will be done.";
      echo "## Combine the DAGChainer synteny pairs into a file to be clustered."
      cat 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' | awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
    else 
      echo "## The ks_block_wgd_cutoff is set to '$ks_block_wgd_cutoff' and ks_pair_cutoff is set to '$ks_pair_cutoff'";
      echo "## Combine the DAGChainer synteny pairs, with additional filtering by pairwise and block Ks thresholds, into a file to be clustered."
      # Combine the synteny pairs into a file to be clustered
      cat 05_kaksout/*.rptout | 
        awk -v PAIR_CUTOFF=$ks_pair_cutoff -v BLOCK_CUTOFF=$ks_block_wgd_cutoff '
            NF>0 && $1!~/^#/ && ($5<=PAIR_CUTOFF || $7<=BLOCK_CUTOFF) {print $1 "\t" $2}' |
          sort -u > 05_filtered_pairs.tsv
    fi
  else 
    echo "## Combine the homology pairs (filtered by quota if provided) into a file to be clustered."
    cat 04_dag/*_matches.tsv 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' | 
      awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
  fi

  printf "\nDo Markov clustering with inflation parameter $mcl_inflation and ${NPROC} threads\n"
  echo "MCL COMMAND: mcl 05_filtered_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv"
  mcl 05_filtered_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv \
    1>/dev/null
 
  echo "  Add cluster IDs"
  awk -v PRE=$consen_prefix '
    {padnum=sprintf("%05d", NR); print PRE padnum "\t" $0}' tmp.syn_pan.clust.tsv > 06_syn_pan.clust.tsv
  rm tmp.syn_pan.clust.tsv

  echo "  Move singleton and doubleton clusters into a list of leftovers"
  cat 06_syn_pan.clust.tsv | awk 'NF==2 {print $2} NF==3 {print $2; print $3}' > lists/lis.06_syn_pan_1s_2s
  cat 06_syn_pan.clust.tsv | awk 'NF>3 {print $0}' > 06_syn_pan_ge3.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 06_syn_pan_ge3.clust.tsv > 06_syn_pan_ge3.hsh.tsv
}

##########
run_consense() {
  echo; 
  echo "Add previously unclustered sequences into an \"augmented\" pan-gene set, by homology."
  cd "${WORK_DIR}"
  if [ -d 07_pan_fasta ]; then rm -rf 07_pan_fasta; fi
  mkdir -p 07_pan_fasta lists


  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  echo "    Fasta file:" "${protein_files[@]}"
  get_fasta_from_family_file.pl "${protein_files[@]}" -fam 06_syn_pan_ge3.clust.tsv -out 07_pan_fasta

  cat /dev/null > 07_pan_fasta_prot.faa
  for path in 07_pan_fasta/*; do
    pan_file=`basename $path`
    cat $path | awk -v panID=$pan_file ' $1~/^>/ {print ">" panID "__" substr($0,2) }
                      $1!~/^>/ {print $1} ' >> 07_pan_fasta_prot.faa 
  done

  echo "  Get sorted list of all genes, from the original fasta files"
  cat_or_zcat "${protein_files[@]}" | awk '/^>/ {print substr($1,2)}' | sort > lists/09_all_genes

  echo "  Get sorted list of all clustered genes"
  awk '$1~/^>/ {print $1}' 07_pan_fasta/* | sed 's/>//' | sort > lists/09_all_clustered_genes

  echo "  Get list of genes not in clusters"
  comm -13 lists/09_all_clustered_genes lists/09_all_genes > lists/09_genes_not_in_clusters

  echo "  Retrieve the non-clustered genes"
  cat_or_zcat "${protein_files[@]}" |
    get_fasta_subset.pl -in /dev/stdin -clobber -lis lists/09_genes_not_in_clusters \
                        -out 09_genes_not_in_clusters.faa 

  echo "  Search non-clustered genes against genes already clustered."

  # Check sequence type (in case this run function is called separately from the usually-prior ones)
  someseq=$(head 07_pan_fasta_prot.faa | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
  SEQTYPE=$(check_seq_type "${someseq}") # 3=nuc; 1=pep
  echo "SEQTYPE is: $SEQTYPE"

  mmseqs easy-search 09_genes_not_in_clusters.faa \
                     07_pan_fasta_prot.faa \
                     10_unclust.x.07_pan_fasta.m8 \
                     03_mmseqs_tmp \
                     --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null 

  echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
  echo "  Use the \"main set\" $clust_iden threshold."
  export FAM_PRE=${consen_prefix}
  top_line.awk 10_unclust.x.07_pan_fasta.m8 | 
    awk -v IDEN=${clust_iden} '$3>=IDEN {print $2 "\t" $1}' | 
    perl -pe '$prefix=$ENV{FAM_PRE}; s/^($prefix\d+)__\S+/$1/' |  # trims gene ID from preceding family ID
    sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk >  11_syn_pan_leftovers.clust.tsv

  echo "  Retrieve sequences for the leftover genes"
  mkdir -p 11_pan_leftovers
  get_fasta_from_family_file.pl "${protein_files[@]}" \
    -fam 11_syn_pan_leftovers.clust.tsv -out 11_pan_leftovers/

  echo "  Make augmented cluster sets"
  cat /dev/null > 12_syn_pan_aug_pre.clust.tsv
  augment_cluster_sets.awk leftovers_dir=11_pan_leftovers 07_pan_fasta/* > 12_syn_pan_aug_pre.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 12_syn_pan_aug_pre.clust.tsv > 12_syn_pan_aug_pre.hsh.tsv
}

##########
run_cluster_rest() {
  echo
  echo; echo "== Retrieve unclustered sequences and cluster those that can be =="
  cd "${WORK_DIR}"

  echo "  Retrieve genes present in the original protein files but absent from 12_syn_pan_aug.hsh"
  cut -f2 12_syn_pan_aug_pre.hsh.tsv | LC_ALL=C sort > lists/lis.12_syn_pan_aug_complement
  get_fasta_subset.pl -in 02_all_main_prot.faa -out 12_syn_pan_aug_complement.faa \
    -lis lists/lis.12_syn_pan_aug_complement -xclude -clobber

  MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
  mkdir -p 03_mmseqs_rest
  complmt_self_compare="12_syn_pan_aug_complement.x.12_syn_pan_aug_complement"
  cat 12_syn_pan_aug_complement.faa |
    mmseqs easy-cluster stdin 03_mmseqs_rest/$complmt_self_compare $MMTEMP \
    --min-seq-id $clust_iden -c $clust_cov --cov-mode 0 --cluster-reassign 1>/dev/null

  echo "  Cluster the remaining sequences that have matches"
  mcl 03_mmseqs_rest/${complmt_self_compare}_cluster.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan_aug_complement.clust.tsv

  echo "  Find number of clusters in initial (06) results"
  clust_count_06=`wc -l 06_syn_pan_ge3.clust.tsv | awk '{print $1}'`
 
  echo "  Add cluster IDs"
  cat tmp.syn_pan_aug_complement.clust.tsv | 
    awk -v START=$clust_count_06 -v PRE=$consen_prefix '
      NF>1 {padnum=sprintf("%05d", NR+START); print PRE padnum "\t" $0}' |
    cat > 12_syn_pan_aug_complement.clust.tsv
  rm tmp.syn_pan_aug_complement.clust.tsv

  echo "  Combine clusters derived from synteny or restricted homology (12_syn_pan_aug.clust.tsv)"
  echo "  and clusters from the complement of that set. Skip singletons (\"clusters\" of size 1)."
  cat 12_syn_pan_aug_pre.clust.tsv 12_syn_pan_aug_complement.clust.tsv | awk 'NF>1' > 12_syn_pan_aug.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 12_syn_pan_aug.clust.tsv > 12_syn_pan_aug.hsh.tsv
}

##########
run_add_extra() {
  echo; echo "== Add extra annotation sets (if provided) to the augmented clusters, by homology =="
  cd "${WORK_DIR}"

  if [ -d 13_extra_out_dir ]; then rm -rf 13_extra_out_dir; fi
  if [ -d 13_pan_aug_fasta ]; then rm -rf 13_pan_aug_fasta; fi
  mkdir -p 13_extra_out_dir 13_pan_aug_fasta

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  get_fasta_from_family_file.pl "${protein_files[@]}" -fam 12_syn_pan_aug.clust.tsv -out 13_pan_aug_fasta
  
  cat /dev/null > 13_pan_aug_fasta.faa
  for path in 13_pan_aug_fasta/*; do
    pan_file=`basename $path`
    cat $path | awk -v panID=$pan_file ' $1~/^>/ {print ">" panID "__" substr($0,2) }
                      $1!~/^>/ {print $1} ' >> 13_pan_aug_fasta.faa 
  done

  if (( ${#protein_files_extra[@]} > 0 ))
  then # handle the "extra" annotation files
    echo "  Search non-clustered genes against pan-gene consensus sequences"
    # Check sequence type (in case this run function is called separately from the usually-prior ones)
    someseq=$(head 07_pan_fasta_prot.faa | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
    SEQTYPE=$(check_seq_type "${someseq}") # 3=nuc; 1=pep
    echo "SEQTYPE is: $SEQTYPE"

    for path in "${protein_files_extra[@]}"; do
      fasta_file=`basename ${path%.*}`
      echo "Extra: $fasta_file"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      mmseqs easy-search "${path}" 13_pan_aug_fasta.faa 13_extra_out_dir/${fasta_file}.x.all_cons.m8 \
                   $MMTEMP --search-type ${SEQTYPE} --cov-mode 5 -c ${clust_cov} 1>/dev/null & # background

      if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    wait # wait for jobs to finish
  
    echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
    echo "  Use identity threshold extr_iden: $extra_iden."
    export FAM_PRE=${consen_prefix}
    top_line.awk 13_extra_out_dir/*.x.all_cons.m8 |
      awk -v IDEN=${extra_iden} '$3>=IDEN {print $2 "\t" $1}' | 
      perl -pe '$prefix=$ENV{FAM_PRE}; s/^($prefix\d+)__\S+/$1/' |  # trims gene ID from preceding family ID
      sort -k1,1 -k2,2 | hash_to_rows_by_1st_col.awk > 14_syn_pan_extra.clust.tsv
  
    echo "  Retrieve sequences for the extra genes"
    if [ -d 16_pan_leftovers_extra ]; then rm -rf 16_pan_leftovers_extra; fi
    mkdir -p 16_pan_leftovers_extra
    get_fasta_from_family_file.pl "${protein_files_extra[@]}" \
       -fam 14_syn_pan_extra.clust.tsv -out 16_pan_leftovers_extra/
  
    echo "  Make augmented cluster sets"
    augment_cluster_sets.awk leftovers_dir=16_pan_leftovers_extra 13_pan_aug_fasta/* |
      cat > 18_syn_pan_aug_extra.clust.tsv

    echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
    perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 18_syn_pan_aug_extra.clust.tsv \
      > 18_syn_pan_aug_extra.hsh.tsv

    echo "  Merge fasta sets (in parallel)"
    if [ -d batches ]; then rm -rf batches; fi
    mkdir -p batches
    ls 13_pan_aug_fasta | split - batches/
    if [ -d 19_pan_aug_leftover_merged_prot ]; then rm -rf 19_pan_aug_leftover_merged_prot; fi
    mkdir -p 19_pan_aug_leftover_merged_prot
    for batch in batches/*; do
      cat $batch | while read -r panID; do
        #echo "$batch: $panID"
        if [[ -f "16_pan_leftovers_extra/$panID" ]]; then
          cat 13_pan_aug_fasta/$panID 16_pan_leftovers_extra/$panID > 19_pan_aug_leftover_merged_prot/$panID
        else
          cp 13_pan_aug_fasta/$panID 19_pan_aug_leftover_merged_prot/
        fi
      done & # Do the merging in parallel because of the number of flies
      if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
    done
    wait 
  
    echo "  Get all protein sequences from files in 19_pan_aug_leftover_merged_prot"
    cat /dev/null > 19_pan_aug_leftover_merged_prot.faa
    for path in 19_pan_aug_leftover_merged_prot/*; do
      pan_file=`basename $path`
      cat $path | awk -v panID=$pan_file ' $1~/^>/ {print ">" panID "__" substr($0,2) }
                        $1!~/^>/ {print $1} ' >> 19_pan_aug_leftover_merged_prot.faa 
    done

  else  
    echo "== No annotations were designated as \"extra\", so just promote the syn_pan_aug files as syn_pan_aug_extra. ==" 
    cp 07_pan_fasta_prot.faa 19_pan_aug_leftover_merged_prot.faa
    cp 12_syn_pan_aug.clust.tsv 18_syn_pan_aug_extra.clust.tsv
    cp 12_syn_pan_aug.hsh.tsv 18_syn_pan_aug_extra.hsh.tsv
  fi
}

##########
run_align() {
  echo; echo "== Retrieve sequences for each family, preparatory to aligning them =="
  cd "${WORK_DIR}"
  if [[ -d 19_pan_aug_leftover_merged_prot ]] && [[ -f 19_pan_aug_leftover_merged_prot/${consen_prefix}00001 ]]; then
    : # do nothing; the directory and file(s) exist
  else 
    mkdir -p 19_pan_aug_leftover_merged_prot
    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    get_fasta_from_family_file.pl "${protein_files[@]}" "${protein_files_extra[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_prot
  fi

  echo; echo "== Move small (<4) and very low-entropy families (sequences are all identical) to the side =="
  mkdir -p 19_pan_aug_small_or_identical
  min_seq_count=4
  for filepath in 19_pan_aug_leftover_merged_prot/*; do
    file=`basename $filepath`
    count=$(awk '$1!~/>/ {print FILENAME "\t" $1}' $filepath | sort -u | wc -l);  
    if [[ $count -lt $min_seq_count ]]; then 
      echo "Set aside small or low-entropy family $file";
      mv $filepath 19_pan_aug_small_or_identical/
    fi; 
  done

  echo; echo "== Align the gene families =="
  mkdir -p 20_aligns
  for filepath in 19_pan_aug_leftover_merged_prot/*; do 
    file=`basename $filepath`;
    echo "  Computing alignment, using program famsa, for file $file"
    famsa -t 2 19_pan_aug_leftover_merged_prot/$file 20_aligns/$file 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_model_and_trim() {
  echo; echo "== Build HMMs =="
  cd "${WORK_DIR}"
  mkdir -p 21_hmm
  for filepath in 20_aligns/*; do 
    file=`basename $filepath`;
    hmmbuild -n $file 21_hmm/$file $filepath 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
 
  echo; echo "== Realign to HMMs =="
  mkdir -p 22_hmmalign
  for filepath in 21_hmm/*; do 
    file=`basename $filepath`;
    hmmalign --trim --outformat A2M --amino -o 22_hmmalign/$file 21_hmm/$file 19_pan_aug_leftover_merged_prot/$file &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Trim HMM alignments to match-states =="
  mkdir -p 23_hmmalign_trim1
  for filepath in 22_hmmalign/*; do 
    file=`basename $filepath`;
    cat $filepath | 
      perl -ne 'if ($_ =~ />/) {print $_} else {$line = $_; $line =~ s/[a-z]//g; print $line}' |
      sed '/^$/d' > 23_hmmalign_trim1/$file &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Filter alignments prior to tree calculation =="
  mkdir -p 23_hmmalign_trim2 23_hmmalign_trim2_log
  min_depth=3
  min_pct_depth=20
  min_pct_aligned=20
  for filepath in 23_hmmalign_trim1/*; do 
    file=`basename $filepath`
    filter_align.pl -in $filepath -out 23_hmmalign_trim2/$file -log 23_hmmalign_trim2_log/$file \
                    -depth $min_depth -pct_depth $min_pct_depth -min_pct_aligned $min_pct_aligned &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_calc_trees() {
  echo; echo "== Calculate trees =="
  cd "${WORK_DIR}"
  
  mkdir -p 24_trees
  
  # By default, FastTreeMP uses all machine cores, e.g. 64 on lathyrus. 
  # It is more efficient to run more jobs on one core each by setting an environment variable.
  OMP_NUM_THREADS=1
  export OMP_NUM_THREADS
  for filepath in 23_hmmalign_trim2/*; do
    file=`basename $filepath`
    echo "  Calculating tree for $file"
    fasttree -quiet $filepath > 24_trees/$file &
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
 
  #param_string="${ST}.id${clust_iden}.cov${clust_cov}.cns${consen_iden}.ext${extra_iden}.I${mcl_inflation}"
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

  for dir in 19_pan_aug_leftover_merged_prot 21_hmm 22_hmmalign 23_hmmalign_trim2 24_trees; do
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
export MMSEQS_NUM_THREADS=${NPROC} # mmseqs otherwise uses all cores by default

# mmseqs uses significant number of threads on its own. Set a maximum, which may be below NPROC.
MMSEQSTHREADS=$(( 4 < ${NPROC} ? 4 : ${NPROC} ))

pandagma_conf_params='clust_iden clust_cov consen_iden extra_iden mcl_inflation strict_synt dagchainer_args 
  out_dir_base pctl_low pctl_med pctl_hi consen_prefix annot_str_regex work_dir'

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
for program in mmseqs dagchainer mcl cons famsa run_DAG_chainer.pl; do
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
commandlist="ingest mmseqs filter dagchainer mcl consense cluster_rest add_extra align model_and_trim calc_trees summarize"

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


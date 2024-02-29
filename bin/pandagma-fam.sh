#!/usr/bin/env bash
#
# Configuration and run script which, with other scripts in this package, generates 
# gene-family orthogroups using the programs mmseqs, dagchainer, and mcl. 
# Authors: Steven Cannon, Hyunoh Lee, Joel Berendzen, Nathan Weeks, 2020-2023
#
scriptname='pandagma fam'

# shellcheck source=/dev/null
. pandagma-common.sh 

define HELP_DOC <<'EOS'
Compute orthogroups using a combination of synteny and homology,
using the programs mmseqs, dagchainer, and mcl, and additional pre- and post-refinement steps.

Usage:
  $scriptname -c CONFIG_FILE [options]

  Required:
           -c (path to the config file)

  Options: -s (subcommand to run. If \"all\" or omitted, all steps will be run; otherwise, run specified step)
           -w (working directory, for temporary and intermediate files [default: './work_pandagma'].)
           -o OUTPUT_DIR (name for the output directory [default: './out_pandagma'].
                Applicable only to "all", "summarize", and "xfr_aligns_trees" steps.)
           -n (number of processors to use. Defaults to number of processors available to the parent process)
           -r (retain. Don't do subcommand \"clean\" after running \"all\".)
           -v (version)
           -h (help)
           -m (more information)

Environment requirements: The following packages need to be available in your PATH:
    mmseqs dagchainer mcl cons famsa hmmalign hmmbuild fasttree

Also, please add the pandagma utility programs in the bin directory adjacent to pandagma-fam.sh, e.g.
    PATH=$PWD/bin:\$PATH

Primary protein sequences and annotation (GFF3 or BED) files must be listed in the
config file, in the arrays annotation_files and protein_files. See example files.

Subcommands (in order they are usually run):
  Run these first (if using the ks_peaks.tsv file; otherwise, run all main steps and
  ks filtering will be done using parameters ks_block_wgd_cutoff and max_pair_ks)
                all - All of the steps below, except for ks_filter, clean
                        (Or equivalently: omit the -s flag; \"all\" is default).
             ingest - Prepare the assembly and annotation files for analysis.
             mmseqs - Run mmseqs to do initial clustering of genes from pairs of assemblies.
             filter - Filter the synteny results for chromosome pairings, returning gene pairs.
         dagchainer - Run DAGchainer to filter for syntenic blocks.
            ks_calc - Calculation of Ks values on gene pairs from DAGchainer output.

  Evaluate the stats/ks_histplots.tsv and stats/ks_peaks_auto.tsv files and
  put ks_peaks.tsv into the \${WORK_DIR}/stats directory, then run the following commands:
          ks_filter - Filtering based on provided ks_peaks.tsv file (assumes prior ks_calc step)
                mcl - Derive clusters, with Markov clustering.
           consense - Calculate a consensus sequences from each pan-gene set, 
                      adding sequences missed in the first clustering round.
       cluster_rest - Retrieve unclustered sequences and cluster those that can be.
          add_extra - Add other gene model sets to the primary clusters. Useful for adding
                      annotation sets that may be of lower or uncertain quality.
         tabularize - Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv
          summarize - Move results into output directory, and report summary statistics.

  If generating alignments, models, and trees, run the following steps:
              align - Align families.
     model_and_trim - Build HMMs and trim the alignments, preparatory to calculating trees.
         calc_trees - Calculate gene trees.

  Run the following subcommands separately if you wish:
              clean - Clean (delete) files in the working directory that are not needed 
                        for later addition of data using add_extra and subsequent run commands.
                        By default, \"clean\" is run as part of \"all\" unless the -r flag is set.
EOS

define MORE_INFO <<'EOS'
Two options are provided for filtering gene homologies based on additional provided information.

1. One option for additional filtering is to provide a file of gene matches for each species. 
In this workflow, this set of values is termed \"expected_quotas\", and can be provided
in the family conf file.
This indicates the number of expected paralogs for a species in a gene family at the desired evolutionary depth,
considering a known history of whole-genome duplications that affect the included species. For example,
for the legume family, originating ~70 mya:

expected_quotas=(
  Arachis     4
  Cercis      1
  Glycine     4
  Phaseolus   2
)

These quotas are used in a regular expression to identify the initial portion of the prefixed seqIDs,
for example, \"Cicer\" and \"cerca\" matching the genus and \"gensp\" matching prefixes for these seqIDs:
  Cicer.pan1.chr06
  cerca.ISC453364.gnm3.Chr03

2. The second option for additional filtering is to provide a file of Ks peak values for synteny-block-median Ks values.
This requires calculating gene-pair and block-median Ks values from DAGChainer output, using the first five steps,
through the \"ks_calc\" . Several ks-related parameters also need to be provided in the configuration file. These
are described briefly below. 

Variables in pandagma config file (Set the config with the CONF environment variable)
         clust_iden - Minimum identity threshold for mmseqs clustering [0.40]
          clust_cov - Minimum coverage for mmseqs clustering [0.40]
         extra_iden - Minimum identity threshold for mmseqs addition of \"extra\" annotations [0.30]
      mcl_inflation - Inflation parameter, for Markov clustering [1.6]
        strict_synt - For clustering of the \"main\" annotations, use only syntenic pairs [1]
                        The alternative (0) is to use all homologous pairs that satisfy expected_quotas
      ks_low_cutoff - For inferring Ks peak per species pair. Don't consider Ks block-median values less than this. [0.5]
       ks_hi_cutoff - For inferring Ks peak per species pair. Don't consider Ks block-median values greater than this. [2.0]
         ks_binsize - For calculating and displaying histograms. [0.05]
ks_block_wgd_cutoff - Fallback, if a ks_peaks.tsv file is not provided. [1.75]
        max_pair_ks - Fallback value for excluding gene pairs, if a ks_peaks.tsv file is not provided. [4.0]

      consen_prefix - Prefix to use in orthogroup names
    annot_str_regex - Regular expression for capturing annotation name from gene ID, e.g. 
                        \"([^.]+\.[^.]+)\..+\"
                          for two dot-separated fields, e.g. vigan.Shumari
                        or \"(\D+\d+\D+)\d+.+\" for Zea assembly+annot string, e.g. Zm00032ab
    preferred_annot - String to match and select an annotation set, from a gene ID.
                        This is used for picking representative IDs+sequence from an orthogroup, when
                        this annotation is among those with the median length for the orthogroup.
                        Otherwise, one is selected at random from those with median length.
    expected_quotas - (Optional) array of seqid prefixes & expected number of
                        paralogs for the species identified by the prefix; e.g.:
                        expected_quotas=(glyma 4 medtr 2)

File sets (arrays):
   annotation_files
      protein_files
          cds_files
annotation_files_extra
   protein_files_extra

Optional files with cutoff values
          ks_peaks.tsv
EOS

########################################
# Helper functions begin here

canonicalize_paths() {
  cd "${DATA_DIR}" || exit
  echo "Entering canonicalize_paths. Annotation files: " "${annotation_files[@]}"

  mapfile -t cds_files < <(realpath --canonicalize-existing "${cds_files[@]}")
  mapfile -t annotation_files < <(realpath --canonicalize-existing "${annotation_files[@]}")
  mapfile -t protein_files < <(realpath --canonicalize-existing "${protein_files[@]}")

  if [[ -v protein_files_extra ]]
  then
    mapfile -t protein_files_extra < <(realpath --canonicalize-existing "${protein_files_extra[@]}")
  fi

  if [[ -v cds_files_extra ]]
  then
    mapfile -t cds_files_extra < <(realpath --canonicalize-existing "${cds_files_extra[@]}")
    mapfile -t annotation_files_extra < <(realpath --canonicalize-existing "${annotation_files_extra[@]}")
  fi

  cd "${OLDPWD}" || exit
  readonly submit_dir=${PWD}

  fasta_file=$(basename "${protein_files[0]}" .gz)
  fa="${fasta_file##*.}"
}

########################################
# run functions 

##########
run_ingest() {
# Add positional information from GFF3 or 4- or 6-column BED to FASTA IDs
# BED start coordinate converted to 1-based
  cd "${WORK_DIR}" || exit
  echo; echo "Run ingest: from fasta and gff or bed data, create fasta with IDs containing positional info."
  
  mkdir -p 02_fasta_nuc 02_fasta_prot 01_posn_hsh stats

    # Prepare the tmp.gene_count_start to be joined, in run_summarize.
    # This is captured from the gene IDs using the annot_str_regex set in the config file.
    cat /dev/null > stats/tmp.gene_count_start
    cat /dev/null > stats/tmp.fasta_seqstats
    start_time=$(date)
    printf "Run started at: %s\n" "$start_time" > stats/tmp.timing

  export ANN_REX=${annot_str_regex}

  echo "  Get position information from the main annotation sets (protein)."
  cat /dev/null > 02_all_main_prot.faa # Collect all starting protein sequences, for later comparisons
  for (( file_num = 0; file_num < ${#protein_files[@]} ; file_num++ )); do
    file_base=$(basename "${protein_files[file_num]%.*}")
    cat_or_zcat "${protein_files[file_num]}" >> 02_all_main_prot.faa # Collect original seqs for later comparisons
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | 
      gff_or_bed_to_hash5.awk > 01_posn_hsh/"$file_base".hsh
    hash_into_fasta_id.pl -nodef -fasta "${protein_files[file_num]}" \
                          -hash 01_posn_hsh/"$file_base".hsh \
                          -out 02_fasta_prot/"$file_base"
    # calc basic sequence stats
    annot_name=$(basename 02_fasta_prot/"$file_base" | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
    printf "  Main:  " >> stats/tmp.fasta_seqstats
    cat_or_zcat "${protein_files[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
  done

  echo "  Get position information from the main annotation sets (cds)."
  cat /dev/null > 02_all_main_cds.fna # Collect all starting cds sequences, for later comparisons
  for (( file_num = 0; file_num < ${#cds_files[@]} ; file_num++ )); do
    file_base=$(basename "${cds_files[file_num]%.*}")
    cat_or_zcat "${cds_files[file_num]}" >> 02_all_main_cds.fna # Collect original seqs for later comparisons
    echo "  Adding positional information to fasta file $file_base"
    cat_or_zcat "${annotation_files[file_num]}" | 
      gff_or_bed_to_hash5.awk > 01_posn_hsh/"$file_base".hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files[file_num]}" \
                          -hash 01_posn_hsh/"$file_base".hsh \
                          -out 02_fasta_nuc/"$file_base"
  done

  echo "  Get position information from the extra annotation sets (protein), if any."
  if [[ -v protein_files_extra ]]
  then
    cat /dev/null > 02_all_extra_protein.faa # Collect all starting sequences, for later comparisons
    for (( file_num = 0; file_num < ${#protein_files_extra[@]} ; file_num++ )); do
      file_base=$(basename "${protein_files_extra[file_num]%.*}")
      cat_or_zcat "${protein_files_extra[file_num]}" >> 02_all_extra_protein.faa  # Collect original seqs for later comparisons
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra[file_num]}" | 
        gff_or_bed_to_hash5.awk > 01_posn_hsh/"$file_base".hsh
      hash_into_fasta_id.pl -nodef -fasta "${protein_files_extra[file_num]}" \
                            -hash 01_posn_hsh/"$file_base".hsh \
                            -out 02_fasta_prot/"$file_base"
      # calc basic sequence stats
      annot_name=$(basename 02_fasta_prot/"$file_base" | perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' )
      echo "  CHECK: report annot_name to stats file? [$annot_name]"
      printf "  Extra: " >> stats/tmp.fasta_seqstats
      cat_or_zcat "${protein_files_extra[file_num]}" | calc_seq_stats >> stats/tmp.fasta_seqstats
    done
  fi

  echo "  Get position information from the extra annotation sets (cds), if any."
  if [[ -v cds_files_extra ]]
  then
    cat /dev/null > 02_all_extra_cds.fna # Collect all starting sequences, for later comparisons
    for (( file_num = 0; file_num < ${#cds_files_extra[@]} ; file_num++ )); do
      file_base=$(basename "${cds_files_extra[file_num]%.*}")
      cat_or_zcat "${cds_files_extra[file_num]}" >> 02_all_extra_cds.fna  # Collect original seqs for later comparisons
      echo "  Adding positional information to extra fasta file $file_base"
      cat_or_zcat "${annotation_files_extra[file_num]}" | 
        gff_or_bed_to_hash5.awk > 01_posn_hsh/"$file_base".hsh
      hash_into_fasta_id.pl -nodef -fasta "${cds_files_extra[file_num]}" \
                            -hash 01_posn_hsh/"$file_base".hsh \
                            -out 02_fasta_nuc/"$file_base"
    done
  fi

  echo "  Count starting sequences, for later comparisons"
  for file in 02_fasta_prot/*."$fa"; do
    awk '$0~/UNDEFINED/ {ct++} 
      END{if (ct>0){print "Warning: " FILENAME " has " ct " genes without position (HASH UNDEFINED)" } }' "$file"
    grep '>' "$file" | perl -pe 's/__/\t/g' | cut -f2 | # extracts the GeneName from combined genspchr__GeneName__start__end__orient
      perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex.+/$1/' |
      grep -v UNDEFINED | sort | uniq -c | awk '{print $2 "\t" $1}' >> stats/tmp.gene_count_start
  done

  sort -o stats/tmp.gene_count_start stats/tmp.gene_count_start
}

##########
run_mmseqs() {
  # Do mmseqs clustering on all pairings of the main annotation sets (not the extra ones though)
  cd "${WORK_DIR}" || exit
  echo; echo "Run mmseqs -- at ${clust_iden} percent identity and minimum of ${clust_cov}% coverage."
  #

  if [ -d 03_mmseqs ]; then rm -rf 03_mmseqs ; fi
  mkdir -p 03_mmseqs 03_mmseqs_tmp
  SEQTYPE=1;  # 1=pep; 3=nuc
  for (( file1_num = 0; file1_num < ${#protein_files[@]} ; file1_num++ )); do
    qry_base=$(basename "${protein_files[file1_num]%.*}" ."$fa")
    for (( file2_num = file1_num; file2_num < ${#protein_files[@]} ; file2_num++ )); do  # file2_num = $file1_num includes self-comparisons
      sbj_base=$(basename "${protein_files[file2_num]%.*}" ."$fa")
      echo "  Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)

      if [[ ${qry_base} == "${sbj_base}" ]]; then # self-comparison, so use flag --add-self-matches
        mmseqs easy-search \
          02_fasta_prot/"$qry_base"."$fa" 02_fasta_prot/"$sbj_base"."$fa" 03_mmseqs/"${qry_base}".x."${sbj_base}".m8 "$MMTEMP" \
          --add-self-matches --search-type ${SEQTYPE} --cov-mode 0 -c "${clust_cov}" --min-seq-id "${clust_iden}" \
          --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln" 1>/dev/null 
      else # not a self-comparison, do omit flag --add-self-matches
        mmseqs easy-search \
          02_fasta_prot/"$qry_base"."$fa" 02_fasta_prot/"$sbj_base"."$fa" 03_mmseqs/"${qry_base}".x."${sbj_base}".m8 "$MMTEMP" \
          --search-type ${SEQTYPE} --cov-mode 0 -c "${clust_cov}" --min-seq-id "${clust_iden}" \
          --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln" 1>/dev/null 
      fi
    done
    echo
  done
  wait # wait for last jobs to finish
}

##########
run_filter() {
  echo; echo "From homology (mmseqs -m8) output, split out the following fields: molecule, gene, start, stop."
  cd "${WORK_DIR}" || exit
  if [ -d 04_dag ]; then rm -rf 04_dag ; fi
  mkdir -p 04_dag
  if [[ -v expected_quotas ]]; then  # filter based on list of match quotas if provided
    echo "Filtering on quotas from expected_quotas"
    for mmseqs_path in 03_mmseqs/*.m8; do
      outfilebase=$(basename "$mmseqs_path" .m8)
      echo "  Filtering $outfilebase based on expected quotas"
      filter_mmseqs_by_quotas.pl -quotas <(printf '%s %s\n' "${expected_quotas[@]}") < "${mmseqs_path}" |
        perl -pe 's/\t[\+-]//g' |  # strip orientation, which isn't used by DAGChainer 
        cat | # Next: for self-comparisons, suppress same chromosome && different gene ID (local dups)
        perl -lane 'print $_ unless($F[0] eq $F[4] && $F[1] ne $F[5])' |
        awk 'NF==8' |  # matches for genes with coordinates. The case of <8 can happen for seqs with UNDEF position.
        cat > 04_dag/"${outfilebase}"_matches.tsv &

      # allow to execute up to $NPROC in parallel
      if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
    done
    wait # wait for last jobs to finish
  else   # don't filter, since quotas file isn't provided; just remove orientation, which isn't used by DAGChainer
    echo "No expected_quotas provided, so proceeding without quota (expected gene-count) filtering."
    for mmseqs_path in 03_mmseqs/*.m8; do
      outfilebase=$(basename "$mmseqs_path" .m8)
      cut -f1,2 "${mmseqs_path}" | 
        perl -pe 's/__/\t/g; s/\t[\+-]//g' | 
        cat | # Next: for self-comparisons, suppress same chromosome && different gene ID (local dups)
        perl -lane 'print $_ unless($F[0] eq $F[4] && $F[1] ne $F[5])' |
        awk 'NF==8' |  # matches for genes with coordinates. The case of <8 can happen for seqs with UNDEF position.
        cat > 04_dag/"${outfilebase}"_matches.tsv 
    done
  fi
}

##########
run_dagchainer() {
  # Identify syntenic blocks, using DAGchainer
  cd "${WORK_DIR}" || exit
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
    ave_gene_gap=$(cat 02_fasta_prot/"$qryfile"."$fa" 02_fasta_prot/"$sbjfile"."$fa" | 
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
      cd "${tmpdir}" || exit
      run_DAG_chainer.pl "$dagchainer_args"  -g "$ave_gene_gap" -D "$max_gene_gap" -i "${OLDPWD}/${match_path}" 1>/dev/null
      rmdir "${tmpdir}"
    ) &
    # allow to execute up to $NPROC in parallel
    if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
  done
  wait # wait for last jobs to finish
}

##########
run_ks_calc() {
  echo; echo "Calculate Ks values from processed DAGChainer output, generating gene-pair and "
  echo       "block-median Ks values for subsequent filtering. For this step, alignments from mmseqs "
  echo       "are converted to berkeleydb key-value hashes, where the key is composed of "
  echo       "\"query subject\" or \"subject query\", in lexical order -- since the output of "
  echo       "DAGChainer has query and subject in lexical order."
  echo "NPROC: ${NPROC}"

  cd "${WORK_DIR}" || exit

  if [ ! -d "$WORK_DIR"/05_kaksout ]; then
    echo "creating output directory $WORK_DIR/05_kaksout"
    mkdir -p "$WORK_DIR"/05_kaksout
  fi
  
  echo "  Create berkeleydb index file (directory.index) for fasta file(s) 02_*_cds.fna"
  perl -MBio::DB::Fasta -e 'Bio::DB::Fasta->new(".", -glob=>"02_*_cds.fna")'
    
  for DAGFILE in 04_dag/*aligncoords; do
    base=$(basename "$DAGFILE" _matches.tsv.aligncoords)
    echo "  Calculate Ks values for $base"

    echo "cat $DAGFILE | calc_ks_from_dag.pl -fasta_db . \\ ";
    echo "   -align_method precalc -match_table 03_mmseqs/$base.db \\ ";
    echo "  -report_out 05_kaksout/$base.rptout";
    < "$DAGFILE" calc_ks_from_dag.pl -fasta_db . \
      -match_table 03_mmseqs/"$base".m8  -align_method precalc \
      -report_out 05_kaksout/"$base".rptout & #2> /dev/null & # discard verbose warnings from LocatableSeq.pm
    echo
    if [[ $(jobs -r -p | wc -l) -ge ${NPROC} ]]; then wait -n; fi
  done
  wait

  echo "  Delete Berkeleydb files (they derive from the .m8 files in 03_mmseqs and 02_*_cds.fna)"
  rm -f 03_mmseqs/*.index directory.index # Bio::DB::Fasta index

  echo "Determine provisional Ks peaks (Ks values and amplitudes) and generate Ks plots."
  cat /dev/null > stats/ks_peaks_auto.tsv
  cat /dev/null > stats/ks_histograms.tsv
  cat /dev/null > stats/ks_histplots.tsv
  for ks_path in 05_kaksout/*.rptout; do
    filebase=$(basename "$ks_path" .rptout)
    qry_ann=$(echo "$filebase" | perl -pe 'BEGIN{$REX=$ENV{"ANN_REX"}}; s/^$REX/$1/')
    sbj_ann=$(echo "$filebase" | perl -pe 'BEGIN{$REX=qr($ENV{"ANN_REX"})}; s/.+\.x\.$REX/$1/')

    < "$ks_path" awk 'NF==7 && $7<=3 {print $7}' | histogram.pl -z -n -s "$ks_binsize" |
      awk -v KSLC="$ks_low_cutoff" -v KSHC="$ks_hi_cutoff" -v QA="$qry_ann" -v SA="$sbj_ann" '
        $1>=KSLC && $1<=KSHC && $2>=maxampl { maxampl=$2; maxbin=$1 } 
        END{ printf("%s\t%s\t%.2f\t%d\n", QA, SA, maxbin, maxampl)}' >> stats/ks_peaks_auto.tsv

    ks_bin=$(< stats/ks_peaks_auto.tsv awk -v QA="$qry_ann" -v SA="$sbj_ann" -v OFS="\t" \
                                            '$1 == QA && $2 == SA {print $3}')
    ks_amplitude=$(< stats/ks_peaks_auto.tsv awk -v QA="$qry_ann" -v SA="$sbj_ann" -v OFS="\t" \
                                            '$1 == QA && $2 == SA {print $4}')

    export ks_amplitude
    ks_amplitude_pct=$(perl -e '$KSA=$ENV{ks_amplitude}; {printf("%.2f", $KSA/100)}')
      
    printf "  Amplitude, ks_peak_bin, qry, sbj:\t%s\t%s\t%s\t%s\n" "$ks_amplitude" "$ks_amplitude_pct" "$qry_ann" "$sbj_ann"

    {
      echo "# $filebase" 
      awk 'NF==7 && $7<=3 {print $7}' "$ks_path" | histogram.pl -z -n -s "$ks_binsize" 
      echo "" 
    } >> stats/ks_histograms.tsv

    {
      echo "# $filebase" 
      echo "# Normalized relative to Ks peak inferred at bin $ks_bin, with amplitude $ks_amplitude_pct" 
      awk 'NF==7 && $7<=3 {print $7}' "$ks_path" | histogram.pl -z -n -s "$ks_binsize" | 
        histplot.pl -d "$ks_amplitude_pct" 
      echo "" 
    } >> stats/ks_histplots.tsv
  done
}

##########
run_ks_filter() {
  echo; 
  echo "From optional provided ${WORK_DIR}/stats/ks_peaks.tsv file and from processed DAGChainer output,i "
  echo "(with gene-pair and block-median Ks values added by calc_ks_from_dag.pl), filter on "
  echo "block-median Ks values, and combine the homology data into a file with gene pairs to be clustered."
  echo;

  cd "${WORK_DIR}" || exit
  if [ -f stats/ks_peaks.tsv ]; then  # filter based on list of Ks values
    echo "Filtering on quotas from expected_quotas and ks_pair_cutoff values provided in ks_peaks.tsv"
    if [ -d 05_kaksout_ks_filtered ]; then rm -rf 05_kaksout_ks_filtered ; fi
    mkdir -p 05_kaksout_ks_filtered
    ks_peaks=stats/ks_peaks.tsv
    for ks_path in 05_kaksout/*.rptout; do
      outfilebase=$(basename "$ks_path" .rptout)
      echo "  Filtering $ks_path based on expected block-median Ks values"
      < "${ks_path}" filter_mmseqs_by_ks.pl \
          -ks_peak "${ks_peaks}" -annot_regex "$annot_str_regex" -max_pair_ks "$max_pair_ks" |
        awk 'NF==7' > 05_kaksout_ks_filtered/"${outfilebase}".rptout 
    done
    cat 05_kaksout_ks_filtered/*.rptout | awk 'NF==7 {print $1 "\t" $2}' | sort -u > 05_filtered_pairs.tsv
  else # ks_peaks.tsv file isn't provided. Warn but proceed with a single Ks value: ks_hi_cutoff
    echo "INFO: No ks_peaks.tsv file was provided. It is recommended to review the provisional stats/ks_peaks_auto.tsv"
    echo "and stats/ks_histplots.tsv files and provide a file stats/ks_peaks.tsv -- either simply copying"
    echo "ks_peaks_auto.tsv to ks_peaks.tsv if the reported Ks peak values (column 3) are acceptable,"
    echo "or revising those values based on interpretation of stats/ks_histplots.tsv."
    echo "If stats/ks_histplots.tsv is not provided, then Ks filtering will be done using values provided"
    echo "in the config file for ks_block_wgd_cutoff and max_pair_ks."; echo

    if [ "$strict_synt" -eq 1 ]; then
      # https://stackoverflow.com/questions/3601515/how-to-check-if-a-variable-is-set-in-bash
      if [ -z ${ks_block_wgd_cutoff+x} ] || [ -z ${max_pair_ks+x} ] || 
          [ ! -d "${WORK_DIR}/05_kaksout" ] && [ ! "$(ls -A 05_kaksout/*.rptout)" ]; then
        echo "Variables ks_block_wgd_cutoff and/or max_pair_ks are unset or 05_kaksout doesn't exist; no Ks filtering will be done.";
        echo "Combine the DAGChainer synteny pairs into a file to be clustered."
        cat 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' | awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
      else 
        echo "The ks_block_wgd_cutoff is set to '$ks_block_wgd_cutoff' and max_pair_ks is set to '$max_pair_ks'";
        echo "Combine DAGChainer synteny pairs, filtering by pairwise and block Ks thresholds, into a file to be clustered."
        # Combine the synteny pairs into a file to be clustered
        cat 05_kaksout/*.rptout | 
          awk -v PAIR_CUTOFF="$max_pair_ks" -v BLOCK_CUTOFF="$ks_block_wgd_cutoff" '
              NF>0 && $1!~/^#/ && ($5<=PAIR_CUTOFF || $7<=BLOCK_CUTOFF) {print $1 "\t" $2}' |
            sort -u > 05_filtered_pairs.tsv
      fi
    else 
      echo "## Combine the homology pairs (filtered by quota if provided) into a file to be clustered."
      cat 04_dag/*_matches.tsv 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' | 
        awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
    fi
  fi
}

##########
run_mcl() {
  # Calculate clusters using Markov clustering
  cd "${WORK_DIR}" || exit

  printf "\nDo Markov clustering with inflation parameter %f and %d threads\n" "$mcl_inflation" "$NPROC"
  echo "MCL COMMAND: mcl 05_filtered_pairs.tsv -I $mcl_inflation -te ${NPROC} --abc -o tmp.syn_pan.clust.tsv"
  mcl 05_filtered_pairs.tsv -I "$mcl_inflation" -te "${NPROC}" --abc -o tmp.syn_pan.clust.tsv \
    1>/dev/null
 
  echo "  Add cluster IDs"
  awk -v PRE="$consen_prefix" '
    {padnum=sprintf("%05d", NR); print PRE padnum "\t" $0}' tmp.syn_pan.clust.tsv > 06_syn_pan.clust.tsv
  rm tmp.syn_pan.clust.tsv

  echo "  Move singleton and doubleton clusters into a list of leftovers"
  mkdir -p lists
  awk 'NF==2 {print $2} NF==3 {print $2; print $3}' 06_syn_pan.clust.tsv > lists/lis.06_syn_pan_1s_2s
  awk 'NF>3 {print $0}' 06_syn_pan.clust.tsv > 06_syn_pan_ge3.clust.tsv

  echo "  Reshape from mcl output format (clustered IDs on one line) to a hash format (clust_ID gene)"
  perl -lane 'for $i (1..scalar(@F)-1){print $F[0], "\t", $F[$i]}' 06_syn_pan_ge3.clust.tsv > 06_syn_pan_ge3.hsh.tsv
}

##########
run_consense() {
  echo; 
  echo "Add previously unclustered sequences into an \"augmented\" pan-gene set, by homology."
  cd "${WORK_DIR}" || exit
  if [ -d 07_pan_fasta ]; then rm -rf 07_pan_fasta; fi
  mkdir -p 07_pan_fasta lists


  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  echo "    Fasta file:" "${protein_files[@]}"
  get_fasta_from_family_file.pl "${protein_files[@]}" -fam 06_syn_pan_ge3.clust.tsv -out 07_pan_fasta

  echo "  Merge fasta files in 07_pan_fasta, prefixing them with panID__"
  merge_files_to_pan_fasta.awk 07_pan_fasta/* > 07_pan_fasta_prot.faa

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
                     --search-type "${SEQTYPE}" --cov-mode 5 -c "${clust_cov}" 1>/dev/null 

  echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
  echo "  Use the \"main set\" $clust_iden threshold."
  export FAM_PRE=${consen_prefix}
  top_line.awk 10_unclust.x.07_pan_fasta.m8 | 
    awk -v IDEN="${clust_iden}" '$3>=IDEN {print $2 "\t" $1}' | 
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
  cd "${WORK_DIR}" || exit

  echo "  Retrieve genes present in the original protein files but absent from 12_syn_pan_aug.hsh"
  cut -f2 12_syn_pan_aug_pre.hsh.tsv | LC_ALL=C sort > lists/lis.12_syn_pan_aug_complement
  get_fasta_subset.pl -in 02_all_main_prot.faa -out 12_syn_pan_aug_complement.faa \
    -lis lists/lis.12_syn_pan_aug_complement -xclude -clobber

  MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
  mkdir -p 03_mmseqs_rest
  complmt_self_compare="12_syn_pan_aug_complement.x.12_syn_pan_aug_complement"
  < 12_syn_pan_aug_complement.faa \
    mmseqs easy-cluster stdin 03_mmseqs_rest/$complmt_self_compare "$MMTEMP" \
    --min-seq-id "$clust_iden" -c "$clust_cov" --cov-mode 0 --cluster-reassign 1>/dev/null

  echo "  Cluster the remaining sequences that have matches"
  mcl 03_mmseqs_rest/${complmt_self_compare}_cluster.tsv -I "$mcl_inflation" -te "${NPROC}" --abc -o tmp.syn_pan_aug_complement.clust.tsv

  echo "  Find number of clusters in initial (06) results"
  clust_count_06=$(wc -l 06_syn_pan_ge3.clust.tsv | awk '{print $1}')
 
  echo "  Add cluster IDs"
  < tmp.syn_pan_aug_complement.clust.tsv awk -v START="$clust_count_06" -v PRE="$consen_prefix" '
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
  cd "${WORK_DIR}" || exit

  if [ -d 13_extra_out_dir ]; then rm -rf 13_extra_out_dir; fi
  if [ -d 13_pan_aug_fasta ]; then rm -rf 13_pan_aug_fasta; fi
  mkdir -p 13_extra_out_dir 13_pan_aug_fasta

  echo "  For each pan-gene set, retrieve sequences into a multifasta file."
  get_fasta_from_family_file.pl "${protein_files[@]}" -fam 12_syn_pan_aug.clust.tsv -out 13_pan_aug_fasta
  
  echo "Merge fasta files in 13_pan_aug_fasta, prefixing IDs with panID__"
  merge_files_to_pan_fasta.awk 13_pan_aug_fasta/* > 13_pan_aug_fasta.faa

  if [[ -v protein_files_extra ]]
  then # handle the "extra" annotation files
    echo "  Search non-clustered genes against pan-gene consensus sequences"
    # Check sequence type (in case this run function is called separately from the usually-prior ones)
    someseq=$(head 07_pan_fasta_prot.faa | grep -v '>' | awk -v ORS="" '{print toupper($1)}')
    SEQTYPE=$(check_seq_type "${someseq}") # 3=nuc; 1=pep
    echo "SEQTYPE is: $SEQTYPE"

    for filepath in "${protein_files_extra[@]}"; do
      fasta_file=$(basename "${filepath%.*}")
      echo "Extra: $fasta_file"
      MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
      mmseqs easy-search "${filepath}" 13_pan_aug_fasta.faa 13_extra_out_dir/"${fasta_file}".x.all_cons.m8 \
                   "$MMTEMP" --search-type "${SEQTYPE}" --cov-mode 5 -c "${clust_cov}" 1>/dev/null & # background

      if [[ $(jobs -r -p | wc -l) -ge ${MMSEQSTHREADS} ]]; then wait -n; fi
    done
    wait # wait for jobs to finish
  
    echo "  Place unclustered genes into their respective pan-gene sets, based on top mmsearch hits."
    echo "  Use identity threshold extr_iden: $extra_iden."
    export FAM_PRE=${consen_prefix}
    top_line.awk 13_extra_out_dir/*.x.all_cons.m8 |
      awk -v IDEN="${extra_iden}" '$3>=IDEN {print $2 "\t" $1}' | 
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

    echo "  For each family set, retrieve sequences into a multifasta file."
    printf "    Fasta file: %s %s\n" "${protein_files[@]}" "${protein_files_extra[@]}"
    if [ -d 19_pan_aug_leftover_merged_prot ]; then rm -rf 19_pan_aug_leftover_merged_prot ; fi
    mkdir -p 19_pan_aug_leftover_merged_prot
    get_fasta_from_family_file.pl "${protein_files[@]}" "${protein_files_extra[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_prot

    echo "  Merge fasta files from 19_pan_aug_leftover_merged_prot, prefixing IDs with panID__"
    merge_files_to_pan_fasta.awk 19_pan_aug_leftover_merged_prot/* > 19_pan_aug_leftover_merged_prot.faa

  else  
    echo "== No annotations were designated as \"extra\", so just promote the syn_pan_aug files as syn_pan_aug_extra. ==" 
    cp 07_pan_fasta_prot.faa 19_pan_aug_leftover_merged_prot.faa
    cp 12_syn_pan_aug.clust.tsv 18_syn_pan_aug_extra.clust.tsv
    cp 12_syn_pan_aug.hsh.tsv 18_syn_pan_aug_extra.hsh.tsv
  fi
}

##########
run_tabularize() {
  echo; echo "== Derive a table-format version of 18_syn_pan_aug_extra.clust.tsv"
  cd "${WORK_DIR}" || exit

  # Get table header
  pangene_tabularize.pl -pan 18_syn_pan_aug_extra.clust.tsv -annot_str_regex "$ANN_REX" > tmp.18_syn_pan_aug_extra.clust.tsv
  head -1 tmp.18_syn_pan_aug_extra.clust.tsv > tmp.table_header

  echo "  Sort, putting header row at top, and don't print pangenes that are all NONE"
    < tmp.18_syn_pan_aug_extra.clust.tsv sort -k1,1 | sed '/^$/d; /^#pangene/d' |
    perl -lane '$ct=0; for $gn (@F){if ($gn=~/NONE/){$ct++}}; if ($ct<(scalar(@F)-1)){print $_}' |
    cat tmp.table_header - > 18_syn_pan_aug_extra.table.tsv

    rm tmp.18_syn_pan_aug_extra.clust.tsv
    rm tmp.table_header
}

##########
run_align() {
  echo; echo "== Retrieve sequences for each family, preparatory to aligning them =="
  cd "${WORK_DIR}" || exit
  if [[ -d 19_pan_aug_leftover_merged_prot ]] && [[ -f 19_pan_aug_leftover_merged_prot/${consen_prefix}00001 ]]; then
    : # do nothing; the directory and file(s) exist
  else 
    mkdir -p 19_pan_aug_leftover_merged_prot
    echo "  For each pan-gene set, retrieve sequences into a multifasta file."
    get_fasta_from_family_file.pl "${protein_files[@]}" "${protein_files_extra[@]}" \
      -fam 18_syn_pan_aug_extra.clust.tsv -out 19_pan_aug_leftover_merged_prot
  fi

  echo; echo "== Align the gene families =="
  mkdir -p 20_aligns
  for filepath in 19_pan_aug_leftover_merged_prot/*; do 
    file=$(basename "$filepath");
    echo "  Computing alignment, using program famsa, for file $file"
    famsa -t 2 19_pan_aug_leftover_merged_prot/"$file" 20_aligns/"$file" 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
}

##########
run_model_and_trim() {
  echo; echo "== Build HMMs =="
  cd "${WORK_DIR}" || exit
  mkdir -p 21_hmm
  for filepath in 20_aligns/*; do 
    file=$(basename "$filepath");
    hmmbuild -n "$file" 21_hmm/"$file" "$filepath" 1>/dev/null &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait

  echo; echo "== Realign to HMMs =="
  mkdir -p 22_hmmalign
  for filepath in 21_hmm/*; do
    file=$(basename "$filepath");
    printf "%s " "$file"
    hmmalign --trim --outformat A2M --amino -o 22_hmmalign/"$file" 21_hmm/"$file" 19_pan_aug_leftover_merged_prot/"$file" &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo

  echo; echo "== Trim HMM alignments to match-states =="
  mkdir -p 23_hmmalign_trim1
  for filepath in 22_hmmalign/*; do 
    file=$(basename "$filepath");
    printf "%s " "$file"
    < "$filepath" perl -ne 'if ($_ =~ />/) {print $_} else {$line = $_; $line =~ s/[a-z]//g; print $line}' |
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
    printf "%s " "$file"
    filter_align.pl -in "$filepath" -out 23_hmmalign_trim2/"$file" -log 23_hmmalign_trim2_log/"$file" \
                    -depth $min_depth -pct_depth $min_pct_depth -min_pct_aligned $min_pct_aligned &
    if [[ $(jobs -r -p | wc -l) -ge $((NPROC/2)) ]]; then wait -n; fi
  done
  wait
  echo
}

##########
run_summarize() {
  echo; echo "Summarize: Move results into output directory, and report some summary statistics"

  cd "${WORK_DIR}" || exit
  echo "  work_dir: $PWD"

  echo "  Calculate matrix of gene counts per orthogroup and annotation set"
  calc_pan_stats.pl -annot_regex "$ANN_REX" -pan 18_syn_pan_aug_extra.clust.tsv -out 18_syn_pan_aug_extra.counts.tsv
 
  full_out_dir="${out_dir}"
  stats_file=${full_out_dir}/stats.txt

  cd "${submit_dir}" || exitd

  if [ ! -d "$full_out_dir" ]; then
      echo "creating output directory \"${full_out_dir}/\""
      mkdir -p "$full_out_dir"
  fi

  for file in 06_syn_pan.clust.tsv 06_syn_pan_ge3.hsh.tsv \
              12_syn_pan_aug.clust.tsv 12_syn_pan_aug.hsh.tsv \
              18_syn_pan_aug_extra.clust.tsv  18_syn_pan_aug_extra.hsh.tsv \
              18_syn_pan_aug_extra.table.tsv 18_syn_pan_aug_extra.counts.tsv; do
    if [ -f "${WORK_DIR}"/$file ]; then
      cp "${WORK_DIR}"/$file "${full_out_dir}"/
    else 
      echo "Warning: couldn't find file ${WORK_DIR}/$file; skipping"
    fi
  done

  echo "Copy manifest file into the output directory"
  if [ -f "${submit_dir}/manifests/MANIFEST_output_fam.yml" ]; then
    cp "${submit_dir}/manifests/MANIFEST_output_fam.yml" "$full_out_dir"/
  else
    echo "Couldn't find file manifests/MANIFEST_output_fam.yml"
  fi

  for dir in 19_pan_aug_leftover_merged_prot 21_hmm 22_hmmalign 23_hmmalign_trim2 24_trees; do
    if [ -d "${WORK_DIR}"/$dir ]; then
      echo "Copying directory $dir to output directory"
      cp -r "${WORK_DIR}"/$dir "${full_out_dir}"/
    else 
      echo "Warning: couldn't find dir ${WORK_DIR}/$dir; skipping"
    fi
  done

  printf "Run of program %s, version %s\n" "$scriptname" "$version" > "${stats_file}"

  end_time=$(date)
  cat "${WORK_DIR}"/stats/tmp.timing >> "${stats_file}"
  printf "Run ended at:   %s\n\n" "$end_time" >> "${stats_file}"

  echo "  Report parameters from config file"
  printf "Parameter  \tvalue\n" >> "${stats_file}"
  for key in ${pandagma_conf_params}; do
    printf '%-15s\t%s\n' "${key}" "${!key}" >> "${stats_file}"
  done

  printf "\nOutput directory for this run:\t%s\n" "$full_out_dir" >> "${stats_file}"

  echo "  Report orthogroup composition statistics for the three main cluster-calculation steps"

  echo "  Print sequence composition statistics for each annotation set"
  
  {
  printf "\n== Sequence stats for protein files\n" 
  printf "  Class:  seqs     min max    N50    ave     annotation_name\n"  
  if [ -f "${WORK_DIR}"/stats/tmp.fasta_seqstats ]; then
    cat "${WORK_DIR}"/stats/tmp.fasta_seqstats 

    printf "\n  Avg:   " >> "${stats_file}" 
    < "${WORK_DIR}"/stats/tmp.fasta_seqstats transpose.pl |
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

  cut -f2 "${WORK_DIR}"/18_syn_pan_aug_extra.hsh.tsv | 
    perl -pe '$ann_rex=qr($ENV{"ANN_REX"}); s/$ann_rex/$1/' |
    sort | uniq -c | awk '{print $2 "\t" $1}' > "${WORK_DIR}"/stats/tmp.gene_count_all_end

  paste "${WORK_DIR}"/stats/tmp.gene_count_start \
        "${WORK_DIR}"/stats/tmp.gene_count_all_end |
    awk 'BEGIN{print "  Start\tEnd\tPct_kept\tAnnotation_name"} 
        { printf "  %i\t%i\t%2.1f\t%s\n", $2, $4, 100*($4/$2), $1 }'  >> "${stats_file}"

  echo "  Print counts per accession"
  if [ -f "${full_out_dir}"/18_syn_pan_aug_extra.counts.tsv ]; then
    printf "\n== For all annotation sets, counts of genes-in-orthogroups and counts of orthogroups-with-genes:\n"
    printf "  gns-in-OGs  OGs-w-gns  OGs-w-gns/gns  pct-non-null-OGs  pct-null-OGs  annot-set\n" 
    transpose.pl "${full_out_dir}"/18_syn_pan_aug_extra.counts.tsv | 
      perl -lane 'next if ($.<=3); 
        $ct=0; $sum=0; $nulls=0; $OGs=0;
        for $i (@F[1..(@F-1)]){ $OGs++; if ($i>0){$ct++; $sum+=$i} if ($i==0){$nulls++} }; 
        printf("  %d\t%d\t%.2f\t%.2f\t%.2f\t%s\n", $sum, $ct, 100*$ct/$sum, 100*($OGs-$nulls)/$OGs, 100*$nulls/$OGs, $F[0])' \
      >> "${stats_file}"
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
  cd "${WORK_DIR}" || exit
  echo "  work_dir: $PWD"
  if [ -d MMTEMP ]; then rm -rf MMTEMP/*; 
  fi
  for dir in 11_pan_leftovers 13_extra_out_dir 16_pan_leftovers_extra 19_pan_aug_leftover_merged_prot; do
    if [ -d $dir ]; then echo "  Removing directory $dir"; rm -rf $dir &
    fi
  done
  #for file in 10* 11* 14* 20* 21* 23* 24* consen*; do
  for file in 10* 11* 14* 20*; do
    if [ -f "$file" ]; then echo "  Removing file $file"; rm "$file"; 
    fi
  done
  wait
  cd "$OLDPWD" || exit
}

########################################
# Main program

pandagma_conf_params='clust_iden clust_cov extra_iden mcl_inflation strict_synt 
ks_low_cutoff ks_hi_cutoff ks_binsize ks_block_wgd_cutoff max_pair_ks 
consen_prefix annot_str_regex'

export commandlist="ingest mmseqs filter dagchainer \
             ks_calc ks_filter \
             mcl consense cluster_rest add_extra tabularize \
             align model_and_trim calc_trees summarize"

export dependencies='dagchainer mcl cons famsa run_DAG_chainer.pl'

declare out_dir version clust_iden clust_cov extra_iden mcl_inflation strict_synt \
        ks_low_cutoff ks_hi_cutoff ks_binsize ks_block_wgd_cutoff max_pair_ks \
        consen_prefix annot_str_regex

main_pan_fam "$@"

#!/usr/bin/env bash
#
# This workflow, pandagma-TEfilter.sh, is for optional pre-filtering of CDS and protein data
# to exclude genes with matches to transposable elements or other sequences that the researcher wishes
# not to include in the pangene calculations.
# Authors: Steven Cannon, Joel Berendzen, Nathan Weeks, 2020-2023
#
scriptname='pandagma TEfilter'

. pandagma-common.sh

define HELP_DOC <<'EOS'
This workflow, pandagma-TEfilter.sh, is for optional pre-filtering of CDS and protein data
to exclude genes with matches to transposable elements or other sequences that the 
researcher wishes not to include in the pangene calculations. A nucleotide fasta file of 
sequences-to-avoid needs to be provided in the initial data directory; and another directory
needs to be indicated to receive the cleaned data.

Usage:
  $scriptname -c CONFIG_FILE [options]

  Required:
           -c (path to the config file)
           -i (path to initial data directory, with CDS, protein, annotation, and transposable element files)
           -o (path to output data directory)

  Options: -s (subcommand to run. If \"all\" or omitted, all steps will be run; otherwise, run specified step)
           -w (working directory, for temporary and intermediate files [default: './pandagma_work'].)
           -o OUTPUT_DIR (name for the output directory [default: './pandagma_out'].
                Applicable only to "all" and "summarize" steps.)
           -n (number of processors to use. Defaults to number of processors available to the parent process)
           -v (version)
           -h (help)
           -m (more information)

Environment requirements: The following packages need to be available in your PATH:
    mmseqs

Also, please add the pandagma utility programs in the bin directory adjacent to pandagma-fam.sh, e.g.
    PATH=$PWD/bin:\$PATH

Subcommands (in order they are usually run):
                all - All of the steps below
                        (Or equivalently: omit the -s flag; \"all\" is default).
         TEfilter - Search and remove matches to transposable elements or other sequences in the provided exclude_TE_match_file
EOS

########################################
# Helper functions begin here

canonicalize_paths() {
  echo "Entering canonicalize_paths. Fasta files: "
  echo "${cds_files[@]}"

  mkdir -p ${DATA_DIR}_TEfilt
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
  echo "CHECK exclude_TE_match_file: $exclude_TE_match_file"
  if [[ -v exclude_TE_match_file ]]
  then
    readonly exclude_TE_match_filepath=${exclude_TE_match_file:+$(realpath "${exclude_TE_match_file}")}
    echo "exclude_TE_match_file: $exclude_TE_match_file"
    echo "exclude_TE_match_filepath: $exclude_TE_match_filepath"
  fi

  cd "${OLDPWD}" || exit
  readonly submit_dir=${PWD}

  cds_fasta_file=$(basename "${cds_files[0]}" .gz)
  fna="${cds_fasta_file##*.}"

  prot_fasta_file=$(basename "${protein_files[0]}" .gz)
  faa="${prot_fasta_file##*.}"
}

########################################
# run functions 

##########
run_TEfilter() {
  # Search CDS files against file of sequences to exclude, e.g. transposable elements and plastid sequences
  cd "${WORK_DIR}" || exit
  echo; echo "Run mmseqs search against ${exclude_TE_match_file} -- at ${TE_match_iden} percent identity & ${clust_cov} query coverage."

  if [[ -v exclude_TE_match_file ]]
  then
    echo "  Copying exclude_TE_match_file ${exclude_TE_match_filepath}"
    mkdir -p 00_exclude_TE_match_file
    exclude_TE_match_base=$(basename "$exclude_TE_match_file" ".$fna.gz")
    cat_or_zcat "${exclude_TE_match_filepath}" > 00_exclude_TE_match_file/"${exclude_TE_match_base}.$fna"
  fi

  mkdir -p 03_mmseqs_tmp
  for dir in 00_TEsearch 00_fasta_nuc_orig 00_fasta_prot_orig ; do
    if [ -d $dir ]; then rm -rf $dir ; fi
    mkdir -p $dir 
  done
  exclude_TE_match_base=$(basename "$exclude_TE_match_file" ".$fna.gz")

  MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)

  ##########
  echo "Handle cds_files and corresponding protein files"
  for (( file1_num = 0; file1_num < ${#cds_files[@]} ; file1_num++ )); do
    cds_base=$(basename "${cds_files[file1_num]%.*}" ".$fna")
    prot_base=$(basename "${protein_files[file1_num]%.*}" ".$faa")
    annot_file=$(basename "${annotation_files[file1_num]}")
    echo; echo "  Running mmseqs on comparison: ${cds_base}.x.${exclude_TE_match_base}"

    cat_or_zcat "${cds_files[file1_num]}" > 00_fasta_nuc_orig/"$cds_base.$fna" # Pull file locally; uncompressed needed later
    cat_or_zcat "${protein_files[file1_num]}" > 00_fasta_prot_orig/"$prot_base.$faa" # Pull file locally; uncompressed needed later
      mmseqs easy-search "00_fasta_nuc_orig/${cds_base}.$fna" \
                 00_exclude_TE_match_file/"$exclude_TE_match_base"."$fna" \
                 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8 \
                 "$MMTEMP" --min-seq-id "$TE_match_iden" --search-type 3 --cov-mode 2 -c "${clust_cov}" 1>/dev/null

    echo "  Get list of top matches for each 00_TEsearch file $cds_base.x.$exclude_TE_match_base.m8"
    top_line.awk 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8 > 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8_top
    cut -f1 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8_top > 00_TEsearch/lis."$cds_base"

    echo "  Generate CDS and protein files, excluding matches from the TEsearch"
    get_fasta_subset.pl -in 00_fasta_nuc_orig/"$cds_base.$fna" -clobber \
                         -xclude -lis 00_TEsearch/lis."$cds_base" -out ${DATA_DIR}_TEfilt/"$cds_base.$fna"
    get_fasta_subset.pl -in 00_fasta_prot_orig/"$prot_base.$faa" -clobber \
                         -xclude -lis 00_TEsearch/lis."$cds_base" -out ${DATA_DIR}_TEfilt/"$prot_base.$faa"
    cp ${DATA_DIR}/$annot_file ${DATA_DIR}_TEfilt/$annot_file
    gzip -f ${DATA_DIR}_TEfilt/"$cds_base.$fna" ${DATA_DIR}_TEfilt/"$prot_base.$faa"
  done

  ##########
  if [[ -v cds_files_extra_constr ]]
  then
    echo "Handle cds_files_extra_constr and corresponding protein files"
    for (( file1_num = 0; file1_num < ${#cds_files_extra_constr[@]} ; file1_num++ )); do
      cds_base=$(basename "${cds_files_extra_constr[file1_num]%.*}" ".$fna")
      prot_base=$(basename "${protein_files_extra_constr[file1_num]%.*}" ".$faa")
      annot_file=$(basename "${annotation_files_extra_constr[file1_num]}")
      echo; echo "  Running mmseqs on comparison: ${cds_base}.x.${exclude_TE_match_base}"

      cat_or_zcat "${cds_files_extra_constr[file1_num]}" > 00_fasta_nuc_orig/"$cds_base.$fna" # Pull uncompressed file locally
      cat_or_zcat "${protein_files_extra_constr[file1_num]}" > 00_fasta_prot_orig/"$prot_base.$faa" # Pull uncompressed file locally
        mmseqs easy-search "00_fasta_nuc_orig/${cds_base}.$fna" \
                   00_exclude_TE_match_file/"$exclude_TE_match_base"."$fna" \
                   00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8 \
                   "$MMTEMP" --min-seq-id "$TE_match_iden" --search-type 3 --cov-mode 2 -c "${clust_cov}" 1>/dev/null

      echo "  Get list of top matches for each 00_TEsearch file $cds_base.x.$exclude_TE_match_base.m8"
      top_line.awk 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8 > 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8_top
      cut -f1 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8_top > 00_TEsearch/lis."$cds_base"

      echo "  Generate CDS and protein files, excluding matches from the TEsearch"
      get_fasta_subset.pl -in 00_fasta_nuc_orig/"$cds_base.$fna" -clobber \
                           -xclude -lis 00_TEsearch/lis."$cds_base" -out ${DATA_DIR}_TEfilt/"$cds_base.$fna"
      get_fasta_subset.pl -in 00_fasta_prot_orig/"$prot_base.$faa" -clobber \
                           -xclude -lis 00_TEsearch/lis."$cds_base" -out ${DATA_DIR}_TEfilt/"$prot_base.$faa"
      cp ${DATA_DIR}/$annot_file ${DATA_DIR}_TEfilt/$annot_file
      gzip -f ${DATA_DIR}_TEfilt/"$cds_base.$fna" ${DATA_DIR}_TEfilt/"$prot_base.$faa"
    done
  fi

  ##########
  if [[ -v cds_files_extra_free ]]
  then
    echo "Handle cds_files_extra_free and corresponding protein files"
    for (( file1_num = 0; file1_num < ${#cds_files_extra_free[@]} ; file1_num++ )); do
      cds_base=$(basename "${cds_files_extra_free[file1_num]%.*}" ".$fna")
      prot_base=$(basename "${protein_files_extra_free[file1_num]%.*}" ".$faa")
      annot_file=$(basename "${annotation_files_extra_free[file1_num]}")
      echo; echo "  Running mmseqs on comparison: ${cds_base}.x.${exclude_TE_match_base}"

      cat_or_zcat "${cds_files_extra_free[file1_num]}" > 00_fasta_nuc_orig/"$cds_base.$fna" # Pull uncompressed file locally
      cat_or_zcat "${protein_files_extra_free[file1_num]}" > 00_fasta_prot_orig/"$prot_base.$faa" # Pull uncompressed file locally
        mmseqs easy-search "00_fasta_nuc_orig/${cds_base}.$fna" \
                   00_exclude_TE_match_file/"$exclude_TE_match_base"."$fna" \
                   00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8 \
                   "$MMTEMP" --min-seq-id "$TE_match_iden" --search-type 3 --cov-mode 2 -c "${clust_cov}" 1>/dev/null

      echo "  Generate CDS and protein files, excluding matches from the TEsearch"
      top_line.awk 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8 > 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8_top
      cut -f1 00_TEsearch/"$cds_base".x."$exclude_TE_match_base".m8_top > 00_TEsearch/lis."$cds_base"

      get_fasta_subset.pl -in 00_fasta_nuc_orig/"$cds_base.$fna" -clobber \
                           -xclude -lis 00_TEsearch/lis."$cds_base" -out ${DATA_DIR}_TEfilt/"$cds_base.$fna"
      get_fasta_subset.pl -in 00_fasta_prot_orig/"$prot_base.$faa" -clobber \
                           -xclude -lis 00_TEsearch/lis."$cds_base" -out ${DATA_DIR}_TEfilt/"$prot_base.$faa"
      cp ${DATA_DIR}/$annot_file ${DATA_DIR}_TEfilt/$annot_file
      gzip -f ${DATA_DIR}_TEfilt/"$cds_base.$fna" ${DATA_DIR}_TEfilt/"$prot_base.$faa"
    done
  fi
}

########################################
# Main program

pandagma_conf_params='clust_cov '

# Run all specified steps 
commandlist="TEfilter"

dependencies='mmseqs'

declare  clust_cov 

main_pan_fam "$@"

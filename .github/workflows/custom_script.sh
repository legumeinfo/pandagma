#!/usr/bin/env bash
#
# Example custom script to supply to GitHub Action workflow:
# .github/workflows/custom_script.yml
#
# See pandagma/README.md "Workflow Automation with GitHub Actions"

set -o errexit -o nounset -o pipefail -o xtrace

# download data
mkdir data

CURL='curl --fail --no-progress-meter -LO'

cd data

${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm1.ann1.nRhs/cicar.CDCFrontier.gnm1.ann1.nRhs.cds.fna.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm1.ann1.nRhs/cicar.CDCFrontier.gnm1.ann1.nRhs.gene_models_main.bed.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm1.ann1.nRhs/cicar.CDCFrontier.gnm1.ann1.nRhs.protein.faa.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm2.ann1.9M1L/cicar.CDCFrontier.gnm2.ann1.9M1L.cds_primary.fna.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm2.ann1.9M1L/cicar.CDCFrontier.gnm2.ann1.9M1L.gene_models_main.bed.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm2.ann1.9M1L/cicar.CDCFrontier.gnm2.ann1.9M1L.protein_primary.faa.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7/cicar.CDCFrontier.gnm3.ann1.NPD7.cds.fna.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7/cicar.CDCFrontier.gnm3.ann1.NPD7.gene_models_main.bed.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7/cicar.CDCFrontier.gnm3.ann1.NPD7.protein.faa.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/ICC4958.gnm2.ann1.LCVX/cicar.ICC4958.gnm2.ann1.LCVX.cds_primary.fna.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/ICC4958.gnm2.ann1.LCVX/cicar.ICC4958.gnm2.ann1.LCVX.gene_models_main.bed.gz
${CURL} https://data.legumeinfo.org/Cicer/arietinum/annotations/ICC4958.gnm2.ann1.LCVX/cicar.ICC4958.gnm2.ann1.LCVX.protein_primary.faa.gz
${CURL} https://data.legumeinfo.org/Cicer/echinospermum/annotations/S2Drd065.gnm1.ann1.YZ9H/cicec.S2Drd065.gnm1.ann1.YZ9H.cds.fna.gz
${CURL} https://data.legumeinfo.org/Cicer/echinospermum/annotations/S2Drd065.gnm1.ann1.YZ9H/cicec.S2Drd065.gnm1.ann1.YZ9H.gene_models_main.bed.gz
${CURL} https://data.legumeinfo.org/Cicer/echinospermum/annotations/S2Drd065.gnm1.ann1.YZ9H/cicec.S2Drd065.gnm1.ann1.YZ9H.protein.faa.gz
${CURL} https://data.legumeinfo.org/Cicer/reticulatum/annotations/Besev079.gnm1.ann1.F01Z/cicre.Besev079.gnm1.ann1.F01Z.cds.fna.gz
${CURL} https://data.legumeinfo.org/Cicer/reticulatum/annotations/Besev079.gnm1.ann1.F01Z/cicre.Besev079.gnm1.ann1.F01Z.cds.fna.gz
${CURL} https://data.legumeinfo.org/Cicer/reticulatum/annotations/Besev079.gnm1.ann1.F01Z/cicre.Besev079.gnm1.ann1.F01Z.gene_models_main.bed.gz
${CURL} https://data.legumeinfo.org/Cicer/reticulatum/annotations/Besev079.gnm1.ann1.F01Z/cicre.Besev079.gnm1.ann1.F01Z.protein.faa.gz

cd ..

# create config

cat > pandagma_custom.conf <<'END'
clust_iden='0.90'
clust_cov="0.50"
extra_iden='0.80'
TE_match_iden='0.40'
mcl_inflation='1.6'
strict_synt="1"
pctl_low="25"
pctl_med="50"
pctl_hi="75"
consen_prefix='Cicer.pan2'
annot_str_regex='([^.]+\.[^.]+\.[^.]+\.[^.]+)\..+'
preferred_annot='CDCFrontier.gnm3.ann1'
order_method="reference"

##### (required) list of GFF & FASTA file paths
# Uncomment add file paths to the the annotation_files and cds_files arrays.
# The nth listed GFF file corresponds to the nth listed FASTA file.

annotation_files=(
  cicar.CDCFrontier.gnm3.ann1.NPD7.gene_models_main.bed.gz
  cicar.ICC4958.gnm2.ann1.LCVX.gene_models_main.bed.gz
  cicec.S2Drd065.gnm1.ann1.YZ9H.gene_models_main.bed.gz
  cicre.Besev079.gnm1.ann1.F01Z.gene_models_main.bed.gz
)

cds_files=(
  cicar.CDCFrontier.gnm3.ann1.NPD7.cds.fna.gz
  cicar.ICC4958.gnm2.ann1.LCVX.cds_primary.fna.gz
  cicec.S2Drd065.gnm1.ann1.YZ9H.cds.fna.gz
  cicre.Besev079.gnm1.ann1.F01Z.cds.fna.gz
)

protein_files=(
  cicar.CDCFrontier.gnm3.ann1.NPD7.protein.faa.gz
  cicar.ICC4958.gnm2.ann1.LCVX.protein_primary.faa.gz
  cicec.S2Drd065.gnm1.ann1.YZ9H.protein.faa.gz
  cicre.Besev079.gnm1.ann1.F01Z.protein.faa.gz
)

### (optional) Extra GFF & FASTA files
annotation_files_extra_constr=(
  cicar.CDCFrontier.gnm1.ann1.nRhs.gene_models_main.bed.gz
  cicar.CDCFrontier.gnm2.ann1.9M1L.gene_models_main.bed.gz
)

cds_files_extra_constr=(
  cicar.CDCFrontier.gnm1.ann1.nRhs.cds.fna.gz
  cicar.CDCFrontier.gnm2.ann1.9M1L.cds_primary.fna.gz
)

protein_files_extra_constr=(
  cicar.CDCFrontier.gnm1.ann1.nRhs.protein.faa.gz
  cicar.CDCFrontier.gnm2.ann1.9M1L.protein_primary.faa.gz
)


#### (optional) expected_chr_matches
expected_chr_matches=(
  1 1
  2 2
  3 3
  4 4
  5 5
  6 6
  7 7
  8 8
)
END

# must be executable to be sourced
chmod +x pandagma_custom.conf

# Run pandagma

pandagma pan -c ./pandagma_custom.conf -d data


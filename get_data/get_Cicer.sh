# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data ]; then mkdir -p data; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org/Cicer"

## data
base_dir=$PWD
cd $base_dir/data/

curl -O $url_base/arietinum/annotations/CDCFrontier.gnm1.ann1.nRhs/cicar.CDCFrontier.gnm1.ann1.nRhs.gene_models_main.bed.gz
curl -O $url_base/arietinum/annotations/CDCFrontier.gnm2.ann1.9M1L/cicar.CDCFrontier.gnm2.ann1.9M1L.gene_models_main.bed.gz
curl -O $url_base/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7/cicar.CDCFrontier.gnm3.ann1.NPD7.gene_models_main.bed.gz
curl -O $url_base/arietinum/annotations/ICC4958.gnm2.ann1.LCVX/cicar.ICC4958.gnm2.ann1.LCVX.gene_models_main.bed.gz
curl -O $url_base/echinospermum/annotations/S2Drd065.gnm1.ann1.YZ9H/cicec.S2Drd065.gnm1.ann1.YZ9H.gene_models_main.bed.gz
curl -O $url_base/reticulatum/annotations/Besev079.gnm1.ann1.F01Z/cicre.Besev079.gnm1.ann1.F01Z.gene_models_main.bed.gz

curl -O $url_base/arietinum/annotations/CDCFrontier.gnm1.ann1.nRhs/cicar.CDCFrontier.gnm1.ann1.nRhs.cds.fna.gz
curl -O $url_base/arietinum/annotations/CDCFrontier.gnm2.ann1.9M1L/cicar.CDCFrontier.gnm2.ann1.9M1L.cds_primary.fna.gz
curl -O $url_base/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7/cicar.CDCFrontier.gnm3.ann1.NPD7.cds.fna.gz
curl -O $url_base/arietinum/annotations/ICC4958.gnm2.ann1.LCVX/cicar.ICC4958.gnm2.ann1.LCVX.cds_primary.fna.gz
curl -O $url_base/echinospermum/annotations/S2Drd065.gnm1.ann1.YZ9H/cicec.S2Drd065.gnm1.ann1.YZ9H.cds.fna.gz
curl -O $url_base/reticulatum/annotations/Besev079.gnm1.ann1.F01Z/cicre.Besev079.gnm1.ann1.F01Z.cds.fna.gz

curl -O $url_base/arietinum/annotations/CDCFrontier.gnm1.ann1.nRhs/cicar.CDCFrontier.gnm1.ann1.nRhs.protein.faa.gz
curl -O $url_base/arietinum/annotations/CDCFrontier.gnm2.ann1.9M1L/cicar.CDCFrontier.gnm2.ann1.9M1L.protein_primary.faa.gz
curl -O $url_base/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7/cicar.CDCFrontier.gnm3.ann1.NPD7.protein.faa.gz
curl -O $url_base/arietinum/annotations/ICC4958.gnm2.ann1.LCVX/cicar.ICC4958.gnm2.ann1.LCVX.protein_primary.faa.gz
curl -O $url_base/echinospermum/annotations/S2Drd065.gnm1.ann1.YZ9H/cicec.S2Drd065.gnm1.ann1.YZ9H.protein.faa.gz
curl -O $url_base/reticulatum/annotations/Besev079.gnm1.ann1.F01Z/cicre.Besev079.gnm1.ann1.F01Z.protein.faa.gz

# Fix the forms of chromosome IDs, if necessary
  # Not necessary for Cicer
#  for file in *; do gunzip $file & done

# Check chromosome names
  for file in *bed.gz; do zcat $file | awk '$1~/chr|Chr|Ca/ {print $1}' | sort | uniq -c ; echo; done

  for file in *bed.gz; do echo $file; zcat $file | awk '$1~/chr|Chr|Ca/ {print $1}' | sort | uniq ; echo; done

# Re-compress the files
#  for file in *bed *fna *faa; do gzip $file &
#  done

cat <<DATA > expected_chr_matches.tsv
# Expected chromosome matches for Cicer
1 1
2 2
3 3
4 4
5 5
6 6
7 7
8 8
DATA

cd $base_dir



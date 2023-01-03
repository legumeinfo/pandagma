# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data ]; then mkdir -p data; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org/Phaseolus"

## data
base_dir=$PWD
cd $base_dir/data/
  curl -O $url_base/acutifolius/annotations/Frijol_Bayo.gnm1.ann1.ML22/phaac.Frijol_Bayo.gnm1.ann1.ML22.cds_primary.fna.gz
  curl -O $url_base/acutifolius/annotations/W6_15578.gnm2.ann1.LVZ1/phaac.W6_15578.gnm2.ann1.LVZ1.cds_primary.fna.gz
  curl -O $url_base/lunatus/annotations/G27455.gnm1.ann1.JD7C/phalu.G27455.gnm1.ann1.JD7C.cds_primary.fna.gz
  curl -O $url_base/vulgaris/annotations/5-593.gnm1.ann1.3FBJ/phavu.5-593.gnm1.ann1.3FBJ.cds_primary.fna.gz
  curl -O $url_base/vulgaris/annotations/G19833.gnm1.ann1.pScz/phavu.G19833.gnm1.ann1.pScz.cds_primary.fna.gz
  curl -O $url_base/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.cds_primary.fna.gz
  curl -O $url_base/vulgaris/annotations/LaborOvalle.gnm1.ann1.L1DY/phavu.LaborOvalle.gnm1.ann1.L1DY.cds_primary.fna.gz
  curl -O $url_base/vulgaris/annotations/UI111.gnm1.ann1.8L4N/phavu.UI111.gnm1.ann1.8L4N.cds_primary.fna.gz

  curl -O $url_base/acutifolius/annotations/Frijol_Bayo.gnm1.ann1.ML22/phaac.Frijol_Bayo.gnm1.ann1.ML22.cds.bed.gz
  curl -O $url_base/acutifolius/annotations/W6_15578.gnm2.ann1.LVZ1/phaac.W6_15578.gnm2.ann1.LVZ1.cds.bed.gz
  curl -O $url_base/lunatus/annotations/G27455.gnm1.ann1.JD7C/phalu.G27455.gnm1.ann1.JD7C.cds.bed.gz
  curl -O $url_base/vulgaris/annotations/5-593.gnm1.ann1.3FBJ/phavu.5-593.gnm1.ann1.3FBJ.cds.bed.gz
  curl -O $url_base/vulgaris/annotations/G19833.gnm1.ann1.pScz/phavu.G19833.gnm1.ann1.pScz.cds.bed.gz
  curl -O $url_base/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.cds.bed.gz
  curl -O $url_base/vulgaris/annotations/LaborOvalle.gnm1.ann1.L1DY/phavu.LaborOvalle.gnm1.ann1.L1DY.cds.bed.gz
  curl -O $url_base/vulgaris/annotations/UI111.gnm1.ann1.8L4N/phavu.UI111.gnm1.ann1.8L4N.cds.bed.gz

  curl -O $url_base/acutifolius/annotations/Frijol_Bayo.gnm1.ann1.ML22/phaac.Frijol_Bayo.gnm1.ann1.ML22.protein_primary.faa.gz
  curl -O $url_base/acutifolius/annotations/W6_15578.gnm2.ann1.LVZ1/phaac.W6_15578.gnm2.ann1.LVZ1.protein_primary.faa.gz
  curl -O $url_base/lunatus/annotations/G27455.gnm1.ann1.JD7C/phalu.G27455.gnm1.ann1.JD7C.protein_primary.faa.gz
  curl -O $url_base/vulgaris/annotations/5-593.gnm1.ann1.3FBJ/phavu.5-593.gnm1.ann1.3FBJ.protein_primary.faa.gz
  curl -O $url_base/vulgaris/annotations/G19833.gnm1.ann1.pScz/phavu.G19833.gnm1.ann1.pScz.protein_primary.faa.gz
  curl -O $url_base/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.protein_primary.faa.gz
  curl -O $url_base/vulgaris/annotations/LaborOvalle.gnm1.ann1.L1DY/phavu.LaborOvalle.gnm1.ann1.L1DY.protein_primary.faa.gz
  curl -O $url_base/vulgaris/annotations/UI111.gnm1.ann1.8L4N/phavu.UI111.gnm1.ann1.8L4N.protein_primary.faa.gz


# Fix the forms of chromosome IDs, if necessary
  # Not necessary for Phaseolus
#  for file in *; do gunzip $file & done

# Check chromosome names
  for file in *bed.gz; do zcat $file | awk '$1~/Chr|Pl/i {print $1}' | sort | uniq -c ; echo; done

  for file in *bed.gz; do echo $file; zcat $file | awk '$1~/Chr|Pl/i {print $1}' | sort | uniq ; echo; done

# Re-compress the files
#  for file in *bed *fna *faa; do gzip $file &
#  done

cat <<DATA > expected_chr_matches.tsv
# Expected chromosome matches for Phaseolus vulgaris
01 01
02 02
03 03
04 04
05 05
06 06
07 07
08 08
09 09
10 10
11 11
DATA

cd $base_dir



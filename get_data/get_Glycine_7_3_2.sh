# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset
 
if [ ! -d data_pan ]; then mkdir -p data_pan; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org/Glycine"

## data
base_dir=$PWD
cd $base_dir/data_pan/

# CDS
curl -f -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/JD17.gnm1.ann1.CLFP/glyma.JD17.gnm1.ann1.CLFP.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.cds.fna.gz
curl -f -O $url_base/max/annotations/PI_398296.gnm1.ann1.B0XR/glyma.PI_398296.gnm1.ann1.B0XR.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/PI_548362.gnm1.ann1.LL84/glyma.PI_548362.gnm1.ann1.LL84.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.cds.fna.gz
curl -f -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.cds_primary.fna.gz
curl -f -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.cds_primary.fna.gz
# 
# BED or GFF
curl -f -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.gene_models_main.gff3.gz
curl -f -O $url_base/max/annotations/JD17.gnm1.ann1.CLFP/glyma.JD17.gnm1.ann1.CLFP.gene_models_main.gff3.gz
curl -f -O $url_base//max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.gene_models_main.gff3.gz
curl -f -O $url_base/max/annotations/PI_398296.gnm1.ann1.B0XR/glyma.PI_398296.gnm1.ann1.B0XR.gene_models_main.gff3.gz
curl -f -O $url_base/max/annotations/PI_548362.gnm1.ann1.LL84/glyma.PI_548362.gnm1.ann1.LL84.gene_models_main.gff3.gz
curl -f -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.gff3.gz
curl -f -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz
curl -f -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.gff3.gz
curl -f -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.gene_models_main.gff3.gz
curl -f -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.gene_models_main.gff3.gz
curl -f -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.gene_models_main.gff3.gz
curl -f -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.gene_models_main.gff3.gz

# Protein
curl -f -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/JD17.gnm1.ann1.CLFP/glyma.JD17.gnm1.ann1.CLFP.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.protein.faa.gz
curl -f -O $url_base/max/annotations/PI_398296.gnm1.ann1.B0XR/glyma.PI_398296.gnm1.ann1.B0XR.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/PI_548362.gnm1.ann1.LL84/glyma.PI_548362.gnm1.ann1.LL84.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.protein.faa.gz
curl -f -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.protein_primary.faa.gz
curl -f -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.protein_primary.faa.gz

# Convert from gff3 to bed7
for file in *gff3.gz; do
  base=`basename $file .gff3.gz`
  echo "Deriving a BED7 file from $base.gff3"
  zcat $file | sort_gff.pl | gff_to_bed7_mRNA.awk > $base.bed
  gzip $base.bed
done

# Patch glyso.W05.gnm1.ann1 annotations: the CDS and protein files have Glysoja.10G027808, but the BED does not.
zcat glyso.W05.gnm1.ann1.T47J.protein_primary.faa.gz | fasta_to_table.awk | awk '$1!~/Glysoja.10G027808/' | 
  awk -v FS="\t" '{print ">" $1; print $2}' > tmp.prot_primary.faa
  mv tmp.prot_primary.faa glyso.W05.gnm1.ann1.T47J.protein_primary.faa
  gzip -f glyso.W05.gnm1.ann1.T47J.protein_primary.faa

zcat glyso.W05.gnm1.ann1.T47J.cds_primary.fna.gz | fasta_to_table.awk | awk '$1!~/Glysoja.10G027808/' |
  awk -v FS="\t" '{print ">" $1; print $2}'  > tmp.cds_primary.fna
  mv tmp.cds_primary.fna glyso.W05.gnm1.ann1.T47J.cds_primary.fna
  gzip -f glyso.W05.gnm1.ann1.T47J.cds_primary.fna
  

# Write file expected_chr_matches.tsv
cat <<DATA > expected_chr_matches.tsv
# Expected chromosome matches for Glycine soja and Glycine max
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
11 13
11 11
12 12
13 13
14 14
15 15
16 16
17 17
18 18
19 19
20 20
DATA

cd $base_dir


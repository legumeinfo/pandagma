# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data ]; then mkdir -p data; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org"

## data
base_dir=$PWD
cd $base_dir/data/

curl -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.cds.bed.gz

curl -O $url_base/Arachis/GENUS/pangenes/Arachis.pan1.4LN9/Arachis.pan1.4LN9.pctl25_named_protein.faa.gz
curl -O $url_base/Cicer/GENUS/pangenes/Cicer.pan1.SV8C/Cicer.pan1.SV8C.pctl25_named_protein.faa.gz
curl -O $url_base/Glycine/GENUS/pangenes/Glycine.pan3.YWTW/Glycine.pan3.YWTW.pctl25_named_protein.faa.gz
curl -O $url_base/Medicago/GENUS/pangenes/Medicago.pan1.XXQ6/Medicago.pan1.XXQ6.pctl25_named_protein.faa.gz
curl -O $url_base/Phaseolus/GENUS/pangenes/Phaseolus.pan1.X2PC/Phaseolus.pan1.X2PC.pctl25_named_protein.faa.gz
curl -O $url_base/Vigna/GENUS/pangenes/Vigna.pan1.X2PC/Vigna.pan1.X2PC.pctl25_named_protein.faa.gz

curl -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.protein_primary.faa.gz

# Rename CDS-based bed files to protein
for file in *.cds.bed.gz; do
  base=`basename $file .cds.bed.gz`
  mv $file $base.protein.bed.gz
done

# From the pan-gene protein files, derive bed files
for file in *_named_protein.faa.gz; do 
  base=`basename $file .pctl25_named_protein.faa.gz`
  zcat $file | awk -v OFS="\t" '$1~/^>/ {print substr($1,2), $2, $3, substr($1,2), 0, $4}' | 
    perl -pe 's/^(\w+\.\w+\.[^_]+)_\S+/$1/' > $base.pctl25_named_protein.bed
done

# Check chromosome names
  for file in *bed.gz; do zcat $file | awk '$1~/chr|Chr|Ca/ {print $1}' | sort | uniq -c ; echo; done

  for file in *bed.gz; do echo $file; zcat $file | awk '$1~/chr|Chr|Ca/ {print $1}' | sort | uniq ; echo; done

# Re-compress the files
for file in *bed; do gzip $file &
done


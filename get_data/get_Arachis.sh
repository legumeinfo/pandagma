# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

mkdir -p data data_orig

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org/Arachis"

## data
base_dir=$PWD
cd $base_dir/data_orig/
curl -O $url_base/duranensis/annotations/V14167.gnm1.ann1.cxSM/aradu.V14167.gnm1.ann1.cxSM.cds.fna.gz
curl -O $url_base/hypogaea/annotations/Tifrunner.gnm1.ann1.CCJH/arahy.Tifrunner.gnm1.ann1.CCJH.cds_primary.fna.gz
curl -O $url_base/hypogaea/annotations/Tifrunner.gnm2.ann1.4K0L/arahy.Tifrunner.gnm2.ann1.4K0L.cds_primary.fna.gz
curl -O $url_base/ipaensis/annotations/K30076.gnm1.ann1.J37m/araip.K30076.gnm1.ann1.J37m.cds.fna.gz

curl -O $url_base/duranensis/annotations/V14167.gnm1.ann1.cxSM/aradu.V14167.gnm1.ann1.cxSM.cds.bed.gz
curl -O $url_base/hypogaea/annotations/Tifrunner.gnm1.ann1.CCJH/arahy.Tifrunner.gnm1.ann1.CCJH.cds.bed.gz
curl -O $url_base/hypogaea/annotations/Tifrunner.gnm2.ann1.4K0L/arahy.Tifrunner.gnm2.ann1.4K0L.cds.bed.gz
curl -O $url_base/ipaensis/annotations/K30076.gnm1.ann1.J37m/araip.K30076.gnm1.ann1.J37m.cds.bed.gz

curl -O $url_base/duranensis/annotations/V14167.gnm1.ann1.cxSM/aradu.V14167.gnm1.ann1.cxSM.protein.faa.gz
curl -O $url_base/hypogaea/annotations/Tifrunner.gnm1.ann1.CCJH/arahy.Tifrunner.gnm1.ann1.CCJH.protein_primary.faa.gz
curl -O $url_base/hypogaea/annotations/Tifrunner.gnm2.ann1.4K0L/arahy.Tifrunner.gnm2.ann1.4K0L.protein_primary.faa.gz
curl -O $url_base/ipaensis/annotations/K30076.gnm1.ann1.J37m/araip.K30076.gnm1.ann1.J37m.protein.faa.gz

# Merge duranensis and ipaensis to give a pseudo-allotetraploid that can be compared with A. hypogaea
zcat aradu.V14167.gnm1.ann1.cxSM.cds.bed.gz araip.K30076.gnm1.ann1.J37m.cds.bed.gz |
  perl -pe 's/aradu.V14167.gnm1.Adur/araduip.V14167K30076.gnm1.scaffA_/; 
            s/araip.K30076.gnm1.Aipa/araduip.V14167K30076.gnm1.scaffB_/; 
            s/aradu.V14167.gnm1.Aradu.A/araduip.V14167K30076.gnm1.chrA/;
            s/^araip.K30076.gnm1.Araip.B0(\d)/araduip.V14167K30076.gnm1.chrB1$1/; 
            s/^araip.K30076.gnm1.Araip.B1(\d)/araduip.V14167K30076.gnm1.chrB2$1/' |
  cat > ../data/araduip.V14167K30076.gnm1.ann1.cxSMJ37m.cds.bed

# Tweak chromosome prefixes in hypogaea to make them distinguishable in regex vs. ipaensis and duranensis
# and tweak the chrom names, in order to track distinct annotations via gnm#.ann#
  for filepath in arahy.Tifrunner.*.bed.gz; do
    base=`basename $filepath .gz`
    echo $base
    zcat $filepath | perl -pe 's/Arahy.(\d\d)/Arahy.Ah$1/' > ../data/$base
  done

# The cds and protein sequence IDs should be OK without modification, but the dur and ipa need to be merged.
  zcat aradu.V14167.gnm1.ann1.cxSM.cds.fna.gz araip.K30076.gnm1.ann1.J37m.cds.fna.gz |
    cat > ../data/araduip.V14167K30076.gnm1.ann1.cxSMJ37m.cds.fna

  zcat aradu.V14167.gnm1.ann1.cxSM.protein.faa.gz araip.K30076.gnm1.ann1.J37m.protein.faa.gz |
    cat > ../data/araduip.V14167K30076.gnm1.ann1.cxSMJ37m.protein.faa

  zcat arahy.Tifrunner.gnm1.ann1.CCJH.cds_primary.fna.gz > ../data/arahy.Tifrunner.gnm1.ann1.CCJH.cds_primary.fna
  zcat arahy.Tifrunner.gnm1.ann1.CCJH.protein_primary.faa.gz > ../data/arahy.Tifrunner.gnm1.ann1.CCJH.protein_primary.faa

  zcat arahy.Tifrunner.gnm2.ann1.4K0L.cds_primary.fna.gz > ../data/arahy.Tifrunner.gnm2.ann1.4K0L.cds_primary.fna
  zcat arahy.Tifrunner.gnm2.ann1.4K0L.protein_primary.faa.gz > ../data/arahy.Tifrunner.gnm2.ann1.4K0L.protein_primary.faa

cd $base_dir/data/

for file in *; do gzip $file &
done

cat <<DATA > expected_chr_matches.tsv
# Expected chromosome matches for Arachis
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



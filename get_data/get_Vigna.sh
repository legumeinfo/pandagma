# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data ]; then mkdir -p data; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://legumeinfo.org/data/v2/Vigna"

## data
base_dir=$PWD
cd $base_dir/data/
  curl -O $url_base/angularis/annotations/Gyeongwon.gnm3.ann1.3Nz5/vigan.Gyeongwon.gnm3.ann1.3Nz5.cds_primary.fna.gz
  curl -O $url_base/angularis/annotations/Shumari.gnm1.ann1.8BRS/vigan.Shumari.gnm1.ann1.8BRS.cds_primary.fna.gz
  curl -O $url_base/radiata/annotations/VC1973A.gnm6.ann1.M1Qs/vigra.VC1973A.gnm6.ann1.M1Qs.cds.fna.gz
  curl -O $url_base/unguiculata/annotations/CB5-2.gnm1.ann1.0GKC/vigun.CB5-2.gnm1.ann1.0GKC.cds_primary.fna.gz
  curl -O $url_base/unguiculata/annotations/IT97K-499-35.gnm1.ann1.zb5D/vigun.IT97K-499-35.gnm1.ann1.zb5D.cds_primary.fna.gz
  curl -O $url_base/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.cds_primary.fna.gz
  curl -O $url_base/unguiculata/annotations/Sanzi.gnm1.ann1.HFH8/vigun.Sanzi.gnm1.ann1.HFH8.cds_primary.fna.gz
  curl -O $url_base/unguiculata/annotations/Suvita2.gnm1.ann1.1PF6/vigun.Suvita2.gnm1.ann1.1PF6.cds_primary.fna.gz
  curl -O $url_base/unguiculata/annotations/TZ30.gnm1.ann2.59NL/vigun.TZ30.gnm1.ann2.59NL.cds_primary.fna.gz
  curl -O $url_base/unguiculata/annotations/UCR779.gnm1.ann1.VF6G/vigun.UCR779.gnm1.ann1.VF6G.cds_primary.fna.gz
  #curl -O $url_base/unguiculata/annotations/Xiabao_II.gnm1.ann1.4JFL/vigun.Xiabao_II.gnm1.ann1.4JFL.cds.fna.gz
  curl -O $url_base/unguiculata/annotations/ZN016.gnm1.ann2.C7YV/vigun.ZN016.gnm1.ann2.C7YV.cds_primary.fna.gz

  curl -O $url_base/angularis/annotations/Gyeongwon.gnm3.ann1.3Nz5/vigan.Gyeongwon.gnm3.ann1.3Nz5.protein_primary.faa.gz
  curl -O $url_base/angularis/annotations/Shumari.gnm1.ann1.8BRS/vigan.Shumari.gnm1.ann1.8BRS.protein_primary.faa.gz
  curl -O $url_base/radiata/annotations/VC1973A.gnm6.ann1.M1Qs/vigra.VC1973A.gnm6.ann1.M1Qs.protein.faa.gz
  curl -O $url_base/unguiculata/annotations/CB5-2.gnm1.ann1.0GKC/vigun.CB5-2.gnm1.ann1.0GKC.protein_primary.faa.gz
  curl -O $url_base/unguiculata/annotations/IT97K-499-35.gnm1.ann1.zb5D/vigun.IT97K-499-35.gnm1.ann1.zb5D.protein_primary.faa.gz
  curl -O $url_base/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.protein_primary.faa.gz
  curl -O $url_base/unguiculata/annotations/Sanzi.gnm1.ann1.HFH8/vigun.Sanzi.gnm1.ann1.HFH8.protein_primary.faa.gz
  curl -O $url_base/unguiculata/annotations/Suvita2.gnm1.ann1.1PF6/vigun.Suvita2.gnm1.ann1.1PF6.protein_primary.faa.gz
  curl -O $url_base/unguiculata/annotations/TZ30.gnm1.ann2.59NL/vigun.TZ30.gnm1.ann2.59NL.protein_primary.faa.gz
  curl -O $url_base/unguiculata/annotations/UCR779.gnm1.ann1.VF6G/vigun.UCR779.gnm1.ann1.VF6G.protein_primary.faa.gz
  #curl -O $url_base/unguiculata/annotations/Xiabao_II.gnm1.ann1.4JFL/vigun.Xiabao_II.gnm1.ann1.4JFL.protein.faa.gz
  curl -O $url_base/unguiculata/annotations/ZN016.gnm1.ann2.C7YV/vigun.ZN016.gnm1.ann2.C7YV.protein_primary.faa.gz

  curl -O $url_base/angularis/annotations/Gyeongwon.gnm3.ann1.3Nz5/vigan.Gyeongwon.gnm3.ann1.3Nz5.cds.bed.gz
  curl -O $url_base/angularis/annotations/Shumari.gnm1.ann1.8BRS/vigan.Shumari.gnm1.ann1.8BRS.cds.bed.gz
  curl -O $url_base/radiata/annotations/VC1973A.gnm6.ann1.M1Qs/vigra.VC1973A.gnm6.ann1.M1Qs.cds.bed.gz
  curl -O $url_base/unguiculata/annotations/CB5-2.gnm1.ann1.0GKC/vigun.CB5-2.gnm1.ann1.0GKC.cds.bed.gz
  curl -O $url_base/unguiculata/annotations/IT97K-499-35.gnm1.ann1.zb5D/vigun.IT97K-499-35.gnm1.ann1.zb5D.cds.bed.gz
  curl -O $url_base/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.cds.bed.gz
  curl -O $url_base/unguiculata/annotations/Sanzi.gnm1.ann1.HFH8/vigun.Sanzi.gnm1.ann1.HFH8.cds.bed.gz
  curl -O $url_base/unguiculata/annotations/Suvita2.gnm1.ann1.1PF6/vigun.Suvita2.gnm1.ann1.1PF6.cds.bed.gz
  curl -O $url_base/unguiculata/annotations/TZ30.gnm1.ann2.59NL/vigun.TZ30.gnm1.ann2.59NL.cds.bed.gz
  curl -O $url_base/unguiculata/annotations/UCR779.gnm1.ann1.VF6G/vigun.UCR779.gnm1.ann1.VF6G.cds.bed.gz
  #curl -O $url_base/unguiculata/annotations/Xiabao_II.gnm1.ann1.4JFL/vigun.Xiabao_II.gnm1.ann1.4JFL.cds.bed.gz
  curl -O $url_base/unguiculata/annotations/ZN016.gnm1.ann2.C7YV/vigun.ZN016.gnm1.ann2.C7YV.cds.bed.gz

# Fix the forms of chromosome IDs
  for file in *; do gunzip $file & done

  perl -pi -e 's/vigan.Gyeongwon.gnm3.Va(\d+)/vigan.Gyeongwon.gnm3.chr$1/' vigan.Gyeongwon.gnm3*.bed &
  perl -pi -e 's/vigan.Shumari.gnm1.Chr(\d+)/vigan.Shumari.gnm1.chr$1/' vigan.Shumari.gnm1*.bed &
  perl -pi -e 's/vigra.VC1973A.gnm6.Vr(\d+)/vigra.VC1973A.gnm6.chr$1/' vigra.VC1973A.gnm6*.bed &
  perl -pi -e 's/vigun.CB5-2.gnm1.chr(\d)\t/vigun.CB5-2.gnm1.chr0$1\t/' vigun.CB5-2.gnm1*.bed &
  perl -pi -e 's/vigun.IT97K-499-35.gnm1.Vu(\d+)/vigun.IT97K-499-35.gnm1.chr$1/' vigun.IT97K-499-35.gnm1.ann1*.bed &
  perl -pi -e 's/vigun.IT97K-499-35.gnm1.Vu(\d+)/vigun.IT97K-499-35.gnm1.chr$1/' vigun.IT97K-499-35.gnm1.ann2*.bed &
  perl -pi -e 's/vigun.Sanzi.gnm1.chr(\d)\t/vigun.Sanzi.gnm1.chr0$1\t/' vigun.Sanzi.gnm1.*.bed &
  perl -pi -e 's/vigun.Suvita2.gnm1.chr(\d)\t/vigun.Suvita2.gnm1.chr0$1\t/' vigun.Suvita2.gnm1*.bed &
  perl -pi -e 's/vigun.TZ30.gnm1.chr(\d)\t/vigun.TZ30.gnm1.chr0$1\t/' vigun.TZ30.gnm1*.bed &
  perl -pi -e 's/vigun.UCR779.gnm1.chr(\d)\t/vigun.UCR779.gnm1.chr0$1\t/' vigun.UCR779.gnm1*.bed &
  #perl -pi -e 's/vigun.Xiabao_II.gnm1.LG(\d)\t/vigun.Xiabao_II.gnm1.chr0$1\t/' vigun.Xiabao_II.gnm1*.bed &
  #perl -pi -e 's/vigun.Xiabao_II.gnm1.LG(\d\d)\t/vigun.Xiabao_II.gnm1.chr$1\t/' vigun.Xiabao_II.gnm1*.bed &
  perl -pi -e 's/vigun.ZN016.gnm1.chr(\d)\t/vigun.ZN016.gnm1.chr0$1\t/' vigun.ZN016.gnm1*.bed &

# Check chromosome names
  for file in *bed; do cat $file | awk '$1~/chr/ {print $1}' | sort | uniq -c ; echo; done

  for file in *bed; do echo $file; cat $file | awk '$1~/chr/ {print $1}' | sort | uniq | wc -l ; echo; done

# Re-compress the files
  for file in *bed *fna *faa; do gzip $file &
  done

cat << 'DATA' > expected_chr_matches.tsv
# Expected chromosome matches for Vigna unguiculata
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



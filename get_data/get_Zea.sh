# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp). 
# Retrieved files are typically nucleotide CDS files and/or protein files, and 
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data ]; then mkdir -p data; fi
if [ ! -d data_orig ]; then mkdir -p data_orig; fi
base_dir=$PWD

# Base URL for remote data repository, where annotation files are found
url_base="https://download.maizegdb.org"

cd data_orig
curl -O $url_base/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.cds.fa.gz
curl -O $url_base/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical.cds.fa.gz
curl -O $url_base/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0_Zm00018ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-CML103-REFERENCE-NAM-1.0/Zm-CML103-REFERENCE-NAM-1.0_Zm00021ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-CML228-REFERENCE-NAM-1.0/Zm-CML228-REFERENCE-NAM-1.0_Zm00022ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0_Zm00023ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0_Zm00024ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-CML322-REFERENCE-NAM-1.0/Zm-CML322-REFERENCE-NAM-1.0_Zm00025ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-CML333-REFERENCE-NAM-1.0/Zm-CML333-REFERENCE-NAM-1.0_Zm00026ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-CML52-REFERENCE-NAM-1.0/Zm-CML52-REFERENCE-NAM-1.0_Zm00019ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-CML69-REFERENCE-NAM-1.0/Zm-CML69-REFERENCE-NAM-1.0_Zm00020ab.1.cds.fa.gz
curl -O $url_base/Zm-HP301-REFERENCE-NAM-1.0/Zm-HP301-REFERENCE-NAM-1.0_Zm00027ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Il14H-REFERENCE-NAM-1.0/Zm-Il14H-REFERENCE-NAM-1.0_Zm00028ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0_Zm00030ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Ki3-REFERENCE-NAM-1.0/Zm-Ki3-REFERENCE-NAM-1.0_Zm00029ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-M162W-REFERENCE-NAM-1.0/Zm-M162W-REFERENCE-NAM-1.0_Zm00033ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-M37W-REFERENCE-NAM-1.0/Zm-M37W-REFERENCE-NAM-1.0_Zm00032ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Mo18W-REFERENCE-NAM-1.0/Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Ms71-REFERENCE-NAM-1.0/Zm-Ms71-REFERENCE-NAM-1.0_Zm00035ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0_Zm00036ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-NC358-REFERENCE-NAM-1.0/Zm-NC358-REFERENCE-NAM-1.0_Zm00037ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0_Zm00039ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0_Zm00038ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0_Zm00040ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0_Zm00008a.1.cds.fa.gz
curl -O $url_base/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0_Zm00041ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0_Zm00042ab.1.canonical.cds.fa.gz
curl -O $url_base/Zm-W22-REFERENCE-NRGENE-2.0/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.cds.fa.gz

curl -O $url_base/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.protein.fa.gz
curl -O $url_base/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa.gz
curl -O $url_base/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0_Zm00018ab.1.protein.fa.gz
curl -O $url_base/Zm-CML103-REFERENCE-NAM-1.0/Zm-CML103-REFERENCE-NAM-1.0_Zm00021ab.1.protein.fa.gz
curl -O $url_base/Zm-CML228-REFERENCE-NAM-1.0/Zm-CML228-REFERENCE-NAM-1.0_Zm00022ab.1.protein.fa.gz
curl -O $url_base/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0_Zm00023ab.1.protein.fa.gz
curl -O $url_base/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0_Zm00024ab.1.protein.fa.gz
curl -O $url_base/Zm-CML322-REFERENCE-NAM-1.0/Zm-CML322-REFERENCE-NAM-1.0_Zm00025ab.1.protein.fa.gz
curl -O $url_base/Zm-CML333-REFERENCE-NAM-1.0/Zm-CML333-REFERENCE-NAM-1.0_Zm00026ab.1.protein.fa.gz
curl -O $url_base/Zm-CML52-REFERENCE-NAM-1.0/Zm-CML52-REFERENCE-NAM-1.0_Zm00019ab.1.protein.fa.gz
curl -O $url_base/Zm-CML69-REFERENCE-NAM-1.0/Zm-CML69-REFERENCE-NAM-1.0_Zm00020ab.1.protein.fa.gz
curl -O $url_base/Zm-HP301-REFERENCE-NAM-1.0/Zm-HP301-REFERENCE-NAM-1.0_Zm00027ab.1.protein.fa.gz
curl -O $url_base/Zm-Il14H-REFERENCE-NAM-1.0/Zm-Il14H-REFERENCE-NAM-1.0_Zm00028ab.1.protein.fa.gz
curl -O $url_base/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0_Zm00030ab.1.protein.fa.gz
curl -O $url_base/Zm-Ki3-REFERENCE-NAM-1.0/Zm-Ki3-REFERENCE-NAM-1.0_Zm00029ab.1.protein.fa.gz
curl -O $url_base/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.protein.fa.gz
curl -O $url_base/Zm-M162W-REFERENCE-NAM-1.0/Zm-M162W-REFERENCE-NAM-1.0_Zm00033ab.1.protein.fa.gz
curl -O $url_base/Zm-M37W-REFERENCE-NAM-1.0/Zm-M37W-REFERENCE-NAM-1.0_Zm00032ab.1.protein.fa.gz
curl -O $url_base/Zm-Mo18W-REFERENCE-NAM-1.0/Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034ab.1.protein.fa.gz
curl -O $url_base/Zm-Ms71-REFERENCE-NAM-1.0/Zm-Ms71-REFERENCE-NAM-1.0_Zm00035ab.1.protein.fa.gz
curl -O $url_base/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0_Zm00036ab.1.protein.fa.gz
curl -O $url_base/Zm-NC358-REFERENCE-NAM-1.0/Zm-NC358-REFERENCE-NAM-1.0_Zm00037ab.1.protein.fa.gz
curl -O $url_base/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0_Zm00039ab.1.protein.fa.gz
curl -O $url_base/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0_Zm00038ab.1.protein.fa.gz
curl -O $url_base/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0_Zm00040ab.1.protein.fa.gz
curl -O $url_base/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0_Zm00008a.1.protein.fa.gz
curl -O $url_base/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0_Zm00041ab.1.protein.fa.gz
curl -O $url_base/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0_Zm00042ab.1.protein.fa.gz
curl -O $url_base/Zm-W22-REFERENCE-NRGENE-2.0/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.protein.fa.gz


curl -O $url_base/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz
curl -O $url_base/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz
curl -O $url_base/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0_Zm00018ab.1.gff3.gz
curl -O $url_base/Zm-CML103-REFERENCE-NAM-1.0/Zm-CML103-REFERENCE-NAM-1.0_Zm00021ab.1.gff3.gz
curl -O $url_base/Zm-CML228-REFERENCE-NAM-1.0/Zm-CML228-REFERENCE-NAM-1.0_Zm00022ab.1.gff3.gz
curl -O $url_base/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0_Zm00023ab.1.gff3.gz
curl -O $url_base/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0_Zm00024ab.1.gff3.gz
curl -O $url_base/Zm-CML322-REFERENCE-NAM-1.0/Zm-CML322-REFERENCE-NAM-1.0_Zm00025ab.1.gff3.gz
curl -O $url_base/Zm-CML333-REFERENCE-NAM-1.0/Zm-CML333-REFERENCE-NAM-1.0_Zm00026ab.1.gff3.gz
curl -O $url_base/Zm-CML52-REFERENCE-NAM-1.0/Zm-CML52-REFERENCE-NAM-1.0_Zm00019ab.1.gff3.gz
curl -O $url_base/Zm-CML69-REFERENCE-NAM-1.0/Zm-CML69-REFERENCE-NAM-1.0_Zm00020ab.1.gff3.gz
curl -O $url_base/Zm-HP301-REFERENCE-NAM-1.0/Zm-HP301-REFERENCE-NAM-1.0_Zm00027ab.1.gff3.gz
curl -O $url_base/Zm-Il14H-REFERENCE-NAM-1.0/Zm-Il14H-REFERENCE-NAM-1.0_Zm00028ab.1.gff3.gz
curl -O $url_base/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0_Zm00030ab.1.gff3.gz
curl -O $url_base/Zm-Ki3-REFERENCE-NAM-1.0/Zm-Ki3-REFERENCE-NAM-1.0_Zm00029ab.1.gff3.gz
curl -O $url_base/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.gff3.gz
curl -O $url_base/Zm-M162W-REFERENCE-NAM-1.0/Zm-M162W-REFERENCE-NAM-1.0_Zm00033ab.1.gff3.gz
curl -O $url_base/Zm-M37W-REFERENCE-NAM-1.0/Zm-M37W-REFERENCE-NAM-1.0_Zm00032ab.1.gff3.gz
curl -O $url_base/Zm-Mo18W-REFERENCE-NAM-1.0/Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034ab.1.gff3.gz
curl -O $url_base/Zm-Ms71-REFERENCE-NAM-1.0/Zm-Ms71-REFERENCE-NAM-1.0_Zm00035ab.1.gff3.gz
curl -O $url_base/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0_Zm00036ab.1.gff3.gz
curl -O $url_base/Zm-NC358-REFERENCE-NAM-1.0/Zm-NC358-REFERENCE-NAM-1.0_Zm00037ab.1.gff3.gz
curl -O $url_base/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0_Zm00039ab.1.gff3.gz
curl -O $url_base/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0_Zm00038ab.1.gff3.gz
curl -O $url_base/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0_Zm00040ab.1.gff3.gz
curl -O $url_base/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0_Zm00008a.1.gff3.gz
curl -O $url_base/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0_Zm00041ab.1.gff3.gz
curl -O $url_base/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0_Zm00042ab.1.gff3.gz
curl -O $url_base/Zm-W22-REFERENCE-NRGENE-2.0/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3.gz


echo "Pick longest CDS for several annotation sets, and exclude provisional"
echo "  For GRAMENE-4.0_Zm00001d, also remove 178 GRMZM5 gene models"
  zcat Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.cds.fa.gz | ../bin/fasta_to_table.awk |
    grep -v GRMZM5 |
    perl -pe 's/_T(\d+) ?/\tT$1\t/' |
    awk -v FS="\t" -v OFS="\t" '$0!~/provisional/ {print $1, length($4), $2, $3, $4}' |
    sort -k1,1 -k2nr,2nr | ../bin/top_line.awk | awk -v FS="\t" '{print ">" $1 "_" $3 " " $4; print $5}' |
    cat > ../data/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.canonical.cds.fa

  zcat Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.cds.fa.gz | ../bin/fasta_to_table.awk |
    perl -pe 's/_T(\d+) ?/\tT$1\t/' |
    awk -v FS="\t" -v OFS="\t" '$0!~/provisional/ {print $1, length($4), $2, $3, $4}' |
    sort -k1,1 -k2nr,2nr | ../bin/top_line.awk | awk -v FS="\t" '{print ">" $1 "_" $3 " " $4; print $5}' |
    cat > ../data/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.canonical.cds.fa

  zcat Zm-CML69-REFERENCE-NAM-1.0_Zm00020ab.1.cds.fa.gz | ../bin/fasta_to_table.awk |
    perl -pe 's/_T(\d+) ?/\tT$1\t/' |
    awk -v FS="\t" -v OFS="\t" '$0!~/provisional/ {print $1, length($4), $2, $3, $4}' |
    sort -k1,1 -k2nr,2nr | ../bin/top_line.awk | awk -v FS="\t" '{print ">" $1 "_" $3 " " $4; print $5}' |
    cat > ../data/Zm-CML69-REFERENCE-NAM-1.0_Zm00020ab.1.canonical.cds.fa

  zcat Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0_Zm00008a.1.cds.fa.gz | ../bin/fasta_to_table.awk |
    perl -pe 's/_T(\d+) ?/\tT$1\t/' |
    awk -v FS="\t" -v OFS="\t" '$0!~/provisional/ {print $1, length($4), $2, $3, $4}' |
    sort -k1,1 -k2nr,2nr | ../bin/top_line.awk | awk -v FS="\t" '{print ">" $1 "_" $3 " " $4; print $5}' |
    cat > ../data/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0_Zm00008a.1.canonical.cds.fa


echo "Pick longest proteins for several annotation sets, and exclude provisional"
echo "  For GRAMENE-4.0_Zm00001d, also remove 178 GRMZM5 gene models"
  zcat Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.protein.fa.gz | ../bin/fasta_to_table.awk |
    grep -v GRMZM5 |
    perl -pe 's/_P(\d+) ?/\tP$1\t/' |
    awk -v FS="\t" -v OFS="\t" '$0!~/provisional/ {print $1, length($4), $2, $3, $4}' |
    sort -k1,1 -k2nr,2nr | ../bin/top_line.awk | awk -v FS="\t" '{print ">" $1 "_" $3 " " $4; print $5}' |
    cat > ../data/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.canonical.protein.fa

  zcat Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.protein.fa.gz | ../bin/fasta_to_table.awk |
    perl -pe 's/_P(\d+) ?/\tP$1\t/' |
    awk -v FS="\t" -v OFS="\t" '$0!~/provisional/ {print $1, length($4), $2, $3, $4}' |
    sort -k1,1 -k2nr,2nr | ../bin/top_line.awk | awk -v FS="\t" '{print ">" $1 "_" $3 " " $4; print $5}' |
    cat > ../data/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.canonical.protein.fa

  zcat Zm-CML69-REFERENCE-NAM-1.0_Zm00020ab.1.protein.fa.gz | ../bin/fasta_to_table.awk |
    perl -pe 's/_P(\d+) ?/\tP$1\t/' |
    awk -v FS="\t" -v OFS="\t" '$0!~/provisional/ {print $1, length($4), $2, $3, $4}' |
    sort -k1,1 -k2nr,2nr | ../bin/top_line.awk | awk -v FS="\t" '{print ">" $1 "_" $3 " " $4; print $5}' |
    cat > ../data/Zm-CML69-REFERENCE-NAM-1.0_Zm00020ab.1.canonical.protein.fa

  zcat Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0_Zm00008a.1.protein.fa.gz | ../bin/fasta_to_table.awk |
    perl -pe 's/_P(\d+) ?/\tP$1\t/' |
    awk -v FS="\t" -v OFS="\t" '$0!~/provisional/ {print $1, length($4), $2, $3, $4}' |
    sort -k1,1 -k2nr,2nr | ../bin/top_line.awk | awk -v FS="\t" '{print ">" $1 "_" $3 " " $4; print $5}' |
    cat > ../data/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0_Zm00008a.1.canonical.protein.fa

  echo "Copy protein files from data_orig to data ..."
  echo "in the process, converting from e.g. Zm00004b003353_P001 to Zm00004b003353_T001,"
  echo "to allow ID matches between CDS and protein sequences."
  for file in *protein.fa.gz; do 
    base=`basename $file .gz`
    echo "  Copying $base, with ID modification"
    zcat $file | perl -pe 's/>(\S+)_P(\d+)/>$1_T$2/' > ../data/$base
  done

echo "Copy existing canonical files from data_orig/ to data/"
  for file in *canonical.*.fa.gz; do 
    gunzip $file &
  done
  wait
  cp *canonical.cds.fa ../data/

echo "Derive BED from GFF. "
echo "Add annotation name (e.g. Zm-W22_NRGENE-2) as prefix to the chromosome/scaffold names."
  for path in *gff3.gz; do
    base=`basename $path .gz`
    echo $base
    export annot_name=$(echo $base | perl -pe 's/(.+)-REFERENCE[-_](.+\d)\.\d_Z\w\d+.+/$1_$2/')
    zcat $path | ../bin/gff_to_bed6_mRNA.awk | 
      perl -pe '$prefix=$ENV{'annot_name'}; s/^(\D+)/$prefix.$1/; s/transcript://' |
       cat > ../data/$base.bed &
  done
  wait

echo "Change chromosome strings in Zm-B73-REFERENCE-GRAMENE-4.0"
  perl -pi -e 's/Chr/chr/' ../data/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.bed

echo "Re-compress files"
  for file in *.fa *.gff3; do
    gzip $file &
  done
  wait

  for file in ../data/*fa ../data/*bed; do
    gzip $file &
  done
  wait

echo "Check chromosome IDs"
  for file in ../data/*bed.gz; do echo $file; zcat $file | cut -f1 | awk '$1!~/ctg|scaf/' | sort | uniq -c; echo; done

cat <<DATA > ../data/expected_chr_matches.tsv
# Expected chromosome matches for Zea
1 1
2 2
3 3
4 4
5 5
6 6
7 7
8 8
9 9
9 10
10 10
DATA

cd $base_dir


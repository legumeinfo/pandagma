# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp). 
# Retrieved files are typically nucleotide CDS files and/or protein files, and 
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
if [ ! -d data ]; then mkdir -p data; fi
url_base="https://data.legumeinfo.org/Medicago"

## data
base_dir=$PWD
cd $base_dir/data
curl -f -O $url_base/sativa/annotations/XinJiangDaYe.gnm1.ann1.RKB9/medsa.XinJiangDaYe.gnm1.ann1.RKB9.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/A17.gnm5.ann1_6.L2RX/medtr.A17.gnm5.ann1_6.L2RX.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM004.gnm1.ann1.2XTB/medtr.HM004.gnm1.ann1.2XTB.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM010.gnm1.ann1.WV9J/medtr.HM010.gnm1.ann1.WV9J.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM022.gnm1.ann1.6C8N/medtr.HM022.gnm1.ann1.6C8N.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM023.gnm1.ann1.WZN8/medtr.HM023.gnm1.ann1.WZN8.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM034.gnm1.ann1.YR6S/medtr.HM034.gnm1.ann1.YR6S.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM050.gnm1.ann1.GWRX/medtr.HM050.gnm1.ann1.GWRX.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM056.gnm1.ann1.CHP6/medtr.HM056.gnm1.ann1.CHP6.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM058.gnm1.ann1.LXPZ/medtr.HM058.gnm1.ann1.LXPZ.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM060.gnm1.ann1.H41P/medtr.HM060.gnm1.ann1.H41P.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM095.gnm1.ann1.55W4/medtr.HM095.gnm1.ann1.55W4.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM125.gnm1.ann1.KY5W/medtr.HM125.gnm1.ann1.KY5W.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM129.gnm1.ann1.7FTD/medtr.HM129.gnm1.ann1.7FTD.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM185.gnm1.ann1.GB3D/medtr.HM185.gnm1.ann1.GB3D.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/HM324.gnm1.ann1.SQH2/medtr.HM324.gnm1.ann1.SQH2.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/R108_HM340.gnm1.ann1.85YW/medtr.R108_HM340.gnm1.ann1.85YW.gene_models_main.bed.gz
curl -f -O $url_base/truncatula/annotations/R108.gnmHiC_1.ann1.Y8NH/medtr.R108.gnmHiC_1.ann1.Y8NH.gene_models_main.bed.gz

curl -f -O $url_base/sativa/annotations/XinJiangDaYe.gnm1.ann1.RKB9/medsa.XinJiangDaYe.gnm1.ann1.RKB9.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.cds_primary.fna.gz
curl -f -O $url_base/truncatula/annotations/A17.gnm5.ann1_6.L2RX/medtr.A17.gnm5.ann1_6.L2RX.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM004.gnm1.ann1.2XTB/medtr.HM004.gnm1.ann1.2XTB.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM010.gnm1.ann1.WV9J/medtr.HM010.gnm1.ann1.WV9J.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM022.gnm1.ann1.6C8N/medtr.HM022.gnm1.ann1.6C8N.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM023.gnm1.ann1.WZN8/medtr.HM023.gnm1.ann1.WZN8.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM034.gnm1.ann1.YR6S/medtr.HM034.gnm1.ann1.YR6S.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM050.gnm1.ann1.GWRX/medtr.HM050.gnm1.ann1.GWRX.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM056.gnm1.ann1.CHP6/medtr.HM056.gnm1.ann1.CHP6.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM058.gnm1.ann1.LXPZ/medtr.HM058.gnm1.ann1.LXPZ.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM060.gnm1.ann1.H41P/medtr.HM060.gnm1.ann1.H41P.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM095.gnm1.ann1.55W4/medtr.HM095.gnm1.ann1.55W4.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM125.gnm1.ann1.KY5W/medtr.HM125.gnm1.ann1.KY5W.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM129.gnm1.ann1.7FTD/medtr.HM129.gnm1.ann1.7FTD.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM185.gnm1.ann1.GB3D/medtr.HM185.gnm1.ann1.GB3D.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/HM324.gnm1.ann1.SQH2/medtr.HM324.gnm1.ann1.SQH2.cds.fna.gz
curl -f -O $url_base/truncatula/annotations/R108_HM340.gnm1.ann1.85YW/medtr.R108_HM340.gnm1.ann1.85YW.cds_primary.fna.gz
curl -f -O $url_base/truncatula/annotations/R108.gnmHiC_1.ann1.Y8NH/medtr.R108.gnmHiC_1.ann1.Y8NH.cds.fna.gz

curl -f -O $url_base/sativa/annotations/XinJiangDaYe.gnm1.ann1.RKB9/medsa.XinJiangDaYe.gnm1.ann1.RKB9.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.protein_primary.faa.gz
curl -f -O $url_base/truncatula/annotations/A17.gnm5.ann1_6.L2RX/medtr.A17.gnm5.ann1_6.L2RX.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM004.gnm1.ann1.2XTB/medtr.HM004.gnm1.ann1.2XTB.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM010.gnm1.ann1.WV9J/medtr.HM010.gnm1.ann1.WV9J.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM022.gnm1.ann1.6C8N/medtr.HM022.gnm1.ann1.6C8N.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM023.gnm1.ann1.WZN8/medtr.HM023.gnm1.ann1.WZN8.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM034.gnm1.ann1.YR6S/medtr.HM034.gnm1.ann1.YR6S.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM050.gnm1.ann1.GWRX/medtr.HM050.gnm1.ann1.GWRX.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM056.gnm1.ann1.CHP6/medtr.HM056.gnm1.ann1.CHP6.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM058.gnm1.ann1.LXPZ/medtr.HM058.gnm1.ann1.LXPZ.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM060.gnm1.ann1.H41P/medtr.HM060.gnm1.ann1.H41P.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM095.gnm1.ann1.55W4/medtr.HM095.gnm1.ann1.55W4.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM125.gnm1.ann1.KY5W/medtr.HM125.gnm1.ann1.KY5W.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM129.gnm1.ann1.7FTD/medtr.HM129.gnm1.ann1.7FTD.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM185.gnm1.ann1.GB3D/medtr.HM185.gnm1.ann1.GB3D.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/HM324.gnm1.ann1.SQH2/medtr.HM324.gnm1.ann1.SQH2.protein.faa.gz
curl -f -O $url_base/truncatula/annotations/R108_HM340.gnm1.ann1.85YW/medtr.R108_HM340.gnm1.ann1.85YW.protein_primary.faa.gz
curl -f -O $url_base/truncatula/annotations/R108.gnmHiC_1.ann1.Y8NH/medtr.R108.gnmHiC_1.ann1.Y8NH.protein.faa.gz


echo "Change chromosome names for some of the files, to allow later filtering by chromosome correspondence"

echo "  Uncompress files"
for file in *; do gunzip $file & 
done

wait

echo "  Tweak chromosome names for medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.bed"
    perl -pi -e 's/MtrunA17Chr0c/scaff_/' medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.bed
    perl -pi -e 's/MtrunA17Chr/chr/; s/MtrunA17CP/CP/; s/MtrunA17MT/MT/' medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.bed

echo "  For sativa, split it into five files: scaffolds in one, and .1 .2 .3 .4 chromosomes into four others."
    cat medsa.XinJiangDaYe.gnm1.ann1.RKB9.gene_models_main.bed | awk '$1!~/chr/' | 
      perl -pe 's/(medsa.XinJiangDaYe.gnm1).(\d+)/$1.scaff$2/' > medsa.XinJiangDaYe_sc.gnm1.ann1.RKB9.gene_models_main.bed
    cat medsa.XinJiangDaYe.gnm1.ann1.RKB9.gene_models_main.bed | awk '$1~/chr[1-8].1/' > medsa.XinJiangDaYe_1.gnm1.ann1.RKB9.gene_models_main.bed
    cat medsa.XinJiangDaYe.gnm1.ann1.RKB9.gene_models_main.bed | awk '$1~/chr[1-8].2/' > medsa.XinJiangDaYe_2.gnm1.ann1.RKB9.gene_models_main.bed
    cat medsa.XinJiangDaYe.gnm1.ann1.RKB9.gene_models_main.bed | awk '$1~/chr[1-8].3/' > medsa.XinJiangDaYe_3.gnm1.ann1.RKB9.gene_models_main.bed
    cat medsa.XinJiangDaYe.gnm1.ann1.RKB9.gene_models_main.bed | awk '$1~/chr[1-8].4/' > medsa.XinJiangDaYe_4.gnm1.ann1.RKB9.gene_models_main.bed
    #rm medsa.XinJiangDaYe.gnm1.ann1.RKB9.gene_models_main.bed

    cat medsa.XinJiangDaYe.gnm1.ann1.RKB9.cds.fna | fasta_to_table.awk > medsa.XJDY.fna1
    cat medsa.XJDY.fna1 | awk '$2!~/chr/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_sc.gnm1.ann1.RKB9.cds.fna
    cat medsa.XJDY.fna1 | awk '$2~/chr[1-8].1/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_1.gnm1.ann1.RKB9.cds.fna
    cat medsa.XJDY.fna1 | awk '$2~/chr[1-8].2/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_2.gnm1.ann1.RKB9.cds.fna
    cat medsa.XJDY.fna1 | awk '$2~/chr[1-8].3/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_3.gnm1.ann1.RKB9.cds.fna
    cat medsa.XJDY.fna1 | awk '$2~/chr[1-8].4/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_4.gnm1.ann1.RKB9.cds.fna
    #rm medsa.XJDY.fna1
    #rm medsa.XinJiangDaYe.gnm1.ann1.RKB9.cds.fna

    cat medsa.XinJiangDaYe.gnm1.ann1.RKB9.protein.faa | perl -pe 's/ \[translate_table: standard\]//' | 
      fasta_to_table.awk > medsa.XJDY.faa1
    cat medsa.XJDY.faa1 | awk '$2!~/chr/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_sc.gnm1.ann1.RKB9.protein.faa
    cat medsa.XJDY.faa1 | awk '$2~/chr[1-8].1/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_1.gnm1.ann1.RKB9.protein.faa
    cat medsa.XJDY.faa1 | awk '$2~/chr[1-8].2/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_2.gnm1.ann1.RKB9.protein.faa
    cat medsa.XJDY.faa1 | awk '$2~/chr[1-8].3/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_3.gnm1.ann1.RKB9.protein.faa
    cat medsa.XJDY.faa1 | awk '$2~/chr[1-8].4/ {print ">" $1 " " $2 "\n" $3}' > medsa.XinJiangDaYe_4.gnm1.ann1.RKB9.protein.faa
    #rm medsa.XJDY.faa1
    #rm medsa.XinJiangDaYe.gnm1.ann1.RKB9.protein.faa

echo "  Then shorten the chr names from e.g. chr8.1 to chr8, because we want them co-equal (chr1.1 can match chr1.4)"
    perl -pi -e 's/(chr\d)\.\d/$1/' medsa.XinJiangDaYe_?.gnm1.ann1.RKB9.gene_models_main.bed

echo "  Name the sativa subgenome prefixes in both the chromosome and gene name fields, in all medsa.XinJiangDaYe files"
  perl -pi -e 's/XinJiangDaYe.gnm1/XinJiangDaYe_1.gnm1/g' medsa.XinJiangDaYe_1.gnm1.ann1.*
  perl -pi -e 's/XinJiangDaYe.gnm1/XinJiangDaYe_2.gnm1/g' medsa.XinJiangDaYe_2.gnm1.ann1.*
  perl -pi -e 's/XinJiangDaYe.gnm1/XinJiangDaYe_3.gnm1/g' medsa.XinJiangDaYe_3.gnm1.ann1.*
  perl -pi -e 's/XinJiangDaYe.gnm1/XinJiangDaYe_4.gnm1/g' medsa.XinJiangDaYe_4.gnm1.ann1.*
  perl -pi -e 's/XinJiangDaYe.gnm1/XinJiangDaYe_sc.gnm1/g' medsa.XinJiangDaYe_sc.gnm1.ann1.*

echo "  Fix a gene-name discrepancy in medtr.A17_HM341.gnm4.ann2"
    # bed file has gene names like this:    medtr.A17_HM341.gnm4.ann2.mRNA:MtrunA17Chr2g0284411
    # fasta file has gene names like this:  medtr.A17_HM341.gnm4.ann2.MtrunA17_Chr2g0284411
    # Use the latter form, by removing "mRNA:"
      perl -pi -e 's/mRNA://' medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.bed
      perl -pi -e 's/MtrunA17/MtrunA17_/' medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.bed
      perl -pi -e 's/chr0c/MtrunA17_Chr0c/' medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.bed

echo "  Fix funky scaffold names in  medtr.A17.gnm5.ann1_6.L2RX.gene_models_main.bed"
      perl -pi -e 's/MtrunA17Chr0c/scaff_/' medtr.A17.gnm5.ann1_6.L2RX.gene_models_main.bed


echo "  Check chromosome names. For some of these tests, the results should be empty."
      echo "Checking chromosome names:"; echo

echo "Check: Which annotations have chromosomes identifiers?"
      for file in *bed; do echo $file; cat $file | awk 'tolower($1) !~ /contig|scaf|sc/' | head -4; echo; done

echo "Check: Which annotations have something weird (not contig|scaf|sc|chr)?"
      for file in *bed; do echo $file; cat $file | awk 'tolower($1) !~ /contig|scaf|sc|chr/' | head -4; echo; done

echo "Check: how many chromosome-based features are in each annotation?"
      for file in *bed; do echo $file; cat $file | awk 'tolower($1) ~ /chr/ {print $1}' | sort | uniq -c ; echo; done

echo "Check: how many CHROMOSOMES are in each annotation?"
      for file in *bed; do echo $file; cat $file | awk 'tolower($1) ~ /chr/ {print $1}' | sort | uniq | wc -l; echo; done 
       # awk -v OFS="" '$1~/^med/ {print $1} $1~/[^123456789]/ {print "\t" $1} NF==0 {print "\n"}'


echo  "  Re-compress the files"
    for file in *bed *fna *.faa ; do gzip $file &
    done

wait

# Skip chromosome correspondences, because most of the assemblies/annotations have no chromosome assignments.
cat<<DATA >expected_chr_matches.tsv
# Expected chromosome matches for Medicago truncatula and Medicago sativa
1 1
2 2
3 3
4 4
5 5
6 6
7 7
8 8
4 8
DATA

cd $base_dir


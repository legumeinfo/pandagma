# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files (with gene positions), in GFF3 or BED format.
# Pandagma operates on seven-column BED files (BED6 + 7th column with gene ID), but these
# can be derived from GFF3 using the gff_to_bed7_mRNA.awk utility.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data_fam ]; then mkdir -p data_fam; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org"

## data
base_dir=$PWD
cd $base_dir/data_fam/

#curl -O $url_base/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds_primary.fna.gz
#curl -O $url_base/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz
#curl -O $url_base/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz
#
#curl -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.cds_primary.fna.gz
#curl -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.gff3.gz
#curl -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.protein_primary.faa.gz
#
#curl -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.cds_primary.fna.gz
#curl -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.gene_models_main.gff3.gz
#curl -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.protein_primary.faa.gz
#
#curl -O $url_base/Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.cds_primary.fna.gz
#curl -O $url_base/Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.gene_models_main.gff3.gz
#curl -O $url_base/Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.protein_primary.faa.gz
#
#curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.cds_primary.fna.gz
#curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.gene_models_main.gff3.gz
#curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.protein_primary.faa.gz
#
#curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.cds.fna.gz
#curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.gene_models_main.gff3.gz
#curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.protein.faa.gz
#
#curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.cds.fna.gz
#curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.gene_models_main.gff3.gz
#curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.protein.faa.gz
#
#
#curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.cds_primary.fna.gz
#curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.gene_models_main.gff3.gz
#curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.protein_primary.faa.gz
#
#curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.cds_primary.fna.gz
#curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.gene_models_main.gff3.gz
#curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.protein_primary.faa.gz
#
#curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.cds_primary.fna.gz
#curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.gene_models_main.gff3.gz
#curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.protein_primary.faa.gz
#
## Convert from gff3 to bed7 
#for file in *gff3.gz; do 
#  base=`basename $file .gff3.gz`
#  echo "Deriving a BED7 file from $base.gff3"
#  zcat $file | sort_gff.pl | gff_to_bed7_mRNA.awk > $base.bed
#  gzip $base.bed
#done

# Write file expected_quotas.tsv
cat <<DATA > expected_quotas.tsv
# Expected quotas
arath 2
glyma 4
medtr 2
phavu 2
pissa 2
prupe 1
sento 2
singl 2
vigun 2
vitvi 1
DATA

cd $base_dir


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

curl -O $url_base/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds_primary.fna.gz
curl -O $url_base/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz
curl -O $url_base/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz

curl -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.cds_primary.fna.gz
curl -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.gff3.gz
curl -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.protein_primary.faa.gz

curl -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.cds_primary.fna.gz
curl -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.gene_models_main.gff3.gz
curl -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.protein_primary.faa.gz

curl -O $url_base/Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.cds_primary.fna.gz
curl -O $url_base/Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.gene_models_main.gff3.gz
curl -O $url_base/Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.protein_primary.faa.gz

curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.cds_primary.fna.gz
curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.gene_models_main.gff3.gz
curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.protein_primary.faa.gz

curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.cds.fna.gz
curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.gene_models_main.gff3.gz
curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.protein.faa.gz

curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.cds.fna.gz
curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.gene_models_main.gff3.gz
curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.protein.faa.gz


curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.cds_primary.fna.gz
curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.gene_models_main.gff3.gz
curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.protein_primary.faa.gz

curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.cds_primary.fna.gz
curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.gene_models_main.gff3.gz
curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.protein_primary.faa.gz

curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.cds_primary.fna.gz
curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.gene_models_main.gff3.gz
curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.protein_primary.faa.gz

# Convert from gff3 to bed7 
for file in *gff3.gz; do 
  base=`basename $file .gff3.gz`
  echo "Deriving a BED7 file from $base.gff3"
  zcat $file | sort_gff.pl | gff_to_bed7_mRNA.awk > $base.bed
  gzip $base.bed
done

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

# Write file ks_peaks.tsv file.
# Note that this is typically derived from stats/ks_peaks_auto.tsv, 
# after user evaluation for events such as private WGDs in a lineage.
cat <<DATA > ks_peaks.tsv
# ks_peaks.tsv
glyma.Wm82	glyma.Wm82	0.65	4430
glyma.Wm82	medtr.A17_HM341	0.85	11573	# 0.60 auto
glyma.Wm82	phavu.G19833	0.70	5347
glyma.Wm82	pissa.Cameor	0.60	4194
glyma.Wm82	sento.Myeongyun	0.75	11567
glyma.Wm82	singl.CAF01	0.80	14427
glyma.Wm82	vigun.IT97K-499-35	0.75	5409
medtr.A17_HM341	medtr.A17_HM341	1.00	342
medtr.A17_HM341	phavu.G19833	0.95	5780	# 0.70 auto
medtr.A17_HM341	pissa.Cameor	1.00	1279	# 0.05 auto
medtr.A17_HM341	sento.Myeongyun	0.90	4583
medtr.A17_HM341	singl.CAF01	1.00	5156
medtr.A17_HM341	vigun.IT97K-499-35	1.00	8361	# 0.70 auto
phavu.G19833	phavu.G19833	0.85	1046
phavu.G19833	pissa.Cameor	0.70	2842
phavu.G19833	sento.Myeongyun	0.80	4929
phavu.G19833	singl.CAF01	0.85	6515
phavu.G19833	vigun.IT97K-499-35	0.80	2627
pissa.Cameor	pissa.Cameor	1.05	81
pissa.Cameor	sento.Myeongyun	1.00	1129
pissa.Cameor	singl.CAF01	1.00	1976
pissa.Cameor	vigun.IT97K-499-35	0.75	2908
sento.Myeongyun	sento.Myeongyun	0.60	922
sento.Myeongyun	singl.CAF01	0.70	7760
sento.Myeongyun	vigun.IT97K-499-35	0.85	6375
singl.CAF01	singl.CAF01	0.60	2983
singl.CAF01	vigun.IT97K-499-35	0.90	6820
vigun.IT97K-499-35	vigun.IT97K-499-35	0.90	1091
DATA

cd $base_dir


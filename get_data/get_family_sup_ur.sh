# This script retrieves data for species files to be placed into gene families based on 
# gene families calculated previously, by placement into best-scoring HMMSearch match.
# Retrieved files are nucleotide and protein files, filtered to the primary splice variants.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data_sup ]; then mkdir -p data_sup; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org"

## data
base_dir=$PWD
cd $base_dir/data_sup/

curl -f -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.protein_primary.faa.gz
curl -f -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.protein_primary.faa.gz
curl -f -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.protein_primary.faa.gz
curl -f -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.protein_primary.faa.gz
curl -f -O $url_base/annex/Bauhinia/variegata/annotations/BV-YZ2020.gnm2.ann1.RJ1G/bauva.BV-YZ2020.gnm2.ann1.RJ1G.protein_primary.faa.gz
curl -f -O $url_base/annex/Chamaecrista/fasciculata/annotations/ISC494698.gnm1.ann1.G7XW/chafa.ISC494698.gnm1.ann1.G7XW.protein_primary.faa.gz
curl -f -O $url_base/annex/Quillaja/saponaria/annotations/S10.gnm1.ann1.RQ4J/quisa.S10.gnm1.ann1.RQ4J.protein_primary.faa.gz
curl -f -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.protein.faa.gz
curl -f -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.protein.faa.gz
curl -f -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.protein_primary.faa.gz

curl -f -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.cds_primary.fna.gz
curl -f -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.cds_primary.fna.gz
curl -f -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.cds_primary.fna.gz
curl -f -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.cds_primary.fna.gz
curl -f -O $url_base/annex/Bauhinia/variegata/annotations/BV-YZ2020.gnm2.ann1.RJ1G/bauva.BV-YZ2020.gnm2.ann1.RJ1G.cds_primary.fna.gz
curl -f -O $url_base/annex/Chamaecrista/fasciculata/annotations/ISC494698.gnm1.ann1.G7XW/chafa.ISC494698.gnm1.ann1.G7XW.cds_primary.fna.gz
curl -f -O $url_base/annex/Quillaja/saponaria/annotations/S10.gnm1.ann1.RQ4J/quisa.S10.gnm1.ann1.RQ4J.cds_primary.fna.gz
curl -f -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.cds.fna.gz
curl -f -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.cds.fna.gz
curl -f -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.cds_primary.fna.gz

curl -f -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.gene_models_main.gff3.gz
curl -f -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.gff3.gz
curl -f -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.gene_models_main.gff3.gz
curl -f -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.gene_models_main.gff3.gz
curl -f -O $url_base/annex/Bauhinia/variegata/annotations/BV-YZ2020.gnm2.ann1.RJ1G/bauva.BV-YZ2020.gnm2.ann1.RJ1G.gene_models_main.gff3.gz
curl -f -O $url_base/annex/Chamaecrista/fasciculata/annotations/ISC494698.gnm1.ann1.G7XW/chafa.ISC494698.gnm1.ann1.G7XW.gene_models_main.gff3.gz
curl -f -O $url_base/annex/Quillaja/saponaria/annotations/S10.gnm1.ann1.RQ4J/quisa.S10.gnm1.ann1.RQ4J.gene_models_main.gff3.gz
curl -f -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.gene_models_main.gff3.gz
curl -f -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.gene_models_main.gff3.gz
curl -f -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.gene_models_main.gff3.gz

cd $base_dir


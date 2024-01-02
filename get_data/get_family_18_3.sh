# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files (with gene positions), in BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data ]; then mkdir -p data; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org"

## data
base_dir=$PWD
cd $base_dir/data/

curl -O $url_base/Aeschynomene/evenia/annotations/CIAT22838.gnm1.ann1.ZM3R/aesev.CIAT22838.gnm1.ann1.ZM3R.cds_primary.fna.gz
curl -O $url_base/Aeschynomene/evenia/annotations/CIAT22838.gnm1.ann1.ZM3R/aesev.CIAT22838.gnm1.ann1.ZM3R.gene_models_main.bed.gz
curl -O $url_base/Aeschynomene/evenia/annotations/CIAT22838.gnm1.ann1.ZM3R/aesev.CIAT22838.gnm1.ann1.ZM3R.protein_primary.faa.gz

curl -O $url_base/Arachis/GENUS/pangenes/Arachis.pan2.5P1K/Arachis.pan2.5P1K.pctl40_named_protein.faa.gz
curl -O $url_base/Arachis/GENUS/pangenes/Arachis.pan2.5P1K/Arachis.pan2.5P1K.pctl40_named_cds.fna.gz

curl -O $url_base/Cicer/GENUS/pangenes/Cicer.pan2.CMWZ/Cicer.pan2.CMWZ.pctl25_named_protein.faa.gz
curl -O $url_base/Cicer/GENUS/pangenes/Cicer.pan2.CMWZ/Cicer.pan2.CMWZ.pctl25_named_cds.fna.gz

curl -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.protein_primary.faa.gz
curl -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.cds_primary.fna.gz
curl -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.gene_models_main.bed.gz

curl -O $url_base/annex/Dalbergia/odorifera/annotations/SKLTGB.gnm1.ann1.R67B/dalod.SKLTGB.gnm1.ann1.R67B.cds.fna.gz
curl -O $url_base/annex/Dalbergia/odorifera/annotations/SKLTGB.gnm1.ann1.R67B/dalod.SKLTGB.gnm1.ann1.R67B.protein.faa.gz
curl -O $url_base/annex/Dalbergia/odorifera/annotations/SKLTGB.gnm1.ann1.R67B/dalod.SKLTGB.gnm1.ann1.R67B.gene_models_main.bed.gz

curl -O $url_base/Glycine/GENUS/pangenes/Glycine.pan4.RK4P/Glycine.pan4.RK4P.pctl25_named_protein.faa.gz
curl -O $url_base/Glycine/GENUS/pangenes/Glycine.pan4.RK4P/Glycine.pan4.RK4P.pctl25_named_cds.fna.gz

curl -O $url_base/Lens/culinaris/annotations/CDC_Redberry.gnm2.ann1.5FB4/lencu.CDC_Redberry.gnm2.ann1.5FB4.gene_models_main.bed.gz
curl -O $url_base/Lens/culinaris/annotations/CDC_Redberry.gnm2.ann1.5FB4/lencu.CDC_Redberry.gnm2.ann1.5FB4.cds.fna.gz
curl -O $url_base/Lens/culinaris/annotations/CDC_Redberry.gnm2.ann1.5FB4/lencu.CDC_Redberry.gnm2.ann1.5FB4.protein.faa.gz

curl -O $url_base/Lotus/japonicus/annotations/MG20.gnm3.ann1.WF9B/lotja.MG20.gnm3.ann1.WF9B.gene_models_main.bed.gz
curl -O $url_base/Lotus/japonicus/annotations/MG20.gnm3.ann1.WF9B/lotja.MG20.gnm3.ann1.WF9B.cds_primary.fna.gz
curl -O $url_base/Lotus/japonicus/annotations/MG20.gnm3.ann1.WF9B/lotja.MG20.gnm3.ann1.WF9B.protein_primary.faa.gz

curl -O $url_base/Lupinus/albus/annotations/Amiga.gnm1.ann1.3GKS/lupal.Amiga.gnm1.ann1.3GKS.gene_models_main.bed.gz
curl -O $url_base/Lupinus/albus/annotations/Amiga.gnm1.ann1.3GKS/lupal.Amiga.gnm1.ann1.3GKS.cds.fna.gz
curl -O $url_base/Lupinus/albus/annotations/Amiga.gnm1.ann1.3GKS/lupal.Amiga.gnm1.ann1.3GKS.protein.faa.gz

curl -O $url_base/Medicago/GENUS/pangenes/Medicago.pan2.451K/Medicago.pan2.451K.pctl25_named_protein.faa.gz
curl -O $url_base/Medicago/GENUS/pangenes/Medicago.pan2.451K/Medicago.pan2.451K.pctl25_named_cds.fna.gz

curl -O $url_base/Phaseolus/GENUS/pangenes/Phaseolus.pan2.G5HV/Phaseolus.pan2.G5HV.pctl25_named_protein.faa.gz
curl -O $url_base/Phaseolus/GENUS/pangenes/Phaseolus.pan2.G5HV/Phaseolus.pan2.G5HV.pctl25_named_cds.fna.gz

curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.gene_models_main.bed.gz
curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.cds_primary.fna.gz
curl -O $url_base/Pisum/sativum/annotations/Cameor.gnm1.ann1.7SZR/pissa.Cameor.gnm1.ann1.7SZR.protein_primary.faa.gz

curl -O $url_base/Vicia/faba/annotations/Hedin2.gnm1.ann1.PTNK/vicfa.Hedin2.gnm1.ann1.PTNK.gene_models_main.bed.gz
curl -O $url_base/Vicia/faba/annotations/Hedin2.gnm1.ann1.PTNK/vicfa.Hedin2.gnm1.ann1.PTNK.cds_primary.fna.gz
curl -O $url_base/Vicia/faba/annotations/Hedin2.gnm1.ann1.PTNK/vicfa.Hedin2.gnm1.ann1.PTNK.protein_primary.faa.gz

curl -O $url_base/Vigna/GENUS/pangenes/Vigna.pan2.MQQM/Vigna.pan2.MQQM.pctl25_named_protein.faa.gz
curl -O $url_base/Vigna/GENUS/pangenes/Vigna.pan2.MQQM/Vigna.pan2.MQQM.pctl25_named_cds.fna.gz

curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.protein_primary.faa.gz
curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.cds_primary.fna.gz
curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.gene_models_main.bed.gz

curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.protein_primary.faa.gz
curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.cds_primary.fna.gz
curl -O $url_base/annex/Prunus/persica/annotations/Lovell.gnm2.ann1.S2ZZ/prupe.Lovell.gnm2.ann1.S2ZZ.gene_models_main.bed.gz

curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.protein_primary.faa.gz
curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.cds_primary.fna.gz
curl -O $url_base/annex/Vitis/vinifera/annotations/PN40024.gnm2.ann1.V31M/vitvi.PN40024.gnm2.ann1.V31M.gene_models_main.bed.gz

curl -O $url_base/annex/Quillaja/saponaria/annotations/S10.gnm1.ann1.RQ4J/quisa.S10.gnm1.ann1.RQ4J.protein_primary.faa.gz
curl -O $url_base/annex/Quillaja/saponaria/annotations/S10.gnm1.ann1.RQ4J/quisa.S10.gnm1.ann1.RQ4J.cds_primary.fna.gz
curl -O $url_base/annex/Quillaja/saponaria/annotations/S10.gnm1.ann1.RQ4J/quisa.S10.gnm1.ann1.RQ4J.gene_models_main.bed.gz

curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.protein.faa.gz
curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.cds.fna.gz
curl -O $url_base/annex/Senna/tora/annotations/Myeongyun.gnm1.ann1.5WXB/sento.Myeongyun.gnm1.ann1.gene_models_main.bed.gz

curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.protein.faa.gz
curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.cds.fna.gz
curl -O $url_base/annex/Sindora/glabra/annotations/CAF01.gnm1.ann1.WFKC/singl.CAF01.gnm1.ann1.WFKC.gene_models_main.bed.gz

curl -O $url_base/annex/Chamaecrista/fasciculata/annotations/ISC494698.gnm1.ann1.G7XW/chafa.ISC494698.gnm1.ann1.G7XW.protein_primary.faa.gz
curl -O $url_base/annex/Chamaecrista/fasciculata/annotations/ISC494698.gnm1.ann1.G7XW/chafa.ISC494698.gnm1.ann1.G7XW.cds_primary.fna.gz
curl -O $url_base/annex/Chamaecrista/fasciculata/annotations/ISC494698.gnm1.ann1.G7XW/chafa.ISC494698.gnm1.ann1.G7XW.gene_models_main.bed.gz

# Rename CDS-based bed files to protein
for file in *.gene_models_main.bed.gz; do
  base=`basename $file .gene_models_main.bed.gz`
  mv $file $base.protein.bed.gz
done

# From the pan-gene protein files, derive bed files
for file in *.pctl25_named_protein.faa.gz; do 
  base=`basename $file .pctl25_named_protein.faa.gz`
  zcat $file | awk -v OFS="\t" '$1~/^>/ {print substr($1,2), $2, $3, substr($1,2), 0, $4}' | 
    perl -pe 's/^(\w+\.\w+\.[^_]+)_\S+/$1/' > $base.pctl25_named_protein.bed
done

# Also for Arachis, where the percentile value has been changed to 40
for file in *.pctl40_named_protein.faa.gz; do 
  base=`basename $file .pctl40_named_protein.faa.gz`
  zcat $file | awk -v OFS="\t" '$1~/^>/ {print substr($1,2), $2, $3, substr($1,2), 0, $4}' | 
    perl -pe 's/^(\w+\.\w+\.[^_]+)_\S+/$1/' > $base.pctl40_named_protein.bed
done

# Patch the lupal gene names - many of which look like this (with mRNA:) lupal.Amiga.gnm1.ann1.mRNA:Lalb_Chr25g0290201.1
  perl -pi -e 's/mRNA:/mRNA_/g' lupal*

# Check chromosome names
  for file in *bed.gz; do zcat $file | awk '$1~/chr|Chr|Ca/ {print $1}' | sort | uniq -c ; echo; done

  for file in *bed.gz; do echo $file; zcat $file | awk '$1~/chr|Chr|Ca/ {print $1}' | sort | uniq ; echo; done

# Re-compress the files
for file in *bed; do gzip $file &
done

wait


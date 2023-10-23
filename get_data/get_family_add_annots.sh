# This script retrieves data for species files to be placed into gene families based on 
# gene families calculated previously, by placement into best-scoring HMMSearch match.
# Retrieved files are nucleotide and protein files, filtered to the primary splice variants.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

if [ ! -d data_add_annots ]; then mkdir -p data_add_annots; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org"

## data
base_dir=$PWD
cd $base_dir/data_add_annots/

curl -O $url_base/Arachis/hypogaea/annotations/Tifrunner.gnm2.ann2.PVFB/arahy.Tifrunner.gnm2.ann2.PVFB.cds.fna.gz
curl -O $url_base/Cicer/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7/cicar.CDCFrontier.gnm3.ann1.NPD7.cds.fna.gz
curl -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.cds_primary.fna.gz
curl -O $url_base/Trifolium/pratense/annotations/MilvusB.gnm2.ann1.DFgp/tripr.MilvusB.gnm2.ann1.DFgp.cds_primary.fna.gz
curl -O $url_base/Lotus/japonicus/annotations/MG20.gnm3.ann1.WF9B/lotja.MG20.gnm3.ann1.WF9B.cds_primary.fna.gz
curl -O $url_base/Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.cds_primary.fna.gz
curl -O $url_base/Vigna/angularis/annotations/Gyeongwon.gnm3.ann1.3Nz5/vigan.Gyeongwon.gnm3.ann1.3Nz5.cds_primary.fna.gz
curl -O $url_base/Vigna/radiata/annotations/VC1973A.gnm7.ann1.RWBG/vigra.VC1973A.gnm7.ann1.RWBG.cds.fna.gz
curl -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.cds_primary.fna.gz
curl -O $url_base/Cajanus/cajan/annotations/ICPL87119.gnm2.ann1.L3ZH/cajca.ICPL87119.gnm2.ann1.L3ZH.cds_primary.fna.gz
curl -O $url_base/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds_primary.fna.gz
curl -O $url_base/Lupinus/albus/annotations/Amiga.gnm1.ann1.3GKS/lupal.Amiga.gnm1.ann1.3GKS.cds.fna.gz
curl -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.cds_primary.fna.gz
curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.cds_primary.fna.gz

curl -O $url_base/Arachis/hypogaea/annotations/Tifrunner.gnm2.ann2.PVFB/arahy.Tifrunner.gnm2.ann2.PVFB.protein.faa.gz
curl -O $url_base/Cicer/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7/cicar.CDCFrontier.gnm3.ann1.NPD7.protein.faa.gz
curl -O $url_base/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.protein_primary.faa.gz
curl -O $url_base/Trifolium/pratense/annotations/MilvusB.gnm2.ann1.DFgp/tripr.MilvusB.gnm2.ann1.DFgp.protein_primary.faa.gz
curl -O $url_base/Lotus/japonicus/annotations/MG20.gnm3.ann1.WF9B/lotja.MG20.gnm3.ann1.WF9B.protein_primary.faa.gz
curl -O $url_base/Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K/vigun.IT97K-499-35.gnm1.ann2.FD7K.protein_primary.faa.gz
curl -O $url_base/Vigna/angularis/annotations/Gyeongwon.gnm3.ann1.3Nz5/vigan.Gyeongwon.gnm3.ann1.3Nz5.protein_primary.faa.gz
curl -O $url_base/Vigna/radiata/annotations/VC1973A.gnm7.ann1.RWBG/vigra.VC1973A.gnm7.ann1.RWBG.protein.faa.gz
curl -O $url_base/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.protein_primary.faa.gz
curl -O $url_base/Cajanus/cajan/annotations/ICPL87119.gnm2.ann1.L3ZH/cajca.ICPL87119.gnm2.ann1.L3ZH.protein_primary.faa.gz
curl -O $url_base/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz
curl -O $url_base/Lupinus/albus/annotations/Amiga.gnm1.ann1.3GKS/lupal.Amiga.gnm1.ann1.3GKS.protein.faa.gz
curl -O $url_base/Cercis/canadensis/annotations/ISC453364.gnm3.ann1.3N1M/cerca.ISC453364.gnm3.ann1.3N1M.protein_primary.faa.gz
curl -O $url_base/annex/Arabidopsis/thaliana/annotations/Col0.gnm9.ann11.KH24/arath.Col0.gnm9.ann11.KH24.protein_primary.faa.gz

# Patch the lupal gene names - many of which look like this (with mRNA:) lupal.Amiga.gnm1.ann1.mRNA:Lalb_Chr25g0290201.1
  gunzip lupal.Amiga.gnm1.ann1.3GKS*

  perl -pi -e 's/mRNA:/mRNA_/g' lupal*

  for file in lupal.Amiga.gnm1.ann1.3GKS*; do gzip $file &
  done

cd $base_dir


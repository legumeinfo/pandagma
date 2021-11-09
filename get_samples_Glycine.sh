# This script retrieves data for a pan-gene calculation consisting of eight primary annotation sets
# and three extra ones. The primary ones are used for the initial clustering. The extra ones are
# added to the primary clusters by homology.

# This is a moderately large computational job. It is recommended to be run on a machine with 
# a significant number of available cores. 

# This version of pandagma has not been tested on a cluster with a job manager. Rather, jobs are run 
# in parallel using an indicated proportion of the available cores, as controlled by configuration 
# of the pandagma.sh script.

set -o errexit
set -o nounset

# Base URL for LegumeInfo/SoyBase Data Store, for genus Glycine
url_base="https://legumeinfo.org/data/v2/Glycine"

## data
base_dir=$PWD
cd $base_dir/data/
#  curl -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.cds_primaryTranscript.fna.gz
#  curl -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.cds.bed.gz
#
#  curl -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.cds_primaryTranscript.fna.gz
#  curl -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.cds.bed.gz
#
#  curl -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.cds_primaryTranscript.fna.gz
#  curl -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.cds.bed.gz
#
#  curl -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds_primaryTranscript.fna.gz
#  curl -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds.bed.gz
#
#  curl -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.cds_primaryTranscript.fna.gz
#  curl -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.cds.bed.gz
#
#  curl -O $url_base/max/annotations/Wm82.gnm2.ann2.BG1Q/glyma.Wm82.gnm2.ann2.BG1Q.cds_primaryTranscript.fna.gz
#  curl -O $url_base/max/annotations/Wm82.gnm2.ann2.BG1Q/glyma.Wm82.gnm2.ann2.BG1Q.cds.bed.gz
#
#  curl -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.cds_primaryTranscript.fna.gz
#  curl -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.cds.bed.gz
#
#  curl -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.cds_primaryTranscript.fna.gz
#  curl -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.cds.bed.gz

# "extra" sets
  curl -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.cds_primaryTranscript.fna.gz
  curl -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.cds.bed.gz

### Exclude glyma.Zh13.gnm2.ann1.FJ3G because gene models are too big by half: N50 2376 vs. expected 1572 (Wm82.gnm4.ann1)
##  curl -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.cds.bed.gz
##  curl -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.cds_primaryTranscript.fna.gz


  curl -O $url_base/max/annotations/Hefeng25_IGA1002.gnm1.ann1.320V/glyma.Hefeng25_IGA1002.gnm1.ann1.320V.cds.fna.gz
  curl -O $url_base/max/annotations/Hefeng25_IGA1002.gnm1.ann1.320V/glyma.Hefeng25_IGA1002.gnm1.ann1.320V.cds.bed.gz

  curl -O $url_base/max/annotations/Huaxia3_IGA1007.gnm1.ann1.LKC7/glyma.Huaxia3_IGA1007.gnm1.ann1.LKC7.cds.fna.gz
  curl -O $url_base/max/annotations/Huaxia3_IGA1007.gnm1.ann1.LKC7/glyma.Huaxia3_IGA1007.gnm1.ann1.LKC7.cds.bed.gz

  curl -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.cds.fna.gz
  curl -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.cds.bed.gz

  curl -O $url_base/max/annotations/Wenfeng7_IGA1001.gnm1.ann1.ZK5W/glyma.Wenfeng7_IGA1001.gnm1.ann1.ZK5W.cds.bed.gz
  curl -O $url_base/max/annotations/Wenfeng7_IGA1001.gnm1.ann1.ZK5W/glyma.Wenfeng7_IGA1001.gnm1.ann1.ZK5W.cds.fna.gz

  curl -O $url_base/max/annotations/Wm82_IGA1008.gnm1.ann1.FGN6/glyma.Wm82_IGA1008.gnm1.ann1.FGN6.cds.bed.gz
  curl -O $url_base/max/annotations/Wm82_IGA1008.gnm1.ann1.FGN6/glyma.Wm82_IGA1008.gnm1.ann1.FGN6.cds.fna.gz

  curl -O $url_base/max/annotations/Zh13_IGA1005.gnm1.ann1.87Z5/glyma.Zh13_IGA1005.gnm1.ann1.87Z5.cds.bed.gz
  curl -O $url_base/max/annotations/Zh13_IGA1005.gnm1.ann1.87Z5/glyma.Zh13_IGA1005.gnm1.ann1.87Z5.cds.fna.gz

  curl -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.cds.bed.gz
  curl -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.cds.fna.gz

  curl -O $url_base/soja/annotations/F_IGA1003.gnm1.ann1.G61B/glyso.F_IGA1003.gnm1.ann1.G61B.cds.bed.gz
  curl -O $url_base/soja/annotations/F_IGA1003.gnm1.ann1.G61B/glyso.F_IGA1003.gnm1.ann1.G61B.cds.fna.gz

# perennial species
#  curl -O $url_base/cyrtoloba/annotations/G1267.gnm1.ann1.HRFD/glycy.G1267.gnm1.ann1.HRFD.cds.fna.gz
#  curl -O $url_base/cyrtoloba/annotations/G1267.gnm1.ann1.HRFD/glycy.G1267.gnm1.ann1.HRFD.cds.bed.gz
#
#  curl -O $url_base/dolichocarpa/annotations/G1134.gnm1.ann1.4BJM/glydo.G1134.gnm1.ann1.4BJM.cds.fna.gz
#  curl -O $url_base/dolichocarpa/annotations/G1134.gnm1.ann1.4BJM/glydo.G1134.gnm1.ann1.4BJM.cds.bed.gz
#
#  curl -O $url_base/falcata/annotations/G1718.gnm1.ann1.2KSV/glyfa.G1718.gnm1.ann1.2KSV.cds.fna.gz
#  curl -O $url_base/falcata/annotations/G1718.gnm1.ann1.2KSV/glyfa.G1718.gnm1.ann1.2KSV.cds.bed.gz
#
#  curl -O $url_base/stenophita/annotations/G1974.gnm1.ann1.F257/glyst.G1974.gnm1.ann1.F257.cds.fna.gz
#  curl -O $url_base/stenophita/annotations/G1974.gnm1.ann1.F257/glyst.G1974.gnm1.ann1.F257.cds.bed.gz
#
#  curl -O $url_base/syndetika/annotations/G1300.gnm1.ann1.RRK6/glysy.G1300.gnm1.ann1.RRK6.cds.fna.gz
#  curl -O $url_base/syndetika/annotations/G1300.gnm1.ann1.RRK6/glysy.G1300.gnm1.ann1.RRK6.cds.bed.gz
#
#  curl -O $url_base/tomentella_D3/annotations/G1403.gnm1.ann1.XNZQ/glytoD3.G1403.gnm1.ann1.XNZQ.cds.fna.gz
#  curl -O $url_base/tomentella_D3/annotations/G1403.gnm1.ann1.XNZQ/glytoD3.G1403.gnm1.ann1.XNZQ.cds.bed.gz


cd $base_dir


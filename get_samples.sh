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
  curl -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.cds_primaryTranscript.fna.gz
  curl -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.cds.bed.gz

  curl -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.cds_primaryTranscript.fna.gz
  curl -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.cds.bed.gz

  curl -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.cds_primaryTranscript.fna.gz
  curl -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.cds.bed.gz

  curl -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds_primaryTranscript.fna.gz
  curl -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds.bed.gz

  curl -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.cds_primaryTranscript.fna.gz
  curl -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.cds.bed.gz

  curl -O $url_base/max/annotations/Wm82.gnm2.ann2.BG1Q/glyma.Wm82.gnm2.ann2.BG1Q.cds_primaryTranscript.fna.gz
  curl -O $url_base/max/annotations/Wm82.gnm2.ann2.BG1Q/glyma.Wm82.gnm2.ann2.BG1Q.cds.bed.gz

  curl -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.cds_primaryTranscript.fna.gz
  curl -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.cds.bed.gz

  curl -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.cds_primaryTranscript.fna.gz
  curl -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.cds.bed.gz

# "extra" sets
  curl -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.cds_primaryTranscript.fna.gz
  curl -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.cds.bed.gz

  curl -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.cds.bed.gz
  curl -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.cds_primaryTranscript.fna.gz

  curl -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.cds.fna.gz
  curl -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.cds.bed.gz

cd $base_dir


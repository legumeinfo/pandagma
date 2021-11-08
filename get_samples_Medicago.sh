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

# Base URL for LegumeInfo/SoyBase Data Store, for genus Medicago
url_base="https://legumeinfo.org/data/v2/Medicago"

## data
base_dir=$PWD
cd $base_dir/data/
  curl -O $url_base/sativa/annotations/XinJiangDaYe.gnm1.ann1.RKB9/medsa.XinJiangDaYe.gnm1.ann1.RKB9.cds.bed.gz
  curl -O $url_base/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM004.gnm1.ann1.2XTB/medtr.HM004.gnm1.ann1.2XTB.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM010.gnm1.ann1.WV9J/medtr.HM010.gnm1.ann1.WV9J.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM022.gnm1.ann1.6C8N/medtr.HM022.gnm1.ann1.6C8N.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM023.gnm1.ann1.WZN8/medtr.HM023.gnm1.ann1.WZN8.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM034.gnm1.ann1.YR6S/medtr.HM034.gnm1.ann1.YR6S.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM050.gnm1.ann1.GWRX/medtr.HM050.gnm1.ann1.GWRX.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM056.gnm1.ann1.CHP6/medtr.HM056.gnm1.ann1.CHP6.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM058.gnm1.ann1.LXPZ/medtr.HM058.gnm1.ann1.LXPZ.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM060.gnm1.ann1.H41P/medtr.HM060.gnm1.ann1.H41P.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM095.gnm1.ann1.55W4/medtr.HM095.gnm1.ann1.55W4.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM125.gnm1.ann1.KY5W/medtr.HM125.gnm1.ann1.KY5W.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM129.gnm1.ann1.7FTD/medtr.HM129.gnm1.ann1.7FTD.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM185.gnm1.ann1.GB3D/medtr.HM185.gnm1.ann1.GB3D.cds.bed.gz
  curl -O $url_base/truncatula/annotations/HM324.gnm1.ann1.SQH2/medtr.HM324.gnm1.ann1.SQH2.cds.bed.gz
  curl -O $url_base/truncatula/annotations/jemalong_A17.gnm5.ann1_6.L2RX/medtr.jemalong_A17.gnm5.ann1_6.L2RX.cds.bed.gz
  curl -O $url_base/truncatula/annotations/R108_HM340.gnm1.ann1.85YW/medtr.R108_HM340.gnm1.ann1.85YW.cds.bed.gz

  curl -O $url_base/sativa/annotations/XinJiangDaYe.gnm1.ann1.RKB9/medsa.XinJiangDaYe.gnm1.ann1.RKB9.cds.fna.gz
  curl -O $url_base/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.cds_primaryTranscript.fna.gz
  curl -O $url_base/truncatula/annotations/HM004.gnm1.ann1.2XTB/medtr.HM004.gnm1.ann1.2XTB.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM010.gnm1.ann1.WV9J/medtr.HM010.gnm1.ann1.WV9J.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM022.gnm1.ann1.6C8N/medtr.HM022.gnm1.ann1.6C8N.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM023.gnm1.ann1.WZN8/medtr.HM023.gnm1.ann1.WZN8.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM034.gnm1.ann1.YR6S/medtr.HM034.gnm1.ann1.YR6S.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM050.gnm1.ann1.GWRX/medtr.HM050.gnm1.ann1.GWRX.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM056.gnm1.ann1.CHP6/medtr.HM056.gnm1.ann1.CHP6.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM058.gnm1.ann1.LXPZ/medtr.HM058.gnm1.ann1.LXPZ.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM060.gnm1.ann1.H41P/medtr.HM060.gnm1.ann1.H41P.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM095.gnm1.ann1.55W4/medtr.HM095.gnm1.ann1.55W4.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM125.gnm1.ann1.KY5W/medtr.HM125.gnm1.ann1.KY5W.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM129.gnm1.ann1.7FTD/medtr.HM129.gnm1.ann1.7FTD.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM185.gnm1.ann1.GB3D/medtr.HM185.gnm1.ann1.GB3D.cds.fna.gz
  curl -O $url_base/truncatula/annotations/HM324.gnm1.ann1.SQH2/medtr.HM324.gnm1.ann1.SQH2.cds.fna.gz
  curl -O $url_base/truncatula/annotations/jemalong_A17.gnm5.ann1_6.L2RX/medtr.jemalong_A17.gnm5.ann1_6.L2RX.cds.fna.gz
  curl -O $url_base/truncatula/annotations/R108_HM340.gnm1.ann1.85YW/medtr.R108_HM340.gnm1.ann1.85YW.cds_primaryTranscript.fna.gz

# Shorten some filenames.
  ../scripts/rename.pl 's/primaryTranscript/primary/' *primaryTranscript*

cd $base_dir


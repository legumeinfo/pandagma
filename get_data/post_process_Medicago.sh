# For Medicago_7_16, the Medicago sativa accession medsa.XinJiangDaYe was split into five subgenomes 
# (chr1 - chr8 and sc for scaffolds), in order to use information from chromosome conservation relative 
# to M. truncatula chromosome assemblies. The post-processing in this file recovers the original gene IDs, e.g.
#   medsa.XinJiangDaYe_1.gnm1.ann1.MS.gene87716.t1 --> medsa.XinJiangDaYe.gnm1.ann1.MS.gene87716.t1
#   medsa.XinJiangDaYe_2.gnm1.ann1.MS.gene27472.t1 --> medsa.XinJiangDaYe.gnm1.ann1.MS.gene27472.t1

# Call this from the main pandagma.sh directory, outside out_Medicago_7_16:
#   ./get_data/post_process_Medicago.sh

perl -pi -e 's/(medsa.XinJiangDaYe)_[^.]+/$1/g' out_Medicago_7_16/*pan*



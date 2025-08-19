FILES ::= \
apiam.LA2127.gnm1.ann1.QKWH.gene_models_main.bed.gz \
apiam.LA2127.gnm1_hap2.ann1.NT2R.gene_models_main.bed.gz \
apipr.MO19963523.gnm1.ann1.P2VZ.gene_models_main.bed.gz \
apipr.MO19963523.gnm1_hap2.ann1.PB8C.gene_models_main.bed.gz \
apiam.LA2127.gnm1.ann1.QKWH.cds_primary.fna.gz \
apiam.LA2127.gnm1_hap2.ann1.NT2R.cds_primary.fna.gz \
apipr.MO19963523.gnm1.ann1.P2VZ.cds_primary.fna.gz \
apipr.MO19963523.gnm1_hap2.ann1.PB8C.cds_primary.fna.gz \
apiam.LA2127.gnm1.ann1.QKWH.protein_primary.faa.gz \
apiam.LA2127.gnm1_hap2.ann1.NT2R.protein_primary.faa.gz \
apipr.MO19963523.gnm1.ann1.P2VZ.protein_primary.faa.gz \
apipr.MO19963523.gnm1_hap2.ann1.PB8C.protein_primary.faa.gz \
legume.TE_lib_2024.rpt.6WVT.fna.gz

include $(dir $(realpath $(lastword $(MAKEFILE_LIST))))/common.mk

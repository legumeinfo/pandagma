FILES ::= \
arath.Col0.gnm9.ann11.KH24.cds_primary.fna.gz \
arath.Col0.gnm9.ann11.KH24.gene_models_main.bed.gz \
arath.Col0.gnm9.ann11.KH24.protein_primary.faa.gz \
glyma.Wm82.gnm4.ann1.T8TQ.cds_primary.fna.gz \
glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.bed.gz \
glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz \
medtr.A17_HM341.gnm4.ann2.G3ZY.cds_primary.fna.gz \
medtr.A17_HM341.gnm4.ann2.G3ZY.gene_models_main.bed.gz \
medtr.A17_HM341.gnm4.ann2.G3ZY.protein_primary.faa.gz \
phavu.G19833.gnm2.ann1.PB8d.cds_primary.fna.gz \
phavu.G19833.gnm2.ann1.PB8d.gene_models_main.bed.gz \
phavu.G19833.gnm2.ann1.PB8d.protein_primary.faa.gz \
pissa.Cameor.gnm1.ann1.7SZR.cds_primary.fna.gz \
pissa.Cameor.gnm1.ann1.7SZR.gene_models_main.bed.gz \
pissa.Cameor.gnm1.ann1.7SZR.protein_primary.faa.gz \
prupe.Lovell.gnm2.ann1.S2ZZ.cds_primary.fna.gz \
prupe.Lovell.gnm2.ann1.S2ZZ.gene_models_main.bed.gz \
prupe.Lovell.gnm2.ann1.S2ZZ.protein_primary.faa.gz \
sento.Myeongyun.gnm1.ann1.cds.fna.gz \
sento.Myeongyun.gnm1.ann1.gene_models_main.bed.gz \
sento.Myeongyun.gnm1.ann1.protein.faa.gz \
singl.CAF01.gnm1.ann1.WFKC.cds.fna.gz \
singl.CAF01.gnm1.ann1.WFKC.gene_models_main.bed.gz \
singl.CAF01.gnm1.ann1.WFKC.protein.faa.gz \
vigun.IT97K-499-35.gnm1.ann2.FD7K.cds_primary.fna.gz \
vigun.IT97K-499-35.gnm1.ann2.FD7K.gene_models_main.bed.gz \
vigun.IT97K-499-35.gnm1.ann2.FD7K.protein_primary.faa.gz \
vitvi.PN40024.gnm2.ann1.V31M.cds_primary.fna.gz \
vitvi.PN40024.gnm2.ann1.V31M.gene_models_main.bed.gz \
vitvi.PN40024.gnm2.ann1.V31M.protein_primary.faa.gz

all: $(FILES) ks_peaks.tsv

ks_peaks.tsv:
	cat <<DATA > ks_peaks.tsv
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


include $(dir $(realpath $(lastword $(MAKEFILE_LIST))))/common.mk

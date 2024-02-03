SHELL = bash
.SHELLFLAGS = -o errexit -o noglob -o nounset -o pipefail -o posix -c

.DELETE_ON_ERROR:
.ONESHELL:
.PHONY: all

all: $(FILES)

araduip.V14167K30076.gnm1.ann1.cxSMJ37m.cds.fna.gz \
araduip.V14167K30076.gnm1.ann1.cxSMJ37m.protein.faa.gz:
	suffix=$@
	suffix=$${suffix#araduip.V14167K30076.gnm1.ann1.cxSMJ37m.} # cds.fna.gz or protein.faa.gz
	{
	  fetch-datastore.sh aradu.V14167.gnm1.ann1.cxSM.$${suffix} | gzip -d
	  fetch-datastore.sh araip.K30076.gnm1.ann1.J37m.$${suffix} | gzip -d
	} | sed -e 's/^>aradu.V14167.gnm1/>araduip.V14167K30076.gnm1/' \
	        -e 's/^>araip.K30076.gnm1/>araduip.V14167K30076.gnm1/' | gzip > "$@"

araduip.V14167K30076.gnm1.ann1.cxSMJ37m.gene_models_main.bed.gz:
	{
	  fetch-datastore.sh aradu.V14167.gnm1.ann1.cxSM.gene_models_main.bed.gz | gzip -d
	  fetch-datastore.sh araip.K30076.gnm1.ann1.J37m.gene_models_main.bed.gz | gzip -d
	} | sed -e 's/aradu.V14167.gnm1.Adur/araduip.V14167K30076.gnm1.scaffA_/' \
	        -e 's/araip.K30076.gnm1.Aipa/araduip.V14167K30076.gnm1.scaffB_/' \
	        -e 's/aradu.V14167.gnm1.Aradu.A/araduip.V14167K30076.gnm1.chrA/' \
	        -e 's/^araip.K30076.gnm1.Araip.B0\([0-9]\)/araduip.V14167K30076.gnm1.chrB1\1/' \
	        -e 's/^araip.K30076.gnm1.Araip.B1\([0-9]\)/araduip.V14167K30076.gnm1.chrB2\1/' \
	        -e 's/aradu.V14167.gnm1/araduip.V14167K30076.gnm1/' \
	        -e 's/araip.K30076.gnm1/araduip.V14167K30076.gnm1/' | gzip > "$@"

arahy.Tifrunner.%.bed.gz:
	fetch-datastore.sh "$@" | gzip -d | sed 's/Arahy.\([0-9][0-9]\)/Arahy.Ah\1/' | gzip > "$@"

arahy.Tifrunner.gnm2.ann1.4K0L.cds_primary.fna.gz \
arahy.Tifrunner.gnm2.ann1.4K0L.protein_primary.faa.gz:
	fetch-datastore.sh arahy.Tifrunner.gnm2.ann1.4K0L.info_unmapped_models.txt.gz |
	  gzip -d |
	    awk 'NR == FNR && /^[^#]/ { sub(/gnm1/, "gnm2", $$1); exclude[$$1] }
	         NR != FNR && /^>/ { seqid=substr($$1,2); sub(/\.[0-9]+$$/, "", seqid); keep = !(seqid in exclude) }
	         NR != FNR && keep' - <(fetch-datastore.sh $@ | gzip -d) |
	      gzip > $@

# Patch the lupal gene names - many of which look like this (with mRNA:) lupal.Amiga.gnm1.ann1.mRNA:Lalb_Chr25g0290201.1
lupal.Amiga.gnm1.ann1.3GKS.%.gz:
	fetch-datastore.sh $@ | gzip -d | sed 's/mRNA:/mRNA_/' | gzip > $@



# Patch glyso.W05.gnm1.ann1 annotations: the CDS and protein files have Glysoja.10G027808, but the BED does not.
glyso.W05.gnm1.ann1.T47J.cds_primary.fna.gz \
glyso.W05.gnm1.ann1.T47J.protein_primary.faa.gz:
	fetch-datastore.sh $@ | gzip -d | awk '/^>/ {keep = !index($$1, "Glysoja.10G027808")} keep' | gzip > $@

# fix unusual scaffold naming
medtr.A17.gnm5.ann1_6.L2RX.gene_models_main.bed.gz:
	fetch-datastore.sh $@ | gzip -d | sed -e 's/MtrunA17Chr0c/scaff_/' | gzip > $@


# Standardize on *.chr\d\d chromosome names in Vigna bed files
vig%.bed.gz:
	fetch-datastore.sh "$@" |
	  gzip -d |
	    awk -v OFS='\t' '
	    $$1 ~ /\.[[:alpha:]]{2,3}[[:digit:]]{1,2}$$/ {
	      chr=ref=$$1;
	      sub(/.*+\.[^[:digit:]]+/, "", chr)
	      sub(/\.[^\.]+$$/, "", ref)
	      $$1 = sprintf("%s.chr%02i", ref, chr)
	    }
	    { print }' | gzip > "$@"


# create bed files by extracting info from defline of corresponding *_named_protein.faa files
%named_protein.bed.gz: %named_protein.faa.gz
	gzip -dc $< |
	  awk -v OFS='\t' '/^>/ { split($$1, seqid, /[>_]/); print seqid[2], $$2, $$3, substr($$1, 2), 0, $$4 }' |
	    gzip > $@

%.gz:
	fetch-datastore.sh "$@" > "$@"

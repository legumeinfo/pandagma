#!/usr/bin/env bash

set -o errexit -o noglob -o nounset -o pipefail -o posix
#set -o xtrace
readonly DATAFILE=${1}

: ${DATASTORE:=https://data.legumeinfo.org}

# adjust URL for collections that are located in the annex
case ${DATAFILE} in
  acacr.Acra3RX.gnm1.ann1.6C0V.*|\
  arahy.Tifrunner.gnm1.ann2.TN8K.*|\
  arath.Col0.gnm9.ann11.KH24.*|\
  bauva.BV-YZ2020.gnm2.ann1.RJ1G.*|\
  dalod.SKLTGB.gnm1.ann1.R67B.*|\
  phach.longxuteng.gnm1.ann1.KGX9.*|\
  prupe.Lovell.gnm2.ann1.S2ZZ.*|\
  quisa.S10.gnm1.ann1.RQ4J.*|\
  sento.Myeongyun.gnm1.ann1.5WXB.*|\
  singl.CAF01.gnm1.ann1.WFKC.*|\
  vitvi.PN40024.gnm2.ann1.V31M.*)
    DATASTORE=${DATASTORE}/annex
esac

genspa=${DATAFILE%%.*}
collection=${DATAFILE#*.}
trim_first=${DATAFILE#*.*}
TE_collection=${trim_first%.*.*}
collection=${collection%.*.*.*}
collection_type=annotations

case ${genspa} in
  [A-Z]*) genus=${genspa} species=GENUS collection_type=pangenes collection=${1%.*.*.*} ;;
  acacr) genus=Acacia species=crassicarpa ;;
  aesev) genus=Aeschynomene species=evenia ;;
  aradu) genus=Arachis species=duranensis ;;
  arahy) genus=Arachis species=hypogaea ;;
  araip) genus=Arachis species=ipaensis ;;
  arast) genus=Arachis species=stenosperma ;;
  arath) genus=Arabidopsis species=thaliana ;;
  bauva) genus=Bauhinia species=variegata ;;
  cajca) genus=Cajanus species=cajan ;;
  cerca) genus=Cercis species=canadensis ;;
  cerch) genus=Cercis species=chinensis ;;
  chafa) genus=Chamaecrista species=fasciculata ;;
  cicar) genus=Cicer species=arietinum ;;
  cicec) genus=Cicer species=echinospermum ;;
  cicre) genus=Cicer species=reticulatum ;;
  dalod) genus=Dalbergia species=odorifera ;;
  faial) genus=Faidherbia species=albida ;;
  glycy) genus=Glycine species=cyrtoloba ;;
  glyd3) genus=Glycine species=D3-tomentella ;;
  glydo) genus=Glycine species=dolichocarpa ;;
  glyfa) genus=Glycine species=falcata ;;
  glyma) genus=Glycine species=max ;;
  glyso) genus=Glycine species=soja ;;
  glyst) genus=Glycine species=stenophita ;;
  glysy) genus=Glycine species=syndetika ;;
  labpu) genus=Lablab species=purpureus ;;
  legume) genus=LEGUMES species=Fabaceae ;;
  lencu) genus=Lens species=culinaris ;;
  lener) genus=Lens species=ervoides ;;
  lotja) genus=Lotus species=japonicus ;;
  lupal) genus=Lupinus species=albus ;;
  lupan) genus=Lupinus species=angustifolius ;;
  medsa) genus=Medicago species=sativa ;;
  medtr) genus=Medicago species=truncatula ;;
  phaac) genus=Phaseolus species=acutifolius ;;
  phach) genus=Phanera species=championii ;;
  phaco) genus=Phaseolus species=coccineus ;;
  phalu) genus=Phaseolus species=lunatus ;;
  phavu) genus=Phaseolus species=vulgaris ;;
  pissa) genus=Pisum species=sativum ;;
  prupe) genus=Prunus species=persica ;;
  quisa) genus=Quillaja species=saponaria ;;
  sento) genus=Senna species=tora ;;
  singl) genus=Sindora species=glabra ;;
  tripr) genus=Trifolium species=pratense ;;
  trisu) genus=Trifolium species=subterraneum ;;
  vicfa) genus=Vicia species=faba ;;
  vigan) genus=Vigna species=angularis ;;
  vigra) genus=Vigna species=radiata ;;
  vigun) genus=Vigna species=unguiculata ;;
  vitvi) genus=Vitis species=vinifera ;;
  *) echo "ERROR: genus/species unknown: ${DATAFILE}" >&2; exit 1
esac

# Special case for file of transposable elements, indicated by string TE_lib
if [[ "$collection" == *"TE_lib"* ]]; then
  genus=LEGUMES
  species=Fabaceae
  collection_type=repeats
  collection=$TE_collection
fi

# Special case for Medicago sativa XinJiangDaYe.gnm1.ann1.RKB9, where the annotations have been split into four haplotypes
if [[ "$collection" == *"XinJiangDaYe"* ]]; then
  collection="XinJiangDaYe.gnm1.ann1.RKB9"
fi

curl --no-progress-meter --fail "${DATASTORE}/${genus}/${species}/${collection_type}/${collection}/${DATAFILE}"

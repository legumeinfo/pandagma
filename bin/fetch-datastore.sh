#!/usr/bin/env bash

set -o errexit -o noglob -o nounset -o pipefail -o posix
#set -o xtrace
readonly DATAFILE=${1}

: ${DATASTORE:=https://data.legumeinfo.org}

# adjust URL for collections that are located in the annex
case ${DATAFILE} in
  arahy.Tifrunner.gnm1.ann2.TN8K.*|\
  arath.Col0.gnm9.ann11.KH24.*|\
  bauva.BV-YZ2020.gnm2.ann1.RJ1G.*|\
  chafa.ISC494698.gnm1.ann1.G7XW.*|\
  dalod.SKLTGB.gnm1.ann1.R67B.*|\
  prupe.Lovell.gnm2.ann1.S2ZZ.*|\
  quisa.S10.gnm1.ann1.RQ4J.*|\
  sento.Myeongyun.gnm1.ann1.*|\
  singl.CAF01.gnm1.ann1.WFKC.*|\
  vitvi.PN40024.gnm2.ann1.V31M.*)
    DATASTORE=${DATASTORE}/annex 
esac

genspa=${DATAFILE%%.*}
collection=${DATAFILE#*.}
collection=${collection%.*.*.*}
collection_type=annotations

case ${genspa} in
  [A-Z]*) genus=${genspa} species=GENUS collection_type=pangenes collection=${1%.*.*.*} ;;
  aesev) genus=Aeschynomene species=evenia ;;
  aradu) genus=Arachis species=duranensis ;;
  arahy) genus=Arachis species=hypogaea ;;
  araip) genus=Arachis species=ipaensis ;;
  arast) genus=Arachis species=stenosperma ;;
  arath) genus=Arabidopsis species=thaliana ;;
  bauva) genus=Bauhinia species=variegata ;;
  cajca) genus=Cajanus species=cajan ;;
  cerca) genus=Cercis species=canadensis ;;
  chafa) genus=Chamaecrista species=fasciculata ;;
  cicar) genus=Cicer species=arietinum ;;
  cicec) genus=Cicer species=echinospermum ;;
  cicre) genus=Cicer species=reticulatum ;;
  dalod) genus=Dalbergia species=odorifera ;;
  glycy) genus=Glycine species=cyrtoloba ;;
  glyd3) genus=Glycine species=D3-tomentella ;;
  glydo) genus=Glycine species=dolichocarpa ;;
  glyfa) genus=Glycine species=falcata ;;
  glyma) genus=Glycine species=max ;;
  glyso) genus=Glycine species=soja ;;
  glyst) genus=Glycine species=stenophita ;;
  glysy) genus=Glycine species=syndetika ;;
  lencu) genus=Lens species=culinaris ;;
  lotja) genus=Lotus species=japonicus ;;
  lupal) genus=Lupinus species=albus ;;
  medsa) genus=Medicago species=sativa ;;
  medtr) genus=Medicago species=truncatula ;;
  phaac) genus=Phaseolus species=acutifolius ;;
  phalu) genus=Phaseolus species=lunatus ;;
  phavu) genus=Phaseolus species=vulgaris ;;
  pissa) genus=Pisum species=sativum ;;
  prupe) genus=Prunus species=persica ;;
  quisa) genus=Quillaja species=saponaria ;;
  sento) genus=Senna species=tora ;;
  singl) genus=Sindora species=glabra ;;
  tripr) genus=Trifolium species=pratense ;;
  vicfa) genus=Vicia species=faba ;;
  vigan) genus=Vigna species=angularis ;;
  vigra) genus=Vigna species=radiata ;;
  vigun) genus=Vigna species=unguiculata ;;
  vitvi) genus=Vitis species=vinifera ;;
  *) echo "ERROR: genus/species unknown: ${DATAFILE}" >&2; exit 1
esac

curl --no-progress-meter --fail "${DATASTORE}/${genus}/${species}/${collection_type}/${collection}/${DATAFILE}"

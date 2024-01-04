# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset
 
if [ ! -d data ]; then mkdir -p data; fi

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://data.legumeinfo.org/Glycine"

## data
base_dir=$PWD
cd $base_dir/data/

# CDS
curl -f -O $url_base/max/annotations/58-161.gnm1.ann1.HJ1K/glyma.58-161.gnm1.ann1.HJ1K.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Amsoy.gnm1.ann1.6S5P/glyma.Amsoy.gnm1.ann1.6S5P.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/DongNongNo_50.gnm1.ann1.QSDB/glyma.DongNongNo_50.gnm1.ann1.QSDB.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/FengDiHuang.gnm1.ann1.P6HL/glyma.FengDiHuang.gnm1.ann1.P6HL.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/HanDouNo_5.gnm1.ann1.ZS7M/glyma.HanDouNo_5.gnm1.ann1.ZS7M.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/HeiHeNo_43.gnm1.ann1.PDXG/glyma.HeiHeNo_43.gnm1.ann1.PDXG.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/JD17.gnm1.ann1.CLFP/glyma.JD17.gnm1.ann1.CLFP.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/JiDouNo_17.gnm1.ann1.X5PX/glyma.JiDouNo_17.gnm1.ann1.X5PX.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/JinDouNo_23.gnm1.ann1.SGJW/glyma.JinDouNo_23.gnm1.ann1.SGJW.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/JuXuanNo_23.gnm1.ann1.H8PW/glyma.JuXuanNo_23.gnm1.ann1.H8PW.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/KeShanNo_1.gnm1.ann1.2YX4/glyma.KeShanNo_1.gnm1.ann1.2YX4.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/PI_398296.gnm1.ann1.B0XR/glyma.PI_398296.gnm1.ann1.B0XR.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/PI_548362.gnm1.ann1.LL84/glyma.PI_548362.gnm1.ann1.LL84.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/QiHuangNo_34.gnm1.ann1.WHRV/glyma.QiHuangNo_34.gnm1.ann1.WHRV.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/ShiShengChangYe.gnm1.ann1.VLGS/glyma.ShiShengChangYe.gnm1.ann1.VLGS.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/TieFengNo_18.gnm1.ann1.7GR4/glyma.TieFengNo_18.gnm1.ann1.7GR4.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/TieJiaSiLiHuang.gnm1.ann1.W70Z/glyma.TieJiaSiLiHuang.gnm1.ann1.W70Z.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/TongShanTianEDan.gnm1.ann1.56XW/glyma.TongShanTianEDan.gnm1.ann1.56XW.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/WanDouNo_28.gnm1.ann1.NLYP/glyma.WanDouNo_28.gnm1.ann1.NLYP.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Wm82_ISU01.gnm2.ann1.FGFB/glyma.Wm82_ISU01.gnm2.ann1.FGFB.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/XuDouNo_1.gnm1.ann1.G2T7/glyma.XuDouNo_1.gnm1.ann1.G2T7.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/YuDouNo_22.gnm1.ann1.HCQ1/glyma.YuDouNo_22.gnm1.ann1.HCQ1.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/ZhangChunManCangJin.gnm1.ann1.7HPB/glyma.ZhangChunManCangJin.gnm1.ann1.7HPB.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Zhutwinning2.gnm1.ann1.ZTTQ/glyma.Zhutwinning2.gnm1.ann1.ZTTQ.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/ZiHuaNo_4.gnm1.ann1.FCFQ/glyma.ZiHuaNo_4.gnm1.ann1.FCFQ.cds_primary.fna.gz
curl -f -O $url_base/soja/annotations/PI_549046.gnm1.ann1.65KD/glyso.PI_549046.gnm1.ann1.65KD.cds_primary.fna.gz
curl -f -O $url_base/soja/annotations/PI_562565.gnm1.ann1.1JD2/glyso.PI_562565.gnm1.ann1.1JD2.cds_primary.fna.gz
curl -f -O $url_base/soja/annotations/PI_578357.gnm1.ann1.0ZKP/glyso.PI_578357.gnm1.ann1.0ZKP.cds_primary.fna.gz
curl -f -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.cds_primary.fna.gz

curl -f -O $url_base/cyrtoloba/annotations/G1267.gnm1.ann1.HRFD/glycy.G1267.gnm1.ann1.HRFD.cds.fna.gz
curl -f -O $url_base/D3-tomentella/annotations/G1403.gnm1.ann1.XNZQ/glyd3.G1403.gnm1.ann1.XNZQ.cds.fna.gz
curl -f -O $url_base/dolichocarpa/annotations/G1134.gnm1.ann1.4BJM/glydo.G1134.gnm1.ann1.4BJM.cds.fna.gz
curl -f -O $url_base/falcata/annotations/G1718.gnm1.ann1.2KSV/glyfa.G1718.gnm1.ann1.2KSV.cds.fna.gz
curl -f -O $url_base/max/annotations/Hefeng25_IGA1002.gnm1.ann1.320V/glyma.Hefeng25_IGA1002.gnm1.ann1.320V.cds.fna.gz
curl -f -O $url_base/max/annotations/Huaxia3_IGA1007.gnm1.ann1.LKC7/glyma.Huaxia3_IGA1007.gnm1.ann1.LKC7.cds.fna.gz
curl -f -O $url_base/max/annotations/Hwangkeum.gnm1.ann1.1G4F/glyma.Hwangkeum.gnm1.ann1.1G4F.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.cds.fna.gz
curl -f -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Lee.gnm2.ann1.1FNT/glyma.Lee.gnm2.ann1.1FNT.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Wenfeng7_IGA1001.gnm1.ann1.ZK5W/glyma.Wenfeng7_IGA1001.gnm1.ann1.ZK5W.cds.fna.gz
curl -f -O $url_base/max/annotations/Wm82_IGA1008.gnm1.ann1.FGN6/glyma.Wm82_IGA1008.gnm1.ann1.FGN6.cds.fna.gz
curl -f -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Zh13_IGA1005.gnm1.ann1.87Z5/glyma.Zh13_IGA1005.gnm1.ann1.87Z5.cds.fna.gz
curl -f -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.cds_primary.fna.gz
curl -f -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.cds.fna.gz
curl -f -O $url_base/soja/annotations/F_IGA1003.gnm1.ann1.G61B/glyso.F_IGA1003.gnm1.ann1.G61B.cds.fna.gz
curl -f -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.cds_primary.fna.gz
curl -f -O $url_base/stenophita/annotations/G1974.gnm1.ann1.F257/glyst.G1974.gnm1.ann1.F257.cds.fna.gz
curl -f -O $url_base/syndetika/annotations/G1300.gnm1.ann1.RRK6/glysy.G1300.gnm1.ann1.RRK6.cds.fna.gz
curl -f -O https://data.legumeinfo.org/annex/Glycine/max/annotations/Lee.gnm3.ann1.ZYY3/glyma.Lee.gnm3.ann1.ZYY3.cds.fna.gz
curl -f -O https://data.legumeinfo.org/annex/Glycine/max/annotations/Wm82.gnm5.ann1.J7HW/glyma.Wm82.gnm5.ann1.J7HW.cds.fna.gz

# BED
curl -f -O $url_base/max/annotations/58-161.gnm1.ann1.HJ1K/glyma.58-161.gnm1.ann1.HJ1K.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Amsoy.gnm1.ann1.6S5P/glyma.Amsoy.gnm1.ann1.6S5P.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/DongNongNo_50.gnm1.ann1.QSDB/glyma.DongNongNo_50.gnm1.ann1.QSDB.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/FengDiHuang.gnm1.ann1.P6HL/glyma.FengDiHuang.gnm1.ann1.P6HL.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/HanDouNo_5.gnm1.ann1.ZS7M/glyma.HanDouNo_5.gnm1.ann1.ZS7M.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/HeiHeNo_43.gnm1.ann1.PDXG/glyma.HeiHeNo_43.gnm1.ann1.PDXG.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/JD17.gnm1.ann1.CLFP/glyma.JD17.gnm1.ann1.CLFP.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/JiDouNo_17.gnm1.ann1.X5PX/glyma.JiDouNo_17.gnm1.ann1.X5PX.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/JinDouNo_23.gnm1.ann1.SGJW/glyma.JinDouNo_23.gnm1.ann1.SGJW.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/JuXuanNo_23.gnm1.ann1.H8PW/glyma.JuXuanNo_23.gnm1.ann1.H8PW.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/KeShanNo_1.gnm1.ann1.2YX4/glyma.KeShanNo_1.gnm1.ann1.2YX4.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/PI_398296.gnm1.ann1.B0XR/glyma.PI_398296.gnm1.ann1.B0XR.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/PI_548362.gnm1.ann1.LL84/glyma.PI_548362.gnm1.ann1.LL84.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/QiHuangNo_34.gnm1.ann1.WHRV/glyma.QiHuangNo_34.gnm1.ann1.WHRV.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/ShiShengChangYe.gnm1.ann1.VLGS/glyma.ShiShengChangYe.gnm1.ann1.VLGS.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/TieFengNo_18.gnm1.ann1.7GR4/glyma.TieFengNo_18.gnm1.ann1.7GR4.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/TieJiaSiLiHuang.gnm1.ann1.W70Z/glyma.TieJiaSiLiHuang.gnm1.ann1.W70Z.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/TongShanTianEDan.gnm1.ann1.56XW/glyma.TongShanTianEDan.gnm1.ann1.56XW.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/WanDouNo_28.gnm1.ann1.NLYP/glyma.WanDouNo_28.gnm1.ann1.NLYP.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Wm82_ISU01.gnm2.ann1.FGFB/glyma.Wm82_ISU01.gnm2.ann1.FGFB.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/XuDouNo_1.gnm1.ann1.G2T7/glyma.XuDouNo_1.gnm1.ann1.G2T7.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/YuDouNo_22.gnm1.ann1.HCQ1/glyma.YuDouNo_22.gnm1.ann1.HCQ1.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/ZhangChunManCangJin.gnm1.ann1.7HPB/glyma.ZhangChunManCangJin.gnm1.ann1.7HPB.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Zhutwinning2.gnm1.ann1.ZTTQ/glyma.Zhutwinning2.gnm1.ann1.ZTTQ.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/ZiHuaNo_4.gnm1.ann1.FCFQ/glyma.ZiHuaNo_4.gnm1.ann1.FCFQ.gene_models_main.bed.gz
curl -f -O $url_base/soja/annotations/PI_549046.gnm1.ann1.65KD/glyso.PI_549046.gnm1.ann1.65KD.gene_models_main.bed.gz
curl -f -O $url_base/soja/annotations/PI_562565.gnm1.ann1.1JD2/glyso.PI_562565.gnm1.ann1.1JD2.gene_models_main.bed.gz
curl -f -O $url_base/soja/annotations/PI_578357.gnm1.ann1.0ZKP/glyso.PI_578357.gnm1.ann1.0ZKP.gene_models_main.bed.gz
curl -f -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.gene_models_main.bed.gz

curl -f -O $url_base/cyrtoloba/annotations/G1267.gnm1.ann1.HRFD/glycy.G1267.gnm1.ann1.HRFD.gene_models_main.bed.gz
curl -f -O $url_base/D3-tomentella/annotations/G1403.gnm1.ann1.XNZQ/glyd3.G1403.gnm1.ann1.XNZQ.gene_models_main.bed.gz
curl -f -O $url_base/dolichocarpa/annotations/G1134.gnm1.ann1.4BJM/glydo.G1134.gnm1.ann1.4BJM.gene_models_main.bed.gz
curl -f -O $url_base/falcata/annotations/G1718.gnm1.ann1.2KSV/glyfa.G1718.gnm1.ann1.2KSV.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Hefeng25_IGA1002.gnm1.ann1.320V/glyma.Hefeng25_IGA1002.gnm1.ann1.320V.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Huaxia3_IGA1007.gnm1.ann1.LKC7/glyma.Huaxia3_IGA1007.gnm1.ann1.LKC7.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Hwangkeum.gnm1.ann1.1G4F/glyma.Hwangkeum.gnm1.ann1.1G4F.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Lee.gnm2.ann1.1FNT/glyma.Lee.gnm2.ann1.1FNT.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Wenfeng7_IGA1001.gnm1.ann1.ZK5W/glyma.Wenfeng7_IGA1001.gnm1.ann1.ZK5W.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Wm82_IGA1008.gnm1.ann1.FGN6/glyma.Wm82_IGA1008.gnm1.ann1.FGN6.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Zh13_IGA1005.gnm1.ann1.87Z5/glyma.Zh13_IGA1005.gnm1.ann1.87Z5.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.gene_models_main.bed.gz
curl -f -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.gene_models_main.bed.gz
curl -f -O $url_base/soja/annotations/F_IGA1003.gnm1.ann1.G61B/glyso.F_IGA1003.gnm1.ann1.G61B.gene_models_main.bed.gz
curl -f -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.gene_models_main.bed.gz
curl -f -O $url_base/stenophita/annotations/G1974.gnm1.ann1.F257/glyst.G1974.gnm1.ann1.F257.gene_models_main.bed.gz
curl -f -O $url_base/syndetika/annotations/G1300.gnm1.ann1.RRK6/glysy.G1300.gnm1.ann1.RRK6.gene_models_main.bed.gz
curl -f -O https://data.legumeinfo.org/annex/Glycine/max/annotations/Lee.gnm3.ann1.ZYY3/glyma.Lee.gnm3.ann1.ZYY3.gene_models_main.bed.gz
curl -f -O https://data.legumeinfo.org/annex/Glycine/max/annotations/Wm82.gnm5.ann1.J7HW/glyma.Wm82.gnm5.ann1.J7HW.gene_models_main.bed.gz

# Protein
curl -f -O $url_base/max/annotations/58-161.gnm1.ann1.HJ1K/glyma.58-161.gnm1.ann1.HJ1K.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Amsoy.gnm1.ann1.6S5P/glyma.Amsoy.gnm1.ann1.6S5P.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/DongNongNo_50.gnm1.ann1.QSDB/glyma.DongNongNo_50.gnm1.ann1.QSDB.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/FengDiHuang.gnm1.ann1.P6HL/glyma.FengDiHuang.gnm1.ann1.P6HL.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/HanDouNo_5.gnm1.ann1.ZS7M/glyma.HanDouNo_5.gnm1.ann1.ZS7M.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/HeiHeNo_43.gnm1.ann1.PDXG/glyma.HeiHeNo_43.gnm1.ann1.PDXG.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/JD17.gnm1.ann1.CLFP/glyma.JD17.gnm1.ann1.CLFP.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/JiDouNo_17.gnm1.ann1.X5PX/glyma.JiDouNo_17.gnm1.ann1.X5PX.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/JinDouNo_23.gnm1.ann1.SGJW/glyma.JinDouNo_23.gnm1.ann1.SGJW.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/JuXuanNo_23.gnm1.ann1.H8PW/glyma.JuXuanNo_23.gnm1.ann1.H8PW.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/KeShanNo_1.gnm1.ann1.2YX4/glyma.KeShanNo_1.gnm1.ann1.2YX4.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/PI_398296.gnm1.ann1.B0XR/glyma.PI_398296.gnm1.ann1.B0XR.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/PI_548362.gnm1.ann1.LL84/glyma.PI_548362.gnm1.ann1.LL84.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/QiHuangNo_34.gnm1.ann1.WHRV/glyma.QiHuangNo_34.gnm1.ann1.WHRV.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/ShiShengChangYe.gnm1.ann1.VLGS/glyma.ShiShengChangYe.gnm1.ann1.VLGS.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/TieFengNo_18.gnm1.ann1.7GR4/glyma.TieFengNo_18.gnm1.ann1.7GR4.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/TieJiaSiLiHuang.gnm1.ann1.W70Z/glyma.TieJiaSiLiHuang.gnm1.ann1.W70Z.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/TongShanTianEDan.gnm1.ann1.56XW/glyma.TongShanTianEDan.gnm1.ann1.56XW.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/WanDouNo_28.gnm1.ann1.NLYP/glyma.WanDouNo_28.gnm1.ann1.NLYP.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Wm82_ISU01.gnm2.ann1.FGFB/glyma.Wm82_ISU01.gnm2.ann1.FGFB.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/XuDouNo_1.gnm1.ann1.G2T7/glyma.XuDouNo_1.gnm1.ann1.G2T7.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/YuDouNo_22.gnm1.ann1.HCQ1/glyma.YuDouNo_22.gnm1.ann1.HCQ1.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/ZhangChunManCangJin.gnm1.ann1.7HPB/glyma.ZhangChunManCangJin.gnm1.ann1.7HPB.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Zhutwinning2.gnm1.ann1.ZTTQ/glyma.Zhutwinning2.gnm1.ann1.ZTTQ.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/ZiHuaNo_4.gnm1.ann1.FCFQ/glyma.ZiHuaNo_4.gnm1.ann1.FCFQ.protein_primary.faa.gz
curl -f -O $url_base/soja/annotations/PI_549046.gnm1.ann1.65KD/glyso.PI_549046.gnm1.ann1.65KD.protein_primary.faa.gz
curl -f -O $url_base/soja/annotations/PI_562565.gnm1.ann1.1JD2/glyso.PI_562565.gnm1.ann1.1JD2.protein_primary.faa.gz
curl -f -O $url_base/soja/annotations/PI_578357.gnm1.ann1.0ZKP/glyso.PI_578357.gnm1.ann1.0ZKP.protein_primary.faa.gz
curl -f -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.protein_primary.faa.gz

curl -f -O $url_base/cyrtoloba/annotations/G1267.gnm1.ann1.HRFD/glycy.G1267.gnm1.ann1.HRFD.protein.faa.gz
curl -f -O $url_base/D3-tomentella/annotations/G1403.gnm1.ann1.XNZQ/glyd3.G1403.gnm1.ann1.XNZQ.protein.faa.gz
curl -f -O $url_base/dolichocarpa/annotations/G1134.gnm1.ann1.4BJM/glydo.G1134.gnm1.ann1.4BJM.protein.faa.gz
curl -f -O $url_base/falcata/annotations/G1718.gnm1.ann1.2KSV/glyfa.G1718.gnm1.ann1.2KSV.protein.faa.gz
curl -f -O $url_base/max/annotations/Hefeng25_IGA1002.gnm1.ann1.320V/glyma.Hefeng25_IGA1002.gnm1.ann1.320V.protein.faa.gz
curl -f -O $url_base/max/annotations/Huaxia3_IGA1007.gnm1.ann1.LKC7/glyma.Huaxia3_IGA1007.gnm1.ann1.LKC7.protein.faa.gz
curl -f -O $url_base/max/annotations/Hwangkeum.gnm1.ann1.1G4F/glyma.Hwangkeum.gnm1.ann1.1G4F.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.protein.faa.gz
curl -f -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Lee.gnm2.ann1.1FNT/glyma.Lee.gnm2.ann1.1FNT.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Wenfeng7_IGA1001.gnm1.ann1.ZK5W/glyma.Wenfeng7_IGA1001.gnm1.ann1.ZK5W.protein.faa.gz
curl -f -O $url_base/max/annotations/Wm82_IGA1008.gnm1.ann1.FGN6/glyma.Wm82_IGA1008.gnm1.ann1.FGN6.protein.faa.gz
curl -f -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Zh13_IGA1005.gnm1.ann1.87Z5/glyma.Zh13_IGA1005.gnm1.ann1.87Z5.protein.faa.gz
curl -f -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.protein_primary.faa.gz
curl -f -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.protein.faa.gz
curl -f -O $url_base/soja/annotations/F_IGA1003.gnm1.ann1.G61B/glyso.F_IGA1003.gnm1.ann1.G61B.protein.faa.gz
curl -f -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.protein_primary.faa.gz
curl -f -O $url_base/stenophita/annotations/G1974.gnm1.ann1.F257/glyst.G1974.gnm1.ann1.F257.protein.faa.gz
curl -f -O $url_base/syndetika/annotations/G1300.gnm1.ann1.RRK6/glysy.G1300.gnm1.ann1.RRK6.protein.faa.gz
curl -f -O https://data.legumeinfo.org/annex/Glycine/max/annotations/Lee.gnm3.ann1.ZYY3/glyma.Lee.gnm3.ann1.ZYY3.protein.faa.gz
curl -f -O https://data.legumeinfo.org/annex/Glycine/max/annotations/Wm82.gnm5.ann1.J7HW/glyma.Wm82.gnm5.ann1.J7HW.protein.faa.gz

## The glyma.Zh13.gnm2.ann1.FJ3G gene models are too big by half: N50 2376 vs. expected 1572 (Wm82.gnm4.ann1)
## Exclude glyma.Wm82.gnm2.ann2.BG1Q because the conversion between GFF and CDS was flawed as of mid-2022. Maybe add later.

cat <<DATA > expected_chr_matches.tsv
# Expected chromosome matches for Glycine soja and Glycine max
01 01
02 02
03 03
04 04
05 05
06 06
07 07
08 08
09 09
10 10
11 13
11 11
12 12
13 13
14 14
15 15
16 16
17 17
18 18
19 19
20 20
DATA

cd $base_dir


# This script retrieves data for a pan-gene calculation, from a remote location (typically https or ftp).
# Retrieved files are typically nucleotide CDS files and/or protein files, and
# corresponding annotation files, in GFF3 or BED format.

# Edit this file to identify the base URL for the remote location, and the files to retrieve.

set -o errexit
set -o nounset

# Base URL for LegumeInfo/SoyBase Data Store, for this genus
url_base="https://legumeinfo.org/data/v2/Glycine"

## data
base_dir=$PWD
cd $base_dir/data/

curl -O $url_base/cyrtoloba/annotations/G1267.gnm1.ann1.HRFD/glycy.G1267.gnm1.ann1.HRFD.protein.faa.gz                        
curl -O $url_base/dolichocarpa/annotations/G1134.gnm1.ann1.4BJM/glydo.G1134.gnm1.ann1.4BJM.protein.faa.gz                     
curl -O $url_base/falcata/annotations/G1718.gnm1.ann1.2KSV/glyfa.G1718.gnm1.ann1.2KSV.protein.faa.gz                          
curl -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.protein_primary.faa.gz  
curl -O $url_base/max/annotations/Hefeng25_IGA1002.gnm1.ann1.320V/glyma.Hefeng25_IGA1002.gnm1.ann1.320V.protein.faa.gz        
curl -O $url_base/max/annotations/Huaxia3_IGA1007.gnm1.ann1.LKC7/glyma.Huaxia3_IGA1007.gnm1.ann1.LKC7.protein.faa.gz          
curl -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.protein.faa.gz          
curl -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.protein_primary.faa.gz                
curl -O $url_base/max/annotations/Wenfeng7_IGA1001.gnm1.ann1.ZK5W/glyma.Wenfeng7_IGA1001.gnm1.ann1.ZK5W.protein.faa.gz        
curl -O $url_base/max/annotations/Wm82_IGA1008.gnm1.ann1.FGN6/glyma.Wm82_IGA1008.gnm1.ann1.FGN6.protein.faa.gz                
curl -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.protein_primary.faa.gz              
curl -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.protein_primary.faa.gz              
curl -O $url_base/max/annotations/Wm82.gnm2.ann2.BG1Q/glyma.Wm82.gnm2.ann2.BG1Q.protein_primary.faa.gz              
curl -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz              
curl -O $url_base/max/annotations/Zh13_IGA1005.gnm1.ann1.87Z5/glyma.Zh13_IGA1005.gnm1.ann1.87Z5.protein.faa.gz                
curl -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.protein_primary.faa.gz              
## curl -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.protein_primary.faa.gz            
curl -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.protein.faa.gz                
curl -O $url_base/soja/annotations/F_IGA1003.gnm1.ann1.G61B/glyso.F_IGA1003.gnm1.ann1.G61B.protein.faa.gz                     
curl -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.protein_primary.faa.gz     
curl -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.protein_primary.faa.gz               
curl -O $url_base/stenophita/annotations/G1974.gnm1.ann1.F257/glyst.G1974.gnm1.ann1.F257.protein.faa.gz                       
curl -O $url_base/syndetika/annotations/G1300.gnm1.ann1.RRK6/glysy.G1300.gnm1.ann1.RRK6.protein.faa.gz                        
curl -O $url_base/tomentella_D3/annotations/G1403.gnm1.ann1.XNZQ/glytoD3.G1403.gnm1.ann1.XNZQ.protein.faa.gz                  
                                           
                                           
curl -O $url_base/cyrtoloba/annotations/G1267.gnm1.ann1.HRFD/glycy.G1267.gnm1.ann1.HRFD.cds.fna.gz
curl -O $url_base/dolichocarpa/annotations/G1134.gnm1.ann1.4BJM/glydo.G1134.gnm1.ann1.4BJM.cds.fna.gz
curl -O $url_base/falcata/annotations/G1718.gnm1.ann1.2KSV/glyfa.G1718.gnm1.ann1.2KSV.cds.fna.gz
curl -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.cds_primary.fna.gz
curl -O $url_base/max/annotations/Hefeng25_IGA1002.gnm1.ann1.320V/glyma.Hefeng25_IGA1002.gnm1.ann1.320V.cds.fna.gz
curl -O $url_base/max/annotations/Huaxia3_IGA1007.gnm1.ann1.LKC7/glyma.Huaxia3_IGA1007.gnm1.ann1.LKC7.cds.fna.gz
curl -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.cds.fna.gz
curl -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.cds_primary.fna.gz
curl -O $url_base/max/annotations/Wenfeng7_IGA1001.gnm1.ann1.ZK5W/glyma.Wenfeng7_IGA1001.gnm1.ann1.ZK5W.cds.fna.gz
curl -O $url_base/max/annotations/Wm82_IGA1008.gnm1.ann1.FGN6/glyma.Wm82_IGA1008.gnm1.ann1.FGN6.cds.fna.gz
curl -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.cds_primary.fna.gz
curl -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.cds_primary.fna.gz
curl -O $url_base/max/annotations/Wm82.gnm2.ann2.BG1Q/glyma.Wm82.gnm2.ann2.BG1Q.cds_primary.fna.gz
curl -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds_primary.fna.gz
curl -O $url_base/max/annotations/Zh13_IGA1005.gnm1.ann1.87Z5/glyma.Zh13_IGA1005.gnm1.ann1.87Z5.cds.fna.gz
curl -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.cds_primary.fna.gz
## curl -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.cds_primary.fna.gz
curl -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.cds.fna.gz
curl -O $url_base/soja/annotations/F_IGA1003.gnm1.ann1.G61B/glyso.F_IGA1003.gnm1.ann1.G61B.cds.fna.gz
curl -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.cds_primary.fna.gz
curl -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.cds_primary.fna.gz
curl -O $url_base/stenophita/annotations/G1974.gnm1.ann1.F257/glyst.G1974.gnm1.ann1.F257.cds.fna.gz
curl -O $url_base/syndetika/annotations/G1300.gnm1.ann1.RRK6/glysy.G1300.gnm1.ann1.RRK6.cds.fna.gz
curl -O $url_base/tomentella_D3/annotations/G1403.gnm1.ann1.XNZQ/glytoD3.G1403.gnm1.ann1.XNZQ.cds.fna.gz


curl -O $url_base/cyrtoloba/annotations/G1267.gnm1.ann1.HRFD/glycy.G1267.gnm1.ann1.HRFD.cds.bed.gz
curl -O $url_base/dolichocarpa/annotations/G1134.gnm1.ann1.4BJM/glydo.G1134.gnm1.ann1.4BJM.cds.bed.gz
curl -O $url_base/falcata/annotations/G1718.gnm1.ann1.2KSV/glyfa.G1718.gnm1.ann1.2KSV.cds.bed.gz
curl -O $url_base/max/annotations/FiskebyIII.gnm1.ann1.SS25/glyma.FiskebyIII.gnm1.ann1.SS25.cds.bed.gz
curl -O $url_base/max/annotations/Hefeng25_IGA1002.gnm1.ann1.320V/glyma.Hefeng25_IGA1002.gnm1.ann1.320V.cds.bed.gz
curl -O $url_base/max/annotations/Huaxia3_IGA1007.gnm1.ann1.LKC7/glyma.Huaxia3_IGA1007.gnm1.ann1.LKC7.cds.bed.gz
curl -O $url_base/max/annotations/Jinyuan_IGA1006.gnm1.ann1.2NNX/glyma.Jinyuan_IGA1006.gnm1.ann1.2NNX.cds.bed.gz
curl -O $url_base/max/annotations/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.cds.bed.gz
curl -O $url_base/max/annotations/Wenfeng7_IGA1001.gnm1.ann1.ZK5W/glyma.Wenfeng7_IGA1001.gnm1.ann1.ZK5W.cds.bed.gz
curl -O $url_base/max/annotations/Wm82_IGA1008.gnm1.ann1.FGN6/glyma.Wm82_IGA1008.gnm1.ann1.FGN6.cds.bed.gz
curl -O $url_base/max/annotations/Wm82.gnm1.ann1.DvBy/glyma.Wm82.gnm1.ann1.DvBy.cds.bed.gz
curl -O $url_base/max/annotations/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.cds.bed.gz
curl -O $url_base/max/annotations/Wm82.gnm2.ann2.BG1Q/glyma.Wm82.gnm2.ann2.BG1Q.cds.bed.gz
curl -O $url_base/max/annotations/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.cds.bed.gz
curl -O $url_base/max/annotations/Zh13_IGA1005.gnm1.ann1.87Z5/glyma.Zh13_IGA1005.gnm1.ann1.87Z5.cds.bed.gz
curl -O $url_base/max/annotations/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.cds.bed.gz
## curl -O $url_base/max/annotations/Zh13.gnm2.ann1.FJ3G/glyma.Zh13.gnm2.ann1.FJ3G.cds.bed.gz
curl -O $url_base/max/annotations/Zh35_IGA1004.gnm1.ann1.RGN6/glyma.Zh35_IGA1004.gnm1.ann1.RGN6.cds.bed.gz
curl -O $url_base/soja/annotations/F_IGA1003.gnm1.ann1.G61B/glyso.F_IGA1003.gnm1.ann1.G61B.cds.bed.gz
curl -O $url_base/soja/annotations/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.cds.bed.gz
curl -O $url_base/soja/annotations/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.cds.bed.gz
curl -O $url_base/stenophita/annotations/G1974.gnm1.ann1.F257/glyst.G1974.gnm1.ann1.F257.cds.bed.gz
curl -O $url_base/syndetika/annotations/G1300.gnm1.ann1.RRK6/glysy.G1300.gnm1.ann1.RRK6.cds.bed.gz
curl -O $url_base/tomentella_D3/annotations/G1403.gnm1.ann1.XNZQ/glytoD3.G1403.gnm1.ann1.XNZQ.cds.bed.gz

## Exclude glyma.Zh13.gnm2.ann1.FJ3G because gene models are too big by half: N50 2376 vs. expected 1572 (Wm82.gnm4.ann1)

cat << DATA > expected_chr_matches.tsv
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


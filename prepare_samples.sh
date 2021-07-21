# Prepare example data for a run of pandagma. The example data are from the 
# LegumeInfo/SoyBase Data Store, for genus Glycine. Run this script after running get_samples.sh .
# Call both scripts from where these scripts sit, along with the data and data_extra directories.

echo "Uncompress files"
for path in data/*gz; do gunzip $path &
done
wait

for path in data_extra/*gz; do gunzip $path &
done
wait

rename 's{(\w+\.\w+\.gnm\d+\.ann\d+)\..+\.fna}{$1.fna}' data*/*fna
rename 's{(\w+\.\w+\.gnm\d+\.ann\d+)\..+\.gff3}{$1.gff3}' data*/*gff3

echo
echo "Files should now be ready for analysis with pandagma."
echo "To do next:"
echo "(1) Make sure that the required tools are installed and on your path (see README):"
echo "      perl-bioperl-core, dagchainer, mcl, vsearch, mmseqs2"
echo "(2) Make sure that you are in a cluster environment with sufficient cores. Recommended: >=24 cores."
echo "(3) Execute the script that calls pandagma.sh:"
echo "    nohup ./run_pandagma.sh &"
echo

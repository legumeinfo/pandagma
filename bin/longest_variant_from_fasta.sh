#!/usr/bin/env sh

if test -t 0; then
cat <<Usage-message
  Usage: cat FILE.fa |  longest_variant_from_fasta.pl

  For fasta files (e.g. CDS, protein, mRNA) with splice variants of the form
    >GENE_ID.1      or >GENE_ID.m1 
  (variant signified by dot-separated digit at the end, optionally with leading alpha characters), 
  this script will return only the longest variant of each gene.
  If there are multiple variants with the same length, it will return the one with the lowest variant number, 
  e.g. .1 < .2 < .10  and .m1 < .m2 < .m10

Usage-message
exit
fi

##################################################
# Functions (just to help document what's happening)

# Print fasta onto one line and report sequence length.
# NOTE: This discards the description if one is present.
fasta_to_table() {
  awk 'BEGIN { ORS="" }
       /^>/ && NR==1 { print substr($1,2) "\t" }
       /^>/ && NR!=1 { print "\n" substr($1,2) "\t" }
       /^[^>]/ { print }
       END { print "\n" }'
}

# Split gene splice_variant. Handle variants such as "m1", "mRNA1", and "1"
# Gene forms such as mikado.chr01G1.1, mRNA38015, IDmodified-mrna-3934
# For a gene of length 999 ...
#   ... output for "GENEID.10"  will be  "GENEID zPLACEHOLDERz 10  999"
#   ... output for "GENEID.m10" will be  "GENEID       m       10  999"
# For gene mRNA8903, the output is: "mRNA8903	zPLACEHOLDERz	zQz"
split_gene_splicevar() {
  perl -ne 'if(/^(\S+)\.(\d+)\t(\S+)/){print "$1\tzPLACEHOLDERz\t$2\t$3\n"};
            if(/^([^.]+)\t(\S+)/){print "$1\tzPLACEHOLDERz\tzQz\t$2\n"};
            if(/^(\S+)\.([^\d+])(\d+)\t(\S+)/){print "$1\t$2\t$3\t$4\n"}' |
  awk -v OFS="\t" '{print $1, $2, $3, length($4), $4}'
}

# Sort by geneID, then length (reverse numerically), then by the splice number.
# The last item (splice number) will break ties among variants of the same length,
# putting 1 before 10 and 8 before 9.
sort_by_length() {
  sort -k1,1 -k4nr,4nr -k3n,3n
}

# Top line
top_line() {
  awk 'BEGIN { MAX = 1 }
       $1 == prev && count < MAX { print; count++ }
       $1 != prev { print; count = 1; prev = $1 }'
}

# Reassemble fasta
table_to_fasta() {
  awk -v OFS="" '{print ">" $1 "." $2 $3 "\n" $5}' |
  perl -pe 's/zPLACEHOLDERz//; s/\.zQz//' |
  fold -w100
}

##################################################
# Handle fasta sequence on stdin
fasta_to_table |
split_gene_splicevar |
sort_by_length |
top_line |
table_to_fasta


##########
# Versions
# 2022-12-14 Add sorting of the splice number (e.g. .1 or .m1) to 
#             break ties among variants of the same length.
# 2022-12-15 Report usage message if script is called without stdin
# 2023-03-15 Rename from pick_longest_variant.sh to longest_variant_from_fasta.pl


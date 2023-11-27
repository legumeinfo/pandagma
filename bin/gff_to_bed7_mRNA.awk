#!/usr/bin/awk -f
#
# NAME
#   gff_to_bed7_mRNA.awk - Generate mRNA-level seven-column BED file from GFF
#     Columns: molecule, feature-start, feature-end, mRNA-ID, score(0), strand, gene-ID
#  
#   NOTE: The seventh column (gene-ID) is nonstandard for BED. Retrieving this field
#   allows the BED file to be used to determine the mapping between mRNA and gene IDs.
# 
# SYNOPSIS
#   cat FILE.[gff3] | gff_to_bed7_mRNA.awk 
#     
# OPERANDS
#     INPUT_FILE
#         A file containing GFF data

BEGIN { FS = OFS = "\t" }
NF == 9 && $3 == "mRNA" {
  match($9, /ID=[^;]+/)
  mRNA_ID=substr($9, RSTART+3, RLENGTH-3)

  match($9, /Parent=[^;]+/)
  Parent=substr($9, RSTART+7, RLENGTH-7)

  transcript_Parent[mRNA_ID] = Parent
}
NF == 9 && $3 == "CDS" { 
  match($9, /Parent=[^;]+/)
  mRNA_ID=substr($9, RSTART+7, RLENGTH-7)
  nf = split(transcript_CDS[mRNA_ID], A, ",")
  if (nf == 0)
    transcript_CDS[mRNA_ID] = $1 "," $4 "," $5 "," $7
  else {
    if ($4 < A[2]) A[2] = $4 # min start pos
    if ($5 > A[3]) A[3] = $5 # max end pos
    transcript_CDS[mRNA_ID] = A[1] "," A[2] "," A[3] "," A[4]
  }
}
END {
  for (mRNA_ID in transcript_CDS) { 
    split(transcript_CDS[mRNA_ID], A, ",")
    print A[1], A[2]-1, A[3], mRNA_ID, 0, A[4], transcript_Parent[mRNA_ID]
  }
}

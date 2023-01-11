#!/usr/bin/awk -f
#
# NAME
#   gff_to_bed6_mRNA.awk - Generate mRNA-level six-column BED file from GFF
#     Columns: molecule, feature-start, feature-end, feature-ID, score(0), strand
# 
# SYNOPSIS
#   cat FILE.[gff3] | gff_to_bed6_mRNA.awk 
#     
# OPERANDS
#     INPUT_FILE
#         A file containing GFF data
BEGIN { FS = OFS = "\t" }
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
        print A[1], A[2]-1, A[3], mRNA_ID, 0, A[4]
    }
}

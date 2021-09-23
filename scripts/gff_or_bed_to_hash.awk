#!/usr/bin/awk -f
#
# NAME
#   gff_or_bed_to_hash.awk - Generate hash of gene ID and positional information
# 
# SYNOPSIS
#   cat FILE.[gff3|bed] | gff_or_bed_to_hash.awk 
#     
# OPERANDS
#     INPUT_FILE
#         A file containing GFF or four-column BED data
BEGIN { FS = OFS = "\t" }
NF == 4 { # 4-column BED
   print $4,  $1 "__" $4 "__" $2 + 1 "__" $3
   next
}
NF == 9 && $3 == "CDS" { # GFF3
    match($9, /Parent=[^;]+/)
    mRNA_ID=substr($9, RSTART+7, RLENGTH-7)
    nf = split(transcript_CDS[mRNA_ID], A, ",")
    if (nf == 0)
        transcript_CDS[mRNA_ID] = $1 "," $4 "," $5
    else {
        if ($4 < A[2]) A[2] = $4 # min start pos
        if ($5 > A[3]) A[3] = $5 # max end pos
        transcript_CDS[mRNA_ID] = A[1] "," A[2] "," A[3]
    }
}
END {
    for (mRNA_ID in transcript_CDS) { # if GFF input
        split(transcript_CDS[mRNA_ID], A, ",")
        print mRNA_ID, A[1] "__" mRNA_ID "__" A[2] "__" A[3]
    }
}

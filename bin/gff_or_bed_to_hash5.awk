#!/usr/bin/awk -f
#
# NAME
#   gff_or_bed_to_hash5.awk - Generate hash of gene ID and positional information.
#   The hash value has five fields, separated by "__": chr, geneID, start, end, orient
# 
# SYNOPSIS
#   cat FILE.[gff3|bed] | gff_or_bed_to_hash5.awk 
#     
# OPERANDS
#     INPUT_FILE
#         A file containing GFF or four- or six-column BED data
BEGIN { FS = OFS = "\t" }
NF == 4 { # 4-column BED. Code orientation as "+" in all cases.
   print $4,  $1 "__" $4 "__" $2 + 1 "__" $3 "__+"
   next
}
NF == 6 { # 6-column BED
   print $4,  $1 "__" $4 "__" $2 + 1 "__" $3 "__" $6
   next
}
NF == 9 && $3 == "CDS" { # GFF3
    match($9, /Parent=[^;]+/)
    mRNA_ID=substr($9, RSTART+7, RLENGTH-7)
    nf = split(transcript_CDS[mRNA_ID], A, ",")
    if (nf == 0)
        transcript_CDS[mRNA_ID] = $1 "," $4 "," $5 "," $7
    else {
        if ($4 < A[2]) A[2] = $4 # min start pos
        if ($5 > A[3]) A[3] = $5 # max end pos
        transcript_CDS[mRNA_ID] = A[1] "," A[2] "," A[3] "," $7
    }
}
END {
    for (mRNA_ID in transcript_CDS) { # if GFF input
        split(transcript_CDS[mRNA_ID], A, ",")
        print mRNA_ID, A[1] "__" mRNA_ID "__" A[2] "__" A[3] "__" A[4]
    }
}

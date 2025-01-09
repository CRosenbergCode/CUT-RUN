Central place for annotations to be used.

`VectorBase-67_AaegyptiLVP_AGWG.sorted.gtf` downloaded from VectorBase (requires registration) and sorted with `bedtools sort`

`VectorBase-67_AaegyptiLVP_AGWG.gff` downloaded from VectorBase (requires registration) on Thu Jan  9 11:14:06 MST 2025. No checksum listed on Vectorbase but file size is: 64996066. Number of lines: 515605.

`VectorBase-67_AaegyptiLVP_AGWG.CDS.gff` CDS entries extracted from the larger file:
`bedtools sort -i <(awk '$3 == "CDS" { print $0 }' VectorBase-67_AaegyptiLVP_AGWG.gff) > VectorBase-67_AaegyptiLVP_AGWG.CDS.gff`
Number of lines: 172905

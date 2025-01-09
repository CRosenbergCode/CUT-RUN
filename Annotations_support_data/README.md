# Annotations

## Gene annotations

### VectorBase-67_AaegyptiLVP_AGWG.gff
Downloaded from VectorBase (requires registration) on Thu Jan  9 11:14:06 MST 2025. No checksum listed on Vectorbase but file size is: 64996066. Number of lines: 515605.

### VectorBase-67_AaegyptiLVP_AGWG.CDS.gff 

CDS entries extracted from the larger file using:

`bedtools sort -i <(awk '$3 == "CDS" { print $0 }' VectorBase-67_AaegyptiLVP_AGWG.gff) > VectorBase-67_AaegyptiLVP_AGWG.CDS.gff`

Number of lines: 172905

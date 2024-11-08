Purpose- Diffbind uses raw counts output from MACS2 to compare treatment groups. Is a wrapper for DESeq2 or edgeR differential expression packages. Identifies differential peaks and overlaps between sample groups. 
Input- MACS2 output peak files, .bam files, bam indexes(.bai), .bed file outputs, Input file metadata.csv input files are downloaded locally, because Diffbind is an R package.


Output- .csv files. Use the adjusted p value, which has applied a multiple testing adjustment.

THe RNAREgulatoryRegionExtract.R file is used to pull promoter regions of interest from a list of GOIs.

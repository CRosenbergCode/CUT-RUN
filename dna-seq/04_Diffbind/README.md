Purpose- Diffbind uses raw counts output from MACS2 to compare treatment groups. Is a wrapper for DESeq2 or edgeR differential expression packages. Identifies differential peaks and overlaps between sample groups. 
Input- all peak files, .bam files, and other input files are downloaded locally, because Diffbind is an R package.
Bed file output from MACS2
Input file metadata in .csv format
Output- .csv files. Use the adjusted p value, which has applied a multiple testing adjustment.

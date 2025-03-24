Purpose- Diffbind uses raw counts output from MACS2 to compare treatment groups. Is a wrapper for DESeq2 or edgeR differential expression packages. Identifies differential peaks and overlaps between sample groups. Input- MACS2 output peak files, .bam files, bam indexes(.bai), .bed file outputs, Input file metadata.csv input files are downloaded locally, because Diffbind is an R package.

Output- .csv files. Use the adjusted p value, which has applied a multiple testing adjustment.

DiffBindFunction.R contains the code for utilizing diffbind as described above.

geneSubsetFunctions.R is used to create a gtf file containing only genes of interest (such as those differentially expressed by RNA-seq) as well as pull promoter regions of interest from a list of GOIs.

PeakAdjustmentFunctions.R is used to create a smaller set of peaks based on applying greylists to diffbind or selecting MACS2 peaks which meet the criteria for p/q value and/or fold change.

GeneDistanceShared.R takes a gtf containing a list of genes (or similar features like exons or promoter regions) and compares it to a list of peaks (generally the Diffbind output) to determine whether there are any within a given distance, or inside of the feature of interest.

runSamtools_indexLoop.sh is a bash script to sort and create index files for bam files (sorted bams are necessary for diffbind).

The RNAGeneSets directory contains example gtf files which contain commonly used sets of features such as all genes, genes as well as pseudogenes and long non-coding RNA, and exons and 3’ UTRs (common binding sites for Ago2 in Aedes aegypti).

The DiffbindSampleSheets directory contains examples of the sample sheets used by Diffbind to determine sample metadata and file locations.


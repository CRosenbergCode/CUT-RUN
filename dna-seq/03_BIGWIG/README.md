This is the only step which will require the deeptools_kernel environment, due to the conflicts deeptools can create with other packages in conda.

This pipeline does not directly utilize Bigwig files, but external software such as Integrated Genome Viewer can be used to visualize the data based on these files.

The runBamCoverage.sh script is used to create the initial bigwig files using bams. The bam files utilized must be in the same directory as the script and the bigwig outputs will be written to a directory corresponding to the bin size used.

The runBWCompare script is used to perform input subtraction to adjust sample bigwig files based on the input, to eliminate problematic genomic regions. It requires bigwigs to have been created previously using the prior script.


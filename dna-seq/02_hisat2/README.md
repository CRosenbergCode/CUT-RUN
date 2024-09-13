Purpose- Hisat2 aligns trimmed fastq reads to a reference genome or transcriptome. In this case the Aedes aegypti is aligned against an input genome.
Input- fastq.gz
Output- SAM, which is then converted to BAM 
Requires the following files: 
If no index: requires a genome fasta file named mergedGenome.fa. This should be placed in the 02_hisat2 directory.
If there are index files: they should be placed in the hisatIndex directory. Example- genomeAedesEColi.ht2.# .

These files can be found in the Rosenberg lab remote storage (RSTOR or NAS) under the "Informatics Resources" directory.
In all cases, trimmed and QCed fastq files should be in the 01_fastp/trimmedFastq directory. They will already be there if you have run the script inside the 01_fastp directory.

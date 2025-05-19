# Diffbind
**Purpose -** Determines differentially-bound sites between ChIP-seq experiments. It is a wrapper for DESeq2 and edgeR differential expression packages.

**Input -** .bed, .bam, and .bai files from MACS2 or SEACR, as well as a .RDS greylist file.

**Output -** .csv files. Be sure to use the adjusted p value.

## Usage
Diffbind is an R package, and thus must be run in R Studio. Connecting to an R Studio server is recommended. (Note for later, write more instructions on how to do this)

1. You will need to put your peak caller (MACS2/SEACR) outputs somewhere the machine running R Studio can reach them.
2. Open the `DiffBindFunction.R` script in R studio.
3. Un-comment lines 8-11 and run them. You can highlight them and click the "Run" button in R to do so. (or Ctrl+Enter) Re-comment when you are done, since they only need to be run once.
4. Change the setwd command on line 13 to point to the relevant directory for accessing your peak caller outputs. (This should be the folder where the outputs are)
5. Change the `greylist` variable to point to your .RDS greylist file.
6. 

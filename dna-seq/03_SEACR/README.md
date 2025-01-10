Purpose- SEACR is an alternative to MACS2, MACS2 is more versatile, while SEACR is more highly selective and specifically used for Cut&Run.
SEACR does not appear to work well for broader histone markers such as H3K9Me3 in mosquitoes and should not be used for this purpose at the moment.

These scripts are designed for usage with slurm on the Riviera HPC.

Before starting, the user should move all bam files they want to use to the same directory as the scripts. No additional directories must be created for the scripts.
While not invoked directly by the user, the SEACR program files (SEACR_1.3.R and SEACR_1.3.sh) must be present in the same directory. These should be manually downloaded from either
the SEACR github or our lab's CUTRUN github.

Overall process
Input- BAM
Output- BED- start/end coordinates of signal/peak

03.1_sortBams
Input- Unsorted (or differently sorted) BAM files
Output- BAM sorted by read name

03.2_
Input- BAM sorted by read name
Output- BEDGRAPH

03.3_runSEACR
Input- BEDGRAPH
Output- BED- start/end coordinates of signal/peak

Please note that SEACR requires a different sorting (by name) compared to other programs (chromosomal coordinates). It (03.1_sortBam.sh) only needs to be ran once, but must be run even if .bam files were sorted in another step.

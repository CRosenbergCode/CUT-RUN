#!/bin/bash
# ran on CCTSI server

PTH='rosenberg_feb_2024_cutnrun_rnaseq'

STAR --runThreadN 8 \
--genomeDir STAR_indices \
--readFilesIn $PTH/RNA-RVFV_Rep2_S22_L006_R1_001.fastq.gz $PTH/RNA-RVFV_Rep2_S22_L006_R2_001.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--readFilesCommand zcat

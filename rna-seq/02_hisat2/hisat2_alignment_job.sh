#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=alignmentJob
#SBATCH --output=%x.%j.out
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pegah.eizadkhah@colostate.edu

module purge
module load anaconda
conda activate cutruntools2.1

#trimmed sequence files will have to be changed for your samples
hisat2 -p 8 --rna-strandness RF --phred33 -x aedesA_tran -1 trimmed_RNA-BF_Rep1_S21_L006_R1_001.fastq.gz -2 trimmed_RNA-BF_Rep1_S21_L006_R2_001.fastq.gz -S realData/alignmentSummary_bf_rep1.sam
hisat2 -p 8 --rna-strandness RF --phred33 -x aedesA_tran -1 trimmed_RNA-BF_Rep2_S20_L006_R1_001.fastq.gz -2 trimmed_RNA-BF_Rep2_S20_L006_R2_001.fastq.gz -S realData/alignmentSummary_bf_rep2.sam
hisat2 -p 8 --rna-strandness RF --phred33 -x aedesA_tran -1 trimmed_RNA-RVFV_Rep1-1_S18_L006_R1_001.fastq.gz -2 trimmed_RNA-RVFV_Rep1-1_S18_L006_R2_001.fastq.gz -S realData/alignmentSummary_rvfv_rep1-1.sam
hisat2 -p 8 --rna-strandness RF --phred33 -x aedesA_tran -1 trimmed_RNA-RVFV_Rep1-2_S19_L006_R1_001.fastq.gz -2 trimmed_RNA-RVFV_Rep1-2_S19_L006_R2_001.fastq.gz -S realData/alignmentSummary_rvfv_rep1-2.sam
hisat2 -p 8 --rna-strandness RF --phred33 -x aedesA_tran -1 trimmed_RNA-RVFV_Rep2_S22_L006_R1_001.fastq.gz -2 trimmed_RNA-RVFV_Rep2_S22_L006_R2_001.fastq.gz -S realData/alignmentSummary_rvfv_rep2.sam

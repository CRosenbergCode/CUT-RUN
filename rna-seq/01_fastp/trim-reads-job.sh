#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=adapter_trim_job
#SBATCH --output=%x.%j.out
#SBATCH --time=6:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pegah.eizadkhah@colostate.edu

module purge
module load anaconda
conda activate cutruntools2.1

#directory for sequenced sample data: /projects/haogg@colostate.edu/RosenbergShared/datasets/CUTRUNResults/raw/20240221_LH00407_0020_B22HHNNLT3/RNA
SHARED_DRIVE='/projects/haogg@colostate.edu/RosenbergShared/datasets/CUTRUNResults/raw/20240221_LH00407_0020_B22HHNNLT3/RNA'

#trim adapters off of sequenced reads before alignment 
fastp -i $SHARED_DRIVE/RNA-BF_Rep1_S21_L006_R1_001.fastq.gz -I $SHARED_DRIVE/RNA-BF_Rep1_S21_L006_R2_001.fastq.gz  -o trimmed_RNA-BF_Rep1_S21_L006_R1_001.fastq.gz -O trimmed_RNA-BF_Rep1_S21_L006_R2_001.fastq.gz
fastp -i $SHARED_DRIVE/RNA-BF_Rep2_S20_L006_R1_001.fastq.gz -I $SHARED_DRIVE/RNA-BF_Rep2_S20_L006_R2_001.fastq.gz  -o trimmed_RNA-BF_Rep2_S20_L006_R1_001.fastq.gz -O trimmed_RNA-BF_Rep2_S20_L006_R2_001.fastq.gz
fastp -i $SHARED_DRIVE/RNA-RVFV_Rep1-1_S18_L006_R1_001.fastq.gz -I $SHARED_DRIVE/RNA-RVFV_Rep1-1_S18_L006_R2_001.fastq.gz  -o trimmed_RNA-RVFV_Rep1-1_S18_L006_R1_001.fastq.gz -O trimmed_RNA-RVFV_Rep1-1_S18_L006_R2_001.fastq.gz
fastp -i $SHARED_DRIVE/RNA-RVFV_Rep1-2_S19_L006_R1_001.fastq.gz -I $SHARED_DRIVE/RNA-RVFV_Rep1-2_S19_L006_R2_001.fastq.gz  -o trimmed_RNA-RVFV_Rep1-2_S19_L006_R1_001.fastq.gz -O trimmed_RNA-RVFV_Rep1-2_S19_L006_R2_001.fastq.gz
fastp -i $SHARED_DRIVE/RNA-RVFV_Rep2_S22_L006_R1_001.fastq.gz -I $SHARED_DRIVE/RNA-RVFV_Rep2_S22_L006_R2_001.fastq.gz  -o trimmed_RNA-RVFV_Rep2_S22_L006_R1_001.fastq.gz -O trimmed_RNA-RVFV_Rep2_S22_L006_R2_001.fastq.gz

#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=star-alignment
#SBATCH --output=%x.%j.out
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pegah.eizadkhah@colostate.edu

module load anaconda
conda activate cutruntools2.1

PTH='/projects/haogg@colostate.edu/RosenbergShared/datasets/CUTRUNResults/raw/20240221_LH00407_0020_B22HHNNLT3/RNA'

STAR --runThreadN 8 --genomeDir /scratch/alpine/c836951585@colostate.edu/DAS/STAR_indices \
--readFilesIn $PTH/RNA-BF_Rep1_S21_L006_R1_001.fastq.gz,$PTH/RNA-BF_Rep2_S20_L006_R1_001.fastq.gz,$PTH/RNA-RVFV_Rep1-1_S18_L006_R1_001.fastq.gz,$PTH/RNA-RVFV_Rep1-2_S19_L006_R1_001.fastq.gz,$PTH/RNA-RVFV_Rep2_S22_L006_R1_001.fastq.gz $PTH/RNA-BF_Rep1_S21_L006_R2_001.fastq.gz,$PTH/RNA-BF_Rep2_S20_L006_R2_001.fastq.gz,$PTH/RNA-RVFV_Rep1-1_S18_L006_R2_001.fastq.gz,$PTH/RNA-RVFV_Rep1-2_S19_L006_R2_001.fastq.gz,$PTH/RNA-RVFV_Rep2_S22_L006_R2_001.fastq.gz \
--readFilesCommand zca
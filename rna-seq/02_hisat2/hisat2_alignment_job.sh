#!/bin/bash

#SBATCH --partition=amem
#SBATCH --job-name=alignmentJob
#SBATCH --output=%x.%j.out
#SBATCH --time=6:00:00
#SBATCH --qos=mem
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pegah.eizadkhah@colostate.edu

module purge
module load anaconda
conda activate cutruntools2.1

#directory for sequenced sample data: /projects/haogg@colostate.edu/RosenbergShared/datasets/CUTRUNResults/raw/20240221_LH00407_0020_B22HHNNLT3/RNA
#sequence files will have to be changed for your samples
hisat2 -p 8 --rna-strandness RF --phred33 -x aedesA_tran -1 /projects/haogg@colostate.edu/RosenbergShared/datasets/CUTRUNResults/raw/20240221_LH00407_0020_B22HHNNLT3/RNA/RNA-BF_Rep2_S20_L006_R1_001.fastq.gz -2 /projects/haogg@colostate.edu/RosenbergShared/datasets/CUTRUNResults/raw/20240221_LH00407_0020_B22HHNNLT3/RNA/RNA-BF_Rep2_S20_L006_R2_001.fastq.gz -S realData/alignmentSummary_bf_rep2.sam

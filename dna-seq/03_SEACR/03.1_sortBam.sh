#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=sortBamSEACR
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
source activate base
conda activate rnaPseudo

#This file sorts bams based on NAME, as required by subsequent steps of the process.
#THIS MUST BE RUN EVEN IF HISAT OUTPUT WAS PREVIOUSLY SORTED, AS THE SORTING METHOD IS DIFFERENT

#Place

for FILE in *.bam
do
  EXT=${FILE::-4}
  echo $EXT
  samtools sort -n -@ 8 $FILE -o $EXT.Read.bam
done

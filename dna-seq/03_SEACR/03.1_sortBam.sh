#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=sortBamSEACR
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
module load anaconda
conda activate rnaPseudo


for FILE in ../02_hisat2/HisatAligned/*.bam
do
  EXT=${FILE:12:-4}
  echo $EXT
  #samtools sort -N -@ 8 $FILE -o $EXT.Read.bam
done
